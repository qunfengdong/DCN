#!/usr/local/bin/Rscript

# Combining the original a1_exposure.r and a2_non_exposed.cox.adjust.r

############# FUNCTION ##############

installLibrary <- function(pkg_name){
  #### Check required package and install if the package is not in the library 
  if (pkg_name %in% rownames(installed.packages()) == FALSE) {
    install.packages(pkg_name, repos='http://cran.us.r-project.org', quiet=TRUE)
  }
  require(pkg_name, character.only=TRUE, quietly=TRUE)
}

FindMatch <- function(case, A, infile, opt) {
  ind <- infile[.(c(case[,V3]-opt$days-1, case[,V3]+opt$days)), 
                which=TRUE, roll=TRUE, mult="last"]
  if (is.na(ind[1])) ind[1] <- 0
  recent_visit <- infile[(ind[1]+1):ind[2]]
  recent_visit[, row_index := (ind[1]+1):ind[2]]
  if (recent_visit[,.N] > 1e4) {
    setkey(recent_visit, V5, V6)
    recent_visit <- recent_visit[.(case[,V5], case[,V6])]
  } else {
    recent_visit <- recent_visit[V5==case[,V5] & V6==case[,V6], ]
  }
  match <- recent_visit[V1!=case[,V1] & V2!=A & abs(V4 - case[,V4])<=opt$year, ]
  if (match[,.N] < 1) return(as.integer(NA))
  return(match[match[, list(tmpidx=.I[sample(.N, 1)]), by=V1][sample(.N, 1), tmpidx], row_index])
  #return(match[match[, .I[sample(.N,1)], by=V1][sample(.N, 1), V1], row_index]) # if want to sample based on patient?
}

FindDischarge <- function(DT, ids, B, type, infile, opt){
  DT[, VB := infile[.(ids, B), V3, on=c("V1", "V2")]]
  DT[, event := !is.na(VB) & VB>V3 & VB-V3<=opt$duration*365.24]
  DT[, time := pmin(ifelse(event, VB-V3, Vmax-V3), round(opt$duration*365.24))]
  DT[, pair := 1:nrow(DT)]
  DT[, c("type", "event", "Vmax", "VB") := list(type, event+0, NULL, NULL)]
  setnames(DT, "V1", "ID")
  setcolorder(DT, c("ID", "type", "time", "event","pair", names(DT)[2:(ncol(DT)-4)]))
  return(DT)
}

generate_cases <- function(i, disease, raw, opt) {
  abort_list <- rep(TRUE, length(disease))
  abort_list[i] <- FALSE
  A <- disease[i]
  nameA <- opt$metafile[match(A,opt$metafile[,1]),2]
  cat("Processing disease ", nameA , ".....\n", sep="")
  
  #### Extract exposed group for disease A based on numbers of cases
  # Add unique if duplicate records occured on same disease for same patient
  exposed_ids <- raw[.(A), V1, on="V2"]
  if (length(exposed_ids) < opt$minn) {
    warning(paste0("Patient with Disease ", nameA, " has less than ", 
                   opt$minn, " entries! Abort!!"))
    return(disease[abort_list])
  } else if (length(exposed_ids) > opt$maxn){
    exposed_ids <- sample(exposed_ids, opt$maxn)
  }
  exposed <- raw[.(exposed_ids, A), on=c("V1", "V2")]
  exposed[, Vmax:=raw[.(exposed_ids), .(tmp=max(V3)), by=.EACHI, on="V1"][,tmp]]
  
  # Extract control group for disease A based for each patient in exposed group
  control <- raw[exposed[, FindMatch(.SD, A=A, infile=raw, opt=opt), 
                         by=seq_len(nrow(exposed))][,V1]]
  control_ids <- control[,V1]
  control[, Vmax:=raw[.(control_ids), .(tmp=max(V3)), by=.EACHI, on="V1"][,tmp]]
  
  # Process whether patient has disease B in exposed and control group
  for (j in c(1:length(disease))[-i]){
    B <- disease[j]
    nameB = opt$metafile[match(B,opt$metafile[,1]),2]
    # cat("Processing disease ", B , " under A .....\n", sep="")
    # Whether patient in exposed group has B, type = 1
    exposed_DT <- exposed[, -c("V2"), with=FALSE]
    exposed_DT <- FindDischarge(exposed_DT, exposed_ids, B, 1, raw, opt)
    
    # Whether patient in control group has B, type = 0
    control_DT <- control[, -c("V2"), with=FALSE]
    control_DT <- FindDischarge(control_DT, control_ids, B, 0, raw, opt)
    
    # Remove NAs in control group and delete matched cases in exposed group
    bool_tag <- !is.na(control_DT[,ID])
    exposed_DT <- exposed_DT[bool_tag, ]
    control_DT <- control_DT[bool_tag, ]
    
    # Check whether exposed group has enough discharge cases
    if (exposed_DT[event==1, .N] < opt$minCaseNumber) {
      warning(paste0("Disease ", nameA, " to ", nameB, " has less than ", 
                     opt$minCaseNumber, " discharge cases! Abort!!"))
      next
    }
    # Output exposed and control group to file
    if (opt$saveT0) {
      write.csv(rbindlist(list(exposed_DT, control_DT)), 
                file.path(opt$outfolder, "T0_include", paste0(A, "_", B, ".csv")),
                row.names=FALSE, quote=FALSE)
    }
    
    
    # Remove duration == 0 from exposed and control group
    
    bool_tag <- control_DT[,time] > 0 & exposed_DT[,time] > 0
    exposed_DT <- exposed_DT[bool_tag, ]
    control_DT <- control_DT[bool_tag, ]
    write.csv(rbindlist(list(exposed_DT, control_DT)), 
              file.path(opt$outfolder, paste0(A, "_", B, ".csv")),
              row.names=FALSE, quote=FALSE)
    abort_list[j] <- FALSE
  }
  return(disease[abort_list])
}
.generate_cases_par <- function(i) generate_cases(i, disease, raw, opt)
######### FUNCTION END ##############

#### Required packages ######
.cran_packages <- c("data.table", "doParallel", "foreach", "optparse")
sapply(c(.cran_packages), installLibrary)

#### Command line options #####
option_list = list(
  make_option(c("-i", "--infile"), type="character", default=NULL, 
              help="Input file name [required]. \n\t\tShould be tab delimited with 6 columns in the order of patient_ID, Disease, Disease_Date ([%Y-%m-%d] or [%Y/%m/%d]), Age, Gender, Race. No header is needed. You could provide more information about the patient (additional columns), but only the measurement at disease A will be recorded", metavar="character"),
  make_option(c("-m", "--metafile"), type="character", default=NULL,
              help="Meta file for each disease name [required]. \n\t\tFirst column is the disease ID used in the input file column Disease, second column is the disease name. All the other columns should be additional attribute of the disease, tab separated. Header is needed", metavar="character"),
  make_option(c("--outfolder"), type="character", default=".",
              help="Intermediate files output folder, default is current directory", metavar="character"),
  make_option(c("--maxn"), type="integer", default=10000,
              help="Maximum number of subjects to be qualified as valid exposed cohort, default is 10000. If the total number of patients exceed this number, a subset of 10000 subjects will be randomly selected.", metavar="integer"),
  make_option(c("--minn"), type="integer", default=100,
              help="Minimum number of subjects to be qualified as valid exposed cohort, default is 100", metavar="integer"),
  make_option(c("-c", "--minCaseNumber"), type="integer", default=20,
              help="Minimum number of A->B cases to keep disease pair, default is 20", metavar="integer"),
  make_option(c("--duration"), type="integer", default=5,
              help="study duration (years), default is 5 years", metavar="integer"),
  make_option(c("-y", "--year"), type="integer", default=5,
              help="Age difference between exposed and non-exposed, default is 5 years", metavar="integer"),
  make_option(c("-d", "--days"), type="integer", default=7,
              help="Visit date difference between exposed and non-exposed, unit: days, default is 7 days", metavar="integer"),
  make_option(c("--saveT0"), type="logical",action = "store_true", default=FALSE,
              help="A flag to save the cohort tables with time-to-event = 0. The output CSVs will be saved to T0_include folder, default is false"),
  make_option(c("-p", "--processors"), type="integer", default=8,
              help="Number of cores/CPUs to use, default is 8", metavar="integer")
)

opt_parser = OptionParser(usage = "usage: Rscript %prog -i <inputFileName> -m <metaFileName> [options]",
                          option_list=option_list,
                          description = "\n\t\t >> This the step 1 of EMR package. << \n\tThis R script will find the exposure population for all diesease pairs appearing in the input file that satisfy disease A -> disease B, and and second, it will find the matching non-exposed population for all diesease pairs with some criteria, such as the matching patient should have the same gender and race as a good match control. Outputs are CSV files with paired exposed and non-exposed subjects information. One folder named T0_include with CSV files with time-to-event = 0 could be produced if you turn on the --saveT0 flag. Those files will be used as inputs for step 2 of the package.\n")
opt = parse_args(opt_parser)

if (any(sapply(c(.cran_packages),installLibrary) == F)) {
  print_help(opt_parser)
  stop("At least one of required R packages is not available. Please install it manually.\n", call.=FALSE)
}

if (is.null(opt$infile)){
  print_help(opt_parser)
  stop("Please supply input file name as -i/--infile <filename>.\n", call.=FALSE)
}

if (is.null(opt$metafile)){
  print_help(opt_parser)
  stop("Please supply meta file name as -m/--metafile <filename>.\n", call.=FALSE)
}

opt$metafile <- read.csv(opt$metafile, sep="\t", stringsAsFactors=FALSE)

if (!is.integer(c(opt$maxn,opt$minn,opt$year,opt$days,opt$minCaseNumber,opt$processors,opt$duration))){
  print_help(opt_parser)
  stop("At least one of the integer argumnets (maxn, minn, year, days, minCaseNumber,duration of study and processors) is not integer, please correct!\n", call.=FALSE)
}

if (!file.exists(opt$outfolder)) {
  dir.create(opt$outfolder)
  cat("Creating output folder:", file.path(getwd(), opt$outfolder), "\n")
}

## Add subfolder for different criteria
if (opt$saveT0 & !file.exists(file.path(opt$outfolder, "T0_include")) ) {
  dir.create(file.path(opt$outfolder, "T0_include"))
  cat("Creating sub-folder to save alternative results:", file.path(getwd(), opt$outfolder, "T0_include"), "\n")
}

# if (!file.exists(file.path(opt$outfolder, "D0_exclude"))) {
#   dir.create(file.path(opt$outfolder, "D0_exclude"))
#   cat("Creating sub-folder:", file.path(getwd(), opt$outfolder, "D0_exclude"), "\n")
# }

######### PARSER END ##############

##### Read in Input ########
raw <- fread(opt$infile, sep="\t", stringsAsFactors=FALSE)
raw[, V3 := as.Date(V3)]
##### Filter patient entry with more than 1 visits ####
raw <- raw[V1 %in% raw[,.(visit=uniqueN(V3)), by=V1][visit>1, V1]]

##### Set up Keys and Indices for binary search #####
setkey(raw, V3)
setindex(raw, V1)
setindex(raw, V2)
setindex(raw, V1, V2)

##### Register how many cores to use #####
registerDoParallel(cores=opt$processors)
cat("Using ", getDoParWorkers(), " processors...\n", sep="")

##### Get all diseases ####
disease <- raw[, unique(V2)]

# abort_list <- foreach(i=1:length(disease), .errorhandling="remove") %do% generate_cases(i, disease, raw, opt)

abort_list <- foreach(i=1:length(disease), .errorhandling="remove") %dopar% generate_cases(i, disease, raw, opt)

##### Get aborted disease list #####
names(abort_list) <- disease

# abort_list <- lapply(1:length(disease), function(i) {
#   if (length(abort_list[[i]]) > 0) {
#     cbind(disease[i], abort_list[[i]])
#   } else {
#     NULL
#   }
# })
# abort_list <- do.call('rbind', abort_list)
# 

abort_list <- unlist(lapply(disease, function(x) {
  frmlab = opt$metafile[match(x,opt$metafile[,1]),2]
  x <- toString(x)
  if (length(abort_list[[x]]) > 0) {
    tolab = opt$metafile[match(abort_list[[x]],opt$metafile[,1]),2]
    paste(frmlab, tolab, sep=" to ")
  } else {
    NULL
  }}))

if (length(abort_list) > 0) {
  cat(length(abort_list), "disease pairs are aborted: \n")
  cat(abort_list, sep=", ")
}

cat("\nJob Completed!\n")

proc.time()
