#!/usr/local/bin/Rscript

# Combining the original a1_exposure.r and a2_non_exposed.cox.adjust.r

############# FUNCTION ##############

installLibrary <- function(pkg_name){
  #### Check required package and install if the package is not in the library 
  if (pkg_name %in% rownames(installed.packages()) == FALSE) {
    install.packages(pkg_name,repos='http://cran.us.r-project.org',quiet = T)
  }
  require(pkg_name,character.only=TRUE,quietly = T)
}

List2Df <- function(x) {
  #### Format list of entries to a data frame 
  do.call(rbind,x)
}

RunCoxPH <- function(RL,names,index,sigcut,expcoefcut,metafile) {
  #### Run Cox-PH regression 
  label <- names[[index]]
  dA <- strsplit(label,".vs.")[[1]][1]
  dB <- strsplit(label,".vs.")[[1]][2]
  dAn = metafile[metafile[,1] == dA,2]
  dBn = metafile[metafile[,1] == dB,2]
  #   dA <- strsplit(label,"-")[[1]][1]
  cat(paste0("Processing ",dA," vs ",dB," .... \n"))
  outdf <- RL[[label]]
  basic_formula <- "Surv(time, event) ~ type + V4 + V6"
  if (nlevels(factor(outdf$V5)) >= 2){
    basic_formula <- paste0(basic_formula, ' + V5')
  }
  formula <- paste(c(basic_formula, colnames(outdf)[-(1:8)]), collapse=" + ")
  tryCatch({
    model <- coxph(as.formula(formula), data=outdf)
    testHR <- cox.zph(model) 
    HRp <- testHR$table["GLOBAL","p"]
    sig <- summary(model)$coef[1,]
    out <- data.frame(from=dA,to=dB,t(sig), nrow(outdf),HRp)
    
    if (abs(sig["exp(coef)"]) >= expcoefcut & sig["Pr(>|z|)"] <= sigcut){
      subobj <- survfit(with(outdf,Surv(time,event)) ~ type,outdf, conf.type = "log-log")
      autoplot(subobj, censor.shape = '*') + theme_light() +
        ylab(paste0(dBn," Survival Prob")) + xlab(paste0("Time (Days)\n Exposure = ",dAn)) +
#        guides(fill=guide_legend(title = "Exposure"))
        ggtitle(paste0(dAn, " vs. ", dBn," Survival Curve")) +
        theme(plot.title = element_text(hjust = 0.5))
      ggsave(paste0(gsub(".vs.","_",label),".png"),width = 6.5,height=6.15)
      
      # png(paste0(label,".png"),width = 960, height = 870,res= 150)
      # plot(subobj,col=c("blue","red"),lwd=3, xlab="Days",ylab=paste0(label," Survival Prob"),
      #      main=paste0(label," Survival Curve"),
      #      sub = paste0("Exposure = Disease ",dA))
      # legend("bottomleft",c("Exposure=0","Exposure=1"),col=c("blue","red"),pch=15)
      # dev.off()
      
    }
  },
  error=function(e) {out <<- data.frame(from=NULL,to=NULL,coef=NULL, exp.coef=NULL,se.coef=NULL,z=NULL, Pr...z..=NULL, nrow.outdf.=NULL)})
  out
}

WriteInter <- function(RL,names=names(RL),index){
  ### write out intermediate data file used for survival analysis
  label <- names[[index]]
  outdf <- RL[[label]]
  
  pairfile <- paste0(label,".csv")
  if (file.exists(pairfile)) {file.remove(pairfile)} 
  write.csv(outdf,pairfile,row.names = F,quote = FALSE)
  
}

FilterExposure <- function(outdf,min) {
  #### Filter exposed entries with a minimum number of cases
  sum(outdf$type == 1 & outdf$event == 1) >= min
  
}

FindMatch <- function(case,infile, diseaseA,dDiff=opt$days,yDiff=opt$year){
  #### Find the non-exposed match for each case
  matchbool <- infile[,1] != case[,1] & infile[,2] != diseaseA & 
    as.character(infile[,6]) == as.character(case[,6]) & abs(infile[,3] - case[,3]) <= dDiff & 
    abs(infile[,4] - case[,4]) <= yDiff & as.character(infile[,5]) == as.character(case[,5])
  sumMbool <- sum(matchbool)
  if (sumMbool >= 1) {
    match <- infile[matchbool,]
    n_exp_patientID <- unique(match$V1)
    if (length(n_exp_patientID) == 1){
      randid = n_exp_patientID
    }else{
      randid <- sample(n_exp_patientID,1)
    }
    ne_idx <- sample(1:sum(match$V1 == randid),1)
    match[match$V1 == randid,][ne_idx,]
  }
  
}

RmExCases <- function(outdf,validMatchIDs){
  outdf[outdf$ID %in% validMatchIDs,]
}

GetDiseaeB <- function(FERL,Fnames,index) {
  label <- Fnames[[index]]
  strsplit(label,".vs.")[[1]][2]
  
}
 
TestDiseaseB <- function(entry,infile,diseaseA,diseaseB,sDu=opt$duration){
  
  control_history <- infile[infile$V1 == entry[,1],]
  outlist <- list()
  
  for (B in diseaseB){
    control_B <- control_history[control_history$V2 == B,]
    eventB_bool <- nrow(control_B) >= 1 && as.numeric(control_B$V3 - entry$V3) > 0
    durationB <- ifelse(eventB_bool, control_B$V3, max(control_history$V3)) - as.numeric(entry$V3)
    if (durationB > sDu * 365.24){
      ctl <- data.frame(ID=entry[,1], type=0, time=round(sDu*365.24), event=0, entry[,-c(1,2)])
    }else{
      ctl <- data.frame(ID=entry[,1], type=0, time=durationB, event=as.numeric(eventB_bool), entry[,-c(1,2)])
    }
    outlist[[paste0(diseaseA,".vs.",B)]] <- ctl
  }
  outlist
}

######### FUNCTION END ##############



#### Required packages ######
.cran_packages <- c("doParallel","foreach","OIsurv","optparse","ggplot2","ggfortify")

sapply(c(.cran_packages),installLibrary)

#### Command line options #####
option_list = list(
  make_option(c("-i", "--infile"), type="character", default=NULL, 
              help="Input file name [required]. \n\t\tShould be tab delimited with 6 columns in the order of patient_ID, Disease, Disease_Date ([%Y-%m-%d] or [%Y/%m/%d]), Age, Gender, Race. No header is needed. You could provide more information about the patient (additional columns), but only the measurement at disease A will be recorded", metavar="character"),
  make_option(c("-m","--metafile"),type="character", default=NULL,
              help="Meta file for each disease name [required]. \n\t\tFirst column is the disease ID used in the input file column Disease, second column is the disease name. All the other columns should be additional attribute of the disease, tab separated. Header is needed",metavar="character"),
  make_option(c("-o","--outfile"),type="character", default="all.edges.csv",
              help="Adjacency matrix output file name, default is all.edges.csv",metavar="character"),
  make_option(c("--outfolder"),type="character", default=".",
              help="Intermediate files and survival PNG output folder, default is current directory",metavar="character"),
  make_option(c("--maxn"),type="integer",default = 5000,
              help="Maximum number of subjects to be qualified as valid exposed cohort, default is 5000. If the total number of patients exceed this number, a subset of 5000 subjects will be randomly selected.",metavar="integer"),
  make_option(c("--minn"),type="integer",default = 100,
              help="Minimum number of subjects to be qualified as valid exposed cohort, default is 100",metavar="integer"),
  make_option(c("-c","--caseNumber"),type="integer",default=20,
              help="Minimum number of A->B cases to keep disease pair, default is 20",metavar="integer"),
  make_option(c("--duration"),type="integer", default=10,
              help="study duration (years), default is 10 years",metavar="integer"),
  make_option(c("-y","--year"),type="integer",default = 5,
              help="Age difference between exposed and non-exposed, default is 5 years", metavar="integer"),
  make_option(c("-d","--days"),type="integer",default = 7,
              help="Visit date difference between exposed and non-exposed, unit: days, default is 7 days", metavar="integer"),
  make_option(c("-s","--sigcut"),type="double", default=0.05,
              help="Significant P-value cutoff. Only p-values that are less than this cutoff will be considered as significant P-values, and kept in the adjacency matrix. Default is 0.05",metavar="double"),
  make_option(c("-e","--exp_coef"),type="double", default=1,
              help="Exponentiated coefficients cutoff used in plotting survival curve. Default is 1",metavar="double"),
  make_option(c("-a","--adjust"),type="character",default = "BH",
              help = "Method of adjusting p-values, choices are holm, hochberg, hommel, bonferroni, BY, and fdr. The default is BH",metavar = "character"),
  make_option(c("-r","--writeintermediates"),default=FALSE,action="store_true",
              help="Write out individual disease pair survival tables, default is FALSE"),
  make_option(c("-p","--processors"),type="integer",default = 8,
              help="Number of cores/CPUs to use, default is 8",metavar="integer")
)

opt_parser = OptionParser(usage = "usage: Rscript %prog -i <inputFileName> -m <metaFileName> [options]",
                          option_list=option_list,
                          description = "\nThis R script will first find the exposure population for all diesease pairs appearing in the input file that satisfy disease A -> disease B, and and second, it will find the matching non-exposed population for all diesease pairs with some criteria, such as the matching patient should have the same gender and race as control. After the non-exposed are matched with the exposed group, Cox-PH regression will be performed taking age, race, gender as covariates (or any other confounding factors supplied by user), and final output is an adjacency matrix with all the significant hazard ratios, and survival curve graphs. Multiple test correction will be applied to the p-values.")
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


if (!(opt$adjust %in% p.adjust.methods)){
  print_help(opt_parser)
  stop("Please supply an existing method of adjusting P-value from holm, hochberg, hommel, bonferroni, BY, BH and fdr.\n", call.=FALSE)
}

if (!is.integer(c(opt$maxn,opt$minn,opt$year,opt$days, opt$caseNumber,opt$processors,opt$duration))){
  print_help(opt_parser)
  stop("At least one of the integer argumnets (maxn, minn, year, days, caseNumber,duration of study and processors) is not integer, please correct!\n", call.=FALSE)
}

if (!is.logical(opt$writeintermediates)){
  print_help(opt_parser)
  stop("Please enter either TRUE or FALSE for -r, --writeintermediates!\n",call. = FALSE)
}


##### Read in Input ########
raw <- read.table(opt$infile,sep="\t")
raw[,3] <- as.Date(raw[,3])

meta <- read.csv(opt$metafile,sep="\t")
##### Filter patient entry with more than 2 visits ####
pcount <- table(raw$V1)
raw <- raw[raw$V1 %in% names(pcount)[pcount > 1],]

if (file.exists(opt$outfolder)){
  setwd(file.path(getwd(), opt$outfolder))
} else {
  cat("Creating output folder!\n")
  dir.create(file.path(getwd(), opt$outfolder))
  setwd(file.path(getwd(), opt$outfolder))
}

##### Register how many cores to use #####
# cl = makeCluster(opt$processors)
# registerDoParallel(cl)
# options(cores=opt$processors)
# registerDoMC(opt$processors)
registerDoParallel(cores=opt$processors)
# DiseaseSubset <- lapply(split(raw,raw$V2),function(x) raw[raw$V1 %in% x$V1,])

##### Collect all diseases ####
disease <- unique(raw$V2)

cat("Using ",getDoParWorkers(), " processors...\n", sep="")

if (file.exists(opt$outfile)) {file.remove(opt$outfile)}
write.table(data.frame("from","to","coef","exp_coef","se_coef","z","Pr","N","HRtest-p"),opt$outfile,sep=",",append = T,row.names = F,col.names = F,quote = F)

foreach (i = 1:length(disease),.errorhandling="pass") %dopar% {
# for (i in 1:length(disease)){
  A <- disease[i]
  cat("Processing disease ", A , ".....\n", sep="")
  tmpall <- raw[raw[,1] %in% raw[,1][raw[,2] == A],]
  
  #### Get all disease after A ####
  nonA <- disease[-i]
  tmpids <- unique(tmpall[,1])
  nnr <- length(tmpids)
  exposedResult <- list()

  #### Test to see how many numbers of cases could find ####
  if (nnr <= opt$maxn & nnr >= opt$minn){
    ids <- tmpids
  }
  if (nnr > opt$maxn){
    ids <- sample(tmpids,opt$maxn)
  }
  if (nnr < opt$minn) stop(paste0("Patient with Disease ",A, " has less than ",opt$minn, " entries! Abort!!"))

  for (id in ids){
    for (B in nonA){
      pre <- tmpall[tmpall[,1] ==id & tmpall[,2] == A,]
      post <- tmpall[tmpall[,1] == id & tmpall[,2] == B,]
      
      eventA_bool <- nrow(post) == 1 && (post[,3] - pre[,3]) > 0
      durationA <- ifelse(eventA_bool, post[,3], max(tmpall[tmpall[,1] == id,3])) - as.numeric(pre[,3])
      if (durationA < 0) next 
      if (durationA > opt$duration * 356.24) {
        case <- data.frame(ID=id,type=1,time=round(opt$duration*365.24),event=0,pre[,3:ncol(pre)])
      }else {
        case <- data.frame(ID=id,type=1,time=durationA,event=as.numeric(eventA_bool),pre[,3:ncol(pre)])
      }
      
      
      exposedResult[[paste0(A,".vs.",B)]][[length(exposedResult[[paste0(A,".vs.",B)]])+1]] <- case
    }
  }
  
  ###### Format all exposed entries and filter rare disease pair #######
  exposedResultList <- lapply(exposedResult,List2Df)
  exposedResultList <- exposedResultList[sapply(exposedResultList,FilterExposure,min=opt$caseNumber)]
    
  ###### Find out real exposed IDs ######  
  if (identical(exposedResultList[[1]][["ID"]],ids)){
    validExIDs <- ids
  }else{
    validExIDs <- exposedResultList[[1]][["ID"]]
  }
    
  validList <- split(tmpall[tmpall[,1] %in% validExIDs & tmpall[,2] == A,],validExIDs)
     
  ###### Find matching nonexposed entry ######
  validMatch <- lapply(validList,FindMatch,diseaseA=A, infile=raw)
  
  validMatchIDs <- validExIDs[sapply(validMatch,function(x) !is.null(x))]
    
  FilteredExposedResultList <- lapply(exposedResultList,RmExCases,validMatchIDs=validMatchIDs)

  validMatch <- Filter(Negate(is.null),validMatch)
    
  sm_nonA <- sapply(seq_along(FilteredExposedResultList),GetDiseaeB,
                    FERL=FilteredExposedResultList,Fname=names(FilteredExposedResultList))

  nonexResult <- lapply(validMatch,TestDiseaseB,infile=raw, diseaseA=A,diseaseB=sm_nonA,sDu=opt$duration)
  
  ##### merge Exposed and nonExposed together #####
  ResultList <- list()
  for (na in names(FilteredExposedResultList)){
      ResultList[[na]] <- rbind(do.call(rbind,lapply(nonexResult,function(x) x[[na]])),FilteredExposedResultList[[na]])
  }
  
  ###### Write out intermediate files #####
  if (opt$writeintermediates){
    lapply(seq_along(ResultList),WriteInter,RL=ResultList,names=names(ResultList))
  }
  ###### Output A 2 NonA to all.edges.csv ####
  outResult <- lapply(seq_along(ResultList),RunCoxPH,RL=ResultList,names=names(ResultList),sigcut=opt$sigcut,expcoefcut=opt$exp_coef,metafile=meta)

  write.table(do.call(rbind,outResult),opt$outfile,sep=",",append = T,row.names = F,col.names = F,quote = F)

}

# stopCluster(cl)
final <- read.csv(opt$outfile)
## remove NAs ##

### multiple test correction ####
final$Pr_adjusted <- p.adjust(final$Pr,method = opt$adjust)
final <- final[final$Pr_adjusted <= opt$sigcut & final$exp_coef >= 1,]

### add meta attribute to output ###
frmatr <- paste("from",colnames(meta)[2:ncol(meta)],sep="_")
toatr <- paste("to",colnames(meta)[2:ncol(meta)],sep="_")

atrmeta <- data.frame(meta[match(final$from,meta[,1]),2:ncol(meta)],meta[match(final$to,meta[,1]),2:ncol(meta)])
colnames(atrmeta) <- c(frmatr,toatr)
write.csv(cbind(final,atrmeta),opt$outfile,row.names = F,quote = F)


cat("Job Completed!")

proc.time()
