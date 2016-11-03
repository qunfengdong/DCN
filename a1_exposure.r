#!/usr/local/bin/Rscript

#### Required packages ######
.cran_packages <- c("doParallel","foreach","optparse")

installLibrary <- function(pkg_name){
  if (pkg_name %in% rownames(installed.packages()) == FALSE) {
    install.packages(pkg_name,repos='http://cran.us.r-project.org',quiet = T)
  }
  require(pkg_name,character.only=TRUE,quietly = T)
}

sapply(c(.cran_packages),installLibrary)

#### Command line options #####
option_list = list(
  make_option(c("-i", "--infile"), type="character", default=NULL, 
              help="Input file name [required]. \n\t\tShould be tab delimited with 6 columns in the order of patient_ID, Disease, Disease_Date ([%Y-%m-%d] or [%Y/%m/%d]), Age, Gender, Race. No header is needed. You could provide more information about the patient (additional columns), but only the measurement at disease A will be recorded", metavar="character"),
  make_option(c("-p","--processors"),type="integer",default = 8,
              help="Number of cores/CPUs to use, default is 8",metavar="integer")
)

opt_parser = OptionParser(usage = "usage: Rscript %prog -i <filename> [options]",
                          option_list=option_list,
                          description = "\nThis R script will find the exposure population for all diesease pairs appearing in the input file that satisfy disease A -> disease B, and outputs CSV files of selected exposed patients for each pair of disease combinations. The output from this step will be the input for next step.")
opt = parse_args(opt_parser)

if (is.null(opt$infile)){
  print_help(opt_parser)
  stop("Please supply input file name as -i/--infile <filename>.\n", call.=FALSE)
}

if (!is.integer(opt$processors)){
  print_help(opt_parser)
  stop("Please correct -p/--processors argument to an integer.\n", call.=FALSE)
}


##### Read in Input ########
raw <- read.table(opt$infile,sep="\t")
raw[,3] <- as.Date(raw[,3])

##### Filter patient entry with more than 2 visits ####
pcount <- table(raw$V1)
raw <- raw[raw$V1 %in% names(pcount)[pcount > 1],]

registerDoParallel(cores=opt$processors)
idx <- combn(unique(raw[,2]),2)
idx <- cbind(idx,idx[c(2,1),])

cat("Using ",opt$processors, " processors...\n", sep="")

foreach (i = 1:ncol(idx)) %dopar% {
  if (file.exists(paste0(idx[1,i],".vs.",idx[2,i],".csv"))) {file.remove(paste0(idx[1,i],".vs.",idx[2,i],".csv"))} 

  tmpall <- raw[raw[,1] %in% raw[,1][raw[,2] == idx[1,i]],]

  ids <- unique(tmpall[,1])
  nnr <- length(ids)
  
  for (d in 1:nnr){
    id = ids[d]
    pre <- tmpall[tmpall[,1] ==id & tmpall[,2] == idx[1,i],]
    post <- tmpall[tmpall[,1] == id & tmpall[,2] == idx[2,i],]
    
    if (nrow(post) == 1 && (post[,3] - pre[,3]) > 0){
      dftmp <- data.frame(ID=id,type=1,time=as.numeric(post[,3] - pre[,3]),event=1,pre[,3:ncol(pre)])
      write.table(dftmp,paste0(idx[1,i],".vs.",idx[2,i],".csv"),sep=",",append = T,row.names = F,col.names = F,quote = FALSE)
    }else {
      ## calculate the longest duration without disease ##
      tm <- max(tmpall[,3][tmpall[,1] == id])
      if (as.numeric(tm - pre[,3]) >0){
        dftmp <- data.frame(ID=id,type=1,time=as.numeric(tm - pre[,3]),event=0,pre[,3:ncol(pre)])
        write.table(dftmp,paste0(idx[1,i],".vs.",idx[2,i],".csv"),sep=",",append = T,row.names = F,col.names = F,quote = FALSE)
      }
    }
  }
}

cat("Job Completed!\n")

proc.time()
