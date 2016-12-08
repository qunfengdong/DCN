#!/usr/local/bin/Rscript

#### Required packages ######
.cran_packages <- c("doParallel","foreach","OIsurv","optparse")

installLibrary <- function(pkg_name){
  if (pkg_name %in% rownames(installed.packages()) == FALSE) {
    install.packages(pkg_name,repos='http://cran.us.r-project.org',quiet = T)
  }
  require(pkg_name,character.only=TRUE,quietly = T)
}

sapply(c(.cran_packages),installLibrary)

##### Command line options #####
option_list = list(
  make_option(c("-i", "--infile"), type="character", default=NULL, 
              help="Input file name [required]\n\t\tShould be tab delimited with 6 columns in the order of patient_ID, Disease, Disease_Date ([%Y-%m-%d] or [%Y/%m/%d]), Age, Gender, Race. No header is needed. You could provide more information about the patient (additional columns), but only the measurement at disease A will be recorded. And these/this additional variable(s) will be used as covariate during Cox-PH regression, but won't be used as matching criteria", metavar="character"),
  make_option(c("-o","--outfile"),type="character", default="all.edges.csv",
              help="Adjacency matrix output file name, default is all.edges.csv",metavar="character"),
  make_option(c("-m","--metafile"),type="character", default=NULL,
              help="Meta file for each disease name. First column is the disease ID used in the input file column Disease, all the other column should be additional attribute of the disease, tab separated. Header is needed",metavar="character"),
  make_option(c("-a","--adjust"),type="character",default = "BH",
              help = "Method of adjusting p-values, choices are holm, hochberg, hommel, bonferroni, BY, and fdr. The default is BH",metavar = "character"),
  make_option(c("-s","--sigcut"),type="double", default=0.05,
              help="Significant P-value cutoff. Only p-values that are less than this cutoff will be considered as significant P-values, and kept in the adjacency matrix. Default is 0.05",metavar="double"),
  make_option(c("-y","--year"),type="integer",default = 5,
              help="Age difference between exposed and non-exposed, default is 5 years", metavar="integer"),
  make_option(c("-d","--days"),type="integer",default = 7,
              help="Visit date difference between exposed and non-exposed, unit days, default is 7 days", metavar="integer"),
  make_option(c("-p","--processors"),type="integer",default = 8,
              help="Number of cores/CPUs to use, default is 8",metavar="integer")
)

opt_parser = OptionParser(usage = "usage: Rscript %prog -i <filename> [options]",
                          option_list=option_list,
                          description = "\nThis R script will find the matching non-exposed population for all diesease pairs appearing in the output file from last step with some criteria, such as the matching patient should have the same gender and race as the exposed subject. After the non-exposed are matched with the exposed, Cox-PH regression will be performed taking age, race, gender as covariates (or any other confounding factors supplied by user), and final output is an adjacency matrix with all the adjusted significant hazard ratios, and survival curve graphs.")
opt = parse_args(opt_parser)

#### Check the commandline input are valid #####
if (is.null(opt$infile)){
  print_help(opt_parser)
  stop("Please supply input file name as -i/--infile <filename>.\n", call.=FALSE)
}

if (!is.integer(c(opt$year,opt$days, opt$processors))){
  print_help(opt_parser)
  stop("At least one of the integer argumnets (such as year, days and processors) is not integer, please correct.\n", call.=FALSE)
}

if (!(opt$adjust %in% p.adjust.methods)){
  print_help(opt_parser)
  stop("Please supply an existing method of adjusting P-value from holm, hochberg, hommel, bonferroni, BY, BH and fdr.\n", call.=FALSE)
}

if (opt$sigcut > 1 | opt$sigcut < 0){
  print_help(opt_parser)
  stop("The p-value cutoff you provided is either too small (<0) or too big > (1). Please correct!\n", call.=FALSE)
}

##### Read in Raw input #####
registerDoParallel(cores=opt$processors)

raw <- read.table(opt$infile,sep="\t")
raw[,3] <- as.Date(raw[,3])

## filter patients with more than 2 visits ##
pcount <- table(raw$V1)
raw <- raw[raw$V1 %in% names(pcount)[pcount > 1],]

##### Find all Files #####
# files <- Sys.glob("*.vs.*[^(case)|^(control)|^(all)].csv")
files <- Sys.glob("*.vs.*.csv")

idx <- data.frame(t(sapply(strsplit(files,"\\."),function(x) x[c(1,3)])))
idx <- data.frame(sapply(idx,function(x) as.numeric(as.vector(x))))

colnames(idx) <- c("start_D","end_D")
idx$file <- files
idx$out <- gsub("csv$","control.csv",idx$file)
idx$clean <- gsub("csv$","case.csv",idx$file)
idx$all <- gsub("csv$","all.csv",idx$file)

# remove output file if it exists #
if (file.exists(opt$outfile)) {file.remove(opt$outfile)}

# write header of the output file #
write.table(data.frame("from","to","coef","exp_coef","se_coef","z","Pr","N"),opt$outfile,sep=",",append = T,row.names = F,col.names = F,quote = F)

cat("Using ",opt$processors, " processors...\n", sep="")

foreach (f = 1:nrow(idx)) %dopar% {
  
  #cat(paste0("Processing file: ",idx$file[f],"\n"))
  
  infile <- read.csv(idx$file[f],header=F,colClasses = c("integer","integer","integer","integer","Date","integer","character","factor"))
  infile[,5] <- as.Date(infile[,5])

  if (file.exists(idx$out[f])) {file.remove(idx$out[f])}
  if (file.exists(idx$clean[f])) {file.remove(idx$clean[f])}
  if (file.exists(idx$all[f])) {file.remove(idx$all[f])}
  
  delt <- numeric()
  for (r in 1:nrow(infile)){
    tt <- infile[r,]
    # include original patient with matching age < 5y,visit time < 1wk, same race and gender
    match <- raw[ raw[,1] != tt[,1] & raw[,2] != idx$start_D[f] & as.character(raw[,6]) == as.character(tt[,8]) & abs(raw[,3] - tt[,5]) <= opt$days & abs(raw[,4] - tt[,6]) <= opt$year & as.character(raw[,5]) == as.character(tt[,7]),]
    
    if (nrow(match) >= 1){
      n_exp_patientID <- unique(match$V1)
      if (length(n_exp_patientID) == 1) {
        ctrl <- match[sample(nrow(match),1),]
      }
      else{
        ctrlpool <- match[match$V1 == sample(n_exp_patientID,1),]
        ctrl <- ctrlpool[sample(nrow(ctrlpool),1),]
      }
      tmpctrl <- raw[raw[,1] == ctrl[,1],]
      tpdis <- tmpctrl[tmpctrl[,2]==idx$end_D[f],]
      if (nrow(tpdis)>=1 && (tpdis[,3] - ctrl[,3])>0){
        dftmp <- data.frame(ID=ctrl[,1],type=0,time=as.numeric(tpdis[,3] - ctrl[,3]),event=1,ctrl[,3:ncol(ctrl)])
        write.table(dftmp,idx$out[f],sep=",",append = T,row.names = F,col.names = F,quote = F)
      }
      else{
        tm <- max(tmpctrl[,3])
        if ((tm-ctrl[,3]) > 0){
          dftmp <- data.frame(ID=ctrl[,1],type=0,time=as.numeric(tm - ctrl[,3]),event=0,ctrl[,3:ncol(ctrl)])
          write.table(dftmp,idx$out[f],sep=",",append = T,row.names = F,col.names = F,quote = F)
          
        }
        else{delt <- c(delt, r)}
      }
    }
    else{
      delt <- c(delt, r)
    }
  }
  
  if (file.exists(idx$out[f])) {
    ctrl <- read.table(idx$out[f],sep=",",colClasses = c("integer","integer","integer","integer","Date","integer","character","factor"))
    
    outfile <- infile[-delt,]
    write.table(outfile,idx$clean[f],row.names=F,col.names = F,sep=",",quote = F)
    
    al <- rbind(outfile,ctrl)
    write.table(al,idx$all[f],row.names=F,col.names = F,sep=",",quote = F)
    #cat(paste0("Output file: ", idx$all[f],"\n"))
    al$V7 <- factor(al$V7)
    
    # perform survival analysis, depends on whether gender catogory length #
    tryCatch({
      if (length(levels(al[,7])) == 1){
        vnam <- c("V2","V6",paste0("V",8:ncol(al)))
        mod <- coxph(as.formula(paste("Surv(V3, V4) ~ ",paste(vnam,collapse = "+"))), al)
        su <- summary(mod)
        sig <- su$coefficients[1,]
        tmp <- cbind(idx$start_D[f],idx$end_D[f],t(data.frame(sig)),num=nrow(al))
        write.table(tmp,opt$outfile,sep=",",append = T,row.names = F,col.names = F,quote = F)
        if (sig["exp(coef)"] >=1 & sig["Pr(>|z|)"] <= opt$sigcut){
          subobj <- survfit(with(al,Surv(V3,V4)) ~ V2,al, conf.type = "log-log")
          png(paste0(idx$start_D[f],"_",idx$end_D[f],".png"),width = 960, height = 870,res= 150)
          plot(subobj,col=c("blue","red"),lwd=3, xlab="Days",ylab=paste0("Disease ",idx$end_D[f]," Survival Prob"),
               main=paste0("Disease ",idx$start_D[f]," to Disease ",idx$end_D[f]," \nSurvival Curve"),
               sub = paste0("Exposure = Disease ",idx$start_D[f]))
          legend("bottomleft",c("Exposure=0","Exposure=1"),col=c("blue","red"),pch=15)
          dev.off()
        }
        
        #cat(paste0("Writing 1 entry to file: ",opt$outfile,"\n"))
      }
      else{
        vnam <- c("V2",paste0("V",6:ncol(al)))
        mod <- coxph(as.formula(paste("Surv(V3, V4) ~ ",paste(vnam,collapse = "+"))), al)
        su <- summary(mod)
        sig <- su$coefficients[1,]
        tmp <- cbind(idx$start_D[f],idx$end_D[f],t(data.frame(sig)),num=nrow(al))
        write.table(tmp,opt$outfile,sep=",",append = T,row.names = F,col.names = F,quote = F)
        if (sig["exp(coef)"] >=1 & sig["Pr(>|z|)"] <= opt$sigcut){
          subobj <- survfit(with(al,Surv(V3,V4)) ~ V2,al, conf.type = "log-log")
          png(paste0(idx$start_D[f],"_",idx$end_D[f],".png"),width = 960, height = 870,res= 150)
          plot(subobj,col=c("blue","red"),lwd=3, xlab="Days",ylab=paste0("Disease ",idx$end_D[f]," Survival Prob"),
               main=paste0("Disease ",idx$start_D[f]," to Disease ",idx$end_D[f]," \nSurvival Curve"),
               sub = paste0("Exposure = Disease ",idx$start_D[f]))
          legend("bottomleft",c("Exposure=0","Exposure=1"),col=c("blue","red"),pch=15)
          dev.off()
        }
        
      }
    }, error=function(e){cat("Warning: Cox-PH can not fit:",conditionMessage(e), "Skipping...\n")}
    )
    cat("Processing ",f,"/",nrow(idx),": ",idx$file[f],"\n",sep="")
    file.remove(idx$out[f])
    file.remove(idx$clean[f])
    file.remove(idx$all[f])
  }
  file.remove(idx$file[f])
}

cat("Survival Analysis finished!\n")

### multiple test correction ####
final <- read.csv(opt$outfile)
final <- final[!is.na(final$Pr),]
final$Pr_adjusted <- p.adjust(final$Pr,method = opt$adjust)
final <- final[final$Pr_adjusted <= opt$sigcut & final$exp_coef >= 1,]

### add meta attribute to output ###
if(is.null(opt$metafile)){
  write.csv(final,opt$outfile,row.names = F,quote = F)
}else{
  cat(opt$metafile, " will be incorporated in the final output file: ",opt$outfile,"\n")
  meta <- read.csv(opt$metafile,sep="\t")
  frmatr <- paste("from",colnames(meta)[2:ncol(meta)],sep="_")
  toatr <- paste("to",colnames(meta)[2:ncol(meta)],sep="_")
  
  atrmeta <- data.frame(meta[match(final$from,meta[,1]),2:ncol(meta)],meta[match(final$to,meta[,1]),2:ncol(meta)])
  colnames(atrmeta) <- c(frmatr,toatr)
  write.csv(cbind(final,atrmeta),opt$outfile,row.names = F,quote = F)
}

cat("Job Completed!\n")

proc.time()
