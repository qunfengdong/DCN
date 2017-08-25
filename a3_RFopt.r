#!/usr/local/bin/Rscript


############# FUNCTION ##############

installLibrary <- function(pkg_name){
  #### Check required package and install if the package is not in the library 
  if (pkg_name %in% rownames(installed.packages()) == FALSE) {
    cat("Installing package: ",pkg_name,"\n")
    install.packages(pkg_name,repos='http://cran.us.r-project.org',quiet = T)
  }
  require(pkg_name,character.only=TRUE,quietly = T)
}

######### FUNCTION END ##############



#### Required packages ######
.cran_packages <- c("randomForestSRC","survival","optparse","reshape2","ggplot2")
sapply(c(.cran_packages),installLibrary)
#### Command line options #####
option_list = list(
  make_option(c("-i", "--infile"), type="character", default=NULL, 
              help="Input file name [required]. \n\t\tShould be an intermediate file from a1_runEMR.r", metavar="character"),
  make_option(c("-n", "--ntree"), type="integer", default=1000, 
              help="Number of trees to run random forest, default is 1000", metavar="integer"),
  make_option(c("-t", "--timepoint"),type="character",default="all",
              help="Timepoint(s) to run Wilcoxon rank sum test on, default is all. Possible choices: quantile1st, median, quantile3rd, endpoint, random, or all to run test on all 5 timepoints. You can choose more than one, separated by comma without space.",metavar = "character"),
  make_option(c("--graph"),default=FALSE,action ="store_true",
              help="Draw a simple TimePoint vs. Survival Curve for each individual based on RF estimattion. Default is false")
  
  
)

opt_parser = OptionParser(usage = "usage: Rscript %prog -i <filename> [options]",
                          option_list=option_list,
                          description = "\nThis R script will run a validation analysis on the selected disease pair using Random Forest survival analysis. Wilcoxon rank sum test will be used to test the different between exposed and non-exposed group on 6 different timepoints.")
opt = parse_args(opt_parser)

if (any(sapply(c(.cran_packages),installLibrary) == F)) {
  print_help(opt_parser)
  stop("At least one of required R packages is not available. Please install it manually.\n", call.=FALSE)
  
}

if (is.null(opt$infile)){
  print_help(opt_parser)
  stop("Please supply input file name as -i/--infile <filename>.\n", call.=FALSE)
}

if (!is.integer(opt$ntree)){
  print_help(opt_parser)
  stop("The number of trees is not integer, please correct!\n", call.=FALSE)
}

timeChoice <- c("quantile1st", "median", "quantile3rd", "endpoint","random","all")

tps <- strsplit(opt$timepoint,",")[[1]]

if (any(tps %in% timeChoice == F)){
  print_help(opt_parser)
  stop(paste0("Timepoint: ",tps[!tps %in% timeChoice]," is a invalid choice, please select from ",paste0(timeChoice,collapse = ", "),".\n"),call. = FALSE)
}

if ("all" %in% tps){
  tps = timeChoice[1:6]
}


##### Executive code #####

dat = read.csv(opt$infile, header = T)
dat$type <- factor(dat$type)


####random forest surivival analyis
basic_formula <- "Surv(time, event) ~ type + V4 + V6"
if (nlevels(factor(dat$V5)) >= 2){
  basic_formula <- paste0(basic_formula, ' + V5')
}
formula <- paste(c(basic_formula, colnames(dat)[-(1:8)]), collapse=" + ")

rfsmod <- rfsrc(as.formula(formula), data = dat, ntree = opt$ntree, tree.err=TRUE)
suvtab <- data.frame(rfsmod$survival)
qua <- quantile(1:ncol(suvtab))
colnames(suvtab)[qua] <- timeChoice[1:4]
ran <- seq(1,ncol(suvtab))[!seq(1,ncol(suvtab)) %in% qua]
colnames(suvtab)[sample(ran,1)] <- timeChoice[5]
tsttab <- suvtab[,tps,drop=F]

if (opt$graph){
  tplt <- melt(t(data.frame(rfsmod$survival)))
  colnames(tplt) <- c("timePoint","Subject","Prob")
  tplt$timePoint <- as.numeric(as.vector(gsub("^X","",tplt$timePoint)))
  exp <- factor(dat$type,labels = c("NonExposed","Exposed"))
  names(exp) <- c(1:nrow(rfsmod$survival))
  tplt$Type <- exp[match(tplt$Subject,names(exp))]
  ggplot(tplt, aes(x=timePoint, y=Prob, col=Type))+ theme_bw() + ylab("Survival Prob") +
    geom_line(aes(group=Subject)) + facet_grid(Type~.)
  ggsave(gsub(".csv",".RFsurv.png",opt$infile),width = 6.55,height=5.05)  
}


cat("\n>>Test Result:\n")
(out <- data.frame(Wil_p = sapply(lapply(tsttab,function(prob) {wilcox.test(prob~dat$type,alternative="g")}),function(x) x$p.value)))

