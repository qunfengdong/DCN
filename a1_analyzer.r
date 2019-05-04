#!/usr/local/bin/Rscript

############# FUNCTION ##############
installLibrary <- function(pkg_name){
  #### Check required package and install if the package is not in the library 
  if (pkg_name %in% rownames(installed.packages()) == FALSE) {
    install.packages(pkg_name, repos='http://cran.us.r-project.org', quiet=TRUE)
  }
  require(pkg_name, character.only=TRUE, quietly=TRUE)
}

file_parser <- function(filenames, ...) {
  pattern <- "(?P<A>^[A-Za-z0-9]+?)_(?P<B>[A-Za-z0-9]+?)\\.csv$"
  m <- regexpr(pattern, filenames, perl=TRUE)
  file_lists <- do.call(rbind, lapply(seq_along(m), function(i) {
    if (m[i] == -1) return(NULL)
    st <- attr(m, "capture.start")[i, ]
    lt <- attr(m, "capture.length")[i, ]
    c(substring(filenames[i], st, st+lt-1))
  }))
  if (is.null(file_lists)) {
    stop("Your input directory does not have any valid files. Please provide the right input folder to --inputfolder argument!\n")
  }else{
    colnames(file_lists) <- attr(m, "capture.names")
    file_lists <- cbind(filename=filenames[m!=-1], file_lists)
    return(file_lists)
  }
  
}

.get_formula <- function(x) {
  basic_formula <- "Surv(time, event) ~ type + V4 + V6"
  if (x[,uniqueN(V5)] > 1) basic_formula <- paste0(basic_formula, ' + V5')
  return(paste(c(basic_formula, colnames(x)[-(1:9)]), collapse=" + "))
}

plotStats <- function(outdf,dAn, dBn, label, opt) {
  ##### study duration #####
  du = round(opt$duration*365.24)
  ##### generate table of frequency #####
  tb <- data.frame(dropoff = sapply(split(outdf,outdf$type),function(z) sum(z$event == 0 & z$time != du)),
                   event = sapply(split(outdf,outdf$type),function(z) sum(z$event == 1 )),
                   survived = sapply(split(outdf,outdf$type),function(z) sum(z$event == 0 & z$time == du)))
  tb <- t(tb)
  rownames(tb) <- c("Total Drop-off Frequency","Event Frequency","Survived Frequency")
  colnames(tb) <- c("Non-Exposed","Exposed")
  
  #### subset the input csv to only dropoff patients ####
  dropset = subset(outdf,event == 0 & time != du)
  
  #### calculate ks p-value ####
  ksp = ks.test(dropset[type==1,time],dropset[type==0,time])$p.value

  dropset$type <- factor(dropset$type,labels = c("Non-Exposed","Exposed"))
  cat("  Save drop-off distribution plot:", paste0(label, ".dropoff.png\n"))
  gg_graph <- ggplot(aes(time),data=dropset) + 
#    geom_histogram(aes(y=..density..),colour="black", fill="white") + 
#    geom_histogram(colour="black", fill="white") +
#    geom_freqpoly() + 
    geom_density(alpha=.2, fill="#FF6666") +
    ylab("Density") + theme_bw() + facet_grid(.~type) + 
    xlab(paste0("Time (Days)\nExposure = ",dAn,"         KS test P-value = ",signif(ksp,3)) )+ 
    ggtitle(paste0("Time to Drop-off Distribution\nOutcome = ",dBn)) + 
#    annotation_custom(tableGrob(tb), xmin=2500, xmax=2501, ymin=75, ymax=76) +
    theme(plot.title=element_text(hjust=0.5)) 
#    geom_text(data = tb,aes(label =type))
  g_table <- tableGrob(tb)
  cbg <- grid.arrange(gg_graph, g_table,
               nrow=2,
               as.table=TRUE,
               heights=c(5,2))
  ggsave(file.path(opt$outputfolder, paste0(label, ".dropoff.png")), cbg,
         width=7.25, height=4.35)
}

RunCoxPH <- function(outdf, A, B, opt) {
  ###### Run CoxPH model for A to B ######
  label <- paste0(A, "_", B)
  cat(paste0("Processing ", label, " .... \n"))
  default_formula = .get_formula(outdf)
  
  # Test if user wants to control for pair variable
  if (opt$pair){
    formula <- paste0(default_formula, " + frailty(pair)")
    coetar <- c(1,2,4,6)
  }else{
    formula <- default_formula
    coetar <- c(1,3:5)
  }
  
  out <- tryCatch({
    model <- coxph(as.formula(formula), data=outdf, ties='breslow')
    testHR <- cox.zph(model)
    HRp <- testHR$table["GLOBAL", "p"]
    sig <- summary(model)$coef[1,coetar]
    expcoe <- exp(sig[1])
    testp <- sig[4]
    if (abs(expcoe) >= opt$exp_coef && testp <= opt$sigcut){
      dAn <- opt$metafile[opt$metafile[,1] == A, 2]
      dBn <- opt$metafile[opt$metafile[,1] == B, 2]
      plotStats(outdf=outdf,dAn=dAn, dBn=dBn, label=label, opt=opt) 
      cat("  Save Cox survival plot:", paste0(label, ".cox.survival.png"), "\n")
      subobj <- survfit(with(outdf, Surv(time, event)) ~ type, outdf, 
                        conf.type="log-log")
      
      autoplot(subobj, censor.shape = '*') + theme_light() +
        ylab(paste0(dBn, " Survival Prob")) +
        xlab(paste0("Time (Days)\n Exposure = ", dAn)) +
        # guides(fill=guide_legend(title="Exposure")) +
        ggtitle(paste0(dAn, " vs. ", dBn, " Cox Survival Curve")) +
        theme(plot.title=element_text(hjust=0.5))
      ggsave(file.path(opt$outputfolder, paste0(label, ".cox.survival.png")), 
             width=6.5, height=6.15)
      
      cat("  Save Cox residual plot:", paste0(label, ".cox.residual.png"), "\n")
      coxsnellres <- outdf$event-resid(model, type="martingale")
      fitres <- survfit(coxph(Surv(coxsnellres, outdf$event) ~ 1, 
                              ties='breslow'), type='aalen')
      
      ggplot(data.frame(x=fitres$time, y=-log(fitres$surv)), aes(x=x, y=y)) + 
        geom_line() + theme_light() + 
        ylab(paste0(dBn, " Estimated Cumulative Hazard Rates")) +
        xlab(paste0("Cox-Snell Residuals\n Exposure = ", dAn)) + 
        geom_abline(intercept=0, slope=1, colour="red", linetype=2) + 
        ggtitle(paste0(dAn, " vs. ", dBn, " Cox-Snell Residual Plot")) + 
        theme(plot.title=element_text(hjust=0.5))
      ggsave(file.path(opt$outputfolder, paste0(label, ".cox.residual.png")), 
             width=6.5, height=6.15)
    }
    

    res <- data.frame(from=A, to=B, sig[1],expcoe, t(sig[2:4]),nrow(outdf), HRp, row.names = NULL)
    colnames(res) <- c("from", "to", "coef", "exp_coef", 
                       "se_coef", "z", "Pr", "N", "HRtest-p")
    res
  },
  error = function(e) {
    return(NULL)
  })
  return(out)
}

RunRF <- function(outdf, A, B, opt) {
  ###### Run RF model for A to B ######
  label <- paste0(A, "_", B)

  formula <- .get_formula(outdf)
  outdf <- as.data.frame(unclass(outdf))
  outdf$type <- factor(outdf$type)
  
  dAn <- opt$metafile[opt$metafile[,1] == A, 2]
  dBn <- opt$metafile[opt$metafile[,1] == B, 2]
  
  out <- tryCatch({
    model <- rfsrc(as.formula(formula), data=outdf, 
                   ntree=opt$ntree, tree.err=TRUE)
    t <- c(model$time.interest[1], diff(model$time.interest))
    suvtab <- data.frame(model$survival)
    
    cat("  Save RF survival plot:", paste0(label, ".rf.survival.png"), "\n")
    dplt <- data.frame(timePoint=model$time.interest,
                       sapply(split(suvtab, outdf$type), colMeans))
    colnames(dplt)[2:3] <- c("NonExposed", "Exposed")
    tplt <- melt(dplt, id.vars=c("timePoint"))
    ggplot(tplt, aes(x=timePoint, y=value, col=variable)) +
      geom_line() + theme_light() +
      ylab(paste0(dBn, " Survival Prob")) +
      xlab(paste0("Time (Days)\n Exposure = ", dAn)) +
      guides(colour=guide_legend(title="Group")) +
      ggtitle(paste0(dAn, " vs. ", dBn, " RF Survival Curve"))
    ggsave(file.path(opt$outputfolder, paste0(label, ".rf.survival.png")),
           width=6.95, height=5.05)
    
    area <- rowSums(suvtab*t)
    tresult <- t.test(area ~ outdf$type, alternative='g')
    data.frame(from=A, to=B, 
               Pr_wil=wilcox.test(area ~ outdf$type, alternative='g')$p.value,
               Pr_t=tresult$p.value, 
               NonExposedMean=tresult$estimate[1],
               exposedMean=tresult$estimate[2], 
               row.names=NULL)
  },
  error = function(e) {
    return(NULL)
  })
  return(out)
}

######### FUNCTION END ##############


#### Required packages ######
.cran_packages <- c("KMsurv","survival", "randomForestSRC", "data.table", 
                    "doParallel", "foreach", "optparse", "ggplot2", 
                    "ggfortify","gridExtra","intcox")
sapply(c(.cran_packages), installLibrary)
if (!"OIsurv" %in% rownames(installed.packages())){
  install.packages('https://cran.r-project.org/src/contrib/Archive/OIsurv/OIsurv_0.2.tar.gz', repos = NULL, type="source")
  
}
if (! "intcox" %in% rownames(installed.packages())){
  install.packages("https://cran.r-project.org/src/contrib/Archive/intcox/intcox_0.9.3.tar.gz",repos=NULL, type="source")
}

#### Command line options #####
option_list = list(
  make_option(c("--inputfolder"), type="character", default='.',
              help="Input folder name. \n\t\tShould be a folder that contains all disease trajectories. Each file is named as NameA_NameB.csv. Default is current directory", metavar="character"),
  make_option(c("-m", "--metafile"), type="character", default=NULL,
              help="Meta file for each disease name [required]. \n\t\tFirst column is the disease ID used in the input file column Disease, second column is the disease name. All the other columns should be additional attribute of the disease, tab separated. Header is needed", metavar="character"),
  make_option(c("--outputfolder"), type="character", default=NULL,
              help="Survival analysis / random forest PNG output folder, default is the <inputfolder>", metavar="character"),
  make_option(c("--method"), type="character", default="CoxPH",
              help="Survival Analysis method, choises are CoxPH, RF. The default is CoxPH", metavar="character"),
  make_option(c("-o", "--outfile"), type="character", default=NULL,
              help="Adjacency matrix output file name, default is <method>.edges.csv", metavar="character"),
  make_option(c("--duration"), type="integer", default=5,
              help="study duration (years), default is 5 years", metavar="integer"),
  make_option(c("-s", "--sigcut"), type="double", default=0.05,
              help="Significant P-value cutoff. Only p-values that are less than this cutoff will be considered as significant P-values, and kept in the adjacency matrix. Default is 0.05", metavar="double"),
  make_option(c("-e", "--exp_coef"), type="double", default=1.0,
              help="Exponentiated coefficients cutoff used in plotting survival curve. Default is 1.0", metavar="double"),
  make_option(c("-a", "--adjust"), type="character", default="BH",
              help = "Method of adjusting p-values, choices are holm, hochberg, hommel, bonferroni, BY, and fdr. The default is BH", metavar="character"),
  make_option(c("-n", "--ntree"), type="integer", default=1000, 
              help="Number of trees to run random forest, default is 1000", metavar="integer"),
  make_option(c("--pair"), type="logical", action = "store_true",default=FALSE, 
              help="A flag to control for matched subject in CoxPH. Default is false"),
  make_option(c("-p", "--processors"), type="integer", default=8,
              help="Number of cores/CPUs to use, default is 8",metavar="integer")
)

opt_parser = OptionParser(usage = "usage: Rscript %prog -i <inputFileName> -m <metaFileName> [options]",
                          option_list=option_list,
                          description = "\n\t\t >> This is the step 2 of EMR package. << \n\tThis R script will perform either Cox-PH regression or random forest survival analysis on the valid files in result folder from step 1. \n\t\tOutputs:\n\t\t > CoxPH: an adjacency matrix with all the significant hazard ratios, hazard ratio test p-values, and survival curve and residual graphs.\n\t\t > RF: Wilcoxon and T test p-values, mean time to event in exposed and non-exposed group, and survival curve graphs. \n\t\t*Multiple test correction will be applied to the p-values in both methods.")
opt = parse_args(opt_parser)

###### Check arguments ######

if (any(sapply(c(.cran_packages), installLibrary) == F)) {
  print_help(opt_parser)
  stop("At least one of required R packages is not available. Please install it manually.\n", call.=FALSE)
}

### check validity of meta file format ###
if (is.null(opt$metafile)){
  print_help(opt_parser)
  stop("Please supply meta file name as -m/--metafile <filename>.\n", call.=FALSE)
}

opt$metafile <- read.csv(opt$metafile, sep="\t", stringsAsFactors=FALSE)

if (ncol(opt$metafile) < 2){
  stop("Your metafile has less than 2 columns. Please correct!\n", call.=FALSE)
}


if (!opt$method %in% c("CoxPH","RF")){
  print_help(opt_parser)
  stop("Please enter the correct method for analysis --method <method>. Choices are CoxPH, RF. Case sensitive.\n", call.=FALSE)
}

if (!(opt$adjust %in% p.adjust.methods)){
  print_help(opt_parser)
  stop("Please supply an existing method of adjusting P-value from holm, hochberg, hommel, bonferroni, BY, BH and fdr.\n", call.=FALSE)
}

if (!is.integer(opt$processors)){
  print_help(opt_parser)
  stop("Number of processors is not integer, please correct!\n", call.=FALSE)
}

if (is.null(opt$outputfolder)){
  opt$outputfolder <- opt$inputfolder
}

if (!file.exists(opt$outputfolder)) {
  dir.create(opt$outputfolder)
  cat("Creating output folder:", file.path(getwd(), opt$outputfolder), "\n")
}

if (is.null(opt$outfile)){
  opt$outfile <- paste0(opt$method, ".edge.csv")
}

##### Read in Input ########
# dir <- "/Users/ruichenrong/Projects/Loyola/Huaiying_codes/test_result"
# test_filenames <- c('A_B.csv', '1_2.csv', 'AAA_12dgf.csv', 'A.vs.B.csv', 'a1_Run.r', 'A__B.csv', 'AA_BB.csv.csv', 'AA_BB_CC.csv.csv', '__abc_abc.csv')
valid_files <- file_parser(list.files(opt$inputfolder))
cat(">> Totally", nrow(valid_files), "valid csv files in input folder: \n", 
    opt$inputfolder, "\n")

registerDoParallel(cores=opt$processors)
cat("Using ", getDoParWorkers(), " processors...\n", sep="")

#final <- foreach(i=1:nrow(valid_files), .errorhandling="remove") %do% {
final <- foreach(i=1:nrow(valid_files), .errorhandling="remove") %dopar% {
  # for (i in 1:nrow(valid_files)) {
  dt <- fread(file.path(opt$inputfolder, valid_files[i,1]), 
              sep=",", stringsAsFactors=FALSE)
  #  dt[,V3 := as.Date(V3)]
  get(paste0("Run", opt$method))(dt, valid_files[i,2], valid_files[i,3], opt)
  
}
final <- do.call(rbind.data.frame, final)

### multiple test correction ####
pvalue_col <- grep("^Pr", colnames(final), value=TRUE)
adjust_p <- sapply(pvalue_col, 
                   function(i) p.adjust(final[,i], method=opt$adjust))
colnames(adjust_p) <- paste(pvalue_col, opt$adjust, sep=".")
final <- cbind(final, adjust_p)
if (opt$method == "CoxPH") {
  final <- final[final[, "Pr"] <= opt$sigcut & final$exp_coef >= opt$exp_coef,]
}

### add metafile attribute to output ###
frmatr <- paste("from", colnames(opt$metafile)[-1], sep="_")
toatr <- paste("to", colnames(opt$metafile)[-1], sep="_")
atrmeta <- cbind(opt$metafile[match(final$from, opt$metafile[,1]), -1], 
                 opt$metafile[match(final$to, opt$metafile[,1]), -1])
colnames(atrmeta) <- c(frmatr, toatr)

write.csv(cbind(final, atrmeta), file.path(opt$outputfolder, opt$outfile), 
          row.names=FALSE, quote=FALSE)

cat("\nJob Completed!\n")

proc.time()