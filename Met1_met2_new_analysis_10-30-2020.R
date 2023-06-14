
# This analysis is for my Met1 and Met2 samples, non-targeted, ran on the Perkin Elmer
# GC-MS in Oct. 2020. As a note, these samples were run with a 'split' (jackie knows more)
# and these samples were derivatized with Methoxyamine in Pyridine and MSTFA + TMS. 

## Additionally, after Jackie looked into the analysis a bit more, we changed some parameters 
## from the first R run. This R code below should be less stringent as it uses a Bonferroni correction 
## versus using an FDR correction in the previous analysis. We are hoping this change will 
## lead to more spectra and therefore, more compounds identified. 

### This analysis was originally completed on Oct. 30, 2020. 


library(xcms)
library(csu.pmf.tools)
library(RAMClustR)

###############################################
## Section A start - XCMS processing
###############################################
## we first choose our project directory, then define our experimental parameters
wd  <-  choose.dir()  # wd <- getwd()
setwd(wd)
dir.create('datasets')
ExpDes  <-  pmfDefineExperiment() 
save(ExpDes, file = "datasets/ExpDes.Rdata")
# load("ExpDes.Rdata")
# load("datasets/ExpDes.Rdata")


# run XCMS - GCMS
#default here is whatever you saved in your experimental design file. 
# this step takes a lot of time

xset  <-  pmfxcms(
  cores=4, 
  minpw=1.5,
  maxpw=15, 
  ExpDes = ExpDes, 
  filetype = "cdf",
  bwpre=3, 
  bwpost=1, 
  minTIC=10, 
  outTIC=TRUE, 
  outPCA=TRUE, 
  snthresh=5, 
  minfrac=0.05,
  seqskip=0, 
  regroup=FALSE
)


## Section A End - XCMS

#What was done:
##peak picking, signal to noise, filters signals, fills in missing signals (because we never have a 0 signal, 
##there is always background noise)



###############################################
## Section B Start - Build RC obect, remove features
## found at levels in blank samples similar to QC.
###############################################

## build empty ramclustObj
# this step takes time
RC <- rc.get.xcms.data(
  xcmsObj = xset,
  ExpDes = ExpDes,
  mzdec = 1
)

save(RC, file = "datasets/RCobject_0_postXCMS.Rdata")

#skipping this because our samples names should be correct
fix.sample.names<-read.csv(file='fix.sample.names.csv', header=TRUE, sep=",")
RC$phenoData$sample.names<-fix.sample.names$name

#  ExpDes$design
## turn sample names into factor names

#skipping 
RC <- rc.expand.sample.names(
  ramclustObj = RC, 
  delim = "-",
  factor.names = c(
    "sn", "species", 
    "tissue", "date",
    "label", "trt", 
    "use")
)

## remove 'system peak' features
## we didn't run any blanks for this run so we skip this
RC <- rc.feature.filter.blanks(
  ramclustObj = RC,
  qc.tag = "QC",
  blank.tag = "blank",
  sn = 3
)

## replace missing values with low level noise (optionally also zero values)
##here we get 0 replacements because we didn't have any 0 in our data.
## for the met1 and met2 samples I changed the replace.zero=FALSE with TRUE. TRUE is defualt anyways and we
#need to replace any 0 with a noise level.
RC <- rc.feature.replace.na(
  ramclustObj = RC,
  replace.int = 0.1,
  replace.noise = 0.1,
  replace.zero = TRUE,
  samp.max.missing = 0.8
)
## section B end 


### The resulting number of features from this section = 415324


###############################################
## Section C: Feature normalization to correct for analytical drift
###############################################

## create 'qc' plots prenormalization for comparison
RC <- rc.qc(
  ramclustObj = RC, 
  outfile.basename = "preNorm",
  qc.tag="QC",
  remove.qc = FALSE
)

#take a look at the plots generated from this code and determine whether we need to normalize the data
#or not. based on this data skew, we need to normalize. 

## then restore 'qc' samples back to dataset for normalization
#RC <- rc.restore.qc.samples(
#  ramclustObj = RC
#)


RC <- rc.feature.normalize.qc(
  ramclustObj = RC, 
  order = order(RC$phenoData$filenames),
  batch = rep(1, nrow(RC$phenoData)),
  qc.tag = "QC",
  rsq.cut = 0.1,
  p.cut = 0.05
)

#here we get R telling us: Features were normalized by linearly regressing run order versus qc feature 
#intensities to account for instrument signal intensity drift. Only features with a regression p-value 
#less than 0.05 and an r-squared greater than 0.1 were corrected.  Of 2119 features, 389 were corrected 
#for run order effects.

##so this means that run order did not affect QC analysis. 
#

## noramlize feature values to total ion signal
## can account for both sensitivity drift and sample 
## concentration differences (i.e. urine dilution)
RC <- rc.normalize.tic(
  ramclustObj = RC
)

## rerun rc.qc to visualize QC variance post-normalization. 
RC <- rc.qc(
  ramclustObj = RC, 
  outfile.basename = "post_QCNorm"
)

#don't worry about the file error you get (eg. "Warning message:
#In dir.create("QC") : 'QC' already exists")
#take a look at the new normalized plots generated and see the variability with the QCs now. 
#for this data, the QCs should be pushed far to the left (less variability). this data isn't perfect, but
#its good enough to continue. 
#(up in the 110% of sd) 

save(RC, file = "datasets/RCobject_1_preClustering.Rdata")

## EVALUATE OUTPUT BEFORE PROCEEDING
## IF DATA QUALITY LOOKS INSUFFICIENT
## FIX BEFORE CLUSTERING AND ANNOTATING

## Section C: End

### The resulting number of features from this section = 2119


###############################################
## Section D: Cluster Features
###############################################

# cluster Features into compounds. 
#uses signal intensity and time to identify which clusters they go with
#each compound has various features and those are clustered into a spectra

RC <- rc.ramclustr(
  ramclustObj = RC, 
  st = NULL,
  sr = NULL,
  maxt = NULL,
  deepSplit = FALSE,
  blocksize = 2000,
  mult = 5,
  hmax = 0.9,
  collapse = TRUE,
  minModuleSize = 2,
  linkage = "average",
  cor.method = "pearson",
  rt.only.low.n = TRUE,
  fftempdir = NULL
)

##results we see from this code: 
#"RAMClust has condensed  features 2119 into 200 spectra 
#collapsing feature into spectral signal intensities "

save(RC, file = "datasets/RCobject_2_postClustering.Rdata")

RC <- rc.qc(
  ramclustObj = RC, 
  outfile.basename = "post_Clustering",
  qc.tag = 'QC',
  remove.qc = TRUE
)

#dont worry about the error regarding dir.create
#take a look at the QC post clustering .pdf to see variability and histograms
#these looked pretty good for this dataset

# load("datasets/RCobject_2_postClustering.Rdata")

# Export spectra
rc.export.msp.rc(ramclustObj = RC, one.file = TRUE)

#after exporting the spectra we can use this in RAMSearch for cluster ID

##now we are ready for annotating
## Section D: End




### The resulting number of features from this section = 2119 features condensed into200 spectra






















#####################

# BREAK in SCRIPT:
# ____ run RAMSearch manually

## import MSFinder results
RC <- impRamSearch(ramclustObj = RC)

RC <- annotate(
  ramclustObj = RC,
  standardize.names = FALSE,
  min.msms.score = 6.5,
  database.priority = NULL,
  citation.score = TRUE,
  find.inchikey = FALSE,
  taxonomy.inchi = NULL,
  reset = TRUE)


## Error with this command <- Error in `[.data.frame`(tmp, , 1:5) : undefined columns selected
## use pubchem rest API to retrieve a great deal of data on each annoation
# RC <- rc.cmpd.get.pubchem(ramclustObj = RC)

# The below didn't give an error but likely didn't work based upon the fact that the pubchem command didn't work
## use classyFire to retreive compounds classification for each annotation
# RC <- getClassyFire(ramclustObj = RC)

RC <- degolm(ramclustObj = RC)

save(RC, file = "datasets/RCobject_3_postAnnotation.Rdata")

## Section E: End



###############################################
## Section F: Statistical Analysis
###############################################

# Examine factors to guide analysis
des <- getData(ramclustObj = RC)[[1]]
head(des, n = 3)



#  PCA
# factors were misnamed: sn=sn, species = label, tissue = species, date = tissue, label=date, trt = trt, use = use
RC <- pmfpca(ramclustObj = RC,  
             which.factors = c("date","trt","label"), 
             npc = "auto", num.factors = NULL)


#  ANOVA
RC <- pmfanova(ramclustObj = RC, anova.call = "label", subset=c("date", "Leaf", "use", "Change over time"))
RC <- pmfanova(ramclustObj = RC, anova.call = "trt*date", subset=c( "use", "Metabolomics at Harvest", "label", "March8"))
RC <- pmfanova(ramclustObj = RC, anova.call = "trt*date")

save(RC, file = "datasets/RCobject_4_postStats.Rdata")


## Section F: End



###############################################
## Section G: Reporting and file export
###############################################

## Export SpecAbund Dataset - contains quantitative signal intensity values
exportDataset(ramclustObj=RC, which.data="SpecAbund")

## Export Annotation Summary file (spectra directory)
annotation.summary(ramclustObj=RC)

## Export methods summary for all R based post-XCMS steps
write.methods(ramclustObj = RC)

## print sample names to console, return summarized values in RC object
RC <- reportSampleNames(ramclustObj = RC)

## zip files for sharing
make.zip.files(do.raw = FALSE)

## write a text file with descriptions of output files.  
write.file.summary()

## Section G: End


