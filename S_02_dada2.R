# This script was modified to use for the analysis of bacterial amplicon data
# See: https://benjjneb.github.io/dada2/tutorial.html
# It is not optimized for:
#   amplicons which are sequenced into the rev primer at the 5' end of the fwd read, and vice versa (e.g. ITS)
#   faster computing times at slightly lower output quality (e.g. big data)

# Christiane Hassenruecks own optimization was implemented, snakemake code can be found here
#https://git.io-warnemuende.de/bio_inf/workflow_templates/src/branch/master/Amplicon_dada2_MiSeq
 # differences at merging step to not loose organisms (e.g. large SOX bacteria) with large 16S gene which cannot be merged
 # need reference for mapping for this!

## This script was used for analysis of the bacterial sequences of the Novogene libraries 14,15,16, 46, 56,57,58
# Potter Cove temperature SIP experiment original and deeper re-sequencing
# same samples are in multiple libraries for deeper sequencing 
# samples are pooled after seq table creation

## run on Server

# load module on server:  

# open in separate screen, especially for large jobs:
#screen -S R
# activate environment needed for further mapping below before opening R
#module load bbmap/38.86
#conda activate metawrap-env
#module load R/4.4.1 # version number might be different, use autocomplete

## make sure that all needed directories are in working directory:
# QualityProfiles Logfiles ErrorProfiles Output

# open R while in the correct working directory: 


# open R by typing R

## load packages ####
require(dada2)
require(ShortRead)
require(gridExtra)
require(tidyverse)

packageVersion("dada2")
# 1.32.0

# save and load workspace
## path to working directory
wd <- "WORKINGDIRECTORY/"
setwd(wd)

## workspace file
tf <- "TEMPORARYFILE.Rdata"
# save.image(tf)
# load(tf) # to load workspace again

## path to scripts to load
sloc <- "SCRIPTLOCATION/Scripts_to_source/"

## path to database for mapping
dblm <- "DATABASE_LOCATION_MAP/"

## path to database for classification
dblc <- "DATABASE_LOCATION_CLAS/"

## project name to name files
prj <- "PROJECT"


# specify path to input fastq files
## this file has the names of the flowcells with the lane afterwards in each line, no header
 # created automatically in script dada2_seqprep.bash, can also be easily made by hand
seq_lanes <- read.table(paste0(wd, "flowcell_IDs.txt"),
                        col.names = "flowcell", sep = "\t")
rownames(seq_lanes) <- seq_lanes$flowcell

## here replace for Arc if analysing archaea libraries
TARGET <- "Bac"

extract_sample_paths <- function(flowcell_ID) {
path <- paste0(wd, "seq_by_flowcell/", flowcell_ID, "/Clipped_", TARGET)
fnFR_R1 <- sort(list.files(path, pattern="clip_fr_R1.fastq", full.names = TRUE))
fnFR_R2 <- sort(list.files(path, pattern="clip_fr_R2.fastq", full.names = TRUE))
fnRF_R1 <- sort(list.files(path, pattern="clip_rf_R1.fastq", full.names = TRUE))
fnRF_R2 <- sort(list.files(path, pattern="clip_rf_R2.fastq", full.names = TRUE))
# Extract sample names
sample.names <- sapply(strsplit(basename(fnFR_R1), "_"), `[`, 1)
flowcell_ID <- setNames(list(fnFR_R1, fnFR_R2, fnRF_R1, fnRF_R2, sample.names), 
                        c("fnFR_R1", "fnFR_R2", "fnRF_R1", "fnRF_R2", "sample.names"))
return(flowcell_ID)
}

### create list with each seq_lane as one element
sample.paths <- apply(seq_lanes, 1, extract_sample_paths)

## quality check and trimming ####
 # skip this, use the ones from previous analysis. just 6 samples or so added, should not influence profiles
 # proceed to trimming
 # copy quality profiles and parameter optimization from directory new_analysis_2024_2
source(paste0(sloc, "dada2_quality_check.R"))
mapply(function(x, i) {
  quality_check(
    c(x[['fnFR_R1']], x[['fnRF_R1']]),              # need '' inside [[]] to access list elements by name
    c(x[['fnFR_R2']], x[['fnRF_R2']]),
    file_base = paste0("QualityProfiles/", i,"_QualityProfile_separate") # need to create folder before!
  )},
  x = sample.paths, i = names(sample.paths),        # need mapply instead of lapply to access the names of the list
  SIMPLIFY = F
)


### check pdfs:
# Considerations for trimming:
# expected max length: 252bp (?)
# min overlap: 30bp  --> in total R1 + R2 seq should have 285-290 bp total length
# reads should be truncated so that rev primer is not included at end of fwd reads
# It is recommended to trim to just enough for the required length for sufficient overlap
# Caution: don't remove too much

## Run parameter optimization
# Define ranges for maxEE (normally same for all)
range_maxEE <- matrix(
  c(1, 1,
    2, 1,
    2, 2,
    3, 2,
    3, 3),
  nrow = 5,
  ncol = 2,
  byrow = T
)

# add to list according parameters for testing truncation length and maxEE
# check quality profile pdfs for this:
### seq lane 1/6 H7LM7DRX2_L2
  ## forward (R1)
  # looks very good
  # potential cut-off: 180-190 bp
  ## reverse (R2)
  # also good, in rf some bad peaks at 175-180 bp
  # potential cut-off: 175 bp
sample.paths[['H7LM7DRX2_L2']][['range_truncLen']] <- matrix(
    c(130, 160, # first number for R1 read (= forward), second number for R2 read (=reverse)
      115, 175,
      120, 170,
      100, 190),
  nrow = 4, ncol = 2, byrow = T)
sample.paths[['H7LM7DRX2_L2']][['range_maxEE']] <- range_maxEE

### seq lane 2/6 H7W7NDRX2_L2
## forward (R1)
# potential cut-off: 100 or 125
## reverse (R2)
# cut-off: 145 or 175
sample.paths[['H7W7NDRX2_L2']][['range_truncLen']] <- matrix(
  c(145, 145, # first number for R1 read (= forward), second number for R2 read (=reverse)
    115, 175,
    110, 180,
    105, 185,
    100, 190),
  nrow = 5, ncol = 2, byrow = T)
sample.paths[['H7W7NDRX2_L2']][['range_maxEE']] <- range_maxEE

### seq lane 3/6 H7WKNDRX2_L2
## forward (R1)
# potential cut-off: 105 or 145
## reverse (R2)
# cut-off: 145 or 175-180
sample.paths[['H7WKNDRX2_L2']][['range_truncLen']] <- matrix(
  c(145, 145, # first number for R1 read (= forward), second number for R2 read (=reverse)
    115, 175,
    110, 180,
    105, 185),
  nrow = 4, ncol = 2, byrow = T)
sample.paths[['H7WKNDRX2_L2']][['range_maxEE']] <- range_maxEE

### seq lane 4/6 HC2JWDRX2_L1
## forward (R1)
# potential cut-off: 105 or 120
## reverse (R2)
# cut-off: 145 or 175. fr really good bad rf bad peaks
sample.paths[['HC2JWDRX2_L1']][['range_truncLen']] <- matrix(
  c(145, 145, # first number for R1 read (= forward), second number for R2 read (=reverse)
    115, 175,
    110, 180,
    120, 170),
  nrow = 4, ncol = 2, byrow = T)
sample.paths[['HC2JWDRX2_L1']][['range_maxEE']] <- range_maxEE

### seq lane 5/6 HJW2YDRXX_L1
## forward (R1)
# look really bad
# potential cut-off: as short as possible
## reverse (R2)
# fr drops past 170, rf bad peaks already at 145, 180
# cut-off: 145 or 180
sample.paths[['HJW2YDRXX_L1']][['range_truncLen']] <- matrix(
  c(145, 145, # first number for R1 read (= forward), second number for R2 read (=reverse)
    120, 170,
    110, 180,
    100, 190,
    95, 195),
  nrow = 5, ncol = 2, byrow = T)
sample.paths[['HJW2YDRXX_L1']][['range_maxEE']] <- range_maxEE

### seq lane 6/6 HJW2YDRXX_L2
## forward (R1)
# look really bad
# potential cut-off: as short as possible
## reverse (R2)
# fr drops past 170, rf bad peaks already at 145, 180
# cut-off: 145 or 180
sample.paths[['HJW2YDRXX_L2']][['range_truncLen']] <- matrix(
  c(145, 145, # first number for R1 read (= forward), second number for R2 read (=reverse)
    120, 170,
    110, 180,
    100, 190,
    95, 195),
  nrow = 5, ncol = 2, byrow = T)
sample.paths[['HJW2YDRXX_L2']][['range_maxEE']] <- range_maxEE


# Prepare directories to place filtered files in Filtered_TARGET/ subdirectory in folder seq_by_flowcell individual flowcell folders
sample.paths <- mapply(function(x, i) {
  x[['filtFR_R1']] <- file.path(paste0(wd, "/seq_by_flowcell/", i, "/Filtered_", TARGET), 
                                paste0(x[['sample.names']], "_FR_R1_filt.fastq"))
  x[['filtFR_R2']] <- file.path(paste0(wd, "/seq_by_flowcell/", i, "/Filtered_", TARGET), 
                                paste0(x[['sample.names']], "_FR_R2_filt.fastq"))
  x[['filtRF_R1']] <- file.path(paste0(wd, "/seq_by_flowcell/", i, "/Filtered_", TARGET), 
                                paste0(x[['sample.names']], "_RF_R1_filt.fastq"))
  x[['filtRF_R2']] <- file.path(paste0(wd, "/seq_by_flowcell/", i, "/Filtered_", TARGET), 
                                paste0(x[['sample.names']], "_RF_R2_filt.fastq"))
  names(x[['filtFR_R1']]) <- x[['sample.names']]
  names(x[['filtFR_R2']]) <- x[['sample.names']]
  names(x[['filtRF_R1']]) <- x[['sample.names']]
  names(x[['filtRF_R2']]) <- x[['sample.names']]
  x
},
x = sample.paths, i = names(sample.paths),
SIMPLIFY = F
)

# Run parameter optimization
# This is quite time consuming and should only be attempted on a server with as many cores as you have samples (or at least 20)
source(paste0(sloc, "dada2_screen_settings.R"))
threads <- 150
sample.paths <- mapply(function(x, i) {
  x[['screen_filt_settings_fr']] <- 
    screen_settings(sample.names = x[['sample.names']], fnFs = x[['fnFR_R1']], fnRs = x[['fnFR_R2']], 
                    range_maxEE = x[['range_maxEE']], range_truncLen = x[['range_truncLen']], cpus = threads)
  x[['screen_filt_settings_rf']] <- 
    screen_settings(sample.names = x[['sample.names']], fnFs = x[['fnRF_R1']], fnRs = x[['fnRF_R2']], 
                    range_maxEE = x[['range_maxEE']], range_truncLen = x[['range_truncLen']], cpus = threads)
  x
},
  x = sample.paths, i = names(sample.paths),
  SIMPLIFY = F
)

save.image(tf)



## make plots
mapply(function(x, i) {
pdf(paste0("QualityProfiles/Parameter_screening_", i, ".pdf"), width = 7, height = 7)
plot(
  x[['screen_filt_settings_fr']][, "prop.total"],
  x[['screen_filt_settings_fr']][, "q90"] - x[['screen_filt_settings_fr']][, "q10"],
  col = rep(rainbow(nrow(x[['range_maxEE']])), nrow(x[['range_truncLen']])),
  pch = 16
)
text(
  x[['screen_filt_settings_fr']][, "prop.total"],
  x[['screen_filt_settings_fr']][, "q90"] - x[['screen_filt_settings_fr']][, "q10"],
  pos = 2,
  col = adjustcolor("black", alpha.f = 0.5)
)
plot(
  x[['screen_filt_settings_rf']][, "prop.total"],
  x[['screen_filt_settings_rf']][, "q90"] - x[['screen_filt_settings_rf']][, "q10"],
  col = rep(rainbow(nrow(x[['range_maxEE']])), nrow(x[['range_truncLen']])),
  pch = 16
)
text(
  x[['screen_filt_settings_rf']][, "prop.total"],
  x[['screen_filt_settings_rf']][, "q90"] - x[['screen_filt_settings_rf']][, "q10"],
  pos = 2,
  col = adjustcolor("black", alpha.f = 0.5)
)
dev.off()
}, 
x = sample.paths, i = names(sample.paths),
SIMPLIFY = F
)


# This is just a gut feeling, but I would optimize for the following criteria:
#   small difference between 10 and 90 percentile of retained reads
#   high total proportion of retained reads
#   most stringent maxEE that does not result in severe loss of reads

# look at numbers in table:
sample.paths[['flowcell']][['screen_filt_settings_fr']]
sample.paths[['flowcell']][['screen_filt_settings_rf']]

## Run trimming with optimal parameters

# define optimal trunc length, these parameters will be needed again further down

# seq lane 1/6 H7LM7DRX2_L2 ##
 # truncLen = c(100,190), maxEE = c(2,2);
 # index 18, fr: prop.total 0.9775; q10 0.9723; q90 0.9822; rf: prop.total 0.9714; q10 0.9683; q90 0.9750;
sample.paths[['H7LM7DRX2_L2']][['truncLen_R1']] <- 100
sample.paths[['H7LM7DRX2_L2']][['truncLen_R2']] <- 190
sample.paths[['H7LM7DRX2_L2']][['error_R1']] <- 2
sample.paths[['H7LM7DRX2_L2']][['error_R2']] <- 2

# seq lane 2/6 H7W7NDRX2_L2 ##
# truncLen = c(110,180), maxEE = c(2,2);
# index 8, fr: prop.total 0.9423; q10 0.9326; q90 0.9481; rf: prop.total 0.9295; q10 0.9205; q90 0.9351;
sample.paths[['H7W7NDRX2_L2']][['truncLen_R1']] <- 110
sample.paths[['H7W7NDRX2_L2']][['truncLen_R2']] <- 180
sample.paths[['H7W7NDRX2_L2']][['error_R1']] <- 2
sample.paths[['H7W7NDRX2_L2']][['error_R2']] <- 2

# seq lane 3/6 H7WKNDRX2_L2 ##
# truncLen = c(105,185), maxEE = c(2,2);
# index 18, fr: prop.total 0.9698; q10 0.9641; q90 0.9758; rf: prop.total 0.9623; q10 0.9580; q90 0.9666;
sample.paths[['H7WKNDRX2_L2']][['truncLen_R1']] <- 105
sample.paths[['H7WKNDRX2_L2']][['truncLen_R2']] <- 185
sample.paths[['H7WKNDRX2_L2']][['error_R1']] <- 2
sample.paths[['H7WKNDRX2_L2']][['error_R2']] <- 2

# seq lane 4/6 HC2JWDRX2_L1 ##
# truncLen = c(110,180), maxEE = c(2,2);
# index 13, fr: prop.total 0.9361; q10 0.9326; q90 0.9420; rf: prop.total 0.9277; q10 0.9254; q90 0.9350;
sample.paths[['HC2JWDRX2_L1']][['truncLen_R1']] <- 110
sample.paths[['HC2JWDRX2_L1']][['truncLen_R2']] <- 180
sample.paths[['HC2JWDRX2_L1']][['error_R1']] <- 2
sample.paths[['HC2JWDRX2_L1']][['error_R2']] <- 2

# seq lane 5/6 HJW2YDRXX_L1 ##
# truncLen = c(100,190), maxEE = c(2,2);
# index 18, fr: prop.total 0.9254; q10 0.9065; q90 0.9393; rf: prop.total 0.9134; q10 0.8918; q90 0.9282;
sample.paths[['HJW2YDRXX_L1']][['truncLen_R1']] <- 100
sample.paths[['HJW2YDRXX_L1']][['truncLen_R2']] <- 190
sample.paths[['HJW2YDRXX_L1']][['error_R1']] <- 2
sample.paths[['HJW2YDRXX_L1']][['error_R2']] <- 2

# seq lane 6/6 HJW2YDRXX_L2 ##
# truncLen = c(110,180), maxEE = c(2,2);
# index 13, fr: prop.total 0.9121; q10 0.8906; q90 0.9265; rf: prop.total 0.8998; q10 0.8794; q90 0.9152;
sample.paths[['HJW2YDRXX_L2']][['truncLen_R1']] <- 110
sample.paths[['HJW2YDRXX_L2']][['truncLen_R2']] <- 180
sample.paths[['HJW2YDRXX_L2']][['error_R1']] <- 2
sample.paths[['HJW2YDRXX_L2']][['error_R2']] <- 2


# run trimming
threads <- 150

sample.paths <- mapply(function(x, i) {
  x[['filt_FR.out']] <- filterAndTrim(
  fwd = x[['fnFR_R1']], 
  filt = x[['filtFR_R1']], 
  rev = x[['fnFR_R2']], 
  filt.rev = x[['filtFR_R2']],
  truncLen = c(x[['truncLen_R1']], x[['truncLen_R2']]),
  maxN = 0,
  minQ = 2,
  maxEE = c(x[['error_R1']], x[['error_R2']]), 
  truncQ = 0, 
  rm.phix = TRUE,
  compress = F,
  multithread = threads
  )

  x[['filt_RF.out']] <- filterAndTrim(
  fwd = x[['fnRF_R1']], 
  filt = x[['filtRF_R1']], 
  rev = x[['fnRF_R2']], 
  filt.rev = x[['filtRF_R2']],
  truncLen = c(x[['truncLen_R1']], x[['truncLen_R2']]),
  maxN = 0,
  minQ = 2,
  maxEE = c(x[['error_R1']], x[['error_R2']]), 
  truncQ = 0, 
  rm.phix = TRUE,
  compress = F,
  multithread = threads
  )
  x
}, 
x = sample.paths, i = names(sample.paths),
SIMPLIFY = F
)

## Repeat quality check after trimming
source(paste0(sloc, "dada2_quality_check.R"))
mapply(function(x, i) {
quality_check(
  c(x[['filtFR_R1']], x[['filtRF_R1']]),
  c(x[['filtFR_R2']], x[['filtRF_R2']]),
  file_base = paste0("QualityProfiles/", i, "_QualityProfileFiltered_", x[['truncLen_R1']], "_", x[['truncLen_R2']])
)
}, 
x = sample.paths, i = names(sample.paths),
SIMPLIFY = F
)

save.image(tf)



### Learn error rates ####
# It is generally not necessary to increase the number of nbases used for the error estimation
# It is possible that with 10 rounds (MAX_CONSIST), the algorithm for learning the errors won't converge
# Increasing MAX_CONSIST will lead to longer run times, and may only marginally improve error estimation
# I would not recommend setting MAX_CONSIST higher than 15

## correct error estimates because of binned quality scores?
# https://github.com/benjjneb/dada2/issues/791
# https://github.com/benjjneb/dada2/issues/938
# The dada2 team is working on making some good recommendations for binned quality scores
# For now, there are 2 options:
#   1) run error learning with modified loess function (maybe more elegant)
#   Hack the loessErrfun() of dada2 package (used in both learnErrors and dada): 
#   mod.lo <- loess(rlogp ~ q, df, weights = log10(tot), span = 2)
#   2) coerce any value lower than the Q40 probability to be the Q40 value in the learnErrors() output
#   We will do this here to avoid re-running the error learning


## Use alternative Loess function, this one always gave best results from binned quality scores
source(paste0(sloc, "loessErrfun2.R"))
threads <- 150

sample.paths <- mapply(function(x, i) {
  
sink(paste0("Logfiles/", i, "_", TARGET, "_log_errorEstimation.log"))
  
x[['errFR_R1']] <- learnErrors(x[['filtFR_R1']], errorEstimationFunction = loessErrfun2, multithread = threads, randomize = TRUE, verbose = 1, MAX_CONSIST = 15)
x[['errFR_R2']] <- learnErrors(x[['filtFR_R2']], errorEstimationFunction = loessErrfun2, multithread = threads, randomize = TRUE, verbose = 1, MAX_CONSIST = 15)
x[['errRF_R1']] <- learnErrors(x[['filtRF_R1']], errorEstimationFunction = loessErrfun2, multithread = threads, randomize = TRUE, verbose = 1, MAX_CONSIST = 15)
x[['errRF_R2']] <- learnErrors(x[['filtRF_R2']], errorEstimationFunction = loessErrfun2, multithread = threads, randomize = TRUE, verbose = 1, MAX_CONSIST = 15)

sink()
x

}, 
x = sample.paths, i = names(sample.paths),
SIMPLIFY = F
)

# check logfiles in bash console:
## samples H7LM7DRX2_L2 errFR_R2 did not reach convergence
# rest of samples reached convergence after 7-12 rounds

# if convergence was not reached for all samples, error learning can be re-run only for certain samples
# change sequence lane ID and exact name (L1 or L2) manually in code below:
sink(paste0("Logfiles/H7LM7DRX2_L2_", TARGET, "_log_errorEstimation.log"), append = T)
sample.paths[['H7LM7DRX2_L2']][['errFR_R2']] <- learnErrors(sample.paths[['H7LM7DRX2_L2']][['filtFR_R2']], errorEstimationFunction = loessErrfun2, multithread = threads, randomize = TRUE, verbose = 1, MAX_CONSIST = 20)
sink()
# reached convergence after 8 rounds


# plot error profiles
## do plots manually
x <- "H7LM7DRX2_L2"
pdf(paste0("ErrorProfiles/", x, "_ErrorProfiles_separate.pdf"))
barplot(log10(dada2:::checkConvergence(sample.paths[[x]][['errFR_R1']]) + 1), main = "Convergence_fwd")
barplot(log10(dada2:::checkConvergence(sample.paths[[x]][['errFR_R2']]) + 1), main = "Convergence_rev")
barplot(log10(dada2:::checkConvergence(sample.paths[[x]][['errRF_R1']]) + 1), main = "Convergence_fwd")
barplot(log10(dada2:::checkConvergence(sample.paths[[x]][['errRF_R2']]) + 1), main = "Convergence_rev")
plotErrors(sample.paths[[x]][['errFR_R1']], nominalQ = TRUE)
plotErrors(sample.paths[[x]][['errFR_R2']], nominalQ = TRUE)
plotErrors(sample.paths[[x]][['errRF_R1']], nominalQ = TRUE)
plotErrors(sample.paths[[x]][['errRF_R2']], nominalQ = TRUE)
dev.off()

x <- "H7W7NDRX2_L2"
pdf(paste0("ErrorProfiles/", x, "_ErrorProfiles_separate.pdf"))
barplot(log10(dada2:::checkConvergence(sample.paths[[x]][['errFR_R1']]) + 1), main = "Convergence_fwd")
barplot(log10(dada2:::checkConvergence(sample.paths[[x]][['errFR_R2']]) + 1), main = "Convergence_rev")
barplot(log10(dada2:::checkConvergence(sample.paths[[x]][['errRF_R1']]) + 1), main = "Convergence_fwd")
barplot(log10(dada2:::checkConvergence(sample.paths[[x]][['errRF_R2']]) + 1), main = "Convergence_rev")
plotErrors(sample.paths[[x]][['errFR_R1']], nominalQ = TRUE)
plotErrors(sample.paths[[x]][['errFR_R2']], nominalQ = TRUE)
plotErrors(sample.paths[[x]][['errRF_R1']], nominalQ = TRUE)
plotErrors(sample.paths[[x]][['errRF_R2']], nominalQ = TRUE)
dev.off()

x <- "H7WKNDRX2_L2"
pdf(paste0("ErrorProfiles/", x, "_ErrorProfiles_separate.pdf"))
barplot(log10(dada2:::checkConvergence(sample.paths[[x]][['errFR_R1']]) + 1), main = "Convergence_fwd")
barplot(log10(dada2:::checkConvergence(sample.paths[[x]][['errFR_R2']]) + 1), main = "Convergence_rev")
barplot(log10(dada2:::checkConvergence(sample.paths[[x]][['errRF_R1']]) + 1), main = "Convergence_fwd")
barplot(log10(dada2:::checkConvergence(sample.paths[[x]][['errRF_R2']]) + 1), main = "Convergence_rev")
plotErrors(sample.paths[[x]][['errFR_R1']], nominalQ = TRUE)
plotErrors(sample.paths[[x]][['errFR_R2']], nominalQ = TRUE)
plotErrors(sample.paths[[x]][['errRF_R1']], nominalQ = TRUE)
plotErrors(sample.paths[[x]][['errRF_R2']], nominalQ = TRUE)
dev.off()

x <- "HC2JWDRX2_L1"
pdf(paste0("ErrorProfiles/", x, "_ErrorProfiles_separate.pdf"))
barplot(log10(dada2:::checkConvergence(sample.paths[[x]][['errFR_R1']]) + 1), main = "Convergence_fwd")
barplot(log10(dada2:::checkConvergence(sample.paths[[x]][['errFR_R2']]) + 1), main = "Convergence_rev")
barplot(log10(dada2:::checkConvergence(sample.paths[[x]][['errRF_R1']]) + 1), main = "Convergence_fwd")
barplot(log10(dada2:::checkConvergence(sample.paths[[x]][['errRF_R2']]) + 1), main = "Convergence_rev")
plotErrors(sample.paths[[x]][['errFR_R1']], nominalQ = TRUE)
plotErrors(sample.paths[[x]][['errFR_R2']], nominalQ = TRUE)
plotErrors(sample.paths[[x]][['errRF_R1']], nominalQ = TRUE)
plotErrors(sample.paths[[x]][['errRF_R2']], nominalQ = TRUE)
dev.off()

x <- "HJW2YDRXX_L1"
pdf(paste0("ErrorProfiles/", x, "_ErrorProfiles_separate.pdf"))
barplot(log10(dada2:::checkConvergence(sample.paths[[x]][['errFR_R1']]) + 1), main = "Convergence_fwd")
barplot(log10(dada2:::checkConvergence(sample.paths[[x]][['errFR_R2']]) + 1), main = "Convergence_rev")
barplot(log10(dada2:::checkConvergence(sample.paths[[x]][['errRF_R1']]) + 1), main = "Convergence_fwd")
barplot(log10(dada2:::checkConvergence(sample.paths[[x]][['errRF_R2']]) + 1), main = "Convergence_rev")
plotErrors(sample.paths[[x]][['errFR_R1']], nominalQ = TRUE)
plotErrors(sample.paths[[x]][['errFR_R2']], nominalQ = TRUE)
plotErrors(sample.paths[[x]][['errRF_R1']], nominalQ = TRUE)
plotErrors(sample.paths[[x]][['errRF_R2']], nominalQ = TRUE)
dev.off()

x <- "HJW2YDRXX_L2"
pdf(paste0("ErrorProfiles/", x, "_ErrorProfiles_separate.pdf"))
barplot(log10(dada2:::checkConvergence(sample.paths[[x]][['errFR_R1']]) + 1), main = "Convergence_fwd")
barplot(log10(dada2:::checkConvergence(sample.paths[[x]][['errFR_R2']]) + 1), main = "Convergence_rev")
barplot(log10(dada2:::checkConvergence(sample.paths[[x]][['errRF_R1']]) + 1), main = "Convergence_fwd")
barplot(log10(dada2:::checkConvergence(sample.paths[[x]][['errRF_R2']]) + 1), main = "Convergence_rev")
plotErrors(sample.paths[[x]][['errFR_R1']], nominalQ = TRUE)
plotErrors(sample.paths[[x]][['errFR_R2']], nominalQ = TRUE)
plotErrors(sample.paths[[x]][['errRF_R1']], nominalQ = TRUE)
plotErrors(sample.paths[[x]][['errRF_R2']], nominalQ = TRUE)
dev.off()



### Dereplicate and denoise samples ####
threads <- 150
sample.paths <- mapply(function(x, i) {
  
sink(paste0("Logfiles/", i, "_", TARGET, "_derep_denoise.log"))

x[['dadaFR_R1']] <- dada(x[['filtFR_R1']], err = x[['errFR_R1']], errorEstimationFunction = loessErrfun2, multithread = threads, pool = TRUE)
x[['dadaFR_R2']] <- dada(x[['filtFR_R2']], err = x[['errFR_R2']], errorEstimationFunction = loessErrfun2, multithread = threads, pool = TRUE)
x[['dadaRF_R1']] <- dada(x[['filtRF_R1']], err = x[['errRF_R1']], errorEstimationFunction = loessErrfun2, multithread = threads, pool = TRUE)
x[['dadaRF_R2']] <- dada(x[['filtRF_R2']], err = x[['errRF_R2']], errorEstimationFunction = loessErrfun2, multithread = threads, pool = TRUE)

sink()

x
}, 
x = sample.paths, i = names(sample.paths),
SIMPLIFY = F
)

# it is a good idea to save your workspace here
save.image(tf)


## Merge reads ####
# mismatches are 0 by default
# verbose = TRUE prints merging pair numbers to screen
source(paste0(sloc, "extract_unmerged.R"))

sample.paths <- mapply(function(x, i) {
x[['mergers_FR0']] <- mergePairs(
  x[['dadaFR_R1']],
  x[['filtFR_R1']], 
  x[['dadaFR_R2']], 
  x[['filtFR_R2']], 
  minOverlap = 10,
  verbose = TRUE,
  returnRejects = TRUE
)
x[['mergers_RF0']] <- mergePairs(
  x[['dadaRF_R1']],
  x[['filtRF_R1']], 
  x[['dadaRF_R2']], 
  x[['filtRF_R2']], 
  minOverlap = 10,
  verbose = TRUE,
  returnRejects = TRUE
)


### rescue unmerged reads ####

#### extract them and save in fasta file
x[['unmerged_FR']] <- extract_unmerged(x[['dadaFR_R1']], x[['dadaFR_R2']], x[['mergers_FR0']])
x[['unmerged_RF']] <- extract_unmerged(x[['dadaRF_R1']], x[['dadaRF_R2']], x[['mergers_RF0']])

writeFasta(x[['unmerged_FR']][[1]], file = paste0(i, "_", TARGET, "_unmerged_FR_R1.fasta"))
writeFasta(x[['unmerged_FR']][[2]], file = paste0(i, "_", TARGET, "_unmerged_FR_R2.fasta"))
writeFasta(x[['unmerged_RF']][[1]], file = paste0(i, "_", TARGET, "_unmerged_RF_R1.fasta"))
writeFasta(x[['unmerged_RF']][[2]], file = paste0(i, "_", TARGET, "_unmerged_RF_R2.fasta"))

x
}, 
x = sample.paths, i = names(sample.paths),
SIMPLIFY = F
)


#### outside of R reads will be mapped with bbmap.sh and insert size extracted
   # can access outside with system(), conda environment has to be active when opening R!!!
threads <- 150
map_ref <- paste0(dblm, "SILVA_138.1_SSURef_NR99_tax_silva_trunc.fasta")

sample.paths <- mapply(function(x, i) {
  
system(paste0(
  "bbmap.sh threads=", threads, " ref=", map_ref, " in=", i, "_", TARGET, "_unmerged_FR_R1.fasta in2=", 
  i, "_", TARGET, "_unmerged_FR_R2.fasta out=", i, "_", TARGET, "_unmerged_FR.bam"))

system(paste0(
  "bbmap.sh threads=", threads, " ref=", map_ref, " in=", i, "_", TARGET, "_unmerged_RF_R1.fasta in2=", 
  i, "_", TARGET, "_unmerged_RF_R2.fasta out=", i, "_", TARGET, "_unmerged_RF.bam"))


##### extract insert size
system(paste0(
  "samtools view -F2304 -f66 -m50 ",
  i, "_", TARGET, "_unmerged_FR.bam | cut -f1,9 > ", i, "_", TARGET, "_unmerged_is_FR.txt"))

system(paste0(
  "samtools view -F2304 -f66 -m50 ",
  i, "_", TARGET, "_unmerged_RF.bam | cut -f1,9 > ", i, "_", TARGET, "_unmerged_is_RF.txt"))


#### read insert size back into R
x[['is.fr']] <- read.table(paste0(i, "_", TARGET, "_unmerged_is_FR.txt"), h = F, sep = "\t", col.names = c("seqID", "insert"))
x[['is.rf']] <- read.table(paste0(i, "_", TARGET, "_unmerged_is_RF.txt"), h = F, sep = "\t", col.names = c("seqID", "insert"))


#### filter to insert sizes that exceed maximum length of merged sequences
x[['is_long.fr']] <- x[['is.fr']][x[['is.fr']][['insert']] > (x[['truncLen_R1']] + x[['truncLen_R2']] - 10), ] %>% 
  separate(seqID, into = c("sample_index", "row_index"), sep = "_", remove = F, convert = T)

x[['is_long.rf']] <- x[['is.rf']][x[['is.rf']][['insert']] > (x[['truncLen_R1']] + x[['truncLen_R2']] - 10), ] %>% 
  separate(seqID, into = c("sample_index", "row_index"), sep = "_", remove = F, convert = T)

x
}, 
x = sample.paths, i = names(sample.paths),
SIMPLIFY = F
)

save.image(tf)


#### retrieve and concatenate sequence
source(paste0(sloc, "retr_concat_unmerged.R"))

sample.paths <- mapply(function(x, i) {

x[['mergers_FR']] <- retr_concat_unmerged(x[['mergers_FR0']], x[['is_long.fr']], x[['unmerged_FR']])
x[['mergers_RF']] <- retr_concat_unmerged(x[['mergers_RF0']], x[['is_long.rf']], x[['unmerged_RF']])

x
}, 
x = sample.paths, i = names(sample.paths),
SIMPLIFY = F
)


# Create sequence table ####
# Create sequence table with actual sequences as column names
sample.paths <- mapply(function(x, i) {
  
x[['seqtab_FR']] <- makeSequenceTable(x[['mergers_FR']])
x[['seqtab_RF']] <- makeSequenceTable(x[['mergers_RF']])

x
}, 
x = sample.paths, i = names(sample.paths),
SIMPLIFY = F
)


mapply(function(x, i) {
# for dimension first number samples, second number sequences (=ASVs)
dim(x[['seqtab_FR']])
dim(x[['seqtab_RF']])
}, 
x = sample.paths, i = names(sample.paths),
SIMPLIFY = F
)

# $H7LM7DRX2_L2
# [1]    70 13381
# 
# $H7W7NDRX2_L2
# [1]   70 4598
# 
# $H7WKNDRX2_L2
# [1]    55 11046
# 
# $HC2JWDRX2_L1
# [1]   15 3925
# 
# $HJW2YDRXX_L1
# [1]  115 3231
# 
# $HJW2YDRXX_L2
# [1]   86 4088

# As with the tryRC option of mergeSequenceTables only those sequences which are duplicated
# will be turned, manually turn sequences of RF table
# seqtab_RF_rc

# Generate reverse complement of rf
sample.paths <- mapply(function(x, i) {
  
x[['seqtab_RF_rc']] <- x[['seqtab_RF']]
colnames(x[['seqtab_RF_rc']]) <- rc(colnames(x[['seqtab_RF']]))

x
}, 
x = sample.paths, i = names(sample.paths),
SIMPLIFY = F
)

# Merge sequence tables fr and reverse rf sequences in new list
seq.tables <-  mapply(function(x, i) {
  x[['seqtab']] <- mergeSequenceTables(
    x[['seqtab_FR']],
    x[['seqtab_RF_rc']],
    repeats = "sum" # samples with the same name (so sequence) will be summed together
  )
  # will give this message if samples are summed:
  # Duplicated sample names detected in the sequence table row names.
}, 
x = sample.paths, i = names(sample.paths),
SIMPLIFY = F
)


# look at dimensions of all sequence tables
mapply(function(x, i) {
  dim(x)
}, 
x = seq.tables, i = names(seq.tables),
SIMPLIFY = F
)

# $H7LM7DRX2_L2
# [1]    70 19919
# 
# $H7W7NDRX2_L2
# [1]   70 6502
# 
# $H7WKNDRX2_L2
# [1]    55 15594
# 
# $HC2JWDRX2_L1
# [1]   15 5461
# 
# $HJW2YDRXX_L1
# [1]  115 4826
# 
# $HJW2YDRXX_L2
# [1]   86 5560

# merge sequence tables of different sequencing runs
seqtab <- mergeSequenceTables(tables = seq.tables, repeats = "sum")

dim(seqtab)
#  210 32013


# Remove chimeras ####
# This may remove quite a bit of ASVs, but only a small fraction of your total sequences
threads <- 150
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = threads, verbose = TRUE)
# Identified 12108 bimeras out of 32013 input sequences.
ncol(seqtab.nochim)/ncol(seqtab)
# 0.6217787
summary(rowSums(seqtab.nochim)/rowSums(seqtab))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.9344  0.9740  0.9793  0.9784  0.9843  0.9968
## 63% of ASVs and 98% of reads kept


# Inspect ASV length distribution
table(nchar(colnames(seqtab.nochim)))
ASV_length <- as.data.frame(table(nchar(colnames(seqtab.nochim))))
write_tsv(ASV_length,
          file = paste0(wd, "Output/ASV_length_distribution.txt"))
ASV_length %>% 
  arrange(desc(Freq)) %>% 
  head(n = 10)
# Var1  Freq
# 1   251 17628
# 2   252  1308
# 3   250   451
# 4   300   128
# 5   253    95
# 6   190    49
# 7   204    32
# 8   254    27
# 9   249    24
# 10  185    22


# read numbers per length
table(rep(nchar(colnames(seqtab.nochim)), colSums(seqtab.nochim))) 
read_length <- as.data.frame(table(rep(nchar(colnames(seqtab.nochim)), colSums(seqtab.nochim))))
write_tsv(read_length,
          file = paste0(wd, "Output/reads_length_distribution.txt"))
read_length %>% 
  arrange(desc(Freq)) %>% 
  head(n = 10)
# Var1     Freq
# 1   251 12445952
# 2   250   564128
# 3   252   243301
# 4   253     3636
# 5   249     1061
# 6   254      930
# 7   300      624
# 8   204      254
# 9   190      213
# 10  213      137



# Check unusual sequence lengths
### Check unusual sequence lengths
# subset sequence table to most common sequences of unusual length
seqtab.unus.length <- seqtab.nochim %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "seq") %>% 
  mutate(seq_length = nchar(seq)) %>% 
  filter(seq_length < 249 | seq_length > 254) %>% 
  mutate(seq_count = rowSums(.[rownames(seqtab.nochim)]),
         seq_name = paste0("sq", 1:nrow(.))) %>% 
  mutate(seq_name = paste0(seq_name, ";size=", seq_count)) %>% 
  slice_max(order_by = seq_count, n = 50)

# simple function to write fasta file from data.frame
writeFasta <- function(data, colname_seqname, colname_seq, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum, colname_seqname], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum, colname_seq]))
  }
  fileConn<-file(filename)
  writeLines(fastaLines, fileConn)
  close(fileConn)
}

# save fasta file
writeFasta(seqtab.unus.length, colname_seqname = "seq_name", colname_seq = "seq",
           filename = paste0(wd, "Output/check_", prj, ".fasta"))

### do not blast new sequences but use results from previous analysis
## long sequences (300 bp) associated with Sulfurospirillum -> keep
## short sequences (213 bp, 204 bp, 191 bp) distantly associated with some Desulfobacterium, but has so few reads so don't keep


# Remove potential junk sequences and singletons
# dada does not generate singletons, any singletons are introduced in the merging step
# Adjust range of sequence lengths based on expected length of marker gene fragment and read distribution (here 252 bp)
# to save unmerged reads (which were merged based on alignment) estimate upper threshold of junk:
# shortest trimming - merging overlap: 
seqtab.nochim2 <- seqtab.nochim[, colSums(seqtab.nochim) > 1 & ((nchar(colnames(seqtab.nochim)) >= 249 & nchar(colnames(seqtab.nochim)) <= 254) | nchar(colnames(seqtab.nochim)) == 300)]
dim(seqtab.nochim2) #  210 17232
ncol(seqtab.nochim2)/ncol(seqtab)
# 0.5382813
summary(rowSums(seqtab.nochim2)/rowSums(seqtab))
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.9329  0.9733  0.9788  0.9779  0.9840  0.9967


# Taxonomic classification ####
## making function for flexible assignment choosing different databases and outputting different thresholds
own_fct_assignTax <- function(threads, seqtable, database_loc, database_name, bootstrap_for_testing, path_output) {
  
  ## change option setting to not print in scientific annotation for easier reading
  # get previous option to reset on exit
  default_scipen <- options("scipen")
  on.exit(options(scipen = default_scipen)) # this is run on exiting the function for any reason
  # setting to not displaying scientific format
  options(scipen = 999)
  
  ## check if taxonomy output already exists
  if(exists(paste0(database_name, "_tax"), envir = .GlobalEnv)) {
    ## if it exists (because function was run previously and you just want to test different bootstraps)
    # it simply takes specific tax object and renames it 'tax' in here
    print(paste0("Taxonomy table ", database_name, "_tax already exits, using this."))
    tax <- get(paste0(database_name, "_tax"), envir = .GlobalEnv)     # need the get() function here to use paste for calling of object
  } else {
    ## if it does not exist, run taxonomic assignment
    # assign taxonomy with bootstrap value 0
    print(paste("Assigning taxonomy with ", database_name))
    tax <- assignTaxonomy(
      seqtable,
      database_loc,
      multithread = threads,
      minBoot = 0,
      outputBootstraps = T,
      tryRC = T,
      taxLevels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
    )
    ### return tax object so running the code again without the need of doing the assignment again is possible
    assign(paste0(database_name, "_tax"), tax, envir = .GlobalEnv)
    
  }
  
  ## test different bootstrap values: 70, 80, 90
  ### summary of bootstrap values
  sum_bstr <- apply(tax$boot, 2, summary)
  
  ### function to enter different cutoffs
  bstr_cutoff <- function(cutoff) {
    # save taxonomy table separately
    tax.filt <- tax$tax
    # set all assigned taxonomy NA if bootstrap below cutoff
    tax.filt[tax$boot < cutoff] <- NA
    
    # table how many taxa are unassigned on each level
    taxa_unassigned <- base::as.data.frame(apply(tax.filt, 2, function(x) sum(is.na(x)) )) %>% 
      rownames_to_column()
    colnames(taxa_unassigned) <- c("rank", "ASV")
    
    # table how many reads are unassigned on each level
    read_counts_unassigned <- data.frame()
    for(i in 1:ncol(tax.filt)) {
      output = c(colnames(tax.filt)[i], sum(seqtable[, is.na(tax.filt[, i])]))
      read_counts_unassigned <- rbind(read_counts_unassigned, output)
    }
    colnames(read_counts_unassigned) <- c("rank", "reads_no")
    
    # table percentage of reads unassigned on each level
    read_perc_unassigned <- data.frame()
    for(i in 1:ncol(tax.filt)) {
      output = c(colnames(tax.filt)[i], sum(seqtable[, is.na(tax.filt[, i])])/sum(seqtable))
      read_perc_unassigned <- rbind(read_perc_unassigned, output)
    }
    colnames(read_perc_unassigned) <- c("rank", "reads_perc")
    
    output <- full_join(taxa_unassigned, read_counts_unassigned, by = "rank") %>% 
      full_join(read_perc_unassigned, by = "rank") %>% 
      mutate(across(c(reads_no, reads_perc), as.numeric)) %>% 
      column_to_rownames("rank") %>% 
      rename_with(., ~ paste0(.x, "_bstr", cutoff, recycle0 = T))
    return(output)
  }
  
  # use function created above to test different bootstrap thresholds, 
  # create table as output which lists number of unassigned ASVs and reads on different rank levels
  bstr_test_output <- as.data.frame(sapply(bootstrap_for_testing, bstr_cutoff, simplify = F)) %>% 
    relocate(starts_with("reads_perc")) %>% 
    relocate(starts_with("reads_no")) %>% 
    relocate(starts_with("ASV"))
  
  # print the output to console
  print("ASVs and reads which are unassigned using different bootstrap values")
  print(bstr_test_output)
  
  # request user input which bootstrap to use and control input
  ## need a function for this
  cutoff_input <- function() {
    
    us_in <- NA            # need this to start with for the while loop
    
    # repeat this following loop prompting user input as long as 'us_in' is either NA or not between 1 and 100
    while(is.na(us_in) || !((1 < us_in) & (us_in < 100))) {
      # request user input which bootstrap to use
      user_input_no <- readline("Which bootstrap value should be used? Enter a number between 1 and 100. If pipeline does not continue no correct input was entered. Try again. ")
      us_in <- as.numeric(user_input_no)   # readline input is always text, need to make this numeric. wrapping around readline() does not work!
      
      # check if selected bootstrap was actually tested, if yes ask user if to continue or enter a different value
      if(!(us_in %in% bootstrap_for_testing)) {
        checkContinue <- function() {
          user_input <- readline(paste("Selected bootstrap ", us_in, " was not within the tested bootstraps. Do you want to enter a different bootstrap value? If not, pipeline will continue with given value. (y/n) "))
          if(user_input == "y") {
            # set us_in NA so while loop repeats
            print("Enter new bootstrap value.")
            us_in <- NA
            # need to create some output, just us_in <- NA won't manipulate us_in output
            us_in
          } else if(user_input == "n") {
            us_in
          }
        }
        us_in <- checkContinue()
      }
    }
    
    print(paste(us_in, "is used as bootstrap value."))
    us_in
  }
  
  ## use created function, 'bstr_cutoff_used' is bootstrap cut-off number
  bstr_cutoff_used <- cutoff_input()
  
  # write output for different bootstrap values tested
  ## include actual bootstrap used and format table a bit
  bstr_test_output <- bstr_test_output %>% 
    rownames_to_column(var = "rank") %>% 
    mutate(bstr_used = bstr_cutoff_used)
  write_tsv(x = bstr_test_output, file = paste0(path_output, "tax_assignment_", database_name, "_bootstraps.txt"))
  
  # making final taxonomy and ASV table
  ## save taxonomy table separately
  tax.final <- tax$tax
  ## set all assigned taxonomy NA if bootstrap below cutoff
  tax.final[tax$boot < bstr_cutoff_used] <- NA
  return(tax.final)
  
}

### silva release 138
tax.out_silva <- own_fct_assignTax(threads = 150, seqtable = seqtab.nochim2, database_loc = paste0(dblc, "silva_nr99_v138.1_train_set.fa.gz"),
                                   database_name = "silva", bootstrap_for_testing = list(70, 80, 90), path_output = paste0(wd, "Output/"))
# use 80 as bootstrap value

save.image(tf)


# Write output ####
## string that will be used in all written output, so it has not to be changed everywhere manually
output_name <- prj
output_path <- paste0(wd, "Output/")

## ASV sequences as fasta file
uniquesToFasta(seqtab.nochim2, paste0(output_path, "ASVs_dada2_", output_name, ".fasta")) # sequences of asv

# read written file in again, to get correct ids for sequences
seq_unique <- readDNAStringSet(paste0(output_path, "ASVs_dada2_", output_name, ".fasta"))
ASV_id_seq <- data.frame(ASV_ID = gsub(";.*", "", names(seq_unique)), 
                         seq = as.character(seq_unique))

## ASV table with ASVs in rows and samples in columns, sequence and sequence ID, no taxonomy because there are multiple ones
asv.print <- seqtab.nochim2 %>% 
  t() %>%                                                  # transpose table
  as.data.frame() %>% 
  tibble::rownames_to_column() %>%                         # make sequences as own column
  full_join(ASV_id_seq, by = c("rowname" = "seq")) %>%     # merge ASV ID by sequence
  dplyr::relocate(ASV_ID, rowname) %>%                     # order dataframe
  dplyr::rename("seq" = "rowname")

# full table with sequence and sequence ID
write.table(asv.print, paste0(output_path, "asv_table_with_seq_", output_name, ".txt"), quote = F, sep = "\t", row.names = F)

# table with only ASV IDs, could be imported with ASV IDs as rownames
asv.print %>% 
  select(-c(seq)) %>% 
  write.table(paste0(output_path, "asv_table_", output_name, ".txt"), quote = F, sep = "\t", row.names = F)

# Write asv table for metacoder #
# need to transform sequences to ASV IDs or loading into R will take super long
# seq_table with ASVs in columns and samples in rows
asv.print %>% 
  select(-c(seq)) %>% 
  column_to_rownames("ASV_ID") %>% 
  t() %>% 
  as.data.frame() %>% 
  write.table(paste0(output_path, "seq_table_for_metacoder_", output_name, ".txt"), quote = F, sep = "\t", row.names = T)


## join ASV IDs with different taxonomies
### make function to export taxonomy for different databases
export_taxonomy <- function(database_name) {
  # get taxonomy
  tax.out <- get(paste0("tax.out_", database_name), envir = .GlobalEnv)
  
  # join with ASV ID
  tax.print <- tax.out %>% 
    as.data.frame() %>% 
    tibble::rownames_to_column() %>%                         # make sequences as own column
    full_join(ASV_id_seq, by = c("rowname" = "seq")) %>%     # join with taxonomy
    select(-rowname) %>%                                     # get rid of actual sequence
    dplyr::relocate(ASV_ID)
  
  # write taxonomy output
  ## taxonomy table with ASV IDs, could be imported with ASV IDs as rownames
  write.table(tax.print, paste0(output_path, database_name, "_tax_table_", output_name, ".txt"), quote = F, sep = "\t", row.names = F)
  
  ## taxonomy table for metacoder
  # tax_table with ASVs in rows and ranks of taxonomy in separate columns, similar to above but with ASV IDs as rownames
  # could also use tax_table from above and import with one column used as rownames
  tax.print %>% 
    column_to_rownames("ASV_ID") %>% 
    write.table(paste0(output_path, database_name, "_tax_table_for_metacoder_", output_name, ".txt"), quote = F, sep = "\t", row.names = T)
}

### run function for all taxonomy
for(database in gsub("tax.out_", "", ls(pattern = "tax.out_"))) export_taxonomy(database)




# Get nSeqs summary ####
# source functions for creating nSeq files
source(paste0(sloc, "dada2_nSeq_creation_for_list_notax.R"))

# lib_no needs to contain everything after the default file name which is library specific, so here '_libx' including the '_'
# if the file names are the original ones as only one library is analysed use lib_no = NULL instead
nSeq_all <- mapply(function(x, i) {
  i <- create_nSeq(list_input = x, nSeq_file_input = paste0(wd, "nSeq/nSeqs_", i, "_Bac.txt"))
  
}, 
x = sample.paths, i = names(sample.paths),
SIMPLIFY = F
)

# Warning message 
#1: Expected 1 pieces. Additional pieces discarded in 99 rows [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, ...].
#2: Expected 1 pieces. Additional pieces discarded in 99 rows [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, ...].
#3: Expected 1 pieces. Additional pieces discarded in 99 rows [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, ...].
# is ok and does not influence outcome


# feed into function
nSeq <- amend_nSeq(lib_list_input = nSeq_all)

## calculate how many reads were lost on each step
nSeq_perc <- perc_nSeq(nSeq_input = nSeq)


# Write output ####
## string that will be used in all written output, so it has not to be changed everywhere manually
output_name <- prj

## statistic
write.table(nSeq_indiv_seqlanes, paste0(wd, "Output/nSeq_indiv_lanes_dada2_", output_name, ".txt"), row.names = F, quote = F, sep = "\t")
write.table(nSeq, paste0(wd, "Output/nSeq_dada2_", output_name, ".txt"), row.names = F, quote = F, sep = "\t")
write.table(nSeq_perc, paste0(wd, "Output/nSeq_perc_dada2_", output_name, ".txt"), row.names = F, quote = F, sep = "\t")

# end of script ####
save.image(tf)

# proceed to metacoder import R script
