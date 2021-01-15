# Run IQTmosq Project

# load libraries, functions, and themes
source("setup.R")

# load data
source("load_data.R")

# run data analysis scripts
source("run_prepData.R")                    # merges data into one dataframe
## slow
source("run_meltcurve_analysis.R")          # genotype/haplotype construction
#browser()
source("run_Haplotype_Imputation.R")        # create obj for plot_kdrHaps.R
source("run_repeatedGtest_2013_V1016I.R")   # report gtest.2013
source("run_repeatedGtest_2014_V1016I.R")   # report gtest.2014
source("run_kdr_testHWE.R")                 # HWE analysis  - errors due to fixation of one allele in a given month

# create tables and figures for manuscript
# Table 1 - no script
# Table 2 
source("run_LD_V410L-V1016I.R") # linkage disequilibrium analysis
# Figure 1 - no script
# Figure 2 - no script
# Figure 3 
source("plot_kdrHaps.R")
# Figure 4 - no script
# Figure 5 & Sup Figure 1
source("./WFABC/run_pipeline.R", chdir=T)
# Figure 6 - 
# ***add Michael's script here***
# Figure 7 - 
source('plot_kdrZones.R')
# Sup Table 1 - no script
# Sup Figure 2 
# ***add Michael's script here***
