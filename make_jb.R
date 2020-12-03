# Run IQTmosq Project

# load libraries, functions, and themes
source("setup_jb.R")

# load data
source("load_data_jb.R")

# run data analysis scripts
source("meltcurve_analysis.R")  # genotype/haplotype construction
source("run_Haplotype_Imputation.R")  # create obj for plot_kdrHaps.R
source("run_repeatedGtest_2013_V1016I.R") # repeated g-test for 2013
source("run_repeatedGtest_2014_V1016I.R") # repeated g-test for 2014
source("run_kdr_testHWE.R")     # HWE analysis 

# create tables and figures for manuscript
# Table 1 - no script
# Table 2 
source("run_LD_V410L-V1016I.R") # linkage disequilibrium analysis
# Figure 1 - no script
# Figure 2 - no script
# Figure 3 
source("plot_kdrHaps.R")
# Figure 4 - no script
# Figure 5
source("./WFABC/run_WFABC_pipeline.sh") # ***needs updating***
# Figure 6 - 
# ***add Michael's script here***
# Figure 7 - ***replace with xian's updated script***
source("plot_kdrZones13.R")
source("plot_kdrZones14.R")
# Sup Table 1 - no script
# Sup Figure 1
source("./WFABC/run_WFABC_pipeline.sh") # ***needs updating***
# Sup Figure 2 
# ***add Michael's script here***
