# This script will run data through the WFABC pipeline
## Call w with chdir=T
## all paths relative to current dir (WFABC)
# source("./scripts/function_createDF.R")
source("./scripts/function_createDF_temp.R")
# source("./scripts/function_convertPosteriors.R")
source('./scripts/1_createDF.R')
system('./scripts/2_runWFABC.sh')
# source('./scripts/4_analyze_byMo_sel_converted.R')
source("./scripts/4_analyze_temp.R")
