# Set up workspace for analysis

# Clear working environment -----
rm(list = ls())

# Load libraries -----
library(gdata)
options(gsubfn.engine = "R") # prevents R from stalling while loading sqldf library
library(sqldf)
library(dplyr)
library(reshape2) # for haplotype imputation
library(ggplot2)
# library(RColorBrewer)
library(ggthemes)
library(RVAideMemoire) # for repeated g-test
library(genetics) # for LD analysis


# Source functions -----
source("loop_MeltCurve_catFiles.R") # concatenates various meltcurve rep files
source("function_convertMergedMeltCurve.R") # converts meltcurve output to genotypes
# Source functions to create dataframe of genotype/haplotype counts
source("function_mc.1016.R")
source("function_mc.1534.R")
source("function_mc.410.R")
source("function_mc.haps.R")

# Source ggplot themes -----
source("ggplotTheme_Vella.R")

