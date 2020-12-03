# Set up workspace for analysis

# Clear working environment -----
rm(list = ls())

# Load libraries -----
library(gdata)
options(gsubfn.engine = "R") # prevents R from stalling while loading sqldf library
library(sqldf)
library(dplyr)
library(ggplot2)
library(ggthemes)

# Source functions -----
source("./loop_MeltCurve_catFiles.R")
# Load function to convert melt curve output
source("./function_convertMergedMeltCurve.R")
# Source functions to create dataframe of genotype/haplotype counts
source("./function_mc.1016.R")
source("./function_mc.1534.R")
source("./function_mc.410.R")
source("./function_mc.haps.R")

# Source ggplot themes -----
source("./ggplotTheme_Vella.R")
