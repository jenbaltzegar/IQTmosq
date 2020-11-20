#!/usr/bin/env Rscript

# This code will create the input file to run WFABC model 

##########################################################
### input file format
# The first line you define the number of loci & number of time points.
# The second line indicates the time of the 12 points. 
#   The important thing here is that the difference between 
#   the values corresponds to the number of generations between 
#   the time points. In this example the first time point is labeled 
#   as generation “1”, and the second time point is 9 generations 
#   later (“10”). Adding or subtracting a constant to all the values 
#   won’t change the results. They can also be negative!
# The rest of the file has to contain a pair of lines for each locus 
#   included in the analysis. For each locus the first line corresponds to the 
#   sample size at each time point (in number of chromosomes, so twice 
#   the number of individuals for diploids), and the second line to the 
#   number of A alleles at each time point.
##########################################################

### Set up the working space -----
# clear working environment
rm(list = ls())

# Set working directory
setwd("~/Documents/jen.temp/jfbaltz_kdr/WFABC")

# Load data
mc.1016.yr <- read.csv("./data/mc.1016.yr_reduced.csv")
mc.1534.yr <- read.csv("./data/mc.1534.yr_reduced.csv")
mc.1534.byMo <- read.csv("./data/mc.1534.MoYr_reduced_byMonth.csv")
mc.1016.byMo <- read.csv("./data/mc.1016.MoYr_reduced_byMonth.csv")

# Restrict dfs for period of selection
sel.1016 <- mc.1016.yr[11:15,]
sel.1534 <- mc.1534.yr[3:9,]
mc.1016.byMo.sel <- mc.1016.byMo[124:174,]   # 124:174 = Apr 2010 - Jun 2014 
mc.1534.byMo.sel <- mc.1534.byMo[34:120,]    # 34:120 = Oct 2002 - Dec 2009

### Source the functions -----
source("./scripts/function_createDF.R")

# ### Create dfs for analysis -----
# create.WFABC.df(mc.1016.yr, mc.1534.yr, n_loci = 2, filename = "./multiple_loci.txt")
create.WFABC.df(sel.1016, n_loci = 1, filename = "./V1016I_underSelection.txt")
create.WFABC.df(sel.1534, n_loci = 1, filename = "./F1534C_underSelection.txt")
# create.WFABC.df.byMo(mc.1016.byMo, filename = "./V1016I_byMo.txt")
# create.WFABC.df.byMo(mc.1534.byMo, filename = "./F1534C_byMo.txt")
create.WFABC.df.byMo(mc.1016.byMo.sel, filename = "./V1016I_byMo_sel.txt")
create.WFABC.df.byMo(mc.1534.byMo.sel, filename = "./F1534C_byMo_sel.txt")
# create.WFABC.df.byMo(mc.1016.byMo, mc.1534.byMo, n_loci = 2, filename = "./multiple_loci_byMo.txt")
# create.WFABC.df.byMo(mc.1016.byMo.sel, mc.1534.byMo.sel, n_loci = 2, filename = "./multiple_loci_byMo_sel.txt")


# ### Remove months with no samples - functions sourced from function_WFABC_createDF.R -----
# remove.noSamples.twolocus("./multiple_loci.txt")
remove.noSamples.onelocus("./V1016I_underSelection.txt")
remove.noSamples.onelocus("./F1534C_underSelection.txt")
# remove.noSamples.onelocus("./V1016I_byMo.txt")
# remove.noSamples.onelocus("./F1534C_byMo.txt")
remove.noSamples.onelocus("./V1016I_byMo_sel.txt")
remove.noSamples.onelocus("./F1534C_byMo_sel.txt")
# remove.noSamples.twolocus("./multiple_loci_byMo.txt")
# remove.noSamples.twolocus("./multiple_loci_byMo_sel.txt")

