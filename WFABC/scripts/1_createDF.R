#!/usr/bin/env Rscript

# This code will create the input file to run WFABC model 

##########################################################
### input file format - copied from WFABC manual
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

# load data
mc.1016.byMo <- read.csv("../data/mc.1016.byMo.csv", header = TRUE)
mc.1534.byMo <- read.csv("../data/mc.1534.byMo.csv", header = TRUE)

# source functions
source("./scripts/function_createDF.R")

##! subset by date
# Restrict dfs for period of selection
mc.1016.byMo.sel <- mc.1016.byMo[124:174,]   # 124:174 = Apr 2010 - Jun 2014 
mc.1534.byMo.sel <- mc.1534.byMo[34:120,]    # 34:120 = Oct 2002 - Dec 2009

# ### Create dfs for analysis -----
create.WFABC.df.byMo(mc.1016.byMo.sel, filename = "./V1016I_byMo_sel.txt")
create.WFABC.df.byMo(mc.1534.byMo.sel, filename = "./F1534C_byMo_sel.txt")

# ### Remove months with no samples - functions sourced from function_createDF.R -----
remove.noSamples.onelocus("./V1016I_byMo_sel.txt")
remove.noSamples.onelocus("./F1534C_byMo_sel.txt")

