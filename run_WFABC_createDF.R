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
### Set up the working space ----------------------------
# # clear working environment
# rm(list = ls())

# Set working directory
setwd("/Users/jenbaltz/Dropbox/GouldLab/Project_Mosquito/Database")

# Load data
mc.1016.yr <- read.csv("mc.1016.yr_reduced.csv")
mc.1534.yr <- read.csv("mc.1534.yr_reduced.csv")
sel.1016 <- mc.1016.yr[11:15,]
sel.1534 <- mc.1534.yr[3:9,]
mc.1534.byMo <- read.csv("mc.1534.MoYr_reduced_byMonth.csv")
mc.1016.byMo <- read.csv("mc.1016.MoYr_reduced_byMonth.csv")

# restrict byMo dfs for period of selection
mc.1534.byMo.sel <- mc.1534.byMo[34:120,]
mc.1016.byMo.sel <- mc.1016.byMo[34:120,]

### Function to create df for one or two loci by year
create.WFABC.df <- function(df1, df2 = NULL, n_loci, filename, gens_per_yr = 12) {
  if(!is.null(df2)) df1 <- df1 + df2
  
  # convert certain columns to numeric
  df1$RR <- as.numeric(as.character(df1$RR))
  df1$SR <- as.numeric(as.character(df1$SR))
  df2$RR <- as.numeric(as.character(df2$RR))
  df2$SR <- as.numeric(as.character(df2$SR))
  
  # calculate variables
  n_timepts <- nrow(df1)
  list_timepts <- c(1:n_timepts)
  gens_timepts <- gens_per_yr*list_timepts
  n_chr_timepts <- 2*df1$n
  n_R_timepts <- ((2*df1$RR) + df1$SR)
  n_chr_timepts2 <- 2*df2$n
  n_R_timepts2 <- ((2*df2$RR) + df2$SR)
  
  # write txt file
  sink(filename)
  cat(n_loci, n_timepts, "\n")
  cat(gens_timepts, "\n")
  cat(n_chr_timepts, "\n")
  cat(n_R_timepts, "\n")
  cat(n_chr_timepts2, "\n")
  cat(n_R_timepts2)
  sink()
}

### Function to create df for one or two loci by month
create.WFABC.df.byMo <- function(df1, df2 = NULL, n_loci, filename, gens_per_yr = 12) {
  if(!is.null(df2)) df1 <- df1 + df2
  
  # calculate variables
  n_timepts <- nrow(df1)
  list_timepts <- c(1:n_timepts)
  n_chr_timepts <- 2*df1$n
  n_R_timepts <- ((2*df1$RR) + df1$SR)
  n_chr_timepts2 <- 2*df2$n
  n_R_timepts2 <- ((2*df2$RR) + df2$SR)
  
  # write txt file
  sink(filename)
  cat(n_loci, n_timepts, "\n")
  cat(list_timepts, "\n")
  cat(n_chr_timepts, "\n")
  cat(n_R_timepts, "\n")
  cat(n_chr_timepts2, "\n")
  cat(n_R_timepts2)
  sink()
}

### run function to create multiple dfs for analysis 
create.WFABC.df(mc.1016.yr, mc.1534.yr, n_loci = 2, filename = "./WFABC/multiple_loci.txt")
# create.WFABC.df(sel.1016, n_loci = 1, filename = "./WFABC/V1016I_underSelection.txt")
# create.WFABC.df(sel.1534, n_loci = 1, filename = "./WFABC/F1534C_underSelection.txt")
# create.WFABC.df.byMo(mc.1016.byMo, n_loci = 1, filename = "./WFABC/V1016I_byMo.txt")
# create.WFABC.df.byMo(mc.1534.byMo, n_loci = 1, filename = "./WFABC/F1534C_byMo.txt")
create.WFABC.df.byMo(mc.1016.byMo, mc.1534.byMo, n_loci = 2, filename = "./WFABC/multiple_loci_byMo.txt")
create.WFABC.df.byMo(mc.1016.byMo.sel, mc.1534.byMo.sel, n_loci = 2, filename = "./WFABC/multiple_loci_byMo_sel.txt")

### Function to remove months where no individuals were sampled
remove.noSamples <- function(filename){
  df <- read.csv(file = filename, header = FALSE, sep = "")
  df <- t(df) # transpose df
  df <- t(df[df[,3] != 0,])
  
  # write txt file
  sink(filename)
  cat("2", ncol(df), "\n")
  cat(df[2,], "\n")
  cat(df[3,], "\n")
  cat(df[4,], "\n")
  cat(df[5,], "\n")
  cat(df[6,], "\n")
  sink()
}

# remove months with no samples
remove.noSamples("./WFABC/multiple_loci_byMo.txt")
remove.noSamples("./WFABC/multiple_loci_byMo_sel.txt")

