# This script performs a "Repeated G-Test" on 2013 V1016I data
# Started 28 Mar 2018
# depends on setup_jb.R and meltcurve_analysis.R

# # Required libraries
# library(RVAideMemoire)
# library(dplyr)

# # Load mc.1016.t13 and mc.1016.b13 data
# mc.1016.t13 <- read.csv("../../mc.1016.t13_reduced.csv")
# mc.1016.b13 <- read.csv("../../mc.1016.b13_reduced_expandedBuffer.csv")
# # mc.1016.b13 <- read.csv("mc.1016.b13_reduced.csv")

# rename objs to avoid overwriting
mc.t13 <- mc.1016.t13
mc.b13 <- mc.1016.b13

# Add column for treatment type
# for spray region
mc.t13 <- cbind(mc.t13, zone = rep("treatment"))
mc.t13

# for buffer region
mc.b13 <- cbind(mc.b13, zone = rep("buffer"))
mc.b13

# Combine the datasets
dat <- rbind(mc.b13, mc.t13)
dat

### G test should not be done on zero values. 
# Set minimum observations
.min.obs <- 5
# Subset mc.t13 to remove months with samples size < 5
dat <- subset(dat, n > .min.obs)
dat

# Calculate number of R and S alleles from genotype counts
dat$R <- dat$SR + (2*dat$RR)
dat$S <- dat$SR + (2*dat$SS)
dat

############################################################
# Notes about choosing expected proportions for G.tests
# Options:
# 1 - equal proportions - 0.5 : 0.5 ratio of R:S alleles
# 2 - initial proportions - 0.765377 : 0.234623 ratio of R:S alleles (as of 2/9/18)
# 3 - 2013 mean proportions - 0.7612376 : 0.2387624 ratio of R:S alleles (as of 8/2/17)

# Choice:
# Go with option #2 because we care about how the frequencies 
# are changing over time relative to the beginning of the sampling period.
# This choice is in line three of the Fun.xx functions as ", p=c(meanFreqR, meanFreqS)"
# Code to automatically update initial proportions, these are also the expected proportions for 
# the pooled G-test
# meanFreqR = freq R in April buffer zone + freq R in April spray zone / 2

meanFreqR = (dat$freqR[dat$month==4 & dat$zone=="buffer"] 
             + dat$freqR[dat$month==4 & dat$zone=="treatment"])/2
meanFreqS = 1 - meanFreqR

############################################################
# Individual G-tests
# Ho: Numbers within each expt fit expectations
# Ex: There is a specified proportion (i.e. 1:1 or 3:1) of R:S alleles within each month*trt group
# Use Bonferroni correction for this test

# Functions to calculate individual G's, df's, and p-values
Fun.G = function (Q){             
  G.test(x=c(Q["R"], Q["S"])
         # , p=c(0.765377, 0.234623)
         , p=c(meanFreqR, meanFreqS)
  )$statistic                    
}

Fun.df = function (Q){
  G.test(x=c(Q["R"], Q["S"])
         # , p=c(0.765377, 0.234623)
         , p=c(meanFreqR, meanFreqS)
  )$parameter
}

Fun.p = function (Q){
  G.test(x=c(Q["R"], Q["S"])
         # , p=c(0.765377, 0.234623)
         , p=c(meanFreqR, meanFreqS)
  )$p.value
}


# Calculate proportion of R allele
dat =
  mutate(dat,
         Prop.R = R / (R + S),                         
         G =       apply(dat[c("R", "S")], 1, Fun.G),
         df =      apply(dat[c("R", "S")], 1, Fun.df),
         p.Value = apply(dat[c("R", "S")], 1, Fun.p)
  )
# View data
dat

############################################################
# Heterogeneity G-test
# Note: This is the test to report for these data
# Ho: Relative proportions are equal across different experiments
# Ex: The proportion of R alleles is the same in different month*trt groups

# Create a data matrix to run G-test for heterogeneity
Data.matrix = as.matrix(dat[c("R", "S")])      
Data.matrix                                     

# Heterogeneity
gtest.2013 <- G.test(Data.matrix)  

####### Report this #######
gtest.2013    
###########################

############################################################
# Pooled G-test
# Ho: Pooled data fit expectations
# Ex: The number of R and S alles summed across month*trt group is equal to the expected proportions

# Set up data for pooled G-test
Total.R = sum(dat$R)                           
Total.S = sum(dat$S)                           

observed = c(Total.R, Total.S)
expected = c(meanFreqR, meanFreqS)

G.test(x=observed, p=expected)

############################################################
# Total G-test
# Ho: Data from individual experiments fit expectations

# Set up data for total G-test
Total.G  = sum(dat$G)                          
Total.df = sum(dat$df)

# Run 
Total.G                                       
Total.df

pchisq(Total.G,
       df= Total.df,
       lower.tail=FALSE)






