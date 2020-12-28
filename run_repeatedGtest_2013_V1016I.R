# This script performs a "Repeated G-Test" on 2013 V1016I data
# Started 28 Mar 2018
# depends on setup_jb.R and meltcurve_analysis.R

dat <- mc.1016.zone %>%
  filter(year==2013)

### G test should not be done on zero values. 
# Set minimum observations
.min.obs <- 5
# Subset mc.t13 to remove months with samples size < 5
dat <- subset(dat, n > .min.obs)

# Calculate number of R and S alleles from genotype counts
dat$R <- dat$SR + (2*dat$RR)
dat$S <- dat$SR + (2*dat$SS)

############################################################
# Calculate mean initial frequencies for R and S
# meanFreqR = freq R in April buffer zone + freq R in April spray zone / 2

meanFreqR = (dat$freqR[dat$month=="04" & dat$zone=="buffer"] 
             + dat$freqR[dat$month=="04" & dat$zone=="treatment"])/2
meanFreqS = 1 - meanFreqR

############################################################
# Individual G-tests
# Ho: Numbers within each expt fit expectations
# Ex: There is a specified proportion (i.e. 1:1 or 3:1) of R:S alleles within each month*trt group
# Use Bonferroni correction for this test

# Functions to calculate individual G's, df's, and p-values
Fun.G = function (Q){             
  G.test(x=c(Q["R"], Q["S"])
         , p=c(meanFreqR, meanFreqS)
  )$statistic                    
}

Fun.df = function (Q){
  G.test(x=c(Q["R"], Q["S"])
         , p=c(meanFreqR, meanFreqS)
  )$parameter
}

Fun.p = function (Q){
  G.test(x=c(Q["R"], Q["S"])
         , p=c(meanFreqR, meanFreqS)
  )$p.value
}


# Calculate proportion of R allele
dat =
  mutate(as.data.frame(dat),
         Prop.R = R / (R + S),                         
         G =       apply(dat[c("R", "S")], 1, Fun.G),
         df =      apply(dat[c("R", "S")], 1, Fun.df),
         p.Value = apply(dat[c("R", "S")], 1, Fun.p)
  )

############################################################
# Heterogeneity G-test
# Ho: Relative proportions are equal across different experiments
# Ex: The proportion of R alleles is the same in different month*trt groups

# Create a data matrix to run G-test for heterogeneity
Data.matrix = as.matrix(dat[c("R", "S")])      

# Heterogeneity
gtest.2013 <- G.test(Data.matrix)  

####### Report this #######
print(gtest.2013)    
###########################

