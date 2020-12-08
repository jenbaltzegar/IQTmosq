# Create matrix of haplotype probabilities
# Started: 8/1/17
# depends on setup_jb.R and meltcurve_analysis.R

# This matrix assumes that RS does not exist in the population
# Locus 1016 is always written first

# # Required libraries
# library(reshape2)

### Create matrix of haplotype probabilities ------------------------------------------------------------
genos <- c("SSSS", "SSSR", "SSRR", "SRSS", "SRSR", "SRRR", "RRSS", "RRSR", "RRRR")
haps <- c("SS", "SR", "RS", "RR")
HapProbs <- matrix(data = c(1,0,0,0,0.5,0.5,0,0,0,1,0,0,1,0,0,0,0.33,0.33,0,0.33,0,0.5,0,0.5,0,0,0,0,0,0,0,1,0,0,0,1)
       , nrow = 9, ncol = 4, byrow = T, dimnames = list(genos, haps))

### Impute haplotype numbers ------------------------------------------------------------
# 2000
df <- data.frame()
# For each year, multiply the genotype number by the haplotype probability
for(i in 1:9){
  df_New <- mc.haps.yr[1,i] * HapProbs[i,]
  yr <- mc.haps.yr$year[1]
  df.yr <- c(df_New, yr)
  df <- rbind(df, df.yr)
  names(df) <- c("SS", "SR", "RS", "RR", "Year")
}
# df
# Create vector of haplotype sums for year
sum2000 <- colSums(df[, c(1:4)])
# Rename df
haps2000 <- df


# 2001
df <- data.frame()
# For each year, multiply the genotype number by the haplotype probability
for(i in 1:9){
  df_New <- mc.haps.yr[2,i] * HapProbs[i,]
  yr <- mc.haps.yr$year[2]
  df.yr <- c(df_New, yr)
  df <- rbind(df, df.yr)
  names(df) <- c("SS", "SR", "RS", "RR", "Year")
}
# df
# Create vector of haplotype sums for year
sum2001 <- colSums(df[, c(1:4)])
# Rename df
haps2001 <- df


# 2002
df <- data.frame()
# For each year, multiply the genotype number by the haplotype probability
for(i in 1:9){
  df_New <- mc.haps.yr[3,i] * HapProbs[i,]
  yr <- mc.haps.yr$year[3]
  df.yr <- c(df_New, yr)
  df <- rbind(df, df.yr)
  names(df) <- c("SS", "SR", "RS", "RR", "Year")
}
# df
# Create vector of haplotype sums for year
sum2002 <- colSums(df[, c(1:4)])
# Rename df
haps2002 <- df

# 2003
df <- data.frame()
# For each year, multiply the genotype number by the haplotype probability
for(i in 1:9){
  df_New <- mc.haps.yr[4,i] * HapProbs[i,]
  yr <- mc.haps.yr$year[4]
  df.yr <- c(df_New, yr)
  df <- rbind(df, df.yr)
  names(df) <- c("SS", "SR", "RS", "RR", "Year")
}
# df
# Create vector of haplotype sums for year
sum2003 <- colSums(df[, c(1:4)])
# Rename df
haps2003 <- df

# 2004
df <- data.frame()
# For each year, multiply the genotype number by the haplotype probability
for(i in 1:9){
  df_New <- mc.haps.yr[5,i] * HapProbs[i,]
  yr <- mc.haps.yr$year[5]
  df.yr <- c(df_New, yr)
  df <- rbind(df, df.yr)
  names(df) <- c("SS", "SR", "RS", "RR", "Year")
}
# df
# Create vector of haplotype sums for year
sum2004 <- colSums(df[, c(1:4)])
# Rename df
haps2004 <- df


# 2005
df <- data.frame()
# For each year, multiply the genotype number by the haplotype probability
for(i in 1:9){
  df_New <- mc.haps.yr[6,i] * HapProbs[i,]
  yr <- mc.haps.yr$year[6]
  df.yr <- c(df_New, yr)
  df <- rbind(df, df.yr)
  names(df) <- c("SS", "SR", "RS", "RR", "Year")
}
# df
# Create vector of haplotype sums for year
sum2005 <- colSums(df[, c(1:4)])
# Rename df
haps2005 <- df

# 2006
df <- data.frame()
# For each year, multiply the genotype number by the haplotype probability
for(i in 1:9){
  df_New <- mc.haps.yr[7,i] * HapProbs[i,]
  yr <- mc.haps.yr$year[7]
  df.yr <- c(df_New, yr)
  df <- rbind(df, df.yr)
  names(df) <- c("SS", "SR", "RS", "RR", "Year")
}
# df
# Create vector of haplotype sums for year
sum2006 <- colSums(df[, c(1:4)])
# Rename df
haps2006 <- df


# 2007
df <- data.frame()
# For each year, multiply the genotype number by the haplotype probability
for(i in 1:9){
  df_New <- mc.haps.yr[8,i] * HapProbs[i,]
  yr <- mc.haps.yr$year[8]
  df.yr <- c(df_New, yr)
  df <- rbind(df, df.yr)
  names(df) <- c("SS", "SR", "RS", "RR", "Year")
}
# df
# Create vector of haplotype sums for year
sum2007 <- colSums(df[, c(1:4)])
# Rename df
haps2007 <- df


# 2008
df <- data.frame()
# For each year, multiply the genotype number by the haplotype probability
for(i in 1:9){
  df_New <- mc.haps.yr[9,i] * HapProbs[i,]
  yr <- mc.haps.yr$year[9]
  df.yr <- c(df_New, yr)
  df <- rbind(df, df.yr)
  names(df) <- c("SS", "SR", "RS", "RR", "Year")
}
# df
# Create vector of haplotype sums for year
sum2008 <- colSums(df[, c(1:4)])
# Rename df
haps2008 <- df


# 2009
df <- data.frame()
# For each year, multiply the genotype number by the haplotype probability
for(i in 1:9){
  df_New <- mc.haps.yr[10,i] * HapProbs[i,]
  yr <- mc.haps.yr$year[10]
  df.yr <- c(df_New, yr)
  df <- rbind(df, df.yr)
  names(df) <- c("SS", "SR", "RS", "RR", "Year")
}
# df
# Create vector of haplotype sums for year
sum2009 <- colSums(df[, c(1:4)])
# Rename df
haps2009 <- df

# 2010
df <- data.frame()
# For each year, multiply the genotype number by the haplotype probability
for(i in 1:9){
  df_New <- mc.haps.yr[11,i] * HapProbs[i,]
  yr <- mc.haps.yr$year[11]
  df.yr <- c(df_New, yr)
  df <- rbind(df, df.yr)
  names(df) <- c("SS", "SR", "RS", "RR", "Year")
}
# df
# Create vector of haplotype sums for year
sum2010 <- colSums(df[, c(1:4)])
# Rename df
haps2010 <- df


# 2011
df <- data.frame()
# For each year, multiply the genotype number by the haplotype probability
for(i in 1:9){
  df_New <- mc.haps.yr[12,i] * HapProbs[i,]
  yr <- mc.haps.yr$year[12]
  df.yr <- c(df_New, yr)
  df <- rbind(df, df.yr)
  names(df) <- c("SS", "SR", "RS", "RR", "Year")
}
# df
# Create vector of haplotype sums for year
sum2011 <- colSums(df[, c(1:4)])
# Rename df
haps2011 <- df

# 2012
df <- data.frame()
# For each year, multiply the genotype number by the haplotype probability
for(i in 1:9){
  df_New <- mc.haps.yr[13,i] * HapProbs[i,]
  yr <- mc.haps.yr$year[13]
  df.yr <- c(df_New, yr)
  df <- rbind(df, df.yr)
  names(df) <- c("SS", "SR", "RS", "RR", "Year")
}
# df
# Create vector of haplotype sums for year
sum2012 <- colSums(df[, c(1:4)])
# Rename df
haps2012 <- df

# 2013
df <- data.frame()
# For each year, multiply the genotype number by the haplotype probability
for(i in 1:9){
  df_New <- mc.haps.yr[14,i] * HapProbs[i,]
  yr <- mc.haps.yr$year[14]
  df.yr <- c(df_New, yr)
  df <- rbind(df, df.yr)
  names(df) <- c("SS", "SR", "RS", "RR", "Year")
}
# df
# Create vector of haplotype sums for year
sum2013 <- colSums(df[, c(1:4)])
# Rename df
haps2013 <- df

# 2014
df <- data.frame()
# For each year, multiply the genotype number by the haplotype probability
for(i in 1:9){
  df_New <- mc.haps.yr[15,i] * HapProbs[i,]
  yr <- mc.haps.yr$year[15]
  df.yr <- c(df_New, yr)
  df <- rbind(df, df.yr)
  names(df) <- c("SS", "SR", "RS", "RR", "Year")
}
# df
# Create vector of haplotype sums for year
sum2014 <- colSums(df[, c(1:4)])
# Rename df
haps2014 <- df


# 2015
df <- data.frame()
# For each year, multiply the genotype number by the haplotype probability
for(i in 1:9){
  df_New <- mc.haps.yr[16,i] * HapProbs[i,]
  yr <- mc.haps.yr$year[16]
  df.yr <- c(df_New, yr)
  df <- rbind(df, df.yr)
  names(df) <- c("SS", "SR", "RS", "RR", "Year")
}
# df
# Create vector of haplotype sums for year
sum2015 <- colSums(df[, c(1:4)])
# Rename df
haps2015 <- df


# 2016
df <- data.frame()
# For each year, multiply the genotype number by the haplotype probability
for(i in 1:9){
  df_New <- mc.haps.yr[17,i] * HapProbs[i,]
  yr <- mc.haps.yr$year[17]
  df.yr <- c(df_New, yr)
  df <- rbind(df, df.yr)
  names(df) <- c("SS", "SR", "RS", "RR", "Year")
}
# df
# Create vector of haplotype sums for year
sum2016 <- colSums(df[, c(1:4)])
# Rename df
haps2016 <- df

# 2017
df <- data.frame()
# For each year, multiply the genotype number by the haplotype probability
for(i in 1:9){
  df_New <- mc.haps.yr[18,i] * HapProbs[i,]
  yr <- mc.haps.yr$year[18]
  df.yr <- c(df_New, yr)
  df <- rbind(df, df.yr)
  names(df) <- c("SS", "SR", "RS", "RR", "Year")
}
# df
# Create vector of haplotype sums for year
sum2017 <- colSums(df[, c(1:4)])
# Rename df
haps2017 <- df


### Create table of sums for each haplotype per year ------------------------------------------------------------
sumAll <- rbind(sum2000, sum2001, sum2002, sum2003, sum2004, sum2005, sum2006, sum2007
                , sum2008, sum2009, sum2010, sum2011, sum2012, sum2013, sum2014, sum2015, sum2016, sum2017)

row.names(sumAll) <- (c("2000", "2001", "2002", "2003", "2004", "2005", "2006", "2007", "2008"
                        , "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017"))


### Get the frequency of haps for each year ------------------------------------------------------------
freq2000 <- sum2000/mc.haps.yr$n[1]
freq2001 <- sum2001/mc.haps.yr$n[2]
freq2002 <- sum2002/mc.haps.yr$n[3]
freq2003 <- sum2003/mc.haps.yr$n[4]
freq2004 <- sum2004/mc.haps.yr$n[5]
freq2005 <- sum2005/mc.haps.yr$n[6]
freq2006 <- sum2006/mc.haps.yr$n[7]
freq2007 <- sum2007/mc.haps.yr$n[8]
freq2008 <- sum2008/mc.haps.yr$n[9]
freq2009 <- sum2009/mc.haps.yr$n[10]
freq2010 <- sum2010/mc.haps.yr$n[11]
freq2011 <- sum2011/mc.haps.yr$n[12]
freq2012 <- sum2012/mc.haps.yr$n[13]
freq2013 <- sum2013/mc.haps.yr$n[14]
freq2014 <- sum2014/mc.haps.yr$n[15]
freq2015 <- sum2015/mc.haps.yr$n[16]
freq2016 <- sum2016/mc.haps.yr$n[17]
freq2017 <- sum2017/mc.haps.yr$n[18]

### Create table for frequency of haps for each year ------------------------------------------------------------
freqAll <- rbind(freq2000, freq2001, freq2002, freq2003, freq2004, freq2005, freq2006, freq2007
                , freq2008, freq2009, freq2010, freq2011, freq2012, freq2013, freq2014, freq2015, freq2016, freq2017)

row.names(freqAll) <- (c("2000", "2001", "2002", "2003", "2004", "2005", "2006", "2007", "2008"
                        , "2009", "2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017"))
freqAll <- as.data.frame(freqAll)
year <- 2000:2017
freqAll <- cbind(freqAll, year)

### Calculate 95% Confidence Interval (not for ribbon CIs) ------------------------------------------------------------
## for column 1: SS
df <- data.frame()
for(i in 1:18){
  ci = 1.96 * (sqrt((freqAll[i,1]*(1-freqAll[i,1]))/(2*mc.haps.yr$n[i])))
  df = rbind(df, ci)
}
colnames(df) <- "CI_95"
CI_95.SS <- df
# CI_95.SS

## for column 2: SR
df <- data.frame()
for(i in 1:18){
  ci = 1.96 * (sqrt((freqAll[i,2]*(1-freqAll[i,2]))/(2*mc.haps.yr$n[i])))
  df = rbind(df, ci)
}
colnames(df) <- "CI_95"
CI_95.SR <- df
# CI_95.SR

## for column 3: RS
df <- data.frame()
for(i in 1:18){
  ci = 1.96 * (sqrt((freqAll[i,3]*(1-freqAll[i,3]))/(2*mc.haps.yr$n[i])))
  df = rbind(df, ci)
}
colnames(df) <- "CI_95"
CI_95.RS <- df
# CI_95.RS

## for column 4: RR
df <- data.frame()
for(i in 1:18){
  ci = 1.96 * (sqrt((freqAll[i,4]*(1-freqAll[i,4]))/(2*mc.haps.yr$n[i])))
  df = rbind(df, ci)
}
colnames(df) <- "CI_95"
CI_95.RR <- df
# CI_95.RR

CI_95 <- rbind(CI_95.SS, CI_95.SR, CI_95.RS, CI_95.RR)


### Converting from wide to long format ------------------------------------------------------------
freqAll_long <- melt(freqAll, id.vars="year", variable.name = "Haplotype", value.name = "Frequency")
freqAll_long <- cbind(freqAll_long, CI_95)

# # Subset data based on haplotype
# freqSS <- freqAll_long[freqAll_long$Haplotype=="SS", ]
# freqSR <- freqAll_long[freqAll_long$Haplotype=="SR", ]
# freqRS <- freqAll_long[freqAll_long$Haplotype=="RS", ]
# freqRR <- freqAll_long[freqAll_long$Haplotype=="RR", ]

