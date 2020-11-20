# Script from Michael to compare two datasets

# load libraries
library(dplyr)
library(reshape2)

# set working directory
setwd("~/Dropbox/GouldLab/Project_Mosquito/Database")

# Michael's code from email received 5/28/2020
genotype_names <- c("RR","SR","SS")
dat <- read.csv('mc.1534.MoYr_reduced_byMonth.csv')
full_dat <- read.csv('kdrData_reduced_2020-05-27.csv')

#parse to add month and year columns
full_dat <- cbind(full_dat,colsplit(as.character(full_dat$Date),"/",c("month","day","year")))
full_dat$year <- full_dat$year + 2000
head(full_dat)

month_years <- expand.grid("month"=seq(1,12),"year"=seq(2000,2017))
#inefficient aggregation of genotypes by month
genodat <- bind_rows(apply(month_years,1,function(x){
  monthly <-  subset(full_dat,month==x[1] & year==x[2],select="F1534C_converted")
  data.frame("RR"=sum(monthly$F1534C_converted=="RR"),
             "SR"=sum(monthly$F1534C_converted=="SR"),
             "SS"=sum(monthly$F1534C_converted=="SS"))
})
)
full_dat_month <- cbind(genodat,month_years)
#discrepencies?:
cbind(month_years,full_dat_month[,genotype_names] - dat[,genotype_names])
