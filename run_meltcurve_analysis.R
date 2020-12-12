# This script will apply meltcurve functions to kdrData and create objects required for downstream analysis
# requires run_prepData.R to create kdrData object

### Parse data for allele frequency analysis -----
# Parse data by year and/or month to use for allele frequency analysis

# Select based on year -----
dat.mosq <- sqldf("Select *, strftime('%Y', newDate) as year from kdrData where newDate between '2000-01-01' and '2017-12-31'")
# For haplotypes at all years - requires function_mc.haps.R -----
mc.haps.yr <- (
    dat.mosq
    ##! note: year column is now first (instead of last)
    %>% group_by(year)
    %>% group_modify(~mc.haps(.x))
)

# Select based on month from each year -----
# Select based on month from all 2000 data
mc.dat <- sqldf("Select *, strftime('%Y-%m', newDate) as YrMon from kdrData where newDate between '2000-01-01' and '2017-12-31'") 
## make sure all yr.mon combos are present
yr.mon <- with(
    expand.grid(mon=1:12, yr=2000:2017), 
    data.frame(YrMon = sprintf('%d-%02d', yr, mon), month=mon, year=yr)
)
mc.dat <- (
    merge(yr.mon, mc.dat, all=T, by='YrMon')
    %>% group_by(YrMon)
)    

# For 1016 at all month+year combos -----
# mc.1016.byMo - For 1016 locus by month & year -----
mc.1016.byMo <- group_modify(mc.dat, ~mc.1016(.x))
# write df to csv
#write.csv(mc.1016.byMo, file = "./data/mc.1016.byMo.csv", row.names = FALSE)
# For 1534 at  all month+year combos -----
# mc.1534.byMo - For 1534 locus by month & year -----
mc.1534.byMo <- group_modify(mc.dat, ~mc.1534(.x))
# write df to csv
#write.csv(mc.1534.byMo, file = "./data/mc.1534.byMo.csv", row.names = FALSE)

# Select all from treatment zone -----
trt <- sqldf("Select * from kdrData where project_code is 'treatment'")
# Select based on 2013 month from treatment zone
jan2013t <- sqldf("Select * from trt where newDate between '2013-01-01' and '2013-01-31'")
feb2013t <- sqldf("Select * from trt where newDate between '2013-02-01' and '2013-02-31'")
mar2013t <- sqldf("Select * from trt where newDate between '2013-03-01' and '2013-03-31'")
apr2013t <- sqldf("Select * from trt where newDate between '2013-04-01' and '2013-04-31'")
may2013t <- sqldf("Select * from trt where newDate between '2013-05-01' and '2013-05-31'")
jun2013t <- sqldf("Select * from trt where newDate between '2013-06-01' and '2013-06-31'")
jul2013t <- sqldf("Select * from trt where newDate between '2013-07-01' and '2013-07-31'")
aug2013t <- sqldf("Select * from trt where newDate between '2013-08-01' and '2013-08-31'")
sep2013t <- sqldf("Select * from trt where newDate between '2013-09-01' and '2013-09-31'")
oct2013t <- sqldf("Select * from trt where newDate between '2013-10-01' and '2013-10-31'")
# Select based on 2014 month from treatment zone
jan2014t <- sqldf("Select * from trt where newDate between '2014-01-01' and '2014-01-31'")
feb2014t <- sqldf("Select * from trt where newDate between '2014-02-01' and '2014-02-31'")
mar2014t <- sqldf("Select * from trt where newDate between '2014-03-01' and '2014-03-31'")
apr2014t <- sqldf("Select * from trt where newDate between '2014-04-01' and '2014-04-31'")
may2014t <- sqldf("Select * from trt where newDate between '2014-05-01' and '2014-05-31'")
jun2014t <- sqldf("Select * from trt where newDate between '2014-06-01' and '2014-06-31'")
jul2014t <- sqldf("Select * from trt where newDate between '2014-07-01' and '2014-07-31'")
aug2014t <- sqldf("Select * from trt where newDate between '2014-08-01' and '2014-08-31'")
sep2014t <- sqldf("Select * from trt where newDate between '2014-09-01' and '2014-09-31'")
oct2014t <- sqldf("Select * from trt where newDate between '2014-10-01' and '2014-10-31'")

# Select all from buffer zone -----
buff <- sqldf("Select * from kdrData where project_code is 'buffer'") # for literal buffer zone
buff.x <- kdrData[which(kdrData$project_code != "treatment" & !is.na(kdrData$project_code)),] # for expanded buffer

# Select based on 2013 month from buffer zone expanded
jan2013b.x <- buff.x[which(buff.x$newDate >= "2013-01-01" & buff.x$newDate <= "2013-01-31" & !is.na(buff.x$project_code)),]
feb2013b.x <- buff.x[which(buff.x$newDate >= "2013-02-01" & buff.x$newDate <= "2013-02-31" & !is.na(buff.x$project_code)),]
mar2013b.x <- buff.x[which(buff.x$newDate >= "2013-03-01" & buff.x$newDate <= "2013-03-31" & !is.na(buff.x$project_code)),]
apr2013b.x <- buff.x[which(buff.x$newDate >= "2013-04-01" & buff.x$newDate <= "2013-04-31" & !is.na(buff.x$project_code)),]
may2013b.x <- buff.x[which(buff.x$newDate >= "2013-05-01" & buff.x$newDate <= "2013-05-31" & !is.na(buff.x$project_code)),]
jun2013b.x <- buff.x[which(buff.x$newDate >= "2013-06-01" & buff.x$newDate <= "2013-06-31" & !is.na(buff.x$project_code)),]
jul2013b.x <- buff.x[which(buff.x$newDate >= "2013-07-01" & buff.x$newDate <= "2013-07-31" & !is.na(buff.x$project_code)),]
aug2013b.x <- buff.x[which(buff.x$newDate >= "2013-08-01" & buff.x$newDate <= "2013-08-31" & !is.na(buff.x$project_code)),]
sep2013b.x <- buff.x[which(buff.x$newDate >= "2013-09-01" & buff.x$newDate <= "2013-09-31" & !is.na(buff.x$project_code)),]
oct2013b.x <- buff.x[which(buff.x$newDate >= "2013-10-01" & buff.x$newDate <= "2013-10-31" & !is.na(buff.x$project_code)),]
# Select based on 2014 month from buffer zone
jan2014b <- sqldf("Select * from buff where newDate between '2014-01-01' and '2014-01-31'")
feb2014b <- sqldf("Select * from buff where newDate between '2014-02-01' and '2014-02-31'")
mar2014b <- sqldf("Select * from buff where newDate between '2014-03-01' and '2014-03-31'")
apr2014b <- sqldf("Select * from buff where newDate between '2014-04-01' and '2014-04-31'")
may2014b <- sqldf("Select * from buff where newDate between '2014-05-01' and '2014-05-31'")
jun2014b <- sqldf("Select * from buff where newDate between '2014-06-01' and '2014-06-31'")
jul2014b <- sqldf("Select * from buff where newDate between '2014-07-01' and '2014-07-31'")
aug2014b <- sqldf("Select * from buff where newDate between '2014-08-01' and '2014-08-31'")
sep2014b <- sqldf("Select * from buff where newDate between '2014-09-01' and '2014-09-31'")
oct2014b <- sqldf("Select * from buff where newDate between '2014-10-01' and '2014-10-31'")



### Run functions across all years  -----


# For 1016 - treatment zone 2013 & 2014 -----
# For 1016 locus at all months in 2013 - treatment zone
tJan13 <- mc.1016(jan2013t)
tFeb13 <- mc.1016(feb2013t)
tMar13 <- mc.1016(mar2013t)
tApr13 <- mc.1016(apr2013t)
tMay13 <- mc.1016(may2013t)
tJun13 <- mc.1016(jun2013t)
tJul13 <- mc.1016(jul2013t)
tAug13 <- mc.1016(aug2013t)
tSep13 <- mc.1016(sep2013t)
tOct13 <- mc.1016(oct2013t)

# For 1016 locus at all months in 2014 - treatment zone
tJan <- mc.1016(jan2014t)
tFeb <- mc.1016(feb2014t)
tMar <- mc.1016(mar2014t)
tApr <- mc.1016(apr2014t)
tMay <- mc.1016(may2014t)
tJun <- mc.1016(jun2014t)
tJul <- mc.1016(jul2014t)
tAug <- mc.1016(aug2014t)
tSep <- mc.1016(sep2014t)
tOct <- mc.1016(oct2014t)

# For 1016 - buffer zone 2013 & 2014 -----
# For 1016 locus at all months in 2013 - buffer zone expanded
bJan13.x <- mc.1016(jan2013b.x)
bFeb13.x <- mc.1016(feb2013b.x)
bMar13.x <- mc.1016(mar2013b.x)
bApr13.x <- mc.1016(apr2013b.x)
bMay13.x <- mc.1016(may2013b.x)
bJun13.x <- mc.1016(jun2013b.x)
bJul13.x <- mc.1016(jul2013b.x)
bAug13.x <- mc.1016(aug2013b.x)
bSep13.x <- mc.1016(sep2013b.x)
bOct13.x <- mc.1016(oct2013b.x)

# For 1016 locus at all months in 2014 - buffer zone
bJan <- mc.1016(jan2014b)
bFeb <- mc.1016(feb2014b)
bMar <- mc.1016(mar2014b)
bApr <- mc.1016(apr2014b)
bMay <- mc.1016(may2014b)
bJun <- mc.1016(jun2014b)
bJul <- mc.1016(jul2014b)
bAug <- mc.1016(aug2014b)
bSep <- mc.1016(sep2014b)
bOct <- mc.1016(oct2014b)




### Create dataframes for use downstream -----
# Create objects for dfs -----
# Create list of years included in dataframe
year <- c(2000:2017)
# Create list of months included in dataframe - as numeric
month <- c(1:10)

# mc.1016.t13 - For 1016 locus at all months in 2013 - treatment zone -----
dftrt13 <- rbind(tJan13, tFeb13, tMar13, tApr13, tMay13, tJun13, tJul13, tAug13, tSep13, tOct13)
# Add year ID to rows in df rename
mc.1016.t13 <- cbind(month, dftrt13)
# mc.1016t14 - For 1016 locus at all months in 2014 - treatment zone -----
dftrt <- rbind(tJan, tFeb, tMar, tApr, tMay, tJun, tJul, tAug, tSep, tOct)
# Add year ID to rows in df rename
mc.1016.t14 <- cbind(month, dftrt)

# mc.1016.b13.expanded - For 1016 locus at all months in 2013 - buffer zone expanded -----
dfbuff13.x <- rbind(bJan13.x, bFeb13.x, bMar13.x, bApr13.x, bMay13.x, bJun13.x, bJul13.x, bAug13.x, bSep13.x, bOct13.x)
# Add year ID to rows in df rename
mc.1016.b13.expanded <- cbind(month, dfbuff13.x)
# mc.1016.b14 - For 1016 locus at all months in 2014 - buffer zone -----
dfbuff <- rbind(bJan, bFeb, bMar, bApr, bMay, bJun, bJul, bAug, bSep, bOct)
# Add year ID to rows in df rename
mc.1016.b14 <- cbind(month, dfbuff)

# View dataframes -----
# mc.haps.yr
# mc.1016.t13
# mc.1016.t14
# mc.1016.b13.expanded
# mc.1016.b14
# mc.1016.byMo
# mc.1534.byMo
