# This script will apply meltcurve functions to kdrData and create objects required for downstream analysis
# requires run_prepData.R to create kdrData object

### Parse data for allele frequency analysis -----
# Parse data by year and/or month to use for allele frequency analysis

# Select based on year -----
mosq2000 <- sqldf("Select * from kdrData where newDate between '2000-01-01' and '2000-12-31'")
mosq2001 <- sqldf("Select * from kdrData where newDate between '2001-01-01' and '2001-12-31'")
mosq2002 <- sqldf("Select * from kdrData where newDate between '2002-01-01' and '2002-12-31'")
mosq2003 <- sqldf("Select * from kdrData where newDate between '2003-01-01' and '2003-12-31'")
mosq2004 <- sqldf("Select * from kdrData where newDate between '2004-01-01' and '2004-12-31'")
mosq2005 <- sqldf("Select * from kdrData where newDate between '2005-01-01' and '2005-12-31'")
mosq2006 <- sqldf("Select * from kdrData where newDate between '2006-01-01' and '2006-12-31'")
mosq2007 <- sqldf("Select * from kdrData where newDate between '2007-01-01' and '2007-12-31'")
mosq2008 <- sqldf("Select * from kdrData where newDate between '2008-01-01' and '2008-12-31'")
mosq2009 <- sqldf("Select * from kdrData where newDate between '2009-01-01' and '2009-12-31'")
mosq2010 <- sqldf("Select * from kdrData where newDate between '2010-01-01' and '2010-12-31'")
mosq2011 <- sqldf("Select * from kdrData where newDate between '2011-01-01' and '2011-12-31'")
mosq2012 <- sqldf("Select * from kdrData where newDate between '2012-01-01' and '2012-12-31'")
mosq2013 <- sqldf("Select * from kdrData where newDate between '2013-01-01' and '2013-12-31'")
mosq2014 <- sqldf("Select * from kdrData where newDate between '2014-01-01' and '2014-12-31'")
mosq2015 <- sqldf("Select * from kdrData where newDate between '2015-01-01' and '2015-12-31'")
mosq2016 <- sqldf("Select * from kdrData where newDate between '2016-01-01' and '2016-12-31'")
mosq2017 <- sqldf("Select * from kdrData where newDate between '2017-01-01' and '2017-12-31'")

# Select based on month from each year -----
# Select based on month from all 2000 data
mc.dat <- sqldf("Select *, strftime('%Y-%m', newDate) as YrMon from kdrData where newDate between '2000-01-01' and '2017-12-31'") 
## make sure all yr.mon combos are present
yr.mon <- with(
    expand.grid(mon=1:12, yr=2000:2017), 
    data.frame(YrMon = sprintf('%d-%02d', yr, mon), month=mon, year=yr)
)
mc.dat <- merge(yr.mon, mc.dat, all=T, by='YrMon')

# For 1016 at all month+year combos -----
# mc.1016.byMo - For 1016 locus by month & year -----
mc.1016.byMo <- ddply(mc.dat, .variables='YrMon', mc.1016)
# write df to csv
#write.csv(mc.1016.byMo, file = "./data/mc.1016.byMo.csv", row.names = FALSE)
# For 1534 at  all month+year combos -----
# mc.1534.byMo - For 1534 locus by month & year -----
mc.1534.byMo <- ddply(mc.dat, 'YrMon', mc.1534)
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

# For haplotypes at all years - requires function_mc.haps.R -----
h2000 = mc.haps(mosq2000)
h2001 = mc.haps(mosq2001)
h2002 = mc.haps(mosq2002)
h2003 = mc.haps(mosq2003)
h2004 = mc.haps(mosq2004)
h2005 = mc.haps(mosq2005)
h2006 = mc.haps(mosq2006)
h2007 = mc.haps(mosq2007)
h2008 = mc.haps(mosq2008)
h2009 = mc.haps(mosq2009)
h2010 = mc.haps(mosq2010)
h2011 = mc.haps(mosq2011)
h2012 = mc.haps(mosq2012)
h2013 = mc.haps(mosq2013)
h2014 = mc.haps(mosq2014)
h2015 = mc.haps(mosq2015)
h2016 = mc.haps(mosq2016)
h2017 = mc.haps(mosq2017)


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

MonthYear <- c("1-2000", "2-2000", "3-2000", "4-2000", "5-2000", "6-2000", "7-2000", "8-2000", "9-2000", "10-2000", "11-2000", "12-2000"
               , "1-2001", "2-2001", "3-2001", "4-2001", "5-2001", "6-2001", "7-2001", "8-2001", "9-2001", "10-2001", "11-2001", "12-2001"
               , "1-2002", "2-2002", "3-2002", "4-2002", "5-2002", "6-2002", "7-2002", "8-2002", "9-2002", "10-2002", "11-2002", "12-2002"
               , "1-2003", "2-2003", "3-2003", "4-2003", "5-2003", "6-2003", "7-2003", "8-2003", "9-2003", "10-2003", "11-2003", "12-2003"
               , "1-2004", "2-2004", "3-2004", "4-2004", "5-2004", "6-2004", "7-2004", "8-2004", "9-2004", "10-2004", "11-2004", "12-2004"
               , "1-2005", "2-2005", "3-2005", "4-2005", "5-2005", "6-2005", "7-2005", "8-2005", "9-2005", "10-2005", "11-2005", "12-2005"
               , "1-2006", "2-2006", "3-2006", "4-2006", "5-2006", "6-2006", "7-2006", "8-2006", "9-2006", "10-2006", "11-2006", "12-2006"
               , "1-2007", "2-2007", "3-2007", "4-2007", "5-2007", "6-2007", "7-2007", "8-2007", "9-2007", "10-2007", "11-2007", "12-2007"
               , "1-2008", "2-2008", "3-2008", "4-2008", "5-2008", "6-2008", "7-2008", "8-2008", "9-2008", "10-2008", "11-2008", "12-2008"
               , "1-2009", "2-2009", "3-2009", "4-2009", "5-2009", "6-2009", "7-2009", "8-2009", "9-2009", "10-2009", "11-2009", "12-2009"
               , "1-2010", "2-2010", "3-2010", "4-2010", "5-2010", "6-2010", "7-2010", "8-2010", "9-2010", "10-2010", "11-2010", "12-2010"
               , "1-2011", "2-2011", "3-2011", "4-2011", "5-2011", "6-2011", "7-2011", "8-2011", "9-2011", "10-2011", "11-2011", "12-2011"
               , "1-2012", "2-2012", "3-2012", "4-2012", "5-2012", "6-2012", "7-2012", "8-2012", "9-2012", "10-2012", "11-2012", "12-2012"
               , "1-2013", "2-2013", "3-2013", "4-2013", "5-2013", "6-2013", "7-2013", "8-2013", "9-2013", "10-2013", "11-2013", "12-2013"
               , "1-2014", "2-2014", "3-2014", "4-2014", "5-2014", "6-2014", "7-2014", "8-2014", "9-2014", "10-2014", "11-2014", "12-2014"
               , "1-2015", "2-2015", "3-2015", "4-2015", "5-2015", "6-2015", "7-2015", "8-2015", "9-2015", "10-2015", "11-2015", "12-2015"
               , "1-2016", "2-2016", "3-2016", "4-2016", "5-2016", "6-2016", "7-2016", "8-2016", "9-2016", "10-2016", "11-2016", "12-2016"
               , "1-2017", "2-2017", "3-2017", "4-2017", "5-2017", "6-2017", "7-2017", "8-2017", "9-2017", "10-2017", "11-2017", "12-2017")


# mc.haps.yr - For haplotypes -----
dfHaps <- rbind(h2000, h2001, h2002, h2003, h2004, h2005, h2006
                , h2007, h2008, h2009, h2010, h2011, h2012, h2013
                , h2014, h2015, h2016, h2017)
# Add year ID to rows in df rename
mc.haps.yr <- cbind(dfHaps, year)

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
