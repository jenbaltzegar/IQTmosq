### Melt Curve Analysis
### File Started: 31 May 2017
### Major updates: 2 Dec 2020

# Prep data ---------------------------------------------------------------
# Convert Date to readable form for R
kdrData$newDate <- as.character(as.Date(as.character(kdrData$Date), format = "%m/%d/%Y"))
# Note: If year starts with 00 after conversion, then go back and save IsolationLog_IQT-Mosq.csv
# with Date column in the correct format (MM/DD/YYYY). Weird formating issues happen with .csv files

# Left join masterdec GIS and neighborhood info onto kdrData
kdrData <- merge(x=kdrData, y = masterdec[, c("X", "Y", "NEIGHBORHO", "LOC_CODE")]
              , by.x = "Location_Code", by.y = "LOC_CODE" , all.x = TRUE)


# Merge melt curve data with isolation log to create df "kdrData"
# Rename column Genotype in each df
colnames(merged1016_rep1)[4] <- "V1016I_rep1"
colnames(merged1016_rep2)[4] <- "V1016I_rep2"
colnames(merged1534_rep1)[4] <- "F1534C_rep1"
colnames(merged1534_rep2)[4] <- "F1534C_rep2"
colnames(merged410_rep1)[4] <- "V410L_rep1"
colnames(merged410_rep2)[4] <- "V410L_rep2"

# Left join merged data to kdrData
kdrData <- merge(x = kdrData, y = merged1016_rep1[, c("mosquito_id", "V1016I_rep1")]
                 , by = "mosquito_id", all.x = TRUE)
kdrData <- merge(x = kdrData, y = merged1016_rep2[, c("mosquito_id", "V1016I_rep2")]
                 , by = "mosquito_id", all.x = TRUE)
kdrData <- merge(x = kdrData, y = merged1534_rep1[, c("mosquito_id", "F1534C_rep1")]
                 , by = "mosquito_id", all.x = TRUE)
kdrData <- merge(x = kdrData, y = merged1534_rep2[, c("mosquito_id", "F1534C_rep2")]
                 , by = "mosquito_id", all.x = TRUE)
kdrData <- merge(x = kdrData, y = merged410_rep1[, c("mosquito_id", "V410L_rep1")]
                 , by = "mosquito_id", all.x = TRUE)
kdrData <- merge(x = kdrData, y = merged410_rep2[, c("mosquito_id", "V410L_rep2")]
                 , by = "mosquito_id", all.x = TRUE)


# Add columns with verified replicated genotype or put error 
kdrData$V1016I <- ifelse(kdrData$V1016I_rep1 == kdrData$V1016I_rep2
                       , as.character(kdrData$V1016I_rep1)
                       , "error")
kdrData$F1534C <- ifelse(kdrData$F1534C_rep1 == kdrData$F1534C_rep2
                         , as.character(kdrData$F1534C_rep1)
                         , "error")
kdrData$V410L <- ifelse(kdrData$V410L_rep1 == kdrData$V410L_rep2
                         , as.character(kdrData$V410L_rep1)
                         , "error")


# Left join Zone information to kdrDate
kdrData <- merge(x = kdrData, y = exptZone, by.x = "Location_Code", by.y = "location_code", all.x = TRUE)


# Create and save files with missing project_code
# Those files with no Zone (i.e. project code) specified
noZone <- kdrData[is.na(kdrData$project_code),]
# write.csv(noZone, "noZone.csv", row.names = FALSE)

# Those files with no Zone but with known location codes
noZoneNoUnknowns <- noZone[noZone$Location_Code != "unknown",]
# write.csv(noZoneNoUnknowns, "noZoneNoUnknowns.csv", row.names = FALSE)

####
# Key Code: Remove unknown location codes from kdrData and save file
####
# "Unknown" location codes are typically from mosquitoes collected by Steve Stoddard outside of IQT
kdrData <- kdrData[kdrData$Location_Code != "unknown", ]
write.csv(kdrData,"~/Dropbox/GouldLab/Project_Mosquito/Database/kdrData_all_reduced.csv", row.names = F)


# Convert melt curve output to readable genotype and haplotype
# Uses function_convertMergedMeltCurve.R

# Run function for each column
kdrData$V1016I_converted <- convert1016(kdrData)
kdrData$F1534C_converted <- convert1534(kdrData)
kdrData$V410L_converted <- convert410(kdrData)

# Create haplotype for loci 1016 and 1534 - 410 is not included here
kdrData$haplotype <- paste0(kdrData$V1016I_converted, kdrData$F1534C_converted)

# ### Save file
# # Write to file
# write.csv(kdrData, paste0("~/Dropbox/GouldLab/Project_Mosquito/Database/kdrData_reduced_archived/kdrData_reduced_", Sys.Date(), ".csv"), row.names = F)
# write.csv(kdrData, "~/Dropbox/GouldLab/Project_Mosquito/Database/kdrData_reduced.csv", row.names = F)



# Parse data for allele frequency analysis --------------------------------
### Parse data by year and/or month to use for allele frequency analysis
# Keep only those rows with project code

# Select based on year
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

# Select based on month from all 2013 data
jan2013 <- sqldf("Select * from mosq2013 where newDate between '2013-01-01' and '2013-01-31'")
feb2013 <- sqldf("Select * from mosq2013 where newDate between '2013-02-01' and '2013-02-31'")
mar2013 <- sqldf("Select * from mosq2013 where newDate between '2013-03-01' and '2013-03-31'")
apr2013 <- sqldf("Select * from mosq2013 where newDate between '2013-04-01' and '2013-04-31'")
may2013 <- sqldf("Select * from mosq2013 where newDate between '2013-05-01' and '2013-05-31'")
jun2013 <- sqldf("Select * from mosq2013 where newDate between '2013-06-01' and '2013-06-31'")
jul2013 <- sqldf("Select * from mosq2013 where newDate between '2013-07-01' and '2013-07-31'")
aug2013 <- sqldf("Select * from mosq2013 where newDate between '2013-08-01' and '2013-08-31'")
sep2013 <- sqldf("Select * from mosq2013 where newDate between '2013-09-01' and '2013-09-31'")
oct2013 <- sqldf("Select * from mosq2013 where newDate between '2013-10-01' and '2013-10-31'")

# Select based on month from all 2014 data
jan2014 <- sqldf("Select * from mosq2014 where newDate between '2014-01-01' and '2014-01-31'")
feb2014 <- sqldf("Select * from mosq2014 where newDate between '2014-02-01' and '2014-02-31'")
mar2014 <- sqldf("Select * from mosq2014 where newDate between '2014-03-01' and '2014-03-31'")
apr2014 <- sqldf("Select * from mosq2014 where newDate between '2014-04-01' and '2014-04-31'")
may2014 <- sqldf("Select * from mosq2014 where newDate between '2014-05-01' and '2014-05-31'")
jun2014 <- sqldf("Select * from mosq2014 where newDate between '2014-06-01' and '2014-06-31'")
jul2014 <- sqldf("Select * from mosq2014 where newDate between '2014-07-01' and '2014-07-31'")
aug2014 <- sqldf("Select * from mosq2014 where newDate between '2014-08-01' and '2014-08-31'")
sep2014 <- sqldf("Select * from mosq2014 where newDate between '2014-09-01' and '2014-09-31'")
oct2014 <- sqldf("Select * from mosq2014 where newDate between '2014-10-01' and '2014-10-31'")


# Select all from treatment zone
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

# Select all from buffer zone
buff <- sqldf("Select * from kdrData where project_code is 'buffer'")
# Select based on 2013 month from buffer zone
jan2013b <- sqldf("Select * from buff where newDate between '2013-01-01' and '2013-01-31'")
feb2013b <- sqldf("Select * from buff where newDate between '2013-02-01' and '2013-02-31'")
mar2013b <- sqldf("Select * from buff where newDate between '2013-03-01' and '2013-03-31'")
apr2013b <- sqldf("Select * from buff where newDate between '2013-04-01' and '2013-04-31'")
may2013b <- sqldf("Select * from buff where newDate between '2013-05-01' and '2013-05-31'")
jun2013b <- sqldf("Select * from buff where newDate between '2013-06-01' and '2013-06-31'")
jul2013b <- sqldf("Select * from buff where newDate between '2013-07-01' and '2013-07-31'")
aug2013b <- sqldf("Select * from buff where newDate between '2013-08-01' and '2013-08-31'")
sep2013b <- sqldf("Select * from buff where newDate between '2013-09-01' and '2013-09-31'")
oct2013b <- sqldf("Select * from buff where newDate between '2013-10-01' and '2013-10-31'")
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



# Run functions across all years ------------------------------------------
# For 1016 locus at all years - requires function_mc.1016.R
m2000 = mc.1016(mosq2000)
m2001 = mc.1016(mosq2001)
m2002 = mc.1016(mosq2002)
m2003 = mc.1016(mosq2003)
m2004 = mc.1016(mosq2004)
m2005 = mc.1016(mosq2005)
m2006 = mc.1016(mosq2006)
m2007 = mc.1016(mosq2007)
m2008 = mc.1016(mosq2008)
m2009 = mc.1016(mosq2009)
m2010 = mc.1016(mosq2010)
m2011 = mc.1016(mosq2011)
m2012 = mc.1016(mosq2012)
m2013 = mc.1016(mosq2013)
m2014 = mc.1016(mosq2014)
m2015 = mc.1016(mosq2015)
m2016 = mc.1016(mosq2016)
m2017 = mc.1016(mosq2017)

# For 1534 locus at all years - requires function_mc.1534.R
s2000 = mc.1534(mosq2000)
s2001 = mc.1534(mosq2001)
s2002 = mc.1534(mosq2002)
s2003 = mc.1534(mosq2003)
s2004 = mc.1534(mosq2004)
s2005 = mc.1534(mosq2005)
s2006 = mc.1534(mosq2006)
s2007 = mc.1534(mosq2007)
s2008 = mc.1534(mosq2008)
s2009 = mc.1534(mosq2009)
s2010 = mc.1534(mosq2010)
s2011 = mc.1534(mosq2011)
s2012 = mc.1534(mosq2012)
s2013 = mc.1534(mosq2013)
s2014 = mc.1534(mosq2014)
s2015 = mc.1534(mosq2015)
s2016 = mc.1534(mosq2016)
s2017 = mc.1534(mosq2017)

# For 410 locus at all years - requires function_mc.410.R
n2000 = mc.410(mosq2000)
n2001 = mc.410(mosq2001)
n2002 = mc.410(mosq2002)
n2003 = mc.410(mosq2003)
n2004 = mc.410(mosq2004)
n2005 = mc.410(mosq2005)
n2006 = mc.410(mosq2006)
n2007 = mc.410(mosq2007)
n2008 = mc.410(mosq2008)
n2009 = mc.410(mosq2009)
n2010 = mc.410(mosq2010)
n2011 = mc.410(mosq2011)
n2012 = mc.410(mosq2012)
n2013 = mc.410(mosq2013)
n2014 = mc.410(mosq2014)
n2015 = mc.410(mosq2015)
n2016 = mc.410(mosq2016)
n2017 = mc.410(mosq2017)

# For haplotypes at all years - requires function_mc.haps.R
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

# For 1016 locus at all months in 2013
mJan13 <- mc.1016(jan2013)
mFeb13 <- mc.1016(feb2013)
mMar13 <- mc.1016(mar2013)
mApr13 <- mc.1016(apr2013)
mMay13 <- mc.1016(may2013)
mJun13 <- mc.1016(jun2013)
mJul13 <- mc.1016(jul2013)
mAug13 <- mc.1016(aug2013)
mSep13 <- mc.1016(sep2013)
mOct13 <- mc.1016(oct2013)

# For 1016 locus at all months in 2014
mJan <- mc.1016(jan2014)
mFeb <- mc.1016(feb2014)
mMar <- mc.1016(mar2014)
mApr <- mc.1016(apr2014)
mMay <- mc.1016(may2014)
mJun <- mc.1016(jun2014)
mJul <- mc.1016(jul2014)
mAug <- mc.1016(aug2014)
mSep <- mc.1016(sep2014)
mOct <- mc.1016(oct2014)

# For 1534 locus at all months in 2013
sJan13 <- mc.1534(jan2013)
sFeb13 <- mc.1534(feb2013)
sMar13 <- mc.1534(mar2013)
sApr13 <- mc.1534(apr2013)
sMay13 <- mc.1534(may2013)
sJun13 <- mc.1534(jun2013)
sJul13 <- mc.1534(jul2013)
sAug13 <- mc.1534(aug2013)
sSep13 <- mc.1534(sep2013)
sOct13 <- mc.1534(oct2013)

# For 1534 locus at all months in 2014
sJan <- mc.1534(jan2014)
sFeb <- mc.1534(feb2014)
sMar <- mc.1534(mar2014)
sApr <- mc.1534(apr2014)
sMay <- mc.1534(may2014)
sJun <- mc.1534(jun2014)
sJul <- mc.1534(jul2014)
sAug <- mc.1534(aug2014)
sSep <- mc.1534(sep2014)
sOct <- mc.1534(oct2014)

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

# For 1016 locus at all months in 2013 - buffer zone
bJan13 <- mc.1016(jan2013b)
bFeb13 <- mc.1016(feb2013b)
bMar13 <- mc.1016(mar2013b)
bApr13 <- mc.1016(apr2013b)
bMay13 <- mc.1016(may2013b)
bJun13 <- mc.1016(jun2013b)
bJul13 <- mc.1016(jul2013b)
bAug13 <- mc.1016(aug2013b)
bSep13 <- mc.1016(sep2013b)
bOct13 <- mc.1016(oct2013b)

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




# Create dataframes -------------------------------------------------------
# Create list of years included in dataframe
year <- c(2000:2017)
# Create list of months included in dataframe - as numeric
month <- c(1:10)
# Create list of months included in dataframe - as objects
#month <- c("January", "February", "March", "April", "May", "June", "July", "August", "September", "October")

# For 1016 locus
df1016 <- rbind(m2000, m2001, m2002, m2003, m2004, m2005, m2006
               , m2007, m2008, m2009, m2010, m2011, m2012, m2013
               , m2014, m2015, m2016, m2017)
# Add year ID to rows in df rename
mc.1016.yr <- cbind(year, df1016)

# For 1534 locus
df1534 <- rbind(s2000, s2001, s2002, s2003, s2004, s2005, s2006
               , s2007, s2008, s2009, s2010, s2011, s2012, s2013
               , s2014, s2015, s2016, s2017)
# Add year ID to rows in df rename
mc.1534.yr <- cbind(year, df1534)

# For 410 locus
df410 <- rbind(n2000, n2001, n2002, n2003, n2004, n2005, n2006
                , n2007, n2008, n2009, n2010, n2011, n2012, n2013
                , n2014, n2015, n2016, n2017)
# Add year ID to rows in df rename
mc.410.yr <- cbind(year, df410)

# For haplotypes
dfHaps <- rbind(h2000, h2001, h2002, h2003, h2004, h2005, h2006
                , h2007, h2008, h2009, h2010, h2011, h2012, h2013
                , h2014, h2015, h2016, h2017)
# Add year ID to rows in df rename
mc.haps.yr <- cbind(dfHaps, year)

# For 1016 locus at all months in 2013
df13.1016 <- rbind(mJan13, mFeb13, mMar13, mApr13, mMay13, mJun13, mJul13, mAug13, mSep13, mOct13)
# Add year ID to rows in df rename
mc.1016.mo13 <- cbind(month, df13.1016)
# For 1016 locus at all months in 2014
df14.1016 <- rbind(mJan, mFeb, mMar, mApr, mMay, mJun, mJul, mAug, mSep, mOct)
# Add year ID to rows in df rename
mc.1016.mo <- cbind(month, df14.1016)

# For 1534 locus at all months in 2013
df13.1534 <- rbind(sJan13, sFeb13, sMar13, sApr13, sMay13, sJun13, sJul13, sAug13, sSep13, sOct13)
# Add year ID to rows in df rename
mc.1534.mo13 <- cbind(month, df13.1534)
# For 1534 locus at all months in 2014
df14.1534 <- rbind(sJan, sFeb, sMar, sApr, sMay, sJun, sJul, sAug, sSep, sOct)
# Add year ID to rows in df rename
mc.1534.mo <- cbind(month, df14.1534)

# For 1016 locus at all months in 2013 - treatment zone
dftrt13 <- rbind(tJan13, tFeb13, tMar13, tApr13, tMay13, tJun13, tJul13, tAug13, tSep13, tOct13)
# Add year ID to rows in df rename
mc.1016.t13 <- cbind(month, dftrt13)
# For 1016 locus at all months in 2014 - treatment zone
dftrt <- rbind(tJan, tFeb, tMar, tApr, tMay, tJun, tJul, tAug, tSep, tOct)
# Add year ID to rows in df rename
mc.1016.t <- cbind(month, dftrt)

# For 1016 locus at all months in 2013 - buffer zone
dfbuff13 <- rbind(bJan13, bFeb13, bMar13, bApr13, bMay13, bJun13, bJul13, bAug13, bSep13, bOct13)
# Add year ID to rows in df rename
mc.1016.b13 <- cbind(month, dfbuff13)
# For 1016 locus at all months in 2014 - buffer zone
dfbuff <- rbind(bJan, bFeb, bMar, bApr, bMay, bJun, bJul, bAug, bSep, bOct)
# Add year ID to rows in df rename
mc.1016.b <- cbind(month, dfbuff)

# ### To view dataframes
# mc.1016.yr
# mc.1534.yr
# mc.410.yr
# mc.haps.yr
# mc.1016.mo13
# mc.1016.mo
# mc.1534.mo13
# mc.1534.mo
# mc.1016.t13
# mc.1016.t
# mc.1016.b13
# mc.1016.b


# # Save dataframes ---------------------------------------------------------
# # These are required for plots, selection coefficient, and other analyses
# write.csv(mosq2014, file = "/Users/jenbaltz/Dropbox/GouldLab/Project_Mosquito/Database/mosq2014.csv", row.names = F)
# write.csv(mosq2013, file = "/Users/jenbaltz/Dropbox/GouldLab/Project_Mosquito/Database/mosq2013.csv", row.names = F)
# write.csv(mc.1016.yr, file = "/Users/jenbaltz/Dropbox/GouldLab/Project_Mosquito/Database/mc.1016.yr_reduced.csv", row.names = F)
# write.csv(mc.1534.yr, file = "/Users/jenbaltz/Dropbox/GouldLab/Project_Mosquito/Database/mc.1534.yr_reduced.csv", row.names = F)
# write.csv(mc.410.yr, file = "/Users/jenbaltz/Dropbox/GouldLab/Project_Mosquito/Database/mc.410.yr_reduced.csv", row.names = F)
# write.csv(mc.haps.yr, file = "/Users/jenbaltz/Dropbox/GouldLab/Project_Mosquito/Database/mc.haps.yr_reduced.csv", row.names = F)
# write.csv(mc.1016.mo13, file = "/Users/jenbaltz/Dropbox/GouldLab/Project_Mosquito/Database/mc.1016.mo13_reduced.csv", row.names = F)
# write.csv(mc.1016.mo, file = "/Users/jenbaltz/Dropbox/GouldLab/Project_Mosquito/Database/mc.1016.mo_reduced.csv", row.names = F)
# write.csv(mc.1534.mo13, file = "/Users/jenbaltz/Dropbox/GouldLab/Project_Mosquito/Database/mc.1534.mo13_reduced.csv", row.names = F)
# write.csv(mc.1534.mo, file = "/Users/jenbaltz/Dropbox/GouldLab/Project_Mosquito/Database/mc.1534.mo_reduced.csv", row.names = F)
# write.csv(mc.1016.t13, file = "/Users/jenbaltz/Dropbox/GouldLab/Project_Mosquito/Database/mc.1016.t13_reduced.csv", row.names = F)
# write.csv(mc.1016.t, file = "/Users/jenbaltz/Dropbox/GouldLab/Project_Mosquito/Database/mc.1016.t_reduced.csv", row.names = F)
# write.csv(mc.1016.b13, file = "/Users/jenbaltz/Dropbox/GouldLab/Project_Mosquito/Database/mc.1016.b13_reduced.csv", row.names = F)
# write.csv(mc.1016.b, file = "/Users/jenbaltz/Dropbox/GouldLab/Project_Mosquito/Database/mc.1016.b_reduced.csv", row.names = F)


