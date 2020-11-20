### Melt Curve Analysis
### File Started: 31 May 2017

# Prepare working environment ---------------------------------------------
# set working directory
setwd("~/Dropbox/GouldLab/Project_Mosquito/Database")

# clear working environment
rm(list = ls())

# Load libraries
library(gdata)
options(gsubfn.engine = "R") # prevents R from stalling while loading sqldf library
library(sqldf)
library(dplyr) # required for allele frequency calculations
library(ggplot2) # required to produce plots



# Concatenate MeltCurve data files ----------------------------------------
# Source script
source("R_Scripts/IQTmosq/loop_MeltCurve_catFiles.R")
# If you recieve the following error -- 
# "Error in file(file, "rt") : cannot open the connection' 
# -- make sure all analysis files are available and named properly in MeltCurve folder"

# set working directory back to Database folder
setwd("~/Dropbox/GouldLab/Project_Mosquito/Database")



# Load data -------------------------------------------------------
# Isolation log contains all samples isolated
# Note: make sure that IsolationLog_IQT-Mosq.csv file has the date column in the MM-DD-YYYY format
kdrData <- read.csv("IsolationLog_IQT-Mosq.csv")

# masterdec.csv contains GPS and neighborhood info
masterdec <- read.csv("../QGIS_Files/Routput/masterdec.csv")

# These files contain genotype information
merged1016_rep1 <- read.csv("MeltCurve_1016_rep1.csv")
merged1016_rep2 <- read.csv("MeltCurve_1016_rep2.csv")
merged1534_rep1 <- read.csv("MeltCurve_1534_rep1.csv")
merged1534_rep2 <- read.csv("MeltCurve_1534_rep2.csv")
merged410_rep1 <- read.csv("MeltCurve_410_rep1.csv")
merged410_rep2 <- read.csv("MeltCurve_410_rep2.csv")

# exptZone contains zone information for 2013 & 2014 experiments
exptZone <- read.csv("Location_Zone.csv")



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
write.csv(noZone, "noZone.csv", row.names = FALSE)

# Those files with no Zone but with known location codes
noZoneNoUnknowns <- noZone[noZone$Location_Code != "unknown",]
write.csv(noZoneNoUnknowns, "noZoneNoUnknowns.csv", row.names = FALSE)

####
# Key Code: Remove unknown location codes from kdrData and save file
####
# "Unknown" location codes are typically from mosquitoes collected by Steve Stoddard outside of IQT
kdrData <- kdrData[kdrData$Location_Code != "unknown", ]
write.csv(kdrData,"~/Dropbox/GouldLab/Project_Mosquito/Database/kdrData_all_reduced.csv", row.names = F)


# Convert melt curve output to readable genotype and haplotype
# Load function to convert melt curve output
source("R_Scripts/IQTmosq/function_convertMergedMeltCurve.R")

# Run function for each column
kdrData$V1016I_converted <- convert1016(kdrData)
kdrData$F1534C_converted <- convert1534(kdrData)
kdrData$V410L_converted <- convert410(kdrData)

# Create haplotype for loci 1016 and 1534 - 410 is not included here
kdrData$haplotype <- paste0(kdrData$V1016I_converted, kdrData$F1534C_converted)

### Save file
# Write to file
write.csv(kdrData, paste0("~/Dropbox/GouldLab/Project_Mosquito/Database/kdrData_reduced_archived/kdrData_reduced_", Sys.Date(), ".csv"), row.names = F)
write.csv(kdrData, "~/Dropbox/GouldLab/Project_Mosquito/Database/kdrData_reduced.csv", row.names = F)



# Parse data for allele frequency analysis
# Select based on month from each year---------------------
# Select based on month from all 2000 data
jan2000 <- sqldf("Select * from kdrData where newDate between '2000-01-01' and '2000-01-31'")
feb2000 <- sqldf("Select * from kdrData where newDate between '2000-02-01' and '2000-02-31'")
mar2000 <- sqldf("Select * from kdrData where newDate between '2000-03-01' and '2000-03-31'")
apr2000 <- sqldf("Select * from kdrData where newDate between '2000-04-01' and '2000-04-31'")
may2000 <- sqldf("Select * from kdrData where newDate between '2000-05-01' and '2000-05-31'")
jun2000 <- sqldf("Select * from kdrData where newDate between '2000-06-01' and '2000-06-31'")
jul2000 <- sqldf("Select * from kdrData where newDate between '2000-07-01' and '2000-07-31'")
aug2000 <- sqldf("Select * from kdrData where newDate between '2000-08-01' and '2000-08-31'")
sep2000 <- sqldf("Select * from kdrData where newDate between '2000-09-01' and '2000-09-31'")
oct2000 <- sqldf("Select * from kdrData where newDate between '2000-10-01' and '2000-10-31'")
nov2000 <- sqldf("Select * from kdrData where newDate between '2000-11-01' and '2000-11-31'")
dec2000 <- sqldf("Select * from kdrData where newDate between '2000-12-01' and '2000-12-31'")

# Select based on month from all 2001 data
jan2001 <- sqldf("Select * from kdrData where newDate between '2001-01-01' and '2001-01-31'")
feb2001 <- sqldf("Select * from kdrData where newDate between '2001-02-01' and '2001-02-31'")
mar2001 <- sqldf("Select * from kdrData where newDate between '2001-03-01' and '2001-03-31'")
apr2001 <- sqldf("Select * from kdrData where newDate between '2001-04-01' and '2001-04-31'")
may2001 <- sqldf("Select * from kdrData where newDate between '2001-05-01' and '2001-05-31'")
jun2001 <- sqldf("Select * from kdrData where newDate between '2001-06-01' and '2001-06-31'")
jul2001 <- sqldf("Select * from kdrData where newDate between '2001-07-01' and '2001-07-31'")
aug2001 <- sqldf("Select * from kdrData where newDate between '2001-08-01' and '2001-08-31'")
sep2001 <- sqldf("Select * from kdrData where newDate between '2001-09-01' and '2001-09-31'")
oct2001 <- sqldf("Select * from kdrData where newDate between '2001-10-01' and '2001-10-31'")
nov2001 <- sqldf("Select * from kdrData where newDate between '2001-11-01' and '2001-11-31'")
dec2001 <- sqldf("Select * from kdrData where newDate between '2001-12-01' and '2001-12-31'")

# Select based on month from all 2002 data
jan2002 <- sqldf("Select * from kdrData where newDate between '2002-01-01' and '2002-01-31'")
feb2002 <- sqldf("Select * from kdrData where newDate between '2002-02-01' and '2002-02-31'")
mar2002 <- sqldf("Select * from kdrData where newDate between '2002-03-01' and '2002-03-31'")
apr2002 <- sqldf("Select * from kdrData where newDate between '2002-04-01' and '2002-04-31'")
may2002 <- sqldf("Select * from kdrData where newDate between '2002-05-01' and '2002-05-31'")
jun2002 <- sqldf("Select * from kdrData where newDate between '2002-06-01' and '2002-06-31'")
jul2002 <- sqldf("Select * from kdrData where newDate between '2002-07-01' and '2002-07-31'")
aug2002 <- sqldf("Select * from kdrData where newDate between '2002-08-01' and '2002-08-31'")
sep2002 <- sqldf("Select * from kdrData where newDate between '2002-09-01' and '2002-09-31'")
oct2002 <- sqldf("Select * from kdrData where newDate between '2002-10-01' and '2002-10-31'")
nov2002 <- sqldf("Select * from kdrData where newDate between '2002-11-01' and '2002-11-31'")
dec2002 <- sqldf("Select * from kdrData where newDate between '2002-12-01' and '2002-12-31'")

# Select based on month from all 2003 data
jan2003 <- sqldf("Select * from kdrData where newDate between '2003-01-01' and '2003-01-31'")
feb2003 <- sqldf("Select * from kdrData where newDate between '2003-02-01' and '2003-02-31'")
mar2003 <- sqldf("Select * from kdrData where newDate between '2003-03-01' and '2003-03-31'")
apr2003 <- sqldf("Select * from kdrData where newDate between '2003-04-01' and '2003-04-31'")
may2003 <- sqldf("Select * from kdrData where newDate between '2003-05-01' and '2003-05-31'")
jun2003 <- sqldf("Select * from kdrData where newDate between '2003-06-01' and '2003-06-31'")
jul2003 <- sqldf("Select * from kdrData where newDate between '2003-07-01' and '2003-07-31'")
aug2003 <- sqldf("Select * from kdrData where newDate between '2003-08-01' and '2003-08-31'")
sep2003 <- sqldf("Select * from kdrData where newDate between '2003-09-01' and '2003-09-31'")
oct2003 <- sqldf("Select * from kdrData where newDate between '2003-10-01' and '2003-10-31'")
nov2003 <- sqldf("Select * from kdrData where newDate between '2003-11-01' and '2003-11-31'")
dec2003 <- sqldf("Select * from kdrData where newDate between '2003-12-01' and '2003-12-31'")

# Select based on month from all 2004 data
jan2004 <- sqldf("Select * from kdrData where newDate between '2004-01-01' and '2004-01-31'")
feb2004 <- sqldf("Select * from kdrData where newDate between '2004-02-01' and '2004-02-31'")
mar2004 <- sqldf("Select * from kdrData where newDate between '2004-03-01' and '2004-03-31'")
apr2004 <- sqldf("Select * from kdrData where newDate between '2004-04-01' and '2004-04-31'")
may2004 <- sqldf("Select * from kdrData where newDate between '2004-05-01' and '2004-05-31'")
jun2004 <- sqldf("Select * from kdrData where newDate between '2004-06-01' and '2004-06-31'")
jul2004 <- sqldf("Select * from kdrData where newDate between '2004-07-01' and '2004-07-31'")
aug2004 <- sqldf("Select * from kdrData where newDate between '2004-08-01' and '2004-08-31'")
sep2004 <- sqldf("Select * from kdrData where newDate between '2004-09-01' and '2004-09-31'")
oct2004 <- sqldf("Select * from kdrData where newDate between '2004-10-01' and '2004-10-31'")
nov2004 <- sqldf("Select * from kdrData where newDate between '2004-11-01' and '2004-11-31'")
dec2004 <- sqldf("Select * from kdrData where newDate between '2004-12-01' and '2004-12-31'")

# Select based on month from all 2005 data
jan2005 <- sqldf("Select * from kdrData where newDate between '2005-01-01' and '2005-01-31'")
feb2005 <- sqldf("Select * from kdrData where newDate between '2005-02-01' and '2005-02-31'")
mar2005 <- sqldf("Select * from kdrData where newDate between '2005-03-01' and '2005-03-31'")
apr2005 <- sqldf("Select * from kdrData where newDate between '2005-04-01' and '2005-04-31'")
may2005 <- sqldf("Select * from kdrData where newDate between '2005-05-01' and '2005-05-31'")
jun2005 <- sqldf("Select * from kdrData where newDate between '2005-06-01' and '2005-06-31'")
jul2005 <- sqldf("Select * from kdrData where newDate between '2005-07-01' and '2005-07-31'")
aug2005 <- sqldf("Select * from kdrData where newDate between '2005-08-01' and '2005-08-31'")
sep2005 <- sqldf("Select * from kdrData where newDate between '2005-09-01' and '2005-09-31'")
oct2005 <- sqldf("Select * from kdrData where newDate between '2005-10-01' and '2005-10-31'")
nov2005 <- sqldf("Select * from kdrData where newDate between '2005-11-01' and '2005-11-31'")
dec2005 <- sqldf("Select * from kdrData where newDate between '2005-12-01' and '2005-12-31'")

# Select based on month from all 2006 data
jan2006 <- sqldf("Select * from kdrData where newDate between '2006-01-01' and '2006-01-31'")
feb2006 <- sqldf("Select * from kdrData where newDate between '2006-02-01' and '2006-02-31'")
mar2006 <- sqldf("Select * from kdrData where newDate between '2006-03-01' and '2006-03-31'")
apr2006 <- sqldf("Select * from kdrData where newDate between '2006-04-01' and '2006-04-31'")
may2006 <- sqldf("Select * from kdrData where newDate between '2006-05-01' and '2006-05-31'")
jun2006 <- sqldf("Select * from kdrData where newDate between '2006-06-01' and '2006-06-31'")
jul2006 <- sqldf("Select * from kdrData where newDate between '2006-07-01' and '2006-07-31'")
aug2006 <- sqldf("Select * from kdrData where newDate between '2006-08-01' and '2006-08-31'")
sep2006 <- sqldf("Select * from kdrData where newDate between '2006-09-01' and '2006-09-31'")
oct2006 <- sqldf("Select * from kdrData where newDate between '2006-10-01' and '2006-10-31'")
nov2006 <- sqldf("Select * from kdrData where newDate between '2006-11-01' and '2006-11-31'")
dec2006 <- sqldf("Select * from kdrData where newDate between '2006-12-01' and '2006-12-31'")

# Select based on month from all 2007 data
jan2007 <- sqldf("Select * from kdrData where newDate between '2007-01-01' and '2007-01-31'")
feb2007 <- sqldf("Select * from kdrData where newDate between '2007-02-01' and '2007-02-31'")
mar2007 <- sqldf("Select * from kdrData where newDate between '2007-03-01' and '2007-03-31'")
apr2007 <- sqldf("Select * from kdrData where newDate between '2007-04-01' and '2007-04-31'")
may2007 <- sqldf("Select * from kdrData where newDate between '2007-05-01' and '2007-05-31'")
jun2007 <- sqldf("Select * from kdrData where newDate between '2007-06-01' and '2007-06-31'")
jul2007 <- sqldf("Select * from kdrData where newDate between '2007-07-01' and '2007-07-31'")
aug2007 <- sqldf("Select * from kdrData where newDate between '2007-08-01' and '2007-08-31'")
sep2007 <- sqldf("Select * from kdrData where newDate between '2007-09-01' and '2007-09-31'")
oct2007 <- sqldf("Select * from kdrData where newDate between '2007-10-01' and '2007-10-31'")
nov2007 <- sqldf("Select * from kdrData where newDate between '2007-11-01' and '2007-11-31'")
dec2007 <- sqldf("Select * from kdrData where newDate between '2007-12-01' and '2007-12-31'")

# Select based on month from all 2008 data
jan2008 <- sqldf("Select * from kdrData where newDate between '2008-01-01' and '2008-01-31'")
feb2008 <- sqldf("Select * from kdrData where newDate between '2008-02-01' and '2008-02-31'")
mar2008 <- sqldf("Select * from kdrData where newDate between '2008-03-01' and '2008-03-31'")
apr2008 <- sqldf("Select * from kdrData where newDate between '2008-04-01' and '2008-04-31'")
may2008 <- sqldf("Select * from kdrData where newDate between '2008-05-01' and '2008-05-31'")
jun2008 <- sqldf("Select * from kdrData where newDate between '2008-06-01' and '2008-06-31'")
jul2008 <- sqldf("Select * from kdrData where newDate between '2008-07-01' and '2008-07-31'")
aug2008 <- sqldf("Select * from kdrData where newDate between '2008-08-01' and '2008-08-31'")
sep2008 <- sqldf("Select * from kdrData where newDate between '2008-09-01' and '2008-09-31'")
oct2008 <- sqldf("Select * from kdrData where newDate between '2008-10-01' and '2008-10-31'")
nov2008 <- sqldf("Select * from kdrData where newDate between '2008-11-01' and '2008-11-31'")
dec2008 <- sqldf("Select * from kdrData where newDate between '2008-12-01' and '2008-12-31'")

# Select based on month from all 2009 data
jan2009 <- sqldf("Select * from kdrData where newDate between '2009-01-01' and '2009-01-31'")
feb2009 <- sqldf("Select * from kdrData where newDate between '2009-02-01' and '2009-02-31'")
mar2009 <- sqldf("Select * from kdrData where newDate between '2009-03-01' and '2009-03-31'")
apr2009 <- sqldf("Select * from kdrData where newDate between '2009-04-01' and '2009-04-31'")
may2009 <- sqldf("Select * from kdrData where newDate between '2009-05-01' and '2009-05-31'")
jun2009 <- sqldf("Select * from kdrData where newDate between '2009-06-01' and '2009-06-31'")
jul2009 <- sqldf("Select * from kdrData where newDate between '2009-07-01' and '2009-07-31'")
aug2009 <- sqldf("Select * from kdrData where newDate between '2009-08-01' and '2009-08-31'")
sep2009 <- sqldf("Select * from kdrData where newDate between '2009-09-01' and '2009-09-31'")
oct2009 <- sqldf("Select * from kdrData where newDate between '2009-10-01' and '2009-10-31'")
nov2009 <- sqldf("Select * from kdrData where newDate between '2009-11-01' and '2009-11-31'")
dec2009 <- sqldf("Select * from kdrData where newDate between '2009-12-01' and '2009-12-31'")

# Select based on month from all 2010 data
jan2010 <- sqldf("Select * from kdrData where newDate between '2010-01-01' and '2010-01-31'")
feb2010 <- sqldf("Select * from kdrData where newDate between '2010-02-01' and '2010-02-31'")
mar2010 <- sqldf("Select * from kdrData where newDate between '2010-03-01' and '2010-03-31'")
apr2010 <- sqldf("Select * from kdrData where newDate between '2010-04-01' and '2010-04-31'")
may2010 <- sqldf("Select * from kdrData where newDate between '2010-05-01' and '2010-05-31'")
jun2010 <- sqldf("Select * from kdrData where newDate between '2010-06-01' and '2010-06-31'")
jul2010 <- sqldf("Select * from kdrData where newDate between '2010-07-01' and '2010-07-31'")
aug2010 <- sqldf("Select * from kdrData where newDate between '2010-08-01' and '2010-08-31'")
sep2010 <- sqldf("Select * from kdrData where newDate between '2010-09-01' and '2010-09-31'")
oct2010 <- sqldf("Select * from kdrData where newDate between '2010-10-01' and '2010-10-31'")
nov2010 <- sqldf("Select * from kdrData where newDate between '2010-11-01' and '2010-11-31'")
dec2010 <- sqldf("Select * from kdrData where newDate between '2010-12-01' and '2010-12-31'")

# Select based on month from all 2011 data
jan2011 <- sqldf("Select * from kdrData where newDate between '2011-01-01' and '2011-01-31'")
feb2011 <- sqldf("Select * from kdrData where newDate between '2011-02-01' and '2011-02-31'")
mar2011 <- sqldf("Select * from kdrData where newDate between '2011-03-01' and '2011-03-31'")
apr2011 <- sqldf("Select * from kdrData where newDate between '2011-04-01' and '2011-04-31'")
may2011 <- sqldf("Select * from kdrData where newDate between '2011-05-01' and '2011-05-31'")
jun2011 <- sqldf("Select * from kdrData where newDate between '2011-06-01' and '2011-06-31'")
jul2011 <- sqldf("Select * from kdrData where newDate between '2011-07-01' and '2011-07-31'")
aug2011 <- sqldf("Select * from kdrData where newDate between '2011-08-01' and '2011-08-31'")
sep2011 <- sqldf("Select * from kdrData where newDate between '2011-09-01' and '2011-09-31'")
oct2011 <- sqldf("Select * from kdrData where newDate between '2011-10-01' and '2011-10-31'")
nov2011 <- sqldf("Select * from kdrData where newDate between '2011-11-01' and '2011-11-31'")
dec2011 <- sqldf("Select * from kdrData where newDate between '2011-12-01' and '2011-12-31'")

# Select based on month from all 2012 data
jan2012 <- sqldf("Select * from kdrData where newDate between '2012-01-01' and '2012-01-31'")
feb2012 <- sqldf("Select * from kdrData where newDate between '2012-02-01' and '2012-02-31'")
mar2012 <- sqldf("Select * from kdrData where newDate between '2012-03-01' and '2012-03-31'")
apr2012 <- sqldf("Select * from kdrData where newDate between '2012-04-01' and '2012-04-31'")
may2012 <- sqldf("Select * from kdrData where newDate between '2012-05-01' and '2012-05-31'")
jun2012 <- sqldf("Select * from kdrData where newDate between '2012-06-01' and '2012-06-31'")
jul2012 <- sqldf("Select * from kdrData where newDate between '2012-07-01' and '2012-07-31'")
aug2012 <- sqldf("Select * from kdrData where newDate between '2012-08-01' and '2012-08-31'")
sep2012 <- sqldf("Select * from kdrData where newDate between '2012-09-01' and '2012-09-31'")
oct2012 <- sqldf("Select * from kdrData where newDate between '2012-10-01' and '2012-10-31'")
nov2012 <- sqldf("Select * from kdrData where newDate between '2012-11-01' and '2012-11-31'")
dec2012 <- sqldf("Select * from kdrData where newDate between '2012-12-01' and '2012-12-31'")

# Select based on month from all 2013 data
jan2013 <- sqldf("Select * from kdrData where newDate between '2013-01-01' and '2013-01-31'")
feb2013 <- sqldf("Select * from kdrData where newDate between '2013-02-01' and '2013-02-31'")
mar2013 <- sqldf("Select * from kdrData where newDate between '2013-03-01' and '2013-03-31'")
apr2013 <- sqldf("Select * from kdrData where newDate between '2013-04-01' and '2013-04-31'")
may2013 <- sqldf("Select * from kdrData where newDate between '2013-05-01' and '2013-05-31'")
jun2013 <- sqldf("Select * from kdrData where newDate between '2013-06-01' and '2013-06-31'")
jul2013 <- sqldf("Select * from kdrData where newDate between '2013-07-01' and '2013-07-31'")
aug2013 <- sqldf("Select * from kdrData where newDate between '2013-08-01' and '2013-08-31'")
sep2013 <- sqldf("Select * from kdrData where newDate between '2013-09-01' and '2013-09-31'")
oct2013 <- sqldf("Select * from kdrData where newDate between '2013-10-01' and '2013-10-31'")
nov2013 <- sqldf("Select * from kdrData where newDate between '2013-11-01' and '2013-11-31'")
dec2013 <- sqldf("Select * from kdrData where newDate between '2013-12-01' and '2013-12-31'")

# Select based on month from all 2014 data
jan2014 <- sqldf("Select * from kdrData where newDate between '2014-01-01' and '2014-01-31'")
feb2014 <- sqldf("Select * from kdrData where newDate between '2014-02-01' and '2014-02-31'")
mar2014 <- sqldf("Select * from kdrData where newDate between '2014-03-01' and '2014-03-31'")
apr2014 <- sqldf("Select * from kdrData where newDate between '2014-04-01' and '2014-04-31'")
may2014 <- sqldf("Select * from kdrData where newDate between '2014-05-01' and '2014-05-31'")
jun2014 <- sqldf("Select * from kdrData where newDate between '2014-06-01' and '2014-06-31'")
jul2014 <- sqldf("Select * from kdrData where newDate between '2014-07-01' and '2014-07-31'")
aug2014 <- sqldf("Select * from kdrData where newDate between '2014-08-01' and '2014-08-31'")
sep2014 <- sqldf("Select * from kdrData where newDate between '2014-09-01' and '2014-09-31'")
oct2014 <- sqldf("Select * from kdrData where newDate between '2014-10-01' and '2014-10-31'")
nov2014 <- sqldf("Select * from kdrData where newDate between '2014-11-01' and '2014-11-31'")
dec2014 <- sqldf("Select * from kdrData where newDate between '2014-12-01' and '2014-12-31'")

# Select based on month from all 2015 data
jan2015 <- sqldf("Select * from kdrData where newDate between '2015-01-01' and '2015-01-31'")
feb2015 <- sqldf("Select * from kdrData where newDate between '2015-02-01' and '2015-02-31'")
mar2015 <- sqldf("Select * from kdrData where newDate between '2015-03-01' and '2015-03-31'")
apr2015 <- sqldf("Select * from kdrData where newDate between '2015-04-01' and '2015-04-31'")
may2015 <- sqldf("Select * from kdrData where newDate between '2015-05-01' and '2015-05-31'")
jun2015 <- sqldf("Select * from kdrData where newDate between '2015-06-01' and '2015-06-31'")
jul2015 <- sqldf("Select * from kdrData where newDate between '2015-07-01' and '2015-07-31'")
aug2015 <- sqldf("Select * from kdrData where newDate between '2015-08-01' and '2015-08-31'")
sep2015 <- sqldf("Select * from kdrData where newDate between '2015-09-01' and '2015-09-31'")
oct2015 <- sqldf("Select * from kdrData where newDate between '2015-10-01' and '2015-10-31'")
nov2015 <- sqldf("Select * from kdrData where newDate between '2015-11-01' and '2015-11-31'")
dec2015 <- sqldf("Select * from kdrData where newDate between '2015-12-01' and '2015-12-31'")

# Select based on month from all 2016 data
jan2016 <- sqldf("Select * from kdrData where newDate between '2016-01-01' and '2016-01-31'")
feb2016 <- sqldf("Select * from kdrData where newDate between '2016-02-01' and '2016-02-31'")
mar2016 <- sqldf("Select * from kdrData where newDate between '2016-03-01' and '2016-03-31'")
apr2016 <- sqldf("Select * from kdrData where newDate between '2016-04-01' and '2016-04-31'")
may2016 <- sqldf("Select * from kdrData where newDate between '2016-05-01' and '2016-05-31'")
jun2016 <- sqldf("Select * from kdrData where newDate between '2016-06-01' and '2016-06-31'")
jul2016 <- sqldf("Select * from kdrData where newDate between '2016-07-01' and '2016-07-31'")
aug2016 <- sqldf("Select * from kdrData where newDate between '2016-08-01' and '2016-08-31'")
sep2016 <- sqldf("Select * from kdrData where newDate between '2016-09-01' and '2016-09-31'")
oct2016 <- sqldf("Select * from kdrData where newDate between '2016-10-01' and '2016-10-31'")
nov2016 <- sqldf("Select * from kdrData where newDate between '2016-11-01' and '2016-11-31'")
dec2016 <- sqldf("Select * from kdrData where newDate between '2016-12-01' and '2016-12-31'")

# Select based on month from all 2017 data
jan2017 <- sqldf("Select * from kdrData where newDate between '2017-01-01' and '2017-01-31'")
feb2017 <- sqldf("Select * from kdrData where newDate between '2017-02-01' and '2017-02-31'")
mar2017 <- sqldf("Select * from kdrData where newDate between '2017-03-01' and '2017-03-31'")
apr2017 <- sqldf("Select * from kdrData where newDate between '2017-04-01' and '2017-04-31'")
may2017 <- sqldf("Select * from kdrData where newDate between '2017-05-01' and '2017-05-31'")
jun2017 <- sqldf("Select * from kdrData where newDate between '2017-06-01' and '2017-06-31'")
jul2017 <- sqldf("Select * from kdrData where newDate between '2017-07-01' and '2017-07-31'")
aug2017 <- sqldf("Select * from kdrData where newDate between '2017-08-01' and '2017-08-31'")
sep2017 <- sqldf("Select * from kdrData where newDate between '2017-09-01' and '2017-09-31'")
oct2017 <- sqldf("Select * from kdrData where newDate between '2017-10-01' and '2017-10-31'")
nov2017 <- sqldf("Select * from kdrData where newDate between '2017-11-01' and '2017-11-31'")
dec2017 <- sqldf("Select * from kdrData where newDate between '2017-12-01' and '2017-12-31'")

# Source functions to create dataframe of genotype counts, allele  --------
# Source functions --------------
source("R_Scripts/IQTmosq/function_mc.1016.R")
source("R_Scripts/IQTmosq/function_mc.1534.R")
source("R_Scripts/IQTmosq/function_mc.410.R")
source("R_Scripts/IQTmosq/function_mc.haps.R")

# Run mc.1016() for all month+year combos-------------
# 2000 
mJan00 <- mc.1016(jan2000)
mFeb00 <- mc.1016(feb2000)
mMar00 <- mc.1016(mar2000)
mApr00 <- mc.1016(apr2000)
mMay00 <- mc.1016(may2000)
mJun00 <- mc.1016(jun2000)
mJul00 <- mc.1016(jul2000)
mAug00 <- mc.1016(aug2000)
mSep00 <- mc.1016(sep2000)
mOct00 <- mc.1016(oct2000)
mNov00 <- mc.1016(nov2000)
mDec00 <- mc.1016(dec2000)

# 2001 
mJan01 <- mc.1016(jan2001)
mFeb01 <- mc.1016(feb2001)
mMar01 <- mc.1016(mar2001)
mApr01 <- mc.1016(apr2001)
mMay01 <- mc.1016(may2001)
mJun01 <- mc.1016(jun2001)
mJul01 <- mc.1016(jul2001)
mAug01 <- mc.1016(aug2001)
mSep01 <- mc.1016(sep2001)
mOct01 <- mc.1016(oct2001)
mNov01 <- mc.1016(nov2001)
mDec01 <- mc.1016(dec2001)

# 2002
mJan02 <- mc.1016(jan2002)
mFeb02 <- mc.1016(feb2002)
mMar02 <- mc.1016(mar2002)
mApr02 <- mc.1016(apr2002)
mMay02 <- mc.1016(may2002)
mJun02 <- mc.1016(jun2002)
mJul02 <- mc.1016(jul2002)
mAug02 <- mc.1016(aug2002)
mSep02 <- mc.1016(sep2002)
mOct02 <- mc.1016(oct2002)
mNov02 <- mc.1016(nov2002)
mDec02 <- mc.1016(dec2002)

# 2003 
mJan03 <- mc.1016(jan2003)
mFeb03 <- mc.1016(feb2003)
mMar03 <- mc.1016(mar2003)
mApr03 <- mc.1016(apr2003)
mMay03 <- mc.1016(may2003)
mJun03 <- mc.1016(jun2003)
mJul03 <- mc.1016(jul2003)
mAug03 <- mc.1016(aug2003)
mSep03 <- mc.1016(sep2003)
mOct03 <- mc.1016(oct2003)
mNov03 <- mc.1016(nov2003)
mDec03 <- mc.1016(dec2003)

# 2004 
mJan04 <- mc.1016(jan2004)
mFeb04 <- mc.1016(feb2004)
mMar04 <- mc.1016(mar2004)
mApr04 <- mc.1016(apr2004)
mMay04 <- mc.1016(may2004)
mJun04 <- mc.1016(jun2004)
mJul04 <- mc.1016(jul2004)
mAug04 <- mc.1016(aug2004)
mSep04 <- mc.1016(sep2004)
mOct04 <- mc.1016(oct2004)
mNov04 <- mc.1016(nov2004)
mDec04 <- mc.1016(dec2004)

# 2005 
mJan05 <- mc.1016(jan2005)
mFeb05 <- mc.1016(feb2005)
mMar05 <- mc.1016(mar2005)
mApr05 <- mc.1016(apr2005)
mMay05 <- mc.1016(may2005)
mJun05 <- mc.1016(jun2005)
mJul05 <- mc.1016(jul2005)
mAug05 <- mc.1016(aug2005)
mSep05 <- mc.1016(sep2005)
mOct05 <- mc.1016(oct2005)
mNov05 <- mc.1016(nov2005)
mDec05 <- mc.1016(dec2005)

# 2006 
mJan06 <- mc.1016(jan2006)
mFeb06 <- mc.1016(feb2006)
mMar06 <- mc.1016(mar2006)
mApr06 <- mc.1016(apr2006)
mMay06 <- mc.1016(may2006)
mJun06 <- mc.1016(jun2006)
mJul06 <- mc.1016(jul2006)
mAug06 <- mc.1016(aug2006)
mSep06 <- mc.1016(sep2006)
mOct06 <- mc.1016(oct2006)
mNov06 <- mc.1016(nov2006)
mDec06 <- mc.1016(dec2006)

# 2007 
mJan07 <- mc.1016(jan2007)
mFeb07 <- mc.1016(feb2007)
mMar07 <- mc.1016(mar2007)
mApr07 <- mc.1016(apr2007)
mMay07 <- mc.1016(may2007)
mJun07 <- mc.1016(jun2007)
mJul07 <- mc.1016(jul2007)
mAug07 <- mc.1016(aug2007)
mSep07 <- mc.1016(sep2007)
mOct07 <- mc.1016(oct2007)
mNov07 <- mc.1016(nov2007)
mDec07 <- mc.1016(dec2007)

# 2008 
mJan08 <- mc.1016(jan2008)
mFeb08 <- mc.1016(feb2008)
mMar08 <- mc.1016(mar2008)
mApr08 <- mc.1016(apr2008)
mMay08 <- mc.1016(may2008)
mJun08 <- mc.1016(jun2008)
mJul08 <- mc.1016(jul2008)
mAug08 <- mc.1016(aug2008)
mSep08 <- mc.1016(sep2008)
mOct08 <- mc.1016(oct2008)
mNov08 <- mc.1016(nov2008)
mDec08 <- mc.1016(dec2008)

# 2009 
mJan09 <- mc.1016(jan2009)
mFeb09 <- mc.1016(feb2009)
mMar09 <- mc.1016(mar2009)
mApr09 <- mc.1016(apr2009)
mMay09 <- mc.1016(may2009)
mJun09 <- mc.1016(jun2009)
mJul09 <- mc.1016(jul2009)
mAug09 <- mc.1016(aug2009)
mSep09 <- mc.1016(sep2009)
mOct09 <- mc.1016(oct2009)
mNov09 <- mc.1016(nov2009)
mDec09 <- mc.1016(dec2009)

# 2010 
mJan10 <- mc.1016(jan2010)
mFeb10 <- mc.1016(feb2010)
mMar10 <- mc.1016(mar2010)
mApr10 <- mc.1016(apr2010)
mMay10 <- mc.1016(may2010)
mJun10 <- mc.1016(jun2010)
mJul10 <- mc.1016(jul2010)
mAug10 <- mc.1016(aug2010)
mSep10 <- mc.1016(sep2010)
mOct10 <- mc.1016(oct2010)
mNov10 <- mc.1016(nov2010)
mDec10 <- mc.1016(dec2010)

# 2011 
mJan11 <- mc.1016(jan2011)
mFeb11 <- mc.1016(feb2011)
mMar11 <- mc.1016(mar2011)
mApr11 <- mc.1016(apr2011)
mMay11 <- mc.1016(may2011)
mJun11 <- mc.1016(jun2011)
mJul11 <- mc.1016(jul2011)
mAug11 <- mc.1016(aug2011)
mSep11 <- mc.1016(sep2011)
mOct11 <- mc.1016(oct2011)
mNov11 <- mc.1016(nov2011)
mDec11 <- mc.1016(dec2011)

# 2012 
mJan12 <- mc.1016(jan2012)
mFeb12 <- mc.1016(feb2012)
mMar12 <- mc.1016(mar2012)
mApr12 <- mc.1016(apr2012)
mMay12 <- mc.1016(may2012)
mJun12 <- mc.1016(jun2012)
mJul12 <- mc.1016(jul2012)
mAug12 <- mc.1016(aug2012)
mSep12 <- mc.1016(sep2012)
mOct12 <- mc.1016(oct2012)
mNov12 <- mc.1016(nov2012)
mDec12 <- mc.1016(dec2012)

# 2013 
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
mNov13 <- mc.1016(nov2013)
mDec13 <- mc.1016(dec2013)

# 2014 
mJan14 <- mc.1016(jan2014)
mFeb14 <- mc.1016(feb2014)
mMar14 <- mc.1016(mar2014)
mApr14 <- mc.1016(apr2014)
mMay14 <- mc.1016(may2014)
mJun14 <- mc.1016(jun2014)
mJul14 <- mc.1016(jul2014)
mAug14 <- mc.1016(aug2014)
mSep14 <- mc.1016(sep2014)
mOct14 <- mc.1016(oct2014)
mNov14 <- mc.1016(nov2014)
mDec14 <- mc.1016(dec2014)

# 2015 
mJan15 <- mc.1016(jan2015)
mFeb15 <- mc.1016(feb2015)
mMar15 <- mc.1016(mar2015)
mApr15 <- mc.1016(apr2015)
mMay15 <- mc.1016(may2015)
mJun15 <- mc.1016(jun2015)
mJul15 <- mc.1016(jul2015)
mAug15 <- mc.1016(aug2015)
mSep15 <- mc.1016(sep2015)
mOct15 <- mc.1016(oct2015)
mNov15 <- mc.1016(nov2015)
mDec15 <- mc.1016(dec2015)

# 2016 
mJan16 <- mc.1016(jan2016)
mFeb16 <- mc.1016(feb2016)
mMar16 <- mc.1016(mar2016)
mApr16 <- mc.1016(apr2016)
mMay16 <- mc.1016(may2016)
mJun16 <- mc.1016(jun2016)
mJul16 <- mc.1016(jul2016)
mAug16 <- mc.1016(aug2016)
mSep16 <- mc.1016(sep2016)
mOct16 <- mc.1016(oct2016)
mNov16 <- mc.1016(nov2016)
mDec16 <- mc.1016(dec2016)

# 2017 
mJan17 <- mc.1016(jan2017)
mFeb17 <- mc.1016(feb2017)
mMar17 <- mc.1016(mar2017)
mApr17 <- mc.1016(apr2017)
mMay17 <- mc.1016(may2017)
mJun17 <- mc.1016(jun2017)
mJul17 <- mc.1016(jul2017)
mAug17 <- mc.1016(aug2017)
mSep17 <- mc.1016(sep2017)
mOct17 <- mc.1016(oct2017)
mNov17 <- mc.1016(nov2017)
mDec17 <- mc.1016(dec2017)

# Run mc.1534() for all month+year combos-------------
# 2000 
sJan00 <- mc.1534(jan2000)
sFeb00 <- mc.1534(feb2000)
sMar00 <- mc.1534(mar2000)
sApr00 <- mc.1534(apr2000)
sMay00 <- mc.1534(may2000)
sJun00 <- mc.1534(jun2000)
sJul00 <- mc.1534(jul2000)
sAug00 <- mc.1534(aug2000)
sSep00 <- mc.1534(sep2000)
sOct00 <- mc.1534(oct2000)
sNov00 <- mc.1534(nov2000)
sDec00 <- mc.1534(dec2000)

# 2001 
sJan01 <- mc.1534(jan2001)
sFeb01 <- mc.1534(feb2001)
sMar01 <- mc.1534(mar2001)
sApr01 <- mc.1534(apr2001)
sMay01 <- mc.1534(may2001)
sJun01 <- mc.1534(jun2001)
sJul01 <- mc.1534(jul2001)
sAug01 <- mc.1534(aug2001)
sSep01 <- mc.1534(sep2001)
sOct01 <- mc.1534(oct2001)
sNov01 <- mc.1534(nov2001)
sDec01 <- mc.1534(dec2001)

# 2002
sJan02 <- mc.1534(jan2002)
sFeb02 <- mc.1534(feb2002)
sMar02 <- mc.1534(mar2002)
sApr02 <- mc.1534(apr2002)
sMay02 <- mc.1534(may2002)
sJun02 <- mc.1534(jun2002)
sJul02 <- mc.1534(jul2002)
sAug02 <- mc.1534(aug2002)
sSep02 <- mc.1534(sep2002)
sOct02 <- mc.1534(oct2002)
sNov02 <- mc.1534(nov2002)
sDec02 <- mc.1534(dec2002)

# 2003 
sJan03 <- mc.1534(jan2003)
sFeb03 <- mc.1534(feb2003)
sMar03 <- mc.1534(mar2003)
sApr03 <- mc.1534(apr2003)
sMay03 <- mc.1534(may2003)
sJun03 <- mc.1534(jun2003)
sJul03 <- mc.1534(jul2003)
sAug03 <- mc.1534(aug2003)
sSep03 <- mc.1534(sep2003)
sOct03 <- mc.1534(oct2003)
sNov03 <- mc.1534(nov2003)
sDec03 <- mc.1534(dec2003)

# 2004 
sJan04 <- mc.1534(jan2004)
sFeb04 <- mc.1534(feb2004)
sMar04 <- mc.1534(mar2004)
sApr04 <- mc.1534(apr2004)
sMay04 <- mc.1534(may2004)
sJun04 <- mc.1534(jun2004)
sJul04 <- mc.1534(jul2004)
sAug04 <- mc.1534(aug2004)
sSep04 <- mc.1534(sep2004)
sOct04 <- mc.1534(oct2004)
sNov04 <- mc.1534(nov2004)
sDec04 <- mc.1534(dec2004)

# 2005 
sJan05 <- mc.1534(jan2005)
sFeb05 <- mc.1534(feb2005)
sMar05 <- mc.1534(mar2005)
sApr05 <- mc.1534(apr2005)
sMay05 <- mc.1534(may2005)
sJun05 <- mc.1534(jun2005)
sJul05 <- mc.1534(jul2005)
sAug05 <- mc.1534(aug2005)
sSep05 <- mc.1534(sep2005)
sOct05 <- mc.1534(oct2005)
sNov05 <- mc.1534(nov2005)
sDec05 <- mc.1534(dec2005)

# 2006 
sJan06 <- mc.1534(jan2006)
sFeb06 <- mc.1534(feb2006)
sMar06 <- mc.1534(mar2006)
sApr06 <- mc.1534(apr2006)
sMay06 <- mc.1534(may2006)
sJun06 <- mc.1534(jun2006)
sJul06 <- mc.1534(jul2006)
sAug06 <- mc.1534(aug2006)
sSep06 <- mc.1534(sep2006)
sOct06 <- mc.1534(oct2006)
sNov06 <- mc.1534(nov2006)
sDec06 <- mc.1534(dec2006)

# 2007 
sJan07 <- mc.1534(jan2007)
sFeb07 <- mc.1534(feb2007)
sMar07 <- mc.1534(mar2007)
sApr07 <- mc.1534(apr2007)
sMay07 <- mc.1534(may2007)
sJun07 <- mc.1534(jun2007)
sJul07 <- mc.1534(jul2007)
sAug07 <- mc.1534(aug2007)
sSep07 <- mc.1534(sep2007)
sOct07 <- mc.1534(oct2007)
sNov07 <- mc.1534(nov2007)
sDec07 <- mc.1534(dec2007)

# 2008 
sJan08 <- mc.1534(jan2008)
sFeb08 <- mc.1534(feb2008)
sMar08 <- mc.1534(mar2008)
sApr08 <- mc.1534(apr2008)
sMay08 <- mc.1534(may2008)
sJun08 <- mc.1534(jun2008)
sJul08 <- mc.1534(jul2008)
sAug08 <- mc.1534(aug2008)
sSep08 <- mc.1534(sep2008)
sOct08 <- mc.1534(oct2008)
sNov08 <- mc.1534(nov2008)
sDec08 <- mc.1534(dec2008)

# 2009 
sJan09 <- mc.1534(jan2009)
sFeb09 <- mc.1534(feb2009)
sMar09 <- mc.1534(mar2009)
sApr09 <- mc.1534(apr2009)
sMay09 <- mc.1534(may2009)
sJun09 <- mc.1534(jun2009)
sJul09 <- mc.1534(jul2009)
sAug09 <- mc.1534(aug2009)
sSep09 <- mc.1534(sep2009)
sOct09 <- mc.1534(oct2009)
sNov09 <- mc.1534(nov2009)
sDec09 <- mc.1534(dec2009)

# 2010 
sJan10 <- mc.1534(jan2010)
sFeb10 <- mc.1534(feb2010)
sMar10 <- mc.1534(mar2010)
sApr10 <- mc.1534(apr2010)
sMay10 <- mc.1534(may2010)
sJun10 <- mc.1534(jun2010)
sJul10 <- mc.1534(jul2010)
sAug10 <- mc.1534(aug2010)
sSep10 <- mc.1534(sep2010)
sOct10 <- mc.1534(oct2010)
sNov10 <- mc.1534(nov2010)
sDec10 <- mc.1534(dec2010)

# 2011 
sJan11 <- mc.1534(jan2011)
sFeb11 <- mc.1534(feb2011)
sMar11 <- mc.1534(mar2011)
sApr11 <- mc.1534(apr2011)
sMay11 <- mc.1534(may2011)
sJun11 <- mc.1534(jun2011)
sJul11 <- mc.1534(jul2011)
sAug11 <- mc.1534(aug2011)
sSep11 <- mc.1534(sep2011)
sOct11 <- mc.1534(oct2011)
sNov11 <- mc.1534(nov2011)
sDec11 <- mc.1534(dec2011)

# 2012 
sJan12 <- mc.1534(jan2012)
sFeb12 <- mc.1534(feb2012)
sMar12 <- mc.1534(mar2012)
sApr12 <- mc.1534(apr2012)
sMay12 <- mc.1534(may2012)
sJun12 <- mc.1534(jun2012)
sJul12 <- mc.1534(jul2012)
sAug12 <- mc.1534(aug2012)
sSep12 <- mc.1534(sep2012)
sOct12 <- mc.1534(oct2012)
sNov12 <- mc.1534(nov2012)
sDec12 <- mc.1534(dec2012)

# 2013 
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
sNov13 <- mc.1534(nov2013)
sDec13 <- mc.1534(dec2013)

# 2014 
sJan14 <- mc.1534(jan2014)
sFeb14 <- mc.1534(feb2014)
sMar14 <- mc.1534(mar2014)
sApr14 <- mc.1534(apr2014)
sMay14 <- mc.1534(may2014)
sJun14 <- mc.1534(jun2014)
sJul14 <- mc.1534(jul2014)
sAug14 <- mc.1534(aug2014)
sSep14 <- mc.1534(sep2014)
sOct14 <- mc.1534(oct2014)
sNov14 <- mc.1534(nov2014)
sDec14 <- mc.1534(dec2014)

# 2015 
sJan15 <- mc.1534(jan2015)
sFeb15 <- mc.1534(feb2015)
sMar15 <- mc.1534(mar2015)
sApr15 <- mc.1534(apr2015)
sMay15 <- mc.1534(may2015)
sJun15 <- mc.1534(jun2015)
sJul15 <- mc.1534(jul2015)
sAug15 <- mc.1534(aug2015)
sSep15 <- mc.1534(sep2015)
sOct15 <- mc.1534(oct2015)
sNov15 <- mc.1534(nov2015)
sDec15 <- mc.1534(dec2015)

# 2016 
sJan16 <- mc.1534(jan2016)
sFeb16 <- mc.1534(feb2016)
sMar16 <- mc.1534(mar2016)
sApr16 <- mc.1534(apr2016)
sMay16 <- mc.1534(may2016)
sJun16 <- mc.1534(jun2016)
sJul16 <- mc.1534(jul2016)
sAug16 <- mc.1534(aug2016)
sSep16 <- mc.1534(sep2016)
sOct16 <- mc.1534(oct2016)
sNov16 <- mc.1534(nov2016)
sDec16 <- mc.1534(dec2016)

# 2017 
sJan17 <- mc.1534(jan2017)
sFeb17 <- mc.1534(feb2017)
sMar17 <- mc.1534(mar2017)
sApr17 <- mc.1534(apr2017)
sMay17 <- mc.1534(may2017)
sJun17 <- mc.1534(jun2017)
sJul17 <- mc.1534(jul2017)
sAug17 <- mc.1534(aug2017)
sSep17 <- mc.1534(sep2017)
sOct17 <- mc.1534(oct2017)
sNov17 <- mc.1534(nov2017)
sDec17 <- mc.1534(dec2017)

# Create dataframes -------------------------------------------------------
# Create list of years included in dataframe
year <- c(2000:2017)
# Create list of months included in dataframe - as numeric
month <- c(1:12)

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




# For 1016 locus
df1016 <- rbind(mJan00, mFeb00, mMar00, mApr00, mMay00, mJun00, mJul00, mAug00, mSep00, mOct00, mNov00, mDec00
                , mJan01, mFeb01, mMar01, mApr01, mMay01, mJun01, mJul01, mAug01, mSep01, mOct01, mNov01, mDec01
                , mJan02, mFeb02, mMar02, mApr02, mMay02, mJun02, mJul02, mAug02, mSep02, mOct02, mNov02, mDec02
                , mJan03, mFeb03, mMar03, mApr03, mMay03, mJun03, mJul03, mAug03, mSep03, mOct03, mNov03, mDec03
                , mJan04, mFeb04, mMar04, mApr04, mMay04, mJun04, mJul04, mAug04, mSep04, mOct04, mNov04, mDec04
                , mJan05, mFeb05, mMar05, mApr05, mMay05, mJun05, mJul05, mAug05, mSep05, mOct05, mNov05, mDec05
                , mJan06, mFeb06, mMar06, mApr06, mMay06, mJun06, mJul06, mAug06, mSep06, mOct06, mNov06, mDec06
                , mJan07, mFeb07, mMar07, mApr07, mMay07, mJun07, mJul07, mAug07, mSep07, mOct07, mNov07, mDec07
                , mJan08, mFeb08, mMar08, mApr08, mMay08, mJun08, mJul08, mAug08, mSep08, mOct08, mNov08, mDec08
                , mJan09, mFeb09, mMar09, mApr09, mMay09, mJun09, mJul09, mAug09, mSep09, mOct09, mNov09, mDec09
                , mJan10, mFeb10, mMar10, mApr10, mMay10, mJun10, mJul10, mAug10, mSep10, mOct10, mNov10, mDec10
                , mJan11, mFeb11, mMar11, mApr11, mMay11, mJun11, mJul11, mAug11, mSep11, mOct11, mNov11, mDec11
                , mJan12, mFeb12, mMar12, mApr12, mMay12, mJun12, mJul12, mAug12, mSep12, mOct12, mNov12, mDec12
                , mJan13, mFeb13, mMar13, mApr13, mMay13, mJun13, mJul13, mAug13, mSep13, mOct13, mNov13, mDec13
                , mJan14, mFeb14, mMar14, mApr14, mMay14, mJun14, mJul14, mAug14, mSep14, mOct14, mNov14, mDec14
                , mJan15, mFeb15, mMar15, mApr15, mMay15, mJun15, mJul15, mAug15, mSep15, mOct15, mNov15, mDec15
                , mJan16, mFeb16, mMar16, mApr16, mMay16, mJun16, mJul16, mAug16, mSep16, mOct16, mNov16, mDec16
                , mJan17, mFeb17, mMar17, mApr17, mMay17, mJun17, mJul17, mAug17, mSep17, mOct17, mNov17, mDec17)
# Add year ID to rows in df rename
mc.1016.moyr <- cbind(MonthYear, df1016)

# For 1534 locus
df1534 <- rbind(sJan00, sFeb00, sMar00, sApr00, sMay00, sJun00, sJul00, sAug00, sSep00, sOct00, sNov00, sDec00
                , sJan01, sFeb01, sMar01, sApr01, sMay01, sJun01, sJul01, sAug01, sSep01, sOct01, sNov01, sDec01
                , sJan02, sFeb02, sMar02, sApr02, sMay02, sJun02, sJul02, sAug02, sSep02, sOct02, sNov02, sDec02
                , sJan03, sFeb03, sMar03, sApr03, sMay03, sJun03, sJul03, sAug03, sSep03, sOct03, sNov03, sDec03
                , sJan04, sFeb04, sMar04, sApr04, sMay04, sJun04, sJul04, sAug04, sSep04, sOct04, sNov04, sDec04
                , sJan05, sFeb05, sMar05, sApr05, sMay05, sJun05, sJul05, sAug05, sSep05, sOct05, sNov05, sDec05
                , sJan06, sFeb06, sMar06, sApr06, sMay06, sJun06, sJul06, sAug06, sSep06, sOct06, sNov06, sDec06
                , sJan07, sFeb07, sMar07, sApr07, sMay07, sJun07, sJul07, sAug07, sSep07, sOct07, sNov07, sDec07
                , sJan08, sFeb08, sMar08, sApr08, sMay08, sJun08, sJul08, sAug08, sSep08, sOct08, sNov08, sDec08
                , sJan09, sFeb09, sMar09, sApr09, sMay09, sJun09, sJul09, sAug09, sSep09, sOct09, sNov09, sDec09
                , sJan10, sFeb10, sMar10, sApr10, sMay10, sJun10, sJul10, sAug10, sSep10, sOct10, sNov10, sDec10
                , sJan11, sFeb11, sMar11, sApr11, sMay11, sJun11, sJul11, sAug11, sSep11, sOct11, sNov11, sDec11
                , sJan12, sFeb12, sMar12, sApr12, sMay12, sJun12, sJul12, sAug12, sSep12, sOct12, sNov12, sDec12
                , sJan13, sFeb13, sMar13, sApr13, sMay13, sJun13, sJul13, sAug13, sSep13, sOct13, sNov13, sDec13
                , sJan14, sFeb14, sMar14, sApr14, sMay14, sJun14, sJul14, sAug14, sSep14, sOct14, sNov14, sDec14
                , sJan15, sFeb15, sMar15, sApr15, sMay15, sJun15, sJul15, sAug15, sSep15, sOct15, sNov15, sDec15
                , sJan16, sFeb16, sMar16, sApr16, sMay16, sJun16, sJul16, sAug16, sSep16, sOct16, sNov16, sDec16
                , sJan17, sFeb17, sMar17, sApr17, sMay17, sJun17, sJul17, sAug17, sSep17, sOct17, sNov17, sDec17)
# Add year ID to rows in df rename
mc.1534.moyr <- cbind(MonthYear, df1534)


# ### To view dataframes
# mc.1016.moyr
# mc.1534.moyr


# Save dataframes ---------------------------------------------------------
# These are required for plots, selection coefficient, and other analyses---------
write.csv(mc.1016.moyr, file = "/Users/jenbaltz/Dropbox/GouldLab/Project_Mosquito/Database/mc.1016.MoYr_reduced_byMonth.csv", row.names = F)
write.csv(mc.1534.moyr, file = "/Users/jenbaltz/Dropbox/GouldLab/Project_Mosquito/Database/mc.1534.MoYr_reduced_byMonth.csv", row.names = F)


# scratch work
library(lubridate)
head(parse_date_time(mc.1534.moyr$MonthYear, "m-y"))

plot(parse_date_time(mc.1534.moyr$MonthYear, "m-y"), mc.1534.moyr$freqR, type = "o")


# pull out samples from october 2002
temp <- kdrData[kdrData$newDate >= "2002-10-01" & kdrData$newDate <= "2002-10-31",]
temp
nrow(temp)
temp <- temp[, c(1:2, 11:13, 25:26, 28)]
nrow(temp)
temp <- temp[complete.cases(temp), ]
nrow(temp)
temp
