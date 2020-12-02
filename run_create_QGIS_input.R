###############################################################  
# This file is to create a subset of the masterdec that 
# contains all of the houses from where we have samples
# Started 30 Aug 2017
############################################################### 

##################### 
### clear working environment
rm(list = ls())

### set working directory
setwd("~/Dropbox/GouldLab/Project_Mosquito/Database")

##################### 
### Load datasets
# kdrData_reduced_[date].csv is the file that contains all samples with kdr genotypes and known location codes
kdr <- read.csv("kdrData_reduced.csv")
# lib5 platemap for NGS work
lib5 <- read.table("popmap.lib5.tsv", header = F, sep = "\t")


##################### 
### Prep dataframes for merging
# Create new column with only first 6 characters from Location_Code
kdr$LOC_CODE <- substr(kdr$Location_Code, 1, 6)

##################### 
### Merge dataframes
# Merge kdr with lib5 (aka - left join)
masterdec_ngs <- merge(lib5, kdr, by.x=c("V1"), by.y = c("mosquito_id"), all.x=TRUE)
# Subset on specific columns
masterdec_ngs <- masterdec_ngs[, c("V1", "V2", "Location_Code", "X", "Y", "newDate"
                                   , "F1534C_converted", "V1016I_converted", "haplotype")]
# Rename first two columns
colnames(masterdec_ngs)[1:2] <- c("mosquito_id", "population")
  
##################### 
### Prep dataframes for export to QGIS
# Remove rows with NA and errors in haplotype column
keep <- c("SSSS", "RRRR", "SSRR", "SSSR", "SRRR", "SRSS", "SRSR", "RRSS")
masterdec_kdr_all <- kdr[kdr$haplotype %in% keep,]

# Reduce dataset for easier viewing and use in QGIS
masterdec_kdr_reduced <- masterdec_kdr_all[c("LOC_CODE","NEIGHBORHO", "X","Y","mosquito_id","newDate","haplotype")]

##################### 
### Create .csv files for use in QGIS
write.csv(masterdec_kdr_reduced, "../QGIS_Files/Routput/masterdec_kdr.csv", row.names = F)
write.csv(masterdec_ngs, "../QGIS_Files/Routput/masterdec_ngs.csv", row.names = F)

##################### 
### Subsetting masterdec_kdr_reduced by year and saving for use in QGIS
# Subset masterdec_kdr_reduced$newDate on year
masterdec_2000 <- masterdec_kdr_reduced[format.Date(masterdec_kdr_reduced$newDate, "%Y")==2000,]
masterdec_2001 <- masterdec_kdr_reduced[format.Date(masterdec_kdr_reduced$newDate, "%Y")==2001,]
masterdec_2002 <- masterdec_kdr_reduced[format.Date(masterdec_kdr_reduced$newDate, "%Y")==2002,]
masterdec_2003 <- masterdec_kdr_reduced[format.Date(masterdec_kdr_reduced$newDate, "%Y")==2003,]
masterdec_2004 <- masterdec_kdr_reduced[format.Date(masterdec_kdr_reduced$newDate, "%Y")==2004,]
masterdec_2005 <- masterdec_kdr_reduced[format.Date(masterdec_kdr_reduced$newDate, "%Y")==2005,]
masterdec_2006 <- masterdec_kdr_reduced[format.Date(masterdec_kdr_reduced$newDate, "%Y")==2006,]
masterdec_2007 <- masterdec_kdr_reduced[format.Date(masterdec_kdr_reduced$newDate, "%Y")==2007,]
masterdec_2008 <- masterdec_kdr_reduced[format.Date(masterdec_kdr_reduced$newDate, "%Y")==2008,]
masterdec_2009 <- masterdec_kdr_reduced[format.Date(masterdec_kdr_reduced$newDate, "%Y")==2009,]
masterdec_2010 <- masterdec_kdr_reduced[format.Date(masterdec_kdr_reduced$newDate, "%Y")==2010,]
masterdec_2011 <- masterdec_kdr_reduced[format.Date(masterdec_kdr_reduced$newDate, "%Y")==2011,]
masterdec_2012 <- masterdec_kdr_reduced[format.Date(masterdec_kdr_reduced$newDate, "%Y")==2012,]
masterdec_2013 <- masterdec_kdr_reduced[format.Date(masterdec_kdr_reduced$newDate, "%Y")==2013,]
masterdec_2014 <- masterdec_kdr_reduced[format.Date(masterdec_kdr_reduced$newDate, "%Y")==2014,]
masterdec_2015 <- masterdec_kdr_reduced[format.Date(masterdec_kdr_reduced$newDate, "%Y")==2015,]
masterdec_2016 <- masterdec_kdr_reduced[format.Date(masterdec_kdr_reduced$newDate, "%Y")==2016,]
masterdec_2017 <- masterdec_kdr_reduced[format.Date(masterdec_kdr_reduced$newDate, "%Y")==2017,]

# Create files for each masterdec_[year] data frame
write.csv(masterdec_2000, "../QGIS_Files/Routput/masterdec_2000.csv", row.names = F)
write.csv(masterdec_2001, "../QGIS_Files/Routput/masterdec_2001.csv", row.names = F)
write.csv(masterdec_2002, "../QGIS_Files/Routput/masterdec_2002.csv", row.names = F)
write.csv(masterdec_2003, "../QGIS_Files/Routput/masterdec_2003.csv", row.names = F)
write.csv(masterdec_2004, "../QGIS_Files/Routput/masterdec_2004.csv", row.names = F)
write.csv(masterdec_2005, "../QGIS_Files/Routput/masterdec_2005.csv", row.names = F)
write.csv(masterdec_2006, "../QGIS_Files/Routput/masterdec_2006.csv", row.names = F)
write.csv(masterdec_2007, "../QGIS_Files/Routput/masterdec_2007.csv", row.names = F)
write.csv(masterdec_2008, "../QGIS_Files/Routput/masterdec_2008.csv", row.names = F)
write.csv(masterdec_2009, "../QGIS_Files/Routput/masterdec_2009.csv", row.names = F)
write.csv(masterdec_2010, "../QGIS_Files/Routput/masterdec_2010.csv", row.names = F)
write.csv(masterdec_2011, "../QGIS_Files/Routput/masterdec_2011.csv", row.names = F)
write.csv(masterdec_2012, "../QGIS_Files/Routput/masterdec_2012.csv", row.names = F)
write.csv(masterdec_2013, "../QGIS_Files/Routput/masterdec_2013.csv", row.names = F)
write.csv(masterdec_2014, "../QGIS_Files/Routput/masterdec_2014.csv", row.names = F)
write.csv(masterdec_2015, "../QGIS_Files/Routput/masterdec_2015.csv", row.names = F)
write.csv(masterdec_2016, "../QGIS_Files/Routput/masterdec_2016.csv", row.names = F)
write.csv(masterdec_2017, "../QGIS_Files/Routput/masterdec_2017.csv", row.names = F)



