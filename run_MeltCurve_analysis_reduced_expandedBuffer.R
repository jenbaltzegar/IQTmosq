### Melt Curve Analysis
### File Started: 31 May 2017

# Parse data for allele frequency analysis --------------------------------
### Parse data by year and/or month to use for allele frequency analysis
# Keep only those rows with project code

# # Select based on year
# mosq2000 <- sqldf("Select * from kdrData where newDate between '2000-01-01' and '2000-12-31'")
# mosq2001 <- sqldf("Select * from kdrData where newDate between '2001-01-01' and '2001-12-31'")
# mosq2002 <- sqldf("Select * from kdrData where newDate between '2002-01-01' and '2002-12-31'")
# mosq2003 <- sqldf("Select * from kdrData where newDate between '2003-01-01' and '2003-12-31'")
# mosq2004 <- sqldf("Select * from kdrData where newDate between '2004-01-01' and '2004-12-31'")
# mosq2005 <- sqldf("Select * from kdrData where newDate between '2005-01-01' and '2005-12-31'")
# mosq2006 <- sqldf("Select * from kdrData where newDate between '2006-01-01' and '2006-12-31'")
# mosq2007 <- sqldf("Select * from kdrData where newDate between '2007-01-01' and '2007-12-31'")
# mosq2008 <- sqldf("Select * from kdrData where newDate between '2008-01-01' and '2008-12-31'")
# mosq2009 <- sqldf("Select * from kdrData where newDate between '2009-01-01' and '2009-12-31'")
# mosq2010 <- sqldf("Select * from kdrData where newDate between '2010-01-01' and '2010-12-31'")
# mosq2011 <- sqldf("Select * from kdrData where newDate between '2011-01-01' and '2011-12-31'")
# mosq2012 <- sqldf("Select * from kdrData where newDate between '2012-01-01' and '2012-12-31'")
# mosq2013 <- sqldf("Select * from kdrData where newDate between '2013-01-01' and '2013-12-31'")
# mosq2014 <- sqldf("Select * from kdrData where newDate between '2014-01-01' and '2014-12-31'")
# mosq2015 <- sqldf("Select * from kdrData where newDate between '2015-01-01' and '2015-12-31'")
# mosq2016 <- sqldf("Select * from kdrData where newDate between '2016-01-01' and '2016-12-31'")
# mosq2017 <- sqldf("Select * from kdrData where newDate between '2017-01-01' and '2017-12-31'")
# 
# # Select based on month from all 2013 data
# jan2013 <- sqldf("Select * from mosq2013 where newDate between '2013-01-01' and '2013-01-31'")
# feb2013 <- sqldf("Select * from mosq2013 where newDate between '2013-02-01' and '2013-02-31'")
# mar2013 <- sqldf("Select * from mosq2013 where newDate between '2013-03-01' and '2013-03-31'")
# apr2013 <- sqldf("Select * from mosq2013 where newDate between '2013-04-01' and '2013-04-31'")
# may2013 <- sqldf("Select * from mosq2013 where newDate between '2013-05-01' and '2013-05-31'")
# jun2013 <- sqldf("Select * from mosq2013 where newDate between '2013-06-01' and '2013-06-31'")
# jul2013 <- sqldf("Select * from mosq2013 where newDate between '2013-07-01' and '2013-07-31'")
# aug2013 <- sqldf("Select * from mosq2013 where newDate between '2013-08-01' and '2013-08-31'")
# sep2013 <- sqldf("Select * from mosq2013 where newDate between '2013-09-01' and '2013-09-31'")
# oct2013 <- sqldf("Select * from mosq2013 where newDate between '2013-10-01' and '2013-10-31'")
# 
# # Select based on month from all 2014 data
# jan2014 <- sqldf("Select * from mosq2014 where newDate between '2014-01-01' and '2014-01-31'")
# feb2014 <- sqldf("Select * from mosq2014 where newDate between '2014-02-01' and '2014-02-31'")
# mar2014 <- sqldf("Select * from mosq2014 where newDate between '2014-03-01' and '2014-03-31'")
# apr2014 <- sqldf("Select * from mosq2014 where newDate between '2014-04-01' and '2014-04-31'")
# may2014 <- sqldf("Select * from mosq2014 where newDate between '2014-05-01' and '2014-05-31'")
# jun2014 <- sqldf("Select * from mosq2014 where newDate between '2014-06-01' and '2014-06-31'")
# jul2014 <- sqldf("Select * from mosq2014 where newDate between '2014-07-01' and '2014-07-31'")
# aug2014 <- sqldf("Select * from mosq2014 where newDate between '2014-08-01' and '2014-08-31'")
# sep2014 <- sqldf("Select * from mosq2014 where newDate between '2014-09-01' and '2014-09-31'")
# oct2014 <- sqldf("Select * from mosq2014 where newDate between '2014-10-01' and '2014-10-31'")
# 

# Select all from buffer zone
buff.x <- sqldf("Select * from kdrData where project_code is not 'treatment'")
buff.x <- buff[!is.na(buff$project_code),]

# Select based on 2013 month from buffer zone
jan2013b.x <- sqldf("Select * from buff.x where newDate between '2013-01-01' and '2013-01-31'")
feb2013b.x <- sqldf("Select * from buff.x where newDate between '2013-02-01' and '2013-02-31'")
mar2013b.x <- sqldf("Select * from buff.x where newDate between '2013-03-01' and '2013-03-31'")
apr2013b.x <- sqldf("Select * from buff.x where newDate between '2013-04-01' and '2013-04-31'")
may2013b.x <- sqldf("Select * from buff.x where newDate between '2013-05-01' and '2013-05-31'")
jun2013b.x <- sqldf("Select * from buff.x where newDate between '2013-06-01' and '2013-06-31'")
jul2013b.x <- sqldf("Select * from buff.x where newDate between '2013-07-01' and '2013-07-31'")
aug2013b.x <- sqldf("Select * from buff.x where newDate between '2013-08-01' and '2013-08-31'")
sep2013b.x <- sqldf("Select * from buff.x where newDate between '2013-09-01' and '2013-09-31'")
oct2013b.x <- sqldf("Select * from buff.x where newDate between '2013-10-01' and '2013-10-31'")


# Run function across all months in 2013 buffer zone ------------------------------------------
# For 1016 locus at all months in 2013 - buffer zone
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


# Create dataframes -------------------------------------------------------
# Create list of years included in dataframe
year <- c(2000:2017)
# Create list of months included in dataframe - as numeric
month <- c(1:10)

# For 1016 locus at all months in 2013 - buffer zone
dfbuff13.x <- rbind(bJan13.x, bFeb13.x, bMar13.x, bApr13.x, bMay13.x, bJun13.x, bJul13.x, bAug13.x, bSep13.x, bOct13.x)
# Add year ID to rows in df rename
mc.1016.b13.expanded <- cbind(month, dfbuff13.x)

# ### To view dataframes
# mc.1016.b13.expanded

# # Save dataframes ---------------------------------------------------------
# # These are required for plots, selection coefficient, and other analyses
# write.csv(mc.1016.b13, file = "./data/mc.1016.b13_reduced_expandedBuffer.csv", row.names = F)

