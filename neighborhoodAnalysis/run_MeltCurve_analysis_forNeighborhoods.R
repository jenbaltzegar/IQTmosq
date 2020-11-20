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


# Load data -------------------------------------------------------
kdrData <- read.csv("kdrData_reduced.csv")


# Subset kdrData by latitude --------------------------------------------
# Remove rows where Y = NA
temp <- subset(kdrData, !is.na(Y))
# Remove rows with NA and errors in haplotype column
keep <- c("SSSS", "RRRR", "SSRR", "SSSR", "SRRR", "SRSS", "SRSR", "RRSS")
temp <- temp[temp$haplotype %in% keep,]
# Coerce date to be a date
temp$newDate <- as.Date(temp$newDate)

# Separate samples based on latitude and/or longitude
top <- temp[temp$Y >= '9585862',]
mid <- temp[temp$Y < '9586780' & temp$Y > '9583180',]
bot <- temp[temp$Y <= '9583180',]
top.L <- temp[temp$Y >= '9586780' & temp$X <= '694100',]
top.R <- temp[temp$Y >= '9586780' & temp$X > '694100',]
no.bot <- temp[temp$Y >= '9583180',] # for michael vella modeling fitting 2020-05-07

# Save list of individuals for plotting
write.csv(top, "../QGIS_Files/Routput/north.csv", row.names = FALSE)
write.csv(top.R, "../QGIS_Files/Routput/northEast.csv", row.names = FALSE)
write.csv(top.L, "../QGIS_Files/Routput/northWest.csv", row.names = FALSE)
write.csv(mid, "../QGIS_Files/Routput/central.csv", row.names = FALSE)
write.csv(bot, "../QGIS_Files/Routput/south.csv", row.names = FALSE)
write.csv(no.bot, "../QGIS_Files/Routput/noSouth.csv", row.names = FALSE) # for michael vella modeling fitting 2020-05-07

# Subset top based on year - "North"
top2000 <- top[top$newDate > '2000-01-01' & top$newDate < '2000-12-31',]
top2001 <- top[top$newDate > '2001-01-01' & top$newDate < '2001-12-31',]
top2002 <- top[top$newDate > '2002-01-01' & top$newDate < '2002-12-31',]
top2003 <- top[top$newDate > '2003-01-01' & top$newDate < '2003-12-31',]
top2004 <- top[top$newDate > '2004-01-01' & top$newDate < '2004-12-31',]
top2005 <- top[top$newDate > '2005-01-01' & top$newDate < '2005-12-31',]
top2006 <- top[top$newDate > '2006-01-01' & top$newDate < '2006-12-31',]
top2007 <- top[top$newDate > '2007-01-01' & top$newDate < '2007-12-31',]
top2008 <- top[top$newDate > '2008-01-01' & top$newDate < '2008-12-31',]
top2009 <- top[top$newDate > '2009-01-01' & top$newDate < '2009-12-31',]
top2010 <- top[top$newDate > '2010-01-01' & top$newDate < '2010-12-31',]
top2011 <- top[top$newDate > '2011-01-01' & top$newDate < '2011-12-31',]
top2012 <- top[top$newDate > '2012-01-01' & top$newDate < '2012-12-31',]
top2013 <- top[top$newDate > '2013-01-01' & top$newDate < '2013-12-31',]
top2014 <- top[top$newDate > '2014-01-01' & top$newDate < '2014-12-31',]
top2015 <- top[top$newDate > '2015-01-01' & top$newDate < '2015-12-31',]
top2016 <- top[top$newDate > '2016-01-01' & top$newDate < '2016-12-31',]
top2017 <- top[top$newDate > '2017-01-01' & top$newDate < '2017-12-31',]

# Subset top.L based on year - "North West"
top.L2000 <- top.L[top.L$newDate > '2000-01-01' & top.L$newDate < '2000-12-31',]
top.L2001 <- top.L[top.L$newDate > '2001-01-01' & top.L$newDate < '2001-12-31',]
top.L2002 <- top.L[top.L$newDate > '2002-01-01' & top.L$newDate < '2002-12-31',]
top.L2003 <- top.L[top.L$newDate > '2003-01-01' & top.L$newDate < '2003-12-31',]
top.L2004 <- top.L[top.L$newDate > '2004-01-01' & top.L$newDate < '2004-12-31',]
top.L2005 <- top.L[top.L$newDate > '2005-01-01' & top.L$newDate < '2005-12-31',]
top.L2006 <- top.L[top.L$newDate > '2006-01-01' & top.L$newDate < '2006-12-31',]
top.L2007 <- top.L[top.L$newDate > '2007-01-01' & top.L$newDate < '2007-12-31',]
top.L2008 <- top.L[top.L$newDate > '2008-01-01' & top.L$newDate < '2008-12-31',]
top.L2009 <- top.L[top.L$newDate > '2009-01-01' & top.L$newDate < '2009-12-31',]
top.L2010 <- top.L[top.L$newDate > '2010-01-01' & top.L$newDate < '2010-12-31',]
top.L2011 <- top.L[top.L$newDate > '2011-01-01' & top.L$newDate < '2011-12-31',]
top.L2012 <- top.L[top.L$newDate > '2012-01-01' & top.L$newDate < '2012-12-31',]
top.L2013 <- top.L[top.L$newDate > '2013-01-01' & top.L$newDate < '2013-12-31',]
top.L2014 <- top.L[top.L$newDate > '2014-01-01' & top.L$newDate < '2014-12-31',]
top.L2015 <- top.L[top.L$newDate > '2015-01-01' & top.L$newDate < '2015-12-31',]
top.L2016 <- top.L[top.L$newDate > '2016-01-01' & top.L$newDate < '2016-12-31',]
top.L2017 <- top.L[top.L$newDate > '2017-01-01' & top.L$newDate < '2017-12-31',]

# Subset top.R based on year - "North East"
top.R2000 <- top.R[top.R$newDate > '2000-01-01' & top.R$newDate < '2000-12-31',]
top.R2001 <- top.R[top.R$newDate > '2001-01-01' & top.R$newDate < '2001-12-31',]
top.R2002 <- top.R[top.R$newDate > '2002-01-01' & top.R$newDate < '2002-12-31',]
top.R2003 <- top.R[top.R$newDate > '2003-01-01' & top.R$newDate < '2003-12-31',]
top.R2004 <- top.R[top.R$newDate > '2004-01-01' & top.R$newDate < '2004-12-31',]
top.R2005 <- top.R[top.R$newDate > '2005-01-01' & top.R$newDate < '2005-12-31',]
top.R2006 <- top.R[top.R$newDate > '2006-01-01' & top.R$newDate < '2006-12-31',]
top.R2007 <- top.R[top.R$newDate > '2007-01-01' & top.R$newDate < '2007-12-31',]
top.R2008 <- top.R[top.R$newDate > '2008-01-01' & top.R$newDate < '2008-12-31',]
top.R2009 <- top.R[top.R$newDate > '2009-01-01' & top.R$newDate < '2009-12-31',]
top.R2010 <- top.R[top.R$newDate > '2010-01-01' & top.R$newDate < '2010-12-31',]
top.R2011 <- top.R[top.R$newDate > '2011-01-01' & top.R$newDate < '2011-12-31',]
top.R2012 <- top.R[top.R$newDate > '2012-01-01' & top.R$newDate < '2012-12-31',]
top.R2013 <- top.R[top.R$newDate > '2013-01-01' & top.R$newDate < '2013-12-31',]
top.R2014 <- top.R[top.R$newDate > '2014-01-01' & top.R$newDate < '2014-12-31',]
top.R2015 <- top.R[top.R$newDate > '2015-01-01' & top.R$newDate < '2015-12-31',]
top.R2016 <- top.R[top.R$newDate > '2016-01-01' & top.R$newDate < '2016-12-31',]
top.R2017 <- top.R[top.R$newDate > '2017-01-01' & top.R$newDate < '2017-12-31',]

# Subset mid based on year - "Central"
mid2000 <- mid[mid$newDate > '2000-01-01' & mid$newDate < '2000-12-31',]
mid2001 <- mid[mid$newDate > '2001-01-01' & mid$newDate < '2001-12-31',]
mid2002 <- mid[mid$newDate > '2002-01-01' & mid$newDate < '2002-12-31',]
mid2003 <- mid[mid$newDate > '2003-01-01' & mid$newDate < '2003-12-31',]
mid2004 <- mid[mid$newDate > '2004-01-01' & mid$newDate < '2004-12-31',]
mid2005 <- mid[mid$newDate > '2005-01-01' & mid$newDate < '2005-12-31',]
mid2006 <- mid[mid$newDate > '2006-01-01' & mid$newDate < '2006-12-31',]
mid2007 <- mid[mid$newDate > '2007-01-01' & mid$newDate < '2007-12-31',]
mid2008 <- mid[mid$newDate > '2008-01-01' & mid$newDate < '2008-12-31',]
mid2009 <- mid[mid$newDate > '2009-01-01' & mid$newDate < '2009-12-31',]
mid2010 <- mid[mid$newDate > '2010-01-01' & mid$newDate < '2010-12-31',]
mid2011 <- mid[mid$newDate > '2011-01-01' & mid$newDate < '2011-12-31',]
mid2012 <- mid[mid$newDate > '2012-01-01' & mid$newDate < '2012-12-31',]
mid2013 <- mid[mid$newDate > '2013-01-01' & mid$newDate < '2013-12-31',]
mid2014 <- mid[mid$newDate > '2014-01-01' & mid$newDate < '2014-12-31',]
mid2015 <- mid[mid$newDate > '2015-01-01' & mid$newDate < '2015-12-31',]
mid2016 <- mid[mid$newDate > '2016-01-01' & mid$newDate < '2016-12-31',]
mid2017 <- mid[mid$newDate > '2017-01-01' & mid$newDate < '2017-12-31',]

# Subset bot based on year - "South"
bot2000 <- bot[bot$newDate > '2000-01-01' & bot$newDate < '2000-12-31',]
bot2001 <- bot[bot$newDate > '2001-01-01' & bot$newDate < '2001-12-31',]
bot2002 <- bot[bot$newDate > '2002-01-01' & bot$newDate < '2002-12-31',]
bot2003 <- bot[bot$newDate > '2003-01-01' & bot$newDate < '2003-12-31',]
bot2004 <- bot[bot$newDate > '2004-01-01' & bot$newDate < '2004-12-31',]
bot2005 <- bot[bot$newDate > '2005-01-01' & bot$newDate < '2005-12-31',]
bot2006 <- bot[bot$newDate > '2006-01-01' & bot$newDate < '2006-12-31',]
bot2007 <- bot[bot$newDate > '2007-01-01' & bot$newDate < '2007-12-31',]
bot2008 <- bot[bot$newDate > '2008-01-01' & bot$newDate < '2008-12-31',]
bot2009 <- bot[bot$newDate > '2009-01-01' & bot$newDate < '2009-12-31',]
bot2010 <- bot[bot$newDate > '2010-01-01' & bot$newDate < '2010-12-31',]
bot2011 <- bot[bot$newDate > '2011-01-01' & bot$newDate < '2011-12-31',]
bot2012 <- bot[bot$newDate > '2012-01-01' & bot$newDate < '2012-12-31',]
bot2013 <- bot[bot$newDate > '2013-01-01' & bot$newDate < '2013-12-31',]
bot2014 <- bot[bot$newDate > '2014-01-01' & bot$newDate < '2014-12-31',]
bot2015 <- bot[bot$newDate > '2015-01-01' & bot$newDate < '2015-12-31',]
bot2016 <- bot[bot$newDate > '2016-01-01' & bot$newDate < '2016-12-31',]
bot2017 <- bot[bot$newDate > '2017-01-01' & bot$newDate < '2017-12-31',]



# Source functions to create dataframe of genotype counts, allele  --------
# Source functions
source("R_Scripts/IQTmosq/function_mc.haps.R")

# Run functions for top, mid, & bot ---------------------------------------
# For haplotypes at all years for top, mid, and bot
top.haps.2000 = mc.haps(top2000)
top.haps.2001 = mc.haps(top2001)
top.haps.2002 = mc.haps(top2002)
top.haps.2003 = mc.haps(top2003)
top.haps.2004 = mc.haps(top2004)
top.haps.2005 = mc.haps(top2005)
top.haps.2006 = mc.haps(top2006)
top.haps.2007 = mc.haps(top2007)
top.haps.2008 = mc.haps(top2008)
top.haps.2009 = mc.haps(top2009)
top.haps.2010 = mc.haps(top2010)
top.haps.2011 = mc.haps(top2011)
top.haps.2012 = mc.haps(top2012)
top.haps.2013 = mc.haps(top2013)
top.haps.2014 = mc.haps(top2014)
top.haps.2015 = mc.haps(top2015)
top.haps.2016 = mc.haps(top2016)
top.haps.2017 = mc.haps(top2017)

top.L.haps.2000 = mc.haps(top.L2000)
top.L.haps.2001 = mc.haps(top.L2001)
top.L.haps.2002 = mc.haps(top.L2002)
top.L.haps.2003 = mc.haps(top.L2003)
top.L.haps.2004 = mc.haps(top.L2004)
top.L.haps.2005 = mc.haps(top.L2005)
top.L.haps.2006 = mc.haps(top.L2006)
top.L.haps.2007 = mc.haps(top.L2007)
top.L.haps.2008 = mc.haps(top.L2008)
top.L.haps.2009 = mc.haps(top.L2009)
top.L.haps.2010 = mc.haps(top.L2010)
top.L.haps.2011 = mc.haps(top.L2011)
top.L.haps.2012 = mc.haps(top.L2012)
top.L.haps.2013 = mc.haps(top.L2013)
top.L.haps.2014 = mc.haps(top.L2014)
top.L.haps.2015 = mc.haps(top.L2015)
top.L.haps.2016 = mc.haps(top.L2016)
top.L.haps.2017 = mc.haps(top.L2017)

top.R.haps.2000 = mc.haps(top.R2000)
top.R.haps.2001 = mc.haps(top.R2001)
top.R.haps.2002 = mc.haps(top.R2002)
top.R.haps.2003 = mc.haps(top.R2003)
top.R.haps.2004 = mc.haps(top.R2004)
top.R.haps.2005 = mc.haps(top.R2005)
top.R.haps.2006 = mc.haps(top.R2006)
top.R.haps.2007 = mc.haps(top.R2007)
top.R.haps.2008 = mc.haps(top.R2008)
top.R.haps.2009 = mc.haps(top.R2009)
top.R.haps.2010 = mc.haps(top.R2010)
top.R.haps.2011 = mc.haps(top.R2011)
top.R.haps.2012 = mc.haps(top.R2012)
top.R.haps.2013 = mc.haps(top.R2013)
top.R.haps.2014 = mc.haps(top.R2014)
top.R.haps.2015 = mc.haps(top.R2015)
top.R.haps.2016 = mc.haps(top.R2016)
top.R.haps.2017 = mc.haps(top.R2017)

mid.haps.2000 = mc.haps(mid2000)
mid.haps.2001 = mc.haps(mid2001)
mid.haps.2002 = mc.haps(mid2002)
mid.haps.2003 = mc.haps(mid2003)
mid.haps.2004 = mc.haps(mid2004)
mid.haps.2005 = mc.haps(mid2005)
mid.haps.2006 = mc.haps(mid2006)
mid.haps.2007 = mc.haps(mid2007)
mid.haps.2008 = mc.haps(mid2008)
mid.haps.2009 = mc.haps(mid2009)
mid.haps.2010 = mc.haps(mid2010)
mid.haps.2011 = mc.haps(mid2011)
mid.haps.2012 = mc.haps(mid2012)
mid.haps.2013 = mc.haps(mid2013)
mid.haps.2014 = mc.haps(mid2014)
mid.haps.2015 = mc.haps(mid2015)
mid.haps.2016 = mc.haps(mid2016)
mid.haps.2017 = mc.haps(mid2017)

bot.haps.2000 = mc.haps(bot2000)
bot.haps.2001 = mc.haps(bot2001)
bot.haps.2002 = mc.haps(bot2002)
bot.haps.2003 = mc.haps(bot2003)
bot.haps.2004 = mc.haps(bot2004)
bot.haps.2005 = mc.haps(bot2005)
bot.haps.2006 = mc.haps(bot2006)
bot.haps.2007 = mc.haps(bot2007)
bot.haps.2008 = mc.haps(bot2008)
bot.haps.2009 = mc.haps(bot2009)
bot.haps.2010 = mc.haps(bot2010)
bot.haps.2011 = mc.haps(bot2011)
bot.haps.2012 = mc.haps(bot2012)
bot.haps.2013 = mc.haps(bot2013)
bot.haps.2014 = mc.haps(bot2014)
bot.haps.2015 = mc.haps(bot2015)
bot.haps.2016 = mc.haps(bot2016)
bot.haps.2017 = mc.haps(bot2017)

# Create dataframes -------------------------------------------------------
# Create list of years included in dataframe
year <- c(2000:2017)

#### haps
# For top
dfTop.haps <- rbind(top.haps.2000, top.haps.2001, top.haps.2002, top.haps.2003, top.haps.2004
                    , top.haps.2005, top.haps.2006, top.haps.2007, top.haps.2008, top.haps.2009
                    , top.haps.2010, top.haps.2011, top.haps.2012, top.haps.2013, top.haps.2014
                    , top.haps.2015, top.haps.2016, top.haps.2017)
mc.haps.top <- cbind(dfTop.haps, year)

# For top.L
dftop.L.haps <- rbind(top.L.haps.2000, top.L.haps.2001, top.L.haps.2002, top.L.haps.2003, top.L.haps.2004
                    , top.L.haps.2005, top.L.haps.2006, top.L.haps.2007, top.L.haps.2008, top.L.haps.2009
                    , top.L.haps.2010, top.L.haps.2011, top.L.haps.2012, top.L.haps.2013, top.L.haps.2014
                    , top.L.haps.2015, top.L.haps.2016, top.L.haps.2017)
mc.haps.top.L <- cbind(dftop.L.haps, year)

# For top.R
dftop.R.haps <- rbind(top.R.haps.2000, top.R.haps.2001, top.R.haps.2002, top.R.haps.2003, top.R.haps.2004
                    , top.R.haps.2005, top.R.haps.2006, top.R.haps.2007, top.R.haps.2008, top.R.haps.2009
                    , top.R.haps.2010, top.R.haps.2011, top.R.haps.2012, top.R.haps.2013, top.R.haps.2014
                    , top.R.haps.2015, top.R.haps.2016, top.R.haps.2017)
mc.haps.top.R <- cbind(dftop.R.haps, year)

# For mid
dfmid.haps <- rbind(mid.haps.2000, mid.haps.2001, mid.haps.2002, mid.haps.2003, mid.haps.2004
                    , mid.haps.2005, mid.haps.2006, mid.haps.2007, mid.haps.2008, mid.haps.2009
                    , mid.haps.2010, mid.haps.2011, mid.haps.2012, mid.haps.2013, mid.haps.2014
                    , mid.haps.2015, mid.haps.2016, mid.haps.2017)
mc.haps.mid <- cbind(dfmid.haps, year)

# For bot
dfbot.haps <- rbind(bot.haps.2000, bot.haps.2001, bot.haps.2002, bot.haps.2003, bot.haps.2004
                    , bot.haps.2005, bot.haps.2006, bot.haps.2007, bot.haps.2008, bot.haps.2009
                    , bot.haps.2010, bot.haps.2011, bot.haps.2012, bot.haps.2013, bot.haps.2014
                    , bot.haps.2015, bot.haps.2016, bot.haps.2017)
mc.haps.bot <- cbind(dfbot.haps, year)

# Save dataframes ---------------------------------------------------------
# These are required for plots, selection coefficient, and other analyses
write.csv(mc.haps.top, file = "/Users/jenbaltz/Dropbox/GouldLab/Project_Mosquito/Database/mc.haps.top.csv", row.names = F)
write.csv(mc.haps.top.L, file = "/Users/jenbaltz/Dropbox/GouldLab/Project_Mosquito/Database/mc.haps.top.L.csv", row.names = F)
write.csv(mc.haps.top.R, file = "/Users/jenbaltz/Dropbox/GouldLab/Project_Mosquito/Database/mc.haps.top.R.csv", row.names = F)
write.csv(mc.haps.mid, file = "/Users/jenbaltz/Dropbox/GouldLab/Project_Mosquito/Database/mc.haps.mid.csv", row.names = F)
write.csv(mc.haps.bot, file = "/Users/jenbaltz/Dropbox/GouldLab/Project_Mosquito/Database/mc.haps.bot.csv", row.names = F)