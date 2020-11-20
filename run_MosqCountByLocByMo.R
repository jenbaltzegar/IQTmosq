# Create spreadsheet with number of individual mosquitoes collected per house per month

# Prepare working environment ---------------------------------------------
# set working directory
setwd("~/Dropbox/GouldLab/Project_Mosquito/Database")

# clear working environment
rm(list = ls())

# load libraries
library(lubridate)
library(dplyr)
library(ggplot2)

# Load & Prep Data -----------------
# load df
df <- read.csv("kdrData_all_reduced.csv", header = TRUE)

# Convert Date to readable form for R
df$newDate <- as.character(as.Date(as.character(df$Date), format = "%m/%d/%Y"))
# Note: If year starts with 00 after conversion, then go back and save IsolationLog_IQT-Mosq.csv
# with Date column in the correct format (MM/DD/YYYY). Weird formating issues happen with .csv files

# reduce df
df <- df[,c(1,2,11)]

# change newDate from character to Date obj
df$newDate <- ymd(df$newDate)
df$by_month <- floor_date(df$newDate, "month")

# restrict df to certain years
df <- df[which(df$by_month > "1999-12-01" & df$by_month < "2007-01-01"),]

# Sort by date
df <- df[order(df$by_month, df$Location_Code),]
head(df)

# Summarize data ---------------
# group data by date and location code and count number of indiv. mosq.
df.count <- df %>% 
  dplyr::group_by(by_month, Location_Code) %>% 
  dplyr::summarise(num_indiv = n_distinct(mosquito_id)) 
# convert back to df
df.count <- as.data.frame(df.count)
# head(df.count)
# tail(df.count)
# df.count[df.count$by_month == "2002-01-01",]

# group data by date and location code and calculate mean number of indiv. mosq / mo
df.mean <- df.count %>%
  dplyr::group_by(by_month) %>%
  dplyr::summarise(n_houses = n_distinct(Location_Code)
    , avg_indiv_mo = mean(num_indiv))
# convert back to df
df.mean <- as.data.frame(df.mean)
# head(df.mean)
# tail(df.mean)
# df.mean[df.mean$by_month == "2002-01-01",]

# write dfs to csv. ---------------
write.csv(df.count, "./MosqCount_PerLocPerMonth.csv", row.names = FALSE)
write.csv(df.mean, "./MosqMean_PerHousePerMonth.csv", row.names = FALSE)
