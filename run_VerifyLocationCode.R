## This program is to analyze data generated from Jungkoo Kang's program
## called BlockLevelSampleCount which is implemented in Ubuntu. Output files are 
## saved in ~Dropbox/GouldLab/Database
## Started: 19 Sep 2017
## Last edit: 20 Sep 2017

# Set working directory
setwd("/Users/jenbaltz/Dropbox/GouldLab/Project_Mosquito/Database")

# Import data
raw <- read.csv("BlockLevelSampleCount_output_Data_Raw.csv", sep = ",", header = T)
summary <- read.csv("BlockLevelSampleCount_output_Data_Summary.csv", sep = ",", header = T)
isol <- read.csv("IsolationLog_IQT-Mosq.csv", sep = ",", header = T)

# View data
head(raw); tail(raw)
head(summary); tail(summary)

# ###############################################################################################
# # Find rows that end in P - this means that they were caught outside (i.e., Peri) of the house
# # view head of isol to make sure that Location Code column is still present
# head(isol)
# 
# # Add column called "collected_outside" if location code ends in "P"
# isol$collected_outside <- ifelse(grepl("*[P]$", isol$Location_Code), "yes", "no")
# 
# # # Verify this was added correctly
# # # return rows where Location Code ends in "P"
# # peri <- isol[grep("*[P]$", isol$Location_Code), ]
# # head(peri)
# 
# # Remove the final P from the end of Location_Code
# # head(isol)
# isol$Location_Code <- gsub(pattern = "*[P]$", replacement = "", isol$Location_Code)
# # head(isol)
# 
# # Save isolation log with new column
# write.csv(isol, "IsolationLog_IQT-Mosq.csv", row.names = F)


# ###############################################################################################
# ### Create and save file containing list of mosquitoes with missing location
# # Subset df "raw" where Block is missing
# miss <- raw[raw$Block=="Missing location",]
# # Subset again where location_code is not unknown
# miss <- miss[miss$Location_Code!="unknow",]
# 
# # subset to remove Location code from isol b/c redundant in miss
# #isol <- isol[-3]
# 
# # left join isolation log info to miss
# miss.isol <- merge(miss, isol, by="mosquito_id", all.x=TRUE)
# 
# # Save file to check location codes in database & tubes
# write.csv(miss.isol, "VerifyLocationCode.csv", row.names = F)


###############################################################################################
### Create and save file containing list of mosquitoes with location codes not in city db
# load data
df <- read.csv("kdrData_reduced.csv")
# subset samples that do not have a longitude GPS code
new_df <- df[is.na(df$X),]
# Save file to check location codes in database and tubes\
write.csv(new_df, paste0("./VerifyLocationCodes_", Sys.Date(), ".csv"), row.names = FALSE)



