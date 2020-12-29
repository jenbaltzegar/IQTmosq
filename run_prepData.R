# This script will read in all data and merge into one dataframe called kdrData

# Prep data ---------------------------------------------------------------
# Convert Date to readable form for R
kdrData$newDate <- as.character(as.Date(as.character(kdrData$Date), format = "%m/%d/%Y"))
# Note: If year starts with 00 after conversion, then go back and save IsolationLog_IQT-Mosq.csv
# with Date column in the correct format (MM/DD/YYYY). Weird formating issues happen with .csv files

kdrData <- (
    kdrData
    # Left join masterdec GIS and neighborhood info onto kdrData
    %.>% merge(
        x=., y=masterdec[, c("X", "Y", "NEIGHBORHO", "LOC_CODE")],
        by.x = "Location_Code", by.y = "LOC_CODE", 
        all.x = TRUE
    )
    # Left join Zone information to kdrDate
    %.>% merge(
        x = ., y = exptZone, 
        by.x = "Location_Code", by.y = "location_code", 
        all.x = TRUE
    )
)


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
kdrData <- within(kdrData, {
    V1016I <- ifelse(V1016I_rep1 == V1016I_rep2, 
        as.character(V1016I_rep1), "error"
    )
    F1534C <- ifelse(F1534C_rep1 == F1534C_rep2,
        as.character(F1534C_rep1), "error"
    )
    V410L <- ifelse(V410L_rep1 == V410L_rep2, 
    as.character(V410L_rep1), "error"
    )
    # Convert melt curve output to readable genotype and haplotype
    # Uses function_convertMergedMeltCurve.R
    # Run function for each column
    V1016I_converted <- convert1016(V1016I)
    F1534C_converted <- convert1534(F1534C)
    V410L_converted <- convert410(V410L)
    # Create haplotype for loci 1016 and 1534 - 410 is not included here
    haplotype <- paste0(V1016I_converted, F1534C_converted)
})

####
# Key Code: Remove unknown location codes from kdrData and save file
####
# "Unknown" location codes are typically from mosquitoes collected by Steve Stoddard outside of IQT
kdrData <- subset(kdrData, 
    Location_Code != "unknown"
    ## only 7
    & !duplicated(mosquito_id)
)
# write.csv(kdrData,"~/Dropbox/GouldLab/Project_Mosquito/Database/kdrData_all_reduced.csv", row.names = F)

### Save file
# Write to file
write.csv(kdrData, "./data/kdrData.csv", row.names = F)
