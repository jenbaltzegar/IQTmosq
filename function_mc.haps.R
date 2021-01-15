###################################################################################################
### Function to analyze melt curve data
### Creates dataframe of haplotype counts, allele frequency of R, and +/- 95% confidence interval
# For 1534 Locus
# Apply to each year 
# Created 6-16-17
# Last edit 23 Oct 2017
###################################################################################################

### Function to create dataframe of genotype counts, allele frequency of R, and +/- 95% confidence interval
# For 1534 Locus
# Apply to each year 

# objectName <- top2000
# unique(objectName$haplotype)

mc.haps <- function(objectName){
  # Remove rows with errors and NAs
  objectName <- subset(objectName,
    !(
        grepl("*error*", haplotype)
        | grepl("*NA*", haplotype)
        | is.na(haplotype)
    )
  )
  # Count genotypes per locus per year
  countHaps <- sqldf("select haplotype, count (haplotype) as countGenos from objectName group by haplotype order by countGenos")

  ## get counts for each type (0 if not found)
  types <- c(
    "SSSS", "SSSR", "SSRR", 
    "SRSS", "SRSR", "SRRR", 
    "RRSS", "RRSR", "RRRR"
  )
  result <- lapply(types, function(name) {
    ## sum treats missing values as zero
    sum(subset(countHaps, haplotype==name)$countGenos)
  })
  ## set column names
  names(result) <- types
  ## totals
  result$n <- sum(countHaps$countGenos)
  return(as.data.frame(result))
}


# ### Run this code for testing output and creating larger dataframe
# #library(sqldf)
# tmp_Result0 = mc.haps(top2000)
# print(tmp_Result0)
# 
# tmp_Result1 = mc.haps(top2001)
# print(tmp_Result1)
# 
# # Create data frame from all runs of function
# dfNew <- rbind(tmp_Result0, tmp_Result1)
# # Create list of years included in dataframe
# year <- c(2000:2001)
# # Add year ID to rows in dfNew
# cbind(year, dfNew)
# 

