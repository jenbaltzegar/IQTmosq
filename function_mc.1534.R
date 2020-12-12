###################################################################################################
### Function to analyze melt curve data
### Creates dataframe of genotype counts, allele frequency of R, and +/- 95% confidence interval
# For 1534 Locus
# Apply to each year 
# Created 6-16-17
# Last edit 6-16-17
###################################################################################################

### Function to create dataframe of genotype counts, allele frequency of R, and +/- 95% confidence interval
# For 1534 Locus
# Apply to each year 

## depends on fun.ci in function_shared.R
mc.1534 <- function(objectName){

  # omit rows with errors or NAs
  objectName <- subset(objectName, 
    !( 
        grepl('error', F1534C_converted) 
        | is.na(F1534C_converted)
    )
  )

  # Remove rows with errors
  #objectName <- sqldf(c("Delete from objectName where F1534C_converted like 'error'", "select * from objectName"))
  
  # Remove rows with NAs
  #objectName <- sqldf("select * from objectName where F1534C_converted is not null")
  
  # Count genotypes per locus per year
  count1534 <- sqldf("select F1534C_converted, count (F1534C_converted) as countGenos from objectName group by F1534C_converted order by countGenos")
  
  # Subset count1534 on RR & SR
  countRR = count1534[count1534$F1534C_converted=="RR", "countGenos"]
  countSR = count1534[count1534$F1534C_converted=="SR", "countGenos"]
  countSS = count1534[count1534$F1534C_converted=="SS", "countGenos"]
  
  # Convert counts to usable numbers
  scalarRR = sum(countRR)
  scalarSR = sum(countSR)
  scalarSS = sum(countSS)
  
  # Calculate N
  n = scalarRR+scalarSR+scalarSS
  
  # Calculate frequency of R allele
  freqR = (2*scalarRR + scalarSR)/(2*n)
  
  # Calculate 95% Confidence Interval
  CI_95 = fun.ci(freqR, n)
  
  # # Create output
  # print("Genotype Summary Table")
  # print(count1534)
  
  # Create data frame of results
  df <- data.frame(scalarRR, scalarSR, scalarSS, n, freqR, CI_95)
  colnames(df) <- c("RR", "SR", "SS", "n", "freqR", "CI_95")
  return(df)
}


# ### Run this code for testing output and creating larger dataframe
# #library(sqldf)
# tmp_Result0 = mc.1534.yr(mosq2000)
# print(tmp_Result0)
# 
# tmp_Result1 = mc.1534.yr(mosq2001)
# print(tmp_Result1)
# 
# # Create data frame from all runs of function
# dfNew <- rbind(tmp_Result0, tmp_Result1)
# # Create list of years included in dataframe
# year <- c(2000:2001)
# # Add year ID to rows in dfNew
# cbind(year, dfNew)


