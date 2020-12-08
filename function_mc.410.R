###################################################################################################
### Function to analyze melt curve data
### Creates dataframe of genotype counts, allele frequency of R, and +/- 95% confidence interval
# For 410 Locus
# Apply to each year 
# Created 23 Jan 2019
###################################################################################################

### Function to create dataframe of genotype counts, allele frequency of R, and +/- 95% confidence interval
# For 410 Locus
# Apply to each year 

mc.410 <- function(objectName){
  # Remove rows with errors
  objectName <- sqldf(c("Delete from objectName where V410L_converted like 'error'", "select * from objectName"))
  
  # Remove rows with NAs
  objectName <- sqldf("select * from objectName where V410L_converted is not null")
  
  # Count genotypes per locus per year
  count410 <- sqldf("select V410L_converted, count (V410L_converted) as countGenos from objectName group by V410L_converted order by countGenos")
  
  # Subset count1016 on RR & SR
  countRR = count410[count410$V410L_converted=="RR", "countGenos"]
  countSR = count410[count410$V410L_converted=="SR", "countGenos"]
  countSS = count410[count410$V410L_converted=="SS", "countGenos"]
  
  # Convert counts to usable numbers
  scalarRR = sum(countRR)
  scalarSR = sum(countSR)
  scalarSS = sum(countSS)
  
  # Calculate N
  n = scalarRR+scalarSR+scalarSS
  
  # Calculate frequency of R allele
  freqR = (2*scalarRR + scalarSR)/(2*n)
  
  # Calculate 95% Confidence Interval
  CI_95 = 1.96 * sqrt((freqR*(1-freqR))/(2*n))
  
  
  # # Create output
  # print("Genotype Summary Table")
  # print(count410)
  
  # Create data frame of results
  df <- data.frame(scalarRR, scalarSR, scalarSS, n, freqR, CI_95)
  colnames(df) <- c("RR", "SR", "SS", "n", "freqR", "CI_95")
  return(df)
}


# ### Run this code for testing output and creating larger dataframe
# #library(sqldf)
# tmp_Result0 = mc.1016.yr(mosq2000)
# print(tmp_Result0)
# 
# tmp_Result1 = mc.1016.yr(mosq2001)
# print(tmp_Result1)
# 
# # Create data frame from all runs of function
# dfNew <- rbind(tmp_Result0, tmp_Result1)
# # Create list of years included in dataframe
# year <- c(2000:2001)
# # Add year ID to rows in dfNew
# cbind(year, dfNew)
