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
  
  # Subset countHaps on haplotypes
  countSRSR = countHaps[countHaps$haplotype=="SRSR", "countGenos"]
  countRRSS = countHaps[countHaps$haplotype=="RRSS", "countGenos"]
  countSRSS = countHaps[countHaps$haplotype=="SRSS", "countGenos"]
  countSSSR = countHaps[countHaps$haplotype=="SSSR", "countGenos"]
  countSSRR = countHaps[countHaps$haplotype=="SSRR", "countGenos"]
  countSRRR = countHaps[countHaps$haplotype=="SRRR", "countGenos"]
  countSSSS = countHaps[countHaps$haplotype=="SSSS", "countGenos"]
  countRRRR = countHaps[countHaps$haplotype=="RRRR", "countGenos"]
  countRRSR = countHaps[countHaps$haplotype=="RRSR", "countGenos"]
  
  # Convert counts to useable numbers
  scalarSRSR = sum(countSRSR)
  scalarRRSS = sum(countRRSS)
  scalarSRSS = sum(countSRSS)
  scalarSSSR = sum(countSSSR)
  scalarSSRR = sum(countSSRR)
  scalarSRRR = sum(countSRRR)
  scalarSSSS = sum(countSSSS)
  scalarRRRR = sum(countRRRR)
  scalarRRSR = sum(countRRSR)
  
  # Calculate N
  n = scalarSSSS + scalarSSSR + scalarSSRR + scalarSRSS + scalarSRSR + scalarSRRR + scalarRRSS + scalarRRSR + scalarRRRR
  
  # Calculate frequency of R allele
  freqSSSS = (scalarSSSS)/n
  freqSSSR = (scalarSSSR)/n
  freqSSRR = (scalarSSRR)/n
  freqSRSS = (scalarSRSS)/n
  freqSRSR = (scalarSRSR)/n
  freqSRRR = (scalarSRRR)/n
  freqRRSS = (scalarRRSS)/n
  freqRRSR = (scalarRRSR)/n
  freqRRRR = (scalarRRRR)/n
  
  # Create vector to calculate 95% CI
  haploFreq <- c(freqSSSS, freqSSSR, freqSSRR, freqSRSS, freqSRSR, freqSRRR, freqRRSS, freqRRSR, freqRRRR)
  
  # Calculate 95% Confidence Intervfor(i in 1:length(haploFreq)){
  CI_95_SSSS = 1.96 * sqrt((freqSSSS*(1-freqSSSS))/(n))
  CI_95_SSSR = 1.96 * sqrt((freqSSSR*(1-freqSSSR))/(n))
  CI_95_SSRR = 1.96 * sqrt((freqSSRR*(1-freqSSRR))/(n))
  CI_95_SRSS = 1.96 * sqrt((freqSRSS*(1-freqSRSS))/(n))
  ## Typo? SSSS in demom
  CI_95_SRSR = 1.96 * sqrt((freqSRSR*(1-freqSSSS))/(n))
  CI_95_SRRR = 1.96 * sqrt((freqSRRR*(1-freqSSSS))/(n))
  CI_95_RRSS = 1.96 * sqrt((freqRRSS*(1-freqRRSS))/(n))
  CI_95_RRSR = 1.96 * sqrt((freqRRSR*(1-freqRRSR))/(n))
  CI_95_RRRR = 1.96 * sqrt((freqRRRR*(1-freqRRRR))/(n))
  
  # Create data frame of results
  df <- data.frame(scalarSSSS, scalarSSSR, scalarSSRR
                   , scalarSRSS, scalarSRSR, scalarSRRR
                   , scalarRRSS, scalarRRSR, scalarRRRR, n)
  colnames(df) <- c("SSSS", "SSSR", "SSRR"
                    , "SRSS", "SRSR", "SRRR"
                    , "RRSS", "RRSR", "RRRR", "n")
  return(df)
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

