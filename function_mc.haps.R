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
  
    # Calculate 95% Confidence Intervfor(i in 1:length(haploFreq)){
    result <- within(list(), {
        SSSS <- fun.ci(freqSSSS,n)
        SSSR <- fun.ci(freqSSSR,n)
        SSRR <- fun.ci(freqSSRR,n)
        SRSS <- fun.ci(freqSRSS,n)
        SRSR <- fun.ci(freqSRSR,n)
        SRRR <- fun.ci(freqSRRR,n)
        RRSS <- fun.ci(freqRRSS,n)
        RRSR <- fun.ci(freqRRSR,n)
        RRRR <- fun.ci(freqRRRR,n)
        n <- n
    })

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

