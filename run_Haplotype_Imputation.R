# Create matrix of haplotype probabilities
# Started: 8/1/17
# depends on setup_jb.R and meltcurve_analysis.R

# This matrix assumes that RS does not exist in the population
# Locus 1016 is always written first

### Create matrix of haplotype probabilities ------------------------------------------------------------
genos <- c("SSSS", "SSSR", "SSRR", "SRSS", "SRSR", "SRRR", "RRSS", "RRSR", "RRRR")
haps <- c("SS", "SR", "RS", "RR")
##? how was this computed?
HapProbs <- matrix(data = c(1,0,0,0,0.5,0.5,0,0,0,1,0,0,1,0,0,0,0.33,0.33,0,0.33,0,0.5,0,0.5,0,0,0,0,0,0,0,1,0,0,0,1)
       , nrow = 9, ncol = 4, byrow = T, dimnames = list(genos, haps))

### Impute haplotype numbers ------------------------------------------------------------
# For each year, multiply the genotype number by the haplotype probability

lfreq <- lapply(1:18, function(yr) {
    dat <- mc.haps.yr[yr,]
    ## initialize empty, add year as rowname
    df <- matrix(c(0,0,0,0), nrow=1, ncol=4, dimnames=list(dat$year, haps))
    for(i in 1:9){
      df <- df + dat[[1+i]] * HapProbs[i,]
    }
    result <- with(dat, 
        cbind(df/n, year=as.numeric(year), n=n)
    )
})

## combine into data.frame
freqAll_long <- (
    do.call(rbind, lfreq)
    %.>% as.data.frame(.)
    %.>% melt(., id.vars=c("year","n"), variable.name = "Haplotype", value.name = "Frequency")
    %.>% mutate(., CI_95 = fun.ci(Frequency, n))
)


#write.csv(freqAll_long, 'data/freqAll_long.csv', quote=F, row.names=F)

# # Subset data based on haplotype
# freqSS <- freqAll_long[freqAll_long$Haplotype=="SS", ]
# freqSR <- freqAll_long[freqAll_long$Haplotype=="SR", ]
# freqRS <- freqAll_long[freqAll_long$Haplotype=="RS", ]
# freqRR <- freqAll_long[freqAll_long$Haplotype=="RR", ]
