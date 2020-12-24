# This script will parse kdr data and analyze for HWE

# create new obj to avoid overwriting
kdr <- kdrData
kdr$newDate <- as.Date(kdr$newDate)

# remove rows where newDate is NA
kdr <- subset(kdr, 
    # restrict df from Oct 2002 - Dec 2005 
    ##? Jul 2014
    !is.na(newDate)
    & newDate >= "2002-10-01" 
    & newDate <= "2004-07-31",
    select = c(newDate, F1534C_converted) 
)

# split kdr into month+year groups
kdr.yrmo <- split(kdr, format(kdr$newDate, "%Y-%m"))

# Helper function 
## generate genotype object, run HWE.test
make.hwe <- function(df){
  g1   <- c(rep("R/R",df$RR),
            rep("R/S",df$SR),
            rep("S/S",df$SS))
  g1<-genotype(g1)
  ## test gives error at fixation
  if (df$freqR == 1) return()
  result <- HWE.test(g1)
  return(result)
}

# calculate genotype frequencies for each Mo-Year
genos <- lapply(kdr.yrmo, mc.1534)
# calculate HWE for each month+year group
hwe.list <- lapply(genos, make.hwe)
## inspect individual results:
# hwe.list[[1]]
## inspect contents of test:
# str(dat.hwe[[1]])

## extract pvalues into data.frame 
hwe.pvals <- ldply(hwe.list, function(x) c(p.val=x$test$p.value), .id='YrMo')

## old code
if (F) {
    hwe.pvals <- c(0.1579, 1, 1, #2002
                   0.02098, 0.3412, 1.562e-05, 1, 1, 0.4, 1, 0.2967, 0.2168, 0.3821, 0.08727, 0.05952, #2003
                   0.2405, 0.3647, 0.1994, 0.1327, 0.007385, 0.2348 #2004
                   )

    months <- c("oct.02", "nov.02", "dec.02",
                "jan.03", "feb.03", "mar.03", "apr.03", "may.03", "jun.03", "jul.03", "aug.03", "sep.03", "oct.03", "nov.03", "dec.03",
                "jan.04", "feb.04", "mar.04", "may.04", "jun.04", "jul.04")

    kdr.pvals <- as.data.frame(cbind(months, hwe.pvals))
    print(kdr.pvals)
}
