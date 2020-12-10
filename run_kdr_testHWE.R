# This script will parse kdr data and analyze for HWE


# rename obj to avoid overwriting
kdr <- kdrData
kdr$newDate <- as.Date(kdr$newDate)

# remove rows where newDate is NA
kdr <- kdr[complete.cases(kdr$newDate),]

# restrict df from Oct 2002 - Dec 2005 
kdr <- kdr[kdr$newDate >= "2002-10-01" & kdr$newDate <= "2005-12-31",]
# # remove unnecessary columns
# kdr <- kdr[, c(1:3,11,14,25:28)]

# split kdr into month+year groups
kdr.moYr <- split(kdr, format(kdr$newDate, "%Y-%m"))

# calculate genotype frequencies for each Mo-Year
genos <- lapply(kdr.moYr, mc.1534)


# Write function to generate genotype object for HWE.test
make.df <- function(df){
  g1   <- c(rep("R/R",df$RR),
            rep("R/S",df$SR),
            rep("S/S",df$SS))
  g1<-genotype(g1)
  return(g1)
}

# generate objects
genos <- lapply(genos, make.df)

# calculate HWE for each month+year group
oct.02 <- HWE.test(genos$`2002-10`)
nov.02 <- HWE.test(genos$`2002-11`)
dec.02 <- HWE.test(genos$`2002-12`)

jan.03 <- HWE.test(genos$`2003-01`)
feb.03 <- HWE.test(genos$`2003-02`)
mar.03 <- HWE.test(genos$`2003-03`)
apr.03 <- HWE.test(genos$`2003-04`)
may.03 <- HWE.test(genos$`2003-05`)
jun.03 <- HWE.test(genos$`2003-06`)
jul.03 <- HWE.test(genos$`2003-07`)
aug.03 <- HWE.test(genos$`2003-08`)
sep.03 <- HWE.test(genos$`2003-09`)
oct.03 <- HWE.test(genos$`2003-10`)
nov.03 <- HWE.test(genos$`2003-11`)
dec.03 <- HWE.test(genos$`2003-12`)

jan.04 <- HWE.test(genos$`2004-01`)
feb.04 <- HWE.test(genos$`2004-02`)
mar.04 <- HWE.test(genos$`2004-03`)
apr.04 <- HWE.test(genos$`2004-04`)
may.04 <- HWE.test(genos$`2004-05`)
jun.04 <- HWE.test(genos$`2004-06`)
jul.04 <- HWE.test(genos$`2004-07`)
aug.04 <- HWE.test(genos$`2004-08`)
sep.04 <- HWE.test(genos$`2004-09`)
oct.04 <- HWE.test(genos$`2004-10`)
nov.04 <- HWE.test(genos$`2004-11`)
dec.04 <- HWE.test(genos$`2004-12`)

jan.05 <- HWE.test(genos$`2005-01`)
feb.05 <- HWE.test(genos$`2005-02`)
mar.05 <- HWE.test(genos$`2005-03`)
apr.05 <- HWE.test(genos$`2005-04`)
may.05 <- HWE.test(genos$`2005-05`)
jun.05 <- HWE.test(genos$`2005-06`)
jul.05 <- HWE.test(genos$`2005-07`)
aug.05 <- HWE.test(genos$`2005-08`)
sep.05 <- HWE.test(genos$`2005-09`)
oct.05 <- HWE.test(genos$`2005-10`)
nov.05 <- HWE.test(genos$`2005-11`)
dec.05 <- HWE.test(genos$`2005-12`)

hwe.pvals <- c(0.1579, 1, 1, #2002
               0.02098, 0.3412, 1.562e-05, 1, 1, 0.4, 1, 0.2967, 0.2168, 0.3821, 0.08727, 0.05952, #2003
               0.2405, 0.3647, 0.1994, 0.1327, 0.007385, 0.2348 #2004
               )

months <- c("oct.02", "nov.02", "dec.02",
            "jan.03", "feb.03", "mar.03", "apr.03", "may.03", "jun.03", "jul.03", "aug.03", "sep.03", "oct.03", "nov.03", "dec.03",
            "jan.04", "feb.04", "mar.04", "may.04", "jun.04", "jul.04")

kdr.pvals <- as.data.frame(cbind(months, hwe.pvals))
print(kdr.pvals)





