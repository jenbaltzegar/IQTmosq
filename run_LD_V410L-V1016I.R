# rename dataframe ---------------------------------------------
kdr <- kdrData
head(kdr)


# create functions ---------------------------------------------
# change converted columns into format readable for 'genetics' package
LDconvert1016 <- function(kdr){
  x <- ifelse(kdr$V1016I_converted == "RR", "R/R"
              , ifelse(kdr$V1016I_converted == "SR", "S/R"
                       , ifelse(kdr$V1016I_converted == "SS", "S/S"
                                , "NA")))
  return(x)
}

LDconvert1534 <- function(kdr){
  x <- ifelse(kdr$F1534C_converted == "RR", "R/R"
              , ifelse(kdr$F1534C_converted == "SR", "S/R"
                       , ifelse(kdr$F1534C_converted == "SS", "S/S"
                                , "NA")))
  return(x)
}

LDconvert410 <- function(kdr){
  x <- ifelse(kdr$V410L_converted == "RR", "R/R"
              , ifelse(kdr$V410L_converted == "SR", "S/R"
                       , ifelse(kdr$V410L_converted == "SS", "S/S"
                                , "NA")))
  return(x)
}

# Create new genotype columns in readable form for 'genetics' package --------------
kdr$V1016I_LD <- LDconvert1016(kdr)
kdr$F1534C_LD <- LDconvert1534(kdr)
kdr$V410L_LD <- LDconvert410(kdr)

head(kdr)

# Subset by year to calculate LD
# Select based on year
mosq2000 <- sqldf("Select * from kdr where newDate between '2000-01-01' and '2000-12-31'")
mosq2001 <- sqldf("Select * from kdr where newDate between '2001-01-01' and '2001-12-31'")
mosq2002 <- sqldf("Select * from kdr where newDate between '2002-01-01' and '2002-12-31'")
mosq2003 <- sqldf("Select * from kdr where newDate between '2003-01-01' and '2003-12-31'")
mosq2004 <- sqldf("Select * from kdr where newDate between '2004-01-01' and '2004-12-31'")
mosq2005 <- sqldf("Select * from kdr where newDate between '2005-01-01' and '2005-12-31'")
mosq2006 <- sqldf("Select * from kdr where newDate between '2006-01-01' and '2006-12-31'")
mosq2007 <- sqldf("Select * from kdr where newDate between '2007-01-01' and '2007-12-31'")
mosq2008 <- sqldf("Select * from kdr where newDate between '2008-01-01' and '2008-12-31'")
mosq2009 <- sqldf("Select * from kdr where newDate between '2009-01-01' and '2009-12-31'")
mosq2010 <- sqldf("Select * from kdr where newDate between '2010-01-01' and '2010-12-31'")
mosq2011 <- sqldf("Select * from kdr where newDate between '2011-01-01' and '2011-12-31'")
mosq2012 <- sqldf("Select * from kdr where newDate between '2012-01-01' and '2012-12-31'")
mosq2013 <- sqldf("Select * from kdr where newDate between '2013-01-01' and '2013-12-31'")
mosq2014 <- sqldf("Select * from kdr where newDate between '2014-01-01' and '2014-12-31'")
mosq2015 <- sqldf("Select * from kdr where newDate between '2015-01-01' and '2015-12-31'")
mosq2016 <- sqldf("Select * from kdr where newDate between '2016-01-01' and '2016-12-31'")
mosq2017 <- sqldf("Select * from kdr where newDate between '2017-01-01' and '2017-12-31'")


# Calculate pairwise LD between SNPs ------------------------------
### 2010 ###
# ID genotypes for each locus
v1016 <- genotype(mosq2010$V1016I_LD)
f1534 <- genotype(mosq2010$F1534C_LD)
v410 <-  genotype(mosq2010$V410L_LD)

ld.2010 <- LD(v1016, v410)

### 2011 ###  - 1534 is fixed
# ID genotypes for each locus
v1016 <- genotype(mosq2011$V1016I_LD)
v410 <-  genotype(mosq2011$V410L_LD)

ld.2011 <- LD(v1016, v410)


# 2012 - 1534 is fixed
# ID genotypes for each locus
v1016 <- genotype(mosq2012$V1016I_LD)
v410 <-  genotype(mosq2012$V410L_LD)
ld.2012 <- LD(v1016, v410)

### 2013 ###
v1016 <- genotype(mosq2013$V1016I_LD)
f1534 <- genotype(mosq2013$F1534C_LD)
v410 <-  genotype(mosq2013$V410L_LD)
ld.2013 <- LD(v1016, v410)

### 2014 ###
v1016 <- genotype(mosq2014$V1016I_LD)
f1534 <- genotype(mosq2014$F1534C_LD)
v410 <-  genotype(mosq2014$V410L_LD)
ld.2014 <- LD(v1016, v410)

### 2015 ### - 1534 is fixed
v1016 <- genotype(mosq2015$V1016I_LD)
v410 <-  genotype(mosq2015$V410L_LD)
ld.2015 <- LD(v1016, v410)

### 2016 ###
v1016 <- genotype(mosq2016$V1016I_LD)
v410 <-  genotype(mosq2016$V410L_LD)
ld.2016 <- LD(v1016, v410)

### 2017 ###
v1016 <- genotype(mosq2017$V1016I_LD)
f1534 <- genotype(mosq2017$F1534C_LD)
v410 <-  genotype(mosq2017$V410L_LD)
ld.2017 <- LD(v1016, v410)


# Create df with these data
year <- c("2010", "2011", "2012", "2013", "2014", "2015", "2016", "2017")
a <- c(ld.2010$`D'`, ld.2010$`R^2`, ld.2010$`X^2`, ld.2010$`P-value`, ld.2010$n)
b <- c(ld.2011$`D'`, ld.2011$`R^2`, ld.2011$`X^2`, ld.2011$`P-value`, ld.2011$n)
c <- c(ld.2012$`D'`, ld.2012$`R^2`, ld.2012$`X^2`, ld.2012$`P-value`, ld.2012$n)
d <- c(ld.2013$`D'`, ld.2013$`R^2`, ld.2013$`X^2`, ld.2013$`P-value`, ld.2013$n)
e <- c(ld.2014$`D'`, ld.2014$`R^2`, ld.2014$`X^2`, ld.2014$`P-value`, ld.2014$n)
f <- c(ld.2015$`D'`, ld.2015$`R^2`, ld.2015$`X^2`, ld.2015$`P-value`, ld.2015$n)
g <- c(ld.2016$`D'`, ld.2016$`R^2`, ld.2016$`X^2`, ld.2016$`P-value`, ld.2016$n)
h <- c(ld.2017$`D'`, ld.2017$`R^2`, ld.2017$`X^2`, ld.2017$`P-value`, ld.2017$n)

LD.V410L <- rbind(a, b, c, d, e, f, g, h)
LD.V410L <- cbind(year, LD.V410L)
colnames(LD.V410L) <- c("year", "D'", "R^2", "X^2", "P-value", "n")
LD.V410L <- as.data.frame(LD.V410L)

print(LD.V410L)
