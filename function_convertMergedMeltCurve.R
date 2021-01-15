### These functions are used to convert the melt curve locus data into genotypes

### For the 1016 Locus
# HomozygousLow melt curve results are homozygous for the resistance allele
# Heterozygous melt curve results are heterozygous for the resistance and susceptible allele
# HomozygousHigh melt curve results are homozygous for the susceptible allele
# All else NA
convert1016 <- function(x){
    x <- factor(x,
        levels=c("HomozygousLow", "Heterozygous", "HomozygousHigh", "error"), 
        labels=c("RR", "SR", "SS", "error")
    )
    x[is.na(x)] <- 'error'
    return(x)
}

### For the 1534 Locus
# HomozygousLow melt curve results are homozygous for the susceptible allele
# Heterozygous melt curve results are heterozygous for the resistance and susceptible allele
# HomozygousHigh melt curve results are homozygous for the resistance allele

convert1534 <- function(x){
    x <- factor(x,
        levels=c("HomozygousLow", "Heterozygous", "HomozygousHigh", "error"), 
        labels=c("SS", "SR", "RR", "error")
    )
    x[is.na(x)] <- 'error'
    return(x)
}


### For the 410 Locus
# HomozygousLow melt curve results are homozygous for the resistance allele
# Heterozygous melt curve results are heterozygous for the resistance and susceptible allele
# HomozygousHigh melt curve results are homozygous for the susceptible allele

convert410 <- function(x){
    x <- factor(x,
        levels=c("HomozygousLow", "Heterozygous", "HomozygousHigh", "error"), 
        labels=c("RR", "SR", "SS", "error")
    )
    x[is.na(x)] <- 'error'
    return(x)
}
