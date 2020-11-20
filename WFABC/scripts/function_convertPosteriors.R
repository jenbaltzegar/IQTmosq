### Convert WFABC posteriors to match Michael Vella's scale

# make a function to convert posteriors for dataframes with two loci
convertPosteriors.twolocus <- function(posteriors.s, posteriors.h){
  # load data
  s <- read.table(posteriors.s)
  h <- read.table(posteriors.h)
  
  # make matrix for locus 1016
  sh.1016 <- cbind(t(s[1,]), t(h[1,]))
  colnames(sh.1016) <- c("s", "h")
  sh.1016 <- as.matrix(sh.1016)
  
  # make matrix for locus 1534
  sh.1534 <- cbind(t(s[2,]), t(h[2,]))
  colnames(sh.1534) <- c("s", "h")
  sh.1534 <- as.matrix(sh.1534)

  # convert posteriors for 1016
  RR <- 1 + sh.1016[,1]
  SR <- 1 + (sh.1016[,1] * sh.1016[,2])
  SS <- 1
  RR.prime <- 1
  SR.prime <- SR/RR
  SS.prime <- SS/RR
  
  cost.SS <- 1 - SS.prime
  cost.SR <- 1 - SR.prime
  
  s.prime.1016 <- cost.SS
  h.prime.1016 <- cost.SR/cost.SS

  # convert posteriors for 1534
  RR <- 1 + sh.1534[,1]
  SR <- 1 + (sh.1534[,1] * sh.1534[,2])
  SS <- 1
  RR.prime <- 1
  SR.prime <- SR/RR
  SS.prime <- SS/RR
  
  cost.SS <- 1 - SS.prime
  cost.SR <- 1 - SR.prime
  
  s.prime.1534 <- cost.SS
  h.prime.1534 <- cost.SR/cost.SS
  
  # make new dfs
  converted.s <- rbind(s.prime.1016, s.prime.1534)
  converted.h <- rbind(h.prime.1016, h.prime.1534)
  write.table(converted.s, file = paste0("converted_", posteriors.s), row.names = FALSE, col.names = FALSE)
  write.table(converted.h, file = paste0("converted_", posteriors.h), row.names = FALSE, col.names = FALSE)
}

# make a function to convert posteriors for dataframes with one locus
convertPosteriors.onelocus <- function(posteriors.s, posteriors.h){
  # load data
  s <- read.table(posteriors.s)
  h <- read.table(posteriors.h)
  
  # make matrix for locus 
  sh <- cbind(t(s[1,]), t(h[1,]))
  colnames(sh) <- c("s", "h")
  sh <- as.matrix(sh)
  
  # convert posteriors 
  RR <- 1 + sh[,1]
  SR <- 1 + (sh[,1] * sh[,2])
  SS <- 1
  RR.prime <- 1
  SR.prime <- SR/RR
  SS.prime <- SS/RR
  
  cost.SS <- 1 - SS.prime
  cost.SR <- 1 - SR.prime
  
  s.prime <- cost.SS
  h.prime <- cost.SR/cost.SS
  
  # make new dfs
  write.table(s.prime, file = paste0("converted_", posteriors.s), row.names = FALSE, col.names = FALSE)
  write.table(h.prime, file = paste0("converted_", posteriors.h), row.names = FALSE, col.names = FALSE)
}

