### These functions either create dataframes for WFABC analysis or remove timepoints with no samples

# make a function to create df for one or two loci by year
create.WFABC.df <- function(df1, df2 = NULL, n_loci, filename, gens_per_yr = 12) {
  if(!is.null(df2)) df1 <- df1 + df2
  
  # convert certain columns to numeric
  df1$RR <- as.numeric(as.character(df1$RR))
  df1$SR <- as.numeric(as.character(df1$SR))
  df2$RR <- as.numeric(as.character(df2$RR))
  df2$SR <- as.numeric(as.character(df2$SR))
  
  # calculate variables
  n_timepts <- nrow(df1)
  list_timepts <- c(1:n_timepts)
  gens_timepts <- gens_per_yr*list_timepts
  n_chr_timepts <- 2*df1$n
  n_R_timepts <- ((2*df1$RR) + df1$SR)
  n_chr_timepts2 <- 2*df2$n
  n_R_timepts2 <- ((2*df2$RR) + df2$SR)
  
  # write txt file
  sink(filename)
  cat(n_loci, n_timepts, "\n")
  cat(gens_timepts, "\n")
  cat(n_chr_timepts, "\n")
  cat(n_R_timepts, "\n")
  cat(n_chr_timepts2, "\n")
  cat(n_R_timepts2)
  sink()
}

# make a function to create df for one locus by month
create.WFABC.df.byMo <- function(df1, filename, gens_per_yr = 12) {

  # calculate variables
  n_timepts <- nrow(df1)
  list_timepts <- c(1:n_timepts)
  n_chr_timepts <- 2*df1$n
  n_R_timepts <- ((2*df1$RR) + df1$SR)

  # write txt file
  sink(filename)
  cat("1", n_timepts, "\n")
  cat(list_timepts, "\n")
  cat(n_chr_timepts, "\n")
  cat(n_R_timepts, "\n")
  sink()
}

# make a function to remove months where no individuals were sampled - two locus
remove.noSamples.twolocus <- function(filename){
  df <- read.csv(file = filename, header = FALSE, sep = "")
  df <- t(df) # transpose df
  df <- t(df[df[,3] != 0,])
  
  # write txt file
  sink(filename)
  cat("2", ncol(df), "\n")
  cat(df[2,], "\n")
  cat(df[3,], "\n")
  cat(df[4,], "\n")
  cat(df[5,], "\n")
  cat(df[6,], "\n")
  sink()
}

# make a function to remove months where no individuals were sampled - one locus
remove.noSamples.onelocus <- function(filename){
  df <- read.csv(file = filename, header = FALSE, sep = "")
  df <- t(df) # transpose df
  df <- t(df[df[,3] != 0,])
  
  # write txt file
  sink(filename)
  cat("1", ncol(df), "\n")
  cat(df[2,], "\n")
  cat(df[3,], "\n")
  cat(df[4,], "\n")
  sink()
}
