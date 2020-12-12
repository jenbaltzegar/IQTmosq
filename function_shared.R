## compute confidence interval
fun.ci <- function(freq, n, z=1.96)
  z * sqrt((freq*(1-freq))/(2*n))
