library(ggplot2)
library(magrittr)
library(reshape2)

# set working directory
setwd("~/Dropbox/GouldLab/Project_Mosquito/Database")
######################################################## 
## set up mathematical model and functions
alleles <- list(c("R","S"))
genotypes <- list(c("R","R"),c("S","R"),c("S","S"))

#turn the list into a vector of strings
genotype_names <- do.call(paste,c(as.data.frame(t(as.data.frame(genotypes))),sep=""))

#setequal is unordered matching.
position_finder <- function(y) Position(function(x) setequal(y,x),x=genotypes)

combos <- list(c(1,1),c(1,2),c(2,1),c(2,2)) #need to do each of these combinations for each locus, and combine all the combinations
ncombos <- length(combos)
ngenos <- length(genotypes)
##
allweights <- seq_len(ngenos) %>% lapply(function(m){
  seq_len(ngenos) %>% lapply(function(n){
    weights <- rep(0,ngenos)
    positions <- lapply(combos,function(y) c(genotypes[[m]][y[1]],genotypes[[n]][y[2]])) %>%
      sapply(position_finder)
    for (pos in 1:ncombos){
      weights[positions[pos]] <- weights[positions[pos]] + 1/ncombos
    }
    weights
  })
})

generate_offprob <- function(){
  off_probs <- array(0,dim=c(rep(ngenos,3)))
  for(m in 1:ngenos){
    for(n in 1:ngenos){
      weights <- allweights[[m]][[n]]
      off_probs[m,n,] <- weights
    }
  }
  nonzero_mate_calc(off_probs)
}

nonzero_mate_calc <- function(off_probs){
  nonzeros <- data.frame(which(off_probs!=0,arr.ind = T))
  names(nonzeros) <- c("fem","male","offs")
  nonzeros$prob <- apply(nonzeros,1,function(x) off_probs[x[1],x[2],x[3]]) #new column with P(ijk)
  mates_list <- lapply(seq(1,ngenos),function(x) nonzeros %>% subset(offs==x))
  lapply(mates_list,as.matrix)
}
## calculate genotype frequencies and normalize
freq_calc <- function(Fem,Male,mates_list){
  freqs <- sapply(mates_list,function(y)
    sum(Fem[y[,1]]*Male[y[,2]]*y[,4])
  )
  freqs
}
## make vectors that give the counts of the H, R, and B alleles for each genotype. to assess fitness
allele_counter <- function(genotypes,allele){genotypes %>% sapply(function(y) sum(y==allele))}
allele_counts <- list("R"="R","S"="S") %>% lapply(function(x) allele_counter(genotypes,x))

## multiplicative fitness costs, constant is automatically vectorized when taking to the power of vector
fitness_assesser <- function(freqs,sSS,sSR,costS=0,type=0){
  if (type == 1){
    newfreqs <- freqs * (1-costS)^allele_counts$S
  }else{
    newfreqs <- freqs
    newfreqs[2] <- freqs[2] * (1-sSR)
    newfreqs[3] <- freqs[3] * (1-sSS)
  }
  newfreqs <- newfreqs / sum(newfreqs)
  newfreqs
}

off_probs=generate_offprob()

full_sim <- function(sSS=0,sSR=0,initRR=0.01,initSR=0.01,off_probs,costS=0,costtype=0){
  
  out <- as.data.frame(matrix(nrow=tend,ncol=3))
  colnames(out) <- genotype_names
  out[1,] <- c(initRR,initSR,1-(initRR+initSR))
  
  for (t in 2:(tend)){
    simM <- array(0,dim=c(12,ngenos))
    simM[1,] <- c(initRR,initSR,1-(initRR+initSR))
    simF <- simM
    for (i in 2:12){
      simFtemp <- simF[i-1,]#fitness_assesser(simF[t-1,],sS=sS,h=h)
      simMtemp <- simM[i-1,]#fitness_assesser(simM[t-1,],sS=sS,h=h)
      nextgen <- freq_calc(simFtemp,simMtemp,off_probs)
      nextgen <- fitness_assesser(nextgen,sSS=sSS,sSR=sSR,costS=costS,type=costtype)
      simF[i,] <- nextgen
      simM[i,] <- nextgen
      
    }
    out[t,] <- colMeans((simM+simF)/2)
    initRR <- out[t,"RR"]
    initSR <- out[t,"SR"]
  }
  out
}

########################################################
## set up parameter estimation functions
## least squares function. or multinomial log likelihood
costfunc <- function(simdat,dat){
  diffs <- (simdat-dat)
  sum(diffs^2)
  # simdat[simdat == 0] <- 1e-9
  # LL <- dat*log(simdat)
  # -sum(LL)
}

## minimizer function
val_func <- function(p,pars,dat){
  p <- transfrom(p)
  pars[parnames] <- p
  simdat <- full_sim(sSS = pars["sSS"],sSR = pars["sSR"],initRR=pars["initRR"],initSR = pars["initSR"],
                     off_probs = off_probs)
  cost <- costfunc(simdat,dat) #calculate least squares
  cost
}

##transform parameters for unconstrained optimization from [0,1]
transto <- function(p) (log(p/(1-p)))
transfrom <- function(p) 1/(1+exp(-p))

########################################################
## estimate parameters 

## read in data from current working directory
dat <- read.csv('mc.1534.yr_reduced.csv')
# dat <- read.csv('mc.1016.yr_reduced.csv')
dat <- dat %>% subset(year <= 2005)
freq_dat <- dat
freq_dat[,c("RR","SR","SS")] <- freq_dat[,c("RR","SR","SS")]/rowSums(freq_dat[,c("RR","SR","SS")])

tend <- nrow(dat)
# choose parameters to estimate. fitness costs, initial proportions of RR and SR
# set parameters to appropriate constant if not in parnames (not going to estimate)

#costtype 0 for separate costs for each genotype: estimate costs sSS and sSR
pars <- c("sSS"=0,"sSR"=0,"initRR"=0,"initSR"=0)
pars[c(3,4)] <- unlist(freq_dat[1,c("RR","SR")])
parnames <- c("sSS","sSR")#,"initRR","initSR")
estim_dat <- freq_dat %>% subset(select=c(RR,SR,SS))

p <- rep(0.3,length(parnames)) # initial parameter values

## parameter estimation
sol <- optim(par=transto(p),fn=val_func,pars=pars,dat=estim_dat,control=list(maxit=2000))
c(sol$convergence,sol$counts)
names(sol$par) <- names(pars[parnames])
paramfit <- transfrom(sol$par)
paramfit

## run simulation with parameters
pars[parnames] <- paramfit
simdat <- full_sim(sSS = pars["sSS"],sSR = pars["sSR"],initRR=pars["initRR"],initSR = pars["initSR"],off_probs = off_probs)
simdat$year <- dat$year

## plot simulation with data
freq_dat %>% subset(select=c(year,RR,SR,SS)) %>%
  melt("year",variable.name="genotype",value.name="frequency") %>%
  ggplot(aes(x=year,y=frequency,color=genotype)) +
  geom_point() +
  geom_line(data=simdat %>% melt("year",variable.name="genotype",value.name="frequency")) +
  ggtitle("Simulated Data") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Year") +
  ylab("Frequency")


####### fitness cost fit grid
parvals <- expand.grid(sSS=seq(0,1,.02),sSR=seq(0,1,.02))#,initRR=c(1e-5,1e-4,1e-3,1e-2),initSR=c(1e-3,1e-2))
parnames <- names(parvals)
parvals$fit <- -10

for (i in 1:nrow(parvals)){
  # p <- parvals[i,1:(ncol(parvals)-1)]
  pars[parnames] <- unlist(parvals[i,1:(ncol(parvals)-1)])
  simdat <- full_sim(sSS = pars["sSS"],sSR = pars["sSR"],initRR=pars["initRR"],initSR = pars["initSR"],
                     off_probs = off_probs)
  parvals[i,"fit"] <- costfunc(simdat,estim_dat) #calculate least squares
}

parvals %>% ggplot() +
  geom_tile(aes(x=sSS,y=sSR,fill=log(fit))) +
  ggtitle("Fitness Cost Fit Grid") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("Cost SS") +
  ylab("Cost SR")
