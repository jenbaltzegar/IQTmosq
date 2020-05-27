
# source this file to load all packages and functions needed
# in model_runfile.R

######################################################## 
# may be able to reduce dependencies
library(reshape2)
library(dplyr)
library(doParallel)
library(ggplot2)
library(gridExtra)
library(ggthemes)
library(viridis)
library(RColorBrewer)

##ggplot theme
my_theme <- function(){
  theme_foundation(base_size=14) + 
    theme(
      plot.title = element_text(face = "bold",size = rel(1.2), hjust = 0.5),
      panel.background = element_rect(colour = NA),
      plot.background = element_rect(colour = NA),
      panel.border = element_rect(colour = NA),
      panel.spacing=unit(1, "lines"),
      axis.title = element_text(face = "bold"),
      axis.title.y = element_text(vjust = 2),
      axis.title.x = element_text(vjust = -0.2),
      axis.line = element_line(colour="black"),
      panel.grid.major = element_line(colour="#f0f0f0"),
      panel.grid.minor = element_blank(),
      plot.margin=unit(c(10,5,5,5),"mm"),
      strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
      strip.text = element_text(face="bold")
    )
}

## set up mathematical model and functions
alleles <- list(c("R","S"))
genotypes <- list(c("R","R"),c("S","R"),c("S","S"))

#turn the list into a vector of strings
genotype_names <- do.call(paste,c(as.data.frame(t(as.data.frame(genotypes))),sep=""))
ngenos <- length(genotypes)

WF_1gen <- function(theta,genos_t,stoch_on=0){
  #alternatively, allele frequencies p^2, 2pq, q^2
  nextgen <- c(
    (genos_t[1]^2 + genos_t[1]*genos_t[2] + 1/4*genos_t[2]^2),
    (1-theta["h"]*theta["sSS"])*(genos_t[1]*genos_t[2] + 2*genos_t[1]*genos_t[3] +
                 1/2*genos_t[2]^2 + genos_t[2]*genos_t[3]),
    (1-theta["sSS"])*(genos_t[3]^2 + genos_t[2]*genos_t[3] + 1/4*genos_t[2]^2)
  )
  names(nextgen) <- genotype_names
  nextgen <- nextgen/sum(nextgen)
  if (stoch_on==1){
    nextgen <- rmultinom(1,theta["popsize"],nextgen)/theta["popsize"]
  }
  nextgen
}

WF_sim <- function(theta,n_gens=12,stoch_on=0){
  # repeat call to WF_1gen to create full simulation
  # alternative parameterization could set initial RR and SR frequencies separately
  # allgens[1,] <- c(theta["RR_0"],theta["SR_0"],1-theta["RR_0"]-theta["SR_0"])
  
  allgens <- matrix(0,ncol=3,nrow=n_gens)
  p1 <- theta["r0"]
  allgens[1,] <- c(p1^2,2*p1*(1-p1),(1-p1)^2) #hardy-weinberg
  for (i in 2:n_gens){
    nextgen <- WF_1gen(theta,genos_t=allgens[i-1,],stoch_on=stoch_on)
    allgens[i,] <- nextgen
  }
  allgens <- as.data.frame(allgens)
  names(allgens) <- genotype_names
  allgens
}

WF_sampler <- function(stoch_dat,Ns){
  # stoch_dat: simulation, Ns:vector of sampling counts
  # produces multinomial samples at each generation of sizes Ns
  samples <- 0 * stoch_dat
  
  for(i in 1:nrow(stoch_dat)){
    if(Ns[i]>0 & sum(stoch_dat[i,genotype_names])>0){
      samples[i,genotype_names] <- rmultinom(1,size=Ns[i],
                                             prob=stoch_dat[i,genotype_names])
    } else{
      samples[i,genotype_names] <- rep(0,3)
    }
  }
  samples
}

## set up parameter estimation functions
costfunc <- function(simdat,dat){
  # multinomial log likelihood
  simdat[dat==0 & simdat == 0] <- 1 #0*log(0)=NaN, but 0*log(1) = 0 gives proper log value
  LL <- dat*log(simdat)
  -sum(LL)
}

val_func <- function(p,parnames,theta,dat){
  # minimizer function
  p <- transfrom(p)
  theta[parnames] <- p
  simdat <- WF_sim(theta=theta,stoch_on=0,n_gens=nrow(dat))
  if(any(simdat<0)){
    return(1e7)
  }
  cost <- costfunc(simdat,dat) #calculates -LL
  
  cost
}

## transform parameters for unconstrained optimization from [0,1]
transto <- function(p) {(log(p/(1-p)))}
transfrom <- function(p) {1/(1+exp(-p))}

get_MLE <- function(dat_in,theta,parnames){
  # given data dat_in, produces MLE for parameters in parnames
  # other parameters set as in theta.
  # typically smooth surface, but need to be careful about local minima
  tend <- nrow(dat_in)
  p <- rep(0.5,length(parnames))
  sol <- optim(par=transto(p),fn=val_func,parnames=parnames,theta=theta,dat=dat_in,control=list(maxit=2000))
  if(sol$convergence!=0){
    print("Solver did not converge.")
  }
  names(sol$par) <- parnames
  paramfit <- transfrom(sol$par)
  theta[parnames] <- paramfit
  theta["LL"] <- -sol$value
  theta
}
