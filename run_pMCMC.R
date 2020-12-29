
library(nimble) #version 0.9.1
library(coda)
library(doParallel)
library(doRNG)

source('model_setup.R')
## read in data for locus 1534 from current working directory
locus <- '1534'
dat <- read.csv(paste0('mc.',locus,'.MoYr_reduced_byMonth.csv'))
dat$gen <- seq_len(nrow(dat))

start_gen <- 34
end_gen <- 12*10
tend <- end_gen - start_gen + 1

#subset data by specified generations 
estim_dat <- dat[dat$gen >= start_gen & dat$gen <= end_gen,genotype_names] #genotype columns only

# choose parameters to estimate. fitness parameters, initial R allele frequency
# set parameters to appropriate constant if not in parnames (not going to estimate)
# popsize is only used in stochastic simulations
my_popsize <- 500
theta1 <- c("sSS"=0.3,"h"=0.05,"r0"=0.01,"popsize"=my_popsize)
parnames <- c("sSS","h","r0")
theta1["r0"] <- dat$freqR[start_gen]

# maximum likelihood estimate from function in model_setup.R
fit <- get_MLE(estim_dat,theta1,parnames)
fit
theta_fit_dat <- fit[names(theta1)]

Nsamples <- rowSums(estim_dat)

##########
## Run pMCMC chains in parallel. Time to run is several hours (depends on machine).
## With nimble, need to build the model separately for every core.
## If using a machine with fewer cores, can instead use runMCMC(nchains=3) 

dir.create("pMCMC_chains")
cores=3
myclust <- makeCluster(cores,outfile="pMCMC_chains/logs.txt")
registerDoParallel(myclust)
myseed=12345
registerDoRNG(myseed)

#set up log file
#t.str <- Sys.time()#strptime(Sys.time(), "%Y-%m-%d %H:%M:%S")
#tfile <- paste0(as.numeric(format(t.str, "%H")),"_",as.numeric(format(t.str, "%M")),"_",as.numeric(format(t.str, "%S")))
#logfile <- paste0("pMCMC_chains/log_",tfile,".txt")
#writeLines(c(""),logfile)

my_niter <- 60000

all_runs <- data.frame(iter=c(1,2,3))
nsims <- 3
numberruns <- 3

allfits <- foreach(x=iter(all_runs,by="row"),
                   .combine=rbind,
                   .packages=c("coda","nimble")) %dopar% {
                     # define Dirichlet-multinomial distribution (given as an example in user manual)
                     #deregisterDistributions('ddirchmulti')
                     ddirchmulti <- nimbleFunction(
                       run = function(x = double(1), Alpha = double(1), size = double(0),
                                      log = integer(0, default = 0)) {
                         returnType(double(0))
                         
                         Alpha <- Alpha + 0.0001
                         
                         logProb <- lgamma(size+1) - sum(lgamma(x+1)) + lgamma(sum(Alpha)) -
                           sum(lgamma(Alpha)) + sum(lgamma(Alpha + x)) - lgamma(sum(Alpha) + size)
                         
                         if(log) return(logProb)
                         else return(exp(logProb))
                       })
                     assign('ddirchmulti', ddirchmulti, envir = .GlobalEnv)
                     rdirchmulti <- nimbleFunction(
                       run = function(n = integer(0), Alpha = double(1), size = double(0)) {
                         returnType(double(1))
                         if(n != 1) print("rdirchmulti only allows n = 1; using n = 1.")
                         p <- rdirch(1, Alpha+0.0001)
                         return(rmulti(1, size = size, prob = p))
                       })
                     assign('rdirchmulti', rdirchmulti, envir = .GlobalEnv)
                     
                     # define model
                     nim_mod <- nimbleCode({
                       Xtemp[1,1] <- r0^2
                       Xtemp[1,2] <- 2*r0*(1-r0)
                       Xtemp[1,3] <- (1-r0)^2
                       Xnorm[1,1:3] <- X[1,1:3]/sum(X[1,1:3])
                       X[1,1:3] ~ dmulti(Xtemp[1,1:3],popsize)
                       
                       Y[1,1:3] ~ ddirchmulti(Alpha = A*Xnorm[1,1:3],size=Ns[1])
                       
                       for(i in 2:tend){
                         Xtemp[i,1] <- Xnorm[i-1,1]^2 + Xnorm[i-1,1]*Xnorm[i-1,2] + 1/4*Xnorm[i-1,2]^2
                         Xtemp[i,2] <- (1-h*sSS)*(Xnorm[i-1,1]*Xnorm[i-1,2] + 2*Xnorm[i-1,1]*Xnorm[i-1,3] +
                                                    1/2*Xnorm[i-1,2]^2 + Xnorm[i-1,2]*Xnorm[i-1,3])
                         Xtemp[i,3] <- (1-sSS)*(Xnorm[i-1,3]^2 + Xnorm[i-1,2]*Xnorm[i-1,3] + 1/4*Xnorm[i-1,2]^2)
                         
                         Xnorm[i,1:3] <- X[i,1:3]/sum(X[i,1:3])
                         
                         X[i,1:3] ~ dmulti(Xtemp[i,1:3],popsize)
                         
                         Y[i,1:3]  ~ ddirchmulti(Alpha = A*Xnorm[i,1:3],size = Ns[i])
                       }
                       
                       # Uninformative priors
                       # A ~ dexp(0.5)
                       A ~ dgamma(0.01,0.01)
                       r0 ~ dbeta(1,1)
                       h ~ dbeta(1,1)
                       sSS ~ dbeta(1,1)
                     })
                     
                     full_model <- nimbleModel(
                       nim_mod,
                       data=list(Y=estim_dat),
                       constants=list(tend=tend,
                                      popsize=theta_fit_dat["popsize"],
                                      Ns=Nsamples
                       )
                     )
                     
                     initsFunction <- function() list(
                       A = 5,
                       sSS = theta_fit_dat["sSS"],
                       r0 = theta_fit_dat["r0"],
                       h = theta_fit_dat["h"]
                     )
                     
                     mcmcConf <- configureMCMC(full_model, nodes = NULL)
                     ## block pMCMC sampler
                     mcmcConf$addSampler(target = c(parnames,'A'),
                                         type = 'RW_PF_block', 
                                         control <- list(
                                           pfType = "bootstrap", 
                                           pfControl = list(saveAll = T, smoothing = T),
                                           pfNparticles =1000,
                                           latents = c("X")))
                     
                     mcmcConf$addMonitors2(c("X"))
                     ## build and compile pMCMC sampler
                     Rmcmc <- buildMCMC(mcmcConf)
                     Cmcmc <- compileNimble(full_model, Rmcmc, resetFunctions = TRUE)

                     samps <- runMCMC(Cmcmc$Rmcmc,
                                      niter=my_niter,
                                      # nburnin=30000,
                                      nchains=1,
                                      inits = initsFunction,
                                      samplesAsCodaMCMC=TRUE)
                     
                     saveRDS(list(samp1=samps$samples,samp2=samps$samples2,
                                  dat=estim_dat,Nsamples=Nsamples,my_popsize=my_popsize,
                                  theta_fit_dat=theta_fit_dat),
                             paste0(getwd(),"/pMCMC_chains/samps",x))
                     return(NULL)
                   }

stopCluster(myclust)


#################### Produce summary statistics and plots
datfolder <- "pMCMC_chains"
flist <- as.list(list.files(paste0(getwd(),"/",datfolder),pattern="samp"))
mcmc_out <- readRDS(paste0(getwd(),"/",datfolder,"/",flist[[1]]))
allchains <- flist %>% lapply(function(fname) {
  mcmc_out <- readRDS(paste0(getwd(),"/",datfolder,"/",fname))
  mcmc_out$samp1
}) %>% mcmc.list()
allchainsX <- flist %>% lapply(function(fname) {
  mcmc_out <- readRDS(paste0(getwd(),"/",datfolder,"/",fname))
  mcmc_out$samp2
}) %>% mcmc.list()
varnames(allchains) <- c("A","h","R0","s")

## remove burn-in
chains <- window(allchains,start=30000)

## plot of chains
# plot(allchains)
# plot(allchains,density=F)
plot(chains)
#plot(chains,density=F)

## diagnostics and parameter estimates
summary(chains)
effectiveSize(chains)
crosscorr(chains)
autocorr.diag(chains)
gelman.diag(chains)

## use thinned chain for plotting
chains <- window(allchains,start=30000,thin=20)
chains_matrix <- as.matrix(chains)
chains_df <- as.data.frame(as.matrix(chains))
ggpairs(chains_df[seq(1,nrow(chains_df),100),]) +my_theme


#################### Plot 95% CI of states

theta_fit <- mcmc_out$theta_fit
yvars <- genotype_names

#get state means and intervals
chainsX <- window(allchainsX,start=30000)

x_summ <- summary(chainsX)
sim_fit <- x_summ$quantiles[,c(1,3,5)]
sim_fit <- rbind(as.data.frame(matrix(sim_fit[,1],ncol=3)),
                 as.data.frame(matrix(sim_fit[,2],ncol=3)),
                 as.data.frame(matrix(sim_fit[,3],ncol=3)))
sim_fit <- sim_fit/theta_fit_dat['popsize']
names(sim_fit) <- genotype_names
sim_fit$gen <- seq(start_gen,end_gen)
sim_fit$quant <- rep(c(2.5,50,97.5),each=tend)

sim_fit_plot <- sim_fit %>% 
  melt(measure.vars=yvars,variable.name="genotype",value.name="frequency")

parnames <- c("sSS","h","r0")

{estim_dat_plot <- estim_dat
  estim_dat_plot$gen <- seq(start_gen,end_gen)
  estim_dat_plot$n <- rowSums(estim_dat_plot[,genotype_names])
  estim_dat_plot$n[estim_dat_plot$n>300] <- 300
  estim_dat_plot[,genotype_names] <- estim_dat_plot[,genotype_names]/rowSums(estim_dat_plot[,genotype_names])
  estim_dat_plot['R allele'] <- estim_dat_plot$RR + estim_dat_plot$SR/2
  estim_dat_plot <- estim_dat_plot %>%
    melt(measure.var=yvars,variable.name="genotype",value.name="frequency")}

# change levels for use as y-axis labels
add_label <- function(text_in){paste0(text_in," frequency")}
levels(estim_dat_plot$genotype) <- add_label(levels(estim_dat_plot$genotype))
levels(sim_fit_plot$genotype) <- add_label(levels(sim_fit_plot$genotype))

# make new data frame for credible interval ribbon
sim_fit_rib1 <- sim_fit_plot %>% subset(quant==2.5)
names(sim_fit_rib1) <- c("gen","quant","genotype","min")
sim_fit_rib2 <- sim_fit_plot %>% subset(quant==97.5)
names(sim_fit_rib2) <- c("gen","quant","genotype","max")
sim_fit_rib <- bind_cols(sim_fit_rib1,sim_fit_rib2)

my_colors <- c("true\n value"="black","PMMH mean\n estimate"="red1","PMMH 95% CI"="red3")
my_breaks <- seq(37,121,12)
my_labels <- seq(2003,2010)

sim_fit_plot %>%
  ggplot(aes(x=gen,y=frequency,group=interaction(genotype))) +
  geom_line(data=sim_fit_plot %>% subset(quant==50),
            aes(color="PMMH mean\n estimate"),
            size=1.2) +
  geom_ribbon(data=sim_fit_rib,
              aes(x=gen,y=NULL,color="PMMH mean\n estimate",fill="PMMH 95% CI",
                  ymin=min,ymax=max),
              size=0.4,alpha=0.2) +
  geom_point(data=estim_dat_plot,
             aes(size=n),
             shape=21,
             color="black",fill="#5A5A5A") +
  scale_size(name="sample\n size",
             labels=c(0,100,200,expression("">=300)))+
  scale_color_manual(name=NULL,values=my_colors)+
  scale_fill_manual(name=NULL,values=my_colors)+
  guides(
    size = guide_legend(order = 1),
    color = guide_legend(order = 2),
    fill = guide_legend(order = 3)
  ) +
  geom_blank(data=data.frame(genotype="RR frequency",size=0,gen=121,frequency=0)) +
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  facet_grid(genotype ~ .,switch="y") +
  scale_x_continuous(name = "Year",breaks = my_breaks,labels = my_labels) +
  ylab(NULL) +
  my_theme +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        strip.text = element_text(size=14))
