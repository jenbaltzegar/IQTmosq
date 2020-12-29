##########################################################
# This code does the following: 
# - load in data and choose which months to use
# - calculate MLE estimate and 95% CI through parametric bootstrapping
# - plot deterministic simulations with different parameter values
# - plot 2-dimensional log likelihood surface
# - plot likelihood profiles
##########################################################

## load packages and functions from estim_setup.R
source('model_setup.R')
## read in data for locus 1534 or 1016 from current working directory
locus <- '1016'
locus <- '1534'
dat <- read.csv(paste0('data/mc.',locus,'.MoYr_reduced_byMonth.csv'))
# make new columns for month and year
dat <- cbind(dat,colsplit(as.character(dat$MonthYear),"-",c("month","year")))
dat$year <- dat$year
dat$gen <- seq_len(nrow(dat))
start_gen <- 34
end_gen <- 120
tend <- end_gen - start_gen + 1

#subset data by specified generations 
estim_dat <- dat[dat$gen >= start_gen & dat$gen <= end_gen,genotype_names] #genotype columns only
subs_dat <- dat[dat$gen >= start_gen & dat$gen <= end_gen,] #all columns
subs_dat2 <- dat[dat$gen <= end_gen,] #all columns

# choose parameters to estimate. fitness parameters, initial R allele frequency
# set parameters to appropriate constant if not in parnames (not going to estimate)
# popsize is only used in stochastic simulations
theta1 <- c("sSS"=0.3,"h"=0.05,"r0"=0.01,"popsize"=500)
parnames <- c("sSS","h")
# parnames <- c("sSS","h","r0")
theta1["r0"] <- dat$freqR[start_gen]
# theta1["r0"] <- 0.01

# maximum likelihood estimate from function in model_setup.R
fit <- get_MLE(estim_dat,theta1,parnames)
fit
theta_fit_dat <- fit[names(theta1)]
# theta_fit_dat['popsize'] <- 200
# run deterministic simulation with MLE parameters
simdat <- WF_sim(theta=theta_fit_dat,n_gens=tend,stoch_on=0) %>%
  as.data.frame()
simdat$gen <- seq(start_gen,end_gen)
simdat$freqR <- simdat$RR + simdat$SR/2
yvars <- genotype_names
simdat_plot <- simdat %>% melt(measure.var=yvars,variable.name="genotype",value.name="frequency")

## parametric bootstrap to produce 95% CI
nboots=500
bootout <- as.data.frame(matrix(-1,ncol=length(parnames),nrow=nboots))
names(bootout) <- parnames
theta_boot <- theta_fit_dat
all_s_dat <- list()
for (i in 1:nrow(bootout)){
  s_dat <- WF_sampler(WF_sim(theta=theta_boot,n_gens=tend,stoch_on=1),
                      rowSums(estim_dat))
  #condition on simulations that R allele doesn't fall out of the population
  #(not necessary if using deterministic simulation)
  # if(s_dat[max(which(rowSums(estim_dat)>0)),"RR"]>0){
    paramfit_temp <- get_MLE(s_dat,theta_boot,parnames)
    bootout[i,parnames] <- paramfit_temp[parnames]
  # }else{
  #   bootout[i,parnames] <- rep(NA,length(parnames))
  # }
  s_dat$simnum <- i
  s_dat$gen <- seq(start_gen,end_gen)
  all_s_dat[[i]] <- s_dat
}
sum(is.na(bootout$h))
newbootout <- bootout[!is.na(bootout$h),]

CIs <- do.call(rbind,
               lapply(parnames,function(x) 
                 data.frame("parameter"=x,"estimate"=theta_boot[x],
                            "CI_L"=quantile(newbootout[,x],0.025),
                            "CI_H"=quantile(newbootout[,x],0.975))
               )
)

CIs #<- CIs
boot_pairplot <- newbootout %>% ggplot(aes(x=sSS,y=h)) +
  geom_point() +
  geom_point(data=data.frame("h"=theta_fit_dat["h"],"sSS"=theta_fit_dat["sSS"]),color="red",size=3) +
  geom_point(data=data.frame("h"=mean(newbootout$h),"sSS"=mean(newbootout$sSS)),color="green",size=2) +
  # coord_fixed(ratio=1)+
  my_theme
boot_pairplot
# pairwise <- ggpairs(newbootout) + my_theme

boot_boxplot <- newbootout %>% 
  melt(measure.var=parnames,value.name="estimate",variable.name="parameter") %>%
  ggplot(aes(x=parameter,y=estimate)) +
  geom_boxplot() +
  geom_point(data=CIs,color="red")+
  my_theme
boot_boxplot

grid.arrange(boot_pairplot,boot_boxplot,nrow=1)

####### produce fits for different fixed parameters
# choose parameters to vary over (h set to vary)
parvals1 <- expand.grid("start_gen" = start_gen,
                        "end_gen" = end_gen,
                        "h"=as.vector(c(seq(0,1,.25))),#,theta_fit_dat["h"])),
                        "sSS" = as.vector(theta_fit_dat["sSS"]),
                        # "sSS" = 0.8,
                        "r0" = 0.02
                        # "r0" = as.vector(theta_fit_dat["r0"])
                        )

# choose parameters to estimate
parnames_temp <- c("sSS")
parnames_temp <- c("r0")
# parnames_temp <- c("sSS","r0")

fits_df <- foreach(x=iter(parvals1,by="row"), .combine=rbind,.packages=c("magrittr")) %do% {
  theta_temp <- c("sSS"=x[["sSS"]],"h"=x[["h"]],"r0"=x[["r0"]],"popsize"=500)
  estim_dat_temp <- dat[dat$gen >= x[["start_gen"]] & dat$gen <= x[["end_gen"]],genotype_names]
  par_fit <- get_MLE(dat=estim_dat_temp,theta=theta_temp,parnames=parnames_temp)
  data.frame(start_gen=x[["start_gen"]],end_gen=x[["end_gen"]],
             sSS = par_fit["sSS"],h=par_fit["h"],r0=par_fit["r0"])
}
fits_df$num <- seq_len(nrow(fits_df))
fits_df$sSS <- 0.8
fits_df$r0 <- 0.02
# run simulations for each row of parameter vales
all_simdat <- foreach(x=iter(fits_df,by="row"),.combine=rbind,.packages=c("magrittr")) %do% {
  tend <- x[["end_gen"]] - x[["start_gen"]] + 1
  theta_temp <- c("sSS"=x[["sSS"]],"h"=x[["h"]],"r0"=x[["r0"]],"popsize"=500)
  simdat <- WF_sim(theta=theta_temp,n_gens=tend,stoch_on=0) %>%
    as.data.frame()
  simdat$gen <- seq(x[["start_gen"]],x[["end_gen"]])
  simdat$h <- x[["h"]]
  simdat$sSS <- x[["sSS"]]
  simdat$r0 <- x[["r0"]]
  simdat$simnum <- x[["num"]]
  simdat
}
#choose whether to plot genotypes or allele frequency
yvars <- genotype_names
yvars <- "freqR"

#set up simulations to plot
all_simdat$freqR <- all_simdat$RR + all_simdat$SR/2
simdat_plot <- all_simdat %>% na.omit() %>%
  melt(measure.var=yvars,variable.name="genotype",value.name="frequency")
simdat_plot$h <- as.factor(round(simdat_plot$h,3))

#set up data to plot
dat_plot <- subs_dat2
dat_plot[,genotype_names] <- dat_plot[,genotype_names]/rowSums(dat_plot[,genotype_names])
max_samp_size_plotted <- 300
dat_plot$n[dat_plot$n>max_samp_size_plotted] <- max_samp_size_plotted
dat_plot <- dat_plot %>% na.omit() %>%
  melt(measure.var=yvars,variable.name="genotype",value.name="frequency")

#set up x-axis breaks
year_breaks <- seq(min(which(dat$year==min(dat[dat$month=="Jan","year"]))),
                   nrow(dat),12)

#set up h lines
#need to modify for in the case that theta_fit_dat["h"] == h value in expand.grid
vir_colors <- viridis(5)
which_fit_val <- which(levels(simdat_plot$h)==round(theta_fit_dat["h"],3)) #make proper line stand out
my_colors <- vir_colors
mysizes <- rep(1.5,5)
myalpha <- rep(0.8,5)
# my_colors <- append(vir_colors,"red",which_fit_val-1)
# mysizes <- append(rep(1.5,5),2,which_fit_val-1)
# myalpha <- append(rep(0.6,5),1,which_fit_val-1)
# mylinetypes <- append(rep("dashed",5),"solid",which_fit_val-1)

#To do: move facet labels to left-hand side: R allele frequency/ RR frequency, SR....
#       ***legend may need adjustment***. needs to move if plotting genotype frequencies
simplot <- simdat_plot %>% #subset(h>0.7) %>%
  ggplot(aes(x=gen,y=frequency)) +
  geom_line(aes(color=h,alpha=h,#linetype=h,
                group=interaction(genotype,simnum)),size=2) +
  scale_color_manual(name="h",values=my_colors)+
  scale_alpha_manual(name="h",values=myalpha)+ 
  # scale_linetype_manual(name="h",values=mylinetypes)+
  # scale_size_manual(name="h",values=mysizes)+ #need to add new scale to add line sizes
  
  geom_point(data=dat_plot,aes(size=n),shape=21,
             color="black",fill="#5A5A5A") +
  scale_size(name="sample\n size", #need to manually change currently
             labels=c(0,100,200,expression("">=300)))+
  
  facet_grid(genotype ~ .)+#,scales = "free_y") +
  
  scale_y_continuous(expand = c(0.03,0))+
  scale_x_continuous(name = "Year",breaks = dat$gen[year_breaks],labels = dat$year[year_breaks]) +
  ylab("R Allele Frequency") +
  # ggtitle("Simulations with different h values") +
  # theme(plot.title = element_text(hjust = 0.5)) +
  my_theme +
  theme(
        legend.position=c(0.85,0.35),
        strip.background = element_blank(),
        strip.text.y = element_blank(),
        legend.background = element_rect(color="black"),
        legend.key = element_rect(color="white"),
        legend.title.align=0.5,
        legend.text.align = 0
        ) +
  guides(color = guide_legend(override.aes = list(size = 3)))
simplot

####### 2-dimensional log likelihood
parvals2 <- expand.grid(h=seq(0,1,0.01),sSS=seq(0,0.99,0.01),
                       # "r0"=dat$freqR[start_gen])
                       "r0"=theta_fit_dat["r0"])
                       # "r0"=c(0.001,0.01,.1,0.2))

pars_set <- names(parvals2)
parvals2$LL <- -1
theta_temp <- theta_fit_dat

#get log likelihood for each row of parameter values
for (i in 1:nrow(parvals2)){
  theta_temp[pars_set] <- unlist(parvals2[i,pars_set])
  simdat_temp <- WF_sim(theta=theta_temp,n_gens=tend,stoch_on=0)
  parvals2[i,"LL"] <- -costfunc(simdat_temp,estim_dat)
}

# produce figure
parvals2_plot <- parvals2
# alternative sSS and h parameterization
# parvals2_plot$alt_sSS <- 1/(1-parvals2_plot$sSS)-1
# parvals2_plot$alt_h <- ((1-parvals2_plot$h*parvals2_plot$sSS)/(1-parvals2_plot$sSS)-1)/
#                       (parvals2_plot$alt_sSS)

# set up breaks (takes some manual configuration currently)
break_vals <- c(min(parvals2_plot$LL),-10000,-5000,-3000,
                max(parvals2_plot$LL)-100,max(parvals2_plot$LL))
# break_vals <- c(-10000,-3000,-2000,-1500,-1000,max(parvals2_plot$LL)-100,max(parvals2_plot$LL))
parvals2_plot$LL_cut <- cut(parvals2$LL,break_vals)
parvals2_plot$LL_cut[is.na(parvals2_plot$LL_cut)] <- levels(parvals2_plot$LL_cut)[1]
my_heat_colors <- brewer.pal(length(break_vals),"YlOrRd")

# TO DO: improve legend. make lines smoother? Add contours? improve point visibility
gridplot <- parvals2_plot %>%
  ggplot(aes(x=sSS,y=h)) +
  geom_tile(aes(fill=LL_cut))+
  scale_fill_manual(name="log likelihood",values=my_heat_colors) +
  
  scale_y_continuous(limits=c(0,1), #name = "Cost SS",
                     expand=c(0,0)) +
  scale_x_continuous(limits=c(0,1), #name = "Dominance",
                     expand=c(0,0)) +
  coord_fixed(ratio = 1) +
  
  geom_point(data=fits_df,size=2.5,shape="square", #outline
             aes(x=sSS,y=h),color="black") +
  geom_point(data=fits_df,size=2,shape="square", #fits from above
             aes(x=sSS,y=h,color=as.factor(h))) +

  scale_color_manual(values = my_colors)+
  guides(color=F)+
  my_theme
gridplot

grid.arrange(simplot,gridplot) #need to properly configure legend if doing this

## asymptotic 95% CI based on MLE assumptions (could also obtain below)
# (need to verify calculation) and use finer grid of estimates
ML_pars <- parvals2[which(parvals2$LL ==max(parvals2$LL)),] #should match theta_fit_dat
ML_CI <- parvals2[which(parvals2$LL > max(parvals2$LL)-0.5*qchisq(p=0.95,df=2)),]
ML_CI[c(1,nrow(ML_CI)),]

############# likelihood profile
parnames <- c("sSS","h","r0")

#choose which parameter to profile over and set values for that parameter
profile_var <- "h"    #seq(0,1,.01)
prof_vals <- as.data.frame(seq(0,1,.01))

profile_var <- "r0"    #seq(0.01,0.99,.01)
prof_vals <- as.data.frame(seq(0.001,0.3,.001))

profile_var <- "sSS" #seq(0,.99,.01)
prof_vals <- as.data.frame(seq(0,0.5,.01))

#set up data frame
names(prof_vals) <- as.character(profile_var)
new_parnames <- parnames[parnames!=profile_var] #maximize over all other parameters
prof_vals[as.character(new_parnames)] <- -1
prof_vals$LL <- -1

theta2 <- theta_fit_dat

for (i in 1:nrow(prof_vals)){
  theta2[profile_var] <- prof_vals[i,profile_var]
  MLE_temp <- get_MLE(estim_dat,theta2,new_parnames)
  prof_vals[i,new_parnames] <- MLE_temp[new_parnames]
  prof_vals[i,"LL"] <- MLE_temp["LL"]
}

head(prof_vals)
prof_vals %>% melt(names(prof_vals[1])) %>%
  ggplot() +
  geom_line(aes_(x=as.name(profile_var),y=quote(value))) +
  facet_grid(variable~.,scale="free_y") +
  scale_x_continuous(limits=c(0,max(prof_vals[profile_var])),
                     name = as.name(profile_var),expand=c(0,0)) +
  my_theme

###############################
