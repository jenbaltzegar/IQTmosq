##########################################################
# This code does the following: 
# - load in data and choose which months to use
# - calculate MLE estimate
# - plot deterministic simulations with different parameter values (Figs 6, S1)
##########################################################

## load packages and functions from estim_setup.R
source('model_setup.R')
## read in data for locus 1534 or 1016 from current working directory
# locus <- '1016' # *only for plotting data*
locus <- '1534'
dat <- read.csv(paste0('mc.',locus,'.MoYr_reduced_byMonth.csv'))
# make new columns for month and year
dat <- cbind(dat,colsplit(as.character(dat$MonthYear),"-",c("month","year")))
dat$year <- dat$year + 2000
dat$gen <- seq_len(nrow(dat))
start_gen <- if(locus=='1534'){34}else{97}
end_gen <- if(locus=='1534'){120}else{216}
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
theta1["r0"] <- dat$freqR[start_gen]
if(locus=='1016'){theta1["r0"] <- 0.0001}
# maximum likelihood estimate from function in model_setup.R
fit <- get_MLE(estim_dat,theta1,parnames)
fit
theta_fit_dat <- fit[names(theta1)]

# run deterministic simulation with MLE parameters
simdat <- WF_sim(theta=theta_fit_dat,n_gens=tend,stoch_on=0) %>%
  as.data.frame()
simdat$gen <- seq(start_gen,end_gen)
simdat$freqR <- simdat$RR + simdat$SR/2
yvars <- genotype_names
simdat_plot <- simdat %>% melt(measure.var=yvars,variable.name="genotype",value.name="frequency")

####### produce fits for different fixed parameters
# choose parameters to vary over (h set to vary)
parvals1 <- expand.grid("start_gen" = start_gen,
                        "end_gen" = end_gen,
                        "h"=as.vector(c(seq(0,1,.25))),#,theta_fit_dat["h"])),
                        "sSS" = as.vector(theta_fit_dat["sSS"]),
                        "r0" = as.vector(theta1["r0"])
                        )

# choose parameters to estimate
parnames_temp <- c("sSS")

fits_df <- foreach(x=iter(parvals1,by="row"), .combine=rbind,.packages=c("magrittr")) %do% {
  theta_temp <- c("sSS"=x[["sSS"]],"h"=x[["h"]],"r0"=x[["r0"]],"popsize"=500)
  estim_dat_temp <- dat[dat$gen >= x[["start_gen"]] & dat$gen <= x[["end_gen"]],genotype_names]
  par_fit <- get_MLE(dat=estim_dat_temp,theta=theta_temp,parnames=parnames_temp)
  data.frame(start_gen=x[["start_gen"]],end_gen=x[["end_gen"]],
             sSS = par_fit["sSS"],h=par_fit["h"],r0=par_fit["r0"])
}
fits_df$num <- seq_len(nrow(fits_df))

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

#set up data to plot
dat_plot <- subs_dat2
dat_plot[,genotype_names] <- dat_plot[,genotype_names]/rowSums(dat_plot[,genotype_names])
max_samp_size_plotted <- 300
dat_plot$n[dat_plot$n>max_samp_size_plotted] <- max_samp_size_plotted

#set up simulations to plot
all_simdat$freqR <- all_simdat$RR + all_simdat$SR/2

#melt by allele frequency 
yvars <- "freqR"
simdat_plot <- all_simdat %>% na.omit() %>%
  melt(measure.var=yvars,variable.name="genotype",value.name="frequency")
simdat_plot$h <- as.factor(round(simdat_plot$h,3))
dat_plot_allele <- dat_plot %>% na.omit() %>%
  melt(measure.var=yvars,variable.name="genotype",value.name="frequency")

#melt by genotype frequency 
yvars <- genotype_names
simdat_plot_genos <- all_simdat %>% na.omit() %>%
  melt(measure.var=yvars,variable.name="genotype",value.name="frequency")
simdat_plot_genos$h <- as.factor(round(simdat_plot_genos$h,3))
dat_plot_genos <- dat_plot %>% na.omit() %>%
  melt(measure.var=yvars,variable.name="genotype",value.name="frequency")

#set up h lines
vir_colors <- viridis(5)
mysizes <- rep(1.5,5)
myalpha <- rep(0.8,5)

# for 1016, produce figures without simulation results
if(locus=="1016"){
  simdat_plot <- simdat_plot %>% subset(gen<0)
  dat_plot_allele <- dat_plot_allele %>% subset(year>=minyear)
  simdat_plot_genos <- simdat_plot_genos %>% subset(gen<0)
  dat_plot_genos <- dat_plot_genos %>% subset(year>=minyear)
  my_legend_position <- c(0.75,0.22)
}else{
  my_legend_position <- c(0.85,0.26)
}

# Plot R allele frequency
simplot <- simdat_plot %>%
  ggplot(aes(x=gen,y=frequency)) +
  geom_line(aes(color=h,alpha=h,
                group=interaction(genotype,simnum)),size=2) +
  scale_color_manual(name="h",values=vir_colors)+
  scale_alpha_manual(name="h",values=myalpha)+ 
  
  geom_point(data=dat_plot_allele,aes(size=n),shape=21,
             color="black",fill="#5A5A5A") +
  scale_size(name="sample\n size",
             labels=c(0,100,200,expression("">=300)))+
  
  scale_y_continuous(expand = c(0.03,0))+
  scale_x_continuous(name = "Year",breaks = seq(1,12*19,12),labels = seq(2000,2018)) +
  # scale_x_continuous(name = "Year",breaks = dat$gen[year_breaks],labels = dat$year[year_breaks]) +
  ylab("R Allele Frequency") +
  
  my_theme +
  theme(
        legend.position=my_legend_position,
        strip.background = element_blank(),
        legend.background = element_rect(color=NA),
        legend.key = element_rect(color="white"),
        legend.title.align=0.5,
        legend.text.align = 0
        ) +
  guides(color = guide_legend(override.aes = list(size = 3)))
simplot
pdf(file = paste0("time_series_Rallele_",locus,".pdf"),width = 9.89,height=7.72)
simplot
dev.off()

# Plot genotype frequencies
levels(simdat_plot_genos$genotype) <- paste(levels(simdat_plot_genos$genotype),"Frequency")
levels(dat_plot_genos$genotype) <- paste(levels(dat_plot_genos$genotype),"Frequency")
simplot_genos <- simdat_plot_genos %>% 
  ggplot(aes(x=gen,y=frequency)) +
  geom_line(aes(color=h,alpha=h,
                group=interaction(genotype,simnum)),size=2) +
  scale_color_manual(name="h",values=vir_colors)+
  scale_alpha_manual(name="h",values=myalpha)+ 
  
  geom_point(data=dat_plot_genos,aes(size=n),shape=21,
             color="black",fill="#5A5A5A") +
  scale_size(name="sample\n size",
             labels=c(0,100,200,expression("">=300)))+
  
  facet_grid(genotype ~ .,switch="y")+
  
  scale_y_continuous(expand = c(0.03,0))+
  scale_x_continuous(name = "Year",breaks = seq(1,12*19,12),labels = seq(2000,2018)) +
  ylab(NULL) +
  
  my_theme +
  theme(
    strip.background = element_blank(),
    strip.text =element_text(size=14),
    strip.placement = "outside",
    legend.background = element_rect(color=NA),
    legend.key = element_rect(color="white"),
    legend.title.align=0.5,
    legend.text.align = 0
  ) +
  guides(color = guide_legend(override.aes = list(size = 3)))
simplot_genos
pdf(file = paste0("time_series_genotype_",locus,".pdf"),width = 9.89,height=9.1)
simplot_genos
dev.off()
