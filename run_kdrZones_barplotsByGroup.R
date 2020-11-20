# Prepare working environment ---------------------------------------------
# set working directory
setwd("~/Dropbox/GouldLab/Project_Mosquito/Database")

# clear working environment
rm(list = ls())

# Load libraries
library(gdata)
options(gsubfn.engine = "R") # prevents R from stalling while loading sqldf library
library(sqldf)
library(dplyr) # required for allele frequency calculations
library(ggplot2) # required to produce plots


# Read in dataframes
mc.1016.t <- read.csv("mc.1016.t.csv")
mc.1016.b <- read.csv("mc.1016.b.csv")

df.t <- mc.1016.t[, 1:5]
df.b <- mc.1016.b[, 1:5]

before.t <- df.t[1, ] + df.t[2, ] + df.t[3, ]
after.t <- df.t[4, ] + df.t[5, ] + df.t[6, ] + df.t[7, ] + df.t[8, ] + df.t[9, ] + df.t[10, ]
before.t; after.t

before.b <- df.b[1, ] + df.b[2, ] + df.b[3, ]
after.b <- df.b[4, ] + df.b[5, ] + df.b[6, ] + df.b[7, ] + df.b[8, ] + df.b[9, ] + df.b[10, ]
before.b; after.b

new.t <- rbind(before.t, after.t)
new.b <- rbind(before.b, after.b)
new.t; new.b

new.t$freqR <- ((2 * new.t$RR) + new.t$SR) / (2 * new.t$n)
new.b$freqR <- ((2 * new.b$RR) + new.b$SR) / (2 * new.b$n)
new.t; new.b

# Calculate 95% Confidence Interval
new.t$CI_95 = 1.96 * sqrt((new.t$freqR*(1-new.t$freqR))/(2*new.t$n))
new.b$CI_95 = 1.96 * sqrt((new.b$freqR*(1-new.b$freqR))/(2*new.b$n))
new.t; new.b

new.t$month <- c(3, 10)
new.b$month <- c(3, 10) 
new.t; new.b

new.t$sampled <- c("before", "post")
new.b$sampled <- c("before", "post")
new.t; new.b

new.t$zone <- "treatment"
new.b$zone <- "buffer"
new.t; new.b

df2 <- rbind(new.t, new.b)
df2


## Bar Graph #########################################
# Basic barplot
p <- ggplot(data=df2, aes(x=sampled, y=freqR, fill=zone)) +
  theme_bw() + #removes grey background 
  theme(plot.background = element_blank()
        , panel.grid.major = element_blank()
        , panel.grid.minor = element_blank()
        , panel.border = element_blank()
        , axis.line = element_line(color = 'dark grey')
        , plot.subtitle=element_text(size=18, hjust=0.5, face="italic")
        , axis.title=element_text(size=14,face="bold")
        , plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
        , axis.text=element_text(size=12)
        , legend.title = element_text(size = 14, face = "bold")
        , legend.text = element_text(size=12)) +
  geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=freqR-CI_95, ymax=freqR+CI_95), width=.2,
                position=position_dodge(.9)) +
  scale_fill_manual("Treatment", values = c("buffer" = "blue", "treatment" = "Red")) +
  scale_x_discrete(breaks=c("before", "post"), labels=c("Pre-Spray", "Post-Spray")) +
  xlab("Groups") +
  ylab("Haplotype Frequency")
p

# Write plot to png
ggsave(filename = paste0("figures/kdrZones/kdrZones_2014/kdrZones_barplot_", Sys.Date(), ".png"), width = 11, height = 8, dpi = 600, units = "in", device='png')



