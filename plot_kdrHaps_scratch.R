# Plot the frequency of kdr haplotypes over time with ggplot()

# Load required libraries
library(ggplot2)
library(RColorBrewer)

# Load data for plotting
freqAll_long <- read.csv("freqAll_long.csv")

# reassign Haplotypes to have specific text for legend
levels(freqAll_long$Haplotype)[match("SS", levels(freqAll_long$Haplotype))] <- "Val1016/Phe1534"
levels(freqAll_long$Haplotype)[match("SR", levels(freqAll_long$Haplotype))] <- "Val1016/Cys1534"
levels(freqAll_long$Haplotype)[match("RS", levels(freqAll_long$Haplotype))] <- "Ile1016/Phe1534"
levels(freqAll_long$Haplotype)[match("RR", levels(freqAll_long$Haplotype))] <- "Ile1016/Cys1534"

# # Load mc.1534.yr and mc.1016.yr for selection coefficient annotation
# mc.1534.yr <- read.csv("mc.1534.yr_withSelCoeff.csv")
# mc.1016.yr <- read.csv("mc.1016.yr_withSelCoeff.csv")
# mc.haps.yr <- read.csv("mc.haps.yr_reduced.csv")

# Create custom color palettes
hapColors <- c("#117733", "#332288", "#CC6677", "#882255")
names(hapColors) <- rev(levels(freqAll_long$Haplotype))
colScale <- scale_colour_manual(name = "Haplotype",values = hapColors)

sprayColors <- c("#A0A0A0", "#D1D1D1","#BBBBBB", "#E6E6E6", "#7F7F7F", "#FFFFFF")


##########################
# Key code for this script
##########################
# Split data at year 2013
# head(freqAll_long)

# To plot-----
kdrHaps <- ggplot(data = freqAll_long[!is.na(freqAll_long$Frequency),], aes(x=year, y=Frequency)) +
  theme_bw() + #removes grey background  
  theme(plot.background = element_blank()
        , panel.grid.major = element_blank()
        , panel.grid.minor = element_blank()
        , panel.border = element_blank()
        , axis.line = element_line(color = 'black')) +
  # xlim(c(2000, 2017)) +
  #ylim(c(-0.1, 1.1)) +
  scale_y_continuous(limits=c(-0.05, 1.1), breaks=c(0.0,0.2, 0.4, 0.6, 0.8, 1.0)) +
  scale_x_continuous(breaks = pretty(freqAll_long$year, n = 18)) +
  # scale_colour_discrete(name ="Haplotype", labels=c("wt1016/wt1534", "wt1016/Cys1534", "Ile1016/wt1534", "Ile1016/Cys1534")) +
  #geom_hline(yintercept = 0) +

  # ### OLD!!! Add background to represent type of insecticides used -----
  # # # Pyrethroids 2002 - 2014
  # # geom_rect(data = NULL, aes(xmin = 2001.5, xmax = 2014.5, ymin = -Inf, ymax = Inf), fill="#DDCC77", alpha = 0.005) +
  # # No pyrethroids
  # geom_rect(data = NULL, aes(xmin = 1999.0, xmax = 2001.5, ymin = -Inf, ymax = Inf)
  #           # , fill="No Pyrethroids"
  #           , alpha = 0.1) +
  # # Delatmethrin 2002 - 2007
  # geom_rect(data = NULL, aes(xmin = 2001.5, xmax = 2007.5, ymin = -Inf, ymax = Inf)
  #           # , fill="Deltamethrin"
  #           , alpha = 0.1) +
  # # Cypermethrin 2008
  # geom_rect(data = NULL, aes(xmin = 2007.5, xmax = 2008.5, ymin = -Inf, ymax = Inf)
  #           # , fill="Cypermethrin"
  #           , alpha = 0.1) +
  # # a-Cypermethrin 2009 - 2011
  # geom_rect(data = NULL, aes(xmin = 2008.5, xmax = 2011.5, ymin = -Inf, ymax = Inf)
  #           # , fill="a-Cypermethrin"
  #           , alpha = 0.1) +
  # # Cypermethrin & a-Cypermethrin 2012
  # geom_rect(data = NULL, aes(xmin = 2011.5, xmax = 2012.5, ymin = -Inf, ymax = Inf)
  #           # , fill="Cypermethrin & a-Cypermethrin"
  #           , alpha = 0.1) +
  # # Lambda & Alpha & Cyper & Alpha+pyriprox in 2013 
  # geom_rect(data = NULL, aes(xmin = 2012.5, xmax = 2013.5, ymin = -Inf, ymax = Inf)
  #           # , fill="Multiple*"
  #           , alpha = 0.1) +
  # # Cypermethrin 2014
  # geom_rect(data = NULL, aes(xmin = 2013.5, xmax = 2014.5, ymin = -Inf, ymax = Inf)
  #           # , fill="Cypermethrin"
  #           , alpha = 0.1) +
  # # No pyrethroids
  # geom_rect(data = NULL, aes(xmin = 2014.5, xmax = 2018.0, ymin = -Inf, ymax = Inf)
  #           # , fill="No Pyrethroids"
  #           , alpha = 0.1) +
  # 
  ### NEW!!! Add background to represent type of insecticides used -----
  # # Pyrethroids 2002 - 2014
  # geom_rect(data = NULL, aes(xmin = 2001.5, xmax = 2014.5, ymin = -Inf, ymax = Inf), fill="#DDCC77", alpha = 0.005) +
  # No pyrethroids
  geom_rect(data = NULL, aes(xmin = 1999.0, xmax = 2001.5, ymin = -Inf, ymax = Inf, fill="No Pyrethroids"), colour=NA, alpha = 1) +
  # Delatmethrin 2002 - 2007
  geom_rect(data = NULL, aes(xmin = 2001.5, xmax = 2007.5, ymin = -Inf, ymax = Inf, fill="Deltamethrin"), colour=NA, alpha = 1) +
  # Cypermethrin 2008
  geom_rect(data = NULL, aes(xmin = 2007.5, xmax = 2008.5, ymin = -Inf, ymax = Inf, fill="Cypermethrin"), colour=NA, alpha = 1) +
  # a-Cypermethrin 2009 - 2011
  geom_rect(data = NULL, aes(xmin = 2008.5, xmax = 2011.5, ymin = -Inf, ymax = Inf, fill="a-Cypermethrin"), colour=NA, alpha = 1) +
  # Cypermethrin & a-Cypermethrin 2012
  geom_rect(data = NULL, aes(xmin = 2011.5, xmax = 2012.5, ymin = -Inf, ymax = Inf, fill="Cypermethrin & a-Cypermethrin"), colour=NA, alpha = 1) +
  # Lambda & Alpha & Cyper & Alpha+pyriprox in 2013 ## "Lambda, a-Cypermethrin, Cypermethrin, & a-Cypermethrin+pyriporx"
  geom_rect(data = NULL, aes(xmin = 2012.5, xmax = 2013.5, ymin = -Inf, ymax = Inf, fill="Multiple*"), colour=NA, alpha = 1) +
  # Cypermethrin 2014
  geom_rect(data = NULL, aes(xmin = 2013.5, xmax = 2014.5, ymin = -Inf, ymax = Inf, fill="Cypermethrin"), colour=NA, alpha = 1) +
  # No pyrethroids
  geom_rect(data = NULL, aes(xmin = 2014.5, xmax = 2018.0, ymin = -Inf, ymax = Inf, fill="No Pyrethroids"), colour=NA, alpha = 1) +

  # ### Add text to represent type of insecticides used -----
  # annotate("text", x= 2008.0, y=1.1, label = "Pyrethroids", size = 5, fontface = 2) +
  # annotate("text", x=2000.25, y=1.1, label = "No \n Pyrethroids", size = 5, fontface = 2) +
  # annotate("text", x=2004.5, y=1.05, label= "Delta", size = 4) +
  # annotate("text", x=2008.0, y=1.05, label= "Cyper", size = 4) +
  # # annotate("text", x=2010, y=1.05, label= expression(alpha), size = 4) +
  # annotate("text", x=2010, y=1.05, label= "Alpha", size = 4) +
  # # annotate("text", x=2012.5, y=1.05, label= expression(alpha), size = 4) +
  # annotate("text", x=2012.5, y=1.05, label= "Alpha & \n Cyper", size = 4) +
  # # annotate("text", x=2012.5, y=1, label= "& \n Cyper", size = 4) +
  # annotate("text", x=2014.0, y=1.05, label = "Cyper", size = 4) +
  # annotate("text", x=2016.0, y=1.1, label = "No \n Pyrethroids", size = 5, fontface = 2) +
  # 
  
  ### Add n to each year
  annotate("text", x=1999.3, y=-0.03, label = "n =", fontface = 2, color = "#332288") +
  annotate("text", x=2000, y=-0.03, label = mc.haps.yr$n[1], fontface = 2, color = "#332288") +
  annotate("text", x=2001, y=-0.03, label = mc.haps.yr$n[2], fontface = 2, color = "#332288") +
  annotate("text", x=2002, y=-0.03, label = mc.haps.yr$n[3], fontface = 2, color = "#332288") +
  annotate("text", x=2003, y=-0.03, label = mc.haps.yr$n[4], fontface = 2, color = "#332288") +
  annotate("text", x=2004, y=-0.03, label = mc.haps.yr$n[5], fontface = 2, color = "#332288") +
  annotate("text", x=2005, y=-0.03, label = mc.haps.yr$n[6], fontface = 2, color = "#332288") +
  annotate("text", x=2006, y=-0.03, label = mc.haps.yr$n[7], fontface = 2, color = "#332288") +
  annotate("text", x=2007, y=-0.03, label = mc.haps.yr$n[8], fontface = 2, color = "#332288") +
  annotate("text", x=2008, y=-0.03, label = mc.haps.yr$n[9], fontface = 2, color = "#332288") +
  annotate("text", x=2009, y=-0.03, label = mc.haps.yr$n[10], fontface = 2, color = "#332288") +
  annotate("text", x=2010, y=-0.03, label = mc.haps.yr$n[11], fontface = 2, color = "#332288") +
  annotate("text", x=2011, y=-0.03, label = mc.haps.yr$n[12], fontface = 2, color = "#332288") +
  annotate("text", x=2012, y=-0.03, label = mc.haps.yr$n[13], fontface = 2, color = "#332288") +
  annotate("text", x=2013, y=-0.03, label = mc.haps.yr$n[14], fontface = 2, color = "#332288") +
  annotate("text", x=2014, y=-0.03, label = mc.haps.yr$n[15], fontface = 2, color = "#332288") +
  annotate("text", x=2015, y=-0.03, label = mc.haps.yr$n[16], fontface = 2, color = "#332288") +
  annotate("text", x=2016, y=-0.03, label = mc.haps.yr$n[17], fontface = 2, color = "#332288") +
  annotate("text", x=2017, y=-0.03, label = mc.haps.yr$n[18], fontface = 2, color = "#332288") +
  
  ### Add data
  geom_line(aes(colour=Haplotype), size = 2) +
  geom_point(aes(colour=Haplotype), size = 3) +
  geom_errorbar(data = freqAll_long, aes(ymin=Frequency-CI_95, ymax=Frequency+CI_95, colour=Haplotype), width=.2, size = 0.7) +
  
  labs(x = "Year", y = "Frequency", title = "kdr Haplotype Frequency"
       #, subtitle = "Error Bars = 95% CI"
       ) +
  colScale +
  
  # Adjusting theme
  theme(plot.subtitle=element_text(size=18, hjust=0.5, face="italic")
        , axis.title=element_text(size=14,face="bold")
        , plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
        , axis.text=element_text(size=12)
  ) +
  theme(legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size=12)
        , legend.title.align=0.5) 


# Print graph to screen
kdrHaps

# # Write plot to png
# ggsave(filename = paste0("figures/kdrHaps/kdrHaps_bars/kdrHaps_", Sys.Date(), ".png"), width = 11, height = 8, dpi = 600, units = "in", device='png')


# # Source function to move x-axis with labels
# source("R_Scripts/function_moveXaxis.R")
# # Plot x-axis moved up
# kdrHaps <- shift_axis(kdrHaps, 0)
# kdrHaps


# scratch to pick best graph
# option 1 - change alpha
kdrHaps

# option 2 - with custom sprayColors
kdrHaps +
  scale_fill_manual(name = 'Pyrethroid Application',
                    values = sprayColors,  
                    guide = guide_legend(override.aes = list(alpha = 1))) 

# option 3 - with gray.colors(6)
kdrHaps +
  scale_fill_manual('Spraying Regime',
                    values = gray.colors(6),  
                    guide = guide_legend(override.aes = list(alpha = 1))) 

kdrHaps +
  scale_fill_grey('Spraying Regime'
                  , guide = guide_legend(override.aes = list(alpha = 1))) 

