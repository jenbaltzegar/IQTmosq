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

# create custom theme - from Michael Vella -----
library(ggthemes)
my_theme <- function(){
  theme_foundation(base_size=14) + 
    theme(
      plot.title = element_text(face = "bold",size = rel(1.2), hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      panel.background = element_rect(colour = NA),
      plot.background = element_rect(colour = NA),
      panel.border = element_rect(colour = NA),
      panel.spacing=unit(1, "lines"),
      axis.title = element_text(face = "bold"),
      axis.title.y = element_text(vjust = 2),
      axis.title.x = element_text(vjust = -0.2),
      axis.line = element_line(colour="black"),
      # panel.grid.major = element_line(colour="#f0f0f0"),
      # panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.margin=unit(c(10,5,5,5),"mm"),
      strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
      strip.text = element_text(face="bold")
    )
}

# Create custom color palettes -----
hapColors <- c("#117733", "#332288", "#CC6677", "#882255")
names(hapColors) <- rev(levels(freqAll_long$Haplotype))
hapScale <- scale_colour_manual(name = "Haplotype",values = hapColors)

# Plot kdr Haplotypes -----
option <- "Option_7"
kdrHaps <- ggplot(data = freqAll_long[!is.na(freqAll_long$Frequency),], aes(x=year, y=Frequency)) +
  my_theme() +
  labs(x = "Year", y = "Frequency", title = "kdr Haplotype Frequency",
       subtitle = option) +
  scale_y_continuous(limits=c(-0.05, 1), breaks=c(0.0,0.2, 0.4, 0.6, 0.8, 1.0)) +
  scale_x_continuous(breaks = pretty(freqAll_long$year, n = 18)) +

  ### Add n to each year
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
  
  ### Add Legends
  hapScale 


# Print graph to screen
kdrHaps

# Write plot to png
ggsave(filename = paste0("~/Desktop/kdrHaps_options/", option, ".png"), width = 11, height = 8, dpi = 600, units = "in", device='png')



# # scratch to pick best graph
# # option 1 - change alpha
# kdrHaps
# 
# # option 2 - with custom sprayColors
# kdrHaps +
#   scale_fill_manual(name = 'Pyrethroid Application',
#                     # values = sprayColors,  
#                     guide = guide_legend(override.aes = list(alpha = 1))) +
#   my_theme()
# 
# # option 3 - with gray.colors(6)
# kdrHaps +
#   scale_fill_manual('Spraying Regime',
#                     values = gray.colors(6),  
#                     guide = guide_legend(override.aes = list(alpha = 1))) 
# 
# kdrHaps +
#   scale_fill_grey('Spraying Regime'
#                   , guide = guide_legend(override.aes = list(alpha = 1))) 
# 
# kdrHaps +
#   my_theme()

