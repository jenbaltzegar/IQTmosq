### Plot Frequency of Resistance Allele at 1016 locus by zone for year 2013
# Plot based off of dataframes mc.1016.t13 and mc.1016.b13
# Created: 09 Mar 2018

# Load Libraries
library(ggplot2)

# Read in dataframes
setwd("/Users/jenbaltz/Dropbox/GouldLab/Project_Mosquito/Database")
mc.1016.t13 <- read.csv("mc.1016.t13_reduced.csv")
mc.1016.b13 <- read.csv("mc.1016.b13_reduced_expandedBuffer.csv")

# 3/28/18 - Subset mc.1016.t13 to remove March & September data because sample size is too low
mc.1016.t13 <- mc.1016.t13[-c(3,9),]
mc.1016.b13 <- mc.1016.b13[-c(3,9),]

# create custom theme - from Michael Vella -----
library(ggthemes)
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

# Plot
kdrZones13 <- ggplot(mc.1016.t13, aes(x=month, y=freqR)) +
  my_theme() +
  # theme_bw() + #removes grey background 
  # theme(plot.background = element_blank()
  #       , panel.grid.major = element_blank()
  #       , panel.grid.minor = element_blank()
  #       , panel.border = element_blank()
  #       , axis.line = element_line(color = 'dark grey')
  #       , plot.subtitle=element_text(size=18, hjust=0.5, face="italic")
  #       , axis.title=element_text(size=14,face="bold")
  #       , plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
  #       , axis.text=element_text(size=12)
  #       , legend.title = element_text(size = 14, face = "bold")
  #       , legend.text = element_text(size=12)) +
  
  labs(x = "Month", y = "Frequency"
       # , title = "Year 2013 \nFrequency of Ile1016/Cys1534 Haplotype By Treatment Group"
       #, subtitle = "Error Bars = 95% CI"
       ) +
  scale_color_manual(name = "Treatment"
                     , values = c("Treatment" = "red", "Buffer" = "blue")
                     , labels = c("Spray Zone", "Buffer Zone")
                     , breaks=c("Treatment","Buffer")) +
  # Add yellow background to represent spray periods
  # Spray date range: 4/29/13 (4.96) - 6/3/13 (6.1)
  geom_rect(data = NULL, aes(xmin = 4.5, xmax = 5.5, ymin = -Inf, ymax = Inf)
            , fill="yellow", alpha = 0.05) +
  annotate("text", x=5, y=1.0, label= "Experimental", size = 5) +
  #Add data
  geom_point(data = mc.1016.t13, aes(color = "Treatment"), size = 5) +
  geom_line(data = mc.1016.t13[!is.na(mc.1016.t13$freqR),], color = "red", size = 2) +
  geom_errorbar(data = mc.1016.t13, aes(ymin=freqR-CI_95, ymax=freqR+CI_95), width=.2
                , color = " red", size = 0.7) +
  geom_point(data = mc.1016.b13, aes(color = "Buffer"), size = 5) +
  geom_line(data = mc.1016.b13[!is.na(mc.1016.b13$freqR),], color = "blue", size = 2) +
  geom_errorbar(data = mc.1016.b13, aes(ymin=freqR-CI_95, ymax=freqR+CI_95), width=.2
                , color = "blue", size = 0.7) +
  #Add n labels for each point
  #annotate("text", x=1, y=0.03, label = mc.1016.t13$n[1], color = "red", fontface = 2) +
  #annotate("text", x=2, y=0.03, label = mc.1016.t13$n[2], color = "red", fontface = 2) +
  #annotate("text", x=3, y=0.03, label = mc.1016.t13$n[3], color = "red", fontface = 2) +
  annotate("text", x=4, y=0.03, label = mc.1016.t13$n[3], color = "red", fontface = 2) +
  annotate("text", x=5, y=0.03, label = mc.1016.t13$n[4], color = "red", fontface = 2) +
  annotate("text", x=6, y=0.03, label = mc.1016.t13$n[5], color = "red", fontface = 2) +
  annotate("text", x=7, y=0.03, label = mc.1016.t13$n[6], color = "red", fontface = 2) +
  annotate("text", x=8, y=0.03, label = mc.1016.t13$n[7], color = "red", fontface = 2) +
  #annotate("text", x=9, y=0.03, label = mc.1016.t13$n[9], color = "red", fontface = 2) +
  #annotate("text", x=10, y=0.03, label = mc.1016.t13$n[10], color = "red", fontface = 2) +
  #annotate("text", x=1, y=0.0, label = mc.1016.b13$n[1], color = "blue", fontface = 2) +
  #annotate("text", x=2, y=0.0, label = mc.1016.b13$n[2], color = "blue", fontface = 2) +
  #annotate("text", x=3, y=0.0, label = mc.1016.b13$n[3], color = "blue", fontface = 2) +
  annotate("text", x=4, y=0.0, label = mc.1016.b13$n[3], color = "blue", fontface = 2) +
  annotate("text", x=5, y=0.0, label = mc.1016.b13$n[4], color = "blue", fontface = 2) +
  annotate("text", x=6, y=0.0, label = mc.1016.b13$n[5], color = "blue", fontface = 2) +
  annotate("text", x=7, y=0.0, label = mc.1016.b13$n[6], color = "blue", fontface = 2) +
  annotate("text", x=8, y=0.0, label = mc.1016.b13$n[7], color = "blue", fontface = 2) +
  #annotate("text", x=9, y=0.0, label = mc.1016.b13$n[9], color = "blue", fontface = 2) +
  #annotate("text", x=10, y=0.0, label = mc.1016.b13$n[10], color = "blue", fontface = 2) +

  ylim(c(0.0,1)) +
  #xlim(c(1,11)) +
  scale_x_continuous(breaks = pretty(mc.1016.t13$month, n = 10)) 


# View plot
kdrZones13
 
# Write plot to png
ggsave(filename = paste0("figures/kdrZones/kdrZones_2013/kdrZones13_", Sys.Date(), ".png"), width = 11, height = 8, dpi = 600, units = "in", device='png')
 
# 
# # Write plot to pdf
# pdf(file = paste("figures/kdrZones/kdrZones_", Sys.Date(), ".pdf", sep = ""), 11, 8.5)
# print(kdrZones)
# dev.off()
