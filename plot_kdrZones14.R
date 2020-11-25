### Plot Frequency of Resistance Allele at 1016 locus by zone for year 2014
# Plot based off of dataframes mc.1016.t and mc.1016.b

#Load Library
library(ggplot2)

# set working directory
# setwd("/Users/jenbaltz/Dropbox/GouldLab/Project_Mosquito/Database")

# Read in dataframes
mc.1016.t <- read.csv("./data/mc.1016.t.csv")
# mc.1016.b <- read.csv("mc.1016.b.csv")

# 3/19/18 - Subset mc.1016.t to remove March data because sample size is too low currently n = 4
mc.1016.t <- mc.1016.t[-3,]

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
kdrZones14 <- ggplot(mc.1016.t, aes(x=month, y=freqR)) +
  my_theme() +
  labs(x = "Month", y = ""
       # , y = "Frequency"
       # , title = "Year 2014 \n Frequency of Ile1016/Cys1534 Haplotype By Treatment Group"
       #, subtitle = "Error Bars = 95% CI"
       ) +
  scale_color_manual(name = "Treatment"
                     , values = c("Treatment" = "red", "Buffer" = "blue")
                     , labels = c("Spray Zone", "Buffer Zone")
                     , breaks=c("Treatment","Buffer")) +
  #Add yellow background to represent spray periods
     # experimental: 4/28/14 - 6/2/14
  geom_rect(data = NULL, aes(xmin = 1.95, xmax = 2.05, ymin = -Inf, ymax = Inf), fill="yellow", alpha = 0.05) +
  annotate("text", x=2, y=1.0, label= "City Wide", size = 5) +
  geom_rect(data = NULL, aes(xmin = 4.6, xmax = 5.4, ymin = -Inf, ymax = Inf), fill="yellow", alpha = 0.05) +
  annotate("text", x=5, y=1.0, label= "Experimental", size = 5) +
  #Add data
  geom_point(data = mc.1016.t, aes(color = "Treatment"), size = 5) +
  geom_line(data = mc.1016.t[!is.na(mc.1016.t$freqR),], color = "red", size = 2) +
  geom_errorbar(data = mc.1016.t, aes(ymin=freqR-CI_95, ymax=freqR+CI_95), width=.2
                , color = " red", size = 0.7) +
  
  geom_point(data = mc.1016.b, aes(color = "Buffer"), size = 5) +
  geom_line(data = mc.1016.b[!is.na(mc.1016.b$freqR),], color = "blue", size = 2) +
  geom_errorbar(data = mc.1016.b, aes(ymin=freqR-CI_95, ymax=freqR+CI_95), width=.2
                , color = "blue", size = 0.7) +
  #geom_smooth(data = mc.1016.t[!is.na(mc.1016.t$freqR),], method='glm', color = "red") +
  #geom_smooth(data = mc.1016.b[!is.na(mc.1016.b$freqR),], method = "lm", color = "blue") +
  
#Add n labels for each point
  annotate("text", x=1, y=0.03, label = mc.1016.t$n[1], color = "red", fontface = 2) +
  annotate("text", x=2, y=0.03, label = mc.1016.t$n[2], color = "red", fontface = 2) +
  annotate("text", x=4, y=0.03, label = mc.1016.t$n[3], color = "red", fontface = 2) +
  annotate("text", x=5, y=0.03, label = mc.1016.t$n[4], color = "red", fontface = 2) +
  annotate("text", x=6, y=0.03, label = mc.1016.t$n[5], color = "red", fontface = 2) +
  annotate("text", x=7, y=0.03, label = mc.1016.t$n[6], color = "red", fontface = 2) +
  annotate("text", x=8, y=0.03, label = mc.1016.t$n[7], color = "red", fontface = 2) +
  annotate("text", x=10, y=0.03, label = mc.1016.t$n[9], color = "red", fontface = 2) +
  annotate("text", x=1, y=0.0, label = mc.1016.b$n[1], color = "blue", fontface = 2) +
  annotate("text", x=2, y=0.0, label = mc.1016.b$n[2], color = "blue", fontface = 2) +
  annotate("text", x=4, y=0.0, label = mc.1016.b$n[4], color = "blue", fontface = 2) +
  annotate("text", x=5, y=0.0, label = mc.1016.b$n[5], color = "blue", fontface = 2) +
  annotate("text", x=6, y=0.0, label = mc.1016.b$n[6], color = "blue", fontface = 2) +
  annotate("text", x=7, y=0.0, label = mc.1016.b$n[7], color = "blue", fontface = 2) +
  annotate("text", x=8, y=0.0, label = mc.1016.b$n[8], color = "blue", fontface = 2) +
  annotate("text", x=10, y=0.0, label = mc.1016.b$n[10], color = "blue", fontface = 2) +

  ylim(c(0.0,1)) +
  #xlim(c(1,11)) +
  scale_x_continuous(breaks = pretty(mc.1016.t$month, n = 10)) 


# View plot
kdrZones14
 
# Write plot to png
# ggsave(filename = paste0("figures/kdrZones/kdrZones_2014/kdrZones14_", Sys.Date(), ".png"), width = 11, height = 8, dpi = 600, units = "in", device='png')
ggsave(filename = paste0("./data/kdrZones14_", Sys.Date(), ".png"), width = 11, height = 8, dpi = 600, units = "in", device='png')

# 
# # Write plot to pdf
# pdf(file = paste("figures/kdrZones/kdrZones_", Sys.Date(), ".pdf", sep = ""), 11, 8.5)
# print(kdrZones)
# dev.off()
