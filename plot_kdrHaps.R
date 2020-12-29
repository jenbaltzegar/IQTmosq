# Plot the frequency of kdr haplotypes over time with ggplot()
# depends on setup_jb.R and run_Haplotype_Imputation.R

## subset and modify
freq_plot <- subset(freqAll_long, !is.na(Frequency))
# reassign Haplotypes to have specific text for legend
freq_plot$Haplotype  <- factor(freq_plot$Haplotype,
    levels = c('SS','SR', 'RS', 'RR'),
    labels = c("Val1016/Phe1534", "Val1016/Cys1534", "Ile1016/Phe1534", "Ile1016/Cys1534")
)
 
# Create custom color palettes -----
hapColors <- c("#117733", "#332288", "#CC6677", "#882255")
#names(hapColors) <- rev(levels(freq_plot$Haplotype))
colScale <- scale_colour_manual(name = "Haplotype",values = hapColors)

sprayColors <- c("#A0A0A0", "#D1D1D1","#BBBBBB", "#E6E6E6", "#7F7F7F", "#FFFFFF")

# Plot kdr Haplotypes -----
kdrHaps <- ggplot(data = freq_plot, aes(x=year, y=Frequency)) +
  my_theme +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + # modify theme
  labs(x = "Year", y = "Frequency"
       #, title = "kdr Haplotype Frequency"
       #, subtitle = option
       ) +
  scale_y_continuous(breaks=c(0.0,0.2, 0.4, 0.6, 0.8, 1.0), expand=expansion(add=0.11)) +
  scale_x_continuous(breaks = pretty(freq_plot$year, n = 18), expand=expansion(0,0)) +
  coord_cartesian(clip = 'off') +
  ### Add background to represent type of insecticides used 
  annotate("rect", xmin = 1999.0, xmax = 2001.5, ymin = -Inf, ymax = Inf, alpha=0.15, fill="grey") + # No pyrethroids
  annotate("rect", xmin = 2001.5, xmax = 2007.5, ymin = -Inf, ymax = Inf, alpha=0.15, fill="red") + # Delatmethrin 2002 - 2007
  annotate("rect", xmin = 2007.5, xmax = 2008.5, ymin = -Inf, ymax = Inf, alpha=0.15, fill="orange") + # Cypermethrin 2008
  annotate("rect", xmin = 2008.5, xmax = 2011.5, ymin = -Inf, ymax = Inf, alpha=0.15, fill="yellow") + # a-Cypermethrin 2009 - 2011
  annotate("rect", xmin = 2011.5, xmax = 2013.5, ymin = -Inf, ymax = Inf, alpha=0.15, fill="salmon4") + # Cypermethrin & a-Cypermethrin 2012 # Lambda, a-Cypermethrin, Cypermethrin, & a-Cypermethrin+pyriporx 2013
  annotate("rect", xmin = 2013.5, xmax = 2014.5, ymin = -Inf, ymax = Inf, alpha=0.15, fill="orange") + # Cypermethrin 2014
  annotate("rect", xmin = 2014.5, xmax = 2017.5, ymin = -Inf, ymax = Inf, alpha=0.15, fill="gray") + # No pyrethroids
  
  ### Add text to represent type of insecticides used
  annotate("text", y=1.1, size=4,
    x=c(2000.25, 2004.5, 2008, 2010, 2012.5, 2014, 2016),
    label=c("No \n Pyrethroids", "Delta",  "Cyper", "Alpha", "Multiple*",  "Cyper",  "No \n Pyrethroids")
  ) +
  annotate("text", y=-0.09, x=c(1999.5, 2000:2017), label = c(
'N=',mc.haps.yr$n), fontface = 2, color = "#332288") +
  
  ### Add data
  geom_line(aes(colour=Haplotype), size = 1.1) +
  geom_point(aes(colour=Haplotype), size = 3) +
  geom_errorbar(data = freq_plot, aes(ymin=Frequency-CI_95, ymax=Frequency+CI_95, colour=Haplotype), width=.14, size = 0.7) +
  ## 2-column legend on top
  theme(
    legend.position='top', 
    legend.direction='horizontal'
  ) + 
  guides(color=guide_legend(ncol=2)) +
  ### Add Legends
  colScale 
  
# Print graph to screen
#print(kdrHaps)
