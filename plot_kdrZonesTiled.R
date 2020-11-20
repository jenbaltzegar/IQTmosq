#Combine kdrZones13 and kdrZones14 plots for publication

# load library
library(gridExtra)

# source objects
source("./R_scripts/IQTmosq/plot_kdrZones13.R")
source("./R_scripts/IQTmosq/plot_kdrZones14.R")

# check plots and remove legend from 2013
kdrZones13 <- kdrZones13 + theme(legend.position = "none")
kdrZones14

# make tiled plot
zones_tiled <- grid.arrange(kdrZones13, kdrZones14, nrow = 1, widths = c(0.9, 1))
ggsave(plot = zones_tiled, filename = paste0("figures/kdrZones/kdrZonesTiled_", Sys.Date(), ".png")
       , width = 11, height = 8, dpi = 600, units = "in"
       , device='png')
