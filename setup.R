# Load Libraries
library(ggplot2)
library(ggthemes)
library(wrapr)

# create custom theme - from Michael Vella -----
my_theme <- (
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
)
