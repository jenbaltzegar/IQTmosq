#!/usr/bin/env Rscript

# This script will analyze & plot output from WFABC posteriors

### Set up the working space ----------------------------
# clear working environment
rm(list = ls())

# Set working directory
setwd("/home/gould/Desktop/WFABC")

# Load libraries
require(MASS)
library(matrixStats)

# Convert the posteriors
source("./scripts/function_WFABC_convertPosteriors.R")
convertPosteriors("multiple_loci_byMo_sel_posterior_s.txt", "multiple_loci_byMo_sel_posterior_h.txt")

# Load data
post_s <- read.table("converted_multiple_loci_byMo_sel_posterior_s.txt")
post_h <- read.table("converted_multiple_loci_byMo_sel_posterior_h.txt")

# Plot Posteriors for Selection
png("./plot_converted_posterior_s_1016_byMo_sel.png")
plot(density(t(post_s[1,])), lwd=2, main="Posteriors for Selection \n Ile1016 under selection", xlab="s")
dev.off()

png("./plot_converted_posterior_s_1534_byMo_sel.png")
plot(density(t(post_s[2,])), lwd=2, main="Posteriors for Selection \n Cys1534 under selection", xlab="s")
dev.off()

# Plot the posteriors for Dominance

png("./plot_converted_posterior_h_1016_byMo_sel.png")
plot(density(t(post_h[1,])), lwd=2 ,main="Posteriors for Dominance \n Ile1016 under selection", xlab="h")
dev.off()

png("./plot_converted_posterior_h_1534_byMo_sel.png")
plot(density(t(post_h[2,])), lwd=2, main="Posteriors Dominance \n Cys1534 under selection", xlab="h")
dev.off()


# Plot the 2D posterior probabilities
z_1016 <- kde2d(t(post_s[1,]),t(post_h[1,]),n=300)
z_1534 <- kde2d(t(post_s[2,]),t(post_h[2,]),n=300)

# Save plots
png("./plot_2D_1016_byMo_sel.png")
image(z_1016, xlab="Selection Coefficient 2002-2009", ylab="Dominance", main="Ile1016 under selection")
dev.off()

png("./plot_2D_1534_byMo_sel.png")
image(z_1534, xlab="Selection Coefficient 2002-2009", ylab="Dominance", main="Cys1534 under selection")
dev.off()

# posterior stats at all loci
mean.post.s <- rowMeans(post_s)
mean.post.h <- rowMeans(post_h)
sds.post.s <- rowSds(as.matrix(post_s))
sds.post.h <- rowSds(as.matrix(post_h))
range.post.s <- rowRanges(as.matrix(post_s))
range.post.h <- rowRanges(as.matrix(post_h))
ci.s.1534 <- quantile(post_s[2,], c(0.025, 0.5, 0.975))
ci.h.1534 <- quantile(post_h[2,], c(0.025, 0.5, 0.975))


df <- rbind(mean.post.s, sds.post.s, mean.post.h, sds.post.h)
colnames(df) <- c("Ile1016", "Cys1534")
write.csv(df, "posteriorMeans.csv")

# Create boxplot for loci
library(ggplot2)
library(reshape2)
library(gridExtra)
library(gridGraphics)

####
# Plotting for by Mo Under Selection -----------
post_s$group <- row.names(post_s)
post_s.m <- melt(post_s, id.vars = "group")
boxplot_s <- ggplot(post_s.m, aes(group, value)) + 
  geom_boxplot() +
  xlab("Locus") +
  ylab("Selection Coefficient") +
  ggtitle("Selection Coefficients 2002-2009") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels=c("1" = "Ile1016", "2" = "Cys1534"))
boxplot_s
ggsave("./boxplot_s.png", dpi = 600)
violin_s <- ggplot(post_s.m, aes(group, value)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  xlab("Locus") +
  ylab("Selection Coefficient") +
  ggtitle("Selection Coefficients 2002-2009") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels=c("1" = "Ile1016", "2" = "Cys1534"))
violin_s
ggsave("./violin_s.png", dpi = 600)

post_h$group <- row.names(post_h)
post_h.m <- melt(post_h, id.vars = "group")
boxplot_h <- ggplot(post_h.m, aes(group, value)) + 
  geom_boxplot() +
  xlab("Locus") +
  ylab("Dominance") +
  ggtitle("Dominance 2002-2009") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels=c("1" = "Ile1016", "2" = "Cys1534"))
boxplot_h
ggsave("./boxplot_h.png", dpi = 600)
violin_h <- ggplot(post_h.m, aes(group, value)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  xlab("Locus") +
  ylab("Dominance") +
  ggtitle("Dominance 2002-2009") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels=c("1" = "Ile1016", "2" = "Cys1534"))
violin_h
ggsave("./violin_h.png", dpi = 600)

boxplot_tile <- grid.arrange(boxplot_s, boxplot_h, nrow = 1)
ggsave(plot = boxplot_tile, "./boxplot_tiled.png")

violin_tile <- grid.arrange(violin_s, violin_h, nrow = 1)
ggsave(plot = violin_tile, "./violin_tiled.png")



