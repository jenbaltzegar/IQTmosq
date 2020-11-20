#!/usr/bin/env Rscript

# This script will analyze & plot output from WFABC posteriors.
# Only the single loci under selection will be analyzed here.

### Set up the working space ----------------------------
# clear working environment
rm(list = ls())

# Set working directory
setwd("/home/gould/Documents/jen.temp/jfbaltz_kdr/WFABC")

# Load libraries
require(MASS)
library(matrixStats)


# Load data ----
post_s.1016 <- as.numeric(as.matrix(read.table("converted_V1016I_byMo_sel_posterior_s.txt")))
post_h.1016 <- as.numeric(as.matrix(read.table("converted_V1016I_byMo_sel_posterior_h.txt")))
post_s.1534 <- as.numeric(as.matrix(read.table("converted_F1534C_byMo_sel_posterior_s.txt")))
post_h.1534 <- as.numeric(as.matrix(read.table("converted_F1534C_byMo_sel_posterior_h.txt")))


# posterior stats at all loci ----
mean.post.s.1016 <- mean(post_s.1016)
mean.post.s.1534 <- mean(post_s.1534)
mean.post.h.1016 <- mean(post_h.1016)
mean.post.h.1534 <- mean(post_h.1534)
ci.s.1016 <- quantile(post_s.1016, c(0.025, 0.5, 0.975))
ci.s.1534 <- quantile(post_s.1534, c(0.025, 0.5, 0.975))
ci.h.1016 <- quantile(post_h.1016, c(0.025, 0.5, 0.975))
ci.h.1534 <- quantile(post_h.1534, c(0.025, 0.5, 0.975))
# view stats
# mean.post.s.1016; mean.post.s.1534; mean.post.h.1016; mean.post.h.1534
# ci.s.1016; ci.s.1534; ci.h.1016; ci.h.1534
df.1016 <- c(mean.post.s.1016, mean.post.h.1016, ci.s.1016[1], ci.s.1016[3], ci.h.1016[1], ci.h.1016[3])
df.1534 <- c(mean.post.s.1534, mean.post.h.1534, ci.s.1534[1], ci.s.1534[3], ci.h.1534[1], ci.h.1534[3])
df <- cbind(df.1016, df.1534)
colnames(df) <- c("Ile1016", "Cys1534")
rownames(df) <- c("s_mean", "h_mean", "s_CI_2.5", "s_CI_97.5", "h_CI_2.5", "h_CI_97.5")
# df
write.csv(df, "./results/converted/posteriorMeans_converted.csv")

# Plot Posteriors for Selection
png("./results/converted/plot_converted_posterior_s_1016_byMo_sel.png")
plot(density(t(post_s.1016)), lwd=2, main="Posteriors for Selection \n Ile1016 under selection", xlab="s")
dev.off()

png("./results/converted/plot_converted_posterior_s_1534_byMo_sel.png")
plot(density(t(post_s.1534)), lwd=2, main="Posteriors for Selection \n Cys1534 under selection", xlab="s")
dev.off()

# Plot the posteriors for Dominance
png("./results/converted/plot_converted_posterior_h_1016_byMo_sel.png")
plot(density(t(post_h.1016)), lwd=2 ,main="Posteriors for Dominance \n Ile1016 under selection", xlab="h")
dev.off()

png("./results/converted/plot_converted_posterior_h_1534_byMo_sel.png")
plot(density(t(post_h.1534)), lwd=2, main="Posteriors for Dominance \n Cys1534 under selection", xlab="h")
dev.off()

# Make temp dfs
post.1016 <- cbind(post_s.1016, post_h.1016)
post.1534 <- cbind(post_s.1534, post_h.1534)

# Plot the 2D posterior probabilities
z_1016 <- kde2d(post.1016[,1], post.1016[,2], n=300)
z_1534 <- kde2d(post.1534[,1], post.1534[,2], n=300)


# Save plots
png("./results/converted/plot_converted_2D_1016_byMo_sel.png")
image(z_1016, xlab="Selection Coefficient", ylab="Dominance", main="Ile1016 under selection"
      # , xlim = c(-1,1), ylim = c(0,1)
      )
dev.off()

png("./results/converted/plot_converted_2D_1534_byMo_sel.png")
image(z_1534, xlab="Selection Coefficient", ylab="Dominance", main="Cys1534 under selection"
      # , xlim = c(-1,1), ylim = c(0,1)
      )
dev.off()



# Create boxplot for loci
library(ggplot2)
library(reshape2)
library(gridExtra)
library(gridGraphics)

####
# Load data - same data, this time without converting to numeric
post_s.1016 <- read.table("converted_V1016I_byMo_sel_posterior_s.txt")
post_h.1016 <- read.table("converted_V1016I_byMo_sel_posterior_h.txt")
post_s.1534 <- read.table("converted_F1534C_byMo_sel_posterior_s.txt")
post_h.1534 <- read.table("converted_F1534C_byMo_sel_posterior_h.txt")

# Plotting for by Mo Under Selection
# post_s.1534 -----
boxplot_s.1534 <- ggplot(post_s.1534, aes("", V1)) +
  geom_boxplot() +
  coord_cartesian(ylim = c(0, 1)) +
  xlab("Locus 1534") +
  ylab("Selection Coefficient") +
  # ggtitle("Selection Coefficients") +
  theme(plot.title = element_text(hjust = 0.5)) 
boxplot_s.1534
ggsave("./results/converted/plot_converted_boxplot_s.1534.png", dpi = 600)

# post_s.1016 -----
boxplot_s.1016 <- ggplot(post_s.1016, aes("", V1)) +
  geom_boxplot() +
  coord_cartesian(ylim = c(0, 1)) +
  xlab("Locus 1016") +
  ylab("Selection Coefficient") +
  # ggtitle("Selection Coefficients") +
  theme(plot.title = element_text(hjust = 0.5)) 
boxplot_s.1016
ggsave("./results/converted/plot_converted_boxplot_s.1016.png", dpi = 600)

# post_h.1534 -----
boxplot_h.1534 <- ggplot(post_h.1534, aes("", V1)) +
  geom_boxplot() +
  coord_cartesian(ylim = c(0, 1)) +
  xlab("Locus 1534") +
  ylab("Dominance") +
  # ggtitle("Dominance") +
  theme(plot.title = element_text(hjust = 0.5)) 
boxplot_h.1534
ggsave("./results/converted/plot_converted_boxplot_h.1534.png", dpi = 600)

# post_h.1016 -----
boxplot_h.1016 <- ggplot(post_h.1016, aes("", V1)) +
  geom_boxplot() +
  coord_cartesian(ylim = c(0, 1)) +
  xlab("Locus 1016") +
  ylab("Dominance") +
  # ggtitle("Dominance") +
  theme(plot.title = element_text(hjust = 0.5)) 
boxplot_h.1016
ggsave("./results/converted/plot_converted_boxplot_h.1016.png", dpi = 600)

# Tiling the plots -----
boxplot_tile <- grid.arrange(boxplot_s.1016, boxplot_h.1016, boxplot_s.1534, boxplot_h.1534, nrow = 2)
ggsave(plot = boxplot_tile, "./results/converted/plot_converted_boxplot_tiled.png")




