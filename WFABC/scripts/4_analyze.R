#!/usr/bin/env Rscript

# This script will analyze output from the run_WFABC.sh program

### Set up the working space ----------------------------
# clear working environment
rm(list = ls())

# Set working directory
setwd("~/Documents/jen.temp/jfbaltz_kdr/WFABC")

# Load libraries
require(MASS)
library(matrixStats)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(gridGraphics)


# Load Data -----
# NOT converted
# s.1016 <- as.matrix(read.table("V1016I_byMo_posterior_s.txt"))
# h.1016 <- as.matrix(read.table("V1016I_byMo_posterior_h.txt"))
# s.1534 <- as.matrix(read.table("F1534C_byMo_posterior_s.txt"))
# h.1534 <- as.matrix(read.table("F1534C_byMo_posterior_h.txt"))
# s.1016 <- as.matrix(read.table("V1016I_byMo_sel_posterior_s.txt"))
# h.1016 <- as.matrix(read.table("V1016I_byMo_sel_posterior_h.txt"))
# s.1534 <- as.matrix(read.table("F1534C_byMo_sel_posterior_s.txt"))
# h.1534 <- as.matrix(read.table("F1534C_byMo_sel_posterior_h.txt"))

# CONVERTED
# s.1016 <- as.matrix(read.table("converted_V1016I_byMo_posterior_s.txt"))
# h.1016 <- as.matrix(read.table("converted_V1016I_byMo_posterior_h.txt"))
# s.1534 <- as.matrix(read.table("converted_F1534C_byMo_posterior_s.txt"))
# h.1534 <- as.matrix(read.table("converted_F1534C_byMo_posterior_h.txt"))
s.1016 <- as.matrix(read.table("converted_V1016I_byMo_sel_posterior_s.txt"))
h.1016 <- as.matrix(read.table("converted_V1016I_byMo_sel_posterior_h.txt"))
s.1534 <- as.matrix(read.table("converted_F1534C_byMo_sel_posterior_s.txt"))
h.1534 <- as.matrix(read.table("converted_F1534C_byMo_sel_posterior_h.txt"))


# Calculate stats of posteriors -----
mean.s.1016 <- mean(s.1016)
mean.s.1534 <- mean(s.1534)
mean.h.1016 <- mean(h.1016)
mean.h.1534 <- mean(h.1534)
ci.s.1016 <- quantile((s.1016), c(0.025, 0.5, 0.975))
ci.s.1534 <- quantile((s.1534), c(0.025, 0.5, 0.975))
ci.h.1016 <- quantile((h.1016), c(0.025, 0.5, 0.975))
ci.h.1534 <- quantile((h.1534), c(0.025, 0.5, 0.975))
# view stats
# mean.s.1016; mean.s.1534; mean.h.1016; mean.h.1534
# ci.s.1016; ci.s.1534; ci.h.1016; ci.h.1534
# create df of results
df.1016 <- c(ci.s.1016[1], mean.s.1016, ci.s.1016[3], ci.h.1016[1], mean.h.1016, ci.h.1016[3])
df.1534 <- c(ci.s.1534[1], mean.s.1534, ci.s.1534[3], ci.h.1534[1], mean.h.1534, ci.h.1534[3])
df <- cbind(df.1016, df.1534)
colnames(df) <- c("Ile1016", "Cys1534")
rownames(df) <- c("s_CI_2.5", "s_mean", "s_CI_97.5", "h_CI_2.5", "h_mean", "h_CI_97.5")
# df
write.csv(df, "posteriorMeans.csv")

# Plot the posteriors for s -----
png("./plot_posterior_s_1016.png")
plot(density(t(s.1016[1,])), lwd=2, main="Posteriors for Selection \n Ile1016", xlab="s")
dev.off()

png("./plot_posterior_s_1534.png")
plot(density(t(s.1534[1,])), lwd=2, main="Posteriors for Selection \n Cys1534", xlab="s")
dev.off()

# Plot the posteriors for h -----
png("./plot_posterior_h_1016.png")
plot(density(t(h.1016[1,])), lwd=2 ,main="Posteriors for Dominance \n Ile1016", xlab="h")
dev.off()

png("./plot_posterior_h_1534.png")
plot(density(t(h.1534[1,])), lwd=2, main="Posteriors Dominance \n Cys1534", xlab="h")
dev.off()


# Plot the 2D posterior probabilities -----
z_1016 <- kde2d(s.1016[1,], h.1016[1,], n=300)
z_1534 <- kde2d(s.1534[1,], h.1534[1,], n=300)

# Save plots
png("./plot_2D_1016.png")
image(z_1016,xlab="Selection Coefficient",ylab="Dominance", main="Ile1016")
dev.off()

png("./plot_2D_1534.png")
image(z_1534,xlab="Selection Coefficient",ylab="Dominance", main="Cys1534")
dev.off()




# Create boxplot for loci -----
# Selection Coefficient plots
cat.s <- as.data.frame(rbind(s.1016, s.1534))
cat.s$group <- row.names(cat.s)
post_s.m <- melt(cat.s, id.vars = "group")

boxplot_s <- ggplot(post_s.m, aes(group, value)) + 
  geom_boxplot() +
  xlab("Locus") +
  ylab("Selection Coefficient") +
  ggtitle("Selection Coefficients") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels=c("1" = "Ile1016", "2" = "Cys1534"))
boxplot_s
ggsave("./boxplot_s.png", dpi = 600)

violin_s <- ggplot(post_s.m, aes(group, value)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  xlab("Locus") +
  ylab("Selection Coefficient") +
  ggtitle("Selection Coefficients") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels=c("1" = "Ile1016", "2" = "Cys1534"))
violin_s
ggsave("./violin_s.png", dpi = 600)

# Dominance plots
cat.h <- as.data.frame(rbind(h.1016, h.1534))
cat.h$group <- row.names(cat.h)
post_h.m <- melt(cat.h, id.vars = "group")

boxplot_h <- ggplot(post_h.m, aes(group, value)) + 
  geom_boxplot() +
  xlab("Locus") +
  ylab("Dominance") +
  ggtitle("Dominance") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels=c("1" = "Ile1016", "2" = "Cys1534"))
boxplot_h
ggsave("./boxplot_h.png", dpi = 600)

violin_h <- ggplot(post_h.m, aes(group, value)) + 
  geom_violin(draw_quantiles = c(0.25, 0.5, 0.75)) +
  xlab("Locus") +
  ylab("Dominance") +
  ggtitle("Dominance") +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels=c("1" = "Ile1016", "2" = "Cys1534"))
violin_h
ggsave("./violin_h.png", dpi = 600)

# Tiled plots
boxplot_tile <- grid.arrange(boxplot_s, boxplot_h, nrow = 1)
ggsave(plot = boxplot_tile, "./boxplot_tiled.png")

violin_tile <- grid.arrange(violin_s, violin_h, nrow = 1)
ggsave(plot = boxplot_tile, "./violin_tiled.png")


