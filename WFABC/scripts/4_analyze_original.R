#!/usr/bin/env Rscript

# This script will analyze output from the run_WFABC.sh program

### Set up the working space ----------------------------
# clear working environment
rm(list = ls())

# Set working directory
setwd("/Users/jenbaltz/Dropbox/GouldLab/Project_Mosquito/Database/WFABC")

# Load libraries
require(MASS)
library(matrixStats)

### Do some stuff ------------------------------------
# Plot the posteriors for s
post_s <- read.table("WFABC_multiple_loci_posterior_s.txt")
post_s_1016 <- read.table("WFABC_V1016I_underSelection_posterior_s.txt")
post_s_1534 <- read.table("WFABC_F1534C_underSelection_posterior_s.txt")

png("./plot_posterior_s_1016.png")
plot(density(t(post_s[1,])), lwd=2, main="Posteriors for Selection \n Ile1016", xlab="s")
dev.off()

png("./plot_posterior_s_1534.png")
plot(density(t(post_s[2,])), lwd=2, main="Posteriors for Selection \n Cys1534", xlab="s")
dev.off()

png("./plot_posterior_s_1016_selection.png")
plot(density(t(post_s_1016[1,])), lwd=2, main="Posteriors for Selection \n Ile1016 under selection", xlab="s")
dev.off()

png("./plot_posterior_s_1534_selection.png")
plot(density(t(post_s_1534[1,])), lwd=2, main="Posteriors for Selection \n Cys1534 under selection", xlab="s")
dev.off()

# Plot the posteriors for h
post_h <- read.table("WFABC_multiple_loci_posterior_h.txt")
post_h_1016 <- read.table("WFABC_V1016I_underSelection_posterior_h.txt")
post_h_1534 <- read.table("WFABC_F1534C_underSelection_posterior_h.txt")

png("./plot_posterior_h_1016.png")
plot(density(t(post_h[1,])), lwd=2 ,main="Posteriors for Dominance \n Ile1016", xlab="h")
dev.off()

png("./plot_posterior_h_1534.png")
plot(density(t(post_h[2,])), lwd=2, main="Posteriors Dominance \n Cys1534", xlab="h")
dev.off()

png("./plot_posterior_h_1016_selection.png")
plot(density(t(post_h_1016[1,])), lwd=2 ,main="Posteriors for Dominance \n Ile1016 under selection", xlab="h")
dev.off()

png("./plot_posterior_h_1534_selection.png")
plot(density(t(post_h_1534[1,])), lwd=2, main="Posteriors Dominance \n Cys1534 under selection", xlab="h")
dev.off()


# Plot the 2D posterior probabilities
z_1016 <- kde2d(t(post_s[1,]),t(post_h[1,]),n=300)
z_1534 <- kde2d(t(post_s[2,]),t(post_h[2,]),n=300)
z_1016_sel <- kde2d(t(post_s_1016[1,]),t(post_h_1016[1,]),n=300)
z_1534_sel <- kde2d(t(post_s_1534[1,]),t(post_h_1534[1,]),n=300)

# Save plots
png("./plot_2D_1016.png")
image(z_1016,xlab="Selection Coefficient",ylab="Dominance", main="Ile1016")
dev.off()

png("./plot_2D_1534.png")
image(z_1534,xlab="Selection Coefficient",ylab="Dominance", main="Cys1534")
dev.off()

png("./plot_2D_1016_selection.png")
image(z_1016_sel,xlab="Selection Coefficient",ylab="Dominance", main="Ile1016 under selection")
dev.off()

png("./plot_2D_1534_selection.png")
image(z_1534_sel,xlab="Selection Coefficient",ylab="Dominance", main="Cys1534 under selection")
dev.off()


# posterior stats at all loci
mean.post.s <- rowMeans(post_s)
mean.post.s.1016.sel <- rowMeans(post_s_1016)
mean.post.s.1534.sel <- rowMeans(post_s_1534)
mean.post.h <- rowMeans(post_h)
mean.post.h.1016.sel <- rowMeans(post_h_1016)
mean.post.h.1534.sel <- rowMeans(post_h_1534)
sds.post.s <- rowSds(as.matrix(post_s))
sds.post.s.1016.sel <- rowSds(as.matrix(post_s_1016))
sds.post.s.1534.sel <- rowSds(as.matrix(post_s_1016))
sds.post.h <- rowSds(as.matrix(post_h))
sds.post.h.1016.sel <- rowSds(as.matrix(post_h_1016))
sds.post.h.1534.sel <- rowSds(as.matrix(post_h_1534))
range.post.s <- rowRanges(as.matrix(post_s))
range.post.s.1016.sel <- rowRanges(as.matrix(post_s_1016))
range.post.s.1534.sel <- rowRanges(as.matrix(post_s_1534))
range.post.h <- rowRanges(as.matrix(post_h))
range.post.h.1016.sel <- rowRanges(as.matrix(post_h_1016))
range.post.h.1534.sel <- rowRanges(as.matrix(post_h_1534))

# Create and save means
mean.post.s.sel <- cbind(mean.post.s.1016.sel, mean.post.s.1534.sel)
mean.post.h.sel <- cbind(mean.post.h.1016.sel, mean.post.h.1534.sel)
sds.post.s.sel <- cbind(sds.post.s.1016.sel, sds.post.s.1534.sel)
sds.post.h.sel <- cbind(sds.post.h.1016.sel, sds.post.h.1534.sel)

df <- rbind(mean.post.s, sds.post.s, mean.post.s.sel, sds.post.s.sel
            , mean.post.h, sds.post.h, mean.post.h.sel, sds.post.h.sel)
# df <- rbind(mean.post.s, sds.post.s, mean.post.h, sds.post.h)
colnames(df) <- c("Ile1016", "Cys1534")
rownames(df) <- c("mean.post.s", "sds.post.s", "mean.post.s.sel"
                  , "sds.post.s.sel", "mean.post.h", "sds.post.h"
                  , "mean.post.h.sel", "sds.post.h.sel")
write.csv(df, "posteriorMeans.csv")

# Create boxplot for loci
library(ggplot2)
library(reshape2)
library(gridExtra)
library(gridGraphics)

post_s$group <- row.names(post_s)
post_s.m <- melt(post_s, id.vars = "group")
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

post_h$group <- row.names(post_h)
post_h.m <- melt(post_h, id.vars = "group")
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

boxplot_tile <- grid.arrange(boxplot_s, boxplot_h, nrow = 1)
ggsave(plot = boxplot_tile, "./boxplot_tiled.png")

violin_tile <- grid.arrange(violin_s, violin_h, nrow = 1)
ggsave(plot = boxplot_tile, "./violin_tiled.png")


