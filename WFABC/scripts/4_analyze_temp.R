#!/usr/bin/env Rscriptp

# This script will analyze & plot output from WFABC posteriors that have been converted.
# Only the single loci under selection will be analyzed here.

# ### Convert WFABC posteriors to match Michael Vella's scale
# convertPosteriors.onelocus("V1016I_byMo_sel_posterior_s.txt", "V1016I_byMo_sel_posterior_h.txt")
# convertPosteriors.onelocus("F1534C_byMo_sel_posterior_s.txt", "F1534C_byMo_sel_posterior_h.txt")

# Load data ----
post_s.1016 <- unlist(read.table("V1016I_byMo_sel_posterior_s.txt"))
post_h.1016 <- unlist(read.table("V1016I_byMo_sel_posterior_h.txt"))
post_s.1534 <- unlist(read.table("F1534C_byMo_sel_posterior_s.txt"))
post_h.1534 <- unlist(read.table("F1534C_byMo_sel_posterior_h.txt"))


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
write.csv(df, "./results/posteriorMeans_converted.csv")

# Plot Posteriors for Selection
#png("./results/plot_posterior_s_1016_byMo_sel.png")
#plot(density(t(post_s.1016)), lwd=2, main="Posteriors for Selection \n Ile1016 under selection", xlab="s")
#dev.off()

#png("./results/plot_posterior_s_1534_byMo_sel.png")
#plot(density(t(post_s.1534)), lwd=2, main="Posteriors for Selection \n Cys1534 under selection", xlab="s")
#dev.off()

# Plot the posteriors for Dominance
#png("./results/plot_posterior_h_1016_byMo_sel.png")
#plot(density(t(post_h.1016)), lwd=2 ,main="Posteriors for Dominance \n Ile1016 under selection", xlab="h")
#dev.off()

#png("./results/plot_posterior_h_1534_byMo_sel.png")
#plot(density(t(post_h.1534)), lwd=2, main="Posteriors for Dominance \n Cys1534 under selection", xlab="h")
#dev.off()

# Make temp dfs
post.1016 <- data.frame(s=post_s.1016, h=post_h.1016)
post.1534 <- data.frame(s=post_s.1534, h=post_h.1534)

dens.1534 <- (
    ggplot(post.1534, aes(x=abs(s), y=h))
    + geom_point(alpha=0.2)
    + geom_density_2d()
    + xlim(0,1.0)
    + ylim(0,1.0)
    + xlab('Selection')
    + ylab('Dominance')
    + my_theme
)

dens.1016 <- dens.1534 %+% post.1016

# Plot the 2D posterior probabilities
#z_1016 <- kde2d(post.1016[,1], post.1016[,2], n=300)
#z_1534 <- kde2d(post.1534[,1], post.1534[,2], n=300)

# Save plots
#png("./results/plot_2D_1016_byMo_sel.png")
#image(z_1016, xlab="Selection Coefficient", ylab="Dominance", main="Ile1016 under selection"
      # , xlim = c(-1,1), ylim = c(0,1)
#      )# path=paste0("../programs/WFABC_v1.1/binaries/", myos, "/")

#dev.off()

#png("./results/plot_2D_1534_byMo_sel.png")
#image(z_1534, xlab="Selection Coefficient", ylab="Dominance", main="Cys1534 under selection"
      # , xlim = c(-1,1), ylim = c(0,1)
#      )
#dev.off()

# Plotting for by Mo Under Selection
# post_s.1534 -----
boxplot_s.1534 <- ggplot(post.1534, aes(y=abs(s), x='')) +
  geom_boxplot() +
  coord_cartesian(ylim = c(0, 1)) +
  #xlab("Locus 1534") +
  xlab("") +
  ylab("Selection") +
  # ggtitle("Selection Coefficients") +
  my_theme

## same for 1016
boxplot_s.1016 <- boxplot_s.1534 %+% post.1016

# post_h.1534 -----
boxplot_h.1534 <- ggplot(post.1534, aes(y=h, x=''))+
  geom_boxplot() +
  coord_cartesian(ylim = c(0, 1)) +
  xlab("") +
  ylab("Dominance") +
  my_theme

boxplot_h.1016 <- boxplot_h.1534 %+% post.1016

# Tiling the plots -----
#boxplot_tile <- grid.arrange(boxplot_s.1016, boxplot_h.1016, boxplot_s.1534, boxplot_h.1534, nrow = 2)
#ggsave(plot = boxplot_tile, "./results/plot_boxplot_tiled.png")




