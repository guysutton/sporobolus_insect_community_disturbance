##########################################################################################
##########################################################################################
##########################################################################################
#####        How do season and disturbance regime affect insect communities?   ###########
#####        Script 03: Model-based analyses (manyglm and LVM's)               ###########
#####        Script author: Guy F. Sutton                                      ###########
#####        Institution  : Department of Zoology and Entomology,              ###########  
#####                       Rhodes University,                                 ###########
#####                       South Africa                                       ###########
#####        Last update  : 09/05/2020                                         ###########
##########################################################################################
##########################################################################################
##########################################################################################

######
### - (1) Model-based hypothesis testing - manyglm 
######

# We need to convert the species abundance matrix to an 'mvabund' object
insect.spp <- mvabund(sim.data)

# Have a look at the data
par(mar=c(2,10,2,2)) # adjusts the margins
boxplot(sim.data,
        horizontal = TRUE,
        las=2, 
        main="Abundance")

# Check the mean-variance relationship
meanvar.plot(insect.spp,
             all.labels = FALSE,
             xlab = "Mean species abundance",
             ylab = "Variance")

# Run global model with poisson distribution
mod1 <- manyglm(insect.spp ~ com_data_lts$disturb_allow_overwinter *
                  com_data_lts$season, 
                family="poisson")

# Check model assumptions
plot(mod1)

# Here, we see a slight fan-shape in the model residuals. 
# - We need to model the variance differently.

# Plot mean-variance relationship directly
# Add a factor variable containing the disturb x season level combos
com_data_lts <- com_data_lts %>%
  unite("tr.block", disturb_allow_overwinter, season, remove = FALSE)

# Make plot
meanvar.plot(insect.spp ~ com_data_lts$tr.block,
             col = as.factor(com_data_lts$disturb_allow_overwinter),
             xlab = " Mean (log scale)",
             ylab = " Variance (log scale)",
             main = " ")

# Add linear increase line
# Poisson model data should fall on this line (linear increase in 
# var with mean)
abline(a=0, 
       b=1, 
       untf = TRUE, 
       lty=3)

# Add legend
legend("bottomright", 
       legend = c("Disturbed", "Undisturbed"), 
       col = c("red", "black"), 
       pch = c(1,1), 
       bty = "n", 
       pt.cex = 1.1, 
       cex = 1.1, 
       text.col = "black", 
       horiz = F)

# Clearly, the variance increases > mean, t.f. poisson model
# likely not appropriate. 
# *** Because our predictor here is orthogonal (factor), this is not NB. 

# Try negative-binomial model
mod2 <- manyglm(insect.spp ~ com_data_lts$disturb_allow_overwinter *
                             com_data_lts$season, 
                  family="negative_binomial")

# Check model assumptions
plot(mod2)

# Much better.

### Now we are going to test the multivariate hypothesis of whether 
### insect species assemblges varied according to disturbance regimes
### and/or season (to account for potential differences between seasons)/ 
### This analysis uses likelihood-ratio tests (LRT) and 
### resampled P-values to test significance.

anova(mod2)
anova(mod2, 
      nBoot=999, 
      test="wald")

### We can see from this table that there is a significant effect of disturbance
### (LRT = 67.56.82, P = 0.001), meaning that the species composition 
### of herbivores clearly differs between disturbance regimes

# Save the output of 'mvabund' global analysis as a word file
capture.output(anova(mod2, 
                     nBoot=999, 
                     test="wald"), 
               file="./results/mvabund_global_model_tests.doc")

### To examine this further, and see which insect species are more likely 
### to be found associated with the different disturbance regimes:
### We run univariate tests for each species separately.
### This is done by using the p.uni="adjusted" argument in the anova function. 
### The "adjusted" part of the argument refers to the resampling method used to compute the p values,
### taking into account the correlation between the response variables. 
### This correlation is often found in ecological systems where different species 
### will interact with each other, competing with or facilitating each others' resource use.
anova(mod2, 
      nBoot=999, 
      p.uni="adjusted")

# Save the output of 'mvabund' univariate analysis as a word file
capture.output(anova(mod2, p.uni="adjusted"), 
               file="./results/mvabund_univariate_tests.doc")

###########################################################
# - Boral R package - bayesian model based analyses
###########################################################

# Model based approach to unconstrained ordination
library(boral)

# Set defaults
example_mcmc_control <- list(n.burnin = 10000, 
                             n.iteration = 40000, 
                             n.thin = 30)

# Drop the dummy species 
sim.data <- sim.data[, 1:6]

# Define predictor variables 
x <- org.data %>%
  dplyr::select(disturb = disturb_allow_overwinter,
                season = season)
x <- as.matrix(x)
x

###
# - Fit poisson model
###

# Poisson model with fixed effects 
ins_p <- boral(sim.data, x = x,
                family = "poisson",
                lv.control = list(num.lv = 2),
                row.eff = "fixed",
                # Sampled 19 sites, each on six occasions (fixed effect)
                row.ids = matrix(rep(1:19, each=6)),
                save.model = TRUE,
                mcmc.control = example_mcmc_control)
summary(ins_p)

# Get the Geweke convergence statistic (z scores)
ins_p$geweke.diag
gew.pvals <- 2*pnorm(abs(unlist(ins_p$geweke.diag[[1]])), 
                     lower.tail = FALSE)
p.adjust(gew.pvals, method = "holm")

# Check model convergnece by MCMC using 'coda' package
# plot(get.mcmcsamples(ins_nb))

# Plot residuals
par(mfrow = c(2,2))
plot(ins_p, est = "median", jitter = FALSE)

# Save to disc
png("./figures/boral_poisson_diagnostics.png")
par(mfrow = c(2,2))
plot(ins_p, est = "median", jitter = TRUE)
dev.off()

###
# - Fit NB model
###

# NB model with fixed effects 
ins_nb <- boral(sim.data, x = x,
                family = "negative.binomial",
                lv.control = list(num.lv = 2),
                row.eff = "fixed",
                # Sampled 19 sites, each on six occasions (fixed effect)
                row.ids = matrix(rep(1:19, each=6)),
                save.model = TRUE,
                mcmc.control = example_mcmc_control)
summary(ins_nb)

# Get the Geweke convergence statistic (z scores)
ins_nb$geweke.diag
gew.pvals <- 2*pnorm(abs(unlist(ins_nb$geweke.diag[[1]])), 
                     lower.tail = FALSE)
p.adjust(gew.pvals, method = "holm")

# Check model convergnece by MCMC using 'coda' package
# plot(get.mcmcsamples(ins_nb))

# Plot residuals
par(mfrow = c(2,2))
plot(ins_nb, est = "mean", jitter = FALSE)

###
# - Plot unconstrained ordination 
###
# Plot ordination (can be compared to nMDS)
lvsplot(ins_nb)
summary(ins_nb)
plot(ins_nb)

# Extract locations of sites in ordination space
xy.bor <- ins_nb$lv.median
xy.bor <- as.data.frame(xy.bor)
xy.bor

# Add the categorical variables
xy.bor$disturb <- org.data$disturb_allow_overwinter
xy.bor$season <- org.data$season
xy.bor$site <- org.data$site_code
xy.bor

# Plot the nMDS with colour-coded points for disturbance regimes, 
# and shapes for season
ggplot(xy.bor, aes(lv1, 
                   lv2, 
                   colour = disturb)) + 
  geom_point(size = 2.5, 
             aes(shape = season)) +
  geom_jitter(width = 0.2, 
              height = 0.2) +
  labs(x = "LV1",
       y = "LV2",
       colour = "Disturbance",
       shape = "Season") +
  theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(colour = "black"),
        axis.title.x = element_text(margin = unit(c(2, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
        legend.position = "right")

ggsave("./figures/fig_boral_plot.png", 
       width = 6, 
       height = 4)

###
# - Plot with convex hulls around groups 
###

head(xy.bor)

# Calculate the hulls for each group
hull_boral <- xy.bor %>%
  group_by(disturb) %>%
  slice(chull(lv1, lv2))
hull_boral

# Plot the nMDS with colour-coded points for disturbance regimes, 
# shapes for season and convex hulls for disturance regimes 
ggplot(xy.bor, aes(lv1, 
                   lv2, 
                   colour = disturb,
                   fill = factor(disturb))) + 
  geom_point(size = 2.5, 
             aes(shape = season)) +
  geom_jitter(width = 0.2, 
              height = 0.2) +
  geom_polygon(data = hull_boral, alpha = 0.1) + 
  labs(x = "LV1",
       y = "LV2",
       colour = "Disturbance",
       shape = "Season") +
  theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(colour = "black"),
        axis.title.x = element_text(margin = unit(c(2, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
        legend.position = "right") +
  guides(fill = "none") +
  theme(legend.key=element_blank(),
        legend.background=element_blank())

ggsave("./figures/fig_boral_plot_with_ellipses.png", 
       width = 6, 
       height = 4)
