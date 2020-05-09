##########################################################################################
##########################################################################################
##########################################################################################
#####        How do season and disturbance regime affect insect communities?   ###########
#####        Script 01: Data import and cleaning                               ###########
#####        Script author: Guy F. Sutton                                      ###########
#####        Institution  : Department of Zoology and Entomology,              ###########  
#####                       Rhodes University,                                 ###########
#####                       South Africa                                       ###########
#####        Last update  : 24/10/2018                                         ###########
##########################################################################################
##########################################################################################
##########################################################################################

##########################################################################################
#####                            Load required packages                             ######
##########################################################################################

library(vegan)
library(MASS)
library(mvabund)
library(tidyr)
library(readr)
library(tidyverse)

##########################################################################################
#####                                 Import data                                   ######
##########################################################################################

# Import raw data 
raw_data <- read_csv2("./data_raw/all_surveys_data.csv")

# Check the data 
head(raw_data)
View(raw_data)

# Drop some columns that we don't need 
raw_data <- raw_data %>%
  dplyr::select(site:insect_abun_total, disturbance_type:disturb_allow_overwinter) %>%
  dplyr::filter(plant_species == c("Sporobolus pyramidalis", "Sporobolus natalensis")) %>%
  mutate(disturb_allow_overwinter = case_when(
    disturb_allow_overwinter == "No" ~ "Disturbed",
    disturb_allow_overwinter == "Yes" ~ "Undisturbed"
  ))
head(raw_data)

# Summarise site x sampling date abundance data
com_data <- raw_data %>%
  group_by(site_code, sampling_date) %>%
  mutate(
    sp_tet1 = sum(tetramesa_sp1_abun, na.rm = TRUE),
    sp_tet2 = sum(tetramesa_sp2_abun, na.rm = TRUE),
    sp_bru1 = sum(bruchophagus_sp1_abun, na.rm = TRUE),
    sp_sht1 = sum(shotholeborer_sp1_abun, na.rm = TRUE),
    sp_chl1 = sum(chloropid_sp1_abun, na.rm = TRUE),
    sp_eur4 = sum(eurytomid_sp4_abun, na.rm = TRUE)
  ) %>%
  dplyr::select(site_code, sampling_date, season, 
                survey_type, plant_species, 
                disturbance_type:disturb_allow_overwinter, sp_tet1:sp_eur4) %>%
  slice(1)
com_data
View(com_data)

##########################################################################################
#####                     Visualise community composition - nMDS                    ######
##########################################################################################

# Keep only the LTS sites because of the differences in sampling effort 
com_data_lts <- com_data %>%
  dplyr::filter(survey_type == "LTS")

# Now we need to extract only the columns with species abundance data
com_data_sp <- com_data_lts %>%
  group_by(site_code, sampling_date, season, disturb_allow_overwinter) %>%
  dplyr::select(starts_with("sp_")) %>%
  ungroup()
com_data_sp

# Store this data.frame so that we can extract factor names from it later with ease 
org.data <- com_data_sp
org.data

# Now keep only the species abundances columns 
com_data_sp <- com_data_sp %>% 
  dplyr::select(-c(site_code, sampling_date, season, disturb_allow_overwinter))
com_data_sp

# Because of issues where rows sum to 0, we must use 
# zero-adjusted Bray-Curtis vals - add dummy species where abundance = 1
com_data_sp$sp_dum1 <- 1
com_data_sp

# Make this into a matrix
com_matrix <- as.matrix(com_data_sp)

# Perform a preliminary nMDS
comm.mds <- metaMDS(comm = com_matrix, # matrix of species abundance data
                    distance = "bray", # dissimilarity metric
                    trace = FALSE, # Suppress R information during processing of data
                    autotransform = FALSE,
                    k = 2, 
                    trymax = 100) # Stops VEGAN package from auto-transforming data

# R will not allow any analyses where row sums for the species x matrix are < 1, if this is the case, it throws:
# Error in cmdscale(dist, k = k) : NA values not allowed in 'd'
# In addition: Warning messages:
# 1: In distfun(comm, method = distance, ...) : you have empty rows: their dissimilarities may be meaningless in method â€œbrayâ€?
# 2: In distfun(comm, method = distance, ...) : missing values in results

# Plot a stressplot (Sheppard plot)
stressplot(comm.mds)

### Large scatter around the line suggests that original dissimilarities are not 
### well preserved in the reduced number of dimensions. 
### In this case, the sheppard plot looks okay (Goodness of fit (R2) > 0.9 = excellent) 

# Plot preliminary nMDS
plot(comm.mds)

### This is not very informative. Let's add some sample details to the plot. 

# Extract the XY co-ordinates from the nMDS plot
xy.mds <- data.frame(comm.mds$points)

# Now add the categorical factors from the original data.frame
# Now add the categorical factors from the original data.frame
xy.mds$disturb <- org.data$disturb_allow_overwinter
xy.mds$season <- org.data$season
xy.mds$site <- org.data$site_code

xy.mds <- xy.mds %>%
  mutate(season = as.factor(season),
         disturb = as.factor(disturb))

# Plot the nMDS with colour-coded points for disturbance regimes, and shapes for season
ggplot(data = xy.mds, aes(x = MDS1, 
                          y = MDS2,
                          shape = as.factor(season),
                          colour = as.factor(disturb))) + 
  geom_point(size = 2.5) +
  geom_jitter(width = 0.2, height = 0.2) +
  labs(x = "nMDS axis 1",
       y = "nMDS axis 2",
       colour = "Disturbance",
       shape = "Season") +
  theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(colour = "black"),
        axis.title.x = element_text(margin = unit(c(2, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
        legend.position = "right")

ggsave("./figures/fig_2_nMDS_plot.png", 
       width = 6, 
       height = 4)

# Extract the stress value 
### Stress provides a measure of the degree to which the distance between samples in 
### reduced dimensional space (usually 2-dimensions) corresponds with the actual 
### multivariate distance between the samples. 
### Lower stress values indicate greater conformity and therefore are desirable. 
### High stress values indicate that there was no 2-dimensional arrangement of your points 
### that reflect their similarities. 
### Clarke 1993 suggests the following guidelines for acceptable stress values: 
### <0.05 = excellent
### <0.10 = good
### <0.20 = usable 
### >0.20 = not acceptable
### Clarke, K. R. (1993). Non-parametric multivariate analysis of changes in community structure. Austral J Ecol 18: 117-143.
comm.mds$stress

# Perform PERMANOVA using the 'adonis' function (species abundances as dependent variable)
# Dependent varaible = species abundance matrix; interaction between season and disturbance
adonis.int <- adonis(com_data_sp ~ season * disturb, 
                    # Constrain permutations to within site (can change p vals)
          
                              strata = xy.mds$site_code, 
                    data = xy.mds,
                    permutations = 999)
adonis.int

# Save the output of PERMANOVA as a word file
capture.output(adonis.int,
               file="./results/PERManova_interaction.doc")

# Now without the interaction term
adonis.add <- adonis(com_data_sp ~ season + disturb, 
                     # Constrain permutations to within site (can change p vals)
                     strata = xy.mds$site_code, 
                     data = xy.mds,
                     permutations = 999)
adonis.add

# Save the output of PERMANOVA as a word file
capture.output(adonis.add,
               file="./results/PERManova_no_interaction.doc")

# Perform SIMPER analysis
# SIMPER assesses the contribution of each species the overall 
# Bray-Curtis dissimilarity index (i.e. zero-inflated B-C index in this case), 
sim.data <- com_data_sp
sim <- simper(sim.data, org.data$disturb_allow_overwinter,
              permutations = 999)
summary(sim)

# Save the output of SIMPER analysis by disturbance as a word file
capture.output(summary(sim), 
               file="./results/lts_simper_disturbance.doc")

sim <- simper(sim.data, org.data$season,
              permutations = 999)
summary(sim)

# Save the output of SIMPER analysis by seasonas a word file
capture.output(summary(sim), 
               file="./results/lts_simper_season.doc")

#################
# - Model-based multivariate analysis
################

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

# Run preliminary model
mod1 <- manyglm(insect.spp ~ com_data_lts$disturb_allow_overwinter *
                             com_data_lts$season, 
                family="poisson")

# Check model assumptions
plot(mod1)

### We see a fan-shape in the model residuals. 
# We need to model the variance differently. 

# Run negative-binomial model
mod2 <- manyglm(insect.spp ~ com_data_lts$disturb_allow_overwinter *
                             com_data_lts$season, 
                family="negative_binomial")

# Check the variance is better
plot(mod2)

### Aweh. Much better. 

### Now we are going to test the multivariate hypothesis of whether 
### insect species assemblges varied according to disturbance regimes. 
### This analysis uses likelihood-ratio tests (LRT) and 
### resampled P-values to test significance.

anova(mod2)
anova(mod2, 
      nBoot=999, 
      test="wald")

### We can see from this table that there is a significant effect of disturbance
### (LRT = 67.56.82, P = 0.001), meaning that the species composition 
### of herbivores clearly differs between disturbance regimes

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

# Save the output of 'mvabund' GLZ analysis as a word file
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

# Model 2 of the Youtube video - latent variable with no covariates
ins_nb <- boral(sim.data, 
                family = "negative.binomial",
                lv.control = list(num.lv = 2),
                row.eff = "fixed",
                save.model = TRUE,
                mcmc.control = example_mcmc_control)
summary(ins_nb)

# Get the Geweke convergence statistic (z scores)
ins_nb$geweke.diag
gew.pvals <- 2*pnorm(abs(unlist(ins_nb$geweke.diag[[1]])), 
                     lower.tail = FALSE)
p.adjust(gew.pvals, method = "holm")

# Check model convergnece by MCMC using 'coda' package
plot(get.mcmcsamples(ins_nb))

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





# Extract coefs
coef.bor <- ins_nb$lv.coefs.mean
coef.bor
coef.bor <- ins_nb$lv.coefs.mean
coef.bor

attributes(ins_nb)
coef.bor <- ins_nb$hpdintervals
coef.bor

############################################################
# - Calculate summary statistics for Tetramesa sp. 1 
############################################################

# Calculate average abundance per site 
abun_tet1 <- com_data_lts %>%
  group_by(season, disturb_allow_overwinter) %>%
  summarise(mean = mean(sp_tet1),
            sd = sd(sp_tet1)) %>%
  mutate(species = "Tetramesa sp. 1")
abun_tet1

# Calculate prop of sites occupied 
occ_tet1 <- com_data_lts %>%
  group_by(site_code) %>%
  mutate(total_abun = cumsum(sp_tet1)) %>%
  mutate(occ = case_when(
    total_abun > 0 ~ 1,
    TRUE ~ 0
  )) %>%
  ungroup() %>%
  group_by(season, disturb_allow_overwinter) %>%
  summarise(occ = mean(occ)) %>%
  mutate(species = "Tetramesa sp. 1")
occ_tet1

# Calculate prop of visits it was recorded  
sites_tet1 <- com_data_lts %>%
  group_by(site_code, sampling_date) %>%
  mutate(total_abun = cumsum(sp_tet1)) %>%
  mutate(occ = case_when(
    total_abun > 0 ~ 1,
    TRUE ~ 0
  )) %>%
  ungroup() %>%
  group_by(season, disturb_allow_overwinter) %>%
  summarise(mean = mean(occ),
            sd = sd(occ)) %>%
  mutate(species = "Tetramesa sp. 1") 
sites_tet1

#sites <- sites_tet1 %>% 
#  ungroup() %>%
#  filter(occ == 1) %>% 
#  filter(disturb_allow_overwinter == "Disturbed") %>%
#  mutate(species = "Tetramesa sp. 1") 
#sites

############################################################
# - Calculate summary statistics for Tetramesa sp. 2 
############################################################

# Calculate average abundance per site 
abun_tet2 <- com_data_lts %>%
  group_by(season, disturb_allow_overwinter) %>%
  summarise(mean = mean(sp_tet2),
            sd = sd(sp_tet2)) %>%
  mutate(species = "Tetramesa sp. 2")
abun_tet2

# Calculate prop of sites occupied 
occ_tet2 <- com_data_lts %>%
  group_by(site_code) %>%
  mutate(total_abun = cumsum(sp_tet2)) %>%
  mutate(occ = case_when(
    total_abun > 0 ~ 1,
    TRUE ~ 0
  )) %>%
  ungroup() %>%
  group_by(season, disturb_allow_overwinter) %>%
  summarise(occ = mean(occ)) %>%
  mutate(species = "Tetramesa sp. 2")
occ_tet2

# Calculate prop of visits it was recorded  
sites_tet2 <- com_data_lts %>%
  group_by(site_code, sampling_date) %>%
  mutate(total_abun = cumsum(sp_tet2)) %>%
  mutate(occ = case_when(
    total_abun > 0 ~ 1,
    TRUE ~ 0
  )) %>%
  ungroup() %>%
  group_by(season, disturb_allow_overwinter) %>%
  summarise(mean = mean(occ),
            sd = sd(occ)) %>%
  mutate(species = "Tetramesa sp. 2") 
sites_tet2

############################################################
# - Calculate summary statistics for Bruchophagus sp. 1 
############################################################

# Calculate average abundance per site 
abun_bru1 <- com_data_lts %>%
  group_by(season, disturb_allow_overwinter) %>%
  summarise(mean = mean(sp_bru1),
            sd = sd(sp_bru1)) %>%
  mutate(species = "Bruchophagus sp. 1")
abun_bru1

# Calculate prop of sites occupied 
occ_bru1 <- com_data_lts %>%
  group_by(site_code) %>%
  mutate(total_abun = cumsum(sp_bru1)) %>%
  mutate(occ = case_when(
    total_abun > 0 ~ 1,
    TRUE ~ 0
  )) %>%
  ungroup() %>%
  group_by(season, disturb_allow_overwinter) %>%
  summarise(occ = mean(occ)) %>%
  mutate(species = "Bruchophagus sp. 1")
occ_bru1

# Calculate prop of visits it was recorded  
sites_bru1 <- com_data_lts %>%
  group_by(site_code, sampling_date) %>%
  mutate(total_abun = cumsum(sp_bru1)) %>%
  mutate(occ = case_when(
    total_abun > 0 ~ 1,
    TRUE ~ 0
  )) %>%
  ungroup() %>%
  group_by(season, disturb_allow_overwinter) %>%
  summarise(mean = mean(occ),
            sd = sd(occ)) %>%
  mutate(species = "Bruchophagus sp. 1") 
sites_bru1

############################################################
# - Calculate summary statistics for Eurytomidae sp. 4 
############################################################

# Calculate average abundance per site 
abun_eur4 <- com_data_lts %>%
  group_by(season, disturb_allow_overwinter) %>%
  summarise(mean = mean(sp_eur4),
            sd = sd(sp_eur4)) %>%
  mutate(species = "Eurytomidae sp. 4")
abun_eur4

# Calculate prop of sites occupied 
occ_bru1 <- com_data_lts %>%
  group_by(site_code) %>%
  mutate(total_abun = cumsum(sp_bru1)) %>%
  mutate(occ = case_when(
    total_abun > 0 ~ 1,
    TRUE ~ 0
  )) %>%
  ungroup() %>%
  group_by(season, disturb_allow_overwinter) %>%
  summarise(occ = mean(occ)) %>%
  mutate(species = "Bruchophagus sp. 1")
occ_bru1

# Calculate prop of visits it was recorded  
sites_eur4 <- com_data_lts %>%
  group_by(site_code, sampling_date) %>%
  mutate(total_abun = cumsum(sp_eur4)) %>%
  mutate(occ = case_when(
    total_abun > 0 ~ 1,
    TRUE ~ 0
  )) %>%
  ungroup() %>%
  group_by(season, disturb_allow_overwinter) %>%
  summarise(mean = mean(occ),
            sd = sd(occ)) %>%
  mutate(species = "Eurytomidae sp. 4") 
sites_eur4

############################################################
# - Calculate summary statistics for Scolytidae sp. 1
############################################################

# Calculate average abundance per site 
abun_sht1 <- com_data_lts %>%
  group_by(season, disturb_allow_overwinter) %>%
  summarise(mean = mean(sp_sht1),
            sd = sd(sp_sht1)) %>%
  mutate(species = "Scolytidae sp. 1")
abun_eur4

# Calculate prop of sites occupied 
occ_bru1 <- com_data_lts %>%
  group_by(site_code) %>%
  mutate(total_abun = cumsum(sp_bru1)) %>%
  mutate(occ = case_when(
    total_abun > 0 ~ 1,
    TRUE ~ 0
  )) %>%
  ungroup() %>%
  group_by(season, disturb_allow_overwinter) %>%
  summarise(occ = mean(occ)) %>%
  mutate(species = "Bruchophagus sp. 1")
occ_bru1

# Calculate prop of visits it was recorded  
sites_sht1 <- com_data_lts %>%
  group_by(site_code, sampling_date) %>%
  mutate(total_abun = cumsum(sp_sht1)) %>%
  mutate(occ = case_when(
    total_abun > 0 ~ 1,
    TRUE ~ 0
  )) %>%
  ungroup() %>%
  group_by(season, disturb_allow_overwinter) %>%
  summarise(mean = mean(occ),
            sd = sd(occ)) %>%
  mutate(species = "Scolytidae sp. 1") 
sites_sht1

############################################################
# - Calculate summary statistics for Chloropidae sp. 1
###########################################################

# Calculate average abundance per site 
abun_chl1 <- com_data_lts %>%
  group_by(season, disturb_allow_overwinter) %>%
  summarise(mean = mean(sp_chl1),
            sd = sd(sp_chl1)) %>%
  mutate(species = "Chloropidae sp. 1")
abun_chl1

# Calculate prop of sites occupied 
occ_bru1 <- com_data_lts %>%
  group_by(site_code) %>%
  mutate(total_abun = cumsum(sp_bru1)) %>%
  mutate(occ = case_when(
    total_abun > 0 ~ 1,
    TRUE ~ 0
  )) %>%
  ungroup() %>%
  group_by(season, disturb_allow_overwinter) %>%
  summarise(occ = mean(occ)) %>%
  mutate(species = "Bruchophagus sp. 1")
occ_bru1

# Calculate prop of visits it was recorded  
sites_chl1 <- com_data_lts %>%
  group_by(site_code, sampling_date) %>%
  mutate(total_abun = cumsum(sp_chl1)) %>%
  mutate(occ = case_when(
    total_abun > 0 ~ 1,
    TRUE ~ 0
  )) %>%
  ungroup() %>%
  group_by(season, disturb_allow_overwinter) %>%
  summarise(mean = mean(occ),
            sd = sd(occ)) %>%
  mutate(species = "Chloropidae sp. 1") 
sites_chl1

# ---------------------
# Plot abundance graphs
# ---------------------

# First, create the dataset we require to plot 
abun_all <- bind_rows(abun_tet1, abun_tet2, abun_bru1,
                      abun_eur4, abun_sht1, abun_chl1)
abun_all <- abun_all %>%
  mutate(species = as.factor(species)) %>%
  mutate(lower = mean - sd,
         upper = mean + sd) %>%
  mutate(lower = case_when(
    lower < 0 ~ 0,
    TRUE ~ lower
  )) %>%
  mutate(disturb_allow_overwinter = fct_relevel(disturb_allow_overwinter,
                                                "Undisturbed",
                                                "Disturbed"),
         species = fct_relevel(species,
                               "Tetramesa sp. 1",
                               "Tetramesa sp. 2",
                               "Bruchophagus sp. 1",
                               "Eurytomidae sp. 4",
                               "Scolytidae sp. 1", 
                               "Chloropidae sp. 1"))
abun_all

# Plot the graph 
ggplot(data = abun_all,aes(x = disturb_allow_overwinter,
                           fill = season)) +
  geom_point(aes(y = mean,
                 colour = season),
             position=position_dodge(1)) + 
  geom_errorbar(aes(ymin = lower, 
                    ymax = upper,
                    colour = season),
                position=position_dodge(1),
                width = 0.7) +
  scale_colour_manual(values=c("black", "gray60")) + 
  labs(x = "Disturbance regime",
       y = "Abundance \n (mean individuals per site)",
       colour = "Season") + 
  geom_vline(xintercept = 1.5, linetype = "dashed",
            colour = "gray60") +
  facet_wrap(~ species, nrow = 2, ncol = 3) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(colour = "black"),
        axis.title.x = element_text(margin = unit(c(2, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
        legend.position = "right") +
  guides(fill = "none")

ggsave("fig_abundance_by_season_disturbance.png", 
       width = 7, 
       height = 4,
       dpi = "print")
ggsave("fig_abundance_by_season_disturbance.pdf", 
       width = 7, 
       height = 4)

# ---------------------
# Plot proportion of surveys recorded graphs
# ---------------------

# First, create the dataset we require to plot 
sites_all <- bind_rows(sites_tet1, sites_tet2, sites_bru1,
                       sites_eur4, sites_sht1, sites_chl1)
sites_all <- sites_all %>%
  mutate(species = as.factor(species)) %>%
  mutate(lower = mean - sd,
         upper = mean + sd) %>%
  mutate(lower = case_when(
    lower < 0 ~ 0,
    TRUE ~ lower
  )) %>%
  ungroup() %>%
  mutate(disturb_allow_overwinter = fct_relevel(disturb_allow_overwinter,
                                                "Undisturbed",
                                                "Disturbed"),
         species = fct_relevel(species,
                               "Tetramesa sp. 1",
                               "Tetramesa sp. 2",
                               "Bruchophagus sp. 1",
                               "Eurytomidae sp. 4",
                               "Scolytidae sp. 1", 
                               "Chloropidae sp. 1"))
sites_all

# Plot the graph 
ggplot(data = sites_all,aes(x = disturb_allow_overwinter, 
                            fill = season)) +
  geom_errorbar(aes(ymin = mean, 
                    ymax = mean,
                    colour = season),
                position=position_dodge(1),
                width = 0.7) +
  labs(x = "Disturbance regime",
       y = "Proportion surveys recorded",
       colour = "Season") + 
  scale_colour_manual(values=c("black", "gray60")) + 
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0), 
                     limits = c(0, 1)) +
  geom_vline(xintercept = 1.5, linetype = "dashed",
             colour = "gray60") +
  facet_wrap(~ species, nrow = 2, ncol = 3) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(colour = "black"),
        axis.title.x = element_text(margin = unit(c(2, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
        legend.position = "right")

ggsave("fig_proportion_by_season_disturbance.png", 
       width = 7, 
       height = 4,
       dpi = "print")
ggsave("fig_proportion_by_season_disturbance.pdf", 
       width = 7, 
       height = 4)

########################################################
# - Does plant architechture change with disturbance? 
########################################################

# Get the right data 
head(raw_data)

# Get in correct format 
raw_data_diam <- raw_data %>%
  mutate(disturb_regime = as.factor(disturb_allow_overwinter)) %>%
  drop_na(tiller_diam) %>%
  drop_na(disturb_regime) %>%
  mutate(tiller_diam = tiller_diam / 10) %>%
  filter(tiller_diam >= 1)
raw_data_diam
View(raw_data_diam)

# Plot histogram
ggplot(data = raw_data_diam, aes(x = tiller_diam,
                            group = disturb_allow_overwinter)) +
  #geom_histogram(aes(colour = disturb_allow_overwinter,
  #                   fill = disturb_allow_overwinter)) +
  geom_density() +
  scale_x_continuous(breaks = c(0, 1, 2, 3, 4, 5, 6, 7), 
                     limits = c(0, 7)) +
  labs(x = "Tiller diameter (mm)",
       y = "Proportion of tillers") + 
  geom_vline(xintercept = 3.0, linetype = "dashed") + 
  annotate("text", 
           x = 6, 
           y = 0.4, 
           label = "Overlap = 42.86%") +
  annotate("text", 
           x = 6.04, 
           y = 0.3, 
           label = "KS test: P < 0.001") +
  facet_wrap(~ disturb_allow_overwinter, nrow = 2) +
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(colour = "black"),
        axis.title.x = element_text(margin = unit(c(2, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
        legend.position = "none")

ggsave("fig_proportion_tillers_density.png", 
       width = 5, 
       height = 4,
       dpi = "print")
ggsave("fig_proportion_tillers_density.pdf", 
       width = 5, 
       height = 4)

# Use the overlap package
library(overlapping)
g1 <- raw_data_diam %>%
  filter(disturb_regime == "Undisturbed") %>%
  filter(tiller_diam >= 1) %>%
  pull(tiller_diam) 
g1 = as.integer((g1))
g2 <- raw_data_diam %>%
  filter(disturb_regime == "Disturbed") %>%
  filter(tiller_diam >= 1) %>%
  pull(tiller_diam)
g2 = as.integer((g2))
dataList <- list( G1 = g1, G2 = g2)
overlap( dataList )$OV * 100

# KS test
ks.test(g1, g2)

