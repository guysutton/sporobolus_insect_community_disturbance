##########################################################################################
##########################################################################################
##########################################################################################
#####        How do season and disturbance regime affect insect communities?   ###########
#####        Script 02: Distance-based analyses (nMDS and PERMANOVA)           ###########
#####        Script author: Guy F. Sutton                                      ###########
#####        Institution  : Department of Zoology and Entomology,              ###########  
#####                       Rhodes University,                                 ###########
#####                       South Africa                                       ###########
#####        Last update  : 09/05/2020                                         ###########
##########################################################################################
##########################################################################################
##########################################################################################

######
### - (1) Visualise community composition - nMDS 
######

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
# Hence, why we used the zero-adjusted Bray-Curtis index. NB!

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
xy.mds$disturb <- org.data$disturb_allow_overwinter
xy.mds$season <- org.data$season
xy.mds$site <- org.data$site_code

# Changes the covariates into factors
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
  scale_colour_manual(values=c("black", "gray60")) + 
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
# - Stress provides a measure of the degree to which the distance between samples in 
#   reduced dimensional space (usually 2-dimensions) corresponds with the actual 
#   multivariate distance between the samples. 
#   - Lower stress values indicate greater conformity and therefore are desirable. 
#   - High stress values indicate that there was no 2-dimensional arrangement of 
#     your points that reflect their similarities. 
#   - Clarke 1993 suggests the following guidelines for acceptable stress values: 
#     - <0.05 = excellent
#     - <0.10 = good
#     - <0.20 = usable 
#     - >0.20 = not acceptable
#     - Clarke, K. R. (1993). Non-parametric multivariate analysis of changes in 
#       community structure. Austral J Ecol 18: 117-143.
comm.mds$stress

# Here, our stress value = 0.12, so our nMDS is okay. 

######
### - (2) PERMANOVA
######

# Perform PERMANOVA using the 'adonis' function 
# - Dependent variable = species abundances matrix
# - Predictors = Disturbance * Season (check whether interaction is required)

# Fit global PERManova (i.e. with interaction term)
adonis.int <- adonis(com_data_sp ~ season * disturb, 
                     # Constrain permutations to within site (can change p vals)
                     strata = xy.mds$site_code, 
                     data = xy.mds,
                     permutations = 999)
adonis.int

# Here, interaction term is not required. Refit model below without interaction. 

# Save the output of PERMANOVA as a word file (to look at later)
capture.output(adonis.int,
               file="./results/PERManova_interaction.doc")

# Fit PERManova without the interaction term
adonis.add <- adonis(com_data_sp ~ season + disturb, 
                     # Constrain permutations to within site (can change p vals)
                     strata = xy.mds$site_code, 
                     data = xy.mds,
                     permutations = 999)
adonis.add

# Save the output of PERMANOVA as a word file
capture.output(adonis.add,
               file="./results/PERManova_no_interaction.doc")

######
### - (3) SIMPER analysis 
######

# SIMPER assesses the contribution of each species the overall 
# - We need the Bray-Curtis dissimilarity index
sim.data <- com_data_sp

# Run SIMPER by disturbance regime
sim <- simper(sim.data, org.data$disturb_allow_overwinter,
              permutations = 999)
summary(sim)

# Save the output of SIMPER analysis by disturbance as a word file
capture.output(summary(sim), 
               file="./results/lts_simper_disturbance.doc")

# Run SIMPER by season
sim <- simper(sim.data, org.data$season,
              permutations = 999)
summary(sim)

# Save the output of SIMPER analysis by seasonas a word file
capture.output(summary(sim), 
               file="./results/lts_simper_season.doc")