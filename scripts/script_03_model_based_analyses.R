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

# Run negative-binomial model
mod2 <- manyglm(insect.spp ~ com_data_lts$disturb_allow_overwinter *
                  com_data_lts$season, 
                family="negative_binomial")

# Check the variance is better
plot(mod2)

### Much better. 

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