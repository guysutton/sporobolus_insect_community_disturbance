##########################################################################################
##########################################################################################
##########################################################################################
#####              Does tiller size influence insect attack rates?                   #####
#####              Script author: Guy F. Sutton                                      #####
#####              Institution  : Department of Zoology and Entomology,              #####  
#####                             Rhodes University,                                 #####
#####                             South Africa                                       #####
#####              Last update  : 24/10/2018                                         #####
##########################################################################################
##########################################################################################
##########################################################################################

##########################################################################################
##########################################################################################
##########################################################################################
#####                          1. Load required packages                             #####
##########################################################################################
##########################################################################################
##########################################################################################

library(broom)
library(MuMIn)
library(tidyverse) 
library(lmtest)
library(readxl)
library(pbkrtest)
library(multcomp)
library(visreg)
library(boot)
library(pbnm)
library(AICcmodavg)
library(dotwhisker)
library(lme4)
library(lmerTest)

##########################################################################################
##########################################################################################
##########################################################################################
#####                          2. Import and clean data                              #####
##########################################################################################
##########################################################################################
##########################################################################################

### Import all survey data
lts.data <- readxl::read_xlsx("./data/all_surveys_data.xlsx")
str(lts.data)

# Rename columns that imported incorrectly
colnames(lts.data)[1] <- "site"
str(lts.data)

# Convert insect pres/absence columns from interger to numeric
field_data <- lts.data %>%
  mutate(tet1_pres = as.factor(tetramesa_sp1_pres),
         tet2_pres = as.factor(tetramesa_sp2_pres),
         bru1_pres = as.factor(bruchophagus_sp1_pres),
         eur1_pres = as.factor(eurytomid_sp4_pres),
         sht1_pres = as.factor(shotholeborer_sp1_pres),
         chl1_pres = as.factor(chloropid_sp1_pres),
         tiller_diam = as.numeric(tiller_diam),
         tiller_height = as.numeric(tiller_height))
str(field_data)

# Keep only the columns we need
field_data <- field_data %>%
  filter(plant_species == c("Sporobolus pyramidalis", "Sporobolus natalensis")) %>%
  dplyr::select(site, sampling_date, plant_species:tiller_height, tet1_pres:chl1_pres)
str(field_data)

##########################################################################################
##########################################################################################
##########################################################################################
#####                           3. Visualise data                                    #####
##########################################################################################
##########################################################################################
##########################################################################################

# Select tiller_diam attacked by Tetramesa_sp1
tet1.diam <- field_data %>%
  filter(tet1_pres == 1)
tet1.diam$species <- "Tetramesa sp. 1"
tet2.diam <- field_data %>%
  filter(tet2_pres == 1)
tet2.diam$species <- "Tetramesa sp. 2"
bru1.diam <- field_data %>%
  filter(bru1_pres == 1)
bru1.diam$species <- "Bruchophagus sp. 1"
eur1.diam <- field_data %>%
  filter(eur1_pres == 1)
eur1.diam$species <- "Eurytomidae sp. 4"
sht1.diam <- field_data %>%
  filter(sht1_pres == 1)
sht1.diam$species <- "Scolytidae sp. 1"
chl1.diam <- field_data %>%
  filter(chl1_pres == 1)
chl1.diam$species <- "Chloropidae sp. 1"

# Combine data
species.diam <- bind_rows(tet1.diam, tet2.diam, bru1.diam,
                          eur1.diam, sht1.diam, chl1.diam)

# Change order that factors will be plotted on the graph 
species.diam <- species.diam %>% 
  mutate(species = fct_relevel(species, 
                               "Tetramesa sp. 1", 
                               "Tetramesa sp. 2", 
                               "Bruchophagus sp. 1",
                               "Eurytomidae sp. 4",
                               "Scolytidae sp. 1",
                               "Chloropidae sp. 1"))

# Create df with background data (all data without 'species' column)
bg.data <- field_data

# Plot graphs 
ggplot(species.diam, aes(tiller_diam, fill = species)) + 
  geom_histogram(data = bg.data, fill = "grey", alpha = 0.5) + 
  geom_histogram() + 
  scale_y_continuous(labels=function(x)x/length(x)) + 
  stat_bin(bins = 30) + 
  facet_wrap(~species, nrow = 2) + 
  theme_bw() + 
  labs(y = "Percentage of tillers attacked (%)",
       x = "Tiller diameter (mm)") + 
  theme(
    axis.ticks = element_blank(),
    axis.text.y = element_text(color = "black", size = 12),
    axis.text.x = element_text(color = "black", size = 12, hjust = 1, vjust = 1),   
    legend.key.size = unit(1.2, 'cm'),
    panel.grid = element_blank(),
    axis.title = element_text(color = 'black'),
    axis.title.y = element_text(margin = margin(t = 0, r = 18, b = 0, l = 0), size = 13),
    axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0), size = 13),
    legend.position = "none",
    legend.box.background = element_rect(),
    legend.box.margin = margin(4, 4, 4, 4))

# Save the plot 
ggsave("tiller_selection_plant_diameter.pdf", 
       width = 6, 
       height = 4, 
       dpi = "print")

ggsave("tiller_selection_plant_diameter.png", 
       width = 6, 
       height = 4)

##########################################################################################
##########################################################################################
##########################################################################################
#####                           4. Model the data                                    #####
##########################################################################################
##########################################################################################
##########################################################################################

##########################################################################################
### 4.1. Determine the best random effects structure for GLMM's 
##########################################################################################

# Select the random effects structure that minimises AICc and maximises -LogLik

# Hold all fixed effects constant and plot varying combinations of random effects,
# including a null model (i.e. do we need to include random effects?)
surv1 <- glmer(tet1_pres ~ 
                + (1 | site), 
               family = binomial(link = "logit"), 
               data = field_data, 
               verbose = FALSE,
               nAGQ = 0)
surv2 <- glmer(tet1_pres ~ 
                 + (1| sampling_date), 
               family = binomial(link = "logit"), 
               data = field_data, 
               verbose = FALSE,
               nAGQ = 0)
surv3 <- glmer(tet1_pres ~ 
                 + (1 | site) + (1| sampling_date), 
               family = binomial(link = "logit"), 
               data = field_data, 
               verbose = FALSE,
               nAGQ = 0)
surv4 <- glmer(tet1_pres ~ 
                 + (1| site/sampling_date), 
               family = binomial(link = "logit"), 
               data = field_data, 
               verbose = FALSE,
               nAGQ = 0)
surv5 <- glm(tet1_pres ~ 1,  
             family = binomial(link = "logit"), 
             data = field_data)

# Print LogLik values for each candidate model
logLik(surv1)
logLik(surv2)
logLik(surv3)
logLik(surv4)
logLik(surv5)

# Print AICc values for each candidate model 
MuMIn::AICc(surv1, surv2, surv3, surv4, surv5)

### Best random effects structure allows for random site and sampling_date intercepts

##########################################################################################
### 4.2. Fit all candidate models and select best performing models 
##########################################################################################

# Fit all candidate models
surv1 <- glmer(tet1_pres ~ tiller_diam + tiller_height + plant_species
               + (1 | site) + (1| sampling_date), 
               family = binomial(link = "logit"), 
               data = field_data, 
               verbose = FALSE,
               nAGQ = 0)
surv2 <- glmer(tet1_pres ~ tiller_diam + tiller_height 
               + (1 | site) + (1| sampling_date), 
               family = binomial(link = "logit"), 
               data = field_data, 
               verbose = FALSE,
               nAGQ = 0)
surv3 <- glmer(tet1_pres ~ tiller_height + plant_species
               + (1 | site) + (1| sampling_date), 
               family = binomial(link = "logit"), 
               data = field_data, 
               verbose = FALSE,
               nAGQ = 0)
surv4 <- glmer(tet1_pres ~ tiller_diam 
               + (1 | site) + (1| sampling_date), 
               family = binomial(link = "logit"), 
               data = field_data, 
               verbose = FALSE,
               nAGQ = 0)
surv5 <- glmer(tet1_pres ~ tiller_height
               + (1 | site) + (1| sampling_date), 
               family = binomial(link = "logit"), 
               data = field_data, 
               verbose = FALSE,
               nAGQ = 0)
surv6 <- glmer(tet1_pres ~ plant_species
               + (1 | site) + (1| sampling_date), 
               family = binomial(link = "logit"), 
               data = field_data, 
               verbose = FALSE,
               nAGQ = 0)
surv7 <- glmer(tet1_pres ~ 
               + (1 | site) + (1| sampling_date), 
               family = binomial(link = "logit"), 
               data = field_data, 
               verbose = FALSE,
               nAGQ = 0)
surv8 <- glm(tet1_pres ~ 1, 
               family = binomial(link = "logit"), 
               data = field_data)

### Get model output table 
# Send all survival models to a list of models (glmer models)
cand.mods <- list(surv1, surv2, surv3, surv4, surv5, surv6, surv7)

# Create a vector of model names to track models in the next tables 
mod.names <- paste("surv_model", 1:length(cand.mods), sep = " ")

## Get AIC table with LogLik 
glmer.mods <- aictab(cand.set = cand.mods, 
                     modnames = mod.names, 
                     sort = TRUE, 
                     digits = 2, 
                     LL = TRUE)

# Send null model to a list of models (glm null models)
cand.mods <- list(surv8)

# Create a vector of model names to track models in the next tables 
mod.names <- paste("null_model", 1:length(cand.mods), sep = " ")

## Get AIC table with LogLik 
null.mod <- aictab(cand.set = cand.mods, 
                   modnames = mod.names, 
                   sort = TRUE, 
                   digits = 2,  
                   LL = TRUE)

# Combine glmer and glm model output tables
model.select <- rbind(glmer.mods, null.mod)
print(surv.mod <- rbind(glmer.mods, null.mod))

# Write table to csv
write.table(model.select, 
            file="surv_model_selection_aic_table.csv", 
            row.names = FALSE, 
            sep=";")

##########################################################################################
### 3.3. Get parameter estimates from top-model  
##########################################################################################

# Print summary 
summary(surv4)

##########################################################################################
### 4.4. Calculate R-squared model fit statistics for candidate model set
##########################################################################################

# model = surv4
mod_ref <- surv4
mod_0 <- surv8
VarF <- var(as.vector(model.matrix(mod_ref) %*% fixef(mod_ref)))
nuN <- 1/attr(VarCorr(mod_0), "sc")^2 # note that glmer report 1/nu not nu as resiudal variance
VarOdN <- 1/nuN # the delta method
VarOlN <- log(1 + 1/nuN) # log-normal approximation
VarOtN <- trigamma(nuN) # trigamma function
c(VarOdN = VarOdN, VarOlN = VarOlN, VarOtN = VarOtN)
nuF <- 1/attr(VarCorr(mod_ref), "sc")^2 # note that glmer report 1/nu not nu as resiudal variance
VarOdF <- 1/nuF # the delta method
VarOlF <- log(1 + 1/nuF) # log-normal approximation
VarOtF <- trigamma(nuF) # trigamma function
c(VarOdF = VarOdF, VarOlF = VarOlF, VarOtF = VarOtF)
marg.3.m <- VarF/(VarF + sum(as.numeric(VarCorr(mod_ref))) + VarOtF)
marg.3.c <- (VarF + sum(as.numeric(VarCorr(mod_ref))))/(VarF +sum(as.numeric(VarCorr(mod_ref)))+VarOtF)
marg.3.m
marg.3.c

# Marginal R2 = Proportion of variance explained by fixed factors only
# Conditional R2 = Proportion of variance explained by fixed and random factors jointly

### Only 14.2% of the variation in Tet1 attack rates is explained by tiller diameter. 
### A further 2.6% is accounted for by tiller_height and plant_species, combined. 





