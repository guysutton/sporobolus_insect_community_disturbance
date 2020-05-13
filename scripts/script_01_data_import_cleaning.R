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

######
### - Load required packages                             
######

library(vegan)
library(MASS)
library(mvabund)
library(tidyr)
library(readr)
library(tidyverse)

######
### - Import data   
######

# Import raw data 
raw_data <- read_csv2("./data_raw/all_surveys_data.csv")

# Check the data 
head(raw_data)

# Drop some columns that we don't need 
raw_data <- raw_data %>%
  dplyr::select(site:insect_abun_total, disturbance_type:disturb_allow_overwinter) %>%
  # Keep plant records for only two target plants
  dplyr::filter(plant_species == c("Sporobolus pyramidalis", "Sporobolus natalensis")) %>%
  # Clean the disturbance variable
  mutate(disturb_allow_overwinter = case_when(
    disturb_allow_overwinter == "No" ~ "Disturbed",
    disturb_allow_overwinter == "Yes" ~ "Undisturbed"))
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
    sp_eur4 = sum(eurytomid_sp4_abun, na.rm = TRUE)) %>%
  dplyr::select(site_code, sampling_date, season, 
                survey_type, plant_species, 
                disturbance_type:disturb_allow_overwinter, sp_tet1:sp_eur4) %>%
  slice(1)
com_data

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

# Calculate species richness per site visit
rich_data <- com_data_lts %>%
  dplyr::select(
    site_code, sampling_date, season, disturb_allow_overwinter,
                sp_tet1:sp_eur4)
rich_data

rich_data$sp_rich <- rowSums(rich_data[, 5:10]!=0)
rich_data <- rich_data %>%
  mutate(disturbance = as.factor(disturb_allow_overwinter),
         season = as.factor(season))

# Plot data 
ggplot(data = rich_data, aes(x = disturb_allow_overwinter,
                             y = sp_rich,
                             colour = season)) +
  geom_jitter(height = 0.1, width = 0.2) +
  scale_colour_manual(values=c("black", "gray60")) + 
  labs(x = "Disturbance regime",
       y = "Species richness",
       colour = "Season") + 
  theme_classic() +
  theme(panel.border = element_rect(colour = "black", fill = NA),
        axis.text = element_text(colour = "black"),
        axis.title.x = element_text(margin = unit(c(2, 0, 0, 0), "mm")),
        axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
        legend.position = "right") +
  guides(fill = "none")