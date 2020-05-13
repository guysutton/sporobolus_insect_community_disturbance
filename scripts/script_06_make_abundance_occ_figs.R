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
ggplot(data = abun_all, aes(x = disturb_allow_overwinter,
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

ggsave("./figures/fig_abundance_by_season_disturbance.png", 
       width = 7, 
       height = 4,
       dpi = "print")
ggsave("./figures/fig_abundance_by_season_disturbance.pdf", 
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

ggsave("./figures/fig_proportion_by_season_disturbance.png", 
       width = 7, 
       height = 4,
       dpi = "print")
ggsave("./figures/fig_proportion_by_season_disturbance.pdf", 
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

ggsave("./figures/fig_proportion_tillers_density.png", 
       width = 5, 
       height = 4,
       dpi = "print")
ggsave("./figures/fig_proportion_tillers_density.pdf", 
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