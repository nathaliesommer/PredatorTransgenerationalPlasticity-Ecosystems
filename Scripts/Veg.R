# Vegetation Allometry

library(dplyr)
library(tidyr)
library(broom)
library(purrr)
library(ggplot2)

## Data Import and Cleaning ----

### Biomass ----
vegmass <- read.csv("Data/Biomass-plots_predators.csv") # biomass plots 

vegmass <- vegmass %>% # unit conversion
  mutate(
    Biomass_Bag_g = Biomass_Bag_kg * 1000,
    Bag_g = Bag_kg * 1000,
    Biomass_g = Biomass_kg * 1000
  )

vegmass <- vegmass %>%
  select(-Biomass_Bag_kg, -Bag_kg, -Biomass_kg) 


### Cover ----
diversity <- read.csv("Data/Vegetation.cover.2023_predators.csv") # plant percent cover 2023


cage_diversity <- diversity %>%
  select(-BARE) %>% # remove bare ground
  rowwise() %>%
  mutate(Cover_Sum = sum(c_across(-c(Sample_ID)))) %>% # sum plant cover 
  ungroup()

cage_diversity <- cage_diversity %>% # normalize cover to 100%
  mutate(across(
    .cols = -c(Sample_ID),
    .fns = ~ round((.x / Cover_Sum) * 100)
  ))

cage_diversity_long <- cage_diversity %>%
  select(-Cover_Sum) %>% 
  pivot_longer(
    cols = -c(Sample_ID),
    names_to = "Species_ID",
    values_to = "Cover"
  )




functional_groups <- cage_diversity_long %>%
  group_by(Sample_ID) %>% 
  summarise(
    SORU = sum(Cover[Species_ID %in% c("SORU2", "SOCA6")], na.rm = TRUE),
    POPRC = sum(Cover[Species_ID == "POPRC"], na.rm = TRUE),
    MISC = sum(Cover[!(Species_ID %in% c("SORU2", "SOCA6", "POPRC"))], na.rm = TRUE)
  ) %>%
  pivot_longer(cols = SORU:MISC, names_to = "Species_ID", values_to = "Cover") %>%
  ungroup() %>% 
  drop_na()



## Allometry ----
# Use biomass plots to estimate the biomass in cages from % cover

# Calculate the biomass for each plant group within biomass plots
vegmass <- vegmass %>%
  mutate(
    SORU_Biomass = (SORU / 100) * Biomass_g,
    POPRC_Biomass = (POPRC / 100) * Biomass_g,
    MISC_Biomass = (MISC / 100) * Biomass_g
  )

# Fit models for each group
SORU_model <- lm(SORU_Biomass ~ SORU, data = vegmass)
POPRC_model <- lm(POPRC_Biomass ~ POPRC, data = vegmass)
MISC_model <- lm(MISC_Biomass ~ MISC, data = vegmass)

# Diagnostic plots for SORU
diagnostic_plot_SORU <- ggplot(vegmass, aes(x = SORU, y = SORU_Biomass)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("SORU Check")

# Diagnostic plots for POPRC
diagnostic_plot_POPRC <- ggplot(vegmass, aes(x = POPRC, y = POPRC_Biomass)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("POPRC Check")

# Diagnostic plots for MISC
diagnostic_plot_MISC <- ggplot(vegmass, aes(x = MISC, y = MISC_Biomass)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle("MISC Check")



functional_groups_wide <- functional_groups %>%
  pivot_wider(names_from = Species_ID, values_from = Cover) %>%
  # Replace missing cover values (e.g., if some species are missing from some plots) with 0
  mutate(
    SORU = as.numeric(coalesce(SORU, 0)),
    POPRC = as.numeric(coalesce(POPRC, 0)),
    MISC = as.numeric(coalesce(MISC, 0))
  )


model_fits <- list(
  SORU_model = lm(SORU_Biomass ~ SORU, data = vegmass),
  POPRC_model = lm(POPRC_Biomass ~ POPRC, data = vegmass),
  MISC_model = lm(MISC_Biomass ~ MISC, data = vegmass)
)

# Apply models row by row
functional_groups_wide <- functional_groups_wide %>%
  rowwise() %>%
  mutate(
    SORU_Biomass = predict(model_fits$SORU_model, newdata = data.frame(SORU = SORU)),
    POPRC_Biomass = predict(model_fits$POPRC_model, newdata = data.frame(POPRC = POPRC)),
    MISC_Biomass = predict(model_fits$MISC_model, newdata = data.frame(MISC = MISC))
  ) %>%
  ungroup()

### Diversity Index Calculation ----

# Function to calculate Shannon-Weiner diversity index
calculate_shannon_weiner <- function(abundances) {
  proportions <- abundances / sum(abundances)
  proportions <- proportions[proportions > 0]  # Remove zero proportions
  -sum(proportions * log(proportions))
}

# Function to calculate diversity indices
calculate_diversity_indices <- function(data) {
  data %>%
    rowwise() %>%
    mutate(
      PlantRichness = sum(c_across(-c(Sample_ID)) > 0),
      PlantDiversity = calculate_shannon_weiner(c_across(-c(Sample_ID)))
    ) %>%
    ungroup() %>%
    select(Sample_ID, PlantRichness, PlantDiversity)
}



# Calculate diversity indices for each dataset
diversity_calc <- calculate_diversity_indices(diversity)
