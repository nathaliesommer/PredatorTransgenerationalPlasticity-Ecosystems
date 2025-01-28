# Vegetation Allometry
## Code by N.R. Sommer

# Object needed from this script: diversity_calc; functional_groups_wide

library(dplyr)
library(tidyr)
library(broom)
library(purrr)
library(ggplot2)

## Data Import and Cleaning ----

### Biomass ----
vegmass <- read.csv("Data/Biomass-plots_predators.csv") # biomass plots 

vegmass <- vegmass %>% 
  mutate(
    Biomass_Bag_g = case_when(
      Year == 2023 ~ Biomass_Bag_kg * 1000,
      TRUE ~ Biomass_Bag_kg 
    ),
    Bag_g = case_when(
      Year == 2023 ~ Bag_kg * 1000,
      TRUE ~ Bag_kg 
    ),
    Biomass_g = case_when(
      Year == 2023 ~ Biomass_kg * 1000,
      TRUE ~ Biomass_kg 
    )
  )

vegmass <- vegmass %>%
  select(-Biomass_Bag_kg, -Bag_kg, -Biomass_kg) 


### Cover ----
diversity <- read.csv("Data/Vegetation.cover_predators.csv") # plant percent cover


cage_diversity <- diversity %>%
  select(-BARE) %>% # remove bare ground
  rowwise() %>%
  mutate(Cover_Sum = sum(c_across(-c(Sample_ID, Year)))) %>% # sum plant cover 
  ungroup()

cage_diversity <- cage_diversity %>% # normalize cover to 100%
  mutate(across(
    .cols = -c(Sample_ID, Year),
    .fns = ~ round((.x / Cover_Sum) * 100)
  ))

cage_diversity_long <- cage_diversity %>%
  select(-Cover_Sum) %>% 
  pivot_longer(
    cols = -c(Sample_ID, Year),
    names_to = "Species_ID",
    values_to = "Cover"
  )




functional_groups <- cage_diversity_long %>%
  group_by(Sample_ID, Year) %>% 
  summarise(
    SORU = sum(Cover[Species_ID %in% c("SORU2", "SOCA6")], na.rm = TRUE),
    POPRC = sum(Cover[Species_ID == "POPRC"], na.rm = TRUE),
    MISC = sum(Cover[!(Species_ID %in% c("SORU2", "SOCA6", "POPRC"))], na.rm = TRUE)
  ) %>%
  pivot_longer(cols = SORU:MISC, names_to = "Species_ID", values_to = "Cover") %>%
  ungroup() %>% 
  drop_na()

# Create functional_groups_wide
functional_groups_wide <- functional_groups %>%
  pivot_wider(
    names_from = Species_ID,
    values_from = Cover
  )



## Allometry ----
# Use biomass plots to estimate the biomass in cages from % cover, separated by Site and Year

# Calculate the biomass for each plant group within biomass plots
vegmass <- vegmass %>%
  mutate(
    SORU_Biomass = (SORU / 100) * Biomass_g,
    POPRC_Biomass = (POPRC / 100) * Biomass_g,
    MISC_Biomass = (MISC / 100) * Biomass_g
  )

# Fit models for each group by Site and Year
model_fits <- vegmass %>%
  group_by(Site, Year) %>%
  group_split() %>%
  map(function(df) {
    list(
      SORU_model = lm(SORU_Biomass ~ SORU, data = df),
      POPRC_model = lm(POPRC_Biomass ~ POPRC, data = df),
      MISC_model = lm(MISC_Biomass ~ MISC, data = df)
    )
  })

# Then apply the models
functional_groups_wide <- functional_groups_wide %>%
  mutate(Site = "UP") %>%  # Add Site if not already present
  group_by(Site, Year) %>%
  group_modify(function(df, group_keys) {
    # Get the corresponding model for this Year
    model_index <- if(group_keys$Year == 2021) 1 else 2
    
    site_year_models <- model_fits[[model_index]]
    
    df %>%
      mutate(
        SORU_Biomass = predict(site_year_models$SORU_model, newdata = data.frame(SORU = SORU)),
        POPRC_Biomass = predict(site_year_models$POPRC_model, newdata = data.frame(POPRC = POPRC)),
        MISC_Biomass = predict(site_year_models$MISC_model, newdata = data.frame(MISC = MISC))
      )
  }) %>%
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
    select(Sample_ID, PlantRichness, PlantDiversity, Year)
}



# Calculate diversity indices for each dataset
diversity_calc <- calculate_diversity_indices(diversity)

