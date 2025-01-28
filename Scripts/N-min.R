### N-mineralization ----
## Code by N.R. Sommer

# Object needed from this script: N_min_calc


## Notes
  # Nitrate = (mg N-NO3/mL)
  # Ammonium = (mg N-NH4/mL)
  # Bulk density = g/cm^3
  # Nitrification was calculated as:
    # NO3 at 28 days minus initial NO3 at time zero; 
    # N mineralization equals NH4 + NO3 at 28 days minus NH4 + NO3 at time zero (Goodale and Aber 2001).
    # This was multiplied by bulk density and the appropriate conversions to obtain results on a g/cm^3 basis Pastor et al 1987: 


library(dplyr)

# Read-in data

N_min <- read.csv("Data/N-min.csv") # NOTE: 2023 N-min data was blank-adjusted by the processing lab
N_min_2021 <- read.csv("Data/2021_N-min.csv") # Baseline data not blank-adjusted, calc below
BulkDensity <- read.csv("Data/BulkDensity.csv")


# Standardize N_min_2021 with blank correction and column names
N_min_2021 <- N_min_2021 %>%
  group_by(Batch) %>%
  mutate(
    N_NO3_blank_avg = mean(`N.NO3..mg.per.mL.`[Sample.ID == "Blank"]),
    N_NH4_blank_avg = mean(`N.NH4..mg.per.mL.`[Sample.ID == "Blank"])
  ) %>%
  mutate(
    NO3 = `N.NO3..mg.per.mL.` - N_NO3_blank_avg,
    NH4 = `N.NH4..mg.per.mL.` - N_NH4_blank_avg
  ) %>%
  mutate(
    Site = "UP"
  ) %>% 
  ungroup() %>%
  select(Sample.ID, Day, NO3, NH4, Site) %>%  # Keep only needed columns
  rename(Sample_ID = Sample.ID)  # Standardize ID column name

# Standardize N_min column names
N_min <- N_min %>%
  rename(
    NO3 = blank.corrected.N.NO3..mg.per.mL.,
    NH4 = blank.corrected.N.NH4..mg.per.mL.
  )

# Get matching samples and add years
matching_samples <- unique(N_min$Sample_ID)

N_min_2021 <- N_min_2021 %>%
  filter(Sample_ID %in% matching_samples) %>%
  mutate(Year = 2021)

N_min <- N_min %>%
  mutate(Year = 2023)

# Combine datasets
N_min_combined <- bind_rows(N_min_2021, N_min)







# Bulk Density calculation
BulkDensity <- BulkDensity %>% 
  mutate(SoilMass = Weight_SoilTray-Weight_Tray) %>% 
  mutate(SoilMassAdj = SoilMass - (Materials_volume_mL_end - Materials_volume_mL_start)) %>% 
  mutate(SoilBulkDensity = (SoilMassAdj / 308.89)) # units = g/cm^3, where 308.89 is the cm^3 of the sampler


# Bulk Density, unit conversion 
BulkDensity_avg <- BulkDensity %>%
  group_by(Site) %>%
  summarize(SoilBulkDensity_avg = mean(SoilBulkDensity))


# Update calculate_rate function to use new column names
calculate_rate <- function(data, soil_depth = 10) {
  data %>%
    left_join(BulkDensity_avg, by = "Site") %>%
    group_by(Sample_ID, Year) %>%
    filter(Day == min(Day) | Day == max(Day)) %>%
    arrange(Day) %>%
    mutate(
      NH4_diff = NH4[2] - NH4[1],
      NO3_diff = NO3[2] - NO3[1],
      NH4_rate = (NH4_diff * 1) * SoilBulkDensity_avg * soil_depth,
      NO3_rate = (NO3_diff * 1) * SoilBulkDensity_avg * soil_depth,
      Overall_mineralization_rate = NH4_rate + NO3_rate
    ) %>%
    ungroup()
}


N_min_calc <- calculate_rate(N_min_combined)

# Clean up the data frame after calculating rates
N_min_calc <- N_min_calc %>%
  group_by(Sample_ID, Year) %>%
  slice(1) %>%  # Keep only the first row for each Sample_ID and year combination
  ungroup() %>%
  select(Sample_ID, NH4_rate, NO3_rate, Overall_mineralization_rate, Year)  # Keep only the specified columns

## Units are (mg N/cm² per month)

