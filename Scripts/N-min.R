### N-mineralization ----

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

N_min <- read.csv("Data/N-min.csv") # NOTE: N-min data was blank-adjusted by the processing lab
BulkDensity <- read.csv("Data/BulkDensity.csv")

# Bulk Density calculation
BulkDensity <- BulkDensity %>% 
  mutate(SoilMass = Weight_SoilTray-Weight_Tray) %>% 
  mutate(SoilMassAdj = SoilMass - (Materials_volume_mL_end - Materials_volume_mL_start)) %>% 
  mutate(SoilBulkDensity = (SoilMassAdj / 308.89)) # units = g/cm^3, where 308.89 is the cm^3 of the sampler



# Bulk Density, unit conversion 
BulkDensity_avg <- BulkDensity %>%
  group_by(Site) %>%
  summarize(SoilBulkDensity_avg = mean(SoilBulkDensity))


calculate_rate <- function(data, year, soil_depth = 10) {
  data %>%
    left_join(BulkDensity_avg, by = "Site") %>%
    group_by(Sample_ID) %>%
    # Filter for the first and last day
    filter(Day == min(Day) | Day == max(Day)) %>%
    arrange(Day) %>%
    # Calculate differences between the last and first day
    mutate(
      NH4_diff = `blank.corrected.N.NH4..mg.per.mL.`[2] - `blank.corrected.N.NH4..mg.per.mL.`[1],
      NO3_diff = `blank.corrected.N.NO3..mg.per.mL.`[2] - `blank.corrected.N.NO3..mg.per.mL.`[1],
      # Calculate rates, scaled to per month and converted to mg/cm³
      NH4_rate = (NH4_diff * 1) * SoilBulkDensity_avg * soil_depth, # mg/cm³ per month
      NO3_rate = (NO3_diff * 1) * SoilBulkDensity_avg * soil_depth, # mg/cm³ per month
      Overall_rate = NH4_rate + NO3_rate
    ) %>%
    ungroup()
}



N_min_calc <- calculate_rate(N_min)

# Clean up the data frame after calculating rates
N_min_calc <- N_min_calc %>%
  filter(Day != 30) %>%  # Drop rows where Day equals 30
  select(Sample_ID, NH4_rate, NO3_rate, Overall_rate)  # Keep only the specified columns

## Units are (mg N/cm² per month)

