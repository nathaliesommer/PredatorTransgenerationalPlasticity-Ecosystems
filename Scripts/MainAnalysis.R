# GenEffects on Ecosystems:

# Objects needed: 
  # N-min.R: N_min_calc
  # SIR.R: SIR_avg
  # Veg.R: diversity_calc; functional_groups_wide

CNdata <- read.csv("Data/CN_Tidy_predators.csv")

CNdata <- CNdata %>%
  dplyr::select(Sample_ID, Population, Site, CageTreatment, PercentN, PercentC, SampleType) %>%
  pivot_wider(names_from = SampleType, values_from = c(PercentN, PercentC))


# Join datasets
combined_data <- full_join(CNdata, N_min_calc, by = "Sample_ID")

combined_data <- full_join(combined_data, SIR_avg, by = "Sample_ID")
combined_data <- full_join(combined_data, diversity_calc, by = "Sample_ID")
combined_data <- full_join(combined_data, functional_groups_wide, by = "Sample_ID")

duplicates(combined_data$Sample_ID)
