# Title: Predator Effects Magnified in Ecosystems
# Script: Analyzes ecosystems data
# Authors: N.R. Sommer


# Data Cleaning ----
# Packages
library(dplyr)

# Objects needed to run this script:
  # N-min.R: N_min_calc
  # SIR.R: SIR_avg
  # Veg.R: diversity_calc; functional_groups_wide

# Read in elemental data
CNdata <- read.csv("Data/CN_Tidy_predators.csv")

# Pivot columns
CNdata <- CNdata %>%
  dplyr::select(Year, Sample_ID, Population, Site, CageTreatment, PercentN, PercentC, SampleType) %>%
  pivot_wider(names_from = SampleType, values_from = c(PercentN, PercentC))

# Join datasets by both Sample_ID and Year
combined_data <- full_join(CNdata, N_min_calc, by = c("Sample_ID", "Year"))
combined_data <- full_join(combined_data, SIR_final, by = c("Sample_ID", "Year"))
combined_data <- full_join(combined_data, diversity_calc, by = c("Sample_ID", "Year"))
combined_data <- full_join(combined_data, functional_groups_wide, by = c("Sample_ID", "Year"))

# Define response variables to spread
response_vars <- c(
  "CO2CperHourperg",
  "Overall_mineralization_rate",
  "SORU_Biomass",
  "POPRC_Biomass",
  "MISC_Biomass",
  "PercentN_SOIL",
  "PercentN_LITTER",
  "PercentN_POPRC",
  "PercentN_MISC",
  "PercentN_SORU",
  "PercentC_SOIL",
  "PercentC_LITTER",
  "PercentC_POPRC",
  "PercentC_MISC",
  "PercentC_SORU",
  "PlantDiversity",
  "PlantRichness"
)
# Create wide format with years as suffixes
combined_data <- combined_data %>%
  pivot_wider(
    id_cols = c(Sample_ID, Population, CageTreatment),
    names_from = Year,
    values_from = all_of(response_vars),
    names_sep = "_"
  )

# Relevel for easier model interpretation 
combined_data$CageTreatment <- as.factor(combined_data$CageTreatment)
combined_data$CageTreatment <- relevel(combined_data$CageTreatment, ref = "Vegetation")




# Linear Mixed-Effects Models ----
# Packages
library(lme4)
library(DHARMa)
library(car)
library(boot)

### SORU Biomass----
SORU_model <- lmer(SORU_Biomass_2023 ~ CageTreatment + PercentN_SOIL_2021 + 
                   POPRC_Biomass_2021 + PercentN_SORU_2021 + SORU_Biomass_2021 + 
                   (1 | Population), data = combined_data)
SORU_model_noRE <- lm(SORU_Biomass_2023 ~ CageTreatment + PercentN_SOIL_2021 + 
                      POPRC_Biomass_2021 + PercentN_SORU_2021 + SORU_Biomass_2021, 
                      data = combined_data)
VarCorr(SORU_model)
anova(SORU_model, SORU_model_noRE)
SORU_model <- SORU_model_noRE
simulation_output <- simulateResiduals(fittedModel = SORU_model)
plot(simulation_output)
testDispersion(simulation_output)
testZeroInflation(simulation_output)


### POPRC Biomass----
POPRC_model <- lmer(POPRC_Biomass_2023 ~ CageTreatment + POPRC_Biomass_2021 + 
                    (1 | Population), data = combined_data)
POPRC_model_noRE <- lm(POPRC_Biomass_2023 ~ CageTreatment + POPRC_Biomass_2021, 
                       data = combined_data)
VarCorr(POPRC_model)
anova(POPRC_model, POPRC_model_noRE)
POPRC_model <- POPRC_model_noRE
simulation_output <- simulateResiduals(fittedModel = POPRC_model)
plot(simulation_output)
testDispersion(simulation_output)
testZeroInflation(simulation_output)


### %N SORU ----
SORUN_model <- lmer(PercentN_SORU_2023 ~ SORU_Biomass_2023 + PlantDiversity_2021 + 
                    PercentN_SORU_2021 + (1 | Population), data = combined_data)
SORUN_model_noRE <- lm(PercentN_SORU_2023 ~ SORU_Biomass_2023 + PlantDiversity_2021 + 
                       PercentN_SORU_2021, data = combined_data)
VarCorr(SORUN_model)
anova(SORUN_model, SORUN_model_noRE)
SORUN_model <- SORUN_model_noRE
simulation_output <- simulateResiduals(fittedModel = SORUN_model)
plot(simulation_output)
testDispersion(simulation_output)
testZeroInflation(simulation_output)


### %C SORU ----
SORUC_model <- lmer(PercentC_SORU_2023 ~ SORU_Biomass_2023 + PlantDiversity_2021 + 
                    PercentC_SORU_2021 + (1 | Population), data = combined_data)
VarCorr(SORUC_model)
SORUC_model_noRE <- lm(PercentC_SORU_2023 ~ SORU_Biomass_2023 + PlantDiversity_2021 + 
                       PercentC_SORU_2021, data = combined_data)
SORUC_model <- SORUC_model_noRE
simulation_output <- simulateResiduals(fittedModel = SORUC_model)
plot(simulation_output)
testDispersion(simulation_output)
testZeroInflation(simulation_output)


### %N POPRC ----
POPRCN_model <- lmer(PercentN_POPRC_2023 ~ POPRC_Biomass_2023 + 
                     Overall_mineralization_rate_2021 + PercentN_POPRC_2021 + 
                     (1 | Population), data = combined_data)
VarCorr(POPRCN_model)
POPRCN_model_noRE <- lm(PercentN_POPRC_2023 ~ POPRC_Biomass_2023 + 
                        Overall_mineralization_rate_2021 + PercentN_POPRC_2021, 
                        data = combined_data)
POPRCN_model <- POPRCN_model_noRE
simulation_output <- simulateResiduals(fittedModel = POPRCN_model)
plot(simulation_output)
testDispersion(simulation_output)
testZeroInflation(simulation_output)


### %C POPRC ----
POPRCC_model <- lmer(PercentC_POPRC_2023 ~ POPRC_Biomass_2023 + SORU_Biomass_2023 + 
                     PlantDiversity_2021 + Overall_mineralization_rate_2021 + 
                     POPRC_Biomass_2021 + PercentC_SOIL_2021 + PercentC_POPRC_2021 + 
                     (1 | Population), data = combined_data)
VarCorr(POPRCC_model)
POPRCC_model_noRE <- lm(PercentC_POPRC_2023 ~ POPRC_Biomass_2023 + SORU_Biomass_2023 + 
                        PlantDiversity_2021 + Overall_mineralization_rate_2021 + 
                        POPRC_Biomass_2021 + PercentC_SOIL_2021 + PercentC_POPRC_2021, 
                        data = combined_data)
POPRCC_model <- POPRCC_model_noRE
vif(POPRCC_model)
simulation_output <- simulateResiduals(fittedModel = POPRCC_model)
plot(simulation_output)
testDispersion(simulation_output)
testZeroInflation(simulation_output)




### %C Litter ----
LitterC_model <- lmer(PercentC_LITTER_2023 ~ PercentC_SORU_2023 + PercentC_POPRC_2023 + 
                      (1 | Population), data = combined_data)
                      
VarCorr(LitterC_model) # dtop random effects

LitterC_model_noRE <- lm(PercentC_LITTER_2023 ~ PercentC_SORU_2023 + PercentC_POPRC_2023, 
                         data = combined_data)

LitterC_model <- LitterC_model_noRE

#### Check assumptions ----

vif(LitterC_model) # good

simulation_output <- simulateResiduals(fittedModel = LitterC_model)
plot(simulation_output) # good
testDispersion(simulation_output) # good
testZeroInflation(simulation_output) # good





### %N Litter ----
LitterN_model <- lmer(PercentN_LITTER_2023 ~ PercentN_SORU_2023 + PercentN_POPRC_2023 + 
                        PercentN_LITTER_2021 + (1 | Population), data = combined_data)
#### Check assumptions ----
# Random effects


VarCorr(LitterN_model)

LitterN_model_noRE <- lm(PercentN_LITTER_2023 ~ PercentN_SORU_2023 + PercentN_POPRC_2023 + 
                         PercentN_LITTER_2021, data = combined_data)

anova(LitterN_model, LitterN_model_noRE) # drop random effects

LitterN_model <- LitterN_model_noRE

vif(LitterN_model) # good

simulation_output <- simulateResiduals(fittedModel = LitterN_model)
plot(simulation_output) # good
testDispersion(simulation_output) # good
testZeroInflation(simulation_output) # good




### %C Soil ----
SoilC_model <- lmer(PercentC_SOIL_2023 ~ CageTreatment + SORU_Biomass_2023 + 
                    CO2CperHourperg_2023 + CO2CperHourperg_2021 + PercentC_POPRC_2021 + 
                    PercentC_SOIL_2021 + (1 | Population), data = combined_data)

VarCorr(SoilC_model) # drop random effects

SoilC_model_noRE <- lm(PercentC_SOIL_2023 ~ CageTreatment + SORU_Biomass_2023 + 
                       CO2CperHourperg_2023 + CO2CperHourperg_2021 + PercentC_POPRC_2021 + 
                       PercentC_SOIL_2021, data = combined_data)

SoilC_model <- SoilC_model_noRE

#### Check assumptions ----

vif(SoilC_model) # good
simulation_output <- simulateResiduals(fittedModel = SoilC_model)


plot(simulation_output)
testDispersion(simulation_output)
testZeroInflation(simulation_output) # assumptions not met

# --> non-parametric bootstrapping

boot_fun <- function(data, indices) {
  d <- data[indices, ]  # resample
  mod <- lm(PercentC_SOIL_2023 ~ 
              CageTreatment +
              SORU_Biomass_2023 +
              CO2CperHourperg_2023 +
              CO2CperHourperg_2021 +
              PercentC_POPRC_2021 +
              PercentC_SOIL_2021,
            data = d)
  return(coef(mod))
}


set.seed(1231) # set seed for reproducibility

boot_model <- boot(combined_data, boot_fun, R = 1000) # perform bootstrapping with 1000 resamples


# extract original confidence intervals for lm model using `coef()`
original_ci <- confint(SoilC_model, parm = names(coef(SoilC_model)), level = 0.95)
original_ci <- as.data.frame(original_ci)
colnames(original_ci) <- c("Original_Lower", "Original_Upper")

# initialize data frame
boot_ci_table <- data.frame(
  Fixed_Effect = names(coef(SoilC_model)),
  Boot_Lower = NA,
  Boot_Upper = NA
)

# iterate over fixed effects to get bootstrapped CIs
for (i in seq_along(coef(SoilC_model))) {
  boot_ci_i <- boot.ci(boot_model, type = "perc", index = i)  # bootstrapped CI for all fixed effects
  if (!is.null(boot_ci_i$percent)) {  # ensure boot.ci returned valid percentiles
    boot_ci_table$Boot_Lower[i] <- boot_ci_i$percent[4] 
    boot_ci_table$Boot_Upper[i] <- boot_ci_i$percent[5]
  } else {
    # error handling
    boot_ci_table$Boot_Lower[i] <- NA
    boot_ci_table$Boot_Upper[i] <- NA
  }
}

if (nrow(original_ci) == nrow(boot_ci_table)) {
  ci_comparison <- cbind(Fixed_Effect = rownames(original_ci), original_ci, boot_ci_table[, -1])
} else {
  warning("Mismatch between the number of fixed effects in original CIs and bootstrapped CIs. Check if all fixed effects are present.")
}


# compare bootstrapped estimates with original fixed effect estimates using `coef()` for lm objects
boot_summary <- data.frame(
  Fixed_Effect = names(coef(SoilC_model)),
  Original_Estimate = coef(SoilC_model),
  Bootstrapped_Median = apply(boot_model$t, 2, median),
  Difference = coef(SoilC_model) - apply(boot_model$t, 2, median)
)

print(ci_comparison) 

print(boot_summary) 

# good; estimates are very stable





### %N Soil ----
# %N Soil 2021 variable excluded for standard deviation = 0; failure in bootstrapping
SoilN_model <- lmer(PercentN_SOIL_2023 ~ CageTreatment + PercentN_LITTER_2023 + 
                    SORU_Biomass_2023 + CO2CperHourperg_2021 + (1 | Population), 
                    data = combined_data)
#### Check assumptions ----
SoilN_model_noRE <- lm(PercentN_SOIL_2023 ~ CageTreatment + PercentN_LITTER_2023 + 
                       SORU_Biomass_2023 + CO2CperHourperg_2021, data = combined_data)

VarCorr(SoilN_model)
anova(SoilN_model, SoilN_model_noRE) #drop random effects

SoilN_model <- SoilN_model_noRE


vif(SoilN_model) # good

simulation_output <- simulateResiduals(fittedModel = SoilN_model)
plot(simulation_output) # bad
testDispersion(simulation_output) # fine
testZeroInflation(simulation_output) # bad



# --> non-parametric bootstrapping

boot_fun <- function(data, indices) {
  d <- data[indices, ]  # resample
  mod <- lm(PercentN_SOIL_2023 ~ CageTreatment + 
              PercentN_LITTER_2023 + 
              SORU_Biomass_2023 +
              CO2CperHourperg_2021,
            data = d)
  return(coef(mod))
}


set.seed(1231) # set seed for reproducibility

boot_model <- boot(combined_data, boot_fun, R = 1000) # perform bootstrapping with 1000 resamples


# extract original confidence intervals for lm model using `coef()`
original_ci <- confint(SoilN_model, parm = names(coef(SoilN_model)), level = 0.95)
original_ci <- as.data.frame(original_ci)
colnames(original_ci) <- c("Original_Lower", "Original_Upper")

# initialize data frame
boot_ci_table <- data.frame(
  Fixed_Effect = names(coef(SoilN_model)),
  Boot_Lower = NA,
  Boot_Upper = NA
)

# iterate over fixed effects to get bootstrapped CIs
for (i in seq_along(coef(SoilN_model))) {
  boot_ci_i <- boot.ci(boot_model, type = "perc", index = i)
  if (!is.null(boot_ci_i$percent)) {
    boot_ci_table$Boot_Lower[i] <- boot_ci_i$percent[4] 
    boot_ci_table$Boot_Upper[i] <- boot_ci_i$percent[5]
  } else {
    # error handling
    boot_ci_table$Boot_Lower[i] <- NA
    boot_ci_table$Boot_Upper[i] <- NA
  }
}

if (nrow(original_ci) == nrow(boot_ci_table)) {
  ci_comparison <- cbind(Fixed_Effect = rownames(original_ci), original_ci, boot_ci_table[, -1])
} else {
  warning("Mismatch between the number of fixed effects in original CIs and bootstrapped CIs. Check if all fixed effects are present.")
}


# compare bootstrapped estimates with original estimates
boot_summary <- data.frame(
  Fixed_Effect = names(coef(SoilN_model)),
  Original_Estimate = coef(SoilN_model),
  Bootstrapped_Median = apply(boot_model$t, 2, median, na.rm = TRUE),
  Difference = coef(SoilN_model) - apply(boot_model$t, 2, median, na.rm = TRUE)
)

print(ci_comparison)
print(boot_summary)



### SIR ----
SIR_model <- lmer(CO2CperHourperg_2023 ~ PercentN_SOIL_2023 + PlantDiversity_2021 + 
                    PercentN_SOIL_2021 + PercentN_POPRC_2021 + POPRC_Biomass_2021 + 
                    POPRC_Biomass_2023 + CO2CperHourperg_2021 + (1 | Population), 
                  data = combined_data)

#### Check assumptions ----
SIR_model_noRE <- lm(CO2CperHourperg_2023 ~ PercentN_SOIL_2023 + PlantDiversity_2021 + 
                       PercentN_SOIL_2021 + PercentN_POPRC_2021 + POPRC_Biomass_2021 + 
                       POPRC_Biomass_2023 + CO2CperHourperg_2021, data = combined_data)

VarCorr(SIR_model) # drop random effects
anova(SIR_model, SIR_model_noRE)

SIR_model <- SIR_model_noRE

vif(SIR_model) # good

simulation_output <- simulateResiduals(fittedModel = SIR_model)
plot(simulation_output) # good
testDispersion(simulation_output) # good
testZeroInflation(simulation_output) # good

### Plant Richness ----
PlantRichness_model <- lmer(PlantRichness_2023 ~ PercentC_SOIL_2023 + SORU_Biomass_2023 + 
                              POPRC_Biomass_2023 + PercentN_SORU_2021 + PlantRichness_2021 + 
                              Overall_mineralization_rate_2021 + PercentN_LITTER_2021 + 
                              (1 | Population), data = combined_data)

#### Check assumptions ----
PlantRichness_model_noRE <- lm(PlantRichness_2023 ~ PercentC_SOIL_2023 + SORU_Biomass_2023 + 
                                 POPRC_Biomass_2023 + PercentN_SORU_2021 + PlantRichness_2021 + 
                                 Overall_mineralization_rate_2021 + PercentN_LITTER_2021, 
                               data = combined_data)

VarCorr(PlantRichness_model)
anova(PlantRichness_model, PlantRichness_model_noRE)

PlantRichness_model <- PlantRichness_model_noRE

vif(PlantRichness_model_noRE) # good

simulation_output <- simulateResiduals(fittedModel = PlantRichness_model_noRE)
plot(simulation_output) # fine
testDispersion(simulation_output) # good
testZeroInflation(simulation_output) # good

### N-min ----
Nmin_model <- lmer(Overall_mineralization_rate_2023 ~ CageTreatment + PercentN_SOIL_2023 + 
                     CO2CperHourperg_2023 + PercentN_POPRC_2021 + 
                     Overall_mineralization_rate_2021 + (1 | Population), data = combined_data)

#### Check assumptions ----
Nmin_model_noRE <- lm(Overall_mineralization_rate_2023 ~ CageTreatment + PercentN_SOIL_2023 + 
                        CO2CperHourperg_2023 + PercentN_POPRC_2021 + 
                        Overall_mineralization_rate_2021, data = combined_data)

VarCorr(Nmin_model)
anova(Nmin_model, Nmin_model_noRE)

Nmin_model <- Nmin_model_noRE

vif(Nmin_model) # good

simulation_output <- simulateResiduals(fittedModel = Nmin_model)
plot(simulation_output) # not great
testDispersion(simulation_output) # good
testZeroInflation(simulation_output) # good

# --> non-parametric bootstrapping
boot_fun <- function(data, indices) {
  d <- data[indices, ]  # resample
  mod <- lm(Overall_mineralization_rate_2023 ~ 
              CageTreatment + 
              PercentN_SOIL_2023 + 
              CO2CperHourperg_2023 + 
              PercentN_POPRC_2021 + 
              Overall_mineralization_rate_2021,
            data = d)
  return(coef(mod))
}

set.seed(1231) # set seed for reproducibility
boot_model <- boot(combined_data, boot_fun, R = 1000)

# Extract original confidence intervals
original_ci <- confint(Nmin_model, parm = names(coef(Nmin_model)), level = 0.95)
original_ci <- as.data.frame(original_ci)
colnames(original_ci) <- c("Original_Lower", "Original_Upper")

# Initialize data frame for bootstrapped CIs
boot_ci_table <- data.frame(
  Fixed_Effect = names(coef(Nmin_model)),
  Boot_Lower = NA,
  Boot_Upper = NA
)

# Get bootstrapped CIs
for (i in seq_along(coef(Nmin_model))) {
  boot_ci_i <- boot.ci(boot_model, type = "perc", index = i)
  if (!is.null(boot_ci_i$percent)) {
    boot_ci_table$Boot_Lower[i] <- boot_ci_i$percent[4] 
    boot_ci_table$Boot_Upper[i] <- boot_ci_i$percent[5]
  }
}

# Combine original and bootstrapped CIs
ci_comparison <- cbind(Fixed_Effect = rownames(original_ci), original_ci, boot_ci_table[, -1])

# Compare bootstrapped estimates with original estimates
boot_summary <- data.frame(
  Fixed_Effect = names(coef(Nmin_model)),
  Original_Estimate = coef(Nmin_model),
  Bootstrapped_Median = apply(boot_model$t, 2, median, na.rm = TRUE),
  Difference = coef(Nmin_model) - apply(boot_model$t, 2, median, na.rm = TRUE)
)

print(ci_comparison)
print(boot_summary)



# SEM ----
library(piecewiseSEM)

sem_model <- psem(
  SORU_model,
  POPRC_model,
  SORUN_model,
  SORUC_model,
  POPRCN_model,
  POPRCC_model,
  LitterC_model,
  LitterN_model,
  SoilC_model,
  SoilN_model,
  SIR_model,
  PlantRichness_model,
  Nmin_model,
  data = combined_data)

summary_sem_model <- summary(sem_model)

options(max.print = 10000)
sink("SEM_Summary_Initial.txt") # saved as SEM_Summary_Initial before incorporating tests of directed separation 
print(summary_sem_model)
sink()

# R2 for POPRCN_model and LitterN_model is <0.05 -- drop from global SEM


reduced_sem_model <- psem(
  SORU_model,
  POPRC_model,
  SORUN_model,
  SORUC_model,
  #POPRCN_model,
  POPRCC_model,
  LitterC_model,
  #LitterN_model,
  SoilC_model,
  SoilN_model,
  SIR_model,
  PlantRichness_model,
  Nmin_model,
  data = combined_data)

summary_sem_model <- summary(reduced_sem_model)

options(max.print = 10000)
sink("SEM-Reduced_Summary_Initial.txt") # saved as SEM_Summary_Initial before incorporating tests of directed separation 
print(summary_sem_model)
sink()




### Indirect path assessment ----
coefficients_table <- summary_sem_model$coefficients

coefficients_table$Estimate <- as.numeric(coefficients_table$Estimate) #NAs warning expected

filtered_coefficients <- coefficients_table[coefficients_table$P.Value < 0.05, ] # filter for significance

result <- filtered_coefficients[, c("Response", "Predictor", "Estimate")]



# Step 1: Identify Intermediate Variables
treatment_categories <- c("CageTreatment = Herbivore", 
                          "CageTreatment = Vegetation",  
                          "CageTreatment = F2_Predator",
                          "CageTreatment = F1_Predator")

# Find direct paths from the treatment categories
direct_paths <- result[result$Predictor %in% treatment_categories, ]

# Function to recursively find indirect paths
find_indirect_paths <- function(current_path, current_estimate, current_intermediate) {
  # Find downstream paths from the current intermediate
  downstream <- result[result$Predictor == current_intermediate, ]
  
  # If no downstream paths, return the current path
  if (nrow(downstream) == 0) {
    return(list(current_path))
  }
  
  # List to store all paths
  all_paths <- list()
  
  # Iterate over each downstream path
  for (i in seq_len(nrow(downstream))) {
    next_intermediate <- downstream$Response[i]
    next_estimate <- downstream$Estimate[i]
    
    # Calculate new indirect effect size
    new_estimate <- current_estimate * next_estimate
    
    # Create new path
    new_path <- append(current_path, list(
      list(
        Intermediate = next_intermediate,
        Estimate = next_estimate,
        Indirect_Effect_Size = new_estimate
      )
    ))
    
    # Recursively find further paths
    all_paths <- c(all_paths, find_indirect_paths(new_path, new_estimate, next_intermediate))
  }
  
  return(all_paths)
}

# Step 2: Map Paths for Indirect Influence
# Create a list to store all indirect paths
all_indirect_paths <- list()

# Iterate over each direct path
for (i in seq_len(nrow(direct_paths))) {
  treatment <- direct_paths$Predictor[i]
  intermediate <- direct_paths$Response[i]
  estimate1 <- direct_paths$Estimate[i]
  
  # Initialize the path
  initial_path <- list(
    list(
      Treatment = treatment,
      Intermediate = intermediate,
      Estimate = estimate1,
      Indirect_Effect_Size = estimate1
    )
  )
  
  # Find all indirect paths starting from this direct path
  paths <- find_indirect_paths(initial_path, estimate1, intermediate)
  
  # Add to the list of all indirect paths
  all_indirect_paths <- c(all_indirect_paths, paths)
}

# Step 3: Convert the list of indirect paths to a data frame
indirect_effects_df <- do.call(rbind, lapply(all_indirect_paths, function(path) {
  # Flatten the path into a single row
  data.frame(
    Treatment = path[[1]]$Treatment,
    Intermediates = paste(sapply(path, function(x) x$Intermediate), collapse = " "),
    Indirect_Effect_Size = path[[length(path)]]$Indirect_Effect_Size
  )
}))

# Print the indirect effects
print(indirect_effects_df)


# Define file path
output_file <- "Indirect_paths.txt"

# Write the summary to a text file
write.table(indirect_effects_df, file = output_file, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)


### Indirect paths by trophic impact ----

# First, calculate the trophic impact for each path
filtered_coefficients_trophic <- filtered_coefficients %>%
  select(c(Response, Predictor, Estimate)) %>%
  filter(!is.na(Estimate)) %>%
  group_by(Response) %>%
  mutate(
    # Subtract vegetation estimate from predator and herbivore estimates
    Vegetation_Estimate = ifelse(
      any(Predictor == "CageTreatment = Vegetation"),
      Estimate[Predictor == "CageTreatment = Vegetation"],
      0  # Default value if Vegetation is not present
    ),
    Trophic_Impact = ifelse(
      Predictor %in% c("CageTreatment = F1_Predator", "CageTreatment = F2_Predator", "CageTreatment = Herbivore"),
      Estimate - Vegetation_Estimate,
      Estimate
    ),
    Predictor = case_when(
      Predictor == "CageTreatment = F1_Predator" ~ "Predator F1",
      Predictor == "CageTreatment = F2_Predator" ~ "Predator F2",
      Predictor == "CageTreatment = Herbivore" ~ "Herbivore",
      TRUE ~ Predictor
    )
  ) %>%
  filter(!Predictor %in% c("CageTreatment = Vegetation")) %>%  # Exclude Vegetation paths
  distinct(Predictor, Response, .keep_all = TRUE)  # Ensure unique paths

# Function to recursively find indirect paths based on trophic impact
find_trophic_impact_paths <- function(current_path, current_estimate, current_intermediate) {
  # Find downstream paths from the current intermediate
  downstream <- filtered_coefficients_trophic[filtered_coefficients_trophic$Predictor == current_intermediate, ]
  
  # If no downstream paths, return the current path
  if (nrow(downstream) == 0) {
    return(list(current_path))
  }
  
  # List to store all paths
  all_paths <- list()
  
  # Iterate over each downstream path
  for (i in seq_len(nrow(downstream))) {
    next_intermediate <- downstream$Response[i]
    next_estimate <- downstream$Trophic_Impact[i]
    
    # Calculate new indirect effect size
    new_estimate <- current_estimate * next_estimate
    
    # Create new path
    new_path <- append(current_path, list(
      list(
        Intermediate = next_intermediate,
        Estimate = next_estimate,
        Indirect_Effect_Size = new_estimate
      )
    ))
    
    # Recursively find further paths
    all_paths <- c(all_paths, find_trophic_impact_paths(new_path, new_estimate, next_intermediate))
  }
  
  return(all_paths)
}

# Initialize list for all indirect paths based on trophic impact
all_trophic_impact_paths <- list()

# Iterate over each direct path
for (i in seq_len(nrow(filtered_coefficients_trophic))) {
  treatment <- filtered_coefficients_trophic$Predictor[i]
  intermediate <- filtered_coefficients_trophic$Response[i]
  estimate1 <- filtered_coefficients_trophic$Trophic_Impact[i]
  
  # Initialize the path
  initial_path <- list(
    list(
      Treatment = treatment,
      Intermediate = intermediate,
      Estimate = estimate1,
      Indirect_Effect_Size = estimate1
    )
  )
  
  # Find all indirect paths starting from this direct path
  paths <- find_trophic_impact_paths(initial_path, estimate1, intermediate)
  
  # Add to the list of all indirect paths
  all_trophic_impact_paths <- c(all_trophic_impact_paths, paths)
}

# Convert the list of indirect paths to a data frame
trophic_impact_effects_df <- do.call(rbind, lapply(all_trophic_impact_paths, function(path) {
  # Flatten the path into a single row
  data.frame(
    Treatment = path[[1]]$Treatment,
    Intermediates = paste(sapply(path, function(x) x$Intermediate), collapse = ", "),
    Indirect_Effect_Size = path[[length(path)]]$Indirect_Effect_Size
  )
}))

# Print the indirect effects based on trophic impact
print(trophic_impact_effects_df)

# Define file path
output_file_trophic <- "Trophic_Impact_Indirect_Paths.txt"

# Write the summary to a text file
write.table(trophic_impact_effects_df, file = output_file_trophic, row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)




# Figures ----
# The first set of figures use the basic DiarammeR format. The second set are updated for publication.

library(DiagrammeR)


#### Predators (Basic) ----

filtered_coefficients_predators <- filtered_coefficients %>%
  select(c(Response, Predictor, Estimate)) %>%
  filter(!is.na(Estimate)) %>%
  group_by(Response) %>%
  # Average the estimates for the predator types
  mutate(
    Estimate = ifelse(
      Predictor %in% c("CageTreatment = F1_Predator", "CageTreatment = F2_Predator"),
      mean(Estimate[Predictor %in% c("CageTreatment = F1_Predator", "CageTreatment = F2_Predator")]),
      Estimate
    ),
    Predictor = ifelse(
      Predictor %in% c("CageTreatment = F1_Predator", "CageTreatment = F2_Predator"),
      "Predator",
      Predictor
    )
  ) %>%
      # Remove duplicate rows after averaging
      distinct()


label_mapping <- c(
  "SORU_Biomass_2023" = "Goldenrod biomass",
  "POPRC_Biomass_2023" = "Grass biomass",
  "PercentN_SORU_2023" = "Goldenrod %N",
  "PercentC_SORU_2023" = "Goldenrod %C",
  "PercentN_POPRC_2023" = "Grass %N",
  "PercentC_POPRC_2023" = "Grass %C",
  "PercentC_SOIL_2023" = "Soil %C",
  "PercentN_SOIL_2023" = "Soil %N",
  "PercentC_LITTER_2023" = "Litter %C",
  "CO2CperHourperg_2023" = "SIR",
  "PlantRichness_2023" = "Plant richness",
  "PercentN_LITTER_2023" = "Litter %N",
  "Overall_mineralization_rate_2023" = "Nitrogen mineralization",
  "CageTreatment = Vegetation" = "Vegetation",
  "CageTreatment = Herbivore" = "Herbivore",
  "Predator" = "Predator",
  "SORU_Biomass_2021" = "Baseline goldenrod biomass",
  "POPRC_Biomass_2021" = "Baseline grass biomass",
  "PercentN_SORU_2021" = "Baseline goldenrod %N",
  "PercentC_SORU_2021" = "Baseline goldenrod %C",
  "PercentN_POPRC_2021" = "Baseline grass %N",
  "PercentC_POPRC_2021" = "Baseline grass %C",
  "PercentC_SOIL_2021" = "Baseline soil %C",
  "PercentN_SOIL_2021" = "Baseline soil %N",
  "CO2CperHourperg_2021" = "Baseline SIR",
  "PlantRichness_2021" = "Baseline plant richness",
  "PercentN_LITTER_2021" = "Baseline litter %N",
  "PercentC_LITTER_2021" = "Baseline litter %C",
  "Overall_mineralization_rate_2021" = "Baseline nitrogen mineralization"
)

# Determine edge color based on Estimate
get_edge_color <- function(estimate) {
  if (estimate < 0) {
    return("red")
  } else {
    return("black")
  }
}

# Label mapping
filtered_coefficients_predators_mapped <- filtered_coefficients_predators %>%
  mutate(
    Predictor_Label = label_mapping[Predictor],
    Response_Label = label_mapping[Response]
  )

# Check for any unmapped labels
unmapped_predictors <- filtered_coefficients_predators$Predictor[is.na(filtered_coefficients_predators_mapped$Predictor_Label)]
unmapped_responses <- filtered_coefficients_predators$Response[is.na(filtered_coefficients_predators_mapped$Response_Label)]

if(length(unmapped_predictors) > 0){
  warning("Unmapped Predictors found: ", paste(unmapped_predictors, collapse = ", "))
}

if(length(unmapped_responses) > 0){
  warning("Unmapped Responses found: ", paste(unmapped_responses, collapse = ", "))
}

filtered_coefficients_predators_mapped <- filtered_coefficients_predators_mapped %>%
  filter(!is.na(Predictor_Label) & !is.na(Response_Label))

# Calculate absolute estimates for scaling
filtered_coefficients_predators_mapped <- filtered_coefficients_predators_mapped %>%
  mutate(
    Absolute_Estimate = abs(Estimate)
  )

min_penwidth <- 2
max_penwidth <- 7.5


max_abs_estimate <- max(filtered_coefficients_predators_mapped$Absolute_Estimate, na.rm = TRUE)

# Scale penwidth based on absolute estimate
filtered_coefficients_predators_mapped <- filtered_coefficients_predators_mapped %>%
  mutate(
    Penwidth = (Absolute_Estimate / max_abs_estimate) * (max_penwidth - min_penwidth) + min_penwidth
  )

# Function to find all downstream nodes and paths
find_downstream_paths <- function(start_nodes, all_paths) {
  visited <- c()
  queue <- start_nodes
  
  while (length(queue) > 0) {
    current_node <- queue[1]
    queue <- queue[-1]
    
    visited <- unique(c(visited, current_node))
    
    downstream_paths <- all_paths %>%
      filter(Predictor_Display_Label == current_node)
    
    new_nodes <- setdiff(downstream_paths$Response_Display_Label, visited)
    queue <- unique(c(queue, new_nodes))
  }
  
  all_paths %>%
    filter(Predictor_Display_Label %in% visited | Response_Display_Label %in% visited)
}


filtered_coefficients_predators_mapped <- filtered_coefficients_predators_mapped %>%
  rename(
    Predictor_Display_Label = Predictor_Label,
    Response_Display_Label = Response_Label
  )

# Filter to include only paths starting from "Herbivore" or "Vegetation" or "Predator" and their downstream paths
filtered_paths <- find_downstream_paths(
  start_nodes = c("Predator", "Vegetation", "Herbivore"),
  all_paths = filtered_coefficients_predators_mapped
)

unique_labels <- unique(c(filtered_paths$Predictor_Display_Label, filtered_paths$Response_Display_Label))

dot_script <- "digraph SEM_PathDiagram { 
  rankdir=TB
  node [shape=rectangle style=filled]
  
  { rank=source \"Predator\" \"Herbivore\" \"Baseline grass biomass\" \"Baseline grass %C\" }
  { rank=same \"Grass biomass\" \"Goldenrod biomass\" \"Soil %N\" \"Soil %C\" }
  { rank=same \"Grass %C\" \"Plant richness\" }
  { rank=sink \"Litter %C\" \"Nitrogen mineralization\" }
"

for(label in unique_labels){
  if(grepl("baseline", label, ignore.case = TRUE)){
    dot_script <- paste0(dot_script, sprintf("  \"%s\" [shape=rectangle, style=filled, fillcolor=gray, fontsize=10];\n", label))
  } else {
    dot_script <- paste0(dot_script, sprintf("  \"%s\" [shape=rectangle, style=filled, fillcolor=white];\n", label))
  }
}

for(i in 1:nrow(filtered_paths)){
  predictor <- filtered_paths$Predictor_Display_Label[i]
  response <- filtered_paths$Response_Display_Label[i]
  estimate <- filtered_paths$Estimate[i]
  color <- get_edge_color(estimate)
  penwidth <- filtered_paths$Penwidth[i]
  label <- sprintf("%.2f", estimate)
  
  dot_script <- paste0(dot_script, sprintf("  \"%s\" -> \"%s\" [color=%s, label=\"%s\", penwidth=%.2f];\n",
                                           predictor, response, color, label, penwidth))
}

dot_script <- paste0(dot_script, "}")

grViz(dot_script)




#### Trophic Impact (Basic) ----

# Calculate adjusted estimates for herbivore and predator paths
filtered_coefficients_trophic <- filtered_coefficients %>%
  select(c(Response, Predictor, Estimate)) %>%
  filter(!is.na(Estimate)) %>%
  group_by(Response) %>%
  mutate(
    # Average the estimates for the predator types
    Predator_Estimate = ifelse(
      Predictor %in% c("CageTreatment = F1_Predator", "CageTreatment = F2_Predator"),
      mean(Estimate[Predictor %in% c("CageTreatment = F1_Predator", "CageTreatment = F2_Predator")]),
      Estimate
    ),
    # Subtract vegetation estimate from herbivore and predator estimates
    Vegetation_Estimate = ifelse(
      any(Predictor == "CageTreatment = Vegetation"),
      Estimate[Predictor == "CageTreatment = Vegetation"],
      0
    ),
    Adjusted_Estimate = ifelse(
      Predictor %in% c("CageTreatment = Herbivore", "CageTreatment = F1_Predator", "CageTreatment = F2_Predator"),
      Predator_Estimate - Vegetation_Estimate,
      Estimate
    ),
    Predictor = ifelse(
      Predictor %in% c("CageTreatment = F1_Predator", "CageTreatment = F2_Predator"),
      "Predator",
      Predictor
    )
  ) %>%
  distinct(Predictor, Response, .keep_all = TRUE)

filtered_coefficients_trophic_mapped <- filtered_coefficients_trophic %>%
  mutate(
    Predictor_Label = label_mapping[Predictor],
    Response_Label = label_mapping[Response]
  ) %>%
  rename(
    Predictor_Display_Label = Predictor_Label,
    Response_Display_Label = Response_Label
  )

filtered_coefficients_trophic_mapped <- filtered_coefficients_trophic_mapped %>%
  filter(Predictor_Display_Label != "Vegetation" & Response_Display_Label != "Vegetation")

unmapped_predictors <- filtered_coefficients_trophic$Predictor[is.na(filtered_coefficients_trophic_mapped$Predictor_Display_Label)]
unmapped_responses <- filtered_coefficients_trophic$Response[is.na(filtered_coefficients_trophic_mapped$Response_Display_Label)]

if(length(unmapped_predictors) > 0){
  warning("Unmapped Predictors found: ", paste(unmapped_predictors, collapse = ", "))
}

if(length(unmapped_responses) > 0){
  warning("Unmapped Responses found: ", paste(unmapped_responses, collapse = ", "))
}

filtered_coefficients_trophic_mapped <- filtered_coefficients_trophic_mapped %>%
  mutate(
    Absolute_Estimate = abs(Adjusted_Estimate)
  )

min_penwidth <- 2
max_penwidth <- 7.5

max_abs_estimate <- max(filtered_coefficients_trophic_mapped$Absolute_Estimate, na.rm = TRUE)

filtered_coefficients_trophic_mapped <- filtered_coefficients_trophic_mapped %>%
  mutate(
    Penwidth = (Absolute_Estimate / max_abs_estimate) * (max_penwidth - min_penwidth) + min_penwidth
  )

filtered_paths <- find_downstream_paths(
  start_nodes = c("Predator", "Herbivore"),  # Removed "Vegetation" from start nodes
  all_paths = filtered_coefficients_trophic_mapped
)

unique_labels <- unique(c(filtered_paths$Predictor_Display_Label, filtered_paths$Response_Display_Label))

dot_script <- "digraph SEM_PathDiagram { 
  rankdir=TB
  node [shape=rectangle style=filled]
  
  { rank=source \"Predator\" \"Herbivore\" \"Baseline grass biomass\" \"Baseline grass %C\" }
  { rank=same \"Grass biomass\" \"Goldenrod biomass\" \"Soil %N\" \"Soil %C\" }
  { rank=same \"Grass %C\" \"Plant richness\" }
  { rank=sink \"Litter %C\" \"Nitrogen mineralization\" }
"

for(label in unique_labels){
  if(grepl("baseline", label, ignore.case = TRUE)){
    dot_script <- paste0(dot_script, sprintf("  \"%s\" [shape=rectangle, style=filled, fillcolor=gray, fontsize=10];\n", label))
  } else {
    dot_script <- paste0(dot_script, sprintf("  \"%s\" [shape=rectangle, style=filled, fillcolor=white];\n", label))
  }
}

for(i in 1:nrow(filtered_paths)){
  predictor <- filtered_paths$Predictor_Display_Label[i]
  response <- filtered_paths$Response_Display_Label[i]
  estimate <- filtered_paths$Adjusted_Estimate[i]
  color <- get_edge_color(estimate)
  penwidth <- filtered_paths$Penwidth[i]
  label <- sprintf("%.2f", estimate)
  
  dot_script <- paste0(dot_script, sprintf("  \"%s\" -> \"%s\" [color=%s, label=\"%s\", penwidth=%.2f];\n",
                                           predictor, response, color, label, penwidth))
}


dot_script <- paste0(dot_script, "}")

grViz(dot_script)





#### Generational Effects (Basic) ----

# Calculate adjusted estimates for F1_Predator and F2_Predator paths
filtered_coefficients_generations <- filtered_coefficients %>%
  select(c(Response, Predictor, Estimate)) %>%
  filter(!is.na(Estimate)) %>%
  group_by(Response) %>%
  mutate(
    # Subtract vegetation estimate from predator estimates
    Vegetation_Estimate = ifelse(
      any(Predictor == "CageTreatment = Vegetation"),
      Estimate[Predictor == "CageTreatment = Vegetation"],
      0  # Default value if Vegetation is not present
    ),
    Adjusted_Estimate = ifelse(
      Predictor %in% c("CageTreatment = F1_Predator", "CageTreatment = F2_Predator"),
      Estimate - Vegetation_Estimate,
      Estimate
    ),
    Predictor = case_when(
      Predictor == "CageTreatment = F1_Predator" ~ "Predator F1",
      Predictor == "CageTreatment = F2_Predator" ~ "Predator F2",
      TRUE ~ Predictor
    )
  ) %>%
  filter(!Predictor %in% c("CageTreatment = Herbivore")) %>%  # Exclude Herbivore paths
  distinct(Predictor, Response, .keep_all = TRUE)  # Ensure unique paths

# New label mapping for generational effects
label_mapping_generations <- c(
  "SORU_Biomass_2023" = "Goldenrod biomass",
  "POPRC_Biomass_2023" = "Grass biomass",
  "PercentN_SORU_2023" = "Goldenrod %N",
  "PercentC_SORU_2023" = "Goldenrod %C",
  "PercentN_POPRC_2023" = "Grass %N",
  "PercentC_POPRC_2023" = "Grass %C",
  "PercentC_SOIL_2023" = "Soil %C",
  "PercentN_SOIL_2023" = "Soil %N",
  "PercentC_LITTER_2023" = "Litter %C",
  "CO2CperHourperg_2023" = "SIR",
  "PlantRichness_2023" = "Plant richness",
  "PercentN_LITTER_2023" = "Litter %N",
  "Overall_mineralization_rate_2023" = "Nitrogen mineralization",
  "CageTreatment = Vegetation" = "Vegetation",
  "Predator F1" = "Predator F1",
  "Predator F2" = "Predator F2",
  "SORU_Biomass_2021" = "Baseline goldenrod biomass",
  "POPRC_Biomass_2021" = "Baseline grass biomass",
  "PercentN_SORU_2021" = "Baseline goldenrod %N",
  "PercentC_SORU_2021" = "Baseline goldenrod %C",
  "PercentN_POPRC_2021" = "Baseline grass %N",
  "PercentC_POPRC_2021" = "Baseline grass %C",
  "PercentC_SOIL_2021" = "Baseline soil %C",
  "PercentN_SOIL_2021" = "Baseline soil %N",
  "CO2CperHourperg_2021" = "Baseline SIR",
  "PlantRichness_2021" = "Baseline plant richness",
  "PercentN_LITTER_2021" = "Baseline litter %N",
  "PercentC_LITTER_2021" = "Baseline litter %C",
  "Overall_mineralization_rate_2021" = "Baseline nitrogen mineralization"
)




# Apply label mapping to predictors and responses
filtered_coefficients_generations_mapped <- filtered_coefficients_generations %>%
  mutate(
    Predictor_Label = label_mapping_generations[Predictor],
    Response_Label = label_mapping_generations[Response]
  ) %>%
  rename(
    Predictor_Display_Label = Predictor_Label,
    Response_Display_Label = Response_Label
  )

# Remove paths where "Vegetation" is a predictor
filtered_coefficients_generations_mapped <- filtered_coefficients_generations_mapped %>%
  filter(Predictor_Display_Label != "Vegetation")

# Ensure no NA values in Adjusted_Estimate
filtered_coefficients_generations_mapped <- filtered_coefficients_generations_mapped %>%
  filter(!is.na(Adjusted_Estimate))

# Calculate absolute estimates for scaling
filtered_coefficients_generations_mapped <- filtered_coefficients_generations_mapped %>%
  mutate(
    Absolute_Estimate = abs(Adjusted_Estimate)
  )

# Define minimum and maximum penwidth
min_penwidth <- 2
max_penwidth <- 7.5

# Calculate scaling factor
max_abs_estimate <- max(filtered_coefficients_generations_mapped$Absolute_Estimate, na.rm = TRUE)

# Scale penwidth based on absolute estimate
filtered_coefficients_generations_mapped <- filtered_coefficients_generations_mapped %>%
  mutate(
    Penwidth = (Absolute_Estimate / max_abs_estimate) * (max_penwidth - min_penwidth) + min_penwidth
  )

# Filter out weak paths before rendering
filtered_paths_generations <- find_downstream_paths(
  start_nodes = c("Predator F1", "Predator F2"),  # Start nodes for the generational effects
  all_paths = filtered_coefficients_generations_mapped
) %>%
  # Add filter for weak paths
  filter(abs(Adjusted_Estimate) >= 0.01)

# Check if the data frame is not empty
if (nrow(filtered_paths_generations) > 0) {
  # Extract unique display labels for nodes involved in these paths
  unique_labels_generations <- unique(c(filtered_paths_generations$Predictor_Display_Label, 
                                      filtered_paths_generations$Response_Display_Label))
  
  # Initialize the DOT script with corrected ranks
  dot_script_generations <- "digraph SEM_PathDiagram { 
    rankdir=TB
    node [shape=rectangle style=filled]
    
    { rank=source \"Predator\" \"Herbivore\" \"Baseline grass biomass\" \"Baseline grass %C\" }
    { rank=same \"Grass biomass\" \"Goldenrod biomass\" \"Soil %N\" \"Soil %C\" }
    { rank=same \"Grass %C\" \"Plant richness\" }
    { rank=sink \"Litter %C\" \"Nitrogen mineralization\" }
  "
  
  # Add nodes with shape customization and background color
  for(label in unique_labels_generations){
    if(grepl("baseline", label, ignore.case = TRUE)){
      # Use rectangle shape with dark gray background for baseline nodes
      dot_script_generations <- paste0(dot_script_generations, sprintf("  \"%s\" [shape=rectangle, style=filled, fillcolor=gray, fontsize=10];\n", label))
    } else {
      # Use default rectangle shape with white background for other nodes
      dot_script_generations <- paste0(dot_script_generations, sprintf("  \"%s\" [shape=rectangle, style=filled, fillcolor=white];\n", label))
    }
  }
  
  # Add edges with labels, colors, and penwidth
  for(i in 1:nrow(filtered_paths_generations)){
    predictor <- filtered_paths_generations$Predictor_Display_Label[i]
    response <- filtered_paths_generations$Response_Display_Label[i]
    estimate <- filtered_paths_generations$Adjusted_Estimate[i]
    color <- get_edge_color(estimate)
    penwidth <- filtered_paths_generations$Penwidth[i]
    label <- sprintf("%.2f", estimate)
    
    dot_script_generations <- paste0(dot_script_generations, sprintf("  \"%s\" -> \"%s\" [color=%s, label=\"%s\", penwidth=%.2f];\n",
                                                                 predictor, response, color, label, penwidth))
  }
  
  # Close the DOT script
  dot_script_generations <- paste0(dot_script_generations, "}")
  
  # Render the graph using grViz
  grViz(dot_script_generations)
} else {
  warning("No paths available to render.")
}


## Predator-Herbivore SEM (Pretty) ----
library(DiagrammeR)

# Define custom colors for treatments
treatment_colors <- c(
  "Herbivore" = "#9BA48C",
  "Predator" = "#B5753C"
)

# Define node colors for treatment plot
treatment_baseline_color <- "#A8A8A8"
treatment_response_color <- "#F0F0F0"

# Define nodes for the treatment SEM plot
treatment_nodes <- data.frame(
  name = c(
    "Predator", "Herbivore",
    "Baseline grass biomass", "Baseline grass %C",
    "Grass biomass", "Goldenrod biomass", "Soil %N", "Soil %C",
    "Grass %C", "Plant richness",
    "Litter %C", "Nitrogen mineralization"
  ),
  fillcolor = c(
    treatment_colors["Predator"], treatment_colors["Herbivore"],
    treatment_baseline_color, treatment_baseline_color,
    rep(treatment_response_color, 8)
  ),
  fontcolor = c(
    "white", "white",  # For treatment nodes
    rep("black", 10)   # For all other nodes
  )
)

# Define edges for the treatment SEM plot
treatment_edges <- data.frame(
  from = c(
    "Baseline grass biomass",
    "Predator", "Predator", "Predator", "Predator", 
    "Herbivore", "Herbivore", "Herbivore", "Herbivore",
    "Baseline grass %C",
    "Grass biomass", "Grass %C",
    "Goldenrod biomass"
  ),
  to = c(
    "Grass biomass",
    "Grass biomass", "Goldenrod biomass", "Soil %N", "Soil %C",
    "Goldenrod biomass", "Soil %N", "Soil %C", "Nitrogen mineralization",
    "Soil %C",
    "Grass %C", "Litter %C",
    "Plant richness"
  ),
  value = c(
    0.39,
    8.64, -5.76, -4.65, 0.27,
    12.52, 0.01, -0.04, 16.96,
    0.10,
    -0.04, 0.46,
    -0.02
  )
)

# Initialize the treatment DOT script
treatment_dot_script <- "digraph SEM_PathDiagram { 
  rankdir=TB
  node [shape=rectangle style=filled fontname=\"Arial\"]
  edge [fontname=\"Arial\"]
  
  { rank=source \"Predator\" \"Herbivore\" \"Baseline grass biomass\" \"Baseline grass %C\" }
  { rank=same \"Grass biomass\" \"Goldenrod biomass\" \"Soil %N\" \"Soil %C\" }
  { rank=same \"Grass %C\" \"Plant richness\" }
  { rank=sink \"Litter %C\" \"Nitrogen mineralization\" }
"

# Add nodes with their properties
for(i in 1:nrow(treatment_nodes)) {
  treatment_dot_script <- paste0(treatment_dot_script,
    sprintf("  \"%s\" [fillcolor=\"%s\" fontcolor=\"%s\"];\n",
            treatment_nodes$name[i], 
            treatment_nodes$fillcolor[i],
            treatment_nodes$fontcolor[i]))
}

# Treatment plot edge creation
for(i in 1:nrow(treatment_edges)) {
  predictor <- treatment_edges$from[i]
  response <- treatment_edges$to[i]
  estimate <- treatment_edges$value[i]
  
  abs_estimate <- abs(estimate)
  penwidth <- 1.5 + (abs_estimate / max(abs(treatment_edges$value))) * 5
  
  if(predictor %in% c("Grass biomass", "Goldenrod biomass", "Grass %C", "Plant richness", "Soil %N", "Soil %C")) {
    # For downstream effects, create two brown lines
    # First line with label
    treatment_dot_script <- paste0(treatment_dot_script, 
      sprintf("  \"%s\" -> \"%s\" [color=\"#A18D6B\", label=\" %.2f\", fontcolor=\"#A18D6B\", penwidth=%.2f, style=\"%s\"];\n",
              predictor, response, estimate, penwidth,
              if(estimate < 0) "dashed" else "solid"))
    # Second line without label
    treatment_dot_script <- paste0(treatment_dot_script, 
      sprintf("  \"%s\" -> \"%s\" [color=\"#A18D6B\", label=\"\", penwidth=%.2f, style=\"%s\"];\n",
              predictor, response, penwidth,
              if(estimate < 0) "dashed" else "solid"))
  } else {
    # For treatment and baseline effects
    edge_color <- if(predictor == "Predator") {
      treatment_colors["Predator"]
    } else if(predictor == "Herbivore") {
      treatment_colors["Herbivore"]
    } else {
      "#000000"  # black for baseline effects
    }
    
    treatment_dot_script <- paste0(treatment_dot_script, 
      sprintf("  \"%s\" -> \"%s\" [color=\"%s\", label=\" %.2f\", fontcolor=\"%s\", penwidth=%.2f, style=\"%s\"];\n",
              predictor, response, edge_color, estimate, edge_color, penwidth,
              if(estimate < 0) "dashed" else "solid"))
  }
}

# Close the treatment DOT script
treatment_dot_script <- paste0(treatment_dot_script, "}")

# Render the treatment plot
treatment_plot <- grViz(treatment_dot_script)

# Display plot
treatment_plot

## Generational Effects (Pretty)----

# Define custom colors for generations
generation_colors <- c(
  "G1" = "#DEB887", 
  "G2" = "#8B4513"
)


# Define node colors for generation plot
generation_baseline_color <- "#A8A8A8"
generation_response_color <- "#F0F0F0"

# Define nodes for the generation SEM plot
generation_nodes <- data.frame(
  name = c(
    "G1", "G2",
    "Baseline grass biomass", "Baseline grass %C",
    "Grass biomass", "Goldenrod biomass", "Soil %N", "Soil %C",
    "Grass %C", "Plant richness",
    "Litter %C"
  ),
  label = c(
    "G<sub>1</sub>", "G<sub>2</sub>",
    "Baseline grass biomass", "Baseline grass %C",
    "Grass biomass", "Goldenrod biomass", "Soil %N", "Soil %C",
    "Grass %C", "Plant richness",
    "Litter %C"
  ),
  fillcolor = c(
    generation_colors["G1"], generation_colors["G2"],
    generation_baseline_color, generation_baseline_color,
    rep(generation_response_color, 7)
  ),
  fontcolor = c(
    "white", "white",  # For generation nodes
    rep("black", 9)    # For all other nodes
  )
)

# Define edges for the generational SEM plot
generation_edges <- data.frame(
  from = c(
    "Baseline grass biomass",
    "Baseline grass biomass",
    "G2", "G2", "G2",
    "G1", "G1", "G1",
    "Baseline grass %C",
    "Grass biomass", "Grass %C",
    "Goldenrod biomass"
  ),
  to = c(
    "Grass biomass",
    "Soil %N",
    "Grass biomass", "Goldenrod biomass", "Soil %C",
    "Grass biomass", "Goldenrod biomass", "Soil %C",
    "Soil %C",
    "Grass %C", "Litter %C",
    "Plant richness"
  ),
  value = c(
    0.39,
    0.02,
    11.36, -9.59, 0.31,
    5.93, -1.93, 0.23,
    0.10,
    -0.04, 0.46,
    -0.02
  )
)

# Initialize the generation DOT script
generation_dot_script <- "digraph SEM_PathDiagram { 
  rankdir=TB
  node [shape=rectangle style=filled fontname=\"Arial\"]
  edge [fontname=\"Arial\"]
  
  { rank=source \"G1\" \"G2\" \"Baseline grass biomass\" \"Baseline grass %C\" }
  { rank=same \"Grass biomass\" \"Goldenrod biomass\" \"Soil %N\" \"Soil %C\" }
  { rank=same \"Grass %C\" \"Plant richness\" }
  { rank=sink \"Litter %C\" }
"

# Add nodes with their properties
for(i in 1:nrow(generation_nodes)) {
  generation_dot_script <- paste0(generation_dot_script,
    sprintf("  \"%s\" [label=<%s>, fillcolor=\"%s\" fontcolor=\"%s\"]\n",
            generation_nodes$name[i], 
            generation_nodes$label[i],
            generation_nodes$fillcolor[i],
            generation_nodes$fontcolor[i]))
}

# Generation plot edge creation
for(i in 1:nrow(generation_edges)) {
  predictor <- generation_edges$from[i]
  response <- generation_edges$to[i]
  estimate <- generation_edges$value[i]
  
  abs_estimate <- abs(estimate)
  penwidth <- 1.5 + (abs_estimate / max(abs(generation_edges$value))) * 5
  
  if(predictor %in% c("Grass biomass", "Goldenrod biomass", "Grass %C", "Plant richness", "Soil %N", "Soil %C")) {
    # For downstream effects, create two brown lines
    # First line with label
    generation_dot_script <- paste0(generation_dot_script, 
      sprintf("  \"%s\" -> \"%s\" [color=\"#B5753C\", label=\" %.2f\", fontcolor=\"#B5753C\", penwidth=%.2f, style=\"%s\"];\n",
              predictor, response, estimate, penwidth,
              if(estimate < 0) "dashed" else "solid"))
    # Second line without label
    generation_dot_script <- paste0(generation_dot_script, 
      sprintf("  \"%s\" -> \"%s\" [color=\"#B5753C\", label=\"\", penwidth=%.2f, style=\"%s\"];\n",
              predictor, response, penwidth,
              if(estimate < 0) "dashed" else "solid"))
  } else {
    # For generation and baseline effects
    edge_color <- if(predictor == "G1") {
      generation_colors["G1"]
    } else if(predictor == "G2") {
      generation_colors["G2"]
    } else {
      "#000000"  # black for baseline effects
    }
    
    generation_dot_script <- paste0(generation_dot_script, 
      sprintf("  \"%s\" -> \"%s\" [color=\"%s\", label=\" %.2f\", fontcolor=\"%s\", penwidth=%.2f, style=\"%s\"];\n",
              predictor, response, edge_color, estimate, edge_color, penwidth,
              if(estimate < 0) "dashed" else "solid"))
  }
}

# Close the generation DOT script
generation_dot_script <- paste0(generation_dot_script, "}")

# Render the generation plot
generation_plot <- grViz(generation_dot_script)

# Display plot
generation_plot

#png("Figures/SEM_treatment.png", width = 2000, height = 2400, res = 300)
#treatment_plot
#dev.off()

#png("Figures/SEM_generation.png", width = 2000, height = 2400, res = 300)
#generation_plot
#dev.off()

## Supporting Information ----
##### 2022 Vegetation Comparison Between Predator F1 and F2 Cages ------

library(ggplot2)
library(dplyr)
library(tidyr)
library(broom)
library(purrr)
library(car)
library(rstatix)


veg_cover_2022 <- read.csv("Data/Vegetation.cover_predators.csv") %>%
  filter(Year == 2022) %>%
  select(-BARE) # Remove bare ground as per Veg.R processing

# Process vegetation data following the same approach as Veg.R
veg_processed_2022 <- veg_cover_2022 %>%
  rowwise() %>%
  mutate(Cover_Sum = sum(c_across(-c(Sample_ID, Year)))) %>% # sum plant cover 
  ungroup()

# Normalize cover to 100%
veg_processed_2022 <- veg_processed_2022 %>% 
  mutate(across(
    .cols = -c(Sample_ID, Year, Cover_Sum),
    .fns = ~ round((.x / Cover_Sum) * 100)
  ))

veg_long_2022 <- veg_processed_2022 %>%
  select(-Cover_Sum) %>% 
  pivot_longer(
    cols = -c(Sample_ID, Year),
    names_to = "Species_ID",
    values_to = "Cover"
  )

# Create functional groups following Veg.R approach
functional_groups_2022 <- veg_long_2022 %>%
  group_by(Sample_ID, Year) %>% 
  summarise(
    SORU = sum(Cover[Species_ID %in% c("SORU2", "SOCA6")], na.rm = TRUE),
    POPRC = sum(Cover[Species_ID == "POPRC"], na.rm = TRUE),
    MISC = sum(Cover[!(Species_ID %in% c("SORU2", "SOCA6", "POPRC"))], na.rm = TRUE)
  ) %>%
  ungroup() %>% 
  drop_na()

# Calculate total vegetation cover and plant richness
# For total cover, use original data (before normalization) to include bare ground
total_cover_richness_2022 <- veg_cover_2022 %>%
  rowwise() %>%
  mutate(
    Total_Cover = sum(c_across(-c(Sample_ID, Year))),
    Plant_Richness = sum(c_across(-c(Sample_ID, Year)) > 0)
  ) %>%
  ungroup() %>%
  select(Sample_ID, Year, Total_Cover, Plant_Richness)

veg_summary_2022 <- functional_groups_2022 %>%
  left_join(total_cover_richness_2022, by = c("Sample_ID", "Year")) %>%
  rename(
    TotalCover = Total_Cover,
    PlantRichness = Plant_Richness
  )

# Using CN data to add treatment column
cage_treatments <- read.csv("Data/CN_Tidy_predators.csv") %>%
  select(Sample_ID, CageTreatment) %>%
  distinct() %>%
  filter(CageTreatment %in% c("F1_Predator", "F2_Predator"))

veg_analysis_2022 <- veg_summary_2022 %>%
  inner_join(cage_treatments, by = "Sample_ID") %>%
  filter(CageTreatment %in% c("F1_Predator", "F2_Predator"))

# Treatment labels for plotting
veg_analysis_2022 <- veg_analysis_2022 %>%
  mutate(
    Treatment = case_when(
      CageTreatment == "F1_Predator" ~ "G1",
      CageTreatment == "F2_Predator" ~ "G2"
    ),
    Treatment = factor(Treatment, levels = c("G1", "G2"))
  )


perform_statistical_test <- function(data, variable) {
  # Check normality
  f1_data <- data %>% filter(Treatment == "G1") %>% pull(!!sym(variable))
  f2_data <- data %>% filter(Treatment == "G2") %>% pull(!!sym(variable))
  
  # Check if data has variation (not all identical values)
  f1_var <- var(f1_data, na.rm = TRUE)
  f2_var <- var(f2_data, na.rm = TRUE)
  
  # If either group has no variation, skip normality test and use Wilcoxon
  if (f1_var == 0 || f2_var == 0 || is.na(f1_var) || is.na(f2_var)) {
    test_result <- wilcox.test(as.formula(paste(variable, "~ Treatment")), 
                              data = data, 
                              exact = FALSE)
    test_type <- "Wilcoxon"
  } else {
    # Shapiro-Wilk test for normality (only if data has variation)
    f1_norm <- shapiro.test(f1_data)
    f2_norm <- shapiro.test(f2_data)
    
    # If both groups are normal, use t-test, otherwise use Wilcoxon
    if (f1_norm$p.value > 0.05 && f2_norm$p.value > 0.05) {
      test_result <- t.test(as.formula(paste(variable, "~ Treatment")), data = data)
      test_type <- "t-test"
    } else {
      # Using exact = FALSE to handle ties properly in Wilcoxon test
      test_result <- wilcox.test(as.formula(paste(variable, "~ Treatment")), 
                                data = data, 
                                exact = FALSE)
      test_type <- "Wilcoxon"
    }
  }
  
  return(list(
    test_type = test_type,
    statistic = test_result$statistic,
    p_value = test_result$p.value,
    f1_mean = mean(f1_data, na.rm = TRUE),
    f2_mean = mean(f2_data, na.rm = TRUE),
    f1_sd = sd(f1_data, na.rm = TRUE),
    f2_sd = sd(f2_data, na.rm = TRUE)
  ))
}

variables_to_test <- c("SORU", "POPRC", "MISC", "TotalCover", "PlantRichness")
statistical_results <- list()

for (var in variables_to_test) {
  statistical_results[[var]] <- perform_statistical_test(veg_analysis_2022, var)
}


generation_colors_supporting <- c(
  "G1" = "#DEB887", 
  "G2" = "#8B4513"
)

# Create individual plots for each variable
create_veg_plot <- function(data, variable, title, y_label) {
  ggplot(data, aes(x = Treatment, y = !!sym(variable), fill = Treatment)) +
    # Box plots
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(aes(color = Treatment), 
                alpha = 0.6, 
                size = 2, 
                width = 0.2,
                height = 0) +
    # Mean points (white diamonds)
    stat_summary(fun = mean, geom = "point", 
                shape = 23, size = 3, fill = "white", color = "black") +
    # Styling
    scale_fill_manual(values = generation_colors_supporting) +
    scale_color_manual(values = generation_colors_supporting) +
    coord_flip() +
    theme_light() +
    theme(
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_line(color = "grey95"),
      legend.position = "none",
      axis.text.y = element_text(size = 10),
      axis.title.y = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12)
    ) +
    labs(
      title = title,
      x = y_label
    )
}

soru_plot <- create_veg_plot(veg_analysis_2022, "SORU", 
                            "SORU Functional Group Cover", "Cover (%)")
poprc_plot <- create_veg_plot(veg_analysis_2022, "POPRC", 
                             "POPRC Functional Group Cover", "Cover (%)")
misc_plot <- create_veg_plot(veg_analysis_2022, "MISC", 
                            "MISC Functional Group Cover", "Cover (%)")
total_cover_plot <- create_veg_plot(veg_analysis_2022, "TotalCover", 
                                   "Total Vegetation Cover", "Cover (%)")
richness_plot <- create_veg_plot(veg_analysis_2022, "PlantRichness", 
                                "Plant Species Richness", "Number of Species")

# Combine plots using patchwork
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  
  # Create a table plot for statistical results
  if (requireNamespace("ggpubr", quietly = TRUE)) {
    library(ggpubr)
    
    # Prepare table data for display
    table_data <- statistical_summary %>%
      mutate(
        G1_Mean_SD = paste0(F1_Mean, " ± ", F1_SD),
        G2_Mean_SD = paste0(F2_Mean, " ± ", F2_SD),
        P_Value_Formatted = ifelse(P_Value < 0.001, "< 0.001", 
                                  ifelse(P_Value < 0.01, sprintf("%.3f", P_Value),
                                         sprintf("%.3f", P_Value)))
      ) %>%
      select(Variable, G1_Mean_SD, G2_Mean_SD, Test_Type, P_Value_Formatted, Significance)
    
    # Create table plot
    results_table <- ggtexttable(table_data, 
                                rows = NULL,
                                theme = ttheme("light", 
                                              base_size = 9,
                                              padding = unit(c(2, 3), "mm"))) +
      theme(plot.margin = margin(5, 5, 5, 5))
    
    # Create combined figure with table
    combined_veg_plot <- (soru_plot + poprc_plot) / 
                        (misc_plot + total_cover_plot) / 
                        (richness_plot + results_table) +
      plot_annotation(
        theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
                     plot.subtitle = element_text(hjust = 0.5, size = 12))
      )
    
  } else {
    # Fallback without table if ggpubr not available
    combined_veg_plot <- (soru_plot + poprc_plot) / 
                        (misc_plot + total_cover_plot) / 
                        (richness_plot + plot_spacer()) +
      plot_annotation(
        theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
                     plot.subtitle = element_text(hjust = 0.5, size = 12))
      )
  }
  
  combined_veg_plot
} else {
  # Display plots individually
  print("SORU Functional Group Cover:")
  print(soru_plot)
  print("POPRC Functional Group Cover:")
  print(poprc_plot)
  print("MISC Functional Group Cover:")
  print(misc_plot)
  print("Total Vegetation Cover:")
  print(total_cover_plot)
  print("Plant Species Richness:")
  print(richness_plot)
}

# Summary table of statistical results
statistical_summary <- data.frame(
  Variable = variables_to_test,
  Test_Type = sapply(statistical_results, function(x) x$test_type),
  F1_Mean = sapply(statistical_results, function(x) round(x$f1_mean, 2)),
  F1_SD = sapply(statistical_results, function(x) round(x$f1_sd, 2)),
  F2_Mean = sapply(statistical_results, function(x) round(x$f2_mean, 2)),
  F2_SD = sapply(statistical_results, function(x) round(x$f2_sd, 2)),
  Statistic = sapply(statistical_results, function(x) round(x$statistic, 3)),
  P_Value = sapply(statistical_results, function(x) round(x$p_value, 4)),
  Significance = sapply(statistical_results, function(x) 
    ifelse(x$p_value < 0.001, "***", 
           ifelse(x$p_value < 0.01, "**", 
                  ifelse(x$p_value < 0.05, "*", "ns"))))
)

print(statistical_summary)



##### 2021 Starting Conditions Comparison Between G1 and G2 Cages ------


# Map cage treatments
cage_treatments_2021 <- read.csv("Data/CN_Tidy_predators.csv") %>%
  filter(Year == 2021) %>%
  select(Sample_ID, CageTreatment) %>%
  distinct()

# Cover comparison for 2021
veg_cover_2021 <- read.csv("Data/Vegetation.cover_predators.csv") %>%
  filter(Year == 2021) %>%
  rowwise() %>%
  mutate(Cover_Sum = sum(c_across(-c(Sample_ID, Year)))) %>% 
  ungroup() %>%
  mutate(across(
    .cols = -c(Sample_ID, Year, Cover_Sum),
    .fns = ~ round((.x / Cover_Sum) * 100)
  )) %>% # normalize cover to 100%
  select(Sample_ID, SORU2, POPRC) %>%
  left_join(cage_treatments_2021, by = "Sample_ID") %>%
  mutate(
    Treatment = case_when(
      CageTreatment == "Vegetation" ~ "Vegetation",
      CageTreatment == "Herbivore" ~ "Herbivore",
      CageTreatment == "F1_Predator" ~ "G1",
      CageTreatment == "F2_Predator" ~ "G2"
    ),
    Treatment = factor(Treatment, levels = c("Vegetation", "Herbivore", "G1", "G2"))
  ) %>%
  drop_na()



# Define colors for all treatments
veg_cover_colors <- c(
  "Vegetation" = "#D3D3D3",  
  "Herbivore" = "#9BA48C", 
  "G1" = "#DEB887",          
  "G2" = "#8B4513"          
)

veg_cover_plot <- ggplot(veg_cover_2021 %>% 
                          pivot_longer(cols = c(SORU2, POPRC), 
                                     names_to = "Species", 
                                     values_to = "Cover") %>%
                          mutate(Species = case_when(Species == "SORU2" ~ "SORU", TRUE ~ Species)),
                        aes(x = Cover, y = Species, fill = Treatment)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = Treatment), 
              alpha = 0.6, 
              size = 2, 
              position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.3, jitter.height = 0.2)) +
  stat_summary(fun = mean, geom = "point", 
              shape = 23, size = 3, fill = "white", color = "black",
              position = position_dodge(width = 0.8)) +
  scale_fill_manual(values = veg_cover_colors) +
  scale_color_manual(values = veg_cover_colors) +
  theme_light() +
  theme(
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_line(color = "grey95"),
    legend.position = "top",
    axis.text.y = element_text(size = 10),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 12)
  ) +
  labs(
    x = "Percent Cover (%)",
    y = "Species",
    fill = "Treatment",
    color = "Treatment"
  )

print(veg_cover_plot)



# Filter combined_data for 2021 baseline variables and G1/G2 treatments
baseline_2021 <- combined_data %>%
  filter(CageTreatment %in% c("F1_Predator", "F2_Predator")) %>%
  select(
    Sample_ID, CageTreatment,
    SORU_Biomass_2021, POPRC_Biomass_2021, MISC_Biomass_2021,
    PlantDiversity_2021, PlantRichness_2021,
    PercentN_SOIL_2021, PercentC_SOIL_2021,
    PercentN_LITTER_2021, PercentC_LITTER_2021,
    PercentN_SORU_2021, PercentC_SORU_2021,
    PercentN_POPRC_2021, PercentC_POPRC_2021,
    PercentN_MISC_2021, PercentC_MISC_2021,
    Overall_mineralization_rate_2021, CO2CperHourperg_2021
  ) %>%
  drop_na()

baseline_2021 <- baseline_2021 %>%
  mutate(
    Treatment = case_when(
      CageTreatment == "F1_Predator" ~ "G1",
      CageTreatment == "F2_Predator" ~ "G2"
    ),
    Treatment = factor(Treatment, levels = c("G1", "G2"))
  )


baseline_variables <- c(
  "SORU_Biomass_2021", "POPRC_Biomass_2021", "MISC_Biomass_2021",
  "PlantDiversity_2021", "PlantRichness_2021",
  "PercentN_SOIL_2021", "PercentC_SOIL_2021",
  "PercentN_LITTER_2021", "PercentC_LITTER_2021",
  "PercentN_SORU_2021", "PercentC_SORU_2021",
  "PercentN_POPRC_2021", "PercentC_POPRC_2021",
  "PercentN_MISC_2021", "PercentC_MISC_2021",
  "Overall_mineralization_rate_2021", "CO2CperHourperg_2021"
)

# Function to screen variables for practical equivalence
screen_variables <- function(data, variables) {
  screening_results <- list()
  
  for (var in variables) {
    f1_data <- data %>% filter(Treatment == "G1") %>% pull(!!sym(var))
    f2_data <- data %>% filter(Treatment == "G2") %>% pull(!!sym(var))
    
    # Calculate basic statistics
    f1_mean <- mean(f1_data, na.rm = TRUE)
    f2_mean <- mean(f2_data, na.rm = TRUE)
    f1_sd <- sd(f1_data, na.rm = TRUE)
    f2_sd <- sd(f2_data, na.rm = TRUE)
    
    # Check for practical equivalence
    mean_diff <- abs(f1_mean - f2_mean)
    pooled_sd <- sqrt((f1_sd^2 + f2_sd^2) / 2)
    
    # Calculate effect size (Cohen's d)
    effect_size <- ifelse(pooled_sd > 0, mean_diff / pooled_sd, 0)
    
    # Determine if practically equivalent
    # Criteria: effect size < 0.1 (very small effect) OR relative difference < 2%
    # This is more conservative - only skip testing if differences are truly minimal
    relative_diff <- ifelse(mean(f1_mean, f2_mean) > 0, 
                           mean_diff / mean(f1_mean, f2_mean) * 100, 0)
    
    is_equivalent <- effect_size < 0.1 || relative_diff < 2
    
    screening_results[[var]] <- list(
      variable = var,
      f1_mean = f1_mean,
      f2_mean = f2_mean,
      f1_sd = f1_sd,
      f2_sd = f2_sd,
      mean_diff = mean_diff,
      effect_size = effect_size,
      relative_diff = relative_diff,
      is_equivalent = is_equivalent,
      should_test = !is_equivalent
    )
  }
  
  return(screening_results)
}

# Screen baseline variables
baseline_screening <- screen_variables(baseline_2021, baseline_variables)

# Display screening results
screening_summary <- data.frame(
  Variable = sapply(baseline_screening, function(x) x$variable),
  G1_Mean = sapply(baseline_screening, function(x) round(x$f1_mean, 3)),
  G2_Mean = sapply(baseline_screening, function(x) round(x$f2_mean, 3)),
  Effect_Size = sapply(baseline_screening, function(x) round(x$effect_size, 3)),
  Relative_Diff_Percent = sapply(baseline_screening, function(x) round(x$relative_diff, 2)),
  Practically_Equivalent = sapply(baseline_screening, function(x) x$is_equivalent),
  Should_Test = sapply(baseline_screening, function(x) x$should_test)
)

print(screening_summary)

# Identify variables that need statistical testing
variables_to_test <- names(baseline_screening)[sapply(baseline_screening, function(x) x$should_test)]
variables_equivalent <- names(baseline_screening)[sapply(baseline_screening, function(x) x$is_equivalent)]

baseline_statistical_results <- list()

# Only perform statistical tests on variables that aren't practically equivalent
for (var in baseline_variables) {
  if (baseline_screening[[var]]$should_test) {
    tryCatch({
      baseline_statistical_results[[var]] <- perform_statistical_test(baseline_2021, var)
    }, error = function(e) {
      # Fallback for any remaining errors
      f1_data <- baseline_2021 %>% filter(Treatment == "G1") %>% pull(!!sym(var))
      f2_data <- baseline_2021 %>% filter(Treatment == "G2") %>% pull(!!sym(var))
      
      test_result <- wilcox.test(as.formula(paste(var, "~ Treatment")), 
                                data = baseline_2021, 
                                exact = FALSE)
      baseline_statistical_results[[var]] <<- list(
        test_type = "Wilcoxon",
        statistic = test_result$statistic,
        p_value = test_result$p.value,
        f1_mean = mean(f1_data, na.rm = TRUE),
        f2_mean = mean(f2_data, na.rm = TRUE),
        f1_sd = sd(f1_data, na.rm = TRUE),
        f2_sd = sd(f2_data, na.rm = TRUE)
      )
    })
  } else {
    # For practically equivalent variables, create results without testing
    screening <- baseline_screening[[var]]
    baseline_statistical_results[[var]] <- list(
      test_type = "Equivalent",
      statistic = NA,
      p_value = NA,
      f1_mean = screening$f1_mean,
      f2_mean = screening$f2_mean,
      f1_sd = screening$f1_sd,
      f2_sd = screening$f2_sd
    )
  }
}


# Create plots for key baseline variables (selecting most important ones for figure)
create_baseline_plot <- function(data, variable, title, y_label) {
  ggplot(data, aes(x = Treatment, y = !!sym(variable), fill = Treatment)) +
    # Box plots
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    # Jittered data points
    geom_jitter(aes(color = Treatment), 
                alpha = 0.6, 
                size = 2, 
                width = 0.2,
                height = 0) +
    # Mean points (white diamonds)
    stat_summary(fun = mean, geom = "point", 
                shape = 23, size = 3, fill = "white", color = "black") +
    # Styling
    scale_fill_manual(values = generation_colors_supporting) +
    scale_color_manual(values = generation_colors_supporting) +
    coord_flip() +
    theme_light() +
    theme(
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_line(color = "grey95"),
      legend.position = "none",
      axis.text.y = element_text(size = 10),
      axis.title.y = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12)
    ) +
    labs(
      title = title,
      x = y_label
    )
}

# Group baseline variables by type for organized plotting
vegetation_vars <- c("SORU_Biomass_2021", "POPRC_Biomass_2021", "MISC_Biomass_2021",
                    "PlantDiversity_2021", "PlantRichness_2021")

soil_chemistry_vars <- c("PercentN_SOIL_2021", "PercentC_SOIL_2021",
                        "PercentN_LITTER_2021", "PercentC_LITTER_2021",
                        "PercentN_SORU_2021", "PercentC_SORU_2021",
                        "PercentN_POPRC_2021", "PercentC_POPRC_2021",
                        "PercentN_MISC_2021", "PercentC_MISC_2021")

ecosystem_function_vars <- c("Overall_mineralization_rate_2021", "CO2CperHourperg_2021")

# Function to create plots for a group of variables
create_variable_group_plots <- function(data, variables, group_name) {
  plots <- list()
  
  for (i in seq_along(variables)) {
    var <- variables[i]
    
    # Create descriptive title and label
    var_name <- gsub("_2021", "", var)
    var_name <- gsub("_", " ", var_name)
    var_name <- tools::toTitleCase(var_name)
    
    # Update specific variable names
    if (var == "Overall_mineralization_rate_2021") {
      var_name <- "N-min"
    } else if (var == "CO2CperHourperg_2021") {
      var_name <- "SIR"
    }
    
    # Determine appropriate y-axis label based on variable type
    if (grepl("Percent", var)) {
      y_label <- "Percentage (%)"
    } else if (grepl("Biomass", var)) {
      y_label <- "Biomass (g)"
    } else if (grepl("Rate", var)) {
      y_label <- "Rate (μg/g/day)"
    } else if (grepl("Diversity", var)) {
      y_label <- "Shannon Index"
    } else if (grepl("Richness", var)) {
      y_label <- "Number of Species"
    } else {
      y_label <- "Value"
    }
    
    plots[[i]] <- create_baseline_plot(data, var, 
                                      var_name, y_label)
  }
  
  return(plots)
}

# Create plots for each group
vegetation_plots <- create_variable_group_plots(baseline_2021, vegetation_vars, "Vegetation")
soil_chemistry_plots <- create_variable_group_plots(baseline_2021, soil_chemistry_vars, "Soil Chemistry")
ecosystem_function_plots <- create_variable_group_plots(baseline_2021, ecosystem_function_vars, "Ecosystem Function")

# Create separate patchwork figures for each variable group with tables
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  
  # Function to create statistical table for a group of variables
  create_group_table <- function(variables, group_name) {
    if (requireNamespace("ggpubr", quietly = TRUE)) {
      library(ggpubr)
      
      # Filter screening results for this group
      group_screening <- baseline_screening[names(baseline_screening) %in% variables]
      
      # Prepare table data for this group
      group_table_data <- data.frame(
        Variable = sapply(group_screening, function(x) {
          var_name <- gsub("_2021", "", x$variable)
          var_name <- gsub("_", " ", var_name)
          tools::toTitleCase(var_name)
        }),
        G1_Mean_SD = sapply(names(group_screening), function(var) {
          result <- baseline_statistical_results[[var]]
          paste0(round(result$f1_mean, 3), " ± ", round(result$f1_sd, 3))
        }),
        G2_Mean_SD = sapply(names(group_screening), function(var) {
          result <- baseline_statistical_results[[var]]
          paste0(round(result$f2_mean, 3), " ± ", round(result$f2_sd, 3))
        }),
        Test_Type = sapply(names(group_screening), function(var) {
          result <- baseline_statistical_results[[var]]
          test_type <- result$test_type
          if (test_type == "Practically equivalent") {
            "Equivalent"
          } else if (test_type == "t-test") {
            "t-test"
          } else if (test_type == "Wilcoxon") {
            "Wilcoxon"
          } else {
            test_type
          }
        }),
        P_Value = sapply(names(group_screening), function(var) {
          result <- baseline_statistical_results[[var]]
          p_val <- result$p_value
          if (is.na(p_val)) {
            "-"
          } else if (p_val < 0.001) {
            "< 0.001"
          } else {
            sprintf("%.3f", p_val)
          }
        }),
        Significance = sapply(names(group_screening), function(var) {
          result <- baseline_statistical_results[[var]]
          p_val <- result$p_value
          if (is.na(p_val)) {
            "-"
          } else if (p_val < 0.001) {
            "***"
          } else if (p_val < 0.01) {
            "**"
          } else if (p_val < 0.05) {
            "*"
          } else {
            "ns"
          }
        })
      )
      
      # Create table plot
      group_table <- ggtexttable(group_table_data, 
                                rows = NULL,
                                theme = ttheme("light", 
                                              base_size = 8,
                                              padding = unit(c(1, 2), "mm"))) +
        theme(plot.margin = margin(5, 5, 5, 5))
      
      return(group_table)
    } else {
      return(plot_spacer())
    }
  }
  
  # Function to create patchwork for a group of plots with table
  create_group_patchwork <- function(plots, variables, group_name) {
    n_plots <- length(plots)
    
    # Create table for this group
    group_table <- create_group_table(variables, group_name)
    
    if (n_plots == 0) {
      return(NULL)
    } else if (n_plots == 1) {
      return(plots[[1]] / group_table)
    } else if (n_plots == 2) {
      return((plots[[1]] + plots[[2]]) / group_table)
    } else if (n_plots == 3) {
      return((plots[[1]] + plots[[2]] + plots[[3]]) / group_table)
    } else if (n_plots == 4) {
      return((plots[[1]] + plots[[2]]) / (plots[[3]] + plots[[4]]) / group_table)
    } else if (n_plots == 5) {
      return((plots[[1]] + plots[[2]] + plots[[3]]) / (plots[[4]] + plots[[5]]) / group_table)
    } else {
      # For 6+ plots, use 3x2 layout with table below
      return((plots[[1]] + plots[[2]] + plots[[3]]) / 
             (plots[[4]] + plots[[5]] + plots[[6]]) / 
             group_table)
    }
  }
  
  # Create separate figures for each group
  vegetation_figure <- create_group_patchwork(vegetation_plots, vegetation_vars, "Vegetation")
  soil_chemistry_figure <- create_group_patchwork(soil_chemistry_plots, soil_chemistry_vars, "Soil Chemistry")
  ecosystem_function_figure <- create_group_patchwork(ecosystem_function_plots, ecosystem_function_vars, "Ecosystem Function")
  
  # Display the figures
  print(vegetation_figure)
  
  print(soil_chemistry_figure)
  
  print(ecosystem_function_figure)
  
} else {
  # Display plots individually if patchwork not available
  for (i in seq_along(vegetation_plots)) {
    print(vegetation_plots[[i]])
  }
  
  for (i in seq_along(soil_chemistry_plots)) {
    print(soil_chemistry_plots[[i]])
  }
  
  for (i in seq_along(ecosystem_function_plots)) {
    print(ecosystem_function_plots[[i]])
  }
}


baseline_statistical_summary <- data.frame(
  Variable = baseline_variables,
  Test_Type = sapply(baseline_statistical_results, function(x) x$test_type),
  G1_Mean = sapply(baseline_statistical_results, function(x) round(x$f1_mean, 3)),
  G1_SD = sapply(baseline_statistical_results, function(x) round(x$f1_sd, 3)),
  G2_Mean = sapply(baseline_statistical_results, function(x) round(x$f2_mean, 3)),
  G2_SD = sapply(baseline_statistical_results, function(x) round(x$f2_sd, 3)),
  Statistic = sapply(baseline_statistical_results, function(x) round(x$statistic, 3)),
  P_Value = sapply(baseline_statistical_results, function(x) round(x$p_value, 4)),
  Significance = sapply(baseline_statistical_results, function(x) 
    ifelse(x$p_value < 0.001, "***", 
           ifelse(x$p_value < 0.01, "**", 
                  ifelse(x$p_value < 0.05, "*", "ns"))))
)


