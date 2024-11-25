# GenEffects on Ecosystems:


# Objects needed: 
  # N-min.R: N_min_calc
  # SIR.R: SIR_avg
  # Veg.R: diversity_calc; functional_groups_wide

CNdata <- read.csv("Data/CN_Tidy_predators.csv")

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


# Relevel for model interpretation (update after filtering)
combined_data$CageTreatment <- as.factor(combined_data$CageTreatment)
combined_data$CageTreatment <- relevel(combined_data$CageTreatment, ref = "Vegetation")

# Packages
library(lme4)
library(DHARMa)
library(car)
library(boot)


# Linear Mixed-Effects Models ----

### SORU Biomass----
SORU_model <- lmer(SORU_Biomass_2023 ~ CageTreatment + PercentN_SOIL_2021 + 
                   POPRC_Biomass_2021 + PercentN_SORU_2021 + SORU_Biomass_2021 + 
                   (1 | Population), data = combined_data)
#### Check assumptions----
   # Random effects
SORU_model_noRE <- lm(SORU_Biomass_2023 ~ CageTreatment + PercentN_SOIL_2021 + 
                      POPRC_Biomass_2021 + PercentN_SORU_2021 + SORU_Biomass_2021, 
                      data = combined_data)
VarCorr(SORU_model)
anova(SORU_model, SORU_model_noRE) # drop random effects

SORU_model <- SORU_model_noRE

simulation_output <- simulateResiduals(fittedModel = SORU_model)
plot(simulation_output) # fine
testDispersion(simulation_output) # good
testZeroInflation(simulation_output) # good




### POPRC Biomass----
POPRC_model <- lmer(POPRC_Biomass_2023 ~ CageTreatment + POPRC_Biomass_2021 + 
                    (1 | Population), data = combined_data)
#### Check assumptions----
   # Random effects
POPRC_model_noRE <- lm(POPRC_Biomass_2023 ~ CageTreatment + POPRC_Biomass_2021, 
                       data = combined_data)
VarCorr(POPRC_model)

anova(POPRC_model, POPRC_model_noRE)


POPRC_model <- POPRC_model_noRE

simulation_output <- simulateResiduals(fittedModel = POPRC_model)
plot(simulation_output) # good
testDispersion(simulation_output) # good
testZeroInflation(simulation_output) # good





### %N SORU ----
SORUN_model <- lmer(PercentN_SORU_2023 ~ SORU_Biomass_2023 + PlantDiversity_2021 + 
                    PercentN_SORU_2021 + (1 | Population), data = combined_data)

#### Check assumptions----
   # Random effects
SORUN_model_noRE <- lm(PercentN_SORU_2023 ~ SORU_Biomass_2023 + PlantDiversity_2021 + 
                       PercentN_SORU_2021, data = combined_data)

VarCorr(SORUN_model)
anova(SORUN_model, SORUN_model_noRE) # drop random effects

SORUN_model <- SORUN_model_noRE

simulation_output <- simulateResiduals(fittedModel = SORUN_model)
plot(simulation_output) # good
testDispersion(simulation_output) # good
testZeroInflation(simulation_output) # good





### %C SORU ----
SORUC_model <- lmer(PercentC_SORU_2023 ~ SORU_Biomass_2023 + PlantDiversity_2021 + 
                    PercentC_SORU_2021 + (1 | Population), data = combined_data)

VarCorr(SORUC_model) # drop random effects 

SORUC_model_noRE <- lm(PercentC_SORU_2023 ~ SORU_Biomass_2023 + PlantDiversity_2021 + 
                       PercentC_SORU_2021, data = combined_data)
SORUC_model <- SORUC_model_noRE

#### Check assumptions ----

simulation_output <- simulateResiduals(fittedModel = SORUC_model)
plot(simulation_output) # good
testDispersion(simulation_output) # good
testZeroInflation(simulation_output) # good





### %N POPRC ----
POPRCN_model <- lmer(PercentN_POPRC_2023 ~ POPRC_Biomass_2023 + 
                     Overall_mineralization_rate_2021 + PercentN_POPRC_2021 + 
                     (1 | Population), data = combined_data)

VarCorr(POPRCN_model) # drop random effects 

POPRCN_model_noRE <- lm(PercentN_POPRC_2023 ~ POPRC_Biomass_2023 + 
                        Overall_mineralization_rate_2021 + PercentN_POPRC_2021, 
                        data = combined_data)
POPRCN_model <- POPRCN_model_noRE

#### Check assumptions ----
simulation_output <- simulateResiduals(fittedModel = POPRCN_model)
plot(simulation_output) # good
testDispersion(simulation_output) # good
testZeroInflation(simulation_output) # good





### %C POPRC ----
POPRCC_model <- lmer(PercentC_POPRC_2023 ~ POPRC_Biomass_2023 + SORU_Biomass_2023 + 
                     PlantDiversity_2021 + Overall_mineralization_rate_2021 + 
                     POPRC_Biomass_2021 + PercentC_SOIL_2021 + PercentC_POPRC_2021 + 
                     (1 | Population), data = combined_data)

VarCorr(POPRCC_model) # drop random effects 

POPRCC_model_noRE <- lm(PercentC_POPRC_2023 ~ POPRC_Biomass_2023 + SORU_Biomass_2023 + 
                        PlantDiversity_2021 + Overall_mineralization_rate_2021 + 
                        POPRC_Biomass_2021 + PercentC_SOIL_2021 + PercentC_POPRC_2021, 
                        data = combined_data)
POPRCC_model <- POPRCC_model_noRE

#### Check assumptions ----

vif(POPRCC_model) # good

simulation_output <- simulateResiduals(fittedModel = POPRCC_model)
plot(simulation_output) # good
testDispersion(simulation_output) # good
testZeroInflation(simulation_output) # good




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
# removed %N Soil 2021 for standard deviation = 0; failure in bootstrapping
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
  Fixed_Effect = names(coef(SoilN_model)),
  Original_Estimate = coef(SoilN_model),
  Bootstrapped_Median = apply(boot_model$t, 2, median),
  Difference = coef(SoilN_model) - apply(boot_model$t, 2, median)
)

print(ci_comparison) 

print(boot_summary) 


# good; estimates are very stable

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


### Plant Diversity ----
PlantDiversity_model <- lmer(PlantDiversity_2023 ~ PercentC_SOIL_2023 + SORU_Biomass_2023 + 
                              POPRC_Biomass_2023 + PercentN_SORU_2021 + PlantRichness_2021 + 
                              Overall_mineralization_rate_2021 + PercentN_LITTER_2021 + 
                              (1 | Population), data = combined_data)


VarCorr(PlantDiversity_model) # drop random effects

#### Check assumptions ----
PlantDiversity_model_noRE <- lm(PlantDiversity_2023 ~ PercentC_SOIL_2023 + SORU_Biomass_2023 + 
                                 POPRC_Biomass_2023 + PercentN_SORU_2021 + PlantRichness_2021 + 
                                 Overall_mineralization_rate_2021 + PercentN_LITTER_2021, 
                               data = combined_data)


PlantDiversity_model <- PlantDiversity_model_noRE

anova(PlantDiversity_model, PlantDiversity_model_noRE)

vif(PlantDiversity_model_noRE) # good

simulation_output <- simulateResiduals(fittedModel = PlantDiversity_model_noRE)
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

# Find direct paths from treatment categories
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
    Intermediates = paste(sapply(path, function(x) x$Intermediate), collapse = " -> "),
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
    Intermediates = paste(sapply(path, function(x) x$Intermediate), collapse = " -> "),
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

library(DiagrammeR)


#### Predators ----

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

# Function to determine edge color based on Estimate
get_edge_color <- function(estimate) {
  if (estimate < 0) {
    return("red")
  } else {
    return("black")
  }
}

# Apply label mapping to predictors and responses
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

# Optionally, remove rows with unmapped labels
filtered_coefficients_predators_mapped <- filtered_coefficients_predators_mapped %>%
  filter(!is.na(Predictor_Label) & !is.na(Response_Label))

# Calculate absolute estimates for scaling
filtered_coefficients_predators_mapped <- filtered_coefficients_predators_mapped %>%
  mutate(
    Absolute_Estimate = abs(Estimate)
  )

# Define minimum and maximum penwidth
min_penwidth <- 2
max_penwidth <- 7.5

# Calculate scaling factor
max_abs_estimate <- max(filtered_coefficients_predators_mapped$Absolute_Estimate, na.rm = TRUE)

# Scale penwidth based on absolute estimate
filtered_coefficients_predators_mapped <- filtered_coefficients_predators_mapped %>%
  mutate(
    Penwidth = (Absolute_Estimate / max_abs_estimate) * (max_penwidth - min_penwidth) + min_penwidth
  )

# Function to find all downstream nodes and paths
find_downstream_paths <- function(start_nodes, all_paths) {
  # Initialize a vector to keep track of visited nodes
  visited <- c()
  # Initialize a queue with the start nodes
  queue <- start_nodes
  
  # While there are nodes to process
  while (length(queue) > 0) {
    # Pop the first node from the queue
    current_node <- queue[1]
    queue <- queue[-1]
    
    # Mark the current node as visited
    visited <- unique(c(visited, current_node))
    
    # Find all paths starting from the current node
    downstream_paths <- all_paths %>%
      filter(Predictor_Display_Label == current_node)
    
    # Add the response nodes to the queue if they haven't been visited
    new_nodes <- setdiff(downstream_paths$Response_Display_Label, visited)
    queue <- unique(c(queue, new_nodes))
  }
  
  # Return all paths that involve the visited nodes
  all_paths %>%
    filter(Predictor_Display_Label %in% visited | Response_Display_Label %in% visited)
}


# Ensure the columns are correctly named
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

# Extract unique display labels for nodes involved in these paths
unique_labels <- unique(c(filtered_paths$Predictor_Display_Label, filtered_paths$Response_Display_Label))

# Initialize the DOT script
dot_script <- "digraph SEM_PathDiagram { 
  rankdir=LR;
  node [shape=rectangle, style=filled, fillcolor=lightblue];
  
"

# Add nodes with shape customization and background color
for(label in unique_labels){
  if(grepl("baseline", label, ignore.case = TRUE)){
    # Use rectangle shape with dark gray background for baseline nodes
    dot_script <- paste0(dot_script, sprintf("  \"%s\" [shape=rectangle, style=filled, fillcolor=gray, fontsize=10];\n", label))
  } else {
    # Use default rectangle shape with white background for other nodes
    dot_script <- paste0(dot_script, sprintf("  \"%s\" [shape=rectangle, style=filled, fillcolor=white];\n", label))
  }
}

# Add edges with labels, colors, and penwidth
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

# Close the DOT script
dot_script <- paste0(dot_script, "}")

# Render the graph using grViz
grViz(dot_script)




#### Trophic Impact ----

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
      0  # Default value if Vegetation is not present
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
  distinct(Predictor, Response, .keep_all = TRUE)  # Ensure unique paths

# Apply label mapping to predictors and responses
filtered_coefficients_trophic_mapped <- filtered_coefficients_trophic %>%
  mutate(
    Predictor_Label = label_mapping[Predictor],
    Response_Label = label_mapping[Response]
  ) %>%
  rename(
    Predictor_Display_Label = Predictor_Label,
    Response_Display_Label = Response_Label
  )

# Remove paths involving "Vegetation"
filtered_coefficients_trophic_mapped <- filtered_coefficients_trophic_mapped %>%
  filter(Predictor_Display_Label != "Vegetation" & Response_Display_Label != "Vegetation")

# Check for any unmapped labels
unmapped_predictors <- filtered_coefficients_trophic$Predictor[is.na(filtered_coefficients_trophic_mapped$Predictor_Display_Label)]
unmapped_responses <- filtered_coefficients_trophic$Response[is.na(filtered_coefficients_trophic_mapped$Response_Display_Label)]

if(length(unmapped_predictors) > 0){
  warning("Unmapped Predictors found: ", paste(unmapped_predictors, collapse = ", "))
}

if(length(unmapped_responses) > 0){
  warning("Unmapped Responses found: ", paste(unmapped_responses, collapse = ", "))
}

# Calculate absolute estimates for scaling
filtered_coefficients_trophic_mapped <- filtered_coefficients_trophic_mapped %>%
  mutate(
    Absolute_Estimate = abs(Adjusted_Estimate)
  )

# Define minimum and maximum penwidth
min_penwidth <- 2
max_penwidth <- 7.5

# Calculate scaling factor
max_abs_estimate <- max(filtered_coefficients_trophic_mapped$Absolute_Estimate, na.rm = TRUE)

# Scale penwidth based on absolute estimate
filtered_coefficients_trophic_mapped <- filtered_coefficients_trophic_mapped %>%
  mutate(
    Penwidth = (Absolute_Estimate / max_abs_estimate) * (max_penwidth - min_penwidth) + min_penwidth
  )

# Render the path diagram
filtered_paths <- find_downstream_paths(
  start_nodes = c("Predator", "Herbivore"),  # Removed "Vegetation" from start nodes
  all_paths = filtered_coefficients_trophic_mapped
)

# Extract unique display labels for nodes involved in these paths
unique_labels <- unique(c(filtered_paths$Predictor_Display_Label, filtered_paths$Response_Display_Label))

# Initialize the DOT script
dot_script <- "digraph SEM_PathDiagram { 
  rankdir=LR;
  node [shape=rectangle, style=filled, fillcolor=lightblue];
  
"

# Add nodes with shape customization and background color
for(label in unique_labels){
  if(grepl("baseline", label, ignore.case = TRUE)){
    # Use rectangle shape with dark gray background for baseline nodes
    dot_script <- paste0(dot_script, sprintf("  \"%s\" [shape=rectangle, style=filled, fillcolor=gray, fontsize=10];\n", label))
  } else {
    # Use default rectangle shape with white background for other nodes
    dot_script <- paste0(dot_script, sprintf("  \"%s\" [shape=rectangle, style=filled, fillcolor=white];\n", label))
  }
}

# Add edges with labels, colors, and penwidth
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

# Close the DOT script
dot_script <- paste0(dot_script, "}")

# Render the graph using grViz
grViz(dot_script)





#### Generational Effects ----

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
  
  # Initialize the DOT script
  dot_script_generations <- "digraph SEM_PathDiagram_Generations { 
    rankdir=LR;
    node [shape=rectangle, style=filled, fillcolor=lightblue];
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
