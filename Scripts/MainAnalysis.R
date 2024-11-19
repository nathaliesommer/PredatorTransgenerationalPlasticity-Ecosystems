# GenEffects on Ecosystems:

# Objects needed: 
  # N-min.R: N_min_calc
  # SIR.R: SIR_avg
  # Veg.R: diversity_calc; functional_groups_wide

# incorporate the baseline variables that were significant from the previous model
# --> goal is to improve the model fit and increase the individual R.squared

CNdata <- read.csv("Data/CN_Tidy_predators.csv")

CNdata <- CNdata %>%
  dplyr::select(Sample_ID, Population, Site, CageTreatment, PercentN, PercentC, SampleType) %>%
  pivot_wider(names_from = SampleType, values_from = c(PercentN, PercentC))


# Join datasets
combined_data <- full_join(CNdata, N_min_calc, by = "Sample_ID")
combined_data <- full_join(combined_data, SIR_avg, by = "Sample_ID")
combined_data <- full_join(combined_data, diversity_calc, by = "Sample_ID")
combined_data <- full_join(combined_data, functional_groups_wide, by = "Sample_ID")

# Relevel for model interpretation
combined_data$CageTreatment <- as.factor(combined_data$CageTreatment)
combined_data$CageTreatment <- relevel(combined_data$CageTreatment, ref = "Vegetation")

# Packages
library(lme4)
library(DHARMa)
library(car)
library(boot)


# Linear Mixed-Effects Models ----

### SORU Biomass----
SORU_model <- lmer(SORU_Biomass ~ CageTreatment + (1 | Population), data = combined_data)

#### Check assumptions----
   # Random effects
SORU_model_noRE <- lm(SORU_Biomass ~ CageTreatment, data = combined_data)
VarCorr(SORU_model)
anova(SORU_model, SORU_model_noRE) # drop random effects

SORU_model <- SORU_model_noRE

simulation_output <- simulateResiduals(fittedModel = SORU_model)
plot(simulation_output) # good
testDispersion(simulation_output) # good
testZeroInflation(simulation_output) # good




### POPRC Biomass----
POPRC_model <- lmer(POPRC_Biomass ~ CageTreatment + (1 | Population), data = combined_data)

#### Check assumptions----
   # Random effects
POPRC_model_noRE <- lm(POPRC_Biomass ~ CageTreatment, data = combined_data)
VarCorr(POPRC_model)
anova(POPRC_model, POPRC_model_noRE) # drop random effects

POPRC_model <- POPRC_model_noRE

simulation_output <- simulateResiduals(fittedModel = POPRC_model)
plot(simulation_output) # good
testDispersion(simulation_output) # good
testZeroInflation(simulation_output) # good





### %N SORU ----
SORUN_model <- lmer(PercentN_SORU ~ SORU_Biomass + (1 | Population), data = combined_data)

#### Check assumptions----
   # Random effects
SORUN_model_noRE <- lm(PercentN_SORU ~ SORU_Biomass, data = combined_data)
VarCorr(SORUN_model)
anova(SORUN_model, SORUN_model_noRE) # drop random effects

SORUN_model <- SORUN_model_noRE

simulation_output <- simulateResiduals(fittedModel = SORUN_model)
plot(simulation_output) # good
testDispersion(simulation_output) # good
testZeroInflation(simulation_output) # good





### %C SORU ----
SORUC_model <- lmer(PercentC_SORU ~ SORU_Biomass + (1 | Population), data = combined_data) #singularity 
VarCorr(SORUC_model) # drop random effects 

SORUC_model_noRE <- lm(PercentC_SORU ~ SORU_Biomass, data = combined_data)

SORUC_model <- SORUC_model_noRE

#### Check assumptions ----

simulation_output <- simulateResiduals(fittedModel = SORUC_model)
plot(simulation_output) # good
testDispersion(simulation_output) # good
testZeroInflation(simulation_output) # good





### %N POPRC ----
POPRCN_model <- lmer(PercentN_POPRC ~ POPRC_Biomass + (1| Population), data = combined_data) #singularity 
VarCorr(POPRCN_model) # drop random effects 

POPRCN_model_noRE <- lm(PercentN_POPRC ~ POPRC_Biomass, data = combined_data)

POPRCN_model <- POPRCN_model_noRE

#### Check assumptions ----
simulation_output <- simulateResiduals(fittedModel = POPRCN_model)
plot(simulation_output) # good
testDispersion(simulation_output) # good
testZeroInflation(simulation_output) # good





### %C POPRC ----
POPRCC_model <- lmer(PercentC_POPRC ~ POPRC_Biomass + SORU_Biomass + (1| Population), data = combined_data) #singularity 
VarCorr(POPRCC_model) # drop random effects 

POPRCC_model_noRE <- lm(PercentC_POPRC ~ POPRC_Biomass + SORU_Biomass, data = combined_data)

POPRCC_model <- POPRCC_model_noRE

#### Check assumptions ----

vif(POPRCC_model) # good

simulation_output <- simulateResiduals(fittedModel = POPRCC_model)
plot(simulation_output) # good
testDispersion(simulation_output) # good
testZeroInflation(simulation_output) # good




### %C Litter ----
LitterC_model <- lmer(PercentC_LITTER ~ PercentC_SORU + PercentC_POPRC + (1 | Population), data = combined_data) #singularity
VarCorr(LitterC_model) # dtop random effects

LitterC_model_noRE <- lm(PercentC_LITTER ~ PercentC_SORU + PercentC_POPRC, data = combined_data)

LitterC_model <- LitterC_model_noRE

#### Check assumptions ----

vif(LitterC_model) # good

simulation_output <- simulateResiduals(fittedModel = LitterC_model)
plot(simulation_output) # good
testDispersion(simulation_output) # good
testZeroInflation(simulation_output) # good





### %N Litter ----
LitterN_model <- lmer(PercentN_LITTER ~ PercentN_SORU + PercentN_POPRC + (1 | Population), data = combined_data)

#### Check assumptions ----
# Random effects
LitterN_model_noRE <- lm(PercentN_LITTER ~ PercentN_SORU + PercentN_POPRC, data = combined_data)
VarCorr(LitterN_model)
anova(LitterN_model, LitterN_model_noRE) # drop random effects

LitterN_model <- LitterN_model_noRE

vif(LitterN_model) # good

simulation_output <- simulateResiduals(fittedModel = LitterN_model)
plot(simulation_output) # good
testDispersion(simulation_output) # good
testZeroInflation(simulation_output) # good



### %C Soil ----
SoilC_model <- lmer(PercentC_SOIL ~ CageTreatment + SORU_Biomass + CO2CperHourperg + (1 | Population), data = combined_data) #singularity
VarCorr(SoilC_model) # drop random effects

SoilC_model_noRE <- lm(PercentC_SOIL ~ CageTreatment + SORU_Biomass + CO2CperHourperg, data = combined_data)

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
  mod <- lm(PercentC_SOIL ~ 
              CageTreatment +
              SORU_Biomass +
              CO2CperHourperg,
            data = d)
  return(coef(mod))  # return coefficients for lm objects
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
SoilN_model <- lmer(PercentN_SOIL ~ CageTreatment + PercentN_LITTER  + (1 | Population), data = combined_data)

#### Check assumptions ----
   # Random effects
SoilN_model_noRE <- lm(PercentN_SOIL ~ CageTreatment + PercentN_LITTER, data = combined_data)
VarCorr(SoilN_model)
anova(SoilN_model, SoilN_model_noRE) # drop random effects

SoilN_model <- SoilN_model_noRE

vif(SoilN_model) # good

simulation_output <- simulateResiduals(fittedModel = SoilN_model)
plot(simulation_output)
testDispersion(simulation_output)
testZeroInflation(simulation_output) # assumptions not met

# --> non-parametric bootstrapping

boot_fun <- function(data, indices) {
  d <- data[indices, ]  # resample
  mod <- lm(PercentN_SOIL ~ 
              CageTreatment +
              PercentN_LITTER +
              PercentC_SOIL,
            data = d)
  return(coef(mod))  # return coefficients for lm objects
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

# good; estimates are stable

### SIR ----
SIR_model <- lmer(CO2CperHourperg ~ PercentN_SOIL + PercentN_POPRC + (1 | Population), data = combined_data)

#### Check assumptions ----
# Random effects
SIR_model_noRE <- lm(CO2CperHourperg ~ PercentN_SOIL  + PercentN_POPRC, data = combined_data)
VarCorr(SIR_model)
anova(SIR_model, SIR_model_noRE) # drop random effects

SIR_model <- SIR_model_noRE

simulation_output <- simulateResiduals(fittedModel = SIR_model)
plot(simulation_output) # assumptions not met
testDispersion(simulation_output)
testZeroInflation(simulation_output)

# --> non-parametric bootstrapping

boot_fun <- function(data, indices) {
  d <- data[indices, ]  # resample
  tryCatch({
    mod <- lm(CO2CperHourperg ~ PercentN_SOIL + PercentN_POPRC, data = d)
    return(coef(mod))  # return coefficients for lm objects
  }, error = function(e) {
    warning("Error in bootstrapping iteration: ", e$message)
    return(rep(NA, length(coef(SIR_model))))  # return NA for all coefficients if an error occurs
  })
}

# Re-run bootstrapping with debugging
set.seed(1231)
boot_model <- boot(combined_data, boot_fun, R = 1000)

# Check for NA values in bootstrapping results
na_indices <- which(is.na(boot_model$t[, which(names(coef(SIR_model)) == "PercentN_SOIL")]))
if (length(na_indices) > 0) {
  warning("NA values found in bootstrapping results for PercentN_SOIL at iterations: ", paste(na_indices, collapse = ", "))
}

# Recalculate bootstrapped medians
boot_medians <- apply(boot_model$t, 2, median, na.rm = TRUE)

# Compare bootstrapped estimates with original fixed effect estimates
boot_summary <- data.frame(
  Fixed_Effect = names(coef(SIR_model)),
  Original_Estimate = coef(SIR_model),
  Bootstrapped_Median = boot_medians,
  Difference = coef(SIR_model) - boot_medians
)

print(boot_summary)

# extract original confidence intervals for lm model using `coef()`
original_ci <- confint(SIR_model, parm = names(coef(SIR_model)), level = 0.95)
original_ci <- as.data.frame(original_ci)
colnames(original_ci) <- c("Original_Lower", "Original_Upper")

# initialize data frame
boot_ci_table <- data.frame(
  Fixed_Effect = names(coef(SIR_model)),
  Boot_Lower = NA,
  Boot_Upper = NA
)

# iterate over fixed effects to get bootstrapped CIs
for (i in seq_along(coef(SIR_model))) {
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

print(ci_comparison) 



### Plant Diversity ----
PlantDiversity_model <- lmer(PlantDiversity ~ PercentC_SOIL  + SORU_Biomass + POPRC_Biomass + (1 | Population), data = combined_data)

#### Check assumptions ----
   # Random effects
PlantDiversity_model_noRE <- lm(PlantDiversity ~ PercentC_SOIL  + SORU_Biomass + POPRC_Biomass, data = combined_data)
VarCorr(PlantDiversity_model)
anova(PlantDiversity_model, PlantDiversity_model_noRE) # drop random effects


PlantDiversity_model <- PlantDiversity_model_noRE

vif(PlantDiversity_model) # good
simulation_output <- simulateResiduals(fittedModel = PlantDiversity_model)
plot(simulation_output)
testDispersion(simulation_output)
testZeroInflation(simulation_output) # assumptions not met

# --> non-parametric bootstrapping

boot_fun <- function(data, indices) {
  d <- data[indices, ]  # resample
  mod <- lm(PlantDiversity ~ 
              PercentC_SOIL +
              SORU_Biomass +
              POPRC_Biomass,
            data = d)
  return(coef(mod))  # return coefficients for lm objects
}


set.seed(1231) # set seed for reproducibility

boot_model <- boot(combined_data, boot_fun, R = 1000) # perform bootstrapping with 1000 resamples


# extract original confidence intervals for lm model using `coef()`
original_ci <- confint(PlantDiversity_model, parm = names(coef(PlantDiversity_model)), level = 0.95)
original_ci <- as.data.frame(original_ci)
colnames(original_ci) <- c("Original_Lower", "Original_Upper")

# initialize data frame
boot_ci_table <- data.frame(
  Fixed_Effect = names(coef(PlantDiversity_model)),
  Boot_Lower = NA,
  Boot_Upper = NA
)

# iterate over fixed effects to get bootstrapped CIs
for (i in seq_along(coef(PlantDiversity_model))) {
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
  Fixed_Effect = names(coef(PlantDiversity_model)),
  Original_Estimate = coef(PlantDiversity_model),
  Bootstrapped_Median = apply(boot_model$t, 2, median),
  Difference = coef(PlantDiversity_model) - apply(boot_model$t, 2, median)
)

print(ci_comparison) 

print(boot_summary) 

# good; estimates are stable


### N-min ----
Nmin_model <- lmer(Overall_rate ~ PercentN_SOIL + CO2CperHourperg + CageTreatment + (1 | Population) , data = combined_data)

#### Check assumptions ----
   # Random effects
Nmin_model_noRE <- lm(Overall_rate ~ PercentN_SOIL + CO2CperHourperg + CageTreatment, data = combined_data)
VarCorr(Nmin_model)
anova(Nmin_model, Nmin_model_noRE) # drop random effects

Nmin_model <- Nmin_model_noRE

vif(Nmin_model) # good
simulation_output <- simulateResiduals(fittedModel = Nmin_model)
plot(simulation_output)
testDispersion(simulation_output)
testZeroInflation(simulation_output) # assumptions not met

# --> non-parametric bootstrapping

boot_fun <- function(data, indices) {
  d <- data[indices, ]  # resample
  tryCatch({
    mod <- lm(Overall_rate ~ PercentN_SOIL + CO2CperHourperg + CageTreatment, data = d)
    return(coef(mod))  # return coefficients for lm objects
  }, error = function(e) {
    warning("Error in bootstrapping iteration: ", e$message)
    return(rep(NA, length(coef(Nmin_model))))  # return NA for all coefficients if an error occurs
  })
}

# Re-run bootstrapping with debugging
set.seed(1231)
boot_model <- boot(combined_data, boot_fun, R = 1000)

# Check for NA values in bootstrapping results
na_indices <- which(is.na(boot_model$t[, which(names(coef(Nmin_model)) == "PercentN_SOIL")]))
if (length(na_indices) > 0) {
  warning("NA values found in bootstrapping results for PercentN_SOIL at iterations: ", paste(na_indices, collapse = ", "))
}

# Recalculate bootstrapped medians
boot_medians <- apply(boot_model$t, 2, median, na.rm = TRUE)

# Compare bootstrapped estimates with original fixed effect estimates
boot_summary <- data.frame(
  Fixed_Effect = names(coef(Nmin_model)),
  Original_Estimate = coef(Nmin_model),
  Bootstrapped_Median = boot_medians,
  Difference = coef(Nmin_model) - boot_medians
)

print(boot_summary)

# extract original confidence intervals for lm model using `coef()`
original_ci <- confint(Nmin_model, parm = names(coef(Nmin_model)), level = 0.95)
original_ci <- as.data.frame(original_ci)
colnames(original_ci) <- c("Original_Lower", "Original_Upper")

# initialize data frame
boot_ci_table <- data.frame(
  Fixed_Effect = names(coef(Nmin_model)),
  Boot_Lower = NA,
  Boot_Upper = NA
)

# iterate over fixed effects to get bootstrapped CIs
for (i in seq_along(coef(Nmin_model))) {
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

print(ci_comparison) 


# SEM ----
library(piecewiseSEM)

# Define the SEM model
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
  PlantDiversity_model,
  Nmin_model,
  data = combined_data)

summary(sem_model)

options(max.print = 10000)
sink("SEM_Summary_Final.txt") # saved as SEM_Summary_Initial before incorporating tests of directed separation 
print(summary(sem_model))
sink()

summary_sem_model <- summary(sem_model)

### Extract standardized coefficients for indirect path assessment ----
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


# Re-define the mapping for labels
label_mapping <- c(
  "SORU_Biomass" = "Goldenrod biomass",
  "POPRC_Biomass" = "Grass biomass",
  "PercentN_SORU" = "Goldenrod %N",
  "PercentC_SORU" = "Goldenrod %C",
  "PercentN_POPRC" = "Grass %N",
  "PercentC_POPRC" = "Grass %C",
  "PercentC_SOIL" = "Soil %C",
  "PercentN_SOIL" = "Soil %N",
  "CO2CperHourperg" = "SIR",
  "PlantDiversity" = "Plant diversity",
  "PercentN_LITTER" = "Litter %N",
  "PercentC_LITTER" = "Litter %C",
  "Overall_rate" = "Nitrogen mineralization",
  "CageTreatment = Vegetation" = "Vegetation",
  "CageTreatment = Herbivore" = "Herbivore",
  "Predator" = "Predator"
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
  "SORU_Biomass" = "Goldenrod biomass",
  "POPRC_Biomass" = "Grass biomass",
  "PercentN_SORU" = "Goldenrod %N",
  "PercentC_SORU" = "Goldenrod %C",
  "PercentN_POPRC" = "Grass %N",
  "PercentC_POPRC" = "Grass %C",
  "PercentC_SOIL" = "Soil %C",
  "PercentN_SOIL" = "Soil %N",
  "CO2CperHourperg" = "SIR",
  "PlantDiversity" = "Plant diversity",
  "PercentN_LITTER" = "Litter %N",
  "PercentC_LITTER" = "Litter %C",
  "Overall_rate" = "Nitrogen mineralization",
  "CageTreatment = Vegetation" = "Vegetation",
  "Predator F1" = "Predator F1",
  "Predator F2" = "Predator F2"
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

# Render the path diagram
filtered_paths_generations <- find_downstream_paths(
  start_nodes = c("Predator F1", "Predator F2"),  # Start nodes for the generational effects
  all_paths = filtered_coefficients_generations_mapped
)

# Check if the data frame is not empty
if (nrow(filtered_paths_generations) > 0) {
  # Extract unique display labels for nodes involved in these paths
  unique_labels_generations <- unique(c(filtered_paths_generations$Predictor_Display_Label, filtered_paths_generations$Response_Display_Label))

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

  