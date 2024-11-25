# G1 and G2 Plasticity Comparison
# Code by Baker and Sommer



# Behavior ----

# Packages 
library(dplyr)
library(lme4)
library(lubridate)
library(ggplot2)
library(tidyr)
library(ggforce)
library(emmeans)
library(ks)
library(gghalves)

### Data Processing ----
Behavior_2023 <- read.csv("Data/Behavior_Raw.csv")
# Raw X, Y data measured in 4cm x 4cm grid coordinates. Converting to cm
Behavior_2023$Y <- Behavior_2023$Y*4
Behavior_2023$X <- Behavior_2023$X*4
Behavior_2023$Distance <- Behavior_2023$Distance*4

# Adjusted start time
Behavior_2023$Time <- hms(Behavior_2023$Time)
behavior_data <- subset(Behavior_2023,Time>(hms("8:40:00")))

behavior_temp <- read.csv("Data/Behavior_Temps_2023_F1_F2_Raw.csv")

behavior_temp <- behavior_temp %>%
  mutate(DateTime = mdy_hm(Date.Time..EDT.))

behavior_data <- behavior_data %>%
  mutate(DateTime = mdy(Date) + hours(hour(Time)) + minutes(minute(Time)))

# Function to find the nearest temperature
find_nearest_temp <- function(obs_time, temp_data) {
  temp_data %>%
    filter(abs(difftime(DateTime, obs_time, units = "mins")) == min(abs(difftime(DateTime, obs_time, units = "mins")))) %>%
    slice(1) %>%
    pull(Ch..1...Temperature.....C.)
}

# Join temperature to behavior_data
behavior_data <- behavior_data %>%
  rowwise() %>%
  mutate(Nearest_Temperature = find_nearest_temp(DateTime, behavior_temp))

# Add a unique id
behavior_data <- behavior_data %>%
  mutate(id = factor(paste(Population, Predation.Treatment, Generation, Day, Terraria, sep = "_")))



# Update any Generation references
behavior_data <- behavior_data %>%
  mutate(Generation = recode(Generation,
                             "F1" = "G1",
                             "F2" = "G2"))

# Create the combined treatment column
behavior_data <- behavior_data %>%
  mutate(Treatment = paste(Generation, Predation.Treatment, sep = "_")) %>%
  mutate(Treatment = factor(Treatment, 
                            levels = c("G1_Herbivore", "G2_Herbivore", 
                                       "G1_Predator", "G2_Predator")))

# Create 'type' column sensu Baker 2022 data

behavior_data <- behavior_data %>%
  mutate(Type = case_when(
    Population %in% c("UP", "MC", "FN", "HF") ~ "Behavior",
    Population %in% c("DC", "SP", "SC", "YF") ~ "Physiology"))

### Analysis ----

# Check population effect
height_model <- lmer(Y ~ Treatment*Type + (1|Population/id), 
                     data = behavior_data)

summary(height_model)

VarCorr(height_model)

# population and individual does explain a significant chunk of variation
# however, given lack of population effect in the ecosystems model, 
# we might not want to model this explicitly.
# --> use population-centered heights for kernel density and BA estimates 
# interpretation is about how F1 and F2 differ in their canopy height, 
# regardless of the absolute height tendencies of different populations

# we also need to use a boundary-corrected kernel density approach;
# limits estimates >0 (i.e., above ground)


# centered population heights
behavior_data <- behavior_data %>%
  group_by(Population) %>%
  mutate(Y_centered = scale(Y, scale = FALSE),
         Y_centered = Y_centered - min(Y_centered)) %>%  # make all values positive
  ungroup()

# data validation
print(range(behavior_data$Y_centered, na.rm = TRUE)) # range of Y_centered
print(sum(!is.finite(behavior_data$Y_centered))) # check for infinite values
print(sum(is.na(behavior_data$Y_centered))) # check for NA

# diagnostics
print(summary(behavior_data$Y_centered))
print(summary(behavior_data$Y_centered[behavior_data$Generation == "G1"]))
print(summary(behavior_data$Y_centered[behavior_data$Generation == "G2"]))

# First, let's check the input data
print("Input data summary:")
print(summary(behavior_data$Y_centered))
print(table(behavior_data$Treatment))

# Let's modify the density estimation approach
density_estimate_ranges <- behavior_data %>%
  group_by(id) %>%
  filter(n() >= 10) %>%  # Ensure enough points for density estimation
  do({
    # Get the data for this group
    y_values <- .$Y_centered
    
    tryCatch({
      # Use base R density estimation which is more robust
      density_estimate <- density(y_values, 
                                  from = 0,  # Force minimum at 0
                                  to = max(y_values),
                                  n = 512,
                                  kernel = "gaussian")
      
      # Calculate cumulative density
      cum_density <- cumsum(density_estimate$y) / sum(density_estimate$y)
      
      # Find the ranges for 50% and 95% density
      core_idx <- which(cum_density <= 0.50)
      broad_idx <- which(cum_density <= 0.95)
      
      if(length(core_idx) > 0 && length(broad_idx) > 0) {
        data.frame(
          core_low = min(density_estimate$x[core_idx]),
          core_high = max(density_estimate$x[core_idx]),
          broad_low = min(density_estimate$x[broad_idx]),
          broad_high = max(density_estimate$x[broad_idx]),
          Treatment = .$Treatment[1]
        )
      } else {
        # Return NA if we couldn't find the ranges
        data.frame(
          core_low = NA_real_,
          core_high = NA_real_,
          broad_low = NA_real_,
          broad_high = NA_real_,
          Treatment = .$Treatment[1]
        )
      }
    }, error = function(e) {
      # Return NA if there was an error
      data.frame(
        core_low = NA_real_,
        core_high = NA_real_,
        broad_low = NA_real_,
        broad_high = NA_real_,
        Treatment = .$Treatment[1]
      )
    })
  }) %>%
  filter(!is.na(core_low)) %>%  # Remove any failed estimations
  group_by(Treatment) %>%
  summarize(
    avg_core_low = mean(core_low),
    avg_core_high = mean(core_high),
    avg_broad_low = mean(broad_low),
    avg_broad_high = mean(broad_high)
  ) %>%
  mutate(Treatment = factor(Treatment,
                            levels = c("G1_Herbivore", "G2_Herbivore", 
                                       "G1_Predator", "G2_Predator")))

# Check the results
print("Updated density estimates:")
print(density_estimate_ranges)

# calculate 1D kernel density estimate for each id
density_estimates <- behavior_data %>%
  group_by(id) %>%
  do({
    density_estimate <- density(.$Y, bw = "nrd0")
    cum_density <- cumsum(density_estimate$y) / sum(density_estimate$y)
    
    # determine 50% and 95% isopleths
    core_threshold <- 0.50
    broad_threshold <- 0.95
    
    core_level <- density_estimate$x[which(cum_density >= core_threshold)[1]]
    broad_level <- density_estimate$x[which(cum_density >= broad_threshold)[1]]
    
    data.frame(core_level = core_level, broad_level = broad_level)
  })

# function to calculate isopleths for each generation
calculate_isopleths <- function(y) {
  density_estimate <- density(y[y > 0], bw = 5)  # restrict to height > 0
  cum_density <- cumsum(density_estimate$y) / sum(density_estimate$y)
  
  core_threshold <- 0.50
  broad_threshold <- 0.95
  
  core_range <- range(density_estimate$x[cum_density <= core_threshold])
  broad_range <- range(density_estimate$x[cum_density <= broad_threshold])
  
  return(list(core_range = core_range, broad_range = broad_range))
}

# calculate isopleths for each generation
isopleths <- behavior_data %>%
  group_by(Generation, Predation.Treatment) %>%
  summarize(core_range = list(calculate_isopleths(Y)$core_range),
            broad_range = list(calculate_isopleths(Y)$broad_range))




#### Bhattacharyya's affinity ----

# Function to calculate Bhattacharyya's affinity
calculate_bhattacharyya_affinity_isopleth <- function(density1, density2, range) {
  # Restrict the density estimates to the specified range
  indices1 <- which(density1$eval.points >= range[1] & density1$eval.points <= range[2])
  indices2 <- which(density2$eval.points >= range[1] & density2$eval.points <= range[2])
  
  # Ensure the restricted grids are the same
  if (length(indices1) != length(indices2) || !all(density1$eval.points[indices1] == density2$eval.points[indices2])) {
    stop("Density grids do not match within the specified range.")
  }
  
  # Normalize the density estimates within the specified range
  density1_norm <- density1$estimate[indices1] / sum(density1$estimate[indices1])
  density2_norm <- density2$estimate[indices2] / sum(density2$estimate[indices2])
  
  # Calculate Bhattacharyya's affinity
  sqrt_density1 <- sqrt(density1_norm)
  sqrt_density2 <- sqrt(density2_norm)
  affinity <- sum(sqrt_density1 * sqrt_density2) * prod(diff(density1$eval.points[indices1]))
  
  return(affinity)
}

# Calculate 1D kernel density estimates for Y data
y_min <- min(behavior_data$Y_centered, na.rm = TRUE)
y_max <- max(behavior_data$Y_centered, na.rm = TRUE)
common_grid <- seq(y_min, y_max, length.out = 100)


# calculate densities for each generation
density_G1P <- stats::density(behavior_data$Y_centered[behavior_data$Treatment == "G1_Predator"],
                              from = y_min, to = y_max, n = 100)
density_G2P <- stats::density(behavior_data$Y_centered[behavior_data$Treatment == "G2_Predator"],
                              from = y_min, to = y_max, n = 100)
density_G1H <- stats::density(behavior_data$Y_centered[behavior_data$Treatment == "G1_Herbivore"],
                              from = y_min, to = y_max, n = 100)
density_G2H <- stats::density(behavior_data$Y_centered[behavior_data$Treatment == "G2_Herbivore"],
                              from = y_min, to = y_max, n = 100)

# Convert to kde format
density_G1P <- list(eval.points = density_G1P$x, estimate = density_G1P$y)
density_G2P <- list(eval.points = density_G2P$x, estimate = density_G2P$y)
density_G1H <- list(eval.points = density_G1H$x, estimate = density_G1H$y)
density_G2H <- list(eval.points = density_G2H$x, estimate = density_G2H$y)

##### 50% isopleth -----
# Calculate ranges for all treatments
cum_density_G1P <- cumsum(density_G1P$estimate) / sum(density_G1P$estimate)
core_range_G1P <- range(density_G1P$eval.points[cum_density_G1P <= 0.5])

cum_density_G2P <- cumsum(density_G2P$estimate) / sum(density_G2P$estimate)
core_range_G2P <- range(density_G2P$eval.points[cum_density_G2P <= 0.5])

cum_density_G1H <- cumsum(density_G1H$estimate) / sum(density_G1H$estimate)
core_range_G1H <- range(density_G1H$eval.points[cum_density_G1H <= 0.5])

cum_density_G2H <- cumsum(density_G2H$estimate) / sum(density_G2H$estimate)
core_range_G2H <- range(density_G2H$eval.points[cum_density_G2H <= 0.5])

##### 95% isopleth -----
# Calculate broad ranges for all treatments
broad_range_G1P <- range(density_G1P$eval.points[cum_density_G1P <= 0.95])
broad_range_G2P <- range(density_G2P$eval.points[cum_density_G2P <= 0.95])
broad_range_G1H <- range(density_G1H$eval.points[cum_density_G1H <= 0.95])
broad_range_G2H <- range(density_G2H$eval.points[cum_density_G2H <= 0.95])

# Calculate BAs for all comparisons
# G1 Predator vs G2 Predator
core_range_P1P2 <- c(max(core_range_G1P[1], core_range_G2P[1]), 
                     min(core_range_G1P[2], core_range_G2P[2]))
broad_range_P1P2 <- c(max(broad_range_G1P[1], broad_range_G2P[1]), 
                      min(broad_range_G1P[2], broad_range_G2P[2]))

BA_P1P2_core <- calculate_bhattacharyya_affinity_isopleth(density_G1P, density_G2P, core_range_P1P2)
BA_P1P2_broad <- calculate_bhattacharyya_affinity_isopleth(density_G1P, density_G2P, broad_range_P1P2)

# G1 Predator vs G1 Herbivore
core_range_P1H1 <- c(max(core_range_G1P[1], core_range_G1H[1]), 
                     min(core_range_G1P[2], core_range_G1H[2]))
broad_range_P1H1 <- c(max(broad_range_G1P[1], broad_range_G1H[1]), 
                      min(broad_range_G1P[2], broad_range_G1H[2]))

BA_P1H1_core <- calculate_bhattacharyya_affinity_isopleth(density_G1P, density_G1H, core_range_P1H1)
BA_P1H1_broad <- calculate_bhattacharyya_affinity_isopleth(density_G1P, density_G1H, broad_range_P1H1)


#### Consolidated stats  ----
# G1 Predator vs G2 Predator
cat("G1 Predator vs G2 Predator:\n")
cat("Core BA:", BA_P1P2_core, "\n")
cat("Broad BA:", BA_P1P2_broad, "\n\n")

# G1 Predator vs G1 Herbivore
cat("G1 Predator vs G1 Herbivore:\n")
cat("Core BA:", BA_P1H1_core, "\n")
cat("Broad BA:", BA_P1H1_broad, "\n\n")



#### Figures ----

# restructure for plotting
density_long <- density_estimate_ranges %>%
  pivot_longer(cols = c(avg_core_low, avg_core_high, avg_broad_low, avg_broad_high),
               names_to = "range_type",
               values_to = "height") %>%
  mutate(type = ifelse(grepl("core", range_type), "core", "broad"))

average_ranges <- density_estimate_ranges


custom_colors <- c(
  "#FED789FF",  # G1 Herbivore (yellow) 50%
  "#FEE4B1FF",  # G1 Herbivore (light yellow) 95% - lighter variant
  "#72874EFF",  # G2 Herbivore (green) 50%
  "#9DAC82FF",  # G2 Herbivore (light green) 95% - lighter variant
  "#70646EFF",  # G1 Predator (purple) 50%
  "#9B939AFF",  # G1 Predator (light purple) 95% - lighter variant
  "#476F84FF",  # G2 Predator (blue) 50%
  "#7B99AAFF"   # G2 Predator (light blue) 95% - lighter variant
)

# Range plot for ALL data
ggplot(average_ranges) +
  # Add half violin plot to the right
  geom_half_violin(data = behavior_data,
                   aes(x = Treatment, y = Y_centered, fill = Treatment),
                   side = "r", alpha = 0.3, position = position_nudge(x = 0.07)) +
  # G1 Herbivore range
  geom_ellipse(data = subset(average_ranges, Treatment == "G1_Herbivore"),
               aes(x0 = 1 - 0.15, 
                   y0 = (avg_broad_low + avg_broad_high) / 2,
                   a = 0.15, b = (avg_broad_high - avg_broad_low) / 2, 
                   angle = 0), 
               fill = custom_colors[2], color = "black", alpha = 0.5) +
  geom_ellipse(data = subset(average_ranges, Treatment == "G1_Herbivore"),
               aes(x0 = 1 - 0.15, 
                   y0 = (avg_core_low + avg_core_high) / 2,
                   a = 0.15, b = (avg_core_high - avg_core_low) / 2, 
                   angle = 0), 
               fill = custom_colors[1], color = "black", alpha = 0.5) +
  # G2 Herbivore range
  geom_ellipse(data = subset(average_ranges, Treatment == "G2_Herbivore"),
               aes(x0 = 2 - 0.15, 
                   y0 = (avg_broad_low + avg_broad_high) / 2,
                   a = 0.15, b = (avg_broad_high - avg_broad_low) / 2, 
                   angle = 0), 
               fill = custom_colors[4], color = "black", alpha = 0.5) +
  geom_ellipse(data = subset(average_ranges, Treatment == "G2_Herbivore"),
               aes(x0 = 2 - 0.15, 
                   y0 = (avg_core_low + avg_core_high) / 2,
                   a = 0.15, b = (avg_core_high - avg_core_low) / 2, 
                   angle = 0), 
               fill = custom_colors[3], color = "black", alpha = 0.5) +
  # G1 Predator range
  geom_ellipse(data = subset(average_ranges, Treatment == "G1_Predator"),
               aes(x0 = 3 - 0.15, 
                   y0 = (avg_broad_low + avg_broad_high) / 2,
                   a = 0.15, b = (avg_broad_high - avg_broad_low) / 2, 
                   angle = 0), 
               fill = custom_colors[6], color = "black", alpha = 0.5) +
  geom_ellipse(data = subset(average_ranges, Treatment == "G1_Predator"),
               aes(x0 = 3 - 0.15, 
                   y0 = (avg_core_low + avg_core_high) / 2,
                   a = 0.15, b = (avg_core_high - avg_core_low) / 2, 
                   angle = 0), 
               fill = custom_colors[5], color = "black", alpha = 0.5) +
  # G2 Predator range
  geom_ellipse(data = subset(average_ranges, Treatment == "G2_Predator"),
               aes(x0 = 4 - 0.15, 
                   y0 = (avg_broad_low + avg_broad_high) / 2,
                   a = 0.15, b = (avg_broad_high - avg_broad_low) / 2, 
                   angle = 0), 
               fill = custom_colors[8], color = "black", alpha = 0.5) +
  geom_ellipse(data = subset(average_ranges, Treatment == "G2_Predator"),
               aes(x0 = 4 - 0.15, 
                   y0 = (avg_core_low + avg_core_high) / 2,
                   a = 0.15, b = (avg_core_high - avg_core_low) / 2, 
                   angle = 0), 
               fill = custom_colors[7], color = "black", alpha = 0.5) +
  # Add mean points and text
  geom_point(data = behavior_data %>% 
               group_by(Treatment) %>% 
               summarize(mean_y = mean(Y_centered)),
             aes(x = as.numeric(Treatment) - 0.15, y = mean_y),
             size = 3, color = "black") +
  geom_text(data = behavior_data %>% 
              group_by(Treatment) %>% 
              summarize(mean_y = mean(Y_centered)),
            aes(x = as.numeric(Treatment) - 0.15, y = mean_y, 
                label = round(mean_y, 1)),
            vjust = -1.5, size = 4) +
  # Aesthetics
  scale_fill_manual(values = c(
    custom_colors[1],  # G1 Herbivore
    custom_colors[3],  # G2 Herbivore
    custom_colors[5],  # G1 Predator
    custom_colors[7]   # G2 Predator
  )) +
  labs(title = "Core (50%) and Broad Ranges (95%) by Treatment", 
       x = "Treatment", 
       y = "Population-Centered Canopy Height") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        text = element_text(size = 12),
        legend.position = "none")

# For the simple violin plot
ggplot(behavior_data, aes(x = Treatment, y = Y_centered, fill = Treatment)) +
  geom_violin(alpha = 0.5) +
  geom_point(data = behavior_data %>% 
               group_by(Treatment) %>% 
               summarize(mean_y = mean(Y_centered)),
             aes(y = mean_y),
             size = 3, color = "black") +
  geom_text(data = behavior_data %>% 
              group_by(Treatment) %>% 
              summarize(mean_y = mean(Y_centered)),
            aes(y = mean_y, label = round(mean_y, 1)),
            vjust = -1.5, size = 4) +
  scale_fill_manual(values = c(
    custom_colors[1],  # G1 Herbivore
    custom_colors[3],  # G2 Herbivore
    custom_colors[5],  # G1 Predator
    custom_colors[7]   # G2 Predator
  )) +
  labs(title = "Height Distribution by Treatment",
       subtitle = "Population-centered mean height",
       x = "Treatment",
       y = "Population-Centered Canopy Height") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 12),
        legend.position = "none")

# Temperature plot
ggplot(behavior_data, aes(x = Nearest_Temperature, y = Y, color = Treatment)) +
  geom_jitter(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = c(
    custom_colors[1],  # G1 Herbivore
    custom_colors[3],  # G2 Herbivore
    custom_colors[5],  # G1 Predator
    custom_colors[7]   # G2 Predator
  )) +
  labs(title = "Relationship between Height and Temperature by Generation",
       x = "Temperature (°C)", 
       y = "Height (cm)") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        text = element_text(size = 12),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom") +
  facet_wrap(~Population)






#### Behavior Type Analysis ----


print("Number of observations by Type:")
print(table(behavior_data$Type))

# Check data structure for behavior type
behavior_subset <- behavior_data %>%
  filter(Type == "Behavior")
print("Number of observations by Treatment for behavior type:")
print(table(behavior_subset$Treatment))

# Check how many observations we have per ID
id_counts <- behavior_subset %>%
  group_by(id) %>%
  summarise(n = n())
print("Observations per ID for behavior type:")
print(summary(id_counts$n))

# Calculate density estimates for behavior type
density_estimate_ranges_behavior <- behavior_data %>%
  filter(Type == "Behavior") %>%
  group_by(id) %>%
  filter(n() >= 10) %>%  # Ensure enough points for density estimation
  do({
    # Get the data for this group
    y_values <- .$Y_centered
    
    tryCatch({
      # Use base R density estimation
      density_estimate <- density(y_values, 
                                  from = 0,  # Force minimum at 0
                                  to = max(y_values),
                                  n = 512,
                                  kernel = "gaussian")
      
      # Calculate cumulative density
      cum_density <- cumsum(density_estimate$y) / sum(density_estimate$y)
      
      # Find the ranges for 50% and 95% density
      core_idx <- which(cum_density <= 0.50)
      broad_idx <- which(cum_density <= 0.95)
      
      if(length(core_idx) > 0 && length(broad_idx) > 0) {
        data.frame(
          core_low = min(density_estimate$x[core_idx]),
          core_high = max(density_estimate$x[core_idx]),
          broad_low = min(density_estimate$x[broad_idx]),
          broad_high = max(density_estimate$x[broad_idx]),
          Treatment = .$Treatment[1]
        )
      } else {
        data.frame(
          core_low = NA_real_,
          core_high = NA_real_,
          broad_low = NA_real_,
          broad_high = NA_real_,
          Treatment = .$Treatment[1]
        )
      }
    }, error = function(e) {
      data.frame(
        core_low = NA_real_,
        core_high = NA_real_,
        broad_low = NA_real_,
        broad_high = NA_real_,
        Treatment = .$Treatment[1]
      )
    })
  }) %>%
  filter(!is.na(core_low)) %>%  # Remove any failed estimations
  group_by(Treatment) %>%
  summarize(
    avg_core_low = mean(core_low),
    avg_core_high = mean(core_high),
    avg_broad_low = mean(broad_low),
    avg_broad_high = mean(broad_high)
  ) %>%
  mutate(Treatment = factor(Treatment,
                            levels = c("G1_Herbivore", "G2_Herbivore", 
                                       "G1_Predator", "G2_Predator")))

# Check the results
print("Behavior type density estimates:")
print(density_estimate_ranges_behavior)


# Run ANOVA on behavior type data only
behavior_aov <- aov(Y_centered ~ Treatment, 
                    data = subset(behavior_data, Type == "Behavior"))

# Run Tukey's HSD
behavior_tukey <- TukeyHSD(behavior_aov)

# Extract only the comparisons we want
selected_comparisons <- behavior_tukey$Treatment[c(
  "G2_Herbivore-G1_Herbivore",  # H1 vs H2
  "G1_Predator-G1_Herbivore",   # H1 vs P1
  "G2_Predator-G1_Predator"     # P1 vs P2
), ]

# Print results
print("Selected Treatment Comparisons for Behavior Type:")
print(selected_comparisons)

# Range Plot for behavior type
ggplot(density_estimate_ranges_behavior) +
  # Add half violin plot to the right
  geom_half_violin(data = subset(behavior_data, Type == "Behavior"),
                   aes(x = Treatment, y = Y_centered, fill = Treatment),
                   side = "r", alpha = 0.3, position = position_nudge(x = 0.07)) +
  # G1 Herbivore range
  geom_ellipse(data = subset(density_estimate_ranges_behavior, Treatment == "G1_Herbivore"),
               aes(x0 = 1 - 0.15, 
                   y0 = (avg_broad_low + avg_broad_high) / 2,
                   a = 0.15, b = (avg_broad_high - avg_broad_low) / 2, 
                   angle = 0), 
               fill = custom_colors[2], color = "black", alpha = 0.5) +
  geom_ellipse(data = subset(density_estimate_ranges_behavior, Treatment == "G1_Herbivore"),
               aes(x0 = 1 - 0.15, 
                   y0 = (avg_core_low + avg_core_high) / 2,
                   a = 0.15, b = (avg_core_high - avg_core_low) / 2, 
                   angle = 0), 
               fill = custom_colors[1], color = "black", alpha = 0.5) +
  # G2 Herbivore range
  geom_ellipse(data = subset(density_estimate_ranges_behavior, Treatment == "G2_Herbivore"),
               aes(x0 = 2 - 0.15, 
                   y0 = (avg_broad_low + avg_broad_high) / 2,
                   a = 0.15, b = (avg_broad_high - avg_broad_low) / 2, 
                   angle = 0), 
               fill = custom_colors[4], color = "black", alpha = 0.5) +
  geom_ellipse(data = subset(density_estimate_ranges_behavior, Treatment == "G2_Herbivore"),
               aes(x0 = 2 - 0.15, 
                   y0 = (avg_core_low + avg_core_high) / 2,
                   a = 0.15, b = (avg_core_high - avg_core_low) / 2, 
                   angle = 0), 
               fill = custom_colors[3], color = "black", alpha = 0.5) +
  # G1 Predator range
  geom_ellipse(data = subset(density_estimate_ranges_behavior, Treatment == "G1_Predator"),
               aes(x0 = 3 - 0.15, 
                   y0 = (avg_broad_low + avg_broad_high) / 2,
                   a = 0.15, b = (avg_broad_high - avg_broad_low) / 2, 
                   angle = 0), 
               fill = custom_colors[6], color = "black", alpha = 0.5) +
  geom_ellipse(data = subset(density_estimate_ranges_behavior, Treatment == "G1_Predator"),
               aes(x0 = 3 - 0.15, 
                   y0 = (avg_core_low + avg_core_high) / 2,
                   a = 0.15, b = (avg_core_high - avg_core_low) / 2, 
                   angle = 0), 
               fill = custom_colors[5], color = "black", alpha = 0.5) +
  # G2 Predator range
  geom_ellipse(data = subset(density_estimate_ranges_behavior, Treatment == "G2_Predator"),
               aes(x0 = 4 - 0.15, 
                   y0 = (avg_broad_low + avg_broad_high) / 2,
                   a = 0.15, b = (avg_broad_high - avg_broad_low) / 2, 
                   angle = 0), 
               fill = custom_colors[8], color = "black", alpha = 0.5) +
  geom_ellipse(data = subset(density_estimate_ranges_behavior, Treatment == "G2_Predator"),
               aes(x0 = 4 - 0.15, 
                   y0 = (avg_core_low + avg_core_high) / 2,
                   a = 0.15, b = (avg_core_high - avg_core_low) / 2, 
                   angle = 0), 
               fill = custom_colors[7], color = "black", alpha = 0.5) +
  # Add mean points and text - adjusted x position to match ellipses
  {
    behavior_means <- subset(behavior_data, Type == "Behavior") %>% 
      group_by(Treatment) %>% 
      summarize(mean_y = mean(Y_centered)) %>%
      mutate(x_pos = as.numeric(Treatment) - 0.15)  # Match ellipse x-position
    
    list(
      geom_point(data = behavior_means,
                 aes(x = x_pos, y = mean_y),
                 size = 3, color = "black"),
      geom_text(data = behavior_means,
                aes(x = x_pos, y = mean_y, 
                    label = round(mean_y, 1)),
                vjust = -1.5, size = 4)
    )
  } +
  # Aesthetics
  scale_fill_manual(values = c(
    custom_colors[1],  # G1 Herbivore
    custom_colors[3],  # G2 Herbivore
    custom_colors[5],  # G1 Predator
    custom_colors[7]   # G2 Predator
  )) +
  scale_x_discrete(labels = c("G1_Herbivore", "G2_Herbivore", 
                              "G1_Predator", "G2_Predator")) +
  labs(title = "Core (50%) and Broad Ranges (95%) by Treatment", 
       subtitle = "Behavior Populations Only",
       x = "Treatment", 
       y = "Population-Centered Canopy Height") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 12),
        legend.position = "none")



#### Physiology Type Analysis ----


print("Number of observations by Type:")
print(table(behavior_data$Type))

# Check data structure for behavior type
physiology_subset <- behavior_data %>%
  filter(Type == "Physiology")
print("Number of observations by Treatment for physiology type:")
print(table(behavior_subset$Treatment))

# Check how many observations we have per ID
id_counts <- physiology_subset %>%
  group_by(id) %>%
  summarise(n = n())
print("Observations per ID for behavior type:")
print(summary(id_counts$n))

# Calculate density estimates for behavior type
density_estimate_ranges_physiology <- behavior_data %>%
  filter(Type == "Physiology") %>%
  group_by(id) %>%
  filter(n() >= 10) %>%  # Ensure enough points for density estimation
  do({
    # Get the data for this group
    y_values <- .$Y_centered
    
    tryCatch({
      # Use base R density estimation
      density_estimate <- density(y_values, 
                                  from = 0,  # Force minimum at 0
                                  to = max(y_values),
                                  n = 512,
                                  kernel = "gaussian")
      
      # Calculate cumulative density
      cum_density <- cumsum(density_estimate$y) / sum(density_estimate$y)
      
      # Find the ranges for 50% and 95% density
      core_idx <- which(cum_density <= 0.50)
      broad_idx <- which(cum_density <= 0.95)
      
      if(length(core_idx) > 0 && length(broad_idx) > 0) {
        data.frame(
          core_low = min(density_estimate$x[core_idx]),
          core_high = max(density_estimate$x[core_idx]),
          broad_low = min(density_estimate$x[broad_idx]),
          broad_high = max(density_estimate$x[broad_idx]),
          Treatment = .$Treatment[1]
        )
      } else {
        data.frame(
          core_low = NA_real_,
          core_high = NA_real_,
          broad_low = NA_real_,
          broad_high = NA_real_,
          Treatment = .$Treatment[1]
        )
      }
    }, error = function(e) {
      data.frame(
        core_low = NA_real_,
        core_high = NA_real_,
        broad_low = NA_real_,
        broad_high = NA_real_,
        Treatment = .$Treatment[1]
      )
    })
  }) %>%
  filter(!is.na(core_low)) %>%  # Remove any failed estimations
  group_by(Treatment) %>%
  summarize(
    avg_core_low = mean(core_low),
    avg_core_high = mean(core_high),
    avg_broad_low = mean(broad_low),
    avg_broad_high = mean(broad_high)
  ) %>%
  mutate(Treatment = factor(Treatment,
                            levels = c("G1_Herbivore", "G2_Herbivore", 
                                       "G1_Predator", "G2_Predator")))

# Check the results
print("Physiology type density estimates:")
print(density_estimate_ranges_physiology)


# Run ANOVA on behavior type data only
physiology_aov <- aov(Y_centered ~ Treatment, 
                      data = subset(behavior_data, Type == "Physiology"))

# Run Tukey's HSD
physiology_tukey <- TukeyHSD(physiology_aov)

# Extract only the comparisons we want
selected_comparisons <- physiology_tukey$Treatment[c(
  "G2_Herbivore-G1_Herbivore",  # H1 vs H2
  "G1_Predator-G1_Herbivore",   # H1 vs P1
  "G2_Predator-G1_Predator"     # P1 vs P2
), ]

# Print results
print("Selected Treatment Comparisons for Physiology Type:")
print(selected_comparisons)

# Range Plot for behavior type
ggplot(density_estimate_ranges_physiology) +
  # Add half violin plot to the right
  geom_half_violin(data = subset(behavior_data, Type == "Physiology"),
                   aes(x = Treatment, y = Y_centered, fill = Treatment),
                   side = "r", alpha = 0.3, position = position_nudge(x = 0.07)) +
  # G1 Herbivore range
  geom_ellipse(data = subset(density_estimate_ranges_physiology, Treatment == "G1_Herbivore"),
               aes(x0 = 1 - 0.15, 
                   y0 = (avg_broad_low + avg_broad_high) / 2,
                   a = 0.15, b = (avg_broad_high - avg_broad_low) / 2, 
                   angle = 0), 
               fill = custom_colors[2], color = "black", alpha = 0.5) +
  geom_ellipse(data = subset(density_estimate_ranges_physiology, Treatment == "G1_Herbivore"),
               aes(x0 = 1 - 0.15, 
                   y0 = (avg_core_low + avg_core_high) / 2,
                   a = 0.15, b = (avg_core_high - avg_core_low) / 2, 
                   angle = 0), 
               fill = custom_colors[1], color = "black", alpha = 0.5) +
  # G2 Herbivore range
  geom_ellipse(data = subset(density_estimate_ranges_physiology, Treatment == "G2_Herbivore"),
               aes(x0 = 2 - 0.15, 
                   y0 = (avg_broad_low + avg_broad_high) / 2,
                   a = 0.15, b = (avg_broad_high - avg_broad_low) / 2, 
                   angle = 0), 
               fill = custom_colors[4], color = "black", alpha = 0.5) +
  geom_ellipse(data = subset(density_estimate_ranges_physiology, Treatment == "G2_Herbivore"),
               aes(x0 = 2 - 0.15, 
                   y0 = (avg_core_low + avg_core_high) / 2,
                   a = 0.15, b = (avg_core_high - avg_core_low) / 2, 
                   angle = 0), 
               fill = custom_colors[3], color = "black", alpha = 0.5) +
  # G1 Predator range
  geom_ellipse(data = subset(density_estimate_ranges_physiology, Treatment == "G1_Predator"),
               aes(x0 = 3 - 0.15, 
                   y0 = (avg_broad_low + avg_broad_high) / 2,
                   a = 0.15, b = (avg_broad_high - avg_broad_low) / 2, 
                   angle = 0), 
               fill = custom_colors[6], color = "black", alpha = 0.5) +
  geom_ellipse(data = subset(density_estimate_ranges_physiology, Treatment == "G1_Predator"),
               aes(x0 = 3 - 0.15, 
                   y0 = (avg_core_low + avg_core_high) / 2,
                   a = 0.15, b = (avg_core_high - avg_core_low) / 2, 
                   angle = 0), 
               fill = custom_colors[5], color = "black", alpha = 0.5) +
  # G2 Predator range
  geom_ellipse(data = subset(density_estimate_ranges_physiology, Treatment == "G2_Predator"),
               aes(x0 = 4 - 0.15, 
                   y0 = (avg_broad_low + avg_broad_high) / 2,
                   a = 0.15, b = (avg_broad_high - avg_broad_low) / 2, 
                   angle = 0), 
               fill = custom_colors[8], color = "black", alpha = 0.5) +
  geom_ellipse(data = subset(density_estimate_ranges_physiology, Treatment == "G2_Predator"),
               aes(x0 = 4 - 0.15, 
                   y0 = (avg_core_low + avg_core_high) / 2,
                   a = 0.15, b = (avg_core_high - avg_core_low) / 2, 
                   angle = 0), 
               fill = custom_colors[7], color = "black", alpha = 0.5) +
  # Add mean points and text - adjusted x position to match ellipses
  {
    physiology_means <- subset(behavior_data, Type == "Physiology") %>% 
      group_by(Treatment) %>% 
      summarize(mean_y = mean(Y_centered)) %>%
      mutate(x_pos = as.numeric(Treatment) - 0.15)  # Match ellipse x-position
    
    list(
      geom_point(data = physiology_means,
                 aes(x = x_pos, y = mean_y),
                 size = 3, color = "black"),
      geom_text(data = physiology_means,
                aes(x = x_pos, y = mean_y, 
                    label = round(mean_y, 1)),
                vjust = -1.5, size = 4)
    )
  } +
  # Aesthetics
  scale_fill_manual(values = c(
    custom_colors[1],  # G1 Herbivore
    custom_colors[3],  # G2 Herbivore
    custom_colors[5],  # G1 Predator
    custom_colors[7]   # G2 Predator
  )) +
  scale_x_discrete(labels = c("G1_Herbivore", "G2_Herbivore", 
                              "G1_Predator", "G2_Predator")) +
  labs(title = "Core (50%) and Broad Ranges (95%) by Treatment", 
       subtitle = "Physiology Populations Only",
       x = "Treatment", 
       y = "Population-Centered Canopy Height") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 12),
        legend.position = "none")





#### Type Comparison Analysis ----

# Analysis for Type differences within each Treatment
# Function to run analysis for each treatment
analyze_type_difference <- function(treatment_level) {
  # Subset data for this treatment
  treatment_data <- subset(behavior_data, Treatment == treatment_level)
  
  # Run ANOVA
  type_aov <- aov(Y_centered ~ Type, data = treatment_data)
  
  # Get summary statistics
  type_means <- treatment_data %>%
    group_by(Type) %>%
    summarize(
      mean = mean(Y_centered),
      sd = sd(Y_centered),
      n = n()
    )
  
  # Get ANOVA results
  aov_results <- summary(type_aov)
  
  # Calculate effect size (Cohen's d)
  means_diff <- diff(tapply(treatment_data$Y_centered, treatment_data$Type, mean))
  pooled_sd <- sqrt(mean(tapply(treatment_data$Y_centered, treatment_data$Type, var)))
  cohens_d <- means_diff / pooled_sd
  
  # Return results
  list(
    treatment = treatment_level,
    aov_results = aov_results,
    means = type_means,
    cohens_d = cohens_d,
    p_value = aov_results[[1]]$`Pr(>F)`[1]
  )
}

# Run analysis for each treatment
treatments <- c("G1_Herbivore", "G2_Herbivore", "G1_Predator", "G2_Predator")
type_analyses <- lapply(treatments, analyze_type_difference)

# Print results
for(analysis in type_analyses) {
  cat("\nResults for", analysis$treatment, ":\n")
  cat("Mean differences between Types:\n")
  print(analysis$means)
  cat("ANOVA p-value:", analysis$p_value, "\n")
  cat("Effect size (Cohen's d):", analysis$cohens_d, "\n")
  cat("----------------------------------------\n")
}


##### BA for Type -----

# Calculate densities for each treatment-type combination
density_H1B <- stats::density(behavior_data$Y_centered[behavior_data$Treatment == "G1_Herbivore" & 
                                                         behavior_data$Type == "Behavior"],
                              from = y_min, to = y_max, n = 100)
density_H1P <- stats::density(behavior_data$Y_centered[behavior_data$Treatment == "G1_Herbivore" & 
                                                         behavior_data$Type == "Physiology"],
                              from = y_min, to = y_max, n = 100)
density_H2B <- stats::density(behavior_data$Y_centered[behavior_data$Treatment == "G2_Herbivore" & 
                                                         behavior_data$Type == "Behavior"],
                              from = y_min, to = y_max, n = 100)
density_H2P <- stats::density(behavior_data$Y_centered[behavior_data$Treatment == "G2_Herbivore" & 
                                                         behavior_data$Type == "Physiology"],
                              from = y_min, to = y_max, n = 100)
density_P1B <- stats::density(behavior_data$Y_centered[behavior_data$Treatment == "G1_Predator" & 
                                                         behavior_data$Type == "Behavior"],
                              from = y_min, to = y_max, n = 100)
density_P1P <- stats::density(behavior_data$Y_centered[behavior_data$Treatment == "G1_Predator" & 
                                                         behavior_data$Type == "Physiology"],
                              from = y_min, to = y_max, n = 100)
density_P2B <- stats::density(behavior_data$Y_centered[behavior_data$Treatment == "G2_Predator" & 
                                                         behavior_data$Type == "Behavior"],
                              from = y_min, to = y_max, n = 100)
density_P2P <- stats::density(behavior_data$Y_centered[behavior_data$Treatment == "G2_Predator" & 
                                                         behavior_data$Type == "Physiology"],
                              from = y_min, to = y_max, n = 100)



# Convert all to kde format
densities <- list(
  H1B = list(eval.points = density_H1B$x, estimate = density_H1B$y),
  H1P = list(eval.points = density_H1P$x, estimate = density_H1P$y),
  H2B = list(eval.points = density_H2B$x, estimate = density_H2B$y),
  H2P = list(eval.points = density_H2P$x, estimate = density_H2P$y),
  P1B = list(eval.points = density_P1B$x, estimate = density_P1B$y),
  P1P = list(eval.points = density_P1P$x, estimate = density_P1P$y),
  P2B = list(eval.points = density_P2B$x, estimate = density_P2B$y),
  P2P = list(eval.points = density_P2P$x, estimate = density_P2P$y)
)

# Calculate cumulative densities and ranges
cum_densities <- lapply(densities, function(x) cumsum(x$estimate) / sum(x$estimate))
core_ranges <- lapply(seq_along(densities), function(i) {
  range(densities[[i]]$eval.points[cum_densities[[i]] <= 0.5])
})
broad_ranges <- lapply(seq_along(densities), function(i) {
  range(densities[[i]]$eval.points[cum_densities[[i]] <= 0.95])
})
names(core_ranges) <- names(densities)
names(broad_ranges) <- names(densities)


# Calculate BAs for type comparisons
# H1 Behavior vs Physiology
core_range_H1 <- c(max(core_ranges$H1B[1], core_ranges$H1P[1]),
                   min(core_ranges$H1B[2], core_ranges$H1P[2]))
broad_range_H1 <- c(max(broad_ranges$H1B[1], broad_ranges$H1P[1]),
                    min(broad_ranges$H1B[2], broad_ranges$H1P[2]))
BA_H1_core <- calculate_bhattacharyya_affinity_isopleth(densities$H1B, densities$H1P, core_range_H1)
BA_H1_broad <- calculate_bhattacharyya_affinity_isopleth(densities$H1B, densities$H1P, broad_range_H1)

# H2 Behavior vs Physiology
core_range_H2 <- c(max(core_ranges$H2B[1], core_ranges$H2P[1]),
                   min(core_ranges$H2B[2], core_ranges$H2P[2]))
broad_range_H2 <- c(max(broad_ranges$H2B[1], broad_ranges$H2P[1]),
                    min(broad_ranges$H2B[2], broad_ranges$H2P[2]))
BA_H2_core <- calculate_bhattacharyya_affinity_isopleth(densities$H2B, densities$H2P, core_range_H2)
BA_H2_broad <- calculate_bhattacharyya_affinity_isopleth(densities$H2B, densities$H2P, broad_range_H2)

# P1 Behavior vs Physiology
core_range_P1 <- c(max(core_ranges$P1B[1], core_ranges$P1P[1]),
                   min(core_ranges$P1B[2], core_ranges$P1P[2]))
broad_range_P1 <- c(max(broad_ranges$P1B[1], broad_ranges$P1P[1]),
                    min(broad_ranges$P1B[2], broad_ranges$P1P[2]))
BA_P1_core <- calculate_bhattacharyya_affinity_isopleth(densities$P1B, densities$P1P, core_range_P1)
BA_P1_broad <- calculate_bhattacharyya_affinity_isopleth(densities$P1B, densities$P1P, broad_range_P1)

# P2 Behavior vs Physiology
core_range_P2 <- c(max(core_ranges$P2B[1], core_ranges$P2P[1]),
                   min(core_ranges$P2B[2], core_ranges$P2P[2]))
broad_range_P2 <- c(max(broad_ranges$P2B[1], broad_ranges$P2P[1]),
                    min(broad_ranges$P2B[2], broad_ranges$P2P[2]))
BA_P2_core <- calculate_bhattacharyya_affinity_isopleth(densities$P2B, densities$P2P, core_range_P2)
BA_P2_broad <- calculate_bhattacharyya_affinity_isopleth(densities$P2B, densities$P2P, broad_range_P2)


# Type comparisons
cat("Type Comparisons:\n")
cat("G1 Herbivore - Behavior vs Physiology:\n")
cat("Core BA:", BA_H1_core, "\n")
cat("Broad BA:", BA_H1_broad, "\n\n")

cat("G2 Herbivore - Behavior vs Physiology:\n")
cat("Core BA:", BA_H2_core, "\n")
cat("Broad BA:", BA_H2_broad, "\n\n")

cat("G1 Predator - Behavior vs Physiology:\n")
cat("Core BA:", BA_P1_core, "\n")
cat("Broad BA:", BA_P1_broad, "\n\n")

cat("G2 Predator - Behavior vs Physiology:\n")
cat("Core BA:", BA_P2_core, "\n")
cat("Broad BA:", BA_P2_broad, "\n\n")

##### Figure -----
type_comparison_plot <- function(treatment_level) {
  # Reuse existing density ranges but filter for specific treatment
  density_ranges_behavior <- density_estimate_ranges_behavior %>%
    filter(Treatment == treatment_level) %>%
    mutate(Type = "Behavior")
  density_ranges_physiology <- density_estimate_ranges_physiology %>%
    filter(Treatment == treatment_level) %>%
    mutate(Type = "Physiology")
  density_ranges <- bind_rows(density_ranges_behavior, density_ranges_physiology)
  
  # Determine colors based on treatment
  if(treatment_level == "G1_Herbivore") {
    colors <- c(custom_colors[1], custom_colors[2])  # Yellow shades
  } else if(treatment_level == "G2_Herbivore") {
    colors <- c(custom_colors[3], custom_colors[4])  # Green shades
  } else if(treatment_level == "G1_Predator") {
    colors <- c(custom_colors[5], custom_colors[6])  # Purple shades
  } else {
    colors <- c(custom_colors[7], custom_colors[8])  # Blue shades
  }
  
  ggplot() +
    # Add half violin plots
    geom_half_violin(data = subset(behavior_data, Treatment == treatment_level & Type == "Behavior"),
                     aes(x = Type, y = Y_centered),
                     side = "l", alpha = 0.3, fill = colors[1],
                     position = position_nudge(x = -0.33)) +
    geom_half_violin(data = subset(behavior_data, Treatment == treatment_level & Type == "Physiology"),
                     aes(x = Type, y = Y_centered),
                     side = "r", alpha = 0.3, fill = colors[1],
                     position = position_nudge(x = .03)) +
    # Add ellipses for each type
    geom_ellipse(data = subset(density_ranges, Type == "Behavior"),
                 aes(x0 = 1 - 0.15,
                     y0 = (avg_broad_low + avg_broad_high) / 2,
                     a = 0.15, b = (avg_broad_high - avg_broad_low) / 2,
                     angle = 0),
                 fill = colors[2], color = "black", alpha = 0.5) +
    geom_ellipse(data = subset(density_ranges, Type == "Behavior"),
                 aes(x0 = 1 - 0.15,
                     y0 = (avg_core_low + avg_core_high) / 2,
                     a = 0.15, b = (avg_core_high - avg_core_low) / 2,
                     angle = 0),
                 fill = colors[1], color = "black", alpha = 0.5) +
    geom_ellipse(data = subset(density_ranges, Type == "Physiology"),
                 aes(x0 = 2 - 0.15,
                     y0 = (avg_broad_low + avg_broad_high) / 2,
                     a = 0.15, b = (avg_broad_high - avg_broad_low) / 2,
                     angle = 0),
                 fill = colors[2], color = "black", alpha = 0.5) +
    geom_ellipse(data = subset(density_ranges, Type == "Physiology"),
                 aes(x0 = 2 - 0.15,
                     y0 = (avg_core_low + avg_core_high) / 2,
                     a = 0.15, b = (avg_core_high - avg_core_low) / 2,
                     angle = 0),
                 fill = colors[1], color = "black", alpha = 0.5) +
    # Add mean points and text
    {
      type_means <- behavior_data %>% 
        filter(Treatment == treatment_level) %>%
        group_by(Type) %>% 
        summarize(mean_y = mean(Y_centered)) %>%
        mutate(x_pos = as.numeric(factor(Type)) - 0.15)
      
      list(
        geom_point(data = type_means,
                   aes(x = x_pos, y = mean_y),
                   size = 2.5, color = "black"),
        geom_text(data = type_means,
                  aes(x = x_pos, y = mean_y, 
                      label = round(mean_y, 1)),
                  vjust = -1.5, size = 3)
      )
    } +
    labs(title = treatment_level,
         x = "Type",
         y = "Population-Centered Canopy Height") +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size = 9),
          legend.position = "none")
}

# Create all four plots
p1 <- type_comparison_plot("G1_Herbivore")
p2 <- type_comparison_plot("G2_Herbivore")
p3 <- type_comparison_plot("G1_Predator")
p4 <- type_comparison_plot("G2_Predator")

# Combine plots using patchwork
library(patchwork)
combined_plot <- (p1 + p2) / (p3 + p4)

# Add overall title
combined_plot + 
  plot_annotation(
    title = "Comparison of Behaviorial and Physiologial Populations by Generation",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))
  )





# Respiration ----

# Packages
library(dplyr)
library(lme4)
library(ggplot2)

### Data Processing ----

Respiration_Full<-read.csv("Data/Respriation_2023_F1_F2_Raw.csv")
Respiration_Full<-Respiration_Full[,-1]

Tempcatlist <- data.frame(matrix(ncol = 1, nrow = 2))
Tempcatlist[,1]<-c("25C","30C")
# For internal consistency with measurements 2022, MB measured these grasshoppers at a 2min dwell time 
# This lead to frequent saturation of the censor making the 35C data patchy. 
# Here MB included only the data for 25C and 30C

Final_output_all_temps <- data.frame(matrix(ncol = 9, nrow = 126))

for (Z in 1:2){
  Final_output_Per_Temp <- data.frame(matrix(ncol = 8, nrow = 63))
  for(R in 1:9){
    Raw<-subset(Respiration_Full, Run==R & Tempcat==Tempcatlist[Z,1])
    
    names(Raw)[1] <- "time"
    names(Raw)[65] <- "respiration"
    output<-subset(Raw,time>18& time<33&respiration>0)
    Logger_Pro_Calcd_Resp<-output[,65]
    Logger_Pro_Calcd_Resp
    
    df <- data.frame(matrix(ncol = 6, nrow = 7))
    
    colnames(df) <- c('delay_low','delay_high','end','Stop_manual','Temp_stop_14mins','Logger_Pro_Calcd')
    df$delay_low<-c(18,20,22,24,26,28,30)
    df$delay_low<-df$delay_low+0.2
    df$delay_high<-c(18,20,22,24,26,28,30)
    df$delay_high<-df$delay_high+1.5
    df$end<-c(18,20,22,24,26,28,30)
    df$end<-df$end+2
    
    names(Raw)[1] <- "time"
    names(Raw)[2] <- "CO2"
    names(Raw)[4] <- "temp"
    names(Raw)[15] <- "flow"
    
    
    for(x in 1:7) {
      
      
      Pulse_and_Background<-subset(Raw,time>df[x,1]&time<=df[x,2])
      Background<-subset(Raw,time>df[x,2]&time<=df[x,3])
      
      X<-(sum(Pulse_and_Background$CO2*(1/60))-(sum(Background$CO2*(1/60))/.5)*1.3)
      df[x,4]<-(((X*1000)/1000000)*mean(Pulse_and_Background$flow))/14
      
      Incubation<-subset(Raw,time>(df[x,3]-16)&time<=(df[x,3]))
      df[x,5]<-mean(Incubation$temp)
    }
    
    df[,6]<-(Logger_Pro_Calcd_Resp)/14
    df<-df[,4:6]
    df$chamber<-c(2,3,4,5,6,7,8)
    
    df$Run<-R
    Final_output_Per_Temp[(((R)*7)-6):((R)*7),]<-df
    
    Final_output_Per_Temp$TempCat<-Tempcatlist[Z,1]
  }
  Final_output_all_temps[((Z*63)-62):((Z*63)),]<-Final_output_Per_Temp[1:63,]
}
Final_output_all_temps<-Final_output_all_temps[,4:9]
colnames(Final_output_all_temps) <- c('Chamber','Run','Manual_Resp_Calc','Incubation_Temp','Logger_Pro_Calcd','Temp_Catagorical')


# Logger Pro automatically calculates the integral for respiration (Logger_Pro_Calcd) 
# MB also calculated the values manually (Manual_Resp_Calc)

Respiration_Mass_2023<-read.csv("Data/Respriation_Mass_2023_F1_F2.csv")
Respiration_Mass_2023<-arrange(Respiration_Mass_2023, Run,Chamber)
Respiration_Mass_2023_times2<-rbind(Respiration_Mass_2023,Respiration_Mass_2023)
respiration_data <-cbind(Final_output_all_temps,Respiration_Mass_2023_times2)

# Fix duplicate columns by selecting the first instance of each column name
respiration_data <- respiration_data[, !duplicated(names(respiration_data))]

# Append gen 1 herbivores:
# Load the new data
Gen1 <- read.csv("Data/UP_Gen_1_Herbivores_2023.csv")

# Check the structure of the new data
str(Gen1)

# Rename columns in the new data to match respiration_data
Gen1 <- Gen1 %>%
  rename(
    Entry_Order = Entry.Order,
    Saturated_at_30 = Saturated
  )

# Ensure the new data has the same columns as respiration_data
# Select only the columns that exist in respiration_data
Gen1 <- Gen1 %>%
  select(names(respiration_data))

# Append the new data to respiration_data
respiration_data <- rbind(respiration_data, Gen1)

# Check the structure of the updated respiration_data
str(respiration_data)




# Cleaning
respiration_data <- respiration_data %>%
  dplyr::mutate(Generation = case_when(
    Generation == 1 ~ "G1",
    Generation == 2 ~ "G2",
    TRUE ~ as.character(Generation)
  )) %>%
  dplyr::filter(
    Exclude != 1
  ) %>%
  dplyr::select(Logger_Pro_Calcd, 
                Temp_Catagorical, 
                Generation, 
                Population, 
                Predation, 
                Individual, 
                Mass) %>%
  dplyr::rename(Temp = Temp_Catagorical)

respiration_data <- respiration_data %>%
  mutate(Treatment = paste(Generation, Predation, sep = "_")) %>%
  mutate(Treatment = factor(Treatment, 
                           levels = c("G1_Herbivore", "G2_Herbivore", 
                                    "G1_Predator", "G2_Predator"))) %>% 
  mutate(Type = case_when(
    Population %in% c("UP", "MC", "FN", "HF") ~ "Behavior",
    Population %in% c("DC", "SP", "SC", "YF") ~ "Physiology"))

# Calculate SMR
respiration_data <- respiration_data %>%
  dplyr::mutate(SMR = Logger_Pro_Calcd / Mass)  # units are uLCO2/min/g

### Analysis ----

# Center temperature to make intercepts more interpretable
respiration_data <- respiration_data %>%
  mutate(Temp_centered = case_when(
    Temp == "25C" ~ -2.5,
    Temp == "30C" ~ 2.5
  ))


# Test for differences between treatments at 25C
respiration_25C <- respiration_data %>%
  filter(Temp == "25C")

aov_25C <- aov(SMR ~ Treatment*Type, data = respiration_25C)
summary(aov_25C)
# Post-hoc if ANOVA is significant
# TukeyHSD(aov_25C)

# Updated plot for 25C
ggplot(respiration_25C, aes(x = Treatment, y = SMR, fill = Treatment)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_fill_manual(values = c(custom_colors[1], custom_colors[3], custom_colors[5], custom_colors[7])) +
  labs(title = "Metabolic Rate by Treatment at 25°C", 
       x = "Treatment", 
       y = expression("Mass-Specific Metabolic Rate (µL CO"[2]*"/min/g)")) +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        text = element_text(size = 12),
        legend.position = "none") +
  facet_wrap(~Type)

# Test for differences between treatments at 30C
respiration_30C <- respiration_data %>%
  filter(Temp == "30C")

# ANOVA instead of t-test for three groups
aov_30C <- aov(SMR ~ Treatment*Type, data = respiration_30C)
summary(aov_30C)
# Post-hoc if ANOVA is significant
#TukeyHSD(aov_30C)

# Updated plot for 30C
ggplot(respiration_30C, aes(x = Treatment, y = SMR, fill = Treatment)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_fill_manual(values = c(custom_colors[1], custom_colors[3], custom_colors[5], custom_colors[7])) +
  labs(title = "Metabolic Rate by Treatment at 30°C", 
       x = "Treatment", 
       y = expression("Mass-Specific Metabolic Rate (µL CO"[2]*"/min/g)")) +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        text = element_text(size = 12),
        legend.position = "none") +
  facet_wrap(~Type)



# Test for differences in thermal reaction norms between Types within each Treatment

# Create separate models for each Treatment
# G1 Predator
g1p_data <- respiration_data %>% 
  filter(Treatment == "G1_Predator")
g1p_model <- lmer(SMR ~ Type * Temp_centered + (1|Individual), data = g1p_data)

summary(g1p_model)
anova(g1p_model)

# G2 Predator
g2p_data <- respiration_data %>% 
  filter(Treatment == "G2_Predator")
g2p_model <- lmer(SMR ~ Type * Temp_centered + (1|Individual), data = g2p_data)
summary(g2p_model)
anova(g2p_model)

# G1 Herbivore
g1h_data <- respiration_data %>% 
  filter(Treatment == "G1_Herbivore")
g1h_model <- lmer(SMR ~ Type * Temp_centered + (1|Individual), data = g1h_data)
summary(g1h_model)
anova(g1h_model)

# G2 Herbivore
# insufficient data, no Type of 'Behavioral'

# Print interpretable results for each comparison
cat("\nG1 Predator - Behavior vs Physiology:\n")
cat("Main effect of Type:", fixef(g1p_model)["TypePhysiology"], "\n")
cat("Interaction (difference in slopes):", fixef(g1p_model)["TypePhysiology:Temp_centered"], "\n")
cat("P-values from ANOVA:\n")
print(anova(g1p_model))

cat("\nG2 Predator - Behavior vs Physiology:\n")
cat("Main effect of Type:", fixef(g2p_model)["TypePhysiology"], "\n")
cat("Interaction (difference in slopes):", fixef(g2p_model)["TypePhysiology:Temp_centered"], "\n")
cat("P-values from ANOVA:\n")
print(anova(g2p_model))

cat("\nG1 Herbivore - Behavior vs Physiology:\n")
cat("Main effect of Type:", fixef(g1h_model)["TypePhysiology"], "\n")
cat("Interaction (difference in slopes):", fixef(g1h_model)["TypePhysiology:Temp_centered"], "\n")
cat("P-values from ANOVA:\n")
print(anova(g1h_model))


### Reaction Norm Plot by Type ----
ggplot(respiration_data, aes(x = Temp, y = SMR, 
                            group = Individual, color = Treatment)) +
  # Add individual reaction norms
  geom_line(alpha = 0.2) +
  geom_point(size = 2, alpha = 0.4) +
  # Add treatment means and SE
  stat_summary(aes(group = Treatment), 
              fun.data = "mean_se",
              geom = "errorbar",
              width = 0.07,
              linewidth = 1,
              alpha = 1) +
  stat_summary(aes(group = Treatment), 
              fun = mean,
              geom = "line",
              linewidth = 1.5,
              alpha = 0.8) +
  stat_summary(aes(group = Treatment), 
              fun = mean,
              geom = "point",
              size = 3,
              alpha = 1) +
  # Aesthetics
  scale_color_manual(values = c(
    "G1_Herbivore" = custom_colors[1],  # Yellow
    "G2_Herbivore" = custom_colors[3],  # Green
    "G1_Predator" = custom_colors[5],   # Purple
    "G2_Predator" = custom_colors[7]    # Blue
  )) +
  labs(title = "Thermal Reaction Norms by Treatment", 
       x = "Temperature (°C)", 
       y = expression("Mass-Specific Metabolic Rate (µL CO"[2]*"/min/g)")) +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        text = element_text(size = 12),
        legend.position = "bottom",
        legend.title = element_blank(),
        strip.text = element_text(face = "bold")) +
  # Add faceting by Type
  facet_wrap(~Type, scales = "free_y", 
             labeller = as_labeller(c(
               Behavior = "Behavior Type Populations",
               Physiology = "Physiology Type Populations"
             )))


# Get sample sizes for each Type*Treatment combination
respiration_data %>%
  group_by(Type, Treatment) %>%
  summarize(
    n_individuals = n_distinct(Individual),
    .groups = 'drop'
  ) %>%
  arrange(Type, Treatment) %>%
  print(n = Inf)  # prints all rows




