# Title: Predator Effects Magnified in Ecosystems
# Script: Analyzes trait data
# Authors: N.Sommer and M.Baker 


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
library(DHARMa)

### Data Processing ----
Behavior_2023 <- read.csv("Data/Behavior.csv")

# Adjusted start time for sparse observations
Behavior_2023$Time <- hms(Behavior_2023$Time)
behavior_data <- subset(Behavior_2023,Time>(hms("8:40:00")))

# Import temperature logger data
behavior_temp <- read.csv("Data/Behavior_Temps.csv")

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

# Update any generation label
behavior_data <- behavior_data %>%
  mutate(Generation = case_when(
    Generation == "F1" ~ "G1",
    Generation == "F2" ~ "G2",
    TRUE ~ as.character(Generation)
  ))

# Create the combined treatment column
behavior_data <- behavior_data %>%
  mutate(Treatment = paste(Generation)) %>%
  mutate(Treatment = factor(Treatment, 
                            levels = c("G1", "G2")))


# Check population effect
height_model <- lmer(Y ~ Treatment + (1|Population), 
                     data = behavior_data)

VarCorr(height_model)

# Center population heights
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

# Check the input data
print("Input data summary:")
print(summary(behavior_data$Y_centered))
print(table(behavior_data$Treatment))

#' Calculate density estimates and ranges for behavioral data
#' @param data Dataframe containing Y_centered values and grouping variables
#' @param min_obs Minimum observations required per group (default 10)
#' @return Dataframe with core and broad ranges for each group
calculate_density_ranges <- function(data, min_obs = 10) {
  data %>%
    group_by(Treatment) %>%
    filter(n() >= min_obs) %>%
    summarize(
      # Calculate quantiles directly
      avg_core_low = quantile(Y_centered, probs = 0.25),
      avg_core_high = quantile(Y_centered, probs = 0.75),
      avg_broad_low = quantile(Y_centered, probs = 0.025),
      avg_broad_high = quantile(Y_centered, probs = 0.975),
      
      # Mean for positioning
      Y_centered = mean(Y_centered),
      .groups = 'drop'
    )
}

#' Print diagnostic information for density ranges
#' @param data The density range data
print_diagnostic <- function(data) {
  data %>%
    group_by(Treatment) %>%
    summarize(
      core_range = paste(round(avg_core_low, 1), "to", round(avg_core_high, 1)),
      core_width = round(avg_core_high - avg_core_low, 1),
      broad_range = paste(round(avg_broad_low, 1), "to", round(avg_broad_high, 1)),
      broad_width = round(avg_broad_high - avg_broad_low, 1),
      mean_height = round(Y_centered, 1)
    ) %>%
    print(n = Inf)
}

# Calculate ranges and print diagnostics
density_estimate_ranges <- calculate_density_ranges(behavior_data)
print("\nDensity Range Diagnostics:")
print_diagnostic(density_estimate_ranges)

### Analysis ----

#### LME repeated measures ----
height_model_rm <- lmer(Y_centered ~ Treatment + (1|id), 
                        data = behavior_data)

summary(height_model_rm)

# emmeans post-hoc contrasts
emmeans(height_model_rm, pairwise ~ Treatment)

#### Distribution test -----

# Wilcoxon rank sum test for G1 vs G2 predators
wilcox_test <- wilcox.test(
  behavior_data$Y_centered[behavior_data$Treatment == "G1"],
  behavior_data$Y_centered[behavior_data$Treatment == "G2"]
)

# Print results
print("Wilcoxon test results for G1 vs G2 predators:")
print(wilcox_test)

## Figures ----

# Define custom colors
custom_colors <- c(
  "G1" = "#C69648",  # G1 - golden
  "G2" = "#A15A29"   # G2 - purple
)

# Define treatment labels for plots
treatment_labels <- c(
  "G1" = "G[1]",
  "G2" = "G[2]"
)

#' Create aggregated panel plot for all data
create_aggregated_plot <- function(data) {
  # Calculate means for the plot
  means <- data %>%
    group_by(Treatment) %>%
    summarise(
      Y_centered = mean(Y_centered),
      .groups = 'drop'
    )
  
  ggplot(data, aes(x = Treatment, y = Y_centered)) +
    # Violin plots
    geom_violin(
      aes(fill = Treatment),
      alpha = 0.5,
      trim = TRUE
    ) +
    # Add jittered points
    geom_jitter(
      aes(color = Treatment),
      alpha = 0.2,
      size = 0.8,
      width = 0.2
    ) +
    # Mean points
    geom_point(
      data = means,
      size = 3,
      color = "black"
    ) +
    # Add mean value labels
    geom_text(
      data = means,
      aes(label = round(Y_centered, 1)),
      vjust = -1,
      size = 3.5
    ) +
    # Styling
    scale_y_continuous(limits = c(0, 80)) +
    scale_x_discrete(
      labels = parse(text = c(
        "G[1]",
        "G[2]"
      ))
    ) +
    scale_fill_manual(values = custom_colors) +
    scale_color_manual(values = custom_colors) +
    theme_light() +
    theme(
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_line(color = "grey95"),
      legend.position = "none",
      axis.text.x = element_text(size = 10)
    ) +
    labs(
      y = "Population-centered canopy height",
      x = NULL
    )
}

# Create all plots
aggregated_plot <- create_aggregated_plot(behavior_data)

# Temperature plot
temperature_plot <- ggplot(behavior_data, aes(x = Nearest_Temperature, y = Y, color = Treatment)) +
  geom_jitter(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  scale_color_manual(values = custom_colors,
                    labels = c(
                      expression(G[1]),
                      expression(G[2])
                    )) +
  labs(x = "Temperature (°C)", 
       y = "Height (cm)") +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        text = element_text(size = 10),
        strip.text = element_text(face = "bold"),
        legend.position = "bottom")

# Display plots
aggregated_plot
temperature_plot


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


# Calculate SMR
respiration_data <- respiration_data %>%
  dplyr::mutate(SMR = Logger_Pro_Calcd / Mass)  # units are uLCO2/min/g


# Extract predator G1-G2 data
respiration_data <- respiration_data %>% 
  dplyr::filter(Predation == "Predator")


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

aov_25C <- aov(SMR ~ Generation, data = respiration_25C)
summary(aov_25C)


# Plot for 25C
ggplot(respiration_25C, aes(x = Generation, y = SMR, fill = Generation)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_fill_manual(values = c(custom_colors[1], custom_colors[2])) +
  labs(title = "Metabolic Rate by Generation at 25°C", 
       x = "Generation", 
       y = expression("Mass-Specific Metabolic Rate (µL CO"[2]*"/min/g)")) +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        text = element_text(size = 12),
        legend.position = "none")

# Test for differences between treatments at 30C
respiration_30C <- respiration_data %>%
  filter(Temp == "30C")

# ANOVA instead of t-test for three groups
aov_30C <- aov(SMR ~ Generation, data = respiration_30C)
summary(aov_30C)

# Plot for 30C
ggplot(respiration_30C, aes(x = Generation, y = SMR, fill = Generation)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  scale_fill_manual(values = c(custom_colors[1], custom_colors[2])) +
  labs(title = "Metabolic Rate by Generation at 30°C", 
       x = "Generation", 
       y = expression("Mass-Specific Metabolic Rate (µL CO"[2]*"/min/g)")) +
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        text = element_text(size = 12),
        legend.position = "none")



# Test for differences in thermal reaction norms

# Separate models for each generation
# G1 Predator
g1p_data <- respiration_data %>% 
  filter(Generation == "G1")
g1p_model <- lmer(SMR ~ Temp_centered + (1 | Individual), data = g1p_data)

summary(g1p_model)

# G2 Predator
g2p_data <- respiration_data %>% 
  filter(Generation == "G2")
g2p_model <- lmer(SMR ~ Temp_centered + (1 | Individual), data = g2p_data)

summary(g2p_model)

g1g2p_model <- lmer(SMR ~ Temp_centered + Generation + (1 | Individual), data = respiration_data)
  
summary(g1g2p_model)

emmeans(g1g2p_model, pairwise ~ Generation)

## Reaction Norm Plot by Type ----
# Reaction Norm Plot
reaction_norm_plot <- ggplot(respiration_data, aes(x = Temp, y = SMR, group = Generation, color = Generation)) +
  geom_point(alpha = 0.2, size = 1) +
  geom_line(aes(group = Individual), alpha = 0.2) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.5) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  scale_color_manual(
    values = setNames(
      c(custom_colors[1], custom_colors[2]),
      unique(respiration_data$Generation)
    ),
    labels = expression(G[1], G[2])
  ) +
  scale_x_discrete(
    labels = c("25°C", "30°C"),
    expand = expansion(mult = 0.2)
  ) +
  theme_light() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(size = 12)
  ) +
  labs(
    x = "Temperature (°C)",
    y = expression(paste("Mass-specific metabolic rate (", mu, "L ", CO[2], "/min/g)")),
    tag = "A"
  )

behavior_plot <- create_aggregated_plot(behavior_data) +
  labs(tag = "B")

# Combine plots using patchwork
library(patchwork)
combined_plot <- reaction_norm_plot + behavior_plot +
  plot_layout(guides = "collect") &  # Collect legends together
  theme(legend.position = "bottom")

# Display combined plot
combined_plot



#  UNUSED CODE ----


# Check sample sizes for each Type*Treatment combination
respiration_data %>%
  group_by(Type, Treatment) %>%
  summarize(
    n_individuals = n_distinct(Individual),
    .groups = 'drop'
  ) %>%
  arrange(Type, Treatment) %>%
  print(n = Inf)  # prints all rows

behavior_data %>%
  group_by(Type, Treatment) %>%
  summarize(
    n_individuals = n_distinct(id),
    .groups = 'drop'
  ) %>%
  arrange(Type, Treatment) %>%
  print(n = Inf)  # prints all rows


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




#### Bhattacharyya's affinity analysis ----

#' Extract density values for a specific treatment
#' @param data Dataframe containing the data
#' @param treatment Treatment to extract
#' @return Vector of Y_centered values for the treatment
calculate_treatment_density <- function(data, treatment) {
  subset(data, Treatment == treatment)$Y_centered
}



#' Calculate Bhattacharyya's affinity between two treatments using core ranges
#' @param density1 First density estimate
#' @param density2 Second density estimate
#' @param range Range of values to consider
#' @return BA value between 0 and 1
calculate_bhattacharyya_affinity <- function(density1, density2, range) {
  # Get core ranges (25th to 75th percentiles)
  core1 <- quantile(density1, probs = c(0.25, 0.75))
  core2 <- quantile(density2, probs = c(0.25, 0.75))
  
  # Calculate overlap within core ranges
  overlap_min <- max(core1[1], core2[1])
  overlap_max <- min(core1[2], core2[2])
  
  if(overlap_max <= overlap_min) {
    return(0)  # No overlap in core ranges
  }
  
  # Calculate BA for core ranges
  x <- seq(overlap_min, overlap_max, length.out = 100)
  f1 <- density(density1, from = overlap_min, to = overlap_max)$y
  f2 <- density(density2, from = overlap_min, to = overlap_max)$y
  
  # Normalize densities
  f1 <- f1 / sum(f1)
  f2 <- f2 / sum(f2)
  
  # Calculate BA
  ba <- sum(sqrt(f1 * f2))
  
  return(ba)
}

# Calculate BA for each combination of treatments and types
calculate_ba_results <- function(data) {
  results <- data.frame(
    Type = character(),
    Comparison = character(),
    BA_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  for(type in unique(data$Type)) {
    type_data <- subset(data, Type == type)
    
    # Calculate densities for each treatment in this type
    g1h_density <- calculate_treatment_density(type_data, "G1_Herbivore")
    g1p_density <- calculate_treatment_density(type_data, "G1_Predator")
    g2h_density <- calculate_treatment_density(type_data, "G2_Herbivore")
    g2p_density <- calculate_treatment_density(type_data, "G2_Predator")
    
    # Calculate range for comparisons
    range <- c(min(data$Y_centered), max(data$Y_centered))
    
    # Within-generation predator effects
    results <- rbind(results, data.frame(
      Type = type,
      Comparison = "G1 Herb vs G1 Pred",
      BA_value = calculate_bhattacharyya_affinity(g1h_density, g1p_density, range)
    ))
    
    # Cross-generation comparisons within treatments
    results <- rbind(results, data.frame(
      Type = type,
      Comparison = "G1 Herb vs G2 Herb",
      BA_value = calculate_bhattacharyya_affinity(g1h_density, g2h_density, range)
    ))
    
    results <- rbind(results, data.frame(
      Type = type,
      Comparison = "G1 Pred vs G2 Pred",
      BA_value = calculate_bhattacharyya_affinity(g1p_density, g2p_density, range)
    ))
  }
  
  return(results)
}

# Calculate and display results
ba_results <- calculate_ba_results(behavior_data)

# Print results in a formatted way
print("Bhattacharyya's Affinity Results:")
print(knitr::kable(ba_results, digits = 3))




