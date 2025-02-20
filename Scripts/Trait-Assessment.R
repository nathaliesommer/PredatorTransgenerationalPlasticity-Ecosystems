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

#### LME repeated measures ---
height_model_rm <- lmer(Y_centered ~ Treatment + (1|id), 
                        data = behavior_data)

summary(height_model_rm)

# emmeans post-hoc contrasts
emmeans(height_model_rm, pairwise ~ Treatment)

#### Distribution test ---

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

## Figures ----
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