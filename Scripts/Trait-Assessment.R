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

Behavior_2023 <- read.csv("Data/Behavior_Pred-G1G2_Herb-G2.csv")
Behavior_2023 <- Behavior_2023 %>%
  filter(Predation.Treatment == "Predator")
Behavior_2023$Time <- hms(Behavior_2023$Time)
behavior_data <- subset(Behavior_2023,Time>(hms("8:40:00")))
behavior_temp <- read.csv("Data/Behavior_Temps.csv")
behavior_temp <- behavior_temp %>%
  mutate(DateTime = mdy_hm(Date.Time..EDT.))
behavior_data <- behavior_data %>%
  mutate(DateTime = mdy(Date) + hours(hour(Time)) + minutes(minute(Time)))
find_nearest_temp <- function(obs_time, temp_data) {
  temp_data %>%
    filter(abs(difftime(DateTime, obs_time, units = "mins")) == min(abs(difftime(DateTime, obs_time, units = "mins")))) %>%
    slice(1) %>%
    pull(Ch..1...Temperature.....C.)
}
behavior_data <- behavior_data %>%
  rowwise() %>%
  mutate(Nearest_Temperature = find_nearest_temp(DateTime, behavior_temp))
behavior_data <- behavior_data %>%
  mutate(Treatment = paste(Generation)) %>%
  mutate(Treatment = factor(Treatment, 
                            levels = c("G1", "G2")))
height_model <- lmer(Y ~ Treatment + (1|Population), 
                     data = behavior_data)
VarCorr(height_model)
behavior_data <- behavior_data %>%
  group_by(Population) %>%
  mutate(Y_centered = scale(Y, scale = FALSE),
         Y_centered = Y_centered - min(Y_centered)) %>%
  ungroup()
print(range(behavior_data$Y_centered, na.rm = TRUE))
print(sum(!is.finite(behavior_data$Y_centered)))
print(sum(is.na(behavior_data$Y_centered)))
print(summary(behavior_data$Y_centered))
print(summary(behavior_data$Y_centered[behavior_data$Generation == "G1"]))
print(summary(behavior_data$Y_centered[behavior_data$Generation == "G2"]))
print("Input data summary:")
print(summary(behavior_data$Y_centered))
print(table(behavior_data$Treatment))
calculate_density_ranges <- function(data, min_obs = 10) {
  data %>%
    group_by(Treatment) %>%
    filter(n() >= min_obs) %>%
    summarize(
      avg_core_low = quantile(Y_centered, probs = 0.25, na.rm = TRUE),
      avg_core_high = quantile(Y_centered, probs = 0.75, na.rm = TRUE),
      avg_broad_low = quantile(Y_centered, probs = 0.025, na.rm = TRUE),
      avg_broad_high = quantile(Y_centered, probs = 0.975, na.rm = TRUE),
      Y_centered = mean(Y_centered, na.rm = TRUE),
      .groups = 'drop'
    )
}
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
density_estimate_ranges <- calculate_density_ranges(behavior_data)
print("\nDensity Range Diagnostics:")
print_diagnostic(density_estimate_ranges)

### Analysis ----

behavior_data <- behavior_data %>%
  mutate(
    id = paste(Population, Predation.Treatment, Generation, Cage.ID, Day, sep = "_"))
height_model_rm <- lmer(Y_centered ~ Treatment + (1|id), 
                        data = behavior_data)
summary(height_model_rm)
emmeans(height_model_rm, pairwise ~ Treatment)
wilcox_test <- wilcox.test(
  behavior_data$Y_centered[behavior_data$Treatment == "G1"],
  behavior_data$Y_centered[behavior_data$Treatment == "G2"]
)
print(wilcox_test)

## Figures ----

custom_colors <- c(
  "G1" = "#C69648",
  "G2" = "#A15A29"
)
treatment_labels <- c(
  "G1" = "G[1]",
  "G2" = "G[2]"
)
create_aggregated_plot <- function(data) {
  means <- data %>%
    group_by(Generation) %>%
    summarise(
      Y_centered = mean(Y_centered),
      .groups = 'drop'
    )
  ggplot(data, aes(x = Generation, y = Y_centered)) +
    geom_violin(
      aes(fill = Generation),
      alpha = 0.5,
      trim = TRUE
    ) +
    geom_jitter(
      aes(color = Generation),
      alpha = 0.2,
      size = 0.8,
      width = 0.2
    ) +
    geom_point(
      data = means,
      size = 3,
      color = "black"
    ) +
    geom_text(
      data = means,
      aes(label = round(Y_centered, 1)),
      vjust = -1,
      size = 3.5
    ) +
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
aggregated_plot <- create_aggregated_plot(behavior_data)
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
aggregated_plot
temperature_plot

# Respiration ----

# Packages
library(dplyr)
library(lme4)
library(ggplot2)

### Data Processing ----
focal_pops <- c("FN", "DC", "UP", "SC")
Respiration_Full<-read.csv("Data/Respriation_2023_F1_F2_Raw.csv")
Respiration_Full<-Respiration_Full[,-1]
Tempcatlist <- data.frame(matrix(ncol = 1, nrow = 2))
Tempcatlist[,1]<-c("25C","30C")
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
Respiration_Mass_2023<-read.csv("Data/Respriation_Mass_2023_F1_F2.csv")
Respiration_Mass_2023<-arrange(Respiration_Mass_2023, Run,Chamber)
Respiration_Mass_2023_times2<-rbind(Respiration_Mass_2023,Respiration_Mass_2023)
respiration_data <-cbind(Final_output_all_temps,Respiration_Mass_2023_times2)
respiration_data <- respiration_data[, !duplicated(names(respiration_data))]
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
  filter(Population %in% focal_pops)
respiration_data <- respiration_data %>%
  dplyr::mutate(SMR = Logger_Pro_Calcd / Mass)
respiration_data_predator <- respiration_data %>% 
  dplyr::filter(Predation == "Predator")
respiration_data_full <- respiration_data
respiration_data_predator <- respiration_data_predator %>%
  mutate(Temp_centered = case_when(
    Temp == "25C" ~ -2.5,
    Temp == "30C" ~ 2.5
  ))
respiration_25C <- respiration_data_predator %>%
  filter(Temp == "25C")
aov_25C <- aov(SMR ~ Generation, data = respiration_25C)
summary(aov_25C)
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
respiration_30C <- respiration_data_predator %>%
  filter(Temp == "30C")
aov_30C <- aov(SMR ~ Generation, data = respiration_30C)
summary(aov_30C)
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
g1p_data <- respiration_data_predator %>% 
  filter(Generation == "G1")
g1p_model <- lmer(SMR ~ Temp_centered + (1 | Individual), data = g1p_data)
summary(g1p_model)
g2p_data <- respiration_data_predator %>% 
  filter(Generation == "G2")
g2p_model <- lmer(SMR ~ Temp_centered + (1 | Individual), data = g2p_data)
summary(g2p_model)
g1g2p_model <- lmer(SMR ~ Temp_centered + Generation + (1 | Individual), data = respiration_data_predator)
summary(g1g2p_model)
emmeans(g1g2p_model, pairwise ~ Generation)
## Figures ----
reaction_norm_plot <- ggplot(respiration_data_predator, aes(x = Temp, y = SMR, group = Generation, color = Generation)) +
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
library(patchwork)
combined_plot <- reaction_norm_plot + behavior_plot +
  plot_layout(guides = "collect") &  
  theme(legend.position = "bottom")
combined_plot

# ---- Supporting Information ----
## Traits ----
behavior_data_supp <- read.csv("Data/Behavior_Pred-G1G2_Herb-G2.csv")
behavior_data_supp <- behavior_data_supp %>%
  mutate(Treatment = case_when(
    Predation.Treatment == "Predator" & Generation == "G2" ~ "G2_Predator",
    Predation.Treatment == "Herbivore" & Generation == "G2" ~ "G2_Herbivore"
  )) %>%
  filter(!is.na(Treatment)) 
behavior_data_supp <- behavior_data_supp %>%
  mutate(id = paste(Population, Predation.Treatment, Generation, Cage.ID, Day, sep = "_"))
behavior_data_supp <- behavior_data_supp %>%
  group_by(Population, Predation.Treatment) %>%
  mutate(Y_centered = scale(Y, scale = FALSE),
         Y_centered = Y_centered - min(Y_centered, na.rm = TRUE)) %>%
  ungroup()
behavior_data_supp$Treatment <- factor(behavior_data_supp$Treatment, levels = c("G2_Predator", "G2_Herbivore"))
supp_colors <- c(
  "G2_Predator" = "#A15A29",  
  "G2_Herbivore" = "#9BA48C"  
)
supp_density_estimate_ranges <- calculate_density_ranges(behavior_data_supp)
print("\nSupplemental Density Range Diagnostics:")
print_diagnostic(supp_density_estimate_ranges)
supp_height_model_rm <- lmer(Y_centered ~ Treatment + (1|id), 
                            data = behavior_data_supp)
summary(supp_height_model_rm)
supp_emmeans <- emmeans(supp_height_model_rm, pairwise ~ Treatment)
print(supp_emmeans)
supp_wilcox_test <- wilcox.test(
  behavior_data_supp$Y_centered[behavior_data_supp$Treatment == "G2_Predator"],
  behavior_data_supp$Y_centered[behavior_data_supp$Treatment == "G2_Herbivore"]
)
print(supp_wilcox_test)
create_supp_aggregated_plot <- function(data) {
  y_range <- range(data$Y_centered, na.rm = TRUE)
  means <- data %>%
    group_by(Treatment) %>%
    summarise(
      Y_centered = mean(Y_centered, na.rm = TRUE),
      .groups = 'drop'
    )
  ggplot(data, aes(x = Treatment, y = Y_centered)) +
    geom_violin(
      aes(fill = Treatment),
      alpha = 0.5,
      trim = TRUE
    ) +
    geom_jitter(
      aes(color = Treatment),
      alpha = 0.2,
      size = 0.8,
      width = 0.2,
      height = 1.5
    ) + 
    geom_point(
      data = means,
      size = 3,
      color = "black"
    ) +
    geom_text(
      data = means,
      aes(label = round(Y_centered, 1)),
      vjust = -1,
      size = 3.5
    ) +
    scale_y_continuous(limits = c(0, 80)) +
    scale_x_discrete(
      labels = c(
        "G2 Predator",
        "G2 Herbivore"
      )
    ) +
    scale_fill_manual(values = supp_colors) +
    scale_color_manual(values = supp_colors) +
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
supp_behavior_plot <- create_supp_aggregated_plot(behavior_data_supp)
g2_pred_resp <- respiration_data_full %>%
  filter(Predation == "Predator", Generation == "G2", Population %in% focal_pops) %>%
  mutate(Treatment = "G2_Predator")
g2_herb_resp <- respiration_data_full %>%
  filter(Predation == "Herbivore", Generation == "G2", Population %in% focal_pops) %>%
  mutate(Treatment = "G2_Herbivore")
supp_resp_data <- bind_rows(
  g2_pred_resp,
  g2_herb_resp
)
supp_resp_data$Treatment <- factor(supp_resp_data$Treatment, levels = c("G2_Predator", "G2_Herbivore"))
supp_resp_data <- supp_resp_data %>%
  mutate(Temp_centered = case_when(
    Temp == "25C" ~ -2.5,
    Temp == "30C" ~ 2.5
  ))
g2p2h_supp_model <- lmer(SMR ~ Temp_centered + Treatment + (1 | Individual), data = supp_resp_data)
  
print(summary(g2p2h_supp_model))
supp_emmeans_resp <- emmeans(g2p2h_supp_model, pairwise ~ Treatment)
print(supp_emmeans_resp)
supp_resp_reaction_norm <- ggplot(supp_resp_data, aes(x = Temp, y = SMR, group = Treatment, color = Treatment)) +
  geom_point(alpha = 0.2, size = 1) +
  geom_line(aes(group = Individual), alpha = 0.2) +
  stat_summary(fun = mean, geom = "line", linewidth = 1.5) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  scale_color_manual(
    values = supp_colors,
    labels = c("G2 Predator", "G2 Herbivore")
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
  )
supp_combined_plot <- supp_resp_reaction_norm + supp_behavior_plot +
  plot_layout(guides = "collect") &  # Collect legends together
  theme(legend.position = "bottom")
print(supp_combined_plot)

