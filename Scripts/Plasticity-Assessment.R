# F1 and F2 Plasticity Comparison
# Code by Baker and Sommer

# Packages 
library(dplyr)
library(lme4)
library(lubridate)
library(sp)
library(ggplot2)
library(tidyr)

# Behavior ----

Behavior_2023 <- read.csv("Data/Behavior_Raw.csv")
# Raw X, Y data in given in 4cm x 4cm grid coordinates. Converting to cm
Behavior_2023$Y <- Behavior_2023$Y*4
Behavior_2023$X <- Behavior_2023$X*4
Behavior_2023$Distance <- Behavior_2023$Distance*4

# Adjusted start time
Behavior_2023$Time <- hms(Behavior_2023$Time)
Behavior_2023 <- subset(Behavior_2023,Time>(hms("8:40:00")))

behavior_data <- Behavior_2023 %>% 
  filter(!(Population %in% c("SP", "HF", "MC"))) %>%  # Exclude specific populations
  filter(Predation.Treatment == "Predator")   # Keep only "Predator" treatment


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

# Add temperature to behavior_data
behavior_data <- behavior_data %>%
  rowwise() %>%
  mutate(Nearest_Temperature = find_nearest_temp(DateTime, behavior_temp))

# Add a unique id
behavior_data <- behavior_data %>%
  mutate(id = factor(paste(Population, Generation, Day, Terraria, sep = "_")))


# Calculate 1D kernel density estimate for each id and determine the range for core and broad levels
density_estimate_ranges <- behavior_data %>%
  group_by(id) %>%
  do({
    density_estimate <- density(.$Y, bw = "nrd0")  # Calculate density for Y within each id
    cum_density <- cumsum(density_estimate$y) / sum(density_estimate$y)
    
    # Determine 50% and 95% isopleths
    core_threshold <- 0.50
    broad_threshold <- 0.95
    
    # Find the range for core and broad levels
    core_range <- range(density_estimate$x[cum_density <= core_threshold])
    broad_range <- range(density_estimate$x[cum_density <= broad_threshold])
    
    # Determine Generation based on id
    generation <- ifelse(grepl("F1", .$id), "F1", "F2")
    
    data.frame(
      core_low = core_range[1],
      core_high = core_range[2],
      broad_low = broad_range[1],
      broad_high = broad_range[2],
      Generation = generation
    )
  })



# Calculate 1D kernel density estimate for each id
density_estimates <- behavior_data %>%
  group_by(id) %>%
  do({
    density_estimate <- density(.$Y, bw = "nrd0")  # Calculate density for Y within each id
    cum_density <- cumsum(density_estimate$y) / sum(density_estimate$y)
    
    # Determine 50% and 95% isopleths
    core_threshold <- 0.50
    broad_threshold <- 0.95
    
    core_level <- density_estimate$x[which(cum_density >= core_threshold)[1]]
    broad_level <- density_estimate$x[which(cum_density >= broad_threshold)[1]]
    
    data.frame(core_level = core_level, broad_level = broad_level)
  })

# Function to calculate isopleths for each generation
calculate_isopleths <- function(y) {
  density_estimate <- density(y[y > 0], bw = 5)  # Restrict to Y > 0
  cum_density <- cumsum(density_estimate$y) / sum(density_estimate$y)
  
  core_threshold <- 0.50
  broad_threshold <- 0.95
  
  core_range <- range(density_estimate$x[cum_density <= core_threshold])
  broad_range <- range(density_estimate$x[cum_density <= broad_threshold])
  
  return(list(core_range = core_range, broad_range = broad_range))
}

# Calculate isopleths for each generation
isopleths <- behavior_data %>%
  group_by(Generation) %>%
  summarize(core_range = list(calculate_isopleths(Y)$core_range),
            broad_range = list(calculate_isopleths(Y)$broad_range))

# Restructure the data for plotting
library(tidyr)

# Convert to long format
density_long <- density_estimate_ranges %>%
  pivot_longer(cols = c(core_low, core_high, broad_low, broad_high),
               names_to = "range_type",
               values_to = "height") %>%
  mutate(type = ifelse(grepl("core", range_type), "core", "broad"))

# Calculate average ranges for each generation
average_ranges <- density_estimate_ranges %>%
  group_by(Generation) %>%
  summarize(
    avg_core_low = mean(core_low),
    avg_core_high = mean(core_high),
    avg_broad_low = mean(broad_low),
    avg_broad_high = mean(broad_high)
  )

# Adjust the ranges to ensure no negative values
average_ranges <- average_ranges %>%
  mutate(
    avg_core_low = pmax(avg_core_low, 0),
    avg_broad_low = pmax(avg_broad_low, 0)
  )

#### Range plots ----

library(ggforce)

# Ensure Generation is a factor with levels F1 and F2
density_estimate_ranges$Generation <- factor(density_estimate_ranges$Generation, levels = c("F1", "F2"))
average_ranges$Generation <- factor(average_ranges$Generation, levels = c("F1", "F2"))



#### Bhattacharyya's affinity ----

library(ks)  # For kernel density estimation

# Function to calculate Bhattacharyya's affinity for a specific range
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
common_grid <- seq(min(behavior_data$Y), max(behavior_data$Y), length.out = 100)

density_F1_Y <- kde(x = behavior_data$Y[behavior_data$Generation == "F1"], eval.points = common_grid)
density_F2_Y <- kde(x = behavior_data$Y[behavior_data$Generation == "F2"], eval.points = common_grid)

##### 50% isopleth -----

# Calculate the 50% isopleth range for F1
cum_density_F1 <- cumsum(density_F1_Y$estimate) / sum(density_F1_Y$estimate)
core_range_F1 <- range(density_F1_Y$eval.points[cum_density_F1 <= 0.5])

# Calculate the 50% isopleth range for F2
cum_density_F2 <- cumsum(density_F2_Y$estimate) / sum(density_F2_Y$estimate)
core_range_F2 <- range(density_F2_Y$eval.points[cum_density_F2 <= 0.5])

# Use the intersection of the core ranges for both generations
core_range <- c(max(core_range_F1[1], core_range_F2[1]), min(core_range_F1[2], core_range_F2[2]))

# Ensure the core range is valid
if (core_range[1] >= core_range[2]) {
  stop("The core ranges do not overlap.")
}

BA_Y_core <- calculate_bhattacharyya_affinity_isopleth(density_F1_Y, density_F2_Y, core_range)
print(paste("Bhattacharyya's Affinity for Y data (50% isopleth):", BA_Y_core))


##### 95% isopleth -----

# Calculate the 95% isopleth range for F1
broad_range_F1 <- range(density_F1_Y$eval.points[cum_density_F1 <= 0.95])

# Calculate the 95% isopleth range for F2
broad_range_F2 <- range(density_F2_Y$eval.points[cum_density_F2 <= 0.95])

# Use the intersection of the broad ranges for both generations
broad_range <- c(max(broad_range_F1[1], broad_range_F2[1]), min(broad_range_F1[2], broad_range_F2[2]))

# Ensure the broad range is valid
if (broad_range[1] >= broad_range[2]) {
  stop("The broad ranges do not overlap.")
}

# Calculate Bhattacharyya's affinity for the 95% isopleth
BA_Y_broad <- calculate_bhattacharyya_affinity_isopleth(density_F1_Y, density_F2_Y, broad_range)
print(paste("Bhattacharyya's Affinity for Y data (95% isopleth):", BA_Y_broad))


# Adjust the ranges to ensure no negative values
average_ranges <- average_ranges %>%
  mutate(
    avg_core_low = pmax(avg_core_low, 0),
    avg_broad_low = pmax(avg_broad_low, 0)
  )

# Calculate mean height for each generation directly from behavior_data
mean_heights <- behavior_data %>%
  group_by(Generation) %>%
  summarize(median_y = median(Y, na.rm = TRUE))

# Ensure Generation is a factor with levels F1 and F2
mean_heights$Generation <- factor(mean_heights$Generation, levels = c("F1", "F2"))

# Average Ranges Figure with Mean Points
ggplot(average_ranges) +

  geom_ellipse(aes(x0 = as.numeric(Generation), 
                   y0 = (avg_broad_low + avg_broad_high) / 2,
                   a = 0.2,  # Same width as core ellipse
                   b = (avg_broad_high - avg_broad_low) / 2, angle = 0, fill = "broad"), 
               color = "black", alpha = 0.5) +
  geom_ellipse(aes(x0 = as.numeric(Generation), 
                   y0 = (avg_core_low + avg_core_high) / 2,
                   a = 0.15,  # Fixed width for both ellipses
                   b = (avg_core_high - avg_core_low) / 2, angle = 0, fill = "core"), 
               color = "black", alpha = 0.5) +
  geom_point(data = mean_heights, aes(x = as.numeric(Generation), y = median_y), 
             color = "red", size = 3) +
  geom_line(data = mean_heights, aes(x = as.numeric(Generation), y = median_y), 
            color = "red", linewidth = 1) +
  scale_fill_manual(values = c("core" = "black", "broad" = "lightgray")) +
  labs(title = "Core (50%) and Broad Ranges (95%) by Generation", x = "Generation", y = "Canopy Height") +
  theme_minimal() +
  scale_x_continuous(breaks = 1:2, labels = c("F1", "F2"))

#### Additional figures ----


# Individual Ranges --> looks ridiculous
ggplot(density_estimate_ranges) +
  geom_ellipse(aes(x0 = as.numeric(Generation), y0 = (core_low + core_high) / 2,
                   a = 0.2, b = (core_high - core_low) / 2, angle = 0, fill = "core"), alpha = 0.3) +
  labs(title = "Individual Core and Broad Ranges by Generation", x = "Generation", y = "Height") +
  theme_minimal() +
  scale_x_continuous(breaks = 1:2, labels = c("F1", "F2"))

### Average height by generation

# Calculate average height (Y) for each id
average_height_data <- behavior_data %>%
  group_by(Generation, id) %>%
  summarize(mean_Y = mean(Y, na.rm = TRUE))

t.test(mean_Y ~ Generation, data = average_height_data)


### Temperature by height

behavior_data$Nearest_Temperature <- as.numeric(behavior_data$Nearest_Temperature)

ggplot(behavior_data, aes(x = Nearest_Temperature, y = Y, color = Generation)) +
  geom_jitter(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Relationship between Height (Y) and Nearest Temperature by Generation",
       x = "Nearest Temperature (°C)", y = "Height (Y)") +
  theme_minimal() +
  scale_color_manual(values = c("F1" = "blue", "F2" = "red")) +
  facet_wrap(~Population)


# Respiration ----

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

# Logger Pro Automatically Calculates the integral for respiration (Logger_Pro_Calcd) 
# MB also calculated the values manually (Manual_Resp_Calc)

Respriation_Mass_2023<-read.csv("Data/Respriation_Mass_2023_F1_F2.csv")
Respriation_Mass_2023<-arrange(Respriation_Mass_2023, Run,Chamber)
Respriation_Mass_2023_times2<-rbind(Respriation_Mass_2023,Respriation_Mass_2023)
respiration_data <-cbind(Final_output_all_temps,Respriation_Mass_2023_times2)

### Analysis ----

# Make a "reaction norm" plot of resp rate across the two temps, for the F1 and F2
# See if there is any difference in slope


