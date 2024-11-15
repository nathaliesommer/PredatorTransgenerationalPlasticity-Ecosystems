# CarbonNitrogen ----

library(tidyverse)
library(ggforce)
library(RColorBrewer)
library(paletteer)

CNdata <- read.csv("Data/Elements/CN_Tidy.csv")

CNdata <- CNdata %>% 
  mutate(
    PercentN = ifelse(Note %in% c("None", "Defoliated"), 0, PercentN),
    PercentC = ifelse(Note %in% c("None", "Defoliated"), 0, PercentC)
  )

CNdata <- CNdata %>%
  rename(Sample_ID = CageID) %>%
  dplyr::select(Sample_ID, Year, Population, Site, PercentN, PercentC, SampleType) %>%
  pivot_wider(names_from = SampleType, values_from = c(PercentN, PercentC))




## Plots for Matthew (works from CNdata load) ----

means_data <- CNdata %>% 
  filter(Year == 2021) %>% 
  filter(SampleLevel == "Cage") %>% 
  filter(Site %in% c("YF", "UP")) %>% 
  filter(SampleType %in% c("SORU", "POPRC")) %>% 
  group_by(Site, SampleType) %>% 
  summarise(mean_PercentN = mean(PercentN, na.rm = TRUE)) %>% 
  ungroup()

CNdata %>% 
  filter(Year == 2021) %>% 
  filter(SampleLevel == "Cage") %>% 
  filter(Site %in% c("YF", "UP")) %>% 
  filter(SampleType %in% c("SORU", "POPRC")) %>% 
  ggplot(aes(x = Site, y = PercentN, fill = Site)) + 
  geom_boxplot(outliers = FALSE) +
  geom_jitter(size = 1, width = .1) +
  theme_classic() +
  scale_fill_paletteer_d("ggsci::alternating_igv") +
  labs(title = "Plant Types %N 2021",
       x = "Type",
       y = "%N") +
  facet_wrap(~SampleType) +  
  geom_text(data = means_data, aes(x = Site, y = 0.5, label = round(mean_PercentN, 1)), 
                                       position = position_dodge(width = 0.75), size = 3, color = "black")



means_data <- CNdata %>% 
  filter(Year == 2021) %>% 
  filter(SampleLevel == "Cage") %>% 
  filter(Site %in% c("YF", "UP")) %>% 
  filter(SampleType %in% c("SORU", "POPRC")) %>% 
  group_by(Site, SampleType) %>% 
  summarise(mean_PercentC = mean(PercentC, na.rm = TRUE)) %>% 
  ungroup()

CNdata %>% 
  filter(Year == 2021) %>% 
  filter(SampleLevel == "Cage") %>% 
  filter(Site %in% c("YF", "UP")) %>% 
  filter(SampleType %in% c("SORU", "POPRC")) %>% 
  filter(PercentC < 60) %>%
  filter(PercentC > 20) %>% 
  ggplot(aes(x = Site, y = PercentC, fill = Site)) + 
  geom_boxplot(outliers = FALSE) +
  geom_jitter(size = 1, width = .1) +
  theme_classic() +
  scale_fill_paletteer_d("ggsci::alternating_igv") +
  labs(title = "Plant Types %C 2021",
       x = "Type",
       y = "%C") +
  facet_wrap(~SampleType) +
  geom_text(data = means_data, aes(x = Site, y = 35, label = round(mean_PercentC, 1)), 
            position = position_dodge(width = 0.75), size = 3, color = "black")

means_data <- CNdata %>% 
  filter(Year == 2021) %>% 
  filter(SampleLevel == "Cage") %>% 
  filter(Site %in% c("YF", "UP")) %>% 
  filter(SampleType %in% c("SORU", "POPRC")) %>% 
  group_by(Site, SampleType) %>% 
  summarise(mean_RatioCN = mean(RatioCN, na.rm = TRUE)) %>% 
  ungroup()

CNdata %>% 
  filter(Year == 2021) %>% 
  filter(SampleLevel == "Cage") %>% 
  filter(Site %in% c("YF", "UP")) %>% 
  filter(SampleType %in% c("SORU", "POPRC")) %>% 
  filter(PercentC < 60) %>%
  filter(PercentC > 20) %>% 
  ggplot(aes(x = Site, y = RatioCN, fill = Site)) + 
  geom_boxplot(outliers = FALSE) +
  geom_jitter(size = 1, width = .1) +
  theme_classic() +
  scale_fill_paletteer_d("ggsci::alternating_igv") +
  labs(title = "Plant Types RatioCN 2021",
       x = "Type",
       y = "RatioCN") +
  facet_wrap(~SampleType) +
  facet_wrap(~SampleType) +
  geom_text(data = means_data, aes(x = Site, y = 7, label = round(mean_RatioCN, 1)), 
            position = position_dodge(width = 0.75), size = 3, color = "black")



means_data <- CNdata %>% 
  filter(Year == 2023) %>% 
  filter(SampleLevel == "Cage") %>% 
  filter(Site %in% c("YF", "UP")) %>% 
  filter(SampleType %in% c("SORU", "POPRC")) %>% 
  group_by(Site, CageTreatment, SampleType) %>% 
  summarise(mean_PercentC = mean(PercentC, na.rm = TRUE)) %>% 
  ungroup()


CNdata %>% 
  filter(Year == 2023) %>% 
  filter(SampleLevel == "Cage") %>% 
  filter(Site %in% c("YF", "UP")) %>% 
  filter(SampleType %in% c("SORU", "POPRC")) %>% 
  filter(PercentC < 80) %>% 
  mutate(CageTreatment = factor(CageTreatment, levels = c("Vegetation", "Herbivore"))) %>%  
  ggplot(aes(x = Site, y = PercentC, fill = CageTreatment)) + 
  geom_boxplot(outliers = FALSE) +
  scale_fill_paletteer_d("ggthemes::wsj_red_green") +
  geom_jitter(size = 1, width = .1) +
  theme_classic() +
  labs(title = "Plant Types %C by Treatment, 2023",
       x = "Type",
       y = "%C") +
  facet_wrap(~SampleType) +
  geom_text(data = means_data, aes(x = Site, y = 37, label = round(mean_PercentC, 1)), 
            position = position_dodge(width = 0.75), size = 3, color = "black")





means_data <- CNdata %>% 
  filter(Year == 2023) %>% 
  filter(SampleLevel == "Cage") %>% 
  filter(Site %in% c("YF", "UP")) %>% 
  filter(SampleType %in% c("SORU", "POPRC")) %>% 
  group_by(Site, CageTreatment, SampleType) %>% 
  summarise(mean_PercentN = mean(PercentN, na.rm = TRUE)) %>% 
  ungroup()


CNdata %>% 
  filter(Year == 2023) %>% 
  filter(SampleLevel == "Cage") %>% 
  filter(Site %in% c("YF", "UP")) %>% 
  filter(SampleType %in% c("SORU", "POPRC")) %>% 
  mutate(CageTreatment = factor(CageTreatment, levels = c("Vegetation", "Herbivore"))) %>%  
  ggplot(aes(x = Site, y = PercentN, fill = CageTreatment)) + 
  geom_boxplot(outliers = FALSE) +
  scale_fill_paletteer_d("ggthemes::wsj_red_green") +
  geom_jitter(size = 1, width = .1) +
  theme_classic() +
  labs(title = "Plant Types %N by Treatment, 2023",
       x = "Type",
       y = "%N") +
  facet_wrap(~SampleType) + 
  geom_text(data = means_data, aes(x = Site, y = .4, label = round(mean_PercentN, 1)), 
            position = position_dodge(width = 0.75), size = 3, color = "black")



means_data <- CNdata %>% 
  filter(Year == 2023) %>% 
  filter(SampleLevel == "Cage") %>% 
  filter(Site %in% c("YF", "UP")) %>% 
  filter(SampleType %in% c("SORU", "POPRC")) %>% 
  filter(RatioCN < 60) %>% 
  group_by(Site, CageTreatment, SampleType) %>% 
  summarise(mean_RatioCN = mean(RatioCN, na.rm = TRUE)) %>% 
  ungroup()

CNdata %>% 
  filter(Year == 2023) %>% 
  filter(SampleLevel == "Cage") %>% 
  filter(Site %in% c("YF", "UP")) %>% 
  filter(SampleType %in% c("SORU", "POPRC")) %>% 
  filter(RatioCN < 60) %>% 
  mutate(CageTreatment = factor(CageTreatment, levels = c("Vegetation", "Herbivore"))) %>%  
  ggplot(aes(x = Site, y = RatioCN, fill = CageTreatment)) + 
  geom_boxplot(outliers = FALSE) +
  scale_fill_paletteer_d("ggthemes::wsj_red_green") +
  geom_jitter(size = 1, width = .1) +
  theme_classic() +
  labs(title = "Plant Types RatioCN by Treatment, 2023",
       x = "Type",
       y = "%RatioCN") +
  facet_wrap(~SampleType) +
  geom_text(data = means_data, aes(x = Site, y = 7, label = round(mean_RatioCN, 1)), 
            position = position_dodge(width = 0.75), size = 3, color = "black")




CNdata %>% 
  filter(Year == 2021) %>% 
  filter(SampleLevel == "Cage") %>% 
  filter(Site %in% c("YF", "UP")) %>% 
  filter(SampleType %in% c("SORU", "POPRC")) %>% 
  group_by(Site) %>% 
  summarise(mean_PercentN = mean(PercentN, na.rm = TRUE))

CNdata %>% 
  filter(Year == 2021) %>% 
  filter(SampleLevel == "Cage") %>% 
  filter(Site %in% c("YF", "UP")) %>% 
  filter(SampleType %in% c("SORU", "POPRC")) %>% 
  group_by(Site, SampleType) %>% 
  summarise(mean_PercentC = mean(PercentC, na.rm = TRUE))

CNdata %>% 
  filter(Year == 2021) %>% 
  filter(SampleLevel == "Cage") %>% 
  filter(Site %in% c("YF", "UP")) %>% 
  filter(SampleType == "Soil") %>% 
  mean(PercentN)






## DataViz ----

### Soil ----
CNdata %>%
  filter(SampleLevel == "Cage", SampleType == "SOIL") %>% 
  ggplot(aes(x = Site, y = RatioCN)) +
  geom_point(size = 3) +
  labs(title = "Soil RatioCN Across Sites",
       x = "Site",
       y = "RatioCN") +
  theme_classic() + 
  coord_flip()

### SORU ----
CNdata %>%
  filter(SampleLevel == "Site", SampleType == "SORU") %>% 
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF", "DC", "SP", "MC", "UP"))) %>% 
  ggplot(aes(x = Site, y = RatioCN)) +
  geom_point(size = 3) +
  labs(title = "SORU RatioCN Across Sites",
       x = "Site",
       y = "RatioCN") +
  theme_classic() + 
  coord_flip()

### POPRC ----
CNdata %>%
  filter(SampleLevel == "Site", SampleType == "POPRC") %>% 
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF", "DC", "SP", "MC", "UP"))) %>% 
  ggplot(aes(x = Site, y = RatioCN)) +
  geom_point(size = 3) +
  labs(title = "POPRC RatioCN Across Sites",
       x = "Site",
       y = "RatioCN") +
  theme_classic() + 
  coord_flip()

### MISC ----
CNdata %>%
  filter(SampleLevel == "Site", SampleType == "MISC") %>% 
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF", "DC", "SP", "MC", "UP"))) %>% 
  ggplot(aes(x = Site, y = RatioCN)) +
  geom_point(size = 3) +
  labs(title = "MISC RatioCN Across Sites",
       x = "Site",
       y = "RatioCN") +
  theme_classic() + 
  coord_flip()

# Variability (2021) ----
# Within versus across site variation in %C and %N, for soil and plants 

## Soil ----

### %C ----
CNdata %>%
  filter(SampleType == "SOIL") %>% 
  filter(!Site %in% c("SP", "HF")) %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "DC", "MC", "UP"))) %>%
  bind_rows(
    CNdata %>%
      filter(SampleLevel == "Cage", SampleType == "SOIL") %>%
      mutate(Site = "All Sites")
  ) %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF", "DC", "SP", "MC", "UP", "All Sites"))) %>%
  ggplot(aes(x = Site, y = PercentC, fill = Site)) +
  geom_violin(trim = FALSE, alpha = 0.7, width = 0.5) + 
  geom_point(position = position_jitter(width = 0.2, height = 0), size = 1, alpha = 0.5) +
  labs(title = "Soil %C Variation: Within and Across Sites",
       x = "Site",
       y = "PercentC") +
  theme_classic() +
  coord_flip()

### %N ----
CNdata %>%
  filter(SampleType == "SOIL") %>% 
  filter(!Site %in% c("SP", "HF")) %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "DC", "MC", "UP"))) %>%
  bind_rows(
    CNdata %>%
      filter(SampleLevel == "Cage", SampleType == "SOIL") %>%
      mutate(Site = "All Sites")
  ) %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF", "DC", "SP", "MC", "UP", "All Sites"))) %>%
  ggplot(aes(x = Site, y = PercentN)) +
  geom_violin(aes(fill = Site), trim = FALSE, alpha = 0.5) +
  labs(title = "Soil %N Variation: Within and Across Sites",
       x = "Site",
       y = "PercentN") +
  theme_classic() +
  coord_flip()

## SORU ----

### %C ----
CNdata %>%
  filter(SampleType == "SORU") %>% 
  filter(!Site %in% c("SP", "HF")) %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "DC", "MC", "UP"))) %>%
  bind_rows(
    CNdata %>%
      filter(SampleLevel == "Cage", SampleType == "SORU") %>%
      mutate(Site = "All Sites")
  ) %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF", "DC", "SP", "MC", "UP", "All Sites"))) %>%
  ggplot(aes(x = Site, y = PercentC)) +
  geom_violin(aes(fill = Site), trim = FALSE, alpha = 0.5) +
  labs(title = "SORU %C Variation: Within and Across Sites",
       x = "Site",
       y = "PercentC") +
  theme_classic() +
  coord_flip()

### %N ----
CNdata %>%
  filter(SampleType == "SORU") %>% 
  filter(!Site %in% c("SP", "HF")) %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "DC", "MC", "UP"))) %>%
  bind_rows(
    CNdata %>%
      filter(SampleLevel == "Cage", SampleType == "SORU") %>%
      mutate(Site = "All Sites")
  ) %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF", "DC", "SP", "MC", "UP", "All Sites"))) %>%
  ggplot(aes(x = Site, y = PercentN)) +
  geom_violin(aes(fill = Site), trim = FALSE, alpha = 0.5) +
  labs(title = "SORU %N Variation: Within and Across Sites",
       x = "Site",
       y = "PercentN") +
  theme_classic() +
  coord_flip()

## POPRC ----

### %C ----
CNdata %>%
  filter(SampleType == "POPRC") %>% 
  filter(!Site %in% c("SP", "HF")) %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "DC", "MC", "UP"))) %>%
  bind_rows(
    CNdata %>%
      filter(SampleLevel == "Cage", SampleType == "POPRC") %>%
      mutate(Site = "All Sites")
  ) %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF", "DC", "SP", "MC", "UP", "All Sites"))) %>%
  ggplot(aes(x = Site, y = PercentC)) +
  geom_violin(aes(fill = Site), trim = FALSE, alpha = 0.5) +
  labs(title = "POPRC %C Variation: Within and Across Sites",
       x = "Site",
       y = "PercentC") +
  theme_classic() +
  coord_flip()

### %N ----
CNdata %>%
  filter(SampleType == "POPRC") %>% 
  filter(!Site %in% c("SP", "HF")) %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "DC", "MC", "UP"))) %>%
  bind_rows(
    CNdata %>%
      filter(SampleLevel == "Cage", SampleType == "POPRC") %>%
      mutate(Site = "All Sites")
  ) %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF", "DC", "SP", "MC", "UP", "All Sites"))) %>%
  ggplot(aes(x = Site, y = PercentN)) +
  geom_violin(aes(fill = Site), trim = FALSE, alpha = 0.5) +
  labs(title = "POPRC %N Variation: Within and Across Sites",
       x = "Site",
       y = "PercentN") +
  theme_classic() +
  coord_flip()

## MISC ----

### %C ----
CNdata %>%
  filter(SampleType == "MISC") %>% 
  filter(!Site %in% c("SP", "HF")) %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "DC", "MC", "UP"))) %>%
  bind_rows(
    CNdata %>%
      filter(SampleLevel == "Cage", SampleType == "MISC") %>%
      mutate(Site = "All Sites")
  ) %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF", "DC", "SP", "MC", "UP", "All Sites"))) %>%
  ggplot(aes(x = Site, y = PercentC)) +
  geom_violin(aes(fill = Site), trim = FALSE, alpha = 0.5) +
  labs(title = "MISC %C Variation: Within and Across Sites",
       x = "Site",
       y = "PercentC") +
  theme_classic() +
  coord_flip()

### %N ----
CNdata %>%
  filter(SampleType == "MISC") %>% 
  filter(!Site %in% c("SP", "HF")) %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "DC", "MC", "UP"))) %>%
  bind_rows(
    CNdata %>%
      filter(SampleLevel == "Cage", SampleType == "MISC") %>%
      mutate(Site = "All Sites")
  ) %>%
  mutate(Site = factor(Site, levels = c("FN", "YF", "SC", "HF", "DC", "SP", "MC", "UP", "All Sites"))) %>%
  ggplot(aes(x = Site, y = PercentN)) +
  geom_violin(aes(fill = Site), trim = FALSE, alpha = 0.5) +
  labs(title = "MISC %N Variation: Within and Across Sites",
       x = "Site",
       y = "PercentN") +
  theme_classic() +
  coord_flip()


# Treatment Effects ----
# Change in %C and change in %N across treatments between 2023 and 2021

PercentC_change <- CNdata %>%
  select(CageID, SampleType, Population, Site, TransplantTreatment, CageTreatment, Year, PercentC) %>%
  group_by(CageID, SampleType) %>% 
  pivot_wider(names_from = Year, values_from = PercentC) %>%  
  mutate(ChangePercentC = `2023` - `2021`)

PercentN_change <- CNdata %>%
  select(CageID, SampleType, Population, Site, TransplantTreatment, CageTreatment, Year, PercentN) %>%
  group_by(CageID, SampleType) %>% 
  pivot_wider(names_from = Year, values_from = PercentN) %>%  
  mutate(ChangePercentN = `2023` - `2021`)

## TROPHIC TYPE ----

### Soil ----

#### %C ----
PercentC_change %>%
  filter(SampleType == "SOIL", !is.na(ChangePercentC)) %>% 
  ggplot(aes(x = CageTreatment, y = ChangePercentC, fill = CageTreatment)) +
  geom_boxplot(alpha = 0.7) + 
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  labs(title = "Change in Soil %C Across Treatment Groups",
       x = "Cage Treatment",
       y = "Change in %C") +
  theme_minimal()

#### %N ----
PercentN_change %>%
  filter(SampleType == "SOIL", !is.na(ChangePercentN)) %>% 
  ggplot(aes(x = CageTreatment, y = ChangePercentN, fill = CageTreatment)) +
  geom_boxplot(alpha = 0.7) + 
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  labs(title = "Change in Soil %N Across Treatment Groups",
       x = "Cage Treatment",
       y = "Change in %N") +
  theme_minimal()

### SORU ----
#### %C ----
PercentC_change %>%
  filter(SampleType == "SORU", !is.na(ChangePercentC)) %>% 
  ggplot(aes(x = CageTreatment, y = ChangePercentC, fill = CageTreatment)) +
  geom_boxplot(alpha = 0.7) + 
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  labs(title = "Change in SORU %C Across Treatment Groups",
       x = "Cage Treatment",
       y = "Change in %C") +
  theme_minimal()

#### %N ----
PercentN_change %>%
  filter(SampleType == "SORU", !is.na(ChangePercentN)) %>% 
  ggplot(aes(x = CageTreatment, y = ChangePercentN, fill = CageTreatment)) +
  geom_boxplot(alpha = 0.7) + 
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  labs(title = "Change in SORU %N Across Treatment Groups",
       x = "Cage Treatment",
       y = "Change in %N") +
  theme_minimal()

### POPRC ----
#### %C ----
PercentC_change %>%
  filter(SampleType == "POPRC", !is.na(ChangePercentC)) %>% 
  ggplot(aes(x = CageTreatment, y = ChangePercentC, fill = CageTreatment)) +
  geom_boxplot(alpha = 0.7) + 
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  labs(title = "Change in POPRC %C Across Treatment Groups",
       x = "Cage Treatment",
       y = "Change in PercentC") +
  theme_minimal()

#### %N ----
PercentN_change %>%
  filter(SampleType == "POPRC", !is.na(ChangePercentN)) %>% 
  ggplot(aes(x = CageTreatment, y = ChangePercentN, fill = CageTreatment)) +
  geom_boxplot(alpha = 0.7) + 
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +
  labs(title = "Change in POPRC %N Across Treatment Groups",
       x = "Cage Treatment",
       y = "Change in %N") +
  theme_minimal()


## TRANSPLANT ----
### Soil ----
#### %C ----
PercentC_change %>%
  filter(SampleType == "SOIL", !is.na(ChangePercentC)) %>%
  mutate(Population = factor(Population, levels = c("FN", "YF", "SC", "HF", "DC", "SP", "MC", "UP"))) %>%
  ggplot(aes(x = Population, y = ChangePercentC, fill = TransplantTreatment)) +
  geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) + 
  labs(title = "Change in Soil %C Across Populations and Transplant Types",
       x = "Population",
       y = "%C") +
  theme_minimal()

#### %N ----

PercentN_change %>%
  filter(SampleType == "SOIL", !is.na(ChangePercentN)) %>%
  mutate(Population = factor(Population, levels = c("FN", "YF", "SC", "HF", "DC", "SP", "MC", "UP"))) %>%
  ggplot(aes(x = Population, y = ChangePercentN, fill = TransplantTreatment)) +
  geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) + 
  labs(title = "Change in Soil %N Across Populations and Transplant Types",
       x = "Population",
       y = "%N") +
  theme_minimal()

### SORU ----
#### %C ----
PercentC_change %>%
  filter(SampleType == "SORU", !is.na(ChangePercentC)) %>%
  mutate(Population = factor(Population, levels = c("FN", "YF", "SC", "HF", "DC", "SP", "MC", "UP"))) %>%
  ggplot(aes(x = Population, y = ChangePercentC, fill = TransplantTreatment)) +
  geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) + 
  labs(title = "Change in SORU %C Across Populations and Transplant Types",
       x = "Population",
       y = "%C") +
  theme_minimal()

#### %N ----
PercentN_change %>%
  filter(SampleType == "SORU", !is.na(ChangePercentN)) %>%
  mutate(Population = factor(Population, levels = c("FN", "YF", "SC", "HF", "DC", "SP", "MC", "UP"))) %>%
  ggplot(aes(x = Population, y = ChangePercentN, fill = TransplantTreatment)) +
  geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) + 
  labs(title = "Change in SORU %N Across Populations and Transplant Types",
       x = "Population",
       y = "%N") +
  theme_minimal()

### POPRC ----
#### %C ----
PercentC_change %>%
  filter(SampleType == "POPRC", !is.na(ChangePercentC)) %>%
  mutate(Population = factor(Population, levels = c("FN", "YF", "SC", "HF", "DC", "SP", "MC", "UP"))) %>%
  ggplot(aes(x = Population, y = ChangePercentC, fill = TransplantTreatment)) +
  geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) + 
  labs(title = "Change in POPRC %C Across Populations and Transplant Types",
       x = "Population",
       y = "%C") +
  theme_minimal()

#### %N ----
PercentN_change %>%
  filter(SampleType == "POPRC", !is.na(ChangePercentN)) %>%
  mutate(Population = factor(Population, levels = c("FN", "YF", "SC", "HF", "DC", "SP", "MC", "UP"))) %>%
  ggplot(aes(x = Population, y = ChangePercentN, fill = TransplantTreatment)) +
  geom_boxplot(alpha = 0.7, position = position_dodge(width = 0.8)) + 
  labs(title = "Change in POPRC %N Across Populations and Transplant Types",
       x = "Population",
       y = "%N") +
  theme_minimal()

