# Ecosystem effects of predators are amplified across generations through prey behavioral plasticity
[Manuscript in review]

Authors: Nathalie R. Sommer, Annise M. Dobson, Matthew S. Baker, Geoffrey C. Trussell, Oswald J. Schmitz

Code by Nathalie R. Sommer



This repository contains all data, scripts, and documentation necessary to reproduce the analyses and figures presented in our manuscript. Upon acceptance, this README will be updated with the manuscript DOI.

## Table of Contents

- [Project Overview](#project-overview)
- [Data](#data)
- [Scripts](#scripts)
- [Figures](#figures)
- [Dependencies](#dependencies)
- [Reproducibility](#reproducibility)

## Project Overview

The research investigates the impact of predators (*P. mira*) on transgenerational differences in prey traits (*M. femurrubrum*) and various ecosystem variables. The data comes from a 3-year common garden experiment, using a BACI (before-after-control-impact) design in field mesocosms, paired with trait assessments. Further details on this experiment can be found in the manuscript.

## Data

All data is contained within the `Data` folder in their raw, original form. Descriptions of each file are contained below, with all data processing and analyses occurring in the [Scripts](#scripts).

### Ecosystems Data
- Data/CN_Tidy_predators.csv: Contains all data on C and N elemental content for plant tissue, litter, and soil
- Data/Vegetation.cover_predators.csv: Contains all data on plant percent cover within each mesocosm
- Data/Biomass-plots_predators.csv: Contains data on plots external to the mesocosms to build allometric equations
- Data/N-min.csv: Contains end-point (2023) data on nitrate and ammonium concentrations from samples on Day 0 and Day 30
- Data/BulkDensity.csv: Contains bulk density estimates from the common garden
- Data/2021_N-min.csv: Contains initial (2021) data on nitrate and ammonium concentrations
- Data/SIR: Folder containing all data pertaining to substrate-induced respiration to estimate microbial biomass
  - SIR Prep_2021.csv: Initial sample processing and calculations of target dry-mass weight for the initial (2021) data
  - SIR Prep_predators.csv: Initial sample processing and calculations of target dry-mass weight for the end-point (2023) data
  - SoilGWC_2021.csv: Gravimetric moisture content for initial (2021) data, by mesocosm ID
  - SoilGWC_predators.csv: Gravimetric moisture content for end-point (2023) data, by mesocosm ID
  - Direct outputs from the IRGA containing data on CO2 accumulation:
    - soil_sir_1.1.csv
    - soil_sir_11.08.csv
    - soil_sir_11.14.csv
    - soil_sir_11.17.csv
    - soil_sir_11.19.csv
    - soil_sir_11.21.csv
    - soil_sir_11.27.csv
    - soil_sir_11.28.csv
    - soil_sir_11.29.csv
    - soil_sir_12.27.csv
    - soil_sir_12.29.csv

### Trait Data
- Data/Behavior_Pred-G1G2_Herb-G2.csv: Contains behavioral data for predators and herbivores (G2 only) across generations
- Data/Behavior_Temps.csv: Contains temperature data associated with behavioral measurements
- Data/Respriation_Mass_2023_F1_F2.csv: Contains respiration data with mass measurements for F1 and F2 generations
- Data/Respriation_2023_F1_F2_Raw.csv: Contains raw respiration data for F1 and F2 generations

## Scripts

All scripts are contained within the `Scripts` folder. Use the code outline panel in RStudio for easiest navigation through this set of scripts.

### Structural Equation Model on Ecosystems Data

- **Data Cleaning**:
  - Scripts/N-min.R: Processes raw data to calculate nitrogen-mineralization rate from ammonium and nitrate data
  - Scripts/SIR.R: Processes raw CO2 values to calculate substrate-induced respiration rates
  - Scripts/Veg.R: Processes plant biomass and plots to calculate allometric relationships, calculates richness from percent cover in mesocosms
- **Analysis**: Scripts/MainAnalysis.R
  - Script for conducting and visualizing the SEM, contains the following sections:
    - Cleans and merges datasets for analysis
    - Builds and evaluates SEM model
    - Generates SEM path diagrams using DiagrammeR

### Trait Assessment
- **Scripts/Trait-Assessment.R**: Processes and analyses trait data, contains the following sections:
  - Behavior:
    - Data Processing
    - Analysis
    - Figures
  - Respiration:
    - Data Processing
    - Analysis
    - Figures


## Dependencies

The analyses were conducted in R (v.4.4.1) with the following packages:

- `dplyr`(v.1.1.4): Data cleaning
- `lme4` (v.1.1-35.5): Linear mixed-effects models
- `DHARma` (v.0.4.6): Evaluate LMM assumptions
- `car` (v.3.1-30): Evaluate LMM assumptions
- `boot` (v.1.3-30): Non-parametric bootstrapping to validate LMMs
- `ggplot2`(v.3.5.1): Data visualization
- `piecewiseSEM` (v.2.3.0.1): SEM analysis on ecosystems data
- `DiagrammeR` (v.1.0.11): Graph and network visualization
- `emmeans` (v.10.5): Post-hoc contrasts on trait data

## Reproducibility

To reproduce the analysis and figures:

1. Clone the repository from GitHub:
   ```bash
   git clone https://github.com/nathaliesommer/PredatorTransgenerationalPlasticity-Ecosystems.git
   ```
2. Set your working directory to the cloned repository.
3. Run the scripts in the order specified:
   - Start with `Scripts/N-min.R`, `Scripts/SIR.R`, and `Scripts/Veg.R` to prepare the data.
   - Proceed with `Scripts/MainAnalysis.R` for SEM analyses and visualization.
   - Use `Scripts/Trait-Assessment.R` for trait analyses and visualization.
4. The results of the SEM analysis will be rendered by the script, documented in:
- **SEM_Summary_Initial.txt**: Initial SEM model results
- **SEM-Reduced_Summary_Initial.txt**: SEM model results, after tests of directed separation
- **Indirect_paths.txt**: Calculated indirect path estimates for the G1 and G2 comparison
- **Trophic_Impact_Indirect_Paths.txt**: Calculated indirect path estimates for the predator and herbivore comparison
5. All other results will print directly in the console.
