# ____________________________________________________________________________________
# BIOINDICATOR–CHEMICAL METRIC + ABIOTIC FACTORS ANALYSIS
# Thesis Project: "Evaluating the Relationship Between Pesticide Pressure and 
# Macroinvertebrate Bioindicators in Swiss Streams"

# Author: Anthony Jeanpier Fow Esteves
# Date: 18.06.2024
# ____________________________________________________________________________________

# Description:
# This script performs the main statistical analyses for the thesis project, integrating synthetic biological,
# chemical, and abiotic datasets to assess relationships and potential drivers of ecological risk.
#
# Workflow summary:
#   1. Initial analyses explore chemical and biological metrics (correlations, patterns, PCA).
#   2. Subsequent steps incorporate env variables (here called abiotic factors) to refine models and 
#   evaluate the influence of bioindicator responses.
#   3. Regression models quantify associations between chemical metrics, other abiotic factors, and biological indicators,
#      with verification of statistical assumptions.
#   4. Cross-validation procedures assess model robustness and predictive power.
#
# IMPORTANT:
# - The datasets used here are entirely synthetic and randomly generated for reproducibility and code demonstration.
# - They closely replicate the structure, columns, and approximate number of observations and streams of the original data,
#   but do NOT contain real observations. Results from this script cannot be interpreted as the actual thesis findings.
# - All data transformations (e.g., log1p, scaling) must be performed as specified in the code.
# - For full variable documentation, see the `/docs` folder.

# Packages ----
Sys.setenv(LANG="EN")
if (!require("forcats")) { install.packages("forcats"); library(forcats) }
if (!require("tidyverse")) { install.packages("tidyverse"); library(tidyverse) }
if (!require("broom")) { install.packages("broom"); library(broom) }
if (!require("car")) { install.packages("car"); library(car) }
if (!require("patchwork")) { install.packages("patchwork"); library(patchwork) }
if (!require("corrplot")) { install.packages("corrplot"); library(corrplot) }
if (!require("lme4")) { install.packages("lme4"); library(lme4) }
if (!require("lmerTest")) { install.packages("lmerTest"); library(lmerTest) }
if (!require("broom.mixed")) { install.packages("broom.mixed"); library(broom.mixed) }
if (!require("performance")) { install.packages("performance"); library(performance) }
if (!require("ggrepel")) { install.packages("ggrepel"); library(ggrepel) }
if (!require("mgcv")) { install.packages("mgcv"); library(mgcv) }
if (!require("ggpubr")) { install.packages("ggpubr"); library(ggpubr) }
if (!require("gridExtra")) { install.packages("gridExtra"); library(gridExtra) }
if (!require("effects")) { install.packages("effects"); library(effects) }
if (!require("viridis")) { install.packages("viridis"); library(viridis) }
if (!require("caret")) { install.packages("caret"); library(caret) }
if (!require("randomForest")) { install.packages("randomForest"); library(randomForest) }
if (!require("lindia")) { install.packages("lindia"); library(lindia) }
if (!require("MuMIn")) { install.packages("MuMIn"); library(MuMIn) }
if (!require("dplyr")) { install.packages("dplyr"); library(dplyr) }
if (!require("tibble")) { install.packages("tibble"); library(tibble) }
if (!require("purrr")) { install.packages("purrr"); library(purrr) }
if (!require("ggeffects")) { install.packages("ggeffects"); library(ggeffects) }
if (!require("stringi")) { install.packages("stringi"); library(stringi) }
if (!require("GGally")) { install.packages("GGally"); library(GGally) }

detach("package:raster", unload = TRUE) # it has conflict with select.

# Activate to check in case any conflict happen, so you can explore it
# detach("package:raster", unload = TRUE) #raster pckg generates problem with "select", be careful
# (.packages())


# Paths to call ----
## user <- "C:\\your user...\\"
user <- "C:\\Users\\ajfe0\\Desktop\\SPEAR_clean25"
file_inputs <- "\\intpus\\" # Here you already have the three data sets

# Paths to save ----
output_graphics_file <- "\\output\\graphics\\" # To save plots

# Call of the necessary data sets from inputs ----
Abiotic_factors_41_streams_ds_prepared <- readRDS(paste0(user, file_inputs, "Abiotic_factors_41_streams_ds_prepared.rds"))
metrics_3.5d_0_R12_pp_ds <- readRDS(paste0(user, file_inputs, "metrics_3.5d_0_R12_pp_ds.rds"))
metrics_14d_0_R12_pp_ds <- readRDS(paste0(user, file_inputs, "metrics_14d_0_R12_pp_ds.rds"))

# First recognize the structure of each data set
str(Abiotic_factors_41_streams_ds_prepared)
str(metrics_3.5d_0_R12_pp_ds)
str(metrics_14d_0_R12_pp_ds)

# 1. Descriptive Analysis ----

## Organization of variables ----
# Metrics:
## Bioindicators used for the analysis
bioindicators <- c("SPEAR", "IBCH_2019", "EPT", "VT_VALUE", "GI_VALUE")

## Chemical metrics (RQs and TUs) under "mixture" and "maximum" approaches and their
## statistical representation (mean, maximum and median)
chemical_metrics_3.5d <- c("ARQmix_mean", "ARQmix_median", "ARQmix_max",
                           "TU_ECmix_mean", "TU_ECmix_median", "TU_ECmix_max",
                           "ARQmax_mean", "ARQmax_median", "ARQmax_max",
                           "TU_ECmax_mean", "TU_ECmax_median", "TU_ECmax_max")
chemical_metrics_14d <-  c("CRQmix_mean", "CRQmix_median", "CRQmix_max",
                           "TU_NOECmix_mean", "TU_NOECmix_median", "TU_NOECmix_max",
                           "CRQmax_mean", "CRQmax_median", "CRQmax_max",
                           "TU_NOECmax_mean", "TU_NOECmax_median", "TU_NOECmax_max")

## Scenarios (time windows in the report); the scenarios set the period or range of time in which 
## a selected group of chemical samples were taken to calculate the chemical metrics before each
## macroinvertebrate monitoring. Each scenario was pre-defined in "data_preparation.R" script.
scenarios_3.5d <- c("3.5_days_1_week", "3.5_days_2_weeks", "3.5_days_1_month", 
                    "3.5_days_2_months", "3.5_days_3_months", "3.5_days_6_months", 
                    "3.5_days_1_year")
scenarios_14d <- c("14_days_1_week", "14_days_2_weeks", "14_days_1_month", 
                    "14_days_2_months", "14_days_3_months", "14_days_6_months", 
                    "14_days_1_year")

# Arrangement of the datasets by scenario (time windows in the report)
metrics_3.5d_0_R12_pp_ds %>% mutate(scenario = fct_relevel(scenario, scenarios_3.5d)) -> metrics_3.5d_0_R12_pp_ds
metrics_14d_0_R12_pp_ds %>% mutate(scenario = fct_relevel(scenario, scenarios_14d)) -> metrics_14d_0_R12_pp_ds

## Frequencies per sub-set of data ----

### Number of observations per year and scenario
metrics_3.5d_0_R12_pp_ds %>% # Here you can replace it by the metrics ds you want to use
  group_by(YEAR,scenario) %>% 
  summarise(n = n()) %>% 
  tidyr::pivot_wider(names_from = YEAR, values_from = n)

### Range of years and number of observations per streams and scenarios (time windows in the report)
metrics_14d_0_R12_pp_ds %>% group_by(Gewasser) %>%
  group_by(Gewasser, scenario) %>%
  summarise(
    year_min = min(YEAR, na.rm = TRUE),
    year_max = max(YEAR, na.rm = TRUE),
    year_range = paste(min(YEAR, na.rm = TRUE), max(YEAR, na.rm = TRUE), sep = "-"),
    n = n(),
    .groups = "drop") %>%
  select(Gewasser, scenario, year_range, n)


## 1.1. Correlation between bioindicators and chemical metrics ----

### Function to run correlation ----

compute_scenario_correlations <- function(
    data, metrics_1, metrics_2,
    transform_1 = "none",  
    transform_2 = "none",  
    cor_method = "spearman" 
) {
  
  # ======================================================================
  # Purpose: Compute correlations between sets of metrics for each scenario,
  #          with optional transformations and tests for normality and 
  #          homoscedasticity.
  #
  # Inputs:
  #   data: Data frame containing the data, with a 'scenario' column
  #   metrics_1: Character vector, names of first set of metrics
  #   metrics_2: Character vector, names of second set of metrics
  #   transform_1: Transformation for metrics_1 ("none", "log", "sqrt")
  #   transform_2: Transformation for metrics_2 ("none", "log", "sqrt")
  #   cor_method: Correlation method ("spearman", "pearson", "kendall")
  #
  # Output: Data frame with correlation results for each scenario and metric pair,
  #         including correlation coefficients, p-values, method, and test results.
  # ======================================================================
  
  # Helper function to apply transformation to a vector
  transform_vec <- function(x, method) {
    if (method == "log") {
      return(log(x + 1))       # Log transform (add 1 to avoid log(0))
    } else if (method == "sqrt") {
      return(sqrt(x))          # Square root transform
    } else {
      return(x)                # No transformation
    }
  }
  
  # Split data by scenario into list of data frames
  data_by_scenario <- data %>%
    group_by(scenario) %>% 
    group_split()
  
  # Function to compute correlations for a single scenario
  compute_corr_scenario <- function(scenario_data, metrics_1, metrics_2) {
    # Generate all combinations of metrics_1 and metrics_2
    combinations <- expand.grid(metrics_1 = metrics_1, metrics_2 = metrics_2, 
                                stringsAsFactors = FALSE)
    
    # Function to compute correlation and test diagnostics for a pair of metrics
    compute_corr <- function(metrics_1, metrics_2) {
      # Apply transformations
      x <- transform_vec(scenario_data[[metrics_1]], transform_1)
      y <- transform_vec(scenario_data[[metrics_2]], transform_2)
      
      # Compute correlation
      cor_test <- cor.test(x, y, method = cor_method)
      cor_res <- broom::tidy(cor_test)
      
      # Test for normality (Shapiro-Wilk)
      norm_x <- tryCatch(shapiro.test(x)$p.value, error = function(e) NA_real_)
      norm_y <- tryCatch(shapiro.test(y)$p.value, error = function(e) NA_real_)
      
      # Test for homoscedasticity (Levene's test, using median split of x)
      group <- as.factor(ifelse(x > median(x, na.rm = TRUE), "high", "low"))
      lev_p <- tryCatch(car::leveneTest(y ~ group)$`Pr(>F)`[1], error = function(e) NA_real_)
      
      # Return results as a tibble
      tibble::tibble(
        estimate = cor_res$estimate,
        p.value = cor_res$p.value,
        method = cor_method,
        norm_x_p = norm_x,
        norm_y_p = norm_y,
        levene_p = lev_p
      )
    }
    
    # Apply compute_corr to all metric combinations and add scenario info
    results <- combinations %>%
      purrr::pmap_dfr(compute_corr) %>%
      mutate(
        metrics_1 = combinations$metrics_1,
        metrics_2 = combinations$metrics_2,
        scenario = unique(scenario_data$scenario),
        transform_1 = transform_1,
        transform_2 = transform_2
      ) %>%
      select(scenario, metrics_1, metrics_2, transform_1, transform_2,
             estimate, p.value, method, norm_x_p, norm_y_p, levene_p)
    
    return(results)
  }
  
  # Process all scenarios and combine results into single dataframe
  all_results <- purrr::map_dfr(
    data_by_scenario, 
    compute_corr_scenario, 
    metrics_1 = metrics_1, 
    metrics_2 = metrics_2
  )
  
  return(all_results)  # Return final results dataframe
}

### Generation of the datasets with the correlation values using "compute_scenario_correlations"

### 3.5-days sampling period
compute_scenario_correlations(
  data = metrics_3.5d_0_R12_pp_ds,
  metrics_1 = chemical_metrics_3.5d,
  metrics_2 = bioindicators,
  transform_1 = "log", 
  transform_2 = "none",
  cor_method = "spearman") -> cor_3.5_bio_chem_ds
### 14-days sampling period
compute_scenario_correlations(
  data = metrics_14d_0_R12_pp_ds,
  metrics_1 = chemical_metrics_14d,
  metrics_2 = bioindicators,
  transform_1 = "log",
  transform_2 = "none",
  cor_method = "spearman") -> cor_14_bio_chem_ds

### Correlograms ----

### 3.5d.plot
ggplot(
  cor_3.5_bio_chem_ds %>% 
    mutate(scenario = fct_relevel(scenario, scenarios_3.5d),
           label = ifelse(p.value > 0.05, "ns", round(estimate, 2)),
           fill_signif = ifelse(p.value > 0.05, NA, estimate)),
  aes(x = metrics_1, y = metrics_2, fill = fill_signif)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), color = "black", size = 2, angle = 45) +
  scale_fill_gradient2(
    low = "#B2182B", mid = "white", high = "#2166AC", # for specific colors "google picker" was used
    midpoint = 0, limit = c(-1, 1),
    na.value = "gray80",
    name = "Correlation") +
  facet_grid(. ~ scenario, scales = "free", space = "free") +
  labs(x="Chemical metrics (log(x+1))", 
       y="Bioindicators")+
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 9)) -> a # Here be careful, I used words like "a", "b" or "p" to momentaneouly save the plots to be exported
### 14d
ggplot(
  cor_14_bio_chem_ds %>% 
    mutate(scenario = fct_relevel(scenario, scenarios_14d),
           label = ifelse(p.value > 0.05, "ns", round(estimate, 2)),
           fill_signif = ifelse(p.value > 0.05, NA, estimate)),
  aes(x = metrics_1, y = metrics_2, fill = fill_signif)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), color = "black", size = 2, angle = 45) +
  scale_fill_gradient2(
    low = "#B2182B", mid = "white", high = "#2166AC",
    midpoint = 0, limit = c(-1, 1),
    na.value = "gray80",
    name = "Correlation") +
  facet_grid(. ~ scenario, scales = "free", space = "free") +
  labs(x="Chemical metrics (log(x+1))", 
       y="Bioindicators")+
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 9)) -> b # Here be careful, I used words like "a", "b" or "p" to momentaneouly save the plots to be exported

### Export the plot:
#png(paste(user, output_graphics_file ,"Bioind_chem_metrics_corr_CORRECTED.png", sep = ""), 
#    width = 40, height = 30, units = "cm", res = 350)
((a / b) +
    plot_annotation(title = "Correlation Heatmap across time window on each sampling period",
                    caption = "Values within each cell represent Spearman’s ρ\n ns: No significant correlation (p-value < 0.05)") +
    plot_layout(guides = "collect")) & 
  theme(
    plot.title        = element_text(face = "bold", size = 14),
    legend.position  = "right",
    legend.title     = element_text(face = "bold", size = 11),
    legend.text      = element_text(size = 10),
    legend.key       = element_rect(fill = "white", colour = "black", linetype = "solid"),
    legend.background = element_rect(fill = NA, colour = "black", linetype = "solid"),
    strip.background = element_rect(fill = "white", colour = "gray12"),
    strip.text       = element_text(face = "bold"),
    axis.text.x      = element_text(angle = 45, 
                                    #vjust = 0.5, 
                                    hjust = 1))

#dev.off()
rm(a, b) #to avoid keep them in the work space

## Pre 1.2. Data-sets preparation ----
### Data set: Metrics 3.5d + static abiotic factors + transformation (e.g. scaling and log1p) -----
left_join(metrics_3.5d_0_R12_pp_ds,
          # Let's use the abiotic ds prepared
          Abiotic_factors_41_streams_ds_prepared %>% 
            # Code to clean the stream (Gewasser name) so that the 
            # the abiotic factor values match with the stream's names
            # from metrics_ds
            mutate(
              Gewasser_clean = Gewasser %>% 
                str_to_lower() %>%
                str_trim() %>%
                stri_trans_general("Latin-ASCII") %>% 
                str_replace_all(" ", "")), 
          by = "Gewasser_clean") %>%
  # Re-name and removal of the Gewasser columns
  mutate(Gewasser = Gewasser.x) %>% 
  select(-Gewasser.y, -Gewasser.x) %>% 
  # Selection of the relevant abiotic factors
  select(
    # Location
    Kanton, Gewasser, scenario, 
    #Categories from the initial ds
    YEAR, Z, Season_MZB, Kategorie_FG_Grosse, 
    #Hydrological variables
    Catch_area_m2, Avg_mod_discharge_m3_s, Proportion_of_wastewater_l_s, FLOZ, Flow_velocity,
    #Morphological variables (the last ecomorphology value used was "Ecomorphology_0_12")
    Ecomorphology_class, Ecomorphology_categories23, Ecomorphology_values23, 
    Stream_bed_construction, Ecomorphology_cont23_Langhans, Ecomorphology_0_12,
    #Land use variables
    Urban_area_frac, Forest_area_frac, Agricultural_area_frac,
    # Modelled water temperature
    Mod_max_temp_summer,
    #Chemical metrics (dynamic)
    all_of(chemical_metrics_3.5d),
    #Precipitation values (dynamic)
    contains("PP"),
    # Bioindicators
    SPEAR, IBCH_2019, EPT, GI_VALUE, VT_VALUE) %>% # str() %>% 
  mutate(
    Gewasser = as.factor(Gewasser),
    YEAR = as.factor(YEAR)) %>% 
  # HERE IS WERE I APPLIED THE SCALE AND LOGP1 TRANSFORMATION
  ## Groupy by scenarios bcs I treat each sub-set of data differently
  group_by(scenario) %>% 
  ## To check the number of point with catchment above 30km2 (so far 6 out of 22)
  ## "30 km2" is a limitation of the SPEAR.
  #dplyr::filter(scenario == "3.5_days_1_year") %>% 
  #mutate(Catch_area_m2 = (Catch_area_m2/1000000)) %>% 
  #dplyr::select(Gewasser, Catch_area_m2) %>% 
  #distinct() %>% summarise(min = min(Catch_area_m2),
  #                         max = max(Catch_area_m2))
  # Scaling of the abiotic factors
  mutate(
    across(
      c(Z, 
        Catch_area_m2, Avg_mod_discharge_m3_s, Proportion_of_wastewater_l_s, Flow_velocity, 
        Urban_area_frac, Forest_area_frac, Agricultural_area_frac,
        Mod_max_temp_summer, Ecomorphology_cont23_Langhans, Ecomorphology_0_12,
        contains("PP_")), scale),
      # Log1p and scaling of the chemical metrics
      across(all_of(chemical_metrics_3.5d), 
           ~ scale(log1p(.)))) %>% 
  # To ensure the loss of unused levels.
  ungroup() %>% droplevels() %>% 
  # Clasification of the various ecomorphology categories tried
  ## Re-group of stream bed type because it has very low observations, and #30-60% is not present
  mutate(
    Stream_bed_construction = case_when(
      Stream_bed_construction %in% "0%, natural" ~ "1: 0%, Completely natural",
      Stream_bed_construction %in% "<10%" ~ "2: <10%, Very low construction",
      Stream_bed_construction %in% c(">60%", "100%, non-natural") ~ "3: >60%, High to non-natural\nconstruction",
      TRUE ~ NA_character_)) %>% 
  mutate(Stream_bed_construction = 
           factor(Stream_bed_construction,
                  levels = c(
                    "1: 0%, Completely natural",
                    "2: <10%, Very low construction",
                    "3: >60%, High to non-natural\nconstruction"))) %>% 
  ## Clasical ecomorphology from 0 - 4
  mutate(
    Ecomorphology_categories23 = case_when(
      Ecomorphology_categories23 %in% "Natürlich, naturnah" ~ "Natural",
      Ecomorphology_categories23 %in% "Wenig beeinträchtigt" ~ "Slightly impaired",
      Ecomorphology_categories23 %in% "Stark beeinträchtigt" ~ "Severely impaired",
      Ecomorphology_categories23 %in% "Naturfremd, künstlich" ~ "Unnatural, artificial",
      Ecomorphology_categories23 %in% "Eingedolt" ~ "Covered",
      Ecomorphology_categories23 %in% "Nicht bestimmt" ~ "Not classified",
      TRUE ~ NA_character_)) %>% 
  ## Clasical ecomorphology from 0 - 4
  mutate(Ecomorphology_categories23 = 
           factor(Ecomorphology_categories23,
                  levels = c(
                    "Natural",
                    "Slightly impaired",
                    "Severely impaired",
                    "Unnatural, artificial",
                    "Covered"))) -> metrics_3.5d_0_R12_abiotic_ds

# Check of missing data and structure
metrics_3.5d_0_R12_abiotic_ds %>% str()
colMeans(is.na(metrics_3.5d_0_R12_abiotic_ds))*100

### Data set: Mean and standard deviation of the raw data (metrics_3.5d_0_R12_abiotic_ds) -----
#____________________________________________________________________________
# "SD_DATASET" IS TO BE USED AS THE REFERENCE
# TO GET THE CHANGE MEANT BY "1-SD" OF INCREASE, AS WELL AS THE
# BASE LINE USED BY EACH REGRESSION MODEL TO GET THE INTERCEPT. HOWEVER, NO
# DIRECT CHANGES IN UNITS MUST BE DONE FOR PLOT BECAUSE THE VALUES CAN BE MISINTERPRETATED
# THIS IS ONLY FOR EXPLORATION
#____________________________________________________________________________

### This is useful to evaluate the baseline values for the interpretation of 
### the regression models
left_join(metrics_3.5d_0_R12_pp_ds,
          # Let's use the abiotic ds prepared
          Abiotic_factors_41_streams_ds_prepared %>% 
            # Code to clean the stream (Gewasser name) so that the 
            # the abiotic factor values match with the stream's names
            # from metrics_ds
            mutate(
              Gewasser_clean = Gewasser %>% 
                str_to_lower() %>%
                str_trim() %>%
                stri_trans_general("Latin-ASCII") %>% 
                str_replace_all(" ", "")), 
          by = "Gewasser_clean") %>%
  mutate(Gewasser = Gewasser.x) %>% 
  select(-Gewasser.y, -Gewasser.x) %>% 
  select(
    scenario, 
    #Hydrological variables
    Catch_area_m2, Avg_mod_discharge_m3_s, Proportion_of_wastewater_l_s, FLOZ, Flow_velocity,
    #Morphological variables
    Ecomorphology_cont23_Langhans, Ecomorphology_0_12,
    #Land use variables
    Urban_area_frac, Forest_area_frac, Agricultural_area_frac,
    # Modeled water temperature
    Mod_max_temp_summer,
    #Chemical metrics (dynamic)
    all_of(chemical_metrics_3.5d),
    #Precipitation values (dynamic)
    contains("PP")) %>% 
  drop_na() %>% 
  # First, let's convert the chemical metrics into log1p
  mutate(
    across(all_of(chemical_metrics_3.5d), ~ log1p(.))) %>%
  # Second, let's wrap the factors from which I need the mean and sd 
  gather(key = term, 
         value = value, Catch_area_m2:ncol(.)) %>% 
  # Finally, calculation of the sd and mean based on scenario (time window) and
  # term. Since, each scenario represent the sub-set of data from which each model
  # is performed. Be careful that the SD and mean is based on log1p
  group_by(scenario, term) %>% 
  summarise(SD = sd(value), mean = mean(value)) -> sd_term_3.5d

### Data set: Metrics 14d + static abiotic factors + transformation (e.g. scaling and log1p) -----
left_join(metrics_14d_0_R12_pp_ds,
          # Let's use the abiotic ds prepared
          Abiotic_factors_41_streams_ds_prepared %>% 
            # Code to clean the stream (Gewasser name) so that the 
            # the abiotic factor values match with the stream's names
            # from metrics_ds
            mutate(
              Gewasser_clean = Gewasser %>% 
                str_to_lower() %>%
                str_trim() %>%
                stri_trans_general("Latin-ASCII") %>% 
                str_replace_all(" ", "")), 
          by = "Gewasser_clean") %>%
  # Re-name and removal of the Gewasser columns
  mutate(Gewasser = Gewasser.x) %>% 
  select(-Gewasser.y, -Gewasser.x) %>% 
  # Selection of the relevant abiotic factors
  select(
    # Location
    Kanton, Gewasser, scenario, 
    #Categories from the initial ds
    YEAR, Z, Season_MZB, Kategorie_FG_Grosse, 
    #Hydrological variables
    Catch_area_m2, Avg_mod_discharge_m3_s, Proportion_of_wastewater_l_s, FLOZ, Flow_velocity,
    #Morphological variables
    Ecomorphology_class, Ecomorphology_categories23, Ecomorphology_values23, 
    Stream_bed_construction, Ecomorphology_cont23_Langhans, Ecomorphology_0_12,
    #Land use variables
    Urban_area_frac, Forest_area_frac, Agricultural_area_frac,
    # Modelled water temperature
    Mod_max_temp_summer,
    #Chemical metrics (dynamic)
    all_of(chemical_metrics_14d),
    #Precipitation values (dynamic)
    contains("PP"),
    # Bioindicators
    SPEAR, IBCH_2019, EPT, GI_VALUE, VT_VALUE) %>% # str() %>% 
  mutate(
    Gewasser = as.factor(Gewasser),
    YEAR = as.factor(YEAR)) %>% 
  # HERE IS WERE I APPLIED THE SCALE AND LOGP1 TRANSF
  ## First remove the "not categorised ecomorphologies
  group_by(scenario) %>% 
  ## To check the number of point with catchment above 30km2 (so far 19 out of 41)
  #dplyr::filter(scenario == "14_days_1_year") %>% 
  #mutate(Catch_area_m2 = (Catch_area_m2/1000000)) %>% 
  #dplyr::select(Gewasser, Catch_area_m2) %>% 
  #distinct() %>% View() summarise(min = min(Catch_area_m2),
  #                         max = max(Catch_area_m2))
  # Scaling of the abiotic factors
  mutate(
    across(
      c(Z, 
        Catch_area_m2, Avg_mod_discharge_m3_s, Proportion_of_wastewater_l_s, Flow_velocity, 
        Urban_area_frac, Forest_area_frac, Agricultural_area_frac,
        Mod_max_temp_summer, Ecomorphology_cont23_Langhans, Ecomorphology_0_12,
        contains("PP_")), scale),
    # Log1p and scaling of the chemical metrics
    across(all_of(chemical_metrics_14d), 
           ~ scale(log1p(.)))) %>% 
  # To ensure the loss of unused levels.
  ungroup() %>% droplevels() %>% 
  # Re-group of stream bed type because it has very low observations, and #30-60% is not present
  mutate(
    Stream_bed_construction = case_when(
      Stream_bed_construction %in% "0%, natural" ~ "0%, Completely natural",
      Stream_bed_construction %in% c("<10%", "10-30%") ~ "<10-30%, Very low to low and\nlow construction",
      Stream_bed_construction %in% c(">60%", "100%, non-natural") ~ ">60%, High to non-natural\nconstruction",
      TRUE ~ NA_character_)) %>% 
  mutate(Stream_bed_construction = 
           factor(Stream_bed_construction,
                  levels = c(
                    "0%, Completely natural",
                    "<10-30%, Very low to low and\nlow construction",
                    ">60%, High to non-natural\nconstruction"))) %>% 
  ## Clasical ecomorphology from 0 - 4
  mutate(
    Ecomorphology_categories23 = case_when(
      Ecomorphology_categories23 %in% "Natürlich, naturnah" ~ "Natural",
      Ecomorphology_categories23 %in% "Wenig beeinträchtigt" ~ "Slightly impaired",
      Ecomorphology_categories23 %in% "Stark beeinträchtigt" ~ "Severely impaired",
      Ecomorphology_categories23 %in% "Naturfremd, künstlich" ~ "Unnatural, artificial",
      Ecomorphology_categories23 %in% "Eingedolt" ~ "Covered",
      Ecomorphology_categories23 %in% "Nicht bestimmt" ~ "Not classified",
      TRUE ~ NA_character_)) %>% 
  ## Clasical ecomorphology from 0 - 4
  mutate(Ecomorphology_categories23 = 
           factor(Ecomorphology_categories23,
                  levels = c(
                    "Natural",
                    "Slightly impaired",
                    "Severely impaired",
                    "Unnatural, artificial",
                    "Covered",
                    "Not classified")))-> metrics_14d_0_R12_abiotic_ds

# Check of missing data and structure
metrics_14d_0_R12_abiotic_ds %>% str()
colMeans(is.na(metrics_14d_0_R12_abiotic_ds))*100


### Data set: Mean and standard deviation of the raw data (metrics_14d_0_R12_abiotic_ds) -----
#____________________________________________________________________________
# "SD_DATASET" IS TO BE USED AS THE REFERENCE
# TO GET THE CHANGE MEANT BY "1-SD" OF INCREASE, AS WELL AS THE
# BASE LINE USED BY EACH REGRESSION MODEL TO GET THE INTERCEPT. HOWEVER, NO
# DIRECT CHANGES IN UNITS MUST BE DONE FOR PLOT BECAUSE THE VALUES CAN BE MISINTERPRETATED
# THIS IS ONLY FOR EXPLORATION
#____________________________________________________________________________
### This is useful to evaluate the baseline values for the interpretation of 
### the regression models
left_join(metrics_14d_0_R12_pp_ds,
          Abiotic_factors_41_streams_ds_prepared %>% 
            mutate(
              Gewasser_clean = Gewasser %>% 
                str_to_lower() %>%
                str_trim() %>%
                stri_trans_general("Latin-ASCII") %>% 
                str_replace_all(" ", "")), 
          by = "Gewasser_clean") %>%
  mutate(Gewasser = Gewasser.x) %>% 
  select(-Gewasser.y, -Gewasser.x) %>% 
  select(
    scenario,
    #Hydrological variables
    Catch_area_m2, Avg_mod_discharge_m3_s, Proportion_of_wastewater_l_s, FLOZ, Flow_velocity,
    #Morphological variables 
    Ecomorphology_cont23_Langhans, Ecomorphology_0_12,
    #Land use variables
    Urban_area_frac, Forest_area_frac, Agricultural_area_frac,
    # Modelled water temperature
    Mod_max_temp_summer,
    all_of(chemical_metrics_14d),
    #Precipitation values (dynamic)
    contains("PP")) %>% 
  drop_na() %>% 
  # First, let's convert the chemical metrics into log1p
  mutate(
    across(all_of(chemical_metrics_14d), ~ log1p(.))) %>% 
  # Second, let's wrap the factors from which I need the mean and sd 
  gather(key = term, 
         value = value, Catch_area_m2:ncol(.)) %>%
  # Finally, calculation of the sd and mean based on scenario (time window) and
  # term. Since, each scenario represent the sub-set of data from which each model
  # is performed. Be careful that the SD and mean is based on log1p
  group_by(scenario, term) %>% 
  summarise(SD = sd(value),
            mean = mean(value)) -> sd_term_14d

## 1.2. Correlation between bioindicators and selected chemical metrics + abiotic factors ----
## Organization of variables ----
desired_order_abiotic <- c(
  # Spatial variables 
  "Z", 
  # Hydrological variables
  "Catch_area_m2", "Avg_mod_discharge_m3_s", "Proportion_of_wastewater_l_s", "Flow_velocity",
  # Land-use variables
  "Urban_area_frac", "Forest_area_frac", "Agricultural_area_frac",
  # Temporal
  "Mod_max_temp_summer",
  "PP_mean_mean",
  # Morphological
  "Ecomorphology_0_12")
desired_order_abiotic_3.5d <- c(
  # Chemical metrics: Some representative (based on 1.1), so that the plot isn't overwhelmed
  "ARQmix_mean", "TU_ECmix_mean", 
  "ARQmix_median", "TU_ECmix_median", 
  "ARQmax_max", "TU_ECmax_max",
  # Other abiotic variables before listed 
  desired_order_abiotic)
desired_order_abiotic_14d <- c(
  # Chemical metrics: Some representative (based on 1.1), so that the plot isn't overwhelmed
  "CRQmix_mean", "TU_NOECmix_mean", 
  "CRQmix_median", "TU_NOECmix_median", 
  "CRQmax_max", "TU_NOECmax_max",
  # Other abiotic variables before listed
  desired_order_abiotic)

# As in the previous correlograms 
bioindicators <- c("SPEAR", "IBCH_2019", "EPT", "VT_VALUE", "GI_VALUE")
  
### Correlograms ----
compute_scenario_correlations(
  data = metrics_3.5d_0_R12_abiotic_ds,
  metrics_1 = bioindicators,
  metrics_2 = desired_order_abiotic_3.5d,
  transform_1 = "none",
  transform_2 = "none",
  cor_method = "spearman") %>% 
  mutate(scenario = fct_relevel(scenario, scenarios_3.5d),
           label = ifelse(p.value > 0.05, "ns", round(estimate, 2)),
           fill_signif = ifelse(p.value > 0.05, NA, estimate),
           metrics_1 = factor(metrics_1, 
                              level = bioindicators), 
           metrics_2 = factor(metrics_2, 
                              level = desired_order_abiotic_3.5d)) %>%
    #filter(abs(estimate) >= 0.3) %>%
    ggplot(aes(x = metrics_1, y = metrics_2, fill = fill_signif)) +
    geom_tile(color = "white") +
    geom_text(aes(label = label), color = "black", size = 3, angle = 45) +
    scale_fill_gradient2(
      low = "#B2182B", mid = "white", high = "#2166AC",
      midpoint = 0, limit = c(-1, 1),
      na.value = "gray80",
      name = "Correlation") +
    facet_grid(. ~ scenario, scales = "free", space = "free") +
    labs(tag = "A",
         subtitle = "3.5-day sampling period",
         x = "Abiotic factors", y = "Abiotic factors")+
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = 9)) -> a
  
compute_scenario_correlations(
  data = metrics_14d_0_R12_abiotic_ds,
  metrics_1 = bioindicators,
  metrics_2 = desired_order_abiotic_14d,
  transform_1 = "none",
  transform_2 = "none",
  cor_method = "spearman") %>% 
  mutate(scenario = fct_relevel(scenario, scenarios_14d),
         label = ifelse(p.value > 0.05, "ns", round(estimate, 2)),
         fill_signif = ifelse(p.value > 0.05, NA, estimate),
         metrics_1 = factor(metrics_1, level = bioindicators),
         metrics_2 = factor(metrics_2, level = desired_order_abiotic_14d)) %>% 
  #filter(abs(estimate) >= 0.3) %>%
  ggplot(aes(x = metrics_1, y = metrics_2, fill = fill_signif)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), color = "black", size = 3, angle = 45) +
  scale_fill_gradient2(
    low = "#B2182B", mid = "white", high = "#2166AC",
    midpoint = 0, limit = c(-1, 1),
    na.value = "gray80",
    name = "Correlation") +
  facet_grid(. ~ scenario, scales = "free", space = "free") +
  labs(tag = "B",
       subtitle = "14-day sampling period",
       x = "Bioindicators", y = "Abiotic factors")+
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 9)) -> b
  
#png(paste(user, output_graphics_file ,"Bioind_abiotic_corr_CORRECTED.png", sep = ""), 
#    width = 40, height = 30, units = "cm", res = 350)
((a / b) +
    plot_annotation(title = "Correlation Heatmap across time window on each sampling period",
                    caption = "Values within each cell represent Spearman’s ρ\n ns: No significant correlation (p-value < 0.05)") +
    plot_layout(guides = "collect")) & 
  theme(
    plot.title        = element_text(face = "bold", size = 14),
    legend.position  = "right",
    legend.title     = element_text(face = "bold", size = 11),
    legend.text      = element_text(size = 10),
    legend.key       = element_rect(fill = "white", colour = "black", linetype = "solid"),
    legend.background = element_rect(fill = NA, colour = "black", linetype = "solid"),
    strip.background = element_rect(fill = "white", colour = "gray12"),
    strip.text       = element_text(face = "bold"),
    axis.text.x      = element_text(angle = 45, 
                                    #vjust = 0.5, 
                                    hjust = 1))
#dev.off()
rm(a,b)
  
## 1.3. Correlation between chemical metrics + abiotic factors ----
### Here we performed the correlation between chemical metrics and abiotic factors in order to 
### assess the degree of relation between factors of the same dimensions or theoretically 
### associated. For this case we used only the sub-set of data "1-year" before the macroinvertebrate
### monitoring. Since the behavior of the correlations are quite similar across time windows.

#png(paste(user, output_graphics_file ,"Abiotc_abiotic_corr_CORRECTED.png", sep = ""), 
#    width = 50, height = 30, units = "cm", res = 350)

# Set the panel  
par(mfrow = c(1, 2))

# 3.5d sampling period
# Classic corrplot
corrplot(
  cor(
    # Only using 1-year ds
    metrics_3.5d_0_R12_abiotic_ds[metrics_3.5d_0_R12_abiotic_ds$scenario == "3.5_days_1_year", desired_order_abiotic_3.5d],
    method = "spearman",
    use = "pairwise.complete.obs"), 
  method = "number", 
  type = "upper",
  sig.level = 0.05,
  insig = "blank",
  addCoef.col = "black",
  col = colorRampPalette(c("#B2182B", "gray", "#2166AC"))(200),
  mar = c(0,0,1,0)) 
title("3.5 days sampling period", line = 1)

# 14d sampling period
# Classic corrplot
corrplot(
  cor(
    # Only using 1-year ds
    metrics_14d_0_R12_abiotic_ds[metrics_14d_0_R12_abiotic_ds$scenario == "14_days_1_year", desired_order_abiotic_14d],
    method = "spearman",
    use = "pairwise.complete.obs"), 
  method = "number", 
  type = "upper",
  sig.level = 0.05,
  insig = "blank",
  addCoef.col = "black",
  col = colorRampPalette(c("#B2182B", "gray", "#2166AC"))(200),
  mar = c(0,0,1,0)) 
title("14 days sampling period", line = 1)

par(mfrow = c(1, 1)) 
#dev.off()

# 2. Environmental Gradient Exploration: PCA Analysis ----

## Definition of the main dimensions for the variables used in the PCA
var_groups_3.5d <- list(
  Chemical_metrics = c("TU_ECmix_median", "ARQmix_median"),
  Land_Use = c("Agricultural_area_frac", "Urban_area_frac", "Forest_area_frac"),
  Spatial = c("Mod_max_temp_summer", "Z"),
  Morphology = "Ecomorphology_0_12",
  Hydrological = c("Flow_velocity", "Avg_mod_discharge_m3_s"),
  Temporal = c(#"PP_mean_mean"
    ))

var_groups_14d <- list(
  Chemical_metrics = c("CRQmix_max", "TU_NOECmix_max"),
  Land_Use = c("Agricultural_area_frac", "Urban_area_frac", "Forest_area_frac"),
  Spatial = c("Mod_max_temp_summer", "Z"),
  Morphology = "Ecomorphology_0_12",
  Hydrological = c("Flow_velocity", "Avg_mod_discharge_m3_s"),
  Temporal = c(#"PP_mean_mean"
  ))

pca_labels <- c(
  "TU_ECmix_median"      = "TU-EC (mixture, median)",
  "ARQmix_median"        = "ARQ (mixture, median)",
  "CRQmix_max"           = "CRQ (mixture, max)",
  "TU_NOECmix_max"       = "TU-NOEC (mixture, max)",
  "Agricultural_area_frac" = "Agricultural proportion",
  "Urban_area_frac"      = "Urban proportion",
  "Forest_area_frac"     = "Forest proportion",
  "Mod_max_temp_summer"  = "Max summer temp.",
  "Z"                    = "Altitude",
  "Ecomorphology_0_12"   = "Ecomorphology",
  "Flow_velocity"        = "Flow velocity",
  "Avg_mod_discharge_m3_s" = "Avg. discharge")


## Function to generate the vector plots per scenario and sampling period
## This includes magnitudes and loadings in the vector plot

# WITH ASSUMPTION PLOTS FIXED
prepare_pca_plot <- function(df, scenario_input, var_groups, arrow_scale = 2) {
  
  # Define a custom color palette for variable groups
  palette_pca_bright <- c("#339CFF", "#FF5C5C", "#444444", "forestgreen", "#FFD84A", "#B97AFF")
  
  # Filter the input data frame to include only rows matching the selected scenario
  df_scen <- df %>% 
    filter(trimws(as.character(scenario)) == trimws(as.character(scenario_input)))
  
  # Print the scenario and number of rows after filtering (for quick feedback)
  print(paste("Scenario:", scenario_input, "n:", nrow(df_scen)))
  
  # Select only the predictor variables (columns) as specified in var_groups
  predictor_vars <- df_scen %>%
    select(all_of(unname(unlist(var_groups))))
  
  # Convert all columns from matrices (if present) to vectors
  predictor_vars <- as.data.frame(lapply(predictor_vars, function(x) as.vector(x)))
  
  # Perform Principal Component Analysis (PCA) on the selected variables
  # Note: centering and scaling are set to FALSE (assumes data is already processed)
  pca_res <- prcomp(predictor_vars, center = FALSE, scale. = FALSE)
  
  # Extract loadings (eigenvectors) for the first two principal components (PC1 and PC2)
  loadings <- as.data.frame(pca_res$rotation[, 1:2])
  loadings$varname_orig <- rownames(loadings)  # Add variable names for labeling
  loadings$varname <- pca_labels[loadings$varname_orig] 
  
  
  # Assign each variable to its group (as defined in var_groups) for coloring
  loadings$Dimension <- purrr::map_chr(loadings$varname_orig, function(v) {
    nm <- names(var_groups)[purrr::map_lgl(var_groups, ~ v %in% .x)]
    if(length(nm) == 0) "Other" else nm
  })
  
  # Calculate the magnitude (length) of each loading vector for alpha (transparency)
  loadings$Magnitude <- sqrt(loadings$PC1^2 + loadings$PC2^2)
  loadings$Magnitude <- (loadings$Magnitude - min(loadings$Magnitude)) / 
    (max(loadings$Magnitude) - min(loadings$Magnitude))
  
  # Scale the loading vectors by arrow_scale to control arrow length in the plot
  loadings <- loadings %>%
    mutate(PC1 = PC1 * arrow_scale,
           PC2 = PC2 * arrow_scale)
  
  # Create the PCA loading plot using ggplot2
  p_loadings <- ggplot() +
    geom_segment(data = loadings,
                 aes(x = 0, y = 0, 
                     xend = PC1, 
                     yend = PC2, 
                     color = Dimension,
                     alpha = Magnitude),
                 arrow = arrow(length = unit(0.2, "cm")),
                 linewidth = 1, 
                 show.legend = TRUE) +
    ggrepel::geom_text_repel(data = loadings,
                             aes(x = PC1, y = PC2, 
                                 label = varname, 
                                 color = Dimension),
                             size = 3, 
                             show.legend = FALSE, 
                             max.overlaps = Inf, 
                             force = 2) +
    labs(title = paste("PCA vectors -", scenario_input),
         x = paste0("PC1 (", round(summary(pca_res)$importance[2,1]*100, 1), "%)"),
         y = paste0("PC2 (", round(summary(pca_res)$importance[2,2]*100, 1), "%)")) +
    scale_color_manual(values = palette_pca_bright) +
    theme_bw() +
    guides(linewidth = "none") +
    theme(
      legend.direction = "vertical",
      legend.position = "right")
  
  # --- Diagnostic plots for PCA assumptions ---
  
  # Pairplot (scatterplot matrix) for linearity and correlation
  p_pairs <- GGally::ggpairs(predictor_vars) +
    labs(title = paste("Pairplot -", scenario_input))
  
  # Histograms for normality
  p_hist <- predictor_vars %>%
    pivot_longer(everything()) %>%
    ggplot(aes(x = value)) +
    geom_histogram(bins = 20) +
    facet_wrap(~name, scales = "free") +
    labs(title = paste("Histograms -", scenario_input))
  
  # Boxplots for outliers
  p_box <- predictor_vars %>%
    pivot_longer(everything()) %>%
    ggplot(aes(y = value)) +
    geom_boxplot() +
    facet_wrap(~name, scales = "free") +
    labs(title = paste("Boxplots -", scenario_input))
  
  # Return a list with all plots
  list(
    loadings_plot = p_loadings,
    pairplot = p_pairs,
    histograms = p_hist,
    boxplots = p_box
  )
}


## Selection of the scenarios based on frequencies
scenarios_PCA3.5d <- c("3.5_days_1_week", "3.5_days_2_months", "3.5_days_1_year")
scenarios_PCA14d <- c("14_days_1_week", "14_days_2_months", "14_days_1_year")

## Application of the function
pca_plots3.5 <- map(scenarios_PCA3.5d,
                    ~ prepare_pca_plot(metrics_3.5d_0_R12_abiotic_ds, .x,
                                       var_groups_3.5d, arrow_scale = 2))
pca_plots14 <- map(scenarios_PCA14d,
                   ~ prepare_pca_plot(metrics_14d_0_R12_abiotic_ds, .x,
                                      var_groups_14d, arrow_scale = 2))

# Extract loading plots from each scenario
loading_plots3.5 <- map(pca_plots3.5, "loadings_plot")
loading_plots14 <- map(pca_plots14, "loadings_plot")

# Arrange and export of the PCA plots
#png(paste(user, output_graphics_file ,"PCA_Comparative2.png", sep = ""), 
#    width = 32, height = 20, units = "cm", res = 300)
((loading_plots3.5[[1]] | loading_plots3.5[[2]] | loading_plots3.5[[3]]) /
    (loading_plots14[[1]] | loading_plots14[[2]] | loading_plots14[[3]])) + 
  plot_annotation(
    title = "Principal Component Analysis (PCA) per sampling period in selected time windows",
    caption = 
    "Each row corresponds the PCA per sampling period. Each column is representative of a distinct time window")+
  plot_layout(guides = "collect")+
  theme(plot.title = element_text(hjust = 0, face = "bold")) -> a
#dev.off()
#ggsave("PCA_improved.png", a, width = 30, height = 20, units = "cm", dpi = 350)
#rm(a)

## 2.1. Checking assumptions of the PCA ----
## Using the same object created before (pca_plots3.5 and pca_plots14)
# For example, to see all diagnostic plots for the first 3.5-day scenario (time window):
pca_plots3.5[[1]]$pairplot
pca_plots3.5[[1]]$histograms
pca_plots3.5[[1]]$boxplots


# 3. Linear Regression in the 3.5-days sampling period -----
## 3.1. LMM 3.5-days "bio-indicator ~ chem metrics + (1|YEAR)" -----

# Purpose of the function: Fits linear mixed models (LMMs) for each combination 
# of scenario and chemical metric, extracts model summaries, and filters results 
# by p-value threshold (0.1)
fit_extract_lmm <- function(df, bioind, chem_metrics, 
                            random, pval_threshold = 0.1) {
  #____________________________________________________________
  # Inputs:
  # df: Data frame containing the data
  # bioind: Name of the biological indicator variable (string)
  # chem_metrics: Vector of chemical metric names (strings)
  # random: Name of the random effect variable (string)
  # pval_threshold: Threshold for filtering significant results (default = 0.1)
  # Outputs: 
  # A data frame of model results, including significance annotations and fit statistics
  #____________________________________________________________
  
  # Generate all combinations of scenario and chemical metric
  combinations <- tidyr::crossing(
    scenario = unique(df$scenario),
    chem_metric = chem_metrics)
  # combinations is a tibble with every scenario paired with every chem_metric
  
  # Internal function for fitting and extracting results for a single scenario and chemical metric
  fit_extract_lmm_inner <- function(scen, chem_metric) {
    # Subset data for the current scenario
    sub_df <- df %>% filter(scenario == scen)
    # Create formula for LMM: bioind ~ chem_metric + (1|random)
    formulation <- as.formula(
      paste(bioind, "~", chem_metric, "+ (1|", random, ")"))
    
    # Fit LMM model with error handling (returns NULL if error occurs)
    mod <- tryCatch(
      lmer(formulation, data = sub_df),
      error = function(e) return(NULL)
    )
    # Skip if model failed to fit
    if (is.null(mod)) return(NULL)
    # Return a tibble with summary statistics and diagnostics
    tibble(
      response = bioind,          # Name of the response variable
      chem_metric = chem_metric,  # Name of the chemical metric
      model = "LMM",              # Type of model
      random = random,            # Name of random effect
      scenario = scen,            # Scenario
      converged = is.null(mod@optinfo$conv$lme4$messages), # Did model converge?
      singular = isSingular(mod), # Is the model singular?
      tidy = list(tidy(mod, effects = "fixed", conf.int=TRUE)), # Tidy model summary
      AIC = AIC(mod),             # Akaike Information Criterion
      BIC = BIC(mod),             # Bayesian Information Criterion
      obs = nobs(mod),            # Number of observations
      R2marg = performance::r2_nakagawa(mod)$R2_marginal,   # Marginal R-squared
      R2cond = performance::r2_nakagawa(mod)$R2_conditional # Conditional R-squared
    )
  }
  
  # Map over all combinations and collect results
  results <- pmap_dfr(
    combinations,
    ~fit_extract_lmm_inner(..1, ..2) # ..1 is scenario, ..2 is chemical_metric
  )
  
  # Filter and annotate significance
  results %>%
    unnest(tidy) %>% # Expand tidy model results
    group_by(scenario, chem_metric) %>%
    filter(all(p.value <= pval_threshold)) %>% # Keep only results below p-value threshold
    mutate(Significance =
             case_when(
               p.value <= 0.05 ~ "Below 0.05",
               p.value > 0.05 & p.value <= 0.1 ~ "Marginal significance",
               TRUE ~ NA_character_
             )) %>%
    ungroup() -> results2
  
  # Print final results for inspection
  print(results2)
}

### 3.1.1.Forest plot: SPEAR; GI; VT; IBCH; EPT - LMM 3.5d -----
# Joining the results from the "fit_extract_lmm" function 
bind_rows(
  fit_extract_lmm(
    df = metrics_3.5d_0_R12_abiotic_ds,
    bioind = "GI_VALUE",
    chem_metrics = chemical_metrics_3.5d,
    random = "YEAR") ,
  fit_extract_lmm(
    df = metrics_3.5d_0_R12_abiotic_ds,
    bioind = "SPEAR",
    chem_metrics = chemical_metrics_3.5d,
    random = "YEAR"),
  fit_extract_lmm(
    df = metrics_3.5d_0_R12_abiotic_ds,
    bioind = "IBCH_2019",
    chem_metrics = chemical_metrics_3.5d,
    random = "YEAR"),
  fit_extract_lmm(
    df = metrics_3.5d_0_R12_abiotic_ds,
    bioind = "EPT",
    chem_metrics = chemical_metrics_3.5d,
    random = "YEAR")) %>% 
  # Generation of labels for the plot
  mutate(
    label = paste(scenario, term, sep = " | "),
    Significance = as.character(Significance),
    singular = as.factor(singular)) %>% 
  # Generation of text for the plot
  mutate(
    AIC = as.numeric(AIC),
    BIC = as.numeric(BIC),
    obs = as.numeric(obs),
    info_label = paste0(
      #"AIC: ", round(AIC, 1),
      #"\nBIC: ", round(BIC, 1),
      "n: ", obs)) %>% 
  # Dataframe edition to remark the comparison between chemical metric family/approaches/representation
  # DATAFRAME EDITION TO REMARK THE COMPARSION BETWEEN CHEMICAL METRICS APPROACHES
  mutate(
    chem_metric_fam = case_when(
      grepl("RQ", chem_metric) ~ "Risk Quotient",
      grepl("TU", chem_metric) ~ "Toxic Unit",
      TRUE ~ NA_character_),
    chem_metric_app = case_when(
      grepl("mix_", chem_metric) ~ "Mixture",
      grepl("max_", chem_metric) ~ "Maximum",
      TRUE ~ NA_character_),
    chem_metric_stat = case_when(
      grepl("_max", chem_metric) ~ "Maximum",
      grepl("_mean", chem_metric) ~ "Mean",
      grepl("_median", chem_metric) ~ "Median",
      TRUE ~ NA_character_)) %>% 
  mutate(
    label2 = paste(scenario, chem_metric_stat, sep = " | ")) %>% 
  # Since we only have one predictor, we can remove the intercept to only plot the slope of the model
  # Also, all the intercepts established their baseline around 20 for SPEAR.
  filter(!term %in% "(Intercept)") %>% 
  
  # Plot
  ggplot(aes( x = estimate,
              y = reorder(label2, estimate),
              xmin = conf.low,
              xmax = conf.high,
              color = Significance,
              shape = as.factor(singular))) +
  geom_point(aes(size = R2marg)) +
  geom_errorbarh(height = 0.2) +
  geom_text(aes(x = Inf, label = info_label),
            hjust = 1.5, vjust = 0.5, size = 2,
            color = "black")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  facet_grid(chem_metric_app~chem_metric_fam+response, 
             scales = "free") +
  scale_color_manual(values = c("Below 0.05" = "navyblue", 
                                "Marginal significance" = "coral", 
                                na.value = "gray70")) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 1)) +
  labs(
    x = "Estimate (with 95% CI)",
    y = "Time Windows | Chemical Metric",
    color = "Significance",
    shape = "Singular fit",
    title = "LMM Effects Across Time windows and Chemical Metrics",
    subtitle = "Risk Quotients & Toxic Units | Bioindicators ",
    tag = "A") +
  theme_bw() +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 8),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold"),
    plot.caption = element_text(face = "italic", hjust = 0, size = 9)) -> a

## 3.2. LM 3.5-days "bio-indicator ~ chem metrics" -----
## Since many singularities were found, the same idea of the LMMs function
## is applied for LMs

# Purpose of the function: Fits linear models (LMs) for each combination of 
# scenario and chemical metric, extracts model summaries, and filters results 
# by p-value threshold (0.1)
fit_extract_lm <- function(df, bioind, chem_metrics, 
                           pval_threshold = 0.1) {
  #_______________________________________________________________________________
  # Inputs:
  #   df: Data frame with data and a 'scenario' column
  #   bioind: Name of the biological indicator variable (string)
  #   chem_metrics: Vector of chemical metric names (strings)
  #   pval_threshold: Threshold for filtering significant results (default = 0.1)
  #
  # Output: Data frame of model results, including significance annotations
  #_______________________________________________________________________________
  
  # Generate all combinations of scenario and chemical metric
  combinations <- tidyr::crossing(
    scenario = unique(df$scenario),
    chem_metric = chem_metrics
  )
  
  # Internal function for fitting and extracting results for a single scenario and chemical metric
  fit_extract_lm_inner <- function(scen, chem_metric) {
    # Subset data for current scenario
    sub_df <- df %>% filter(scenario == scen)
    # Create formula: bioind ~ chem_metric
    formulation <- as.formula(paste(bioind, "~", chem_metric))
    # Fit linear model; return NULL if error
    mod <- tryCatch(
      lm(formulation, data = sub_df),
      error = function(e) return(NULL)
    )
    if (is.null(mod)) return(NULL)
    # Return results as tibble
    tibble(
      response = bioind,
      chem_metric = chem_metric,
      model = "LM",
      random = NA,
      scenario = scen,
      converged = NA,
      singular = NA,
      tidy = list(broom::tidy(mod, conf.int = TRUE)),
      AIC = AIC(mod),
      BIC = BIC(mod),
      obs = nobs(mod),
      adjR2 = summary(mod)$adj.r.squared
    )
  }
  
  # Map over all combinations and collect results
  results <- purrr::pmap_dfr(
    combinations,
    ~fit_extract_lm_inner(..1, ..2)
  )
  
  # Filter and annotate significance
  results %>%
    tidyr::unnest(tidy) %>%                # Unpack model results
    group_by(scenario, chem_metric) %>%    # Group by scenario and metric
    filter(all(p.value <= pval_threshold)) %>% # Filter by p-value threshold
    mutate(Significance =                   # Annotate significance
             dplyr::case_when(
               p.value <= 0.05 ~ "Below 0.05",
               p.value > 0.05 & p.value <= 0.1 ~ "Marginal significance",
               TRUE ~ NA_character_
             )) %>%
    ungroup() -> results2                  # Assign final results
  
  print(results2)  # Print results for inspection
}

### 3.2.1. Forest plot: SPEAR; GI; VT; IBCH; EPT - LM 3.5d -----
# Joining the results from the "fit_extract_lm" function 
bind_rows(
  fit_extract_lm(
    df = metrics_3.5d_0_R12_abiotic_ds,
    bioind = "SPEAR",
    chem_metrics = chemical_metrics_3.5d),
  fit_extract_lm(
    df = metrics_3.5d_0_R12_abiotic_ds,
    bioind = "GI_VALUE",
    chem_metrics = chemical_metrics_3.5d),
  fit_extract_lm(
    df = metrics_3.5d_0_R12_abiotic_ds,
    bioind = "IBCH_2019",
    chem_metrics = chemical_metrics_3.5d),
  fit_extract_lm(
    df = metrics_3.5d_0_R12_abiotic_ds,
    bioind = "EPT",
    chem_metrics = chemical_metrics_3.5d)) %>% 
  # Create a label for y-axis (scenario | term)
  dplyr::mutate(
    label = paste(scenario, term, sep = " | "),
    Significance = as.character(Significance),
    singular = as.factor(singular)) %>% 
  # Generation of labels for the plot
  dplyr::mutate(
    AIC = as.numeric(AIC),
    BIC = as.numeric(BIC),
    obs = as.numeric(obs),
    info_label = paste0(
      #"AIC: ", round(AIC, 1),
      #"\nBIC: ", round(BIC, 1),
      "n: ", obs)) %>% 
  filter(!term %in% "(Intercept)") %>% 
  # DATAFRAME EDITION TO REMARK THE COMPARSION BETWEEN CHEMICAL METRICS APPROACHES
  dplyr::mutate(
    chem_metric_fam = case_when(
      grepl("RQ", chem_metric) ~ "Risk Quotient",
      grepl("TU", chem_metric) ~ "Toxic Unit",
      TRUE ~ NA_character_),
    chem_metric_app = case_when(
      grepl("mix_", chem_metric) ~ "Mixture",
      grepl("max_", chem_metric) ~ "Maximum",
      TRUE ~ NA_character_),
    chem_metric_stat = case_when(
      grepl("_max", chem_metric) ~ "Maximum",
      grepl("_mean", chem_metric) ~ "Mean",
      grepl("_median", chem_metric) ~ "Median",
      TRUE ~ NA_character_)) %>% 
  dplyr::mutate(
    label2 = paste(scenario, chem_metric_stat, sep = " | ")) %>% 
  ggplot(aes( x = estimate,
              y = reorder(label2, estimate),
              xmin = conf.low,
              xmax = conf.high,
              color = Significance)) +
  geom_point(aes(size = adjR2)) +
  geom_errorbarh(height = 0.2) +
  geom_text(aes(x = Inf, label = info_label),
            hjust = 1.5, vjust = 0.5, size = 2,
            color = "black")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  facet_grid(chem_metric_app~chem_metric_fam+response, 
             scales = "free") +
  scale_color_manual(values = c("Below 0.05" = "navyblue", 
                                "Marginal significance" = "coral", 
                                na.value = "gray70")) +
  #scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 1)) +
  labs(
    x = "Estimate (with 95% CI)",
    y = "Time Windows | Chemical Metric",
    color = "Significance",
    size = "adjR2",
    title = "LM Effects Across Time windows and Chemical Metrics",
    subtitle = "Risk Quotients & Toxic Units | Bioindicators ",
    tag = "B") +
  theme_bw() +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 8),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold"),
    plot.caption = element_text(face = "italic", hjust = 0, size = 9)) -> b

## 3.3. Joint LMM and LM forest plots - 3.5d ----
#png(paste(user, output_graphics_file ,"Fig.S2_LMM&LM_3.5d.png", sep = ""), 
#    width = 45, height = 35, units = "cm", res = 350)
(a/b)+
  plot_annotation(
    caption = "A: Linear Mixed Models using YEAR as random effect. B: Linear Model.
    A, B: Each estimate represents a model based on one chemical metric, one bioindicator, and one subset of data. Each chemical metric was used as an individual predictor against each bioindicator across time windows. Each chemical metric is log1p-transformed and scaled.
    Facet columns represent the chemical metric family (Risk Quotient or Toxic Unit), facet rows show the aggregation approach (Mixture or Maximum), and the y-axis displays the statistical representation (Mean, Median, Maximum) of each chemical metric.
    If an indicator, time window, or chemical metric interaction does not appear, it was not at least marginally significant (p < 0.1)."
  )
#dev.off()
rm(a,b)
# NO VT OR EPT WAS AT LEAST MARGINALLY SIGNIFICANT*

## 3.4. Forest plot to join only the most relevant bioindicators 3.5d (SPEAR and GI-value) ----

### Dataset "chemical_metrics_diff_LM_3.5ds" generation
bind_rows(
  fit_extract_lm(
    df = metrics_3.5d_0_R12_abiotic_ds,
    bioind = "SPEAR",
    chem_metrics = chemical_metrics_3.5d),
  fit_extract_lm(
    df = metrics_3.5d_0_R12_abiotic_ds,
    bioind = "GI_VALUE",
    chem_metrics = chemical_metrics_3.5d)) %>% 
  # INCLUSION OF THE SD PER TERM
  left_join(sd_term_3.5d, by = c("scenario", "term")) %>% 
  # THIS IS TO CHECK THE RANGE VALUES FOR THE BIOINDICATORS, TO INTERPRET
  #group_by(response, scenario, term) %>% 
  #summarise(max = max(estimate), min = min(estimate)) %>% View()
  # Create a label for y-axis (scenario | term)
  mutate(
    label = paste(scenario, term, sep = " | "),
    Significance = as.character(Significance),
    singular = as.factor(singular)) %>% 
  mutate(
    AIC = as.numeric(AIC),
    BIC = as.numeric(BIC),
    obs = as.numeric(obs),
    info_label = paste0(
      #"AIC: ", round(AIC, 1),
      #"\nBIC: ", round(BIC, 1),
      "n: ", obs)) %>% 
  # Dataframe edition to remark the comparison between chemical metrics approaches
  mutate(
    chem_metric_fam = case_when(
      grepl("RQ", chem_metric) ~ "Risk Quotient",
      grepl("TU", chem_metric) ~ "Toxic Unit",
      TRUE ~ NA_character_),
    chem_metric_app = case_when(
      grepl("mix_", chem_metric) ~ "Mixture",
      grepl("max_", chem_metric) ~ "Maximum",
      TRUE ~ NA_character_),
    chem_metric_stat = case_when(
      grepl("_max", chem_metric) ~ "Maximum",
      grepl("_mean", chem_metric) ~ "Mean",
      grepl("_median", chem_metric) ~ "Median",
      TRUE ~ NA_character_)) %>% 
  mutate(
    label2 = paste(scenario, chem_metric_stat, sep = " | ")) %>% 
  filter(!term %in% "(Intercept)") %>% 
  filter(Significance %in% c("Below 0.05",
                             "Marginal significance")) %>% 
  
  #_______________________________________________________
  # This section was to explore the interpretation of the slopes
  # and baselines on each relevant LM, so it could be removed
  # in case there is no interest for interpretation at this level.
  
  group_by(response, scenario,
           # since i want to keep all the family and approaches, let's 
           # get the means grouping by all of them seeing the difference 
           # based on stats
           chem_metric_fam, chem_metric_app, chem_metric_stat) %>%
  #summarise(
  # ESTIMATES BY SD INTERPRETATIONS
  #  mean_estimate = mean(estimate),
  # mean_CI_low = mean(conf.low),
  # mean_CI_high = mean(conf.high),
  # ESTIMATES BY RAW UNITS INTERPRETATIONS
  #mean_estimate_raw = mean(estimate/SD),
  #mean_CI_low_raw = mean(conf.low/SD),
  #mean_CI_high_raw = mean(conf.high/SD),
  #R2 DOESNT CHANGE
  #mean_adjR2 = mean(adjR2)) %>%
  
  #________________________________________________________
  
  # Ensure all combinations are present for plotting
  right_join(
    expand_grid(
      response = c("SPEAR", "GI_VALUE"),
      chem_metric_fam = c("Risk Quotient", "Toxic Unit"),
      chem_metric_app = c("Mixture", "Maximum"),
      chem_metric_stat = c("Median", "Maximum", "Mean"),
      scenario = unique(.$scenario)),
    by = c("response", "chem_metric_fam", 
           "chem_metric_app", 
           "chem_metric_stat", "scenario")) -> chemical_metrics_diff_LM_3.5ds

### Counts to report
chemical_metrics_diff_LM_3.5ds %>%
  group_by(response, chem_metric_fam) %>%
  summarise(
    min_estimate = round(min(estimate, na.rm = TRUE),3),
    min_conf.low = round(conf.low[which.min(estimate)],3),
    min_conf.high = round(conf.high[which.min(estimate)],3),
    max_estimate = round(max(estimate, na.rm = TRUE),3),
    max_conf.low = round(conf.low[which.max(estimate)],3),
    max_conf.high = round(conf.high[which.max(estimate)],3),
    ajdR2max = max(adjR2, na.rm = TRUE),
    ajdR2min = min(adjR2, na.rm = TRUE))

### Counts to report
chemical_metrics_diff_LM_3.5ds %>%
  filter(response == "EPT") %>%
  ungroup() %>% droplevels() %>% 
  mutate(range_group = as.numeric(cut_number(adjR2, 4))) %>%
  group_by(range_group) %>%
  summarise(
    min = min(adjR2),
    max = max(adjR2))

### Plot SPEAR
ggplot(
  # Here you can change for the bioindicator you want to explore 
  chemical_metrics_diff_LM_3.5ds %>%
    filter(response %in% "SPEAR") %>% 
    # To remove the long name of each scenario (time window in the report)
    mutate(scenario = str_remove(scenario, "^3\\.5_days[_]?")) %>%  
    # Organization of the time windows for the plot
    mutate(
      scenario = factor(
        scenario,
        levels = c("1_week", "2_weeks", "1_month", "2_months", 
                   "3_months", "6_months", "1_year"), ordered = TRUE)) %>% 
    # First week is removed because "fit" values are good but the visual check show no
    # plausible results
    filter(!scenario == "1_week"),
  aes(x = chem_metric_stat))+
  geom_point(aes(y = estimate, 
                 size = adjR2)) +
  geom_errorbar(aes(ymin = conf.low, 
                    ymax = conf.high), 
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "coral") +
  #scale_color_viridis(option = "viridis", na.value = "gray") +
  facet_grid(
    chem_metric_fam + chem_metric_app ~ scenario , 
    scales = "free") +
  scale_size_continuous(
    name = "Adjusted R²",
    breaks = c(0.03, 0.05, 0.07, 0.15),
    limits = c(0.03,0.15),
    range = c(2, 6))+
  labs(
    x = "Statistical Representation",
    y = "Estimate (Slope)",
    #color = "Mean Estimate",
    size = "Adjusted R²",
    subtitle = "Linear Models: SPEAR ~ Chemical metrics (3.5 days)"
    #title = "Mean Estimate, Confidence Intervals, and Adjusted R² by Grouping, Response, and Category"
  ) +
  theme_bw() +
  theme(
    legend.position = "none",
    legend.direction = "none",
    strip.background = element_rect(fill = "white"),
    axis.text.x = element_text(angle = 45, 
                               hjust = 1),
    strip.text.y = element_text(size = 8)) -> a

### Plot GI_VALUE
ggplot(
  # Here you can change for the bioindicator you want to explore 
  chemical_metrics_diff_LM_3.5ds %>%
    filter(response %in% "GI_VALUE") %>% 
    # To remove the long name of each scenario (time window in the report)
    mutate(scenario = str_remove(scenario, "^3\\.5_days[_]?")) %>% 
    # Organization of the time windows for the plot
    mutate(
      scenario = factor(
        scenario,
        levels = c("1_week", "2_weeks", "1_month", "2_months", 
                   "3_months", "6_months", "1_year"), ordered = TRUE)) %>% 
    # First week is removed because "fit" values are good but the visual check show no
    # plausible results
    filter(!scenario == "1_week"), 
  aes(x = chem_metric_stat))+
  geom_point(aes(y = estimate, 
                 size = adjR2)) +
  geom_errorbar(aes(ymin = conf.low, 
                    ymax = conf.high), 
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "coral") +
  #scale_color_viridis(option = "viridis", na.value = "gray") +
  facet_grid(
    chem_metric_fam + chem_metric_app ~ scenario , 
    scales = "free") +
  scale_size_continuous(
    name = "Adjusted R²",
    breaks = c(0.03, 0.05, 0.07, 0.15),
    limits = c(0.03,0.15),
    range = c(2, 6))+
  labs(
    x = "Statistical Representation",
    y = "Estimate (Slope)",
    #color = "Mean Estimate",
    size = "Adjusted R²",
    subtitle = "Linear Models: GI_VALUE ~ Chemical metrics (3.5 days)"
    #title = "Mean Estimate, Confidence Intervals, and Adjusted R² by Grouping, Response, and Category"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    strip.background = element_rect(fill = "white"),
    axis.text.x = element_text(angle = 45, 
                               hjust = 1),
    strip.text.y = element_text(size = 8)) -> b

### Merge and export the plot
#png(paste(user, output_graphics_file ,"LM_3.5d.png", sep = ""),
#    "Estimates_LM_3.5d.png", width = 35, height = 20, units = "cm", res = 350)
a + b +
  plot_annotation(
    title = "Estimate, Confidence Intervals, and Adjusted R² by grouping, Response, and time windows",
    caption = 
    "Each chemical metric was lop1p-transformed and scaled.
    Facet columns represent the time windows, facet rows show the aggregation approach (Mixture or Maximum) of each chemical metric,
    and the x-axis displays the statistical representation (Mean, Median, Maximum) of each chemical metric.")+
  plot_layout(guides = "collect") &
  theme( legend.position = "bottom",
         legend.direction = "horizontal",
         plot.title = element_text(hjust = 0, face = "bold")) -> c
c
rm(a,b,c)
#ggsave("LM_3.5d_improved.png", c, width = 25, height = 17, units = "cm", dpi = 350)
##

## 3.5. Checking assumptions (only linearity and cook's distance) ----

# Purpose: 
## All the assumptions were first checked using the native plot() function for each model.
## Once linearity and outliers were found to be the least met assumptions, a specific plot was
## generated for reporting purposes.

plot_lm_diagnostics_base <- function(df, response, predictor, main_title = NULL) {
  form <- as.formula(paste(response, "~", predictor))
  modelo <- lm(form, data = df)
  cooks <- cooks.distance(modelo)
  
  # 1. Residuals vs Fitted
  plot(modelo$fitted.values, resid(modelo),
       main = ifelse(is.null(main_title), "Residuals vs Fitted", paste(main_title, "\nResiduals vs Fitted")),
       xlab = "Fitted values", ylab = "Residuals",
       pch = 19, col = "steelblue")
  abline(h = 0, lty = 2, col = "gray40")
  
  # 2. Cook's Distance
  plot(cooks, type = "h",
       main = ifelse(is.null(main_title), "Cook's Distance", paste(main_title, "\nCook's Distance")),
       xlab = "Observation", ylab = "Cook's D",
       col = "pink", lwd = 2)
  abline(h = 0.5, lty = 2, col = "coral")
  abline(h = 1, lty = 2, col = "red")
}

# Now, we can plot the desired response and predictor based on the scenario by filtering the
# dataset
#png(paste(user, output_graphics_file ,"Model_diagnostics_LM_SPEAR_CHEM-3.5d.png", sep = ""), 
#    width = 40, height = 30, units = "cm", res = 350)

par(mfrow = c(4, 6)) # rows, columns

# Row 1, each function gives you 2 plots, so 6 columns here
plot_lm_diagnostics_base(metrics_3.5d_0_R12_abiotic_ds %>% 
                           filter(scenario == "3.5_days_1_week"), 
                         response = "SPEAR", predictor = "TU_ECmix_median", 
                         main_title = "3.5_days_1_week|TU_ECmix_median")
plot_lm_diagnostics_base(metrics_3.5d_0_R12_abiotic_ds %>% 
                           filter(scenario == "3.5_days_2_months"), 
                         response = "SPEAR", predictor = "TU_ECmix_median", 
                         main_title = "3.5_days_2_months|TU_ECmix_median")
plot_lm_diagnostics_base(metrics_3.5d_0_R12_abiotic_ds %>% 
                           filter(scenario == "3.5_days_1_year"), 
                         response = "SPEAR", predictor = "TU_ECmix_median", 
                         main_title = "3.5_days_1_year|TU_ECmix_median")

# Row 2, each function gives you 2 plots, so 6 columns here
plot_lm_diagnostics_base(metrics_3.5d_0_R12_abiotic_ds %>% 
                           filter(scenario == "3.5_days_1_week"), 
                         response = "SPEAR", predictor = "TU_ECmix_max", 
                         main_title = "3.5_days_1_week|TU_ECmix_max")
plot_lm_diagnostics_base(metrics_3.5d_0_R12_abiotic_ds %>% 
                           filter(scenario == "3.5_days_2_months"), 
                         response = "SPEAR", predictor = "TU_ECmix_max", 
                         main_title = "3.5_days_2_months|TU_ECmix_max")
plot_lm_diagnostics_base(metrics_3.5d_0_R12_abiotic_ds %>% 
                           filter(scenario == "3.5_days_1_year"), 
                         response = "SPEAR", predictor = "TU_ECmix_max", 
                         main_title = "3.5_days_1_year|TU_ECmix_max")

# Row 3, each function gives you 2 plots, so 6 columns here
plot_lm_diagnostics_base(metrics_3.5d_0_R12_abiotic_ds %>% 
                           filter(scenario == "3.5_days_1_week"), 
                         response = "SPEAR", predictor = "ARQmix_median", 
                         main_title = "3.5_days_1_week|ARQmix_median")
plot_lm_diagnostics_base(metrics_3.5d_0_R12_abiotic_ds %>% 
                           filter(scenario == "3.5_days_2_months"), 
                         response = "SPEAR", predictor = "ARQmix_median", 
                         main_title = "3.5_days_2_months|ARQmix_median")
plot_lm_diagnostics_base(metrics_3.5d_0_R12_abiotic_ds %>% 
                           filter(scenario == "3.5_days_1_year"), 
                         response = "SPEAR", predictor = "ARQmix_median", 
                         main_title = "3.5_days_1_year|ARQmix_median")

# Row 4, each function gives you 2 plots, so 6 columns here
plot_lm_diagnostics_base(metrics_3.5d_0_R12_abiotic_ds %>% 
                           filter(scenario == "3.5_days_1_week"), 
                         response = "SPEAR", predictor = "ARQmix_max", 
                         main_title = "3.5_days_1_week|ARQmix_max")
plot_lm_diagnostics_base(metrics_3.5d_0_R12_abiotic_ds %>% 
                           filter(scenario == "3.5_days_2_months"), 
                         response = "SPEAR", predictor = "ARQmix_max", 
                         main_title = "3.5_days_2_months|ARQmix_max")
plot_lm_diagnostics_base(metrics_3.5d_0_R12_abiotic_ds %>% 
                           filter(scenario == "3.5_days_1_year"), 
                         response = "SPEAR", predictor = "ARQmix_max", 
                         main_title = "3.5_days_1_year|ARQmix_max")
par(mfrow = c(1, 1))
#dev.off()

# 4. GAM 3.5d + "fit_extract_lm" function adapted to GAM ----
## GAMs are tried as alternative due to the modest  R2 of the LM. Additionally, some assumptions
## were partially met.

# Purpose of the function: Fits generalized additive models (GAMs) for each combination of 
# scenario and chemical metric, extracts model summaries, and filters results 
# by p-value threshold (0.1); (k=6)
fit_extract_gam <- function(df, bioind, chem_metrics, pval_threshold = 0.1) {
  # _______________________________________________________________________
  # Inputs:
  # df: Data frame with data and 'scenario' column
  # bioind: Name of biological indicator (string)
  # chem_metrics: Vector of chemical metric names (strings)
  # pval_threshold: Significance threshold (default = 0.1)
  #
  # Output: Data frame of GAM results with significance annotations
  # _______________________________________________________________________
  
  # Generate all scenario-metric combinations
  combinations <- tidyr::crossing(
    scenario = unique(df$scenario),
    chem_metric = chem_metrics)
  
  # Internal function for GAM fitting and extraction
  fit_extract_gam_inner <- function(scen, chem_metric) {
    # Subset data for current scenario
    sub_df <- df %>% filter(scenario == scen)
    # Create GAM formula with smooth term (k=6)
    formulation <- as.formula(paste(bioind, "~ s(", chem_metric, ", k = 6)", sep = ""))
    
    # Fit GAM model with error handling
    mod <- tryCatch(
      mgcv::gam(formulation, data = sub_df, method = "REML"), 
      error = function(e) return(NULL))
    
    if (is.null(mod)) return(NULL)
    
    # Extract model summary
    gam_sum <- summary(mod)
    smooth_stats <- as.data.frame(gam_sum$s.table)
    
    # Extract parametric coefficients and CIs
    coef_df <- data.frame(
      term = rownames(gam_sum$p.table),
      estimate = gam_sum$p.table[, 1],
      conf.low = gam_sum$p.table[, 1] - 1.96 * gam_sum$p.table[, 2],
      conf.high = gam_sum$p.table[, 1] + 1.96 * gam_sum$p.table[, 2]
    )
    
    # Return results as tibble
    tibble(
      response = bioind,
      chem_metric = chem_metric,
      model = "GAM",
      scenario = scen,
      converged = mod$converged,
      AIC = AIC(mod),
      BIC = BIC(mod),
      obs = nobs(mod),
      adjR2 = gam_sum$r.sq,
      deviance_explained = gam_sum$dev.expl,
      smooth_edf = smooth_stats$edf[1],
      smooth_pvalue = smooth_stats$`p-value`[1],
      tidy = list(coef_df)  # Store coefficients
    )
  }
  
  # Process all combinations
  results <- purrr::pmap_dfr(
    combinations,
    ~fit_extract_gam_inner(..1, ..2)
  )
  
  # Defensive check for empty results
  if (nrow(results) == 0 || !all(c("scenario", "chem_metric") %in% names(results))) {
    warning("No valid models or missing columns.")
    return(invisible(NULL))
  }
  
  # Filter and annotate significance based on smooth term p-value
  results %>% 
    group_by(scenario, chem_metric) %>%
    filter(!is.na(smooth_pvalue) & smooth_pvalue <= pval_threshold) %>%
    dplyr::mutate(Significance =
                    dplyr::case_when(
                      smooth_pvalue <= 0.05 ~ "Below 0.05",
                      smooth_pvalue > 0.05 & smooth_pvalue <= 0.1 ~ "Marginal significance",
                      TRUE ~ NA_character_
                    )) -> results2
  
  print(results2)  # Print final results
}

# Trying the "best" structure observed in LM
gam(SPEAR ~ s(TU_ECmix_max), 
    data = metrics_3.5d_0_R12_abiotic_ds %>% 
      filter(scenario == "3.5_days_3_months"),
    method = "REML") %>% summary()

# 4.1. LMs and GAMs (overlapped) for exploration ----
##____________________________________________________________________________________________________
## Polished exportation of the plots for all the scenarios (time windows) or only the 
## representative time windows. Be careful with the interpretation here, since the script
## fits LM and GAMs without filter if they are significant or not. So, a carefully evaluation
## must be done in line with the forest plot and assumptions analyzed above.

## This plots, ca be replicated for every bioindicator in case the further exploration is needed
## Now it's set for SPEAR, explaration showed as not meaningful results for the rest of bioindicators
##____________________________________________________________________________________________________

# Filter scenarios and select only the representative time windows
df_long_lm_gams <- metrics_3.5d_0_R12_abiotic_ds %>%
  # Here, you can select the scenarios (time windows) of interest, or leave them inactive
  # so that all the iterations are plotted.
  #filter(scenario %in% c("3.5_days_2_weeks",
  #                       "3.5_days_3_months",
  #                       "3.5_days_1_year")) %>% 
  select(scenario, SPEAR, TU_ECmix_median, ARQmix_median, TU_ECmix_max, ARQmix_max)
  
# Long format of the dataset
df_long_lm_gams <- df_long_lm_gams %>%
  pivot_longer(
    # In line with the chemical metrics selected above
    cols = c(TU_ECmix_median, ARQmix_median, TU_ECmix_max, ARQmix_max),
    names_to = "chem_metric",
    values_to = "chem_value")

# Compute model statistics for each panel
stats_lm_gams <- df_long_lm_gams %>%
  group_by(scenario, chem_metric) %>%
  summarise(
    # LM
    adjR2_LM = tryCatch(signif(summary(lm(SPEAR ~ chem_value))$adj.r.squared, 2), error = function(e) NA),
    AIC_LM = tryCatch(signif(AIC(lm(SPEAR ~ chem_value)), 4), error = function(e) NA),
    p_LM = tryCatch(signif(coef(summary(lm(SPEAR ~ chem_value)))[2,4], 2), error = function(e) NA),
    # GAM 
    adjR2_GAM = tryCatch({
      gam_mod <- gam(SPEAR ~ s(chem_value, k = 6), method = "REML")
      signif(summary(gam_mod)$r.sq, 2)
    }, error = function(e) NA),
    AIC_GAM = tryCatch({
      gam_mod <- gam(SPEAR ~ s(chem_value, k = 6), method = "REML")
      signif(AIC(gam_mod), 4)
    }, error = function(e) NA),
    edf_GAM = tryCatch({
      gam_mod <- gam(SPEAR ~ s(chem_value, k = 6), method = "REML")
      signif(summary(gam_mod)$s.table[1, "edf"], 2)
    }, error = function(e) NA),
    p_GAM = tryCatch({
      gam_mod <- gam(SPEAR ~ s(chem_value, k = 6), method = "REML")
      signif(summary(gam_mod)$s.table[1, "p-value"], 2)
    }, error = function(e) NA),
    .groups = "drop") %>%
  mutate(
    stat_text = paste0(
      "LM:\nAdjR² = ", adjR2_LM,
      "\nAIC = ", AIC_LM,
      "\np = ", p_LM,
      "\n\nGAM:\nAdjR² = ", adjR2_GAM,
      "\nAIC = ", AIC_GAM,
      "\nEDF = ", edf_GAM,
      "\np = ", p_GAM))

# Plot
p <- ggplot(df_long_lm_gams, aes(x = chem_value, y = SPEAR)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray35") +
  geom_smooth(method = "lm", se = TRUE, color = "red", fill = "pink", alpha = 0.35,
              linetype = "dashed", linewidth = 1) +
  geom_smooth(method = "gam", formula = y ~ s(x, k = 6), se = TRUE, 
              color = "blue", fill = "lightblue", alpha = 0.35,
              linewidth = 1) +
  facet_grid(chem_metric ~ scenario, scales = "free") +
  geom_text(
    data = stats_lm_gams,
    aes(x = Inf, y = Inf, label = stat_text),
    hjust = 1, vjust = 1.1, size = 3, fontface = "bold", 
    inherit.aes = FALSE) +
  labs(
    # Be careful with the name of the bioindicator
    title = "SPEAR vs Chemical metrics",
    subtitle = "Time windows: All time windows\nLM (red dashed) and GAM (blue) with model statistics",
    x = "scaled(log1p(Chemical metric))",
    # Be careful with the name of the bioindicator
    y = "GI_VALUE",
    # Be careful with the name of the bioindicator in the caption
    caption = 
    "The horizontal line at 0 indicates the lower limit of SPEAR values. Values below 0 are not plausible as SPEAR cannot be negative.
    Model performance is summarized by R2, AIC (Akaike Information Criterion; lower values indicate better fit), EDF (effective degrees of freedom), and p-value (model significance)."
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),  
    strip.text = element_text(face = "bold"))
p

# Export as PNG
#ggsave(
#  filename = file.path(user, output_graphics_file, "Fig_Supp.SPEAR_3.5d_all_time_windows.png"),
#  plot = p,
#  width = 45, height = 25,
#  units = "cm", dpi = 350)
rm(p)

## 4.1.1.Fast check of the assumptions when LM or GAM ----

# Purpose:
# Summarise plot of the most relevant assumptions, including a ggreppel to identify the most
# influential points, so that a sensitivity analysis can be performed.
plot_lm_gam_diagnostics <- function(df, scenario, bioindicator, chemical_metric) {
  # _____________________________________________________________________________________
  # Inputs:
  # df: Data frame containing the data, with columns including Gewasser, YEAR,
  # scenario, bioindicator, and chemical_metric.
  # scenario: String specifying the scenario to filter.
  # bioindicator: String specifying the name of the biological indicator variable.
  # chemical_metric: String specifying the name of the chemical metric variable.
  #
  # Outputs:
  # A single ggplot object (arranged as a grid) containing LM and GAM diagnostic plots,
  # with a main title and detailed caption describing each panel.
  # ______________________________________________________________________________________
  
  # Add point labels for each observation (Gewasser | YEAR) and unique row IDs
  df <- df %>%
    mutate(point_label = paste0(Gewasser, " | ", YEAR),
           row_id = row_number())
  
  # Subset data for the specified scenario and bioindicator
  df_sub <- df %>%
    filter(
      scenario == !!scenario,
      bioindicator == !!bioindicator
    )
  
  # Fit linear model (LM) for the specified variables
  lm_fit <- lm(as.formula(paste(bioindicator, "~", chemical_metric)), data = df_sub)
  # Fit generalized additive model (GAM) for the specified variables
  gam_fit <- gam(as.formula(paste(bioindicator, "~s(", chemical_metric, ",k=6)")), data = df_sub)
  
  # Extract LM statistics: adjusted R-squared, AIC, and p-value
  lm_adjR2 <- signif(summary(lm_fit)$adj.r.squared, 3)
  lm_AIC <- signif(AIC(lm_fit), 3)
  lm_pval <- signif(summary(lm_fit)$coefficients[2,4], 2)
  
  # Extract GAM statistics: adjusted R-squared, AIC, effective degrees of freedom, and p-value
  gam_adjR2 <- signif(summary(gam_fit)$r.sq, 3)
  gam_AIC <- signif(AIC(gam_fit), 3)
  gam_EDF <- signif(summary(gam_fit)$edf[1], 2)
  gam_pval <- signif(summary(gam_fit)$s.table[1,4], 2)
  
  # Prepare LM diagnostic data: fitted values, residuals, standardized residuals, and row IDs
  lm_diag <- data.frame(
    fitted = fitted(lm_fit),
    resid = resid(lm_fit),
    stdresid = rstandard(lm_fit),
    row_id = df_sub$row_id
  )
  
  # LM: Main scatterplot with fitted line and model statistics
  p_lm1 <- ggplot(df_sub, aes_string(x = chemical_metric, y = bioindicator)) +
    geom_point(alpha = 0.5) +
    geom_text_repel(aes(label = paste(row_id, point_label, sep = "-")), size = 2, max.overlaps = 100) +
    geom_smooth(method = "lm", se = TRUE, color = "blue") +
    theme_bw() +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 2,
             label = paste0("AdjR2: ", lm_adjR2,
                            "\nAIC: ", lm_AIC,
                            "\np: ", lm_pval),
             size = 3) +
    labs(title = "Linear Model (LM)")
  
  # LM: Residuals vs Fitted plot
  p_lm2 <- ggplot(lm_diag, aes(x = fitted, y = resid)) +
    geom_point(alpha = 0.5) +
    geom_text_repel(aes(label = row_id), size = 2, max.overlaps = 100) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw() +
    labs(title = "LM: Residuals vs Fitted")
  
  # LM: QQ plot of standardized residuals
  p_lm3 <- ggplot(lm_diag, aes(sample = stdresid)) +
    stat_qq() + stat_qq_line() +
    geom_text(aes(x = quantile(stdresid, probs = ppoints(length(stdresid))), 
                  y = sort(stdresid), label = row_id), 
              size = 2, check_overlap = TRUE, vjust = -0.5) +
    theme_bw() +
    labs(title = "LM: QQ Plot")
  
  # LM: Cook's Distance plot (using custom or external gg_cooksd function)
  p_lm4 <- gg_cooksd(lm_fit, label = TRUE, show.threshold = TRUE, threshold = "convention", scale.factor = 0.7) +
    theme_bw() +
    labs(title = "LM: Cook's Distance")
  
  # Prepare GAM diagnostic data: fitted values, deviance residuals, linear predictor, standardized residuals, and row IDs
  gam_diag <- data.frame(
    fitted = fitted(gam_fit),
    resid = residuals(gam_fit, type = "deviance"),
    linpred = predict(gam_fit, type = "link"),
    stdresid = residuals(gam_fit, type = "pearson"),
    row_id = df_sub$row_id
  )
  
  # GAM: Main scatterplot with fitted smooth and model statistics
  p_gam1 <- ggplot(df_sub, aes_string(x = chemical_metric, y = bioindicator)) +
    geom_point(alpha = 0.5) +
    geom_text_repel(aes(label = paste(row_id, point_label, sep = "-")), size = 2, max.overlaps = 100) +
    stat_smooth(method = "gam", formula = y ~ s(x), se = TRUE, color = "red") +
    theme_bw() +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 2,
             label = paste0("AdjR2: ", gam_adjR2,
                            "\nAIC: ", gam_AIC,
                            "\nEDF: ", gam_EDF,
                            "\np: ", gam_pval),
             size = 3) +
    labs(title = "Generalized Additive Model (GAM)")
  
  # GAM: Deviance residuals vs Fitted plot
  p_gam2 <- ggplot(gam_diag, aes(x = fitted, y = resid)) +
    geom_point(alpha = 0.5) +
    geom_text_repel(aes(label = row_id), size = 2, max.overlaps = 100) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw() +
    labs(title = "GAM: Deviance Residuals vs Fitted")
  
  # GAM: QQ plot of standardized residuals (no labels for clarity)
  p_gam3 <- ggplot(gam_diag, aes(sample = stdresid)) +
    stat_qq() + stat_qq_line() +
    theme_bw() +
    labs(title = "GAM: QQ Plot")
  
  # GAM: Residuals vs Linear Predictor plot
  p_gam4 <- ggplot(gam_diag, aes(x = linpred, y = resid)) +
    geom_point(alpha = 0.5) +
    geom_text_repel(aes(label = row_id), size = 2, max.overlaps = 100) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_bw() +
    labs(title = "GAM: Residuals vs Linear Predictor")
  
  # Arrange LM diagnostic plots in a 2x2 grid
  lm_grid <- ggarrange(p_lm1, p_lm2, p_lm3, p_lm4, ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"))
  # Arrange GAM diagnostic plots in a 2x2 grid
  gam_grid <- ggarrange(p_gam1, p_gam2, p_gam3, p_gam4, ncol = 2, nrow = 2, labels = c("E", "F", "G", "H"))
  
  # Combine LM and GAM grids into a single vertical layout
  final_plot <- ggarrange(lm_grid, gam_grid, ncol = 1, heights = c(1.2, 1))
  
  # Add a main title and a detailed caption to the final plot
  final_plot <- annotate_figure(
    final_plot,
    top = text_grob(
      paste("Model diagnostics for", bioindicator, "vs", chemical_metric, "(", scenario, ")"),
      face = "bold", size = 14),
    bottom = text_grob(" A: Scatterplot of the response variable against the chemical metric, showing the fitted Linear Model (LM) and labeling each point by \"Gewasser | YEAR\". B: LM diagnostics – Residuals vs Fitted values, with each point labeled by its row number to identify potential outliers.
    C: LM diagnostics – QQ plot of standardized residuals, with each point labeled by its row number to assess normality. D: LM diagnostics – Cook's Distance plot, highlighting influential observations according to the conventional threshold.
    E: Scatterplot of the response variable against the chemical metric, showing the fitted Generalized Additive Model (GAM) and labeling each point by \"Gewasser | YEAR\". F: GAM diagnostics – Deviance residuals vs Fitted values, with each point labeled by its row number.
    G: GAM diagnostics – QQ plot of standardized residuals, with each point labeled by its row number. H: GAM diagnostics – Residuals vs Linear Predictor, with each point labeled by its row number.
    Each tag corresponds to the panel in the figure. \"Gewasser | YEAR\" labels allow identification of individual sampling points. Row numbers in diagnostic plots facilitate cross-referencing with the main scatterplots and the dataset. Cook's Distance highlights influential points for the LM fit.",
                       hjust = 0, x = 0, size = 9
    )
  )
  
  # Return the final plot object
  return(final_plot)
}

# Export as PNG, here you can play with the scenario, bioindicator and chemical metric of interes
#png(paste(user, output_graphics_file ,"Model_diagnostics_SPEAR_ARQmix_median-3.5d.png", sep = ""), 
#    width = 25, height = 35, units = "cm", res = 350)
plot_lm_gam_diagnostics(metrics_3.5d_0_R12_abiotic_ds, 
                        scenario = "3.5_days_2_weeks", 
                        bioindicator = "SPEAR", 
                        chemical_metric = "ARQmix_median")
#dev.off()

# 5. GAM 3.5d: SPEAR ~ chem metrics + environmental gradients ----
## Data-set preparation only with the relevant abiotic factors 
metrics_3.5d_0_R12_abiotic_ds %>% 
  select(
    # Organization columns
    Gewasser, scenario, YEAR,
    # Bioindicators
    SPEAR, GI_VALUE, IBCH_2019, VT_VALUE, EPT,
    # Chemical metrics
    contains("max"), contains("mix"), TU_ECmix_median,
    # Land use
    Agricultural_area_frac,
    Forest_area_frac,
    Urban_area_frac,
    # Hydrological
    Flow_velocity,
    Avg_mod_discharge_m3_s,
    # Morphological
    Ecomorphology_0_12,
    # Temporal
    Mod_max_temp_summer,
    PP_mean_mean,
    # Spatial
    Z) %>% na.omit() -> df_metrics_abiotic_3.5d

#_______________________________________________________________________________
# CONCEPTUALIZATION BASED ON RELEVANT FACTORS FOR SPEAR AND OBSERVED 
# ENVIRONMENTAL GRADIENTS IN THE PCA
# Keep chemical and agric, prioritize flow

# Attempts:
# 1. Chemical + Agri + Flow: BEST ONE
# 2. Flow:  It loses explanatory power and interpretability.
# 3. Agri:  It loses explanatory power and interpretability.
# 4. Agri + Flow: SECOND BEST ONE. Agriculture loses interpretability.
#_______________________________________________________________________________

## 5.1. Summaries of the model attempts -----

#____________________________________________________________________________________
## Here you have to change the name of the bioindicator and chemical metric, so it 
## is not necessary to create more objects and you won't end up with too many models.

## Just be careful, because the object "models_gam" is then used for the partial effect 
## plots. To produce a new plot, you need to run a new model from here. So that
## you can check the main values obtained from the new structure. But also for the
## generation of the table with "make_gam_summary_table" function
#____________________________________________________________________________________

## Let's create the empty lists
gam_models <- list()
gam_summaries <- list()
gam_checks <- list()
gam_vif_results <- list()
gam_aic_vals <- c()

## Loop to get fitted GAMs per scenario (time window)
for (scen in scenarios_3.5d) {
  
  # _____________________________________________________________________________
  # Inputs:
  # scenarios_3.5d: Vector of scenario names to loop over.
  # df_metrics_abiotic_3.5d: Data frame containing the data for each scenario,
  # with columns including scenario, SPEAR, TU_ECmix_median, and possibly others.
  # gam_models: List (initialized outside loop) to store GAM model objects.
  # gam_summaries: List (initialized outside loop) to store GAM model summaries.
  # gam_checks: List (initialized outside loop) to store gam.check output.
  # gam_aic_vals: Numeric vector or named list (initialized outside loop) to store AIC values.
  #
  # Outputs:
  # gam_models: List of fitted GAM models, named by scenario.
  # gam_summaries: List of GAM model summaries, named by scenario.
  # gam_checks: List of gam.check outputs, named by scenario.
  # gam_aic_vals: Named vector/list of AIC values for each GAM model.
  # For each scenario, diagnostic plots are displayed.
  # ----------------------------------------------------------------------------
  
  # Subset the data for the current scenario
  df <- df_metrics_abiotic_3.5d %>% 
    filter(scenario == scen)
  
  # Fit a GAM model for SPEAR using the smooth terms you need
  # (other terms are commented out for now)
  model <- gam(SPEAR ~
                 s(TU_ECmix_median, k = 6) +
               #s(Ecomorphology_0_12)+
               s(Agricultural_area_frac, k = 6)+
               #s(Urban_area_frac, k=6) +
               #s(Mod_max_temp_summer)+
               s(Flow_velocity, k = 6),
               data = df, method = "REML")
  
  # Store the GAM model in the models list, using the scenario name as the key
  gam_models[[scen]] <- model
  # Store the summary of the GAM model in the summaries list
  gam_summaries[[scen]] <- summary(model)
  # Capture the output of gam.check (model diagnostics) and store it in gam_checks
  gam_checks[[scen]] <- capture.output(gam.check(model))
  
  # Fit a linear model (LM) for SPEAR using the predictors in the GAM
  # (other terms are commented out for now)
  lm_model <- lm(SPEAR ~ TU_ECmix_median +
                 #Ecomorphology_0_12 +
                 Agricultural_area_frac +
                 Flow_velocity,
                 data = df)
  
  # (Optional: commented out) Store variance inflation factors (VIF) for the LM
  # this must be activated only when two or more predictors are used
  gam_vif_results[[scen]] <- vif(lm_model)
  
  # Store the AIC of the GAM model in aic_vals, using the scenario name as the key
  gam_aic_vals[scen] <- AIC(model)
  
  # Set up a 2x2 plot layout for the next plots
  par(mfrow = c(2,2))
  # Plot the GAM smooths and diagnostics for the current scenario
  plot(model, pages = 1, main = scen)
}
rm(df, model)

# Extract R2, p-values, edf, F as before:
gam_rvals <- sapply(gam_summaries, function(x) x$r.sq)
gam_pvals <- sapply(gam_summaries, function(x) x$s.table[, "p-value"])
gam_edfs  <- sapply(gam_summaries, function(x) x$s.table[, "edf"])
gam_fs    <- sapply(gam_summaries, function(x) x$s.table[, "F"])

# Print or save diagnostics
print(gam_aic_vals)
print(gam_rvals)
print(gam_pvals)
print(gam_edfs)
print(gam_fs)
print(gam_vif_results)

# Plot individuals
print(gam_checks[["3.5_days_6_months"]])
par(mfrow = c(2,2))
plot(gam_models[["3.5_days_1_year"]]) 
par(mfrow = c(1,1))
dev.off()

## Summary table
make_gam_summary_table <- function(summaries, aic_vals, vif_results = NULL, decimals = 3) {
    # _____________________________________________________________________________________________________________
  # Inputs:
  # gam_summaries:   List of GAM model summary objects 
  # gam_aic_vals:    Named vector (or list) of AIC values for each scenario, with names matching those in 'summaries'
  # gam_vif_results: Optional list of VIF (variance inflation factor) results for each scenario (default = NULL)
  # gam_decimals:    Number of decimals for rounding (default = 3)
  #
  # Outputs:
  # A data frame (tibble) with one row per scenario and columns for:
  # - Scenario name
  # - AIC
  # - R-squared (R2)
  # - Mean VIF (if provided)
  # - EDF (effective degrees of freedom) for each smooth term
  # - F-statistic for each smooth term
  # - p-value for each smooth term
  # _____________________________________________________________________________________________________________
  
  # Get all unique smooth term names across all GAM summaries
  all_smooths <- unique(unlist(lapply(summaries, function(x) rownames(x$s.table))))
  
  # For each scenario, create a row with relevant statistics
  rows <- purrr::map2(names(summaries), summaries, function(scen, summ) {
    # Initialize vectors for EDF, F, and p-values for each smooth term
    edfs <- setNames(rep(NA, length(all_smooths)), paste0("edf_", all_smooths))
    Fs   <- setNames(rep(NA, length(all_smooths)), paste0("F_", all_smooths))
    ps   <- setNames(rep(NA, length(all_smooths)), paste0("p_", all_smooths))
    
    # Fill in the values for each smooth term present in the current summary
    for (sm in rownames(summ$s.table)) {
      edfs[paste0("edf_", sm)] <- round(summ$s.table[sm, "edf"], decimals)
      Fs[paste0("F_", sm)]     <- round(summ$s.table[sm, "F"], decimals)
      ps[paste0("p_", sm)]     <- signif(summ$s.table[sm, "p-value"], decimals)  # Use signif for p-values
    }
    
    # Calculate mean VIF for the scenario if VIF results are provided
    vif_mean <- if (!is.null(vif_results)) round(mean(vif_results[[scen]], na.rm = TRUE), decimals) else NA
    
    # Combine all values into a single-row tibble for the scenario
    tibble(
      Scenario = scen,
      AIC = round(aic_vals[scen], decimals),
      R2 = round(summ$r.sq, decimals),
      VIF_mean = vif_mean,
      !!!as.list(edfs),
      !!!as.list(Fs),
      !!!as.list(ps)
    )
  })
  
  # Combine all scenario rows into one data frame and return it
  bind_rows(rows)
}

## USE: BE CAREFUL THAT IS SAVED UNDER THE LAST STRUCTURE IN THE PREVIOUS "gam_model, gam_summaries" LOOP
make_gam_summary_table(gam_summaries, gam_aic_vals, gam_vif_results) %>% View()

## 5.2. Partial effect plots -----

#____________________________________________________________________________________
## Here you have to change set the list of scenarios and abiotic factors to consider 
## matching the ones saved in the "gam_model" performed above in 5.1.
#____________________________________________________________________________________

# List of scenarios (time windows) and predictors to use for the GAMs
scenarios_3.5dGAM <- c("3.5_days_2_weeks", "3.5_days_3_months", "3.5_days_1_year")
predictors_3.5dGAM <- c("TU_ECmix_median",#, "TU_ECmix_median"
                        #"Ecomorphology_0_12",
                        "Agricultural_area_frac", "Flow_velocity"
                        )
gam_predictor_labels <- c(
  "Agricultural_area_frac" = "Agricultural proportion",
  #"Urban_area_frac" = "Urban proportion",
  "Flow_velocity" = "Flow velocity",
  "TU_ECmix_median" = "TU-EC (mixture, median)")

## Purpose: Generation of the dataset of partial effects to plot, based on the 
## gams generated above
get_partial_pred_gam <- function(model, predictor, n = 100) {
  # _________________________________________________________________________________________
  # Inputs:
  # model: A fitted GAM model (mgcv::gam object)
  # predictor: String specifying the name of the predictor variable to generate partial predictions for
  # n: Number of points to use for the predictor sequence (default = 100)
  #
  # Outputs:
  # A tibble with columns:
  # - x: Sequence of predictor values
  # - fit: Predicted smooth effect
  # - se: Standard error of the smooth effect
  # - lower: Lower bound of 95% confidence interval
  # - upper: Upper bound of 95% confidence interval
  # _________________________________________________________________________________________
  
  # Generate a sequence of predictor values covering its range
  x_seq <- seq(min(model$model[[predictor]], na.rm = TRUE),
               max(model$model[[predictor]], na.rm = TRUE), length.out = n)
  
  # Create a new data frame with means of all predictors (except the one of interest)
  # Note: `predictors_3.5dGAM` should be a vector of predictor names used in the model
  newdata <- as.data.frame(lapply(model$model[, predictors_3.5dGAM], mean, na.rm = TRUE))
  newdata <- newdata[rep(1, n), ]  # Replicate the row to match the length of x_seq
  newdata[[predictor]] <- x_seq     # Replace the predictor column with the sequence
  
  # Get predictions for the smooth terms and their standard errors
  pred <- predict(model, newdata = newdata, type = "terms", se.fit = TRUE)
  
  # Check if the desired smooth term is present in the prediction output
  smooth_name <- paste0("s(", predictor, ")")
  if (!(smooth_name %in% colnames(pred$fit))) {
    stop(paste("Smooth term", smooth_name, "not found in prediction output. Available terms:", 
               paste(colnames(pred$fit), collapse = ", ")))
  }
  # Extract the index of the smooth term
  term_index <- which(colnames(pred$fit) == smooth_name)
  
  # Extract the predicted smooth effect and its standard error
  fit_vals <- as.numeric(pred$fit[, term_index])
  se_vals  <- as.numeric(pred$se.fit[, term_index])
  
  # Return a tibble with the predictor sequence, predictions, and confidence intervals
  tibble(
    x = x_seq,
    fit = fit_vals,
    se = se_vals,
    lower = fit_vals - 1.96 * se_vals,
    upper = fit_vals + 1.96 * se_vals
  )
}

## Loop to generate the data frame per each scenario (time window) selected for plotting
## using "get_partial_pred"
gam_plot_data <- purrr::map_dfr(scenarios_3.5dGAM, function(scen) {
  # Get the model for the current scenario
  mod <- gam_models[[scen]]
  # Get the effective degrees of freedom (edf) for each smooth term
  edf_vals <- summary(mod)$s.table[, "edf"]
  # For each predictor, compute partial predictions and add scenario and edf info
  purrr::map2_dfr(predictors_3.5dGAM, edf_vals, function(pred, edf) {
    df_pred <- get_partial_pred_gam(mod, pred)
    df_pred %>%
      mutate(
        predictor = pred,
        scenario = scen,
        edf = edf
      )
  })
})
## Display the structure of the resulting data frame
str(gam_plot_data)

## Extract R2 and AIC per scenario; it uses "gam_models" performed above
gam_model_stats <- tibble(
  scenario = names(gam_models),
  intercept = sapply(gam_models, function(m) signif(summary(m)$p.table["(Intercept)", "Estimate"], 3)),
  r2 = sapply(gam_models, function(m) signif(summary(m)$r.sq, 3)),
  aic = sapply(gam_models, function(m) signif(AIC(m), 5))) %>%
  mutate(
    predictor = predictors_3.5dGAM[1],  # First predictor in the row
    label = paste0("Intercept = ", intercept, "\nR² = ", r2, "\nAIC = ", aic))

## Prepare the EDFs text for each panel in the plot
gam_edf_labels <- gam_plot_data %>%
  distinct(scenario, predictor, edf) %>%
  mutate(label = paste0("EDF = ", signif(edf, 3)))

## Filter gam_plot_data, gam_model_stats and gam_edf_labels, based on the scenarios GAM set before
gam_plot_data_filt <- gam_plot_data %>% filter(scenario %in% scenarios_3.5dGAM)
gam_model_stats_filt <- gam_model_stats %>% filter(scenario %in% scenarios_3.5dGAM)
gam_edf_labels_filt <- gam_edf_labels %>% filter(scenario %in% scenarios_3.5dGAM)

## Partial effect plot
p <- ggplot(gam_plot_data_filt, aes(x = x, y = fit)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.25) +
  geom_line(color = "blue", size = 1) +
  facet_grid(
    factor(scenario, levels = c("3.5_days_2_weeks", "3.5_days_3_months", "3.5_days_1_year")) ~ predictor,
    scales = "free",
    labeller = labeller(
      predictor = gam_predictor_labels)) +
  #facet_grid(scenario ~ predictor, scales = "free_x") +
  #facet_grid(factor(scenario, levels = c("3.5_days_2_weeks", 
  #                                       "3.5_days_3_months", 
  #                                       "3.5_days_1_year")) ~ predictor, scales = "free") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50")+
  
  geom_text(data = gam_edf_labels_filt,
            aes(x = Inf, y = Inf, label = label),
            hjust = 1.1, vjust = 1.5, inherit.aes = FALSE, size = 3, fontface = "bold") +
  geom_text(data = gam_model_stats_filt,
            aes(x = Inf, y = Inf, label = label),
            hjust = 1.1, vjust = 1.5, inherit.aes = FALSE, size = 3.5, fontface = "bold") +
  
  labs(
    title = "Partial Effects of GAMs on SPEAR",
    subtitle = "SPEAR ~ s(TU_ECmix) + s(Agricultural proportion) + s(Flow velocity) for each time window",
    #subtitle = "SPEAR ~ s(TU_ECmix_median) for each time window",
    caption = 
    "Each facet row corresponds to a single GAM including the chemical metric.
    All predictors were scaled (mean = 0, SD = 1); the chemical metric was log1p-transformed and scaled.
    Dashed zero line indicates baseline effect; curve shows bioindicator response change as predictor varies ±1 SD.",
    x = "Predictor value",
    y = "Partial effect on SPEAR"
  )+
  scale_y_continuous(breaks = seq(-20, 20, by = 5))+
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10))

# Export as PNG
#png(paste(user, output_graphics_file ,"GAM_partial_effect_TU_ECmix_matrix.png", sep = ""), 
#    width = 19, height = 22, units = "cm", res = 350)
p
#dev.off()
#ggsave("GAM_PartialPlot_improved.png", p, width = 25, height = 20, units = "cm", dpi = 350)
rm(p)
##

# 6. Linear Regression in the 14-days sampling period ----
## 6.1. LMM 14-days "bio-indicator ~ chem metrics + (1|YEAR)" -----

# Purpose of the function: Fits linear mixed models (LMMs) for each combination 
# of scenario and chemical metric, extracts model summaries, and filters results 
# by p-value threshold (0.1). It was created for 3.5d sampling period, so you can
# find it above.
fit_extract_lmm(
  df = metrics_14d_0_R12_abiotic_ds,
  bioind = "IBCH_2019", #Here change the bioindicator to explore
  chem_metrics = chemical_metrics_14d,
  random = "YEAR")
# VT got only one significant model

## 6.2. Forest plot: SPEAR; GI; VT; IBCH; EPT - LMM 14d ----
# Joining the results from the "fit_extract_lmm" function 
bind_rows(
  bind_rows(
    fit_extract_lmm(
      df = metrics_14d_0_R12_abiotic_ds,
      bioind = "SPEAR",
      chem_metrics = chemical_metrics_14d,
      random = "YEAR"),
    fit_extract_lmm(
      df = metrics_14d_0_R12_abiotic_ds,
      bioind = "GI_VALUE",
      chem_metrics = chemical_metrics_14d,
      random = "YEAR"),
    fit_extract_lmm(
      df = metrics_14d_0_R12_abiotic_ds,
      bioind = "VT_VALUE",
      chem_metrics = chemical_metrics_14d,
      random = "YEAR")),
  fit_extract_lmm(
    df = metrics_14d_0_R12_abiotic_ds,
    bioind = "IBCH_2019",
    chem_metrics = chemical_metrics_14d,
    random = "YEAR"),
  fit_extract_lmm(
    df = metrics_14d_0_R12_abiotic_ds,
    bioind = "EPT",
    chem_metrics = chemical_metrics_14d,
    random = "YEAR")) %>% 
  # Create a label for y-axis (scenario | term)
  mutate(
    label = paste(scenario, term, sep = " | "),
    Significance = as.character(Significance),
    singular = as.factor(singular)) %>% 
  mutate(
    AIC = as.numeric(AIC),
    BIC = as.numeric(BIC),
    obs = as.numeric(obs),
    info_label = paste0(
      #"AIC: ", round(AIC, 1),
      #"\nBIC: ", round(BIC, 1),
      "n: ", obs)) %>% 
  # DATAFRAME EDITION TO REMARK THE COMPARSION BETWEEN CHEMICAL METRICS APPROACHES
  mutate(
    chem_metric_fam = case_when(
      grepl("RQ", chem_metric) ~ "Risk Quotient",
      grepl("TU", chem_metric) ~ "Toxic Unit",
      TRUE ~ NA_character_),
    chem_metric_app = case_when(
      grepl("mix_", chem_metric) ~ "Mixture",
      grepl("max_", chem_metric) ~ "Maximum",
      TRUE ~ NA_character_),
    chem_metric_stat = case_when(
      grepl("_max", chem_metric) ~ "Maximum",
      grepl("_mean", chem_metric) ~ "Average",
      grepl("_median", chem_metric) ~ "Median",
      TRUE ~ NA_character_)) %>% 
  mutate(
    label2 = paste(scenario, chem_metric_stat, sep = " | ")) %>% 
  # Since we only have one predictor, we can remove the intercept to only plot the slope of the model
  # Also, all the intercepts established their baseline around 20 for SPEAR.
  filter(!term %in% "(Intercept)") %>% 
  
  # Plot
  ggplot(aes( x = estimate,
              y = reorder(label2, estimate),
              xmin = conf.low,
              xmax = conf.high,
              color = Significance,
              shape = as.factor(singular))) +
  geom_point(aes(size = R2cond)) +
  geom_errorbarh(height = 0.2) +
  geom_text(aes(x = Inf, label = info_label),
            hjust = 1.5, vjust = 0.5, size = 2,
            color = "black")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  facet_grid(chem_metric_app~chem_metric_fam+response, 
             scales = "free") +
  scale_color_manual(values = c("Below 0.05" = "navyblue", 
                                "Marginal significance" = "coral", 
                                na.value = "gray70")) +
  scale_shape_manual(values = c("FALSE" = 16, "TRUE" = 1)) +
  labs(
    x = "Estimate (with 95% CI)",
    y = "Scenario | Chemical Metric",
    color = "Significance",
    shape = "Singular fit",
    title = "LMM Effects Across Scenarios and Chemical Metrics",
    subtitle = "SPEAR, IBCH, VT, GI & EPT") +
  theme_bw() +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 8))

## 6.3. Forest plot to join only the most relevant bioindicators 14d (SPEAR + GI + IBCH + EPT) ----

### Dataset "chemical_metrics_diff_LM_14ds" generation
bind_rows(
  bind_rows(
    fit_extract_lmm(
      df = metrics_14d_0_R12_abiotic_ds,
      bioind = "SPEAR",
      chem_metrics = chemical_metrics_14d,
      random = "YEAR"),
    fit_extract_lmm(
      df = metrics_14d_0_R12_abiotic_ds,
      bioind = "GI_VALUE",
      chem_metrics = chemical_metrics_14d,
      random = "YEAR"),
    fit_extract_lmm(
      df = metrics_14d_0_R12_abiotic_ds,
      bioind = "VT_VALUE",
      chem_metrics = chemical_metrics_14d,
      random = "YEAR")),
  fit_extract_lmm(
    df = metrics_14d_0_R12_abiotic_ds,
    bioind = "IBCH_2019",
    chem_metrics = chemical_metrics_14d,
    random = "YEAR"),
  fit_extract_lmm(
    df = metrics_14d_0_R12_abiotic_ds,
    bioind = "EPT",
    chem_metrics = chemical_metrics_14d,
    random = "YEAR")) %>% 
  # THIS IS TO CHECK THE RANGE VALUES FOR THE BIOINDICATORS, TO INTERPRET
  #group_by(response, scenario, term) %>% 
  #summarise(max = max(estimate), min = min(estimate)) %>% View()
  
  # Create a label for y-axis (scenario | term)
  mutate(
    label = paste(scenario, term, sep = " | "),
    Significance = as.character(Significance),
    singular = as.factor(singular)) %>% 
  mutate(
    AIC = as.numeric(AIC),
    BIC = as.numeric(BIC),
    obs = as.numeric(obs),
    info_label = paste0(#"AIC: ", round(AIC, 1),
      #"\nBIC: ", round(BIC, 1),
      "n: ", obs)) %>% 
  # Dataframe edition to remark the comparison between chemical metric approaches
  mutate(
    chem_metric_fam = case_when(
      grepl("RQ", chem_metric) ~ "Risk Quotient",
      grepl("TU", chem_metric) ~ "Toxic Unit",
      TRUE ~ NA_character_),
    chem_metric_app = case_when(
      grepl("mix_", chem_metric) ~ "Mixture",
      grepl("max_", chem_metric) ~ "Maximum",
      TRUE ~ NA_character_),
    chem_metric_stat = case_when(
      grepl("_max", chem_metric) ~ "Maximum",
      grepl("_mean", chem_metric) ~ "Mean",
      grepl("_median", chem_metric) ~ "Median",
      TRUE ~ NA_character_)) %>% 
  mutate(
    label2 = paste(scenario, chem_metric_stat, sep = " | ")) %>% 
  filter(!term %in% "(Intercept)") %>% 
  filter(Significance %in% c("Below 0.05",
                             "Marginal significance")) %>% 
  group_by(response, scenario,
           # since i want to keep all the fam and approaches, let 
           # get the means grouping by all of them 
           # seeing the difference based on stats
           chem_metric_fam, chem_metric_app, chem_metric_stat) %>%
  right_join(
    expand_grid(
      response = c("SPEAR", "IBCH_2019", "GI_VALUE", "VT_VALUE", "EPT"),
      chem_metric_fam = c("Risk Quotient", "Toxic Unit"),
      chem_metric_app = c("Mixture", "Maximum"),
      chem_metric_stat = c("Median", "Maximum", "Mean"),
      scenario = unique(.$scenario)),
    by = c("response", "chem_metric_fam", 
           "chem_metric_app", 
           "chem_metric_stat", "scenario")) -> chemical_metrics_diff_LMM_14ds

### Counts to report
chemical_metrics_diff_LMM_14ds %>%
  group_by(response, chem_metric_fam) %>%
  summarise(
    min_estimate = round(min(estimate, na.rm = TRUE),3),
    min_conf.low = round(conf.low[which.min(estimate)],3),
    min_conf.high = round(conf.high[which.min(estimate)],3),
    max_estimate = round(max(estimate, na.rm = TRUE),3),
    max_conf.low = round(conf.low[which.max(estimate)],3),
    max_conf.high = round(conf.high[which.max(estimate)],3))

### Counts to report
chemical_metrics_diff_LMM_14ds %>%
  filter(response == "SPEAR") %>% # Here you can change the bioindicator to explore
  ungroup() %>% droplevels() %>% 
  mutate(range_group = as.numeric(cut_number(R2marg, 4))) %>%
  group_by(range_group) %>%
  summarise(
    min = min(R2marg),
    max = max(R2marg))

### Plot SPEAR
ggplot(
  chemical_metrics_diff_LMM_14ds %>%
    # To remove the long name of each scenario (time window in the report)
    mutate(scenario = str_remove(scenario, "^14\\_days[_]?")) %>%   
    mutate(
      scenario = factor(
        scenario,
        levels = c("1_week", "2_weeks", "1_month", "2_months", 
                   "3_months", "6_months", "1_year"),
        ordered = TRUE)) %>% 
    # Here you can change for the bioindicator you want to explore 
    filter(response %in% "SPEAR"), 
  aes(x = chem_metric_stat))+
  geom_point(aes(y = estimate, 
                 size = R2cond, color = R2marg)) +
  geom_errorbar(aes(ymin = conf.low, 
                    ymax = conf.high), 
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "coral") +
  facet_grid(
    chem_metric_fam + chem_metric_app ~ scenario , 
    scales = "free") +
  scale_color_gradient(
    name = "Marginal R²",
    low = "gray",    
    high = "#08519c", 
    limits = c(0.01, 0.08)) +
  scale_size_continuous(
    name = "Conditional R²",
    breaks = c(0.03, 0.1, 0.3, 0.5, 0.63),
    limits = c(0.03, 0.63),
    range = c(2, 7))+
  labs(
    x = "Statistical representation",
    y = "Estimate (Slope)",
    color = "Marginal R²",
    size = "Conditional R²",
    subtitle = "Linear Mixed Models: SPEAR ~ Chemical metrics (14 days)",
    #title = "Mean Estimate, Confidence Intervals, and Adjusted R² by Grouping, Response, and Category"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, 
                               hjust = 1),
    strip.background = element_rect(fill = "white"),
    strip.text.y = element_text(size = 8)) -> a

### Plot GI_VALUE
ggplot(
  chemical_metrics_diff_LMM_14ds %>% 
    # To remove the long name of each scenario (time window in the report)
    mutate(scenario = str_remove(scenario, "^14\\_days[_]?")) %>%
    mutate(
      scenario = factor(
        scenario,
        levels = c("1_week", "2_weeks", "1_month", "2_months", 
                   "3_months", "6_months", "1_year"),
        ordered = TRUE)) %>% 
    # Here you can change for the bioindicator you want to explore 
    filter(response %in% "GI_VALUE"), 
  aes(x = chem_metric_stat))+
  geom_point(aes(y = estimate, 
                 size = R2cond, color = R2marg)) +
  geom_errorbar(aes(ymin = conf.low, 
                    ymax = conf.high), 
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "coral") +
  facet_grid(
    chem_metric_fam + chem_metric_app ~ scenario , 
    scales = "free") +
  scale_color_gradient(
    name = "Marginal R²",
    low = "gray",    
    high = "#08519c", 
    limits = c(0.01, 0.08)) +
  scale_size_continuous(
    name = "Conditional R²",
    breaks = c(0.03, 0.1, 0.3, 0.5, 0.63),
    limits = c(0.03, 0.63),
    range = c(2, 7))+
  labs(
    x = "Statistical representation",
    y = "Estimate (Slope)",
    color = "Marginal R²",
    size = "Conditional R²",
    subtitle = "Linear Mixed Models: GI_VALUE  ~ Chemical metrics (14 days)",
    #title = "Mean Estimate, Confidence Intervals, and Adjusted R² by Grouping, Response, and Category"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, 
                               hjust = 1),
    strip.background = element_rect(fill = "white"),
    strip.text.y = element_text(size = 8)) -> b

### Plot IBCH
ggplot(
  chemical_metrics_diff_LMM_14ds %>% 
    # To remove the long name of each scenario (time window in the report)
    mutate(scenario = str_remove(scenario, "^14\\_days[_]?")) %>%  
    mutate(
      scenario = factor(
        scenario,
        levels = c("1_week", "2_weeks", "1_month", "2_months", 
                   "3_months", "6_months", "1_year"),
        ordered = TRUE)) %>% 
    # Here you can change for the bioindicator you want to explore 
    filter(response %in% "IBCH_2019"), 
  aes(x = chem_metric_stat))+
  geom_point(aes(y = estimate, 
                 size = R2cond, color = R2marg)) +
  geom_errorbar(aes(ymin = conf.low, 
                    ymax = conf.high), 
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "coral") +
  facet_grid(
    chem_metric_fam + chem_metric_app ~ scenario , 
    scales = "free") +
  scale_color_gradient(
    name = "Marginal R²",
    low = "gray",    
    high = "#08519c", 
    limits = c(0.01, 0.08)) +
  scale_size_continuous(
    name = "Conditional R²",
    breaks = c(0.03, 0.1, 0.3, 0.5, 0.63),
    limits = c(0.03, 0.63),
    range = c(2, 7))+
  labs(
    x = "Statistical representation",
    y = "Estimate (Slope)",
    color = "Marginal R²",
    size = "Conditional R²",
    subtitle = "Linear Mixed Models: IBCH  ~ Chemical metrics (14 days)",
    #title = "Mean Estimate, Confidence Intervals, and Adjusted R² by Grouping, Response, and Category"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, 
                               hjust = 1),
    strip.background = element_rect(fill = "white"),
    strip.text.y = element_text(size = 8)) -> c

### Plot EPT
ggplot(
  chemical_metrics_diff_LMM_14ds %>% 
    # To remove the long name of each scenario (time window in the report)
    mutate(scenario = str_remove(scenario, "^14\\_days[_]?")) %>% 
    mutate(
      scenario = factor(
        scenario,
        levels = c("1_week", "2_weeks", "1_month", "2_months", 
                   "3_months", "6_months", "1_year"),
        ordered = TRUE)) %>% 
    # Here you can change for the bioindicator you want to explore 
    filter(response %in% "EPT"), 
  aes(x = chem_metric_stat))+
  geom_point(aes(y = estimate, 
                 size = R2cond, color = R2marg)) +
  geom_errorbar(aes(ymin = conf.low, 
                    ymax = conf.high), 
                width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "coral") +
  facet_grid(
    chem_metric_fam + chem_metric_app ~ scenario , 
    scales = "free") +
  scale_color_gradient(
    name = "Marginal R²",
    low = "gray",    
    high = "#08519c", 
    limits = c(0.01, 0.08)) +
  scale_size_continuous(
    name = "Conditional R²",
    breaks = c(0.03, 0.1, 0.3, 0.5, 0.63),
    limits = c(0.03, 0.63),
    range = c(2, 7))+
  labs(
    x = "Statistical representation",
    y = "Estimate (Slope)",
    color = "Marginal R²",
    size = "Conditional R²",
    subtitle = "Linear Mixed Models: EPT  ~ Chemical metrics (14 days)",
    #title = "Mean Estimate, Confidence Intervals, and Adjusted R² by Grouping, Response, and Category"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, 
                               hjust = 1),
    strip.background = element_rect(fill = "white"),
    strip.text.y = element_text(size = 8)) -> d

### Export as PNG
#png(paste(user, output_graphics_file ,"Estimates_LM_14d.png", sep = ""), 
#    width = 38, height = 27, units = "cm", res = 350)
(a + b) / (c + d) +
  plot_annotation(
    title = "Estimate, Confidence Intervals, Conditional R² and Marginal R² by grouping, Response, and Time Windows",
    caption = 
    "Each chemical metric was lop1p-transformed and scaled.
    Facet columns represent the time windows, facet rows show the aggregation approach (Mixture or Maximum) of each chemical metric, 
    and the x-axis displays the statistical representation (Mean, Median, Maximum) of each chemical metric.")+
  plot_layout(guides = "collect") &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    plot.title = element_text(hjust = 0, face = "bold")) -> e
e
#dev.off()
rm(a,b,c,d,e)
#ggsave("LMM_14d_improved.png", e, width = 30, height = 26, units = "cm", dpi = 350)
##

## 6.4. Checking assumptions: (only homoscedasticity and normality) ----

# Purpose: 
## All the assumptions were first checked using the native plot() function for each model.
## Once normality and homoscedasticity were found to be the least met assumptions, a specific plot was
## generated for reporting purposes.
plot_lmm_assumptions_base <- function(df, response, predictor, random_effect, main_title = NULL) {
  # It would only fit if it has enough data.
  if(nrow(df) < 4) {
    plot.new(); title(main=paste(main_title, "\nNot enough data"))
    plot.new(); title(main=paste(main_title, "\nNot enough data"))
    return()
  }
  form <- as.formula(paste(response, "~", predictor, "+ (1|", random_effect, ")"))
  modelo <- lmer(form, data = df, REML = FALSE)
  
  residuals <- resid(modelo)
  fitted <- fitted(modelo)
  
  # Here, the linearity check is deactivated for plotting purposes.
  # 1. Residuals vs Fitted
  #plot(fitted, residuals,
  #     main = paste(main_title, "\nResiduals vs Fitted"),
  #     xlab = "Fitted values", ylab = "Residuals",
  #     pch = 21, col = "steelblue", cex.main = 1.2)
  #abline(h = 0, lty = 2, col = "gray40")
  
  # 2. QQ plot residuos
  qqnorm(residuals,
         main = paste(main_title, "\nNormal Q-Q"),
         pch = 16, col = "gray28", cex.main = 1.2)
  qqline(residuals, col = "red3", lwd = 2)
  
  # 3. Scale-Location plot
  sqrt_abs_resid <- sqrt(abs(residuals))
  plot(fitted, sqrt_abs_resid,
       main = paste(main_title, "\nScale-Location"),
       xlab = "Fitted values", ylab = expression(sqrt("|Residuals|")),
       pch = 16, col = "green4", cex.main = 1.2)
  abline(h = mean(sqrt_abs_resid), lty = 2, col = "gray40")
}

## Preparing the loop 
### Variable definition
indices <- c("GI_VALUE", "EPT", "SPEAR")
scenarios <- c("14_days_1_week", "14_days_2_months", "14_days_1_year")
chem_metrics <- c(#TU_NOECmix_median", 
                  "CRQmix_median")
random_effect <- "YEAR"  

# Export as PNG
#png(paste(user, output_graphics_file ,"Model_diagnostics_LMM_GI_EPT_SPEAR.png", 
#          sep = ""), width = 50, height = 30, units = "cm", res = 350)

## Set the panel
par(mfrow = c(3,6))
## Loop based on the pre-defined variables (lists)
for (indice in indices) {
  for (scenario in scenarios) {
    for (chem_metric in chem_metrics) {
      df_sub <- metrics_14d_0_R12_abiotic_ds %>%
        filter(scenario == !!scenario)
      main_title <- paste(indice, "|", scenario, "|", chem_metric, sep ="")
      plot_lmm_assumptions_base(df_sub, response = indice, predictor = chem_metric, 
                                random_effect = random_effect, main_title = main_title)
    }
  }
}

par(mfrow = c(1, 1))
#dev.off()
rm(df_sub, main_title)
##


## 6.5. LMMs and GAMMs (overlapped) for exploration ---- 
##____________________________________________________________________________________________________
## Polished exportation of the plots for all the scenarios (time windows) or only the 
## representative time windows.

## This plots, ca be replicated for every bioindicator in case the further exploration is needed
## Now it's set for SPEAR, exploration showed as not meaningful results for the rest of bioindicators
##____________________________________________________________________________________________________

### Filter scenarios and select only the representative time windows
scenarios_14d_LMM <- c("14_days_2_weeks", "14_days_3_months", "14_days_1_year")

### Long format of the dataset
df_long_lmm_gamms <- metrics_14d_0_R12_abiotic_ds %>%
  filter(scenario %in% c(scenarios_14d_LMM)) %>% #scenarios_14d
  select(scenario, YEAR, SPEAR, 
         #TU_NOECmix_median, CRQmix_median, 
         TU_NOECmix_max, CRQmix_max
         ) %>%
  pivot_longer(
    cols = c(#TU_NOECmix_median, CRQmix_median, 
      TU_NOECmix_max,  CRQmix_max
      ),
    names_to = "chem_metric",
    values_to = "chem_value")

### Compute model statistics for each panel
stats_lmm_gams <- df_long_lmm_gamms %>%
  group_by(scenario, chem_metric) %>%
  # For each group, apply the following function (group_modify expects a function returning a data frame)
  group_modify(~{
    
    # The current group's data frame
    df <- .
    
    # LMM
    lmm <- tryCatch(lmer(SPEAR ~ chem_value + (1|YEAR), data = df), error = function(e) NULL)
    r2_lmm <- tryCatch(r.squaredGLMM(lmm), error = function(e) NA)
    r2_marginal <- tryCatch(signif(r2_lmm[1], 2), error = function(e) NA)
    r2_conditional <- tryCatch(signif(r2_lmm[2], 2), error = function(e) NA)
    p_LMM <- tryCatch(signif(summary(lmm)$coefficients["chem_value", "Pr(>|t|)"], 2), error = function(e) NA)
    AIC_LMM <- tryCatch(signif(AIC(lmm), 4), error = function(e) NA)
    
    # GAMM
    gamm_mod <- tryCatch(gamm(SPEAR ~ s(chem_value, k = 6), random = list(YEAR=~1), data = df), error = function(e) NULL)
    edf_gamm <- tryCatch(signif(summary(gamm_mod$gam)$s.table[1, "edf"], 2), error = function(e) NA)
    p_gamm <- tryCatch(signif(summary(gamm_mod$gam)$s.table[1, "p-value"], 2), error = function(e) NA)
    adjr2_gamm <- tryCatch(signif(summary(gamm_mod$gam)$r.sq, 2), error = function(e) NA)
    AIC_GAMM <- tryCatch(signif(AIC(gamm_mod$lme), 4), error = function(e) NA)
    
    # Return a tibble with all extracted statistics for the current group
    tibble(
      r2_marginal = r2_marginal,
      r2_conditional = r2_conditional,
      p_LMM = p_LMM,
      AIC_LMM = AIC_LMM,
      edf_gamm = edf_gamm,
      p_gamm = p_gamm,
      adjr2_gamm = adjr2_gamm,
      AIC_GAMM = AIC_GAMM)
    }) %>%
  # Remove grouping structure
  ungroup() %>%
  # Create a formatted text summary for each row (scenario and chemical metric)
  mutate(
    stat_text = paste0(
      "LMM:\nMargR² = ", r2_marginal,
      "\nCondR² = ", r2_conditional,
      "\np = ", p_LMM,
      "\nAIC = ", AIC_LMM,
      "\n\nGAMM:\nEDF = ", edf_gamm,
      "\np = ", p_gamm,
      "\nAdjR² = ", adjr2_gamm,
      "\nAIC = ", AIC_GAMM))

### Get predictions and se for LMMs and GAMMs
get_preds_lmm_gamm <- function(df) {
  
  # LMM
  lmm <- tryCatch(lmer(SPEAR ~ chem_value + (1|YEAR), data = df), error = function(e) NULL)
  new_x <- seq(min(df$chem_value, na.rm=TRUE), max(df$chem_value, na.rm=TRUE), length.out=100)
  lmm_pred <- tryCatch({
    pred <- predict(lmm, newdata = data.frame(chem_value = new_x, YEAR = df$YEAR[1]), 
                    re.form = NA, se.fit = TRUE)
    tibble(chem_value = new_x, LMM = pred$fit, LMM_se = pred$se.fit)
  }, error = function(e) tibble(chem_value = new_x, LMM = NA, LMM_se = NA))
  
  # GAMM
  gamm_mod <- tryCatch(gamm(SPEAR ~ s(chem_value, k = 6), random = list(YEAR=~1), data = df), error = function(e) NULL)
  gamm_pred <- tryCatch({
    pred <- predict(gamm_mod$gam, newdata = data.frame(chem_value = new_x), se.fit = TRUE)
    tibble(chem_value = new_x, GAMM = pred$fit, GAMM_se = pred$se.fit)
  }, error = function(e) tibble(chem_value = new_x, GAMM = NA, GAMM_se = NA))
  left_join(lmm_pred, gamm_pred, by = "chem_value")
}

### Apply the function into the df_long_lmm_gamms data set
preds_lmm_gamms <- df_long_lmm_gamms %>%
  group_by(scenario, chem_metric) %>%
  group_modify(~get_preds_lmm_gamm(.)) %>%
  ungroup()

### Plot 
p <- ggplot(df_long_lmm_gamms, aes(x = chem_value, y = SPEAR)) +
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray35") +
  geom_ribbon(data = preds_lmm_gamms, aes(x = chem_value, ymin = LMM - LMM_se, ymax = LMM + LMM_se),
              fill = "pink", alpha = 0.2, inherit.aes = FALSE, na.rm = TRUE) +
  geom_line(data = preds_lmm_gamms, aes(x = chem_value, y = LMM), color = "red", 
            linetype = "dashed", linewidth = 1, na.rm = TRUE) +
  geom_ribbon(data = preds_lmm_gamms, aes(x = chem_value, ymin = GAMM - GAMM_se, ymax = GAMM + GAMM_se), 
              fill = "lightblue", alpha = 0.3, inherit.aes = FALSE, na.rm = TRUE) +
  geom_line(data = preds_lmm_gamms, aes(x = chem_value, y = GAMM), color = "blue", 
            linewidth = 1, na.rm = TRUE) +
  facet_grid(chem_metric ~ scenario, scales = "free") +
  geom_text(
    data = stats_lmm_gams,
    aes(x = Inf, y = Inf, label = stat_text),
    hjust = 1, vjust = 1.1, size = 3, fontface = "bold", 
    inherit.aes = FALSE) +
  labs(
    title = "SPEAR vs Chemical metrics",
    subtitle = paste(
      "Time windows: 2-weeks, 3-months, 1-year\n",
      "LMM (red dashed) and GAMM (blue) with model statistics\n",
      "LMM: Marginal R² (fixed effects), Conditional R² (fixed + random), p-value and AIC.\n",
      "GAMM: EDF and p-value for smooth term, AdjR² for GAM part, and AIC of mixed model."
    ),
    x = "scaled(log1p(Chemical metric))",
    y = "SPEAR",
    caption = 
    "The horizontal line at 0 indicates the lower limit of SPEAR values. Values below 0 are not plausible as SPEAR cannot be negative.
    Model performance is summarized by R2, AIC (Akaike Information Criterion; lower values indicate better fit), EDF (effective degrees of freedom), and p-value (model significance)."
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "white"),  
    strip.text = element_text(face = "bold"))
p

# Export as PNG
#ggsave(
#  filename = file.path(user, output_graphics_file, "SPEAR_14d_LMM_GAMM_R2_pval_AIC_2w_3m_1y.png"),
#  plot = p,
#  width = 28, height = 20,
#  units = "cm", dpi = 350)
rm(p)

##

## 6.6. Intermediate conclusion for LMM ----
# The implementation of GAMMs did not improve the explanatory power or the interpretation of the
# relationship between SPEAR and chemical metrics, either helped with the improvement of the models
# including the other bioindicators, so that LMM is the model selected for further exploration including
# environmental gradients. 


# 7. LMM 14d: SPEAR ~ chem metrics + environmental gradients + (1|YEAR) ----
## Data-set preparation only with the relevant abiotic factors 
metrics_14d_0_R12_abiotic_ds %>% 
  select(
    # Organization columns
    Gewasser, scenario, YEAR,
    # Bioindicators
    SPEAR, GI_VALUE, IBCH_2019, EPT,
    # Chemical metrics
    contains("max"), contains("mix"),
    # Temporal
    PP_mean_mean,
    # Land use
    Agricultural_area_frac,
    Forest_area_frac,
    Urban_area_frac,
    # Hydrological
    Flow_velocity,
    Avg_mod_discharge_m3_s,
    # Morphological
    Ecomorphology_0_12,
    # Temporal
    Mod_max_temp_summer ,
    # Spatial
    Z) %>% 
  na.omit() -> df_metrics_abiotic_14d

#_______________________________________________________________________________
# CONCEPTUALIZATION BASED ON RELEVANT FACTORS FOR SPEAR AND OBSERVED 
# ENVIRONMENTAL GRADIENTS IN THE PCA

# 1. Chemical + Agri + Flow + Ecomorphology
# 1.1. Chemical + Agri + Discharge + Ecomorphology
# 2. Chemical + Agri + Flow + Urban
# 2.1. Chemical + Agri + Discharge + Urban
# 3. Chemical + Agri + Flow
# 3.1. Chemical + Agri + Discharge
# 4. Chemical + Agri
# 5. Chemical + Discharge
# 6. Chemical + Flow

# Final thought: Keep chemical and agriculture, prioritize flow and urban
#_______________________________________________________________________________


## Adaptation of "fit_extract_lmm" to include more predictors ----
### Purpose: Adapt the functions that only includes chemical metrics
### to ble able to include different model structures. So that, we can
### compare the performance of each one.
fit_extract_lmm_abiotic_multi_chem <- function(df, bioind, chem_metrics, 
                                               abiotic_structures, 
                                               random, 
                                               pval_threshold = 0.1,
                                               extra_predictors = NULL) {
  
  # ______________________________________________________________________________________________
  # Inputs:
  # df: Data frame containing the data for analysis
  # bioind: String, name of the biological indicator variable
  # chem_metrics: Vector of strings, names of chemical metrics to include
  # abiotic_structures: Named list, each element is a vector of abiotic predictor names for a 
  #                     model structure
  # random: String, name of the random effect variable (e.g., "YEAR")
  # pval_threshold: Numeric, threshold for filtering significant results (default = 0.1)
  # extra_predictors: Named list, each element is a vector of extra predictors for a model 
  #                   structure (optional)
  #
  # Outputs:
  # A tibble (data frame) with one row per scenario, chemical metric, and abiotic model structure,
  # containing model statistics, fixed effect summaries, and VIF information.
  # Only models with an overall p-value <= pval_threshold are included.
  # Each row also includes a "Significance" annotation for each fixed effect.
  # ______________________________________________________________________________________________
  
  # Create all combinations of scenario, chemical metric, and abiotic model structure
  combinations <- tidyr::crossing(
    scenario = unique(df$scenario),
    chem_metric = chem_metrics,
    abiotic_model = names(abiotic_structures))
  
  # Internal function to fit and extract model info for each combination
  fit_extract_lmm_inner <- function(scen, chem_metric, abiotic_model, random) {
    # Subset data for the current scenario
    sub_df <- df %>% filter(scenario == scen)
    # Get base terms for the current abiotic model
    base_terms <- unlist(abiotic_structures[[abiotic_model]])
    # Add extra predictors if specified for this model
    model_extras <- extra_predictors[[abiotic_model]]
    abiotic_terms <- unique(c(base_terms, model_extras))
    # Construct the right-hand side of the formula
    rhs <- paste(c(chem_metric, abiotic_terms), collapse = " + ")
    # Construct the full formula for the LMM
    formulation <- as.formula(paste(bioind, "~", rhs, "+ (1|", random, ")"))
    
    # Fit the LMM with error handling
    mod <- tryCatch(
      lmer(formulation, data = sub_df),
      error = function(e) return(NULL)
    )
    if (is.null(mod)) return(NULL)
    
    # --- Calculate overall model p-value using likelihood ratio test (LRT) ---
    # Fit a null model (intercept and random effect only)
    null_formula <- as.formula(paste(bioind, "~ 1 + (1|", random, ")"))
    null_mod <- tryCatch(
      lmer(null_formula, data = sub_df, REML = FALSE),
      error = function(e) return(NULL)
    )
    if (is.null(null_mod)) return(NULL)
    # Calculate LRT and extract p-value
    model_pval <- tryCatch({
      lrt <- anova(null_mod, mod, refit = FALSE)
      lrt$`Pr(>Chisq)`[2]  # The p-value is in the second row
    }, error = function(e) NA_real_)
    
    # --- Calculate VIF for fixed effects (for reference, not valid for mixed models) ---
    # Fit a fixed-effects model to calculate VIF
    vif_vals <- tryCatch(
      car::vif(lm(as.formula(paste(bioind, "~", rhs)), data = sub_df)),
      error = function(e) NA
    )
    if (all(is.na(vif_vals))) {
      vif_summary <- NA_character_
      vif_min <- NA_real_
      vif_min_var <- NA_character_
      vif_max <- NA_real_
      vif_max_var <- NA_character_
    } else {
      # Ensure names are present
      if (is.null(names(vif_vals)) || all(names(vif_vals) == "")) {
        names(vif_vals) <- paste0("V", seq_along(vif_vals))
      }
      # Summarize VIF values for output
      vif_summary <- paste(paste(names(vif_vals), round(vif_vals, 2), sep = "="), collapse = "; ")
      vif_min <- min(vif_vals, na.rm = TRUE)
      vif_max <- max(vif_vals, na.rm = TRUE)
      vif_min_var <- names(vif_vals)[which.min(vif_vals)]
      vif_max_var <- names(vif_vals)[which.max(vif_vals)]
    }
    
    # --- Get R2 (marginal and conditional) for the LMM ---
    r2 <- tryCatch(
      performance::r2_nakagawa(mod),
      error = function(e) list(R2_marginal = NA_real_, R2_conditional = NA_real_)
    )
    
    # --- Return results as a tibble row ---
    tibble(
      response = bioind,
      chem_metric = chem_metric,
      abiotic_model = abiotic_model,
      predictors = rhs,
      model = "LMM",
      random = random,
      scenario = scen,
      converged = is.null(mod@optinfo$conv$lme4$messages),
      singular = isSingular(mod),
      tidy = list(broom.mixed::tidy(mod, effects = "fixed", conf.int = TRUE)), # term-level p-values
      AIC = AIC(mod),
      BIC = BIC(mod),
      obs = nobs(mod),
      R2marg = r2$R2_marginal,
      R2cond = r2$R2_conditional,
      model_pval = model_pval,
      vif_summary = vif_summary,
      vif_min = vif_min,
      vif_min_var = vif_min_var,
      vif_max = vif_max,
      vif_max_var = vif_max_var
    )
  }
  
  # Apply the inner function to all combinations using pmap_dfr
  results <- purrr::pmap_dfr(
    list(
      scen = combinations$scenario,
      chem_metric = combinations$chem_metric,
      abiotic_model = combinations$abiotic_model
    ),
    ~ fit_extract_lmm_inner(scen = ..1, chem_metric = ..2, abiotic_model = ..3, random = random)
  )
  
  # If no models were fitted, return an empty tibble and print a message
  if (nrow(results) == 0) {
    message("No models were successfully fitted.")
    return(tibble())
  }
  
  # Unnest term-level results, filter by overall model p-value, and annotate term significance
  results %>%
    tidyr::unnest(tidy) %>%
    group_by(scenario, chem_metric, abiotic_model) %>%
    filter(model_pval <= pval_threshold) %>%
    mutate(Significance =
             dplyr::case_when(
               p.value <= 0.05 ~ "Below 0.05",
               p.value > 0.05 & p.value <= 0.1 ~ "Marginal significance",
               p.value > 0.1 ~ "No significance",
               TRUE ~ NA_character_
             )) -> results2
  
  # Print the final results for inspection
  print(results2)
}

## Application of the function "fit_extract_lmm_abiotic_multi_chem

### Model structures 
mod1 <- c("Agricultural_area_frac", "Flow_velocity", "Ecomorphology_0_12")
mod2 <- c("Agricultural_area_frac", "Flow_velocity", "Urban_area_frac")
mod3 <- c("Agricultural_area_frac", "Flow_velocity")
mod1_1 <- c("Agricultural_area_frac", "Avg_mod_discharge_m3_s", "Ecomorphology_0_12")
mod2_1 <- c("Agricultural_area_frac", "Avg_mod_discharge_m3_s", "Urban_area_frac")
mod3_1 <- c("Agricultural_area_frac", "Avg_mod_discharge_m3_s")
mod4 <- c("Agricultural_area_frac")
mod5 <- c("Avg_mod_discharge_m3_s")
mod6 <- c("Flow_velocity")

### Grouping the model structures
abiotic_structures <- list(mod1 = mod1, mod2 = mod2, mod3 = mod3,
                           mod1_1 = mod1_1, mod2_1 = mod2_1, mod3_1 = mod3_1,
                           mod4 = mod4, mod5 = mod5, mod6 = mod6)

### Attempt
fit_extract_lmm_abiotic_multi_chem(
  df = df_metrics_abiotic_14d, 
  bioind = "SPEAR",
  chem_metrics =  c("CRQmix_max", "Urban_area_frac"), # the best chemical metrics
  abiotic_structures = abiotic_structures, #abiotic_structures
  random = "YEAR",
  pval_threshold = 0.1,
  extra_predictors = NULL) # This was created for the separated implementation of precipitation.
                           # However, no significant model was obtained including precipitation    

## 7.1. Forest plot (Summary): SPEAR ~ Chemical metric + Env gradients (9 models) ----
### This plot also works as a summary for all the models attempted since it shows
### AIC, n, only the ones VIF<5, Conditional R2 and their slope, so no table
### was created for them.
# Joining the results from the "fit_extract_lmm_abiotic_multichem" function 
fit_extract_lmm_abiotic_multi_chem(
  df = df_metrics_abiotic_14d, 
  bioind = "SPEAR",
  chem_metrics =  c("CRQmix_max"),  
  abiotic_structures = abiotic_structures,
  random = "YEAR",
  pval_threshold = 0.1,
  extra_predictors = NULL) %>% 
  mutate(
    label = paste(scenario, term, sep = " | "),
    Significance = as.character(Significance)) %>% 
  mutate(
    AIC = as.numeric(AIC),
    #BIC = as.numeric(BIC),
    obs = as.numeric(obs),
    info_label = paste0("AIC: ", round(AIC, 1), #"\nBIC: ", round(BIC, 1),
                        "\nn: ", obs)) %>% 
  # Filter of intercept and variance inflation factors bellow 5
  filter(!term %in% "(Intercept)", vif_max <=5) %>% 
  mutate(is_chem = ifelse(term == chem_metric, TRUE, FALSE)) %>%
  ggplot(aes( x = estimate,
              y = reorder(term, estimate),
              xmin = conf.low,
              xmax = conf.high,
              color = Significance)) +
  geom_point(aes(size = R2cond, 
                 shape = is_chem)) +
  # To differentiate the chemical metric from the rest of estimates (slopes)
  scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 16)) + # 17 = triangle, 16 = circle
  geom_errorbarh(height = 0.2) +
  geom_text(aes(x = Inf, label = info_label),
            hjust = 1, vjust = 0.5, size = 2,
            color = "black")+
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  scale_color_manual(values = c("Below 0.05" = "navyblue", 
                                "Marginal significance" = "coral", 
                                "No significance" = "gray70")) +
  facet_grid(scenario~abiotic_model, scales = "free")+
  labs(
    x = "Estimate (with 95% CI)",
    y = "Predictors",
    color = "Significance",
    size = "R2cond",
    shape = "Chemical metric",
    title = "LMM Effects Across Time Windows",
    caption = 
    "All predictors were scaled (mean = 0, SD = 1); the chemical metric was log1p-transformed and scaled.
    Facet columns represent each proposed model, facet rows show each time window, and chemical metrics are represented by a triangle.
    Akaike information criterion (AIC) and number of observations (n) are displayed for each subset of data and predictor, same values belong to either the same model or the same sub-set of data.",
    subtitle = paste0(
      "Model structures: All of them include CRQmix; max\n",
      "Mod1: SPEAR ~ Agricultural_area_frac + Flow_velocity + Ecomorphology_0_12 + (1|YEAR)\n",
      "Mod1_1: SPEAR ~ Agricultural_area_frac + Avg_mod_discharge_m3_s + Ecomorphology_0_12 + (1|YEAR)\n",
      "Mod2: SPEAR ~ Agricultural_area_frac + Flow_velocity + Urban_area_frac + (1|YEAR)\n",
      "Mod2_1: SPEAR ~ Agricultural_area_frac + Avg_mod_discharge_m3_s + Urban_area_frac + (1|YEAR)\n",
      "Mod3: SPEAR ~ Agricultural_area_frac + Flow_velocity + (1|YEAR)\n",
      "Mod3_1: SPEAR ~ Agricultural_area_frac + Avg_mod_discharge_m3_s + (1|YEAR)\n",
      "Mod4: SPEAR ~ Agricultural_area_frac + (1|YEAR)\n",
      "Mod5: SPEAR ~ Avg_mod_discharge_m3_s + (1|YEAR)\n",
      "Mod6: SPEAR ~ Flow_velocity + (1|YEAR)")) +
  theme_bw() +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 8))+
  theme(
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10)) -> p

# Export as PNG
#png(paste(user, output_graphics_file ,"Forest_plot_LMM_14d.png", sep = ""), 
#    width = 50, height = 30, units = "cm", res = 350)
p
#dev.off()
rm(p)


## 7.2. Fast check of the assumptions when LMM  ----
lmm_diagnostic_plots <- list()
for (scen in scenarios_14d) {
  df <- df_metrics_abiotic_14d %>% filter(scenario == scen)
  
  model <- lmer(SPEAR ~ 
                  CRQmix_median + 
                  Agricultural_area_frac +
                  Flow_velocity +
                  Urban_area_frac +
                  (1|YEAR),
                data = df)
  
  # Residues and predicted
  res <- resid(model)
  fits <- fitted(model)
  
  # 1. Residues vs predicted
  p1 <- ggplot(data.frame(fits, res), aes(x = fits, y = res)) +
    geom_point() +
    geom_hline(yintercept = 0, 
               linetype = "dashed", 
               color = "red") +
    labs(title = paste("Residues vs Predictes -", scen), 
         x = "Predicted values", y = "Residues")
  
  # 2. QQ plot de residuos
  p2 <- ggplot(data.frame(res), aes(sample = res)) +
    stat_qq() +
    stat_qq_line() +
    labs(title = paste("QQ plot -", scen), 
         x = "Theoretical", 
         y = "Observed")
  
  # 3. Histograma de residuos
  p3 <- ggplot(data.frame(res), aes(x = res)) +
    geom_histogram(bins = 30, fill = "skyblue", color = "black") +
    labs(title = paste("Histogram -", scen), x = "Residues")
  
  # 4. Residuos vs cada predictor (ejemplo con uno)
  p4 <- ggplot(data.frame(df, res), aes(x = CRQmix_median, y = res)) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(title = paste("Residues vs CRQmix_median -", scen), 
         x = "CRQmix_median", y = "Residues")
  
  # Guardar los plots en una lista
  lmm_diagnostic_plots[[scen]] <- list(
    residuos_vs_predichos = p1,
    qqplot_residuos = p2,
    hist_residuos = p3,
    residuos_vs_CRQmix = p4)
}
rm(df, model, res, fits, p1, p2, p3, p4)

# Grid arrange, so that an specific situation can be observed
grid.arrange(
  lmm_diagnostic_plots[["14_days_2_months"]]$residuos_vs_predichos,
  lmm_diagnostic_plots[["14_days_2_months"]]$qqplot_residuos,
  lmm_diagnostic_plots[["14_days_2_months"]]$hist_residuos,
  lmm_diagnostic_plots[["14_days_2_months"]]$residuos_vs_CRQmix, ncol = 2)


# CROSS VALIDATION



## 7.3. Partial effect plots ----

#____________________________________________________________________________________
## Here you have to change the name of the bioindicator and chemical metric, so it 
## is not necessary to create more objects and you won't end up with too many models.

## Just be careful, because the object "models_lmm" is then used for the partial effect 
## plots. To produce a new plot, you need to run a new model from here. So that
## you can check the main values obtained from the new structure.
#____________________________________________________________________________________

scenarios_14d_LMM <- c("14_days_2_weeks", "14_days_3_months", "14_days_1_year")
predictors_14d <- c("Agricultural_area_frac", "Urban_area_frac", "Flow_velocity", "CRQmix_max")

## Let's create the empty lists
lmm_models <- list()
lmm_summaries <- list()
lmm_vif_results <- list()
lmm_aic_vals <- c()
lmm_r2m_vals <- c()
lmm_r2c_vals <- c()
lmm_pvals <- list()
lmm_fstats <- list()

## Loop to get fitted LMMs per scenario (time window)
for (scen in scenarios_14d_LMM) {
  # _______________________________________________________________________________
  # Inputs:
  # scenarios_14d_LMM: Vector of scenario names to loop over.
  # df_metrics_abiotic_14d: Data frame with columns including scenario, SPEAR,
  #                         CRQmix_max, Agricultural_area_frac, Urban_area_frac,
  #                         Flow_velocity, YEAR.
  # predictors_14d: Vector of predictor names (used for partial effects plotting).
  # lmm_models, lmm_summaries, lmm_vif_results, lmm_aic_vals, lmm_r2m_vals,
  # lmm_r2c_vals, lmm_pvals, lmm_fstats: Lists/vectors to store results.
  #
  # Outputs:
  # lmm_models: List of fitted LMM models, named by scenario.
  # lmm_summaries: List of model summaries, named by scenario.
  # lmm_vif_results: List of VIF results for fixed effects, named by scenario.
  # lmm_aic_vals: Named vector of AIC values, by scenario.
  # lmm_r2m_vals: Named vector of marginal R² values, by scenario.
  # lmm_r2c_vals: Named vector of conditional R² values, by scenario.
  # lmm_pvals: List of p-values for fixed effects, by scenario.
  # lmm_fstats: List of F-statistics for fixed effects, by scenario.
  # plot_data: Data frame of partial effects (centered at zero, no intercept),
  #            for each scenario and predictor.
  # _______________________________________________________________________________
  
  # Subset data for the current scenario
  df <- df_metrics_abiotic_14d %>% filter(scenario == scen)
  
  # Fit LMM: SPEAR as response, fixed effects for CRQmix_max, Agricultural_area_frac,
  # Urban_area_frac, Flow_velocity, and random intercept for YEAR
  model <- lmer(SPEAR ~ 
                  CRQmix_max + 
                  Agricultural_area_frac + 
                  Urban_area_frac + 
                  Flow_velocity + 
                  (1|YEAR), 
                data = df)
  
  # Store the fitted model in the list, using scenario name as the key
  lmm_models[[scen]] <- model
  # Store the model summary in the list
  lmm_summaries[[scen]] <- summary(model)
  # Calculate VIF for fixed effects (using linear model for compatibility)
  lmm_vif_results[[scen]] <- vif(lm(SPEAR ~ 
                                      CRQmix_max +
                                      Agricultural_area_frac + 
                                      Urban_area_frac + 
                                      Flow_velocity, 
                                    data = df))
  # Store AIC for the model
  lmm_aic_vals[scen] <- AIC(model)
  
  # Calculate marginal and conditional R²
  lmm_r2s <- performance::r2_nakagawa(model)
  lmm_r2m_vals[scen] <- signif(lmm_r2s$R2_marginal, 3)
  lmm_r2c_vals[scen] <- signif(lmm_r2s$R2_conditional, 3)
  
  # Extract p-values and F-statistics for fixed effects using broom.mixed
  tidy_mod <- broom.mixed::tidy(model, effects = "fixed")
  lmm_pvals[[scen]] <- tidy_mod$p.value
  lmm_fstats[[scen]] <- tidy_mod$statistic
}
# To remove temporary objects that are not useful later
rm(df, model, tidy_mod)

# Partial effects centered at zero (no intercept)
lmm_plot_data <- purrr::map_dfr(scenarios_14d_LMM, function(scen) {
  # Get the model for the current scenario
  mod <- lmm_models[[scen]]
  # Extract the intercept from the fixed effects
  intercept <- fixef(mod)[1]
  # For each predictor, compute partial effects and center at zero
  purrr::map_dfr(predictors_14d, function(pred) {
    # Get predicted effects for the current predictor
    eff <- ggpredict(mod, terms = pred)
    # Return a tibble with centered effects and scenario/predictor info
    tibble(
      x = eff$x,
      fit = eff$predicted - intercept,
      se = eff$std.error,
      lower = eff$conf.low - intercept,
      upper = eff$conf.high - intercept,
      predictor = pred,
      scenario = scen
    )
  })
})
# Force the order of the predictors as factors
lmm_plot_data$predictor <- factor(lmm_plot_data$predictor, levels = predictors_14d)

# Get the estimates (slopes) for each predictor and model
lmm_slopes <- purrr::map_dfr(names(lmm_models), function(scen) {
  mod <- lmm_models[[scen]]
  ests <- broom.mixed::tidy(mod, effects = "fixed") %>%
    filter(term != "(Intercept)") %>%
    select(term, estimate) %>%
    mutate(scenario = scen)
  colnames(ests) <- c("predictor", "slope", "scenario")
  ests
})
# Force predictor order as a factor for consistent plotting
lmm_slopes$predictor <- factor(lmm_slopes$predictor, levels = predictors_14d)

# Create a data frame with model statistics (intercept, R², AIC) for each scenario
lmm_model_stats <- tibble(
  scenario = names(lmm_models),
  intercept = sapply(lmm_models, function(m) signif(fixef(m)[1], 3)),
  r2_marg = lmm_r2m_vals,
  r2_cond = lmm_r2c_vals,
  aic = signif(lmm_aic_vals, 5))

# Join slope data with model statistics for labeling
lmm_label_data <- lmm_slopes %>%
  left_join(lmm_model_stats, by = "scenario") %>%
  mutate(
    # For the panel with "Agricultural_area_frac" (or "CRQmix_max" if changed), show full stats;
    # for other panels, show only the slope
    label = ifelse(
      predictor == "CRQmix_max",  # EDIT: Change to "CRQmix_max" if needed
      paste0(
        "Slope = ", signif(slope, 3),
        "\nIntercept = ", intercept,
        "\nR²_marg = ", r2_marg,
        "\nR²_cond = ", r2_cond,
        "\nAIC = ", aic),
      paste0("Slope = ", signif(slope, 3))))

lmm_predictor_labels <- c(
  "Agricultural_area_frac" = "Agricultural proportion",
  "Urban_area_frac" = "Urban proportion",
  "Flow_velocity" = "Flow velocity",
  "CRQmix_max" = "CRQ (mixture, max)")

# Final plot
p <- ggplot(lmm_plot_data, aes(x = x, y = fit)) +
  # Add confidence interval ribbon
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.25) +
  # Add main effect line
  geom_line(color = "blue", linewidth = 1) +
  # Facet by scenario (rows) and predictor (columns), with free x-axis scales
  #facet_grid(factor(scenario, 
  #                  levels = scenarios_14d_LMM) ~ predictor, scales = "free_x") +
  facet_grid(factor(scenario, levels = scenarios_14d_LMM) ~ predictor, 
             scales = "free_x", 
             labeller = labeller(predictor = lmm_predictor_labels)) +
  # Add reference line at zero
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  # Add model statistics labels to each panel
  geom_text(
    data = lmm_label_data,
    aes(x = Inf, y = Inf, label = label),
    hjust = 1.1, vjust = 1.1, inherit.aes = FALSE, size = 3.3, fontface = "bold"
  ) +
  # Add titles and captions
  labs(
    title = "Partial Effects of LMM on SPEAR",
    subtitle = 
      "Effects centered by subtracting the model intercept.\nSPEAR ~ Agricultural proportion + Urban proportion + Flow velocity + (1|YEAR) for each time window",
    caption = 
      "Each facet row corresponds to a single LMM including each predictor.
      Y-axis shows net effect of each predictor (intercept removed) for comparability with GAMs.
      All predictors were scaled (mean = 0, SD = 1); the chemical metric was log1p-transformed and scaled.
      Dashed zero line indicates baseline effect; curve shows bioindicator response change as predictor varies ±1 SD.",
    x = "Predictor value (scaled)",
    y = "Partial effect on SPEAR (centered)"
  ) +
  # Set y-axis breaks for clarity
  scale_y_continuous(breaks = seq(-20, 20, by = 5)) +
  # Apply theme and styling
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey95"),
    strip.text = element_text(face = "bold"),
    panel.grid.major = element_line(color = "grey90"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10)
  )

# Export as PNG
#png(paste(user, output_graphics_file ,"LMM_partial_effect_NoChem_SPEAR_matrix.png", sep = ""), 
#    width = 25, height = 20, units = "cm", res = 350)
p
#dev.off()
#ggsave("LMM_PartialPlot_improved.png", p, width = 25, height = 20, units = "cm", dpi = 350)
rm(p)



# 8. Cross-Validation (RSME and R2) + Random Forest ----
## 8.1. 3.5d best GAM model in CV + RF ----

### 8.1.1. Data set  ----
ds_3.5d_GAM_CV_RF <- df_metrics_abiotic_3.5d

### Define Bioindicators and Chemical Metrics to include in the model
bioindicators_3.5d_CV <- c("SPEAR")
chemical_metrics_3.5_CV <- c("TU_ECmix_median")

### 8.1.2. Function: Cross-validation of three models (LM, GAM, RF) ----
# Purpose: This function performs k-fold cross-validation for three types of models:
# - Linear Model (LM)
# - Generalized Additive Model (GAM)
# - Random Forest (RF)
# It returns average performance metrics (RMSE and R-squared) and their standard errors.

cv_model_3.5d <- function(data_subset, bioindicator, predictor, k_folds = 5) {
  set.seed(123)  # Set seed for reproducibility of random splits
  
  # Create folds for cross-validation using the bioindicator as grouping variable
  folds <- caret::createFolds(data_subset[[bioindicator]], k = k_folds, list = TRUE)
  
  # Initialize vectors to store performance metrics for each fold
  rmse_lm <- numeric(k_folds)
  r2_lm <- numeric(k_folds)
  rmse_gam <- numeric(k_folds)
  r2_gam <- numeric(k_folds)
  rmse_rf <- numeric(k_folds)
  r2_rf <- numeric(k_folds)
  
  # Loop through each fold
  for (i in 1:k_folds) {
    test_idx <- folds[[i]]                  # Indices for test data in this fold
    train_data <- data_subset[-test_idx, ]  # Training data for this fold
    test_data <- data_subset[test_idx, ]    # Test data for this fold
    
    # --- LINEAR MODEL (LM) ---
    # Define formula for LM using the bioindicator as response and selected predictors (i.e. BEST GAM structure)
    lm_formula <- as.formula(paste(bioindicator, "~", 
                                   predictor, "+ Agricultural_area_frac + Flow_velocity"))
    lm_model <- lm(lm_formula, data = train_data)
    lm_pred <- predict(lm_model, newdata = test_data)
    # Calculate RMSE and R-squared for LM
    rmse_lm[i] <- sqrt(mean((test_data[[bioindicator]] - lm_pred)^2))
    mean_test <- mean(test_data[[bioindicator]])
    r2_lm[i] <- 1 - sum((test_data[[bioindicator]] - lm_pred)^2) / sum((test_data[[bioindicator]] - mean_test)^2)
    
    # --- GENERALIZED ADDITIVE MODEL (GAM) ---
    # Define formula for GAM using smoothed terms for predictors
    gam_formula <- as.formula(paste(bioindicator, "~ s(", 
                                    predictor, ", k=6) + s(Agricultural_area_frac, k=6) + s(Flow_velocity, k=6)"))
    gam_model <- mgcv::gam(gam_formula, data = train_data)
    gam_pred <- predict(gam_model, newdata = test_data)
    # Calculate RMSE and R-squared for GAM
    rmse_gam[i] <- sqrt(mean((test_data[[bioindicator]] - gam_pred)^2))
    r2_gam[i] <- 1 - sum((test_data[[bioindicator]] - gam_pred)^2) / sum((test_data[[bioindicator]] - mean_test)^2)
    
    # --- RANDOM FOREST (RF) ---
    # Define formula for RF using the bioindicator as response and selected predictors
    rf_formula <- as.formula(paste(bioindicator, "~", 
                                   predictor, "+ Agricultural_area_frac + Flow_velocity"))
    rf_model <- randomForest::randomForest(rf_formula, data = train_data, ntree = 500)
    rf_pred <- predict(rf_model, newdata = test_data)
    # Calculate RMSE and R-squared for RF
    rmse_rf[i] <- sqrt(mean((test_data[[bioindicator]] - rf_pred)^2))
    r2_rf[i] <- 1 - sum((test_data[[bioindicator]] - rf_pred)^2) / sum((test_data[[bioindicator]] - mean_test)^2)
  }
  
  # --- CALCULATE STANDARD ERRORS (SE) FOR EACH METRIC ---
  se_rmse_lm  <- sd(rmse_lm)  / sqrt(k_folds)
  se_r2_lm    <- sd(r2_lm)    / sqrt(k_folds)
  se_rmse_gam <- sd(rmse_gam) / sqrt(k_folds)
  se_r2_gam   <- sd(r2_gam)   / sqrt(k_folds)
  se_rmse_rf  <- sd(rmse_rf)  / sqrt(k_folds)
  se_r2_rf    <- sd(r2_rf)    / sqrt(k_folds)
  
  # Return a list with average metrics and their standard errors
  return(list(
    RMSE_LM = mean(rmse_lm),
    SE_RMSE_LM = se_rmse_lm,
    R2_LM = mean(r2_lm),
    SE_R2_LM = se_r2_lm,
    
    RMSE_GAM = mean(rmse_gam),
    SE_RMSE_GAM = se_rmse_gam,
    R2_GAM = mean(r2_gam),
    SE_R2_GAM = se_r2_gam,
    
    RMSE_RF = mean(rmse_rf),
    SE_RMSE_RF = se_rmse_rf,
    R2_RF = mean(r2_rf),
    SE_R2_RF = se_r2_rf
  ))
}

## 8.2. CV pipeline for all scenarios (time windows), bioindicators and chemical metrics selected ----

# Initialize a list to store results for each combination
results_3.5d_CV <- list()
# Loop through each scenario
for (scen in scenarios_3.5d) {
  # Filter data for the current scenario
  subset <- ds_3.5d_GAM_CV_RF %>% filter(scenario == scen)
  
  # Loop through each bioindicator
  for (bio in bioindicators_3.5d_CV) {
    # Loop through each chemical metric
    for (chem in chemical_metrics_3.5_CV) {
      # Perform cross-validation for this combination
      result <- cv_model_3.5d(subset, bio, chem, k_folds = 5)
      # Store results with a unique key (scenario|bioindicator|chemical_metric)
      results_3.5d_CV[[paste(scen, bio, chem, sep = "|")]] <- result
    }
  }
}
rm(subset)

## 8.3. Final 3.5d CV-RF dataset  ----
# Combine all results_3.5d_CV into a single data frame for easy analysis and visualization
results_3.5d_CV <- do.call(rbind, lapply(names(results_3.5d_CV), function(key) {
  # Split the key to extract scenario, bioindicator, and chemical metric
  parts <- strsplit(key, "|", fixed = TRUE)[[1]]
  
  # Create a data frame row for this combination
  data.frame(
    Scenario = parts[1],
    Bioindicator = parts[2],
    Chemical_Metric = parts[3],
    
    RMSE_LM = results_3.5d_CV[[key]]$RMSE_LM,
    SE_RMSE_LM = results_3.5d_CV[[key]]$SE_RMSE_LM,
    
    R2_LM = results_3.5d_CV[[key]]$R2_LM,
    SE_R2_LM = results_3.5d_CV[[key]]$SE_R2_LM,
    
    RMSE_GAM = results_3.5d_CV[[key]]$RMSE_GAM,
    SE_RMSE_GAM = results_3.5d_CV[[key]]$SE_RMSE_GAM,
    
    R2_GAM = results_3.5d_CV[[key]]$R2_GAM,
    SE_R2_GAM = results_3.5d_CV[[key]]$SE_R2_GAM,
    
    RMSE_RF = results_3.5d_CV[[key]]$RMSE_RF,
    SE_RMSE_RF = results_3.5d_CV[[key]]$SE_RMSE_RF,
    
    R2_RF = results_3.5d_CV[[key]]$R2_RF,
    SE_R2_RF = results_3.5d_CV[[key]]$SE_R2_RF)
}))


## 8.4. 14d best LMM model in CV + RF ----

### 8.4.1. Data set  ----
ds_14d_GAM_CV_RF <- df_metrics_abiotic_14d

### Define Bioindicators and Chemical Metrics to include in the model
bioindicators_14d_CV <- c("SPEAR")
chemical_metrics_14d_CV <- c("CRQmix_max")

### 8.4.2. Function: Cross-validation of two models (LMM, RF) ----
# Purpose: This function performs k-fold cross-validation for two types of models:
# - Linear Mixed Model (LM)
# - Random Forest (RF)
# It returns average performance metrics (RMSE and R-squared) and their standard errors.

cv_model_14d <- function(data_subset, bioindicator, predictor, k_folds = 5) {
  set.seed(123)  # Set seed for reproducibility of random splits
  
  # Create folds for cross-validation using the bioindicator as grouping variable
  folds <- caret::createFolds(data_subset[[bioindicator]], k = k_folds, list = TRUE)
  
  # Initialize vectors for performance metrics
  rmse_lmm <- numeric(k_folds)
  r2_lmm <- numeric(k_folds)
  rmse_rf <- numeric(k_folds)
  r2_rf <- numeric(k_folds)
  
  # Loop through each fold
  for (i in 1:k_folds) {
    test_idx <- folds[[i]]
    train_data <- data_subset[-test_idx, ]
    test_data <- data_subset[test_idx, ]
    
    # Linear Mixed Model (LMM)
    lmm_formula <- as.formula(paste(bioindicator, "~", 
                                    predictor, 
                                    "+ Agricultural_area_frac + Urban_area_frac + Flow_velocity + (1|YEAR)"))
    lmm_model <- lmer(lmm_formula, data = train_data)
    lmm_pred <- predict(lmm_model, newdata = test_data)
    rmse_lmm[i] <- sqrt(mean((test_data[[bioindicator]] - lmm_pred)^2))
    mean_test <- mean(test_data[[bioindicator]])
    r2_lmm[i] <- 1 - sum((test_data[[bioindicator]] - lmm_pred)^2) / sum((test_data[[bioindicator]] - mean_test)^2)
    
    # Random Forest (RF)
    rf_formula <- as.formula(paste(bioindicator, "~", 
                                   predictor, 
                                   "+ Agricultural_area_frac + Urban_area_frac + Flow_velocity + YEAR"))
    rf_model <- randomForest::randomForest(rf_formula, data = train_data, ntree = 500)
    rf_pred <- predict(rf_model, newdata = test_data)
    rmse_rf[i] <- sqrt(mean((test_data[[bioindicator]] - rf_pred)^2))
    r2_rf[i] <- 1 - sum((test_data[[bioindicator]] - rf_pred)^2) / sum((test_data[[bioindicator]] - mean_test)^2)
  }
  
  # --- CALCULATE STANDARD ERRORS (SE) FOR EACH METRIC ---
  se_rmse_lmm <- sd(rmse_lmm) / sqrt(k_folds)
  se_r2_lmm <- sd(r2_lmm) / sqrt(k_folds)
  se_rmse_rf <- sd(rmse_rf) / sqrt(k_folds)
  se_r2_rf <- sd(r2_rf) / sqrt(k_folds)
  
  # Return a list with average metrics and their standard errors
  return(list(
    RMSE_LMM = mean(rmse_lmm),
    SE_RMSE_LMM = se_rmse_lmm,
    R2_LMM = mean(r2_lmm),
    SE_R2_LMM = se_r2_lmm,
    
    RMSE_RF = mean(rmse_rf),
    SE_RMSE_RF = se_rmse_rf,
    R2_RF = mean(r2_rf),
    SE_R2_RF = se_r2_rf
  ))
}

## 8.5. CV pipeline for all scenarios (time windows), bioindicators and chemical metrics selected ----

# Initialize a list to store results for each combination
results_14d_CV <- list()

# Loop through each scenario
for (scen in scenarios_14d) {
  # Filter data for the current scenario
  subset <- ds_14d_GAM_CV_RF %>% filter(scenario == scen)
  # Loop through each bioindicator
  for (bio in bioindicators_14d_CV) {
    # Loop through each chem metric
    for (chem in chemical_metrics_14d_CV) {
      # Perform cross-validation for this combination
      result <- cv_model_14d(subset, bio, chem, k_folds = 5)
      # Store results with a unique key (scenario|bioindicator|chemical_metric)
      results_14d_CV[[paste(scen, bio, chem, sep = "|")]] <- result
    }
  }
}
rm(subset)

## 8.6. Final 14d CV-RF dataset  ----
# Combine all results_14d_CV into a single data frame for easy analysis and visualization
results_14d_CV <- do.call(rbind, lapply(names(results_14d_CV), function(key) {
  # Split the key to extract scenario, bioindicator, and chemical metric
  parts <- strsplit(key, "|", fixed = TRUE)[[1]]
  
  # Create a data frame row for this combination
  data.frame(
    Scenario = parts[1],
    Bioindicator = parts[2],
    Chemical_Metric = parts[3],
    
    RMSE_LMM = results_14d_CV[[key]]$RMSE_LMM,
    SE_RMSE_LMM = results_14d_CV[[key]]$SE_RMSE_LMM,
    R2_LMM = results_14d_CV[[key]]$R2_LMM,
    SE_R2_LMM = results_14d_CV[[key]]$SE_R2_LMM,
    
    RMSE_RF = results_14d_CV[[key]]$RMSE_RF,
    SE_RMSE_RF = results_14d_CV[[key]]$SE_RMSE_RF,
    R2_RF = results_14d_CV[[key]]$R2_RF,
    SE_R2_RF = results_14d_CV[[key]]$SE_R2_RF)}))


## 8.7. CV + RF Comparison plots ----

### Add the sampling period to each sub-set of CV data
results_14d_CV$Period <- "14 days"
results_3.5d_CV$Period <- "3.5 days"

## Join of both results
all_results_CV <- bind_rows(results_3.5d_CV, results_14d_CV)

## Filter of the "main data"1_week" time window in the 3.5d sampling period to improve the
## visualization of the plot, also mentioned in the caption
filtered_results_CV <- all_results_CV %>% 
  filter(Scenario != "3.5_days_1_week")

## Re-shaping the data set
long_results_CV <- filtered_results_CV %>%
  select(Scenario, Period, RMSE_LMM, SE_RMSE_LMM, R2_LMM, SE_R2_LMM, RMSE_RF, SE_RMSE_RF, R2_RF, SE_R2_RF,
         RMSE_LM, SE_RMSE_LM, R2_LM, SE_R2_LM, RMSE_GAM, SE_RMSE_GAM, R2_GAM, SE_R2_GAM) %>%
  gather(key = "Metric_Model", value = "Value", 
         RMSE_LMM, R2_LMM, RMSE_RF, R2_RF, RMSE_LM, R2_LM, RMSE_GAM, R2_GAM) %>%
  mutate(
    Validation = ifelse(grepl("RMSE", Metric_Model), "RMSE", "R²"),
    Model = gsub("RMSE_|R2_|_LMM|_RF|_LM|_GAM", "", Metric_Model),
    Model = case_when(
      grepl("LMM", Metric_Model) ~ "LMM",
      grepl("RF", Metric_Model) ~ "RF",
      grepl("LM", Metric_Model) ~ "LM",
      grepl("GAM", Metric_Model) ~ "GAM",
      TRUE ~ Model))

long_results_CV <- long_results_CV %>%
  mutate(
    SE = case_when(
      Metric_Model == "RMSE_LMM" ~ filtered_results_CV$SE_RMSE_LMM[match(Scenario, filtered_results_CV$Scenario)],
      Metric_Model == "R2_LMM" ~ filtered_results_CV$SE_R2_LMM[match(Scenario, filtered_results_CV$Scenario)],
      Metric_Model == "RMSE_RF" ~ filtered_results_CV$SE_RMSE_RF[match(Scenario, filtered_results_CV$Scenario)],
      Metric_Model == "R2_RF" ~ filtered_results_CV$SE_R2_RF[match(Scenario, filtered_results_CV$Scenario)],
      Metric_Model == "RMSE_LM" ~ filtered_results_CV$SE_RMSE_LM[match(Scenario, filtered_results_CV$Scenario)],
      Metric_Model == "R2_LM" ~ filtered_results_CV$SE_R2_LM[match(Scenario, filtered_results_CV$Scenario)],
      Metric_Model == "RMSE_GAM" ~ filtered_results_CV$SE_RMSE_GAM[match(Scenario, filtered_results_CV$Scenario)],
      Metric_Model == "R2_GAM" ~ filtered_results_CV$SE_R2_GAM[match(Scenario, filtered_results_CV$Scenario)],
      TRUE ~ NA_real_))

long_results_CV$Period <- factor(long_results_CV$Period, levels = c("3.5 days", "14 days"))
long_results_CV$Validation <- factor(long_results_CV$Validation, levels = c("RMSE", "R²"))

## Organization of the time windows for the plot
ordered_scenarios_CV <- c(
  "14_days_1_week", "14_days_2_weeks", "14_days_1_month", "14_days_2_months", 
  "14_days_3_months", "14_days_6_months", "14_days_1_year",
  "3.5_days_2_weeks", "3.5_days_1_month", "3.5_days_2_months",
  "3.5_days_3_months", "3.5_days_6_months", "3.5_days_1_year")

## Forcing the order through factor conversion
long_results_CV$Scenario <- factor(long_results_CV$Scenario, levels = ordered_scenarios_CV)

## Specific color pallete
palette_CV <- c("#1f77b4", "#db6363", "gray25", "#59944d")

## Plot
ggplot(long_results_CV, aes(x = Scenario, y = Value, color = Model, shape = Model, group = Model)) +
  geom_point(position = position_dodge(width = 0.5), size = 2) +
  geom_errorbar(aes(ymin = Value - SE, ymax = Value + SE), 
                width = 0.2, position = position_dodge(width = 0.5)) +
  facet_grid(Validation ~ Period, scales = "free") +
  labs(
    title = "Model Validation Metrics Across Sampling Periods and Time Windows",
    subtitle = 
    "3.5-DAYS: SPEAR ~ Agriculture proportion + Flow velocity + TU_ECmix_median\n14-DAYS: SPEAR  ~ Agriculture proportion + Urban proportion + Flow velocity + CRQmix_max",
    y = "Value", x = "Time windows",
    caption = "Note: The time-window '3.5_days_1_week' was excluded due to extreme values that distort visualization.\nRMSE = Root Mean Squared Error; R² = Coefficient of Determination."
  ) +
  scale_color_manual(values = palette_CV) +
  scale_shape_manual(values = c(16, 17, 18, 15, 3)) + 
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.caption = element_text(hjust = 0),
    strip.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold")) -> p

## Export as png
#png(paste(user, output_graphics_file ,"CV_RF_3.5d&14d.png", sep = ""), 
#    width = 25, height = 20, units = "cm", res = 350)
p
#dev.off()
#ggsave("CV_RF_3.5d&14d_improved.png", p, width = 25, height = 20, units = "cm", dpi = 350)
rm(p)


# END :)





