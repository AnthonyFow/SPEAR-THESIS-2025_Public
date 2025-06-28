# _____________________________________________________________
# SYNTHETIC DATASET GENERATION SCRIPT
# Author: Anthony Jeanpier Fow Esteves
# Date: 28.06.2024
# _____________________________________________________________

# Description:
# This script generates synthetic datasets for the thesis analysis, closely mimicking the structure,
# column names, and approximate number of observations and streams of the original data.
# All values are randomly generated and do NOT correspond to real observations.
#
# - You may freely modify this script to change the number of streams, scenarios, or value distributions.
# - The resulting `.rds` files in `/inputs` allow users to run the main statistical analysis script
#   without access to confidential data.
#
# For details on variable meanings and dataset structure, refer to the `/docs` folder.


# Set your seed to make it replicable ----
set.seed(1996)

# 1. Define the name of your 40 streams (you can change the number) ----
n_streams <- 40
stream_names <- paste0("Stream_", seq_len(n_streams))

# 2. Abiotic_factors_41_streams_ds_prepared: one row per stream ----
Abiotic_factors_41_streams_ds_prepared <- tibble(
  Gewasser = stream_names,
  Messstelle_ID = factor(sample(1000:7000, n_streams, replace = FALSE)),
  Catch_area_m2 = round(runif(n_streams, 1e6, 5e6)),
  Avg_mod_discharge_m3_s = round(runif(n_streams, 0.05, 2), 2),
  Proportion_of_wastewater_l_s = round(runif(n_streams, 0, 50), 1),
  FLOZ = sample(1:5, n_streams, replace = TRUE),
  Flow_velocity = round(runif(n_streams, 0.1, 1.5), 3),
  Ecomorphology_class = sample(1:5, n_streams, replace = TRUE),
  Ecomorphology_categories23 = factor(sample(
    c("Eingedolt", "Naturfremd, künstlich", "Stark beeinträchtigt", "Mässig beeinträchtigt", "Wenig beeinträchtigt", "Naturnah"),
    n_streams, replace = TRUE)),
  Ecomorphology_values23 = factor(sample(1:6, n_streams, replace = TRUE)),
  Urban_area_frac = round(runif(n_streams, 0, 50), 1),
  Forest_area_frac = round(runif(n_streams, 0, 80), 1),
  Agricultural_area_frac = round(runif(n_streams, 0, 90), 1),
  Mod_max_temp_summer = round(runif(n_streams, 12, 20), 1),
  Ecomorphology_cont23_Langhans = round(runif(n_streams, 0, 1), 2),
  Ecomorphology_0_12 = round(runif(n_streams, 0, 12), 0),
  Stream_bed_construction = factor(sample(
    c("0%, natural", "20%, partly artificial", "40%, mostly artificial", "60%, artificial", "80%, nearly all artificial", "100%, artificial"), 
    n_streams, replace = TRUE)))

saveRDS(Abiotic_factors_41_streams_ds_prepared, "inputs/Abiotic_factors_41_streams_ds_prepared.rds")

# Scenarios (or time windows, they state that period before each macroinvertebrate 
# monitoring from which the chemical samples were grouped to calculate the chemical
# metrics (e.g., ARQ, TU-NOEC).
scenarios_3.5d <- c("3.5_days_1_week", "3.5_days_2_weeks", "3.5_days_1_month", "3.5_days_2_months", "3.5_days_3_months", "3.5_days_6_months", "3.5_days_1_year")
scenarios_14d  <- c("14_days_1_week", "14_days_2_weeks", "14_days_1_month", "14_days_2_months", "14_days_3_months", "14_days_6_months", "14_days_1_year")

# 3. metrics_14d_0_R12_pp_ds: 1200 rows, it includes all the streams ----
n_14d <- 1200
metrics_14d_0_R12_pp_ds <- tibble(
  Gewasser = sample(stream_names, n_14d, replace = TRUE),
  OBSERVATIONDATE = as.POSIXct("2022-03-02") + sample(0:365, n_14d, replace = TRUE) * 86400,
  CRQmix_max = runif(n_14d, 0, 10),
  CRQmix_max_Probe = factor(sample(4000:5000, n_14d, replace = TRUE)),
  CRQmix_max_DIF_DATE = sample(c("Less than 1 week", "1-2 weeks", "More than 2 weeks"), n_14d, replace = TRUE),
  CRQmix_max_num_Substances = sample(1:70, n_14d, replace = TRUE),
  CRQmix_mean = runif(n_14d, 0, 5),
  CRQmix_median = runif(n_14d, 0, 5),
  CRQmix_n = sample(1:5, n_14d, replace = TRUE),
  CRQmix_mean_num_substances = runif(n_14d, 1, 70),
  CRQmax_max = runif(n_14d, 0, 10),
  CRQmax_max_Probe = factor(sample(4000:5000, n_14d, replace = TRUE)),
  CRQmax_max_DIF_DATE = sample(c("Less than 1 week", "1-2 weeks", "More than 2 weeks"), n_14d, replace = TRUE),
  CRQmax_max_Param = sample(c("Azoxystrobin", "Imidacloprid", "Diuron"), n_14d, replace = TRUE),
  CRQmax_max_Substance = factor(sample(1:3000, n_14d, replace = TRUE)),
  CRQmax_mean = runif(n_14d, 0, 5),
  CRQmax_median = runif(n_14d, 0, 5),
  CRQmax_n = sample(1:5, n_14d, replace = TRUE),
  TU_NOECmix_max = runif(n_14d, 0, 0.02),
  TU_NOECmix_max_Probe = factor(sample(4000:5000, n_14d, replace = TRUE)),
  TU_NOECmix_max_DIF_DATE = sample(c("Less than 1 week", "1-2 weeks", "More than 2 weeks"), n_14d, replace = TRUE),
  TU_NOECmix_max_num_Substances = sample(1:70, n_14d, replace = TRUE),
  TU_NOECmix_mean = runif(n_14d, 0, 0.01),
  TU_NOECmix_median = runif(n_14d, 0, 0.01),
  TU_NOECmix_n = sample(1:5, n_14d, replace = TRUE),
  TU_NOECmix_mean_num_substances = runif(n_14d, 1, 70),
  TU_NOECmax_max = runif(n_14d, 0, 0.02),
  TU_NOECmax_max_Probe = factor(sample(4000:5000, n_14d, replace = TRUE)),
  TU_NOECmax_max_DIF_DATE = sample(c("Less than 1 week", "1-2 weeks", "More than 2 weeks"), n_14d, replace = TRUE),
  TU_NOECmax_max_Param = sample(c("Azoxystrobin", "Imidacloprid", "Diuron"), n_14d, replace = TRUE),
  TU_NOECmax_max_Substance = factor(sample(1:3000, n_14d, replace = TRUE)),
  TU_NOECmax_mean = runif(n_14d, 0, 0.01),
  TU_NOECmax_median = runif(n_14d, 0, 0.01),
  TU_NOECmax_n = sample(1:5, n_14d, replace = TRUE),
  SPEAR = runif(n_14d, 0, 30),
  IBCH_2019 = runif(n_14d, 0, 1),
  GI_VALUE = runif(n_14d, 0, 1),
  VT_VALUE = runif(n_14d, 0, 1),
  EPT = runif(n_14d, 0, 1),
  Kanton = factor(sample(c("AG", "BE", "BL", "BS", "FR", "GE", "GL", "GR", "JU", "LU", "NE", "NW", "OW", "SG", "SH", "SO", "SZ"), n_14d, replace = TRUE)),
  Kategorie_FG_Grosse = factor(sample(c("Grosser Fluss", "Mittlerer Fluss", "Kleiner Fluss", "Bach", "Quellbach"), n_14d, replace = TRUE)),
  YEAR = sample(2018:2022, n_14d, replace = TRUE),
  OID = factor(sample(paste0("CH_", sprintf("%03d", 1:49), "_BS"), n_14d, replace = TRUE)),
  X.y = runif(n_14d, 600000, 800000),
  Y.y = runif(n_14d, 200000, 300000),
  Z = runif(n_14d, 200, 500),
  Season_MZB = factor(sample(c("Spring", "Summer"), n_14d, replace = TRUE)),
  scenario = factor(sample(scenarios_14d, n_14d, replace = TRUE)),
  Gewasser_clean = tolower(Gewasser),
  PP_mean_mean = runif(n_14d, 0, 2))

saveRDS(metrics_14d_0_R12_pp_ds, "inputs/metrics_14d_0_R12_pp_ds.rds")

# 4. metrics_3.5d_0_R12_pp_ds: 300 rows, not all the streams available (as in the thesis) ----
n_3d <- 300
subset_streams <- sample(stream_names, 15) # e.g., 15 out of 40
metrics_3.5d_0_R12_pp_ds <- tibble(
  Gewasser = sample(subset_streams, n_3d, replace = TRUE),
  OBSERVATIONDATE = as.POSIXct("2020-07-27") + sample(0:365, n_3d, replace = TRUE) * 86400,
  ARQmix_max = runif(n_3d, 0, 10),
  ARQmix_max_Probe = factor(sample(4000:5000, n_3d, replace = TRUE)),
  ARQmix_max_DIF_DATE = sample(c("Less than 1 week", "1-2 weeks", "More than 2 weeks"), n_3d, replace = TRUE),
  ARQmix_max_num_Substances = sample(1:70, n_3d, replace = TRUE),
  ARQmix_mean = runif(n_3d, 0, 5),
  ARQmix_median = runif(n_3d, 0, 5),
  ARQmix_n = sample(1:5, n_3d, replace = TRUE),
  ARQmix_mean_num_substances = runif(n_3d, 1, 70),
  ARQmax_max = runif(n_3d, 0, 10),
  ARQmax_max_Probe = factor(sample(4000:5000, n_3d, replace = TRUE)),
  ARQmax_max_DIF_DATE = sample(c("Less than 1 week", "1-2 weeks", "More than 2 weeks"), n_3d, replace = TRUE),
  ARQmax_max_Param = sample(c("Azoxystrobin", "Imidacloprid", "Diuron"), n_3d, replace = TRUE),
  ARQmax_max_Substance = factor(sample(1:3000, n_3d, replace = TRUE)),
  ARQmax_mean = runif(n_3d, 0, 5),
  ARQmax_median = runif(n_3d, 0, 5),
  ARQmax_n = sample(1:5, n_3d, replace = TRUE),
  TU_ECmix_max = runif(n_3d, 0, 0.02),
  TU_ECmix_max_Probe = factor(sample(4000:5000, n_3d, replace = TRUE)),
  TU_ECmix_max_DIF_DATE = sample(c("Less than 1 week", "1-2 weeks", "More than 2 weeks"), n_3d, replace = TRUE),
  TU_ECmix_max_num_Substances = sample(1:70, n_3d, replace = TRUE),
  TU_ECmix_mean = runif(n_3d, 0, 0.01),
  TU_ECmix_median = runif(n_3d, 0, 0.01),
  TU_ECmix_n = sample(1:5, n_3d, replace = TRUE),
  TU_ECmix_mean_num_substances = runif(n_3d, 1, 70),
  TU_ECmax_max = runif(n_3d, 0, 0.02),
  TU_ECmax_max_Probe = factor(sample(4000:5000, n_3d, replace = TRUE)),
  TU_ECmax_max_DIF_DATE = sample(c("Less than 1 week", "1-2 weeks", "More than 2 weeks"), n_3d, replace = TRUE),
  TU_ECmax_max_Param = sample(c("Azoxystrobin", "Imidacloprid", "Diuron"), n_3d, replace = TRUE),
  TU_ECmax_max_Substance = factor(sample(1:3000, n_3d, replace = TRUE)),
  TU_ECmax_mean = runif(n_3d, 0, 0.01),
  TU_ECmax_median = runif(n_3d, 0, 0.01),
  TU_ECmax_n = sample(1:5, n_3d, replace = TRUE),
  SPEAR = runif(n_3d, 0, 30),
  IBCH_2019 = runif(n_3d, 0, 1),
  GI_VALUE = runif(n_3d, 0, 1),
  VT_VALUE = runif(n_3d, 0, 1),
  EPT = runif(n_3d, 0, 1),
  Kanton = factor(sample(c("AG", "BE", "BL", "BS", "FR", "GE", "GL", "GR", "JU", "LU", "NE", "NW", "OW", "SG", "SH", "SO", "SZ"), n_3d, replace = TRUE)),
  Kategorie_FG_Grosse = factor(sample(c("Grosser Fluss", "Mittlerer Fluss", "Kleiner Fluss", "Bach", "Quellbach"), n_3d, replace = TRUE)),
  YEAR = sample(2018:2022, n_3d, replace = TRUE),
  OID = factor(sample(paste0("CH_", sprintf("%03d", 1:49), "_BS"), n_3d, replace = TRUE)),
  X.y = runif(n_3d, 600000, 800000),
  Y.y = runif(n_3d, 200000, 300000),
  Z = runif(n_3d, 200, 500),
  Season_MZB = factor(sample(c("Spring", "Summer"), n_3d, replace = TRUE)),
  scenario = factor(sample(scenarios_3.5d, n_3d, replace = TRUE)),
  Gewasser_clean = tolower(Gewasser),
  PP_mean_mean = runif(n_3d, 0, 2))

saveRDS(metrics_3.5d_0_R12_pp_ds, "inputs/metrics_3.5d_0_R12_pp_ds.rds")

