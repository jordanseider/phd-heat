###
# Identification of reindeer fine-scale foraging behaviour using tri-axial accelerometer data
# Heidi Rautiainen, with acknowledgements to Måns Karlsson for hidden Markov model
# 2020-08-26
# https://doi.org/10.5061/dryad.8sf7m0cs7
#
# Updated code to implement more modern packages for HMM
###

library(tidyverse)
library(data.table)

setwd("C:/Users/jseider.stu/Sync")

## DATA CLEANING
# ACceleration data
acc <- read.csv("./Data/SAMPLE_DATASETS/Accelerometer/acceleration.csv") %>% 
  mutate(Timestamp = as.POSIXct(Timestamp, tz = "UTC"))

# Annotation data (behaviours)
ann <- read.csv("./Data/SAMPLE_DATASETS/Accelerometer/annotations.csv") %>% 
  rename(Observation_ID = Observation.id,
         Behavioral_category = Behavioral.category,
         Behavior_type = Behavior.type,
         Start_s = Start..s., 
         Stop_s = Stop..s., 
         Duration_s = Duration..s., 
         Observation_date = Observation.date) %>% 
  mutate(Observation_date = as.Date(Observation_date, format = "%Y-%m-%d"),
         Start = as.POSIXct(as.numeric(Start_s),
                            origin = "2020-02-27",
                            tz = "UTC"),
         Stop = as.POSIXct(as.numeric(Stop_s),
                           origin = "2020-02-27",
                           tz = "UTC"),
         Modifiers = case_when(is.na(Modifiers) ~ Behavior,
                               .default = Modifiers))

# Add behaviour to timestamped accelerometer data
data <- acc %>% 
  left_join(
    ann %>% select(Start, Stop, Modifiers),
    join_by(Timestamp >= Start, Timestamp <= Stop)
  ) %>% 
  select(-Start, -Stop)

## Accelerometer filtering and data processing
# Vectorized Static & Dynamic Acceleration (5-second window at 10 Hz = 51 samples)
k_win <- 51

data <- data %>%
  mutate(
    stX = frollmean(X, n = k_win, fill = NA, align = "center"),
    stY = frollmean(Y, n = k_win, fill = NA, align = "center"),
    stZ = frollmean(Z, n = k_win, fill = NA, align = "center"),
    
    dyX = X - stX,
    dyY = Y - stY,
    dyZ = Z - stZ,
    
    VeDBA = sqrt(dyX^2 + dyY^2 + dyZ^2), # Crucial for HMM state emissions
    roll  = atan2(stY, stZ),
    pitch = atan2(-stX, sqrt(stY^2 + stZ^2))
    
  ) %>%
  filter(!is.na(stX)) # Drop edge NAs from rolling mean

# Vectorized Coordinate Rotation
data <- data %>%
  mutate(
    rot_X = X,
    rot_Y =  Y * cos(roll) + Z * sin(roll),
    rot_Z = -Y * sin(roll) + Z * cos(roll),
    
    rot_stX = stX,
    rot_stY =  stY * cos(roll) + stZ * sin(roll),
    rot_stZ = -stX * sin(roll) + stZ * cos(roll),
    
    rot_dyX = dyX,
    rot_dyY =  dyY * cos(roll) + dyZ * sin(roll),
    rot_dyZ = -dyX * sin(roll) + dyZ * cos(roll)
  )

## Aggregate data to 5-second window
get_mode <- function(x) {
  ux <- unique(na.omit(x))
  if(length(ux) == 0) return(NA)
  ux[which.max(tabulate(match(x, ux)))]
}

# ==============================================================================
# VARIABLE DEFINITIONS
# ==============================================================================

# See Table A3 from https://pmc.ncbi.nlm.nih.gov/articles/PMC11097004/#ece311380-app-0001 

# --- WINDOW & METADATA ---
# samp      : POSIXct timestamp marking the exact start of the 2-second window
# Timestamp : Formatted character string of the window's starting timestamp
# ID        : Unique Tag ID / Identifier for the individual animal
# Modifiers : Dominant (most frequent) behavioral classification during the window

# --- ROTATED RAW ACCELERATION (m = Median, d = Interquartile Range) ---
# mrot_X    : Median rotated raw acceleration (X-axis / surge / forward-backward)
# drot_X    : IQR of rotated raw acceleration (X-axis)
# mrot_Y    : Median rotated raw acceleration (Y-axis / sway / side-to-side)
# drot_Y    : IQR of rotated raw acceleration (Y-axis)
# mrot_Z    : Median rotated raw acceleration (Z-axis / heave / up-down)
# drot_Z    : IQR of rotated raw acceleration (Z-axis)

# --- ROTATED STATIC ACCELERATION (Gravity component) ---
# mrot_stX  : Median rotated static acceleration (X-axis)
# drot_stX  : IQR of rotated static acceleration (X-axis)
# mrot_stY  : Median rotated static acceleration (Y-axis)
# drot_stY  : IQR of rotated static acceleration (Y-axis)
# mrot_stZ  : Median rotated static acceleration (Z-axis)
# drot_stZ  : IQR of rotated static acceleration (Z-axis)

# --- ROTATED DYNAMIC ACCELERATION (Animal movement component) ---
# mrot_dyX  : Median rotated dynamic acceleration (X-axis)
# drot_dyX  : IQR of rotated dynamic acceleration (X-axis)
# mrot_dyY  : Median rotated dynamic acceleration (Y-axis)
# drot_dyY  : IQR of rotated dynamic acceleration (Y-axis)
# mrot_dyZ  : Median rotated dynamic acceleration (Z-axis)
# drot_dyZ  : IQR of rotated dynamic acceleration (Z-axis)

# --- EXTENDED STATIC STATISTICS ---
# meanrot_stx, meanrot_sty, meanrot_stz : Mean rotated static acceleration (X, Y, Z)
# minrot_stx,  minrot_sty,  minrot_stz  : Minimum rotated static acceleration (X, Y, Z)
# maxrot_stx,  maxrot_sty,  maxrot_stz  : Maximum rotated static acceleration (X, Y, Z)
# sdrot_stx,   sdrot_sty,   sdrot_stz   : Standard deviation of rotated static acceleration (X, Y, Z)

# --- EXTENDED DYNAMIC STATISTICS ---
# meanrot_dyx, meanrot_dyy, meanrot_dyz : Mean rotated dynamic acceleration (X, Y, Z)
# minrot_dyx,  minrot_dyy,  minrot_dyz  : Minimum rotated dynamic acceleration (X, Y, Z)
# maxrot_dyx,  maxrot_dyy,  maxrot_dyz  : Maximum rotated dynamic acceleration (X, Y, Z)
# sdrot_dyx,   sdrot_dyy,   sdrot_dyz   : Standard deviation of rotated dynamic acceleration (X, Y, Z)

# --- POSTURE ANGLES (Derived from static acceleration) ---
# mroll     : Median roll angle (rotation around the X-axis)
# droll     : IQR of roll angle
# mpitch    : Median pitch angle (rotation around the Y-axis)
# dpitch    : IQR of pitch angle

data_5s <- data %>% 
  mutate(window_id = floor_date(Timestamp, "5 seconds")) %>%
  group_by(ID = TagID, window_id) %>%
  filter(n() >= 45) %>% 
  summarize(
    # VeDBA & Pitch features
    mVeDBA = mean(VeDBA, na.rm = TRUE),
    sdVeDBA = sd(VeDBA, na.rm = TRUE),
    mpitch = median(pitch, na.rm = TRUE),
    sdpitch = sd(pitch, na.rm = TRUE),
    mroll  = median(roll, na.rm = TRUE),
    
    # Static & Dynamic Rotated Aggregates
    meanrot_stx = mean(rot_stX, na.rm = TRUE),
    meanrot_sty = mean(rot_stY, na.rm = TRUE),
    meanrot_stz = mean(rot_stZ, na.rm = TRUE),
    sdrot_dyx   = sd(rot_dyX, na.rm = TRUE),
    sdrot_dyy   = sd(rot_dyY, na.rm = TRUE),
    sdrot_dyz   = sd(rot_dyZ, na.rm = TRUE),
    
    # Preserved Label (NA when unobserved)
    Modifier_Label = get_mode(Modifiers),
    .groups = "drop"
  ) %>%
  filter(!Modifier_Label %in% c("Missing data", "High", "While walking")) %>%
  mutate(
    
    Modifier_Label = case_when(
      Modifier_Label == "Agnostic behavior" ~ "Other",
      Modifier_Label == "Hole" ~ "Ground",
      .default = as.character(Modifier_Label)),
    
    behaviour = case_when(
      Modifier_Label %in% c("Standing", "Sleeping", "Lying", "Ruminating", "Resting") ~ "Stationary",
      Modifier_Label %in% c("Grazing", "Ground")              ~ "Grazing",
      Modifier_Label %in% c("Walking", "Trotting", "Running") ~ "Moving",
      is.na(Modifier_Label) ~ NA,
      TRUE ~ NA),
    
    hour_num = hour(window_id) + minute(window_id) / 60 + second(window_id) / 3600,
    cos_tod  = cos(2 * pi * hour_num / 24),
    sin_tod  = sin(2 * pi * hour_num / 24)
    
    )

state_levels <- sort(unique(na.omit(data_5s$behaviour)))
data_5s$known_state <- match(data_5s$behaviour, state_levels)

data_5s$mVeDBA <- ifelse(data_5s$mVeDBA == 0, 1e-6, data_5s$mVeDBA)
data_5s$sdrot_dyx <- ifelse(data_5s$sdrot_dyx == 0, 1e-6, data_5s$sdrot_dyx)

data_hmm <- prepData(
  data = data_5s,
  type = "UTM",
  coordNames = NULL,
  covNames = c("cos_tod", "sin_tod")
)

empirical_stats <- data_5s %>%
  filter(!is.na(known_state)) %>%
  group_by(known_state) %>%
  summarize(
    # Feature 1: mVeDBA (Gamma distribution)
    mean_mVeDBA = mean(mVeDBA, na.rm = TRUE),
    sd_mVeDBA   = sd(mVeDBA, na.rm = TRUE),
    
    # Feature 2: mpitch (Normal distribution)
    mean_mpitch = mean(mpitch, na.rm = TRUE),
    sd_mpitch   = sd(mpitch, na.rm = TRUE),
    
    # Feature 3: sdrot_dyx (Gamma distribution)
    mean_sdrot  = mean(sdrot_dyx, na.rm = TRUE),
    sd_sdrot    = sd(sdrot_dyx, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  arrange(known_state)

Par0 <- list(
  mVeDBA    = c(empirical_stats$mean_mVeDBA,  empirical_stats$sd_mVeDBA),
  mpitch    = c(empirical_stats$mean_mpitch,  empirical_stats$sd_mpitch),
  sdrot_dyx = c(empirical_stats$mean_sdrot,   empirical_stats$sd_sdrot)
)

## 5-fold contiguous block cross-validation
labeled_indices <- which(!is.na(data_5s$known_state))

n_folds <- 5
block_size <- ceiling(length(labeled_indices) / n_folds)
fold_assignments <- rep(1:n_folds, each = block_size)[1:length(labeled_indices)]

data_5s$cv_predicted_state <- NA_character_

for (k in 1:n_folds) {
  # Identify test indices for the current block
  test_idx <- labeled_indices[fold_assignments == k]
  
  # Mask test labels as NA for training
  train_known_states <- data_5s$known_state
  train_known_states[test_idx] <- NA
  
  # Fit PHMM on training sequence (with current block masked)
  phmm_block_fit <- fitHMM(
    data        = data_hmm,
    nbStates    = length(state_levels),
    dist        = list(mVeDBA    = "gamma", 
                       mpitch    = "norm", 
                       sdrot_dyx = "gamma"),
    Par0        = Par0,
    #formula     = ~ cos_tod + sin_tod,
    knownStates = train_known_states,
    stateNames  = state_levels
  )
  
  decoded_path <- viterbi(phmm_block_fit)
  data_5s$cv_predicted_state[test_idx] <- state_levels[decoded_path[test_idx]]
}

test_set <- data_5s %>% filter(!is.na(cv_predicted_state))

real_block_conf_matrix <- confusionMatrix(
  data      = factor(test_set$cv_predicted_state, levels = state_levels),
  reference = factor(test_set$behaviour,          levels = state_levels)
)

print(real_block_conf_matrix)

plot(phmm_block_fit, ask = FALSE)
