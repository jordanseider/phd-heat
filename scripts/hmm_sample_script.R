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
                               Modifiers == "Ruminating" ~ "Resting",
                               Modifiers == "Hole" ~ "Ground",
                               .default = Modifiers))

# Add behaviour to timestamped accelerometer data
data <- acc %>% 
  left_join(
    ann %>% select(Start, Stop, Modifiers),
    join_by(Timestamp >= Start, Timestamp <= Stop)
  ) %>% 
  select(-Start, -Stop) %>% 
  na.omit()

# Summary of data
sum_table <- table(data$Modifiers)
behaviours <- names(sum_table)

r1 <- rle(data$Modifiers)
cr <- cumsum(r1$lengths)


## Accelerometer filtering and data processing
library(chron)

hz <- 10 # sampling frequency
dur <- 5 # window size (seconds) for calculating running mean
AccHz <- hz * dur / 2 # half-window size: see e.g., McClune DW et al. 2014 Fehlmann et al 2017

# Function to calculates a centered moving average for a numeric vector (such as an acceleration axis) 
# using a rolling window of a specified duration
mov.av <- function(obj){
  
  ma <- obj
  for(i in AccHz:length(obj)){
    ma[i] <- mean(obj[(i - AccHz):(i + AccHz)], na.rm = TRUE)
  }

  # Remove edge cases at beginning and end
  for(i in 1:AccHz){
    ma[i] <- NA
  }
  for(i in (length(obj) - AccHz):length(obj)){
    ma[i] <- NA
  }
  
  return(ma)
}

## Calculate static acceleration
data$stX <- mov.av(data$X)
data$stY <- mov.av(data$Y)
data$stZ <- mov.av(data$Z)

## Cacluate dynamic acceleration (DBA) (Dy)
data$dyX <- data$X - data$stX
data$dyY <- data$Y - data$stY
data$dyZ <- data$Z - data$stZ

## Roll
roll <- function(y, z){
  roll = atan2(y, z)
  return(roll)
}

data$roll <- roll(data$stY, data$stZ)

## Pitch
pitch <- function(x, y, z){
  pitch90 = 1 / atan(-x / sqrt(y^2 + z^2))
  return(pitch90)
}

data$pitch <- pitch(data$stX, data$stY, data$stZ)

## Rotation matrix

# This function takes an angel in radians and returns the
# matrix that you multiply a vector with to rotate the vector alpha radians
# around the x-axis.

x_rot <- function(alpha){
  
  R <- matrix(nrow = 3,
              ncol = 3,
              data = c(1, 0, 0,
                       0, cos(alpha), -sin(alpha),
                       0, sin(alpha), cos(alpha)))
  return(R)
  
}

rotation_matrices <- lapply(data$roll, x_rot)

# Compute corrected acceleration values
rotated <- t(sapply(1:length(rotation_matrices), 
                    function(i){
                      return(solve(rotation_matrices[[i]]) %*% c(data$X[i], data$Y[i], data$Z[i]))
                    }))

rotated_static <- t(sapply(1:length(rotation_matrices), 
                           function(i){
                             return(solve(rotation_matrices[[i]]) %*% c(data$stX[i], data$stY[i], data$stZ[i]))
                           }))

rotated_dynamic <- t(sapply(1:length(rotation_matrices), 
                            function(i){
                              return(solve(rotation_matrices[[i]]) %*% c(data$dyX[i], data$dyY[i], data$dyZ[i]))
                            }))

data$rot_X <- rotated[, 1]
data$rot_Y <- rotated[, 2]
data$rot_Z <- rotated[, 3]

data$rot_stX <- rotated_static[, 1]
data$rot_stY <- rotated_static[, 2]
data$rot_stZ <- rotated_static[, 3]

data$rot_dyX <- rotated_dynamic[, 1]
data$rot_dyY <- rotated_dynamic[, 2]
data$rot_dyZ <- rotated_dynamic[, 3]

data <- na.omit(data)

## Aggregate data to 2-second window
mnt <- min(data$Timestamp, na.rm = TRUE)
mxt <- max(data$Timestamp, na.rm = TRUE)

chunks <- seq(from = mnt, to = (mxt - 2), by = 2)

# random selection of chunks
set.seed(1453)

start_idx <- sample(1:2, size = 1) 
samp_times <- chunks[seq(start_idx, length(chunks), by = 1)]

windows <- tibble(
  window_id = seq_along(samp_times),
  samp_start = samp_times,
  samp_end = samp_times + 2
)

get_mode <- function(x) {
  ux <- unique(na.omit(x))
  if(length(ux) == 0) return(NA)
  ux[which.max(tabulate(match(x, ux)))]
}

# ==============================================================================
# DATA_2S VARIABLE DEFINITIONS
# Aggregated 2-second window statistics for tri-axial accelerometer data
# ==============================================================================

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

data_2s <- windows %>% 
  inner_join(data, join_by(samp_start <= Timestamp, samp_end > Timestamp)) %>%
  group_by(window_id, samp = samp_start) %>%
  filter(n() > 19) %>% 
  summarize(
    # Median and IQR (Dispersion)
    mrot_X = median(rot_X, na.rm = TRUE), drot_X = IQR(rot_X, na.rm = TRUE),
    mrot_Y = median(rot_Y, na.rm = TRUE), drot_Y = IQR(rot_Y, na.rm = TRUE),
    mrot_Z = median(rot_Z, na.rm = TRUE), drot_Z = IQR(rot_Z, na.rm = TRUE),
    
    mrot_stX = median(rot_stX, na.rm = TRUE), drot_stX = IQR(rot_stX, na.rm = TRUE),
    mrot_stY = median(rot_stY, na.rm = TRUE), drot_stY = IQR(rot_stY, na.rm = TRUE),
    mrot_stZ = median(rot_stZ, na.rm = TRUE), drot_stZ = IQR(rot_stZ, na.rm = TRUE),
    
    mrot_dyX = median(rot_dyX, na.rm = TRUE), drot_dyX = IQR(rot_dyX, na.rm = TRUE),
    mrot_dyY = median(rot_dyY, na.rm = TRUE), drot_dyY = IQR(rot_dyY, na.rm = TRUE),
    mrot_dyZ = median(rot_dyZ, na.rm = TRUE), drot_dyZ = IQR(rot_dyZ, na.rm = TRUE),
    
    # Means
    meanrot_stx = mean(rot_stX, na.rm = TRUE), meanrot_sty = mean(rot_stY, na.rm = TRUE), meanrot_stz = mean(rot_stZ, na.rm = TRUE),
    meanrot_dyx = mean(rot_dyX, na.rm = TRUE), meanrot_dyy = mean(rot_dyY, na.rm = TRUE), meanrot_dyz = mean(rot_dyZ, na.rm = TRUE),
    
    # Mins and Maxs
    minrot_stx = min(rot_stX, na.rm = TRUE), minrot_sty = min(rot_stY, na.rm = TRUE), minrot_stz = min(rot_stZ, na.rm = TRUE),
    maxrot_stx = max(rot_stX, na.rm = TRUE), maxrot_sty = max(rot_stY, na.rm = TRUE), maxrot_stz = max(rot_stZ, na.rm = TRUE),
    
    minrot_dyx = min(rot_dyX, na.rm = TRUE), minrot_dyy = min(rot_dyY, na.rm = TRUE), minrot_dyz = min(rot_dyZ, na.rm = TRUE),
    maxrot_dyx = max(rot_dyX, na.rm = TRUE), maxrot_dyy = max(rot_dyY, na.rm = TRUE), maxrot_dyz = max(rot_dyZ, na.rm = TRUE),
    
    # Standard Deviations
    sdrot_stx = sd(rot_stX, na.rm = TRUE), sdrot_sty = sd(rot_stY, na.rm = TRUE), sdrot_stz = sd(rot_stZ, na.rm = TRUE),
    sdrot_dyx = sd(rot_dyX, na.rm = TRUE), sdrot_dyy = sd(rot_dyY, na.rm = TRUE), sdrot_dyz = sd(rot_dyZ, na.rm = TRUE),
    
    # Roll and Pitch
    mroll = median(roll, na.rm = TRUE), droll = IQR(roll, na.rm = TRUE),
    mpitch = median(pitch, na.rm = TRUE), dpitch = IQR(pitch, na.rm = TRUE),
    
    # Metadata and Behavioral Modifiers are pulled right here:
    Modifiers = get_mode(Modifiers),
    Timestamp = format(min(Timestamp)),
    ID = as.character(first(TagID)),
    
    .groups = "drop"
  ) %>%
  select(-window_id)

## Training HMM
library(hmmTMB)

# hmmTMB requires the data to be a dataframe with an 'ID' column (to separate individuals) and the data streams.
data_hmm <- data_2s %>%
  dplyr::select(ID, Modifiers, meanrot_stx, sdrot_dyz, drot_Y, 
                droll, maxrot_dyx, minrot_sty, sdrot_stx, mrot_Y) %>%
  mutate(ID = as.factor(ID)) %>%
  na.omit() %>%
  filter(Modifiers != "Missing data") %>%
  mutate(
    Modifiers = case_when(
      Modifiers == "Agnostic behavior" ~ "Other",
      Modifiers == "Standing" ~ "Resting",
      .default = as.character(Modifiers)
    )
  )

# Define Data Streams and Distributions
dists <- list(
  meanrot_stx = "norm",
  sdrot_dyz   = "gamma", # SD is strictly positive, gamma is better than norm
  drot_Y      = "gamma", # IQR is strictly positive
  droll       = "gamma", 
  maxrot_dyx  = "norm",
  minrot_sty  = "norm",
  sdrot_stx   = "gamma",
  mrot_Y      = "norm"
)

# Initialize Parameters using your Known States (Modifiers)
# Because you have a labeled training dataset, you can calculate the exact starting parameters for the model 
# to ensure it converges on your specific behavioral definitions, rather than relying on K-means clustering.
n_states <- length(unique(data_hmm$Modifiers))
state_names <- sort(unique(data_hmm$Modifiers))

# Extract empirical means and standard deviations from your labeled data
initial_pars <- lapply(names(dists), function(stream) {
  
  # Calculate empirical means (m) and standard deviations (s) per state
  m <- tapply(data_hmm[[stream]], data_hmm$Modifiers, mean, na.rm = TRUE)
  s <- tapply(data_hmm[[stream]], data_hmm$Modifiers, sd, na.rm = TRUE)
  
  # Return the properly named numeric vectors based on the distribution
  if (dists[[stream]] == "norm") {
    list(mean = as.numeric(m), sd = as.numeric(s))
  } else if (dists[[stream]] == "gamma") {
    list(shape = as.numeric(m^2 / s^2), scale = as.numeric(s^2 / m))
  }
})

# Assign the stream names (meanrot_stx, sdrot_dyz, etc.) to the list
names(initial_pars) <- names(dists)

# Cosntruct observation model (ties data, distributions, and initial parameters together)
obs <- Observation$new(
  data = data_hmm, 
  dists = dists, 
  par = initial_pars,
  n_states = n_states
)



