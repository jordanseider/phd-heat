###
# Identification of reindeer fine-scale foraging behaviour using tri-axial accelerometer data
# Heidi Rautiainen, with acknowledgements to Måns Karlsson for hidden Markov model
# 2020-08-26
# https://doi.org/10.5061/dryad.8sf7m0cs7
###

library(tidyverse)

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
