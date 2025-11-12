## --------------- HEADER ------------------------------------------------------
## Script name: 1_Clean.R
## Author: David S. Mason
## Affiliation: The Jones Center at Ichauway
## Date Created: 2025-11-11
## Date Last Modified: 2025-11-11
## Copyright (c) David S. Mason, 2025
## Contact: david.mason@jonesctr.org
## Purpose of script: This script cleans the seed rain beneath snags

## --------------- SET—UP WORKSPACE --------------------------------------------
library(tidyverse)
library(lubridate)
library(tidylog)
library(styler)

# Clear the decks
rm(list=ls())

# Bring in the data
seeds <- read.csv("Raw-data/seeds.csv")

# Add zeroes
seeds[1:101, 20:49][is.na(seeds[1:101, 20:49])] <- 0 

# Rename columns
rename <- function(df) {
  start_col <- 34
  end_col <- 49
  num_cols <- end_col - start_col + 1
	new_names <- paste0("sp", 1:num_cols)
  colnames(df)[start_col:end_col] <- new_names
   return(df)
}
seeds <- rename(seeds)

# Fix date
seeds$Date <- mdy(seeds$Date)

# Save matrix for community analysis
write.csv(seeds, "Clean-data/matrix.csv", row.names = FALSE)

# Summarize seeds at traps
seeds$Tot.seeds <- apply(seeds[,20:49], sum, na.rm = TRUE,
												 MARGIN = 1)

# Filter out individual species
seeds <- seeds |>
	dplyr::select(1:19, Tot.seeds)

# Save summarized data for single metric
write.csv(seeds, "Clean-data/summarized.csv", row.names = FALSE)
