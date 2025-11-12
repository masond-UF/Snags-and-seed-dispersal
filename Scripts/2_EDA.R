## --------------- HEADER ------------------------------------------------------
## Script name: 2_EDA.R
## Author: David S. Mason
## Affiliation: The Jones Center at Ichauway
## Date Created: 2025-11-11
## Date Last Modified: 2025-11-11
## Copyright (c) David S. Mason, 2025
## Contact: david.mason@jonesctr.org
## Purpose of script: This script explores the seed rain beneath snags

## --------------- SET—UP WORKSPACE --------------------------------------------
library(tidyverse)
library(lubridate)
library(tidylog)
library(styler)

# Clear the decks
rm(list=ls())

# Bring in the data
seed.mat <- read.csv("Clean-data/matrix.csv")

## --------------- PREPARE DATA ------------------------------------------------

# Summarize seeds at traps
seed.mat$Unknown <- apply(seed.mat[,34:49], sum, na.rm = TRUE,
												 MARGIN = 1)

# Drop individual unknow columns 
seed.mat <- seed.mat |>
	dplyr::select(1:33,50)

# Pivot longer
seeds.lg <- seed.mat |>
	pivot_longer(cols = 20:34, names_to = 'Taxa', values_to = 'Seeds')
str(seeds.lg)

## --------------- TAXA AT SNAGS VS CANOPY -------------------------------------

# summarize total seeds per taxa × trap category
seed.summary <- seeds.lg %>%
  group_by(Taxa, Trap.category) %>%
  summarise(total_seeds = sum(Seeds, na.rm = TRUE), .groups = "drop")

ggplot() +
  geom_jitter(
    data = seeds.lg,
    aes(x = Taxa, y = Seeds, color = Trap.category),
    position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8),
    alpha = 0.5,
    size = 2
  ) +
  geom_col(
    data = seed.summary,
    aes(x = Taxa, y = total_seeds, fill = Trap.category),
    position = position_dodge(width = 0.8),
    alpha = 0.4
  ) +
  theme_bw() +
  labs(
    title = "Seed totals per Taxa",
    y = "Total seeds",
    x = "Taxa"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## --------------- FACETTED BY SITE --------------------------------------------

seed.summary.site <- seeds.lg %>%
  group_by(Site, Taxa, Trap.category) %>%
  summarise(total_seeds = sum(Seeds, na.rm = TRUE), .groups = "drop")

ggplot() +
  geom_jitter(
    data = seeds.lg,
    aes(x = Taxa, y = Seeds, color = Trap.category),
    position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8),
    alpha = 0.5,
    size = 1.5
  ) +
  geom_col(
    data = seed.summary.site,
    aes(x = Taxa, y = total_seeds, fill = Trap.category),
    position = position_dodge(width = 0.8),
    alpha = 0.4
  ) +
  facet_wrap(~ Site, scales = 'free_y') +
  theme_bw() +
  labs(
    title = "Seed totals per Taxa by SITE",
    y = "Total seeds",
    x = "Taxa"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggplot() +
  geom_jitter(
    data = seeds.lg,
    aes(x = Taxa, y = Seeds, color = Trap.category),
    position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8),
    alpha = 0.5,
    size = 1.5
  ) +
  geom_col(
    data = seed.summary.site,
    aes(x = Taxa, y = total_seeds, fill = Trap.category),
    position = position_dodge(width = 0.8),
    alpha = 0.4
  ) +
  facet_wrap(~ Site) +
  theme_bw() +
  labs(
    title = "Seed totals per Taxa by SITE",
    y = "Total seeds",
    x = "Taxa"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
## --------------- FACETTED BY SITE --------------------------------------------

library(tidytext)

trap.totals <- seeds.lg %>%
  group_by(Site, Trap.ID..Snag.ID.No.snag..., Trap.category) %>%
  summarise(total_seeds = sum(Seeds, na.rm = TRUE), .groups = "drop") %>%
  # reverse order: negative total_seeds makes largest value come first
  mutate(trap_ordered = reorder_within(Trap.ID..Snag.ID.No.snag..., -total_seeds, Site))

ggplot(trap.totals, aes(x = trap_ordered, y = total_seeds, fill = Trap.category)) +
  geom_col(alpha = 0.9, width = 0.8) +
  scale_x_reordered() +
  facet_wrap(~ Site, scales = "free_x", nrow = 2) +
  theme_bw() +
  labs(
    title = "Total seeds per trap",
    x = "Trap ID",
    y = "Total seeds",
    fill = "Trap category"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid.major.x = element_blank()
  )

ggplot(trap.totals, aes(x = trap_ordered, y = total_seeds, fill = Trap.category)) +
  geom_col(alpha = 0.9, width = 0.8) +
  scale_x_reordered() +
  theme_bw() +
  labs(
    title = "Total seeds per trap",
    x = "Trap ID",
    y = "Total seeds",
    fill = "Trap category"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    panel.grid.major.x = element_blank()
  )

## --------------- CUMULATIVE TOTALS -------------------------------------------

# Convert Date to date type
seeds.lg <- seeds.lg %>%
  mutate(Date = as.Date(Date))

# cumulative seeds per site through time
cumulative.seeds <- seeds.lg %>%
  group_by(Site, Date) %>%
  summarise(daily_seeds = sum(Seeds, na.rm = TRUE), .groups = "drop") %>%
  group_by(Site) %>%
  arrange(Date) %>%
  mutate(cumulative_seeds = cumsum(daily_seeds))

ggplot(cumulative.seeds, aes(x = Date, y = cumulative_seeds, color = Site,
														 group = Site)) +
  geom_line(size = 1.2) +
  geom_point() +
  theme_bw() +
  labs(
    title = "Cumulative total seeds collected over time",
    y = "Cumulative seeds",
    x = "Date"
  ) +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

cumulative.taxa <- seeds.lg %>%
  group_by(Site, Date, Taxa) %>%
  summarise(present = any(Seeds > 0), .groups = "drop") %>%
  group_by(Site, Date) %>%
  summarise(unique_taxa = sum(present), .groups = "drop") %>%
  group_by(Site) %>%
  arrange(Date) %>%
  mutate(cumulative_taxa = cumsum(unique_taxa))

ggplot(cumulative.taxa, aes(x = Date, y = cumulative_taxa, color = Site,
														group = Site)) +
  geom_line(size = 1.2) +
  geom_point() +
  theme_bw() +
  labs(
    title = "Cumulative number of taxa observed over time",
    y = "Cumulative taxa observed",
    x = "Date"
  ) +
	theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


## --------------- nMDS --------------------------------------------------------

library(vegan)

# Sum all seed species by trap (aggregating across visits)
seed.mat <- seed.mat |>
	dplyr::select(-Unknown) |>
  group_by(Site, Trap.ID..Snag.ID.No.snag..., Trap.category) |>
  summarise(across(NYSY:Euphorbiaceae, sum, na.rm = TRUE)) |>
  ungroup()

# remove empty species
seed.mat <- seed.mat |>
  filter(rowSums(across(where(is.numeric))) > 0)

spec <- seed.mat[,4:17]
pred <- seed.mat[,1:3]

# run NMDS (on taxa only)
nmds <- metaMDS(spec, k = 2, trymax = 100)

# extract scores for plotting
nmds.points <- as.data.frame(scores(nmds, display = "sites"))
nmds.species <- as.data.frame(scores(nmds, display = "species"))

plot.data <- bind_cols(nmds.points, pred)

# plot NMDS ordination
library(ggrepel)
ggplot(data = plot.data, aes(x = NMDS1, y = NMDS2)) +
  # site points
  geom_point(aes(colour = Site, shape = Trap.category), size = 3) +
  # species labels (using their own data frame)
  geom_text_repel(data = nmds.species,
                  aes(x = NMDS1, y = NMDS2, label = rownames(nmds.species)),
                  colour = "black", size = 3) +
  coord_equal() +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank(),
        legend.position = "right") +
  labs(x = "NMDS1", y = "NMDS2",
       colour = "Trap category",
       shape = "Site")
