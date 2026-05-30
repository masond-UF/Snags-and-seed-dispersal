## --------------- HEADER ------------------------------------------------------
## Script name: 3_Analysis.R
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
seed.mat <- read.csv("Clean-data/matrix.csv") |>
	dplyr::select(-SOBI)

# Make the column name easier
colnames(seed.mat)[4] <- "ID"
colnames(seed.mat)[12] <- "Lowest.limb"

# Number of species observed 
sum(colnames(seed.mat) %in% colnames(seed.mat)[
   which(colnames(seed.mat) == "NYSY"):
   which(colnames(seed.mat) == "sp16")
])

# Number of samples
seed.mat |>
	group_by(Site, ID) |>
	tally() |> 
	pull(n) |>
  range()

# Add site classification
nat <- c('Arcadia', 'Greenwood')

seed.mat <- seed.mat |>
	mutate(Understory = ifelse(Site %in% nat, "Natural", "Disturbed"))

rm(nat)
## --------------- TOTAL SEEDS -------------------------------------------------

# Create a dataframe for total seeds
seed.mat.tot <- seed.mat

# Tally seeds at observations
seed.mat.tot$Total.seeds <- rowSums(seed.mat.tot[,20:48], na.rm = TRUE)

# Get total seeds at traps across all observations
seed.mat.tot <- seed.mat.tot |>
	group_by(Site, Understory, ID, Trap.category, Lowest.limb) |>
	summarize(Total.seeds = sum(Total.seeds))

# Total seeds collected
sum(seed.mat.tot$Total.seeds)

seed.mat.tot$Arrival <- ifelse(seed.mat.tot$Total.seeds > 0, 1, 0)

# Model
library(lme4)
tot.mod <- glmer.nb(Total.seeds ~ Understory + Trap.category + (1|Site), 
							 data = seed.mat.tot)
summary(tot.mod)	

library(performance)
check_model(tot.mod, check = c("pp_check", "homogeneity", "outliers", "vif", "overdispersion"),
						residual_type = 'normal')

library(DHARMa)
tot.mod.sim <- simulateResiduals(tot.mod)
plot(tot.mod.sim)

library(emmeans)
emmeans(tot.mod, ~Trap.category, type = 'response') |>
	  pairs(reverse = TRUE)

# Create
# library(multcomp)     
# library(multcompView)
emm.summary <- emmeans(tot.mod, ~Trap.category, type = 'response') |> 
  as.data.frame() |> 
  rename(emmean = response) # Rename for clarity

seed.mat.tot <- seed.mat.tot |> 
  left_join(emm.summary, by = "Trap.category")

tot.plot <- ggplot(seed.mat.tot, aes(x = Trap.category)) +
  geom_jitter(aes(y = Total.seeds, color = Site), 
              width = 0.2, 
              height = 0, 
              alpha = 0.5, 
              size = 2.5) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                width = 0.1,
                size = 0.4) +
  geom_point(aes(y = emmean),
             size = 4,
             color = "black") +
  theme_bw() +
  labs(x = "Trap Category", y = "Predicted Seed Count")+
	annotate("text", x = 0.6, y = Inf, label = "p < 0.001", vjust = 2, hjust = 0)

rm(tot.mod, tot.mod.sim)

## --------------- ARRIVAL -----------------------------------------------------

arrival.mod <- glm(Arrival ~ Understory + Site + Trap.category, 
                   data = seed.mat.tot,
                   family = binomial)

summary(arrival.mod)	

library(performance)
check_model(arrival.mod, 
            check = c("pp_check", "homogeneity", "outliers", "vif"),
            residual_type = 'normal')
library(car)
Anova(arrival.mod)

library(DHARMa)
arrival.mod.sim <- simulateResiduals(arrival.mod)
plot(arrival.mod.sim)

library(emmeans)
emmeans(arrival.mod, ~Trap.category, type = 'response') |>
	  pairs(reverse = TRUE)
emmeans(arrival.mod, ~Site, type = 'response') |>
	  pairs(reverse = TRUE)
emmeans(arrival.mod, ~Understory, type = 'response') |>
	  pairs(reverse = TRUE)

# Create
emm.summary <- emmeans(arrival.mod, ~Trap.category, type = 'response') |> 
  as.data.frame() |> 
  rename(emmean = prob)  # Rename predicted probability for clarity

seed.mat.tot <- seed.mat.tot |> 
  left_join(emm.summary, by = "Trap.category")

arrival.plot <- ggplot(seed.mat.tot, aes(x = Trap.category)) +
  geom_jitter(aes(y = Arrival, color = Site), 
              width = 0.2, 
              height = 0.02, 
              alpha = 0.5, 
              size = 2.5) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                width = 0.1,
                size = 0.4) +
  geom_point(aes(y = emmean),
             size = 4,
             color = "black") +
  theme_bw() +
  labs(x = "Trap Category", y = "Predicted Probability of Seed Arrival") +
	annotate("text", x = 0.6, y = 1.05, label = "p < 0.001", vjust = 1)

rm(seed.mat.tot, arrival.mod, arrival.mod.sim)

## --------------- RICHNESS ----------------------------------------------------

# Create a dataframe for seed richness
seed.mat.rich <- seed.mat

# Getting total seeds for each species
seed.mat.rich <- seed.mat.rich |>
  group_by(ID, Site, Trap.category) |>
	summarise(across(NYSY:sp16, \(x) sum(x, na.rm = TRUE)), .groups = "drop")

# Convert to binary
seed.mat.rich <- seed.mat.rich |>
	mutate(across(4:32, ~ ifelse(.x > 0, 1, 0)))

# Tally species to calculate richness
seed.mat.rich$Richness <- rowSums(seed.mat.rich[,4:32], na.rm = TRUE)

# Model
rich.mod <- glm(Richness ~ Site + Trap.category, 
							 data = seed.mat.rich,
							 family = 'poisson')
summary(rich.mod)
Anova(rich.mod)

# Explained variance
total.dev <- rich.mod$null.deviance
trap.dev.exp <- 16.7129
site.dev.exp <- 3.2434

(trap.dev.exp / total.dev) * 100
(site.dev.exp / total.dev) * 100
(rich.mod$deviance / total.dev) * 100

check_model(rich.mod, check = c("pp_check", "homogeneity", "outliers", "vif", "overdispersion"),
						residual_type = 'normal')

rich.mod.sim <- simulateResiduals(rich.mod)
plot(rich.mod.sim)

emmeans(rich.mod, ~Trap.category, type = 'response') |> 
		 pairs(reverse = TRUE)

emm.summary <- emmeans(rich.mod, ~Trap.category, type = 'response') |> 
  as.data.frame() |> 
  rename(emmean = rate) # Rename for clarity

seed.mat.rich <- seed.mat.rich |> 
  left_join(emm.summary, by = "Trap.category")

rich.plot <- ggplot(seed.mat.rich, aes(x = Trap.category)) +
  geom_jitter(aes(y = Richness, color = Site), 
              width = 0.2, 
              height = 0, 
              alpha = 0.5, 
              size = 2.5) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL),
                width = 0.1,
                size = 0.4) +
  geom_point(aes(y = emmean),
             size = 4,
             color = "black") +
  theme_bw() +
  labs(x = "Trap Category", y = "Seed Rain Richness")+
	annotate("text", x = 0.6, y = Inf, label = "p < 0.001", vjust = 2, hjust = 0)

ggplot(emm.summary, aes(x = Trap.category, y = emmean, fill = Trap.category)) +
  geom_col(color = "black", width = 0.7) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), 
                width = 0.2, 
                size = 1) +
  scale_fill_manual(values = c("Canopy gap" = "#F8766D", "Snag" = "#00BFC4")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none",
    axis.text = element_text(size = 12, color = "black")
  ) +
  labs(x = "Trap category", y = "Predicted seed count")

rm(rich.mod, rich.mod.sim, seed.mat.rich)

## --------------- DIVERSITY ---------------------------------------------------

# Create a dataframe for seed richness
seed.mat.div <- seed.mat

# Getting total seeds for each species
seed.mat.div <- seed.mat.div |>
  group_by(ID, Site, Trap.category) |>
	summarise(across(NYSY:sp16, \(x) sum(x, na.rm = TRUE)), .groups = "drop")

# Calculate diversity
library(vegan)
seed.mat.div <- seed.mat.div |> 
  mutate(Shannon = diversity(dplyr::select(pick(everything()), 4:32), index = "shannon"))

library(fitdistrplus)
descdist(seed.mat.div$Shannon, discrete = FALSE)

div.mod <- lm(Shannon ~ Site + Trap.category, 
							 data = seed.mat.div)
summary(div.mod)

check_model(div.mod, check = c("pp_check", "homogeneity", "outliers", "vif", "overdispersion"),
						residual_type = 'normal')

div.mod.sim <- simulateResiduals(div.mod)
plot(div.mod.sim)

emmeans(div.mod, ~Trap.category, type = 'response') |> pairs(reverse = TRUE)

emm.summary <- emmeans(div.mod, ~Trap.category, type = 'response') |> 
  as.data.frame()

seed.mat.div <- seed.mat.div |> 
  left_join(emm.summary, by = "Trap.category")

div.plot <- ggplot(seed.mat.div, aes(x = Trap.category)) +
  geom_jitter(aes(y = Shannon, color = Site), 
              width = 0.2, 
              height = 0, 
              alpha = 0.5, 
              size = 2.5) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),
                width = 0.1,
                size = 0.4) +
  geom_point(aes(y = emmean),
             size = 4,
             color = "black") +
  theme_bw() +
  labs(x = "Trap Category", y = "Diversity")+
	annotate("text", x = 0.6, y = Inf, label = "p = 0.012", vjust = 2, hjust = 0)


rm(div.mod, div.mod.sim, seed.mat.div, emm.summary)

## --------------- COMPOSITION -------------------------------------------------

# Getting total seeds for each species
seed.mat <- seed.mat |>
  group_by(ID, Site, Understory, Trap.category) |>
	summarise(across(NYSY:sp16, \(x) sum(x, na.rm = TRUE)), .groups = "drop")

# Separate abundance data
spp.data <- seed.mat |> dplyr::select(NYSY:last_col())
env.data <- seed.mat |> dplyr::select(ID, Site, Understory, Trap.category)

# Keep species present in at least 3 samples
common.spp <- spp.data |>
  dplyr::select(where(\(x) sum(x > 0) >= 3))

# Calculate number of seeds from common species
sum(common.spp)
310/422

# Model
library(mvabund)
seed.mv <- mvabund(common.spp)
meanvar.plot(seed.mv)

comm.mod <- manyglm(seed.mv ~ Site + Understory + Trap.category, family = "negative.binomial", data = env.data)

plot(comm.mod)
anova(comm.mod, p.uni = "adjusted")

plot(seed.mv ~ as.factor(env.data$Trap.category), cex.axis = 0.8, cex = 0.8) 
plot(seed.mv ~ as.factor(env.data$Site), cex.axis = 0.8, cex = 0.8) 

rm(comm.mod, seed.mv)
## --------------- MAKE A FIGURE -----------------------------------------------

# Make long for better figure
common.spp <- cbind(env.data, common.spp)
common.spp.lg <- common.spp |>
	pivot_longer(cols = 5:11, names_to = "Taxa", values_to = "Seeds")

seed.summary.site <- common.spp.lg %>%
  group_by(Site, Taxa, Trap.category) %>%
  summarise(mean_seeds = mean(Seeds, na.rm = TRUE), .groups = "drop")

# Fix species names
common.spp.lg <- common.spp.lg |>
  mutate(
    Taxa = dplyr::recode(
      Taxa,
      "QUERCUS.broken" = "Quercus spp.",
      "Pinus" = "Pinus spp."
    )
  )
seed.summary.site <- seed.summary.site |>
  mutate(
    Taxa = case_when(
      Taxa == "QUERCUS.broken" ~ "Quercus spp.",
      Taxa == "Pinus" ~ "Pinus spp.",
      TRUE ~ Taxa
    )
  )
common.spp.lg$Taxa <- factor(common.spp.lg$Taxa)
seed.summary.site$Taxa <- factor(seed.summary.site$Taxa)

ggplot() +
  geom_jitter(
    data = common.spp.lg,
    aes(x = Taxa, y = Seeds, color = Trap.category),
    position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.8),
    alpha = 0.5,
    size = 0.75
  ) +
  geom_col(
    data = seed.summary.site,
    aes(x = Taxa, y = mean_seeds, fill = Trap.category),
    position = position_dodge(width = 0.8),
    alpha = 0.4
  ) +
  facet_wrap(~ Site) +
  scale_y_continuous(
  trans = "log1p",
  breaks = c(0, 1, 5, 10, 25, 50, 100)
	) +
  theme_bw() +
	theme(
  strip.background = element_blank(),
  strip.text = element_text(face = "bold", size = 11),
  aspect.ratio = 0.6,
	) +
  labs(
    y = "Predicted Total Seed Count (log[y + 1])",
    x = "Taxa",
    fill = "Trap Category",
    color = "Trap Category"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
  			axis.title.y = element_text(size = 9, face = "bold"),
  			axis.title.x = element_text(face = "bold")
  			)

rm(common.spp.lg, seed.summary.site)

## --------------- PERMANOVA ---------------------------------------------------

library(vegan)

# Drop empty rows
keep <- rowSums(spp.data) > 0
spp.data <- spp.data[keep, ]
env.data <- env.data[keep, ]

# Create a distance matrix from the species data
dist.bc <- vegdist(sqrt(spp.data), method = "bray")

# No differences in dispersion
bd <- betadisper(dist.bc, env.data$Trap.category)
anova(bd)

bd <- betadisper(dist.bc, env.data$Site)
anova(bd)

# Test significance
adonis2(
  dist.bc ~ Trap.category + Site,
  data = env.data,
  by = "margin"
)

# Percentage of animal-dispersed and wind-dispersed seeds at snags vs. gaps
# vector of animal-dispersed columns
animal.cols <- c("NYSY", "PAQU", "MAGR", "CAAE", "PHAM", "QUERCUS.broken")

results <- common.spp |>
  mutate(animal.total = rowSums(across(all_of(animal.cols))),
         pinus.total = Pinus) |>
  summarise(
    pct.animal.snag =
      sum(animal.total[Trap.category == "Snag"]) /
      sum(animal.total) * 100,
    pct.pinus.gap =
      sum(pinus.total[Trap.category == "Canopy gap"]) /
      sum(pinus.total) * 100
  )

results

