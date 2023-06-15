# Data processing and analysis of YBFMP Beach Seine data -----------------------
# 1/19/2023 TMF

library(tidyverse)
library(lubridate)
library(janitor)
library(vegan)

# Set working directory --------------------------------------------------------
setwd("C:/R/work-projects/")
setwd("./YBFMP-beach-seines")
getwd()

# Set output directory 
output <- "output"

# Clean workspace --------------------------------------------------------------
rm(list=ls())

# Read input data from JT ------------------------------------------------------
df_fish_old < read_csv("data/Beach Seine CPUE All_R_2.csv")
df_fish <- read_csv("data/Beach Seine CPUE All_R_2.csv")

# Clean up column names
df_fish <- df_fish %>% clean_names(case = "big_camel")

# Remove non-fish organisms from dataset ---------------------------------------
unfish <- c("Siberian Prawn",
            "NoCatch",
            "Mississippi Grass Shrimp")

df_fish <- df_fish %>%
  filter(!(CommonName %in% unfish))

# Remove samples from outside AL/BL region -------------------------------------
df_fish <- df_fish %>% filter(Region != "N/A")

# Create yearly averages for Above Lisbon (AL) and Below Lisbon (BL) -----------
df_fish_avg <- df_fish %>%
  group_by(WaterYear, Region, Family, CommonName) %>%
  summarize(Mean.CPUE = mean(Cpue)) %>%
  ungroup()
  
# Take 4th root transform of data ----------------------------------------------
df_fish_avg <- df_fish_avg %>%
  mutate(CPUE.4th = (Mean.CPUE)^(1/4))

# Select and pivot data for NMDS analysis --------------------------------------
# Subset to just cyprinids and centrarchids ------------------------------------
df_fish_CC <- df_fish_avg %>%
  filter(Family %in% c("Cyprinidae", "Centrarchidae"))

# Create NMDS plots ------------------------------------------------------------

df_fish_NMDS <- df_fish_CC %>%
  select(WaterYear, Region, CommonName, CPUE.4th)

genw <- pivot_wider(df_fish_NMDS, 
                    names_from = "CommonName", 
                    values_from = "CPUE.4th",
                    values_fill = 0)

# Calculate the nMDS using vegan 
# A good rule of thumb: stress < 0.05 provides an excellent representation in reduced dimensions,
# < 0.1 is great, < 0.2 is good/ok, and stress < 0.3 provides a poor representation.
df_fish_NMDS <- metaMDS(comm = genw[c(3:21)],
                        distance = "bray",
                        k = 3,
                        trymax = 200)

# Create Shepard plot which shows scatter around the regression between the
# interpoint distances in the final configuration (i.e., the distances between 
# each pair of communities) against their original dissimilarities.
stressplot(df_fish_NMDS)

# Use vegan's scores function to extract the site scores and convert to df
data.scores <- as_tibble(scores(df_fish_NMDS, display = "sites"))

# Combine metadata with NMDS data scores to plot in ggplot
meta <- genw %>% select(1:2)
meta <- cbind(meta, data.scores)

# Read in years as chr so it is not displayed as a gradient
meta$WaterYear <- as.character(meta$WaterYear)

# Plot NMDS values -------------------------------------------------------------
NMDS.Region.plot <- ggplot(meta, aes(x = NMDS1, y = NMDS2, color = Region)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.50) +
  labs(color = "Region") +
  theme_bw()

NMDS.Region.plot

ggsave(path = output,
       filename = "NMDS-all-fish.png", 
       device = "png",
       scale=1.0, 
       units="in",
       height=3,
       width=4, 
       dpi="print")

# Calculate PERMANOVA values for by-region comparison --------------------------
adon.results <- adonis2(genw[c(3:21)] ~ genw$Region, 
                        method = "bray", 
                        perm = 999)

print(adon.results)

anosim(genw[c(3:21)], genw$Region, permutations = 10000, distance = "bray")
