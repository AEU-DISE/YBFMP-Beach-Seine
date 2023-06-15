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
df_fish <- read_csv("data/Bseine_Raw.csv")

# Clean up column names
df_fish <- df_fish %>% 
  clean_names(case = "big_camel") %>%
  rename("CPUE" = "Cpue")

# Remove zeros and NAs
df_fish <- df_fish %>% 
  filter(CPUE != 0)

# Change Inland Silverside to Mississippi Silverside
df_fish <- df_fish %>% 
  mutate(CommonName = str_replace(CommonName, 
                                  "Inland Silverside",
                                  "Mississippi Silverside"))

# Remove samples from before 2011
df_fish <- df_fish %>% 
  filter(WaterYear > 2010)

# Remove non-fish organisms from dataset ---------------------------------------

# Figure out which taxa aren't fish
sort(unique(df_fish$CommonName))

# Siberian Prawn and Mississippi Grass Shrimp are the only two
unfish <- c("Siberian Prawn", "Mississippi Grass Shrimp")

df_fish <- df_fish %>% filter(!(CommonName %in% unfish))

# Remove duplicate samples ------------------------------------------

# Find samples with more than one entry for the same taxon
dups <- df_fish %>%
  dplyr::group_by(SampleDate, WaterYear, WaterYearType, StationCode, Region, CommonName) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n > 1L) 

# write_csv(dups, file = "dups.csv)

# Filter out duplicate samples from AEU resample on 6-27-2011 due to gear issue
df_fish <- df_fish %>% 
  filter(!(SampleDate == "2011-06-27" & CommonName == "Bigscale Logperch" & StationCode == "BL5" & Count == 1)) %>% 
  filter(!(SampleDate == "2011-06-27" & CommonName == "Golden Shiner" & StationCode == "BL5" & Count == 1)) %>% 
  filter(!(SampleDate == "2011-06-27" & CommonName == "Mississippi Silverside" & StationCode == "BL5" & Count == 1)) %>% 
  filter(!(SampleDate == "2011-06-27" & CommonName == "Sacramento Pikeminnow" & StationCode == "BL5" & Count == 1)) %>% 
  filter(!(SampleDate == "2011-06-27" & CommonName == "Shimofuri Goby" & StationCode == "BL5" & Count == 1)) %>% 
  filter(!(SampleDate == "2011-06-27" & CommonName == "Splittail" & StationCode == "BL5" & Count == 6)) %>% 
  filter(!(SampleDate == "2011-06-27" & CommonName == "Striped Bass" & StationCode == "BL5" & Count == 1)) %>% 
  filter(!(SampleDate == "2011-06-27" & CommonName == "Tule Perch" & StationCode == "BL5" & Count == 1))

# Remove duplicate entry on 2/3/2016
df_fish <- df_fish %>% 
  filter(!(SampleDate == "2016-02-03" & CommonName == "Bluegill" & StationCode == "BL1" & Count == 1))

# Delete CHN recorded with count of 3 on 4/27/2017
df_fish <- df_fish %>% 
  filter(!(SampleDate == "2017-04-27" & CommonName == "Chinook Salmon" & StationCode == "BL5" & Count == 3))

# Remove other duplicate samples
df_fish <- distinct(df_fish)

# Create yearly averages for Above Lisbon (AL) and Below Lisbon (BL) -----------

# Add zeros when certain fish weren't caught
df_fish_w <- df_fish %>% 
  pivot_wider(names_from = "CommonName", 
              values_from = "CPUE",
              values_fill = 0) %>% 
  ungroup()

df_fish_avg <- df_fish %>%
  group_by(WaterYear, Region, Family, CommonName) %>%
  summarize(Mean.CPUE = mean(CPUE)) %>%
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
df_fish_NMDS <- metaMDS(comm = genw[c(3:20)],
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
NMDS.Region.plot <- ggplot(meta, aes(x = NMDS1, 
                                     y = NMDS2, 
                                     color = Region)) +
  geom_point(size = 3) +
#  geom_text(nudge_x = 0.04,
#            color = "black") +
  stat_ellipse(level = 0.50) +
  labs(color = "Region") +
  theme_bw()

NMDS.Region.plot +
  labs(x = NULL,
       y = NULL,
       title = "NMDS - Centrarchids and Cyprinids - 2011-2019") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank())

ggsave(path = output,
       filename = "NMDS-CC-only.png", 
       device = "png",
       scale=1.0, 
       units="in",
       height=3,
       width=4, 
       dpi="print")

# Calculate PERMANOVA values for by-region comparison --------------------------
adon.results <- adonis2(genw[c(3:20)] ~ genw$Region, 
                        method = "bray", 
                        perm = 999)

print(adon.results)

anosim(genw[c(3:20)], genw$Region, permutations = 10000, distance = "bray")
