library(factoextra)
library(ppclust)
library(dplyr)
library(ggsignif)
library(ggrepel)
library(performance)  # for model checks
library(lme4)
library(lmerTest) # For p-values
library(MASS)
library(readxl)
library(emmeans)

morph_data <- read.csv("/Users/alexlawson/Masters-Data-Final/Morphoglia/PVN_all/Morphoglia/MorphoGlia/Data/Morphology_HDBSCAN_30_0.1_150_5.csv", stringsAsFactors = FALSE)
View(morph_data)

# Assuming your dataframe is called df

# Count number of unique groups
num_unique_groups <- morph_data %>% 
  distinct(group) %>% 
  nrow()
cat("Number of unique groups:", num_unique_groups, "\n")

# Count number of unique regions within each group
regions_per_group <- morph_data %>%
  group_by(group) %>%
  summarise(unique_regions = n_distinct(region))

View(regions_per_group)

color_map <- c(
  "0" = "184,138,255",   # Soft pink (bubblegum pink)
  "1" = "240,200,200",   # Lavender (pale purple)
  "2" = "231,221,119",   # Pastel teal (light aqua-blue)
  "3" = "59,235,255",    # Lemon yellow (bright pastel yellow)
  "4" = "119,221,119"    # Pastel green (minty green)
)

# Add RGB color column based on 'clusters'
morph_data <- morph_data %>%
  mutate(Cluster_Color = color_map[as.character(Clusters)])

morph_data  <- morph_data [morph_data $Clusters != -1,]



write_csv(morph_data, "/Users/alexlawson/Masters-Data-Final/Morphoglia/PVN_all/Morphoglia/MorphoGlia/Data/Morphology_HDBSCAN_30_0.1_150_5-v2.csv")

data <- read.csv("/Users/alexlawson/Masters-Data-Final/Morphoglia/PVN_all/Morphoglia/MorphoGlia/Data/Morphology_HDBSCAN_15_0.05_150_5-v1.csv")

df = data %>% filter(tissue=="M")

# Step 1: Compute average cell area per mouse
mouse_avg_area <- df %>%
  group_by(group, treatment) %>%
  summarise(mean_area = mean(`Cell_eccentricity`, na.rm = TRUE),
            .groups = 'drop')

# Step 2: Plot the average cell area per mouse by treatment
ggplot(mouse_avg_area, aes(x = treatment, y = mean_area)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
  labs(title = "Average Microglia Cell Area per Mouse by Treatment",
       x = "Treatment Group",
       y = "Mean Cell Area") +
  theme_minimal()

# Step 3: Wilcoxon rank-sum test to compare treatment groups
wilcox_result <- wilcox.test(mean_area ~ treatment, data = mouse_avg_area)
print(wilcox_result)
