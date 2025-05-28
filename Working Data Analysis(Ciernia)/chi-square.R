library(lme4)
library(car)
library(dplyr)

cluster_data <- read.csv("/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/Data Files/cluster-percentage-data.csv")

# Optional: scale to proportion if necessary
if (max(cluster_data$percentage, na.rm = TRUE) > 1) {
  cluster_data$percentage <- cluster_data$percentage / 100
}

brain_regions <- unique(cluster_data$BrainRegion)
wald_results <- list()

for (region in brain_regions) {
  cat("\nAnalyzing region:", region, "\n")
  
  region_data <- cluster_data %>% filter(BrainRegion == region)
  
  # Linear Mixed Model instead of GLMM
  model <- lmer(percentage ~ Cluster * Treatment + (1 | MouseID),
                data = region_data %>% filter(Sex=="M"))
  
  # Type II Wald test
  wald_test <- Anova(model, type = 2)
  print(wald_test)
  
  wald_results[[region]] <- wald_test
}


library(ggplot2)
library(dplyr)

# Filter to PVN only
pvn_data <- cluster_data %>% filter(BrainRegion == "HYPO"&Sex=="M")

# Summarize mean ± SEM per Cluster × Treatment
summary_data <- pvn_data %>%
  group_by(Cluster, Treatment) %>%
  summarise(
    mean_percent = mean(percentage),
    sem = sd(percentage) / sqrt(n())
  )

pvn_data$Treatment <- factor(pvn_data$Treatment, levels = c("CON", "COLD"))
summary_data$Treatment <- factor(summary_data$Treatment, levels = c("CON", "COLD"))


# Plot
ggplot(summary_data, aes(x = Cluster, y = mean_percent, fill = Treatment)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(aes(ymin = mean_percent - sem, ymax = mean_percent + sem),
                position = position_dodge(width = 0.8), width = 0.2) +
  labs(title = "PVN: Cluster Proportions by Treatment",
       x = "Cluster", y = "Mean Percentage of Cells") +
  scale_fill_manual(values = c("CON" = "steelblue", "COLD" = "tomato")) +
  theme_minimal(base_size = 14)
