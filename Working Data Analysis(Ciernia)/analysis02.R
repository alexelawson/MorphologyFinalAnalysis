#LOADING LIBRARIES AND SETTING WD
#clean
library(MicrogliaMorphologyR)
library(factoextra)
library(ppclust)
library(dplyr)
library(ggsignif)
library(ggrepel)
set.seed(1234)
setwd("~/Documents/GitHub/MorphologyAnalysis01")

#reading in the data 
raw_data01 <- read.csv("/Users/alexlawson/Desktop/GliaData/workingdata-final-05.csv", stringsAsFactors = FALSE)

raw_data_section <- read.csv("/Users/alexlawson/Desktop/Masters-Work/data-frames/cluster-percentage-data-w-sections.csv", stringsAsFactors = FALSE)
raw_data01_log <- transform_log(raw_data01,1, start=8, end=34) 
pca_data01_log <- pcadata(raw_data01_log, featurestart=8, featureend=34,
                          pc.start=1, pc.end=10)
pca_data01_log_transformed <- transform_scale(pca_data01_log, start=1, end=3) # scale pca data as input for k-means clustering
kmeans_input01_log <- pca_data01_log_transformed[1:3]
data_kmeans03_log <- kmeans(kmeans_input01_log, centers=3)
pca_kmeans03_log <- cbind(pca_data01_log[1:3], raw_data01, as.data.frame(data_kmeans03_log$cluster)) %>%
  rename(Cluster=`data_kmeans03_log$cluster`) 
clusterfeatures(pca_kmeans03_log, featurestart=11, featureend=37)
#1 - Ameboids
#2 - Rod-Like
#3 - Ramified

write.csv(pca_data01_log, "/Users/alexlawson/Desktop/Masters-Work/statistical-results/pca-data.csv", row.names = FALSE)


plot_log_3 <- clusterplots(pca_kmeans03_log, "PC1", "PC2")
cp_log_3 <- clusterpercentage(pca_kmeans03_log, "Cluster", Treatment, Sex, BrainRegion, MouseID)
plot <- clusterplots(pca_kmeans03_log, "PC1", "PC2")
cp_log_3 <- cp_log_3 %>% mutate(Cluster = 
                                  case_when(Cluster=="1" ~ "Ameboid",
                                            Cluster=="2" ~ "Rod-Like",
                                            Cluster=="3" ~ "Ramified"))




cp_log_3_withslides <- clusterpercentage(pca_kmeans03_log, "Cluster", Treatment, Sex, BrainRegion, MouseID, SlideNumber)
cp_log_3_withslides <- cp_log_3_withslides %>%
  mutate(SlideNumber = case_when(
    SlideNumber == "s07-1" ~ "s07",   
    TRUE ~ SlideNumber              
  ))
cp_log_3_withslides <- cp_log_3_withslides %>% mutate(Cluster = 
                                                        case_when(Cluster=="1" ~ "Ameboid",
                                                                  Cluster=="2" ~ "Rod-Like",
                                                                  Cluster=="3" ~ "Ramified"))
write.csv(cp_log_3_withslides, "/Users/alexlawson/Desktop/Masters-Work/statistical-results/cluster-percentage-data-w-sections.csv", row.names = FALSE)

View(cp_log_3)

# ---- PLOTS ----     
#cold vs. control cluster percentages: split by brainregion and sex

coldvcontrol_all <- cp_log_3 %>% 
  ggplot(aes(x = Cluster, y = percentage, group = interaction(Cluster, Treatment))) +
  facet_grid(~ BrainRegion) + 
  geom_boxplot(aes(group = interaction(Cluster, Treatment), fill = Treatment), outlier.shape = NA) +  # Remove outliers
  geom_point(aes(color = Sex, group = interaction(Cluster, Treatment)), 
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75), 
             size = 2, alpha = 0.5) +  # Points jittered and aligned with boxplots
  scale_fill_manual(values = c("white", "lightgray")) +
  scale_color_manual(values = c("M" = "blue", "F" = "pink")) +  # Colors for Sex
  ggtitle("Cold vs. Control Cluster Percentages Split by Brain Region") +
  labs(x = "Cluster", y = "Percentage", fill = "Treatment") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

coldvcontrol_all


coldvcontrol <- cp_log_3 %>% 
  ggplot(aes(x = Cluster, y = percentage, group = interaction(Cluster, Treatment))) +
  facet_grid(Sex ~ BrainRegion) + 
  geom_boxplot(aes(group = interaction(Cluster, Treatment), fill = Treatment), outlier.shape = NA) +  # Remove outliers
  scale_fill_manual(values = c("white", "lightgray")) +
  ggtitle("Cold vs. Control Cluster Percentages Split by Brain Region and Sex") +
  labs(x = "Cluster", y = "Percentage", fill = "Treatment") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))
coldvcontrol


#control m vs. female split by brain region, with dots 
cp_log_3 %>%
  filter(Treatment == "CON") %>%  # Filter for Control Males and Females
  ggplot(aes(x = Cluster, y = percentage, group = interaction(Cluster, Sex))) +
  facet_grid(. ~ BrainRegion) +  # Columns = Brain Region
  geom_boxplot(aes(fill = Sex), position = position_dodge(width = 0.8), outlier.shape = NA) +  # Compare sexes
  scale_fill_manual(values = c("pink", "lightblue")) +
  geom_point(aes(color = Sex),
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), 
             size = 1.5, alpha = 0.7) +  # Add jittered points for clarity
  ggtitle("Comparison of Control Males vs Control Females Across Brain Regions") +
  labs(x = "Cluster", y = "Percentage",fill = "Sex") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

#control m vs. female, with dots 
stats_input %>%
  filter(Treatment == "COLD") %>%  # Filter for Control Males and Females
  ggplot(aes(x = Cluster, y = percentage, group = interaction(Cluster, Sex))) +
  geom_boxplot(aes(fill = Sex), position = position_dodge(width = 0.8), outlier.shape = NA) +  # Compare sexes
  scale_fill_manual(values = c("pink", "lightblue")) +
  geom_point(aes(color = Sex),
             position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2), 
             size = 1.5, alpha = 0.7) +  # Add jittered points for clarity
  ggtitle("Comparison of Cold Stressed Males vs. Cold Stressed Females") +
  labs(x = "Cluster", y = "Percentage",fill = "Sex") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

#control male vs. female in the different brain regions, without dots
control_braincomparison <- cp_log_3 %>%
  filter(Treatment == "CON", Sex %in% c("M", "F")) %>%
  ggplot(aes(x = Cluster, y = percentage, group = interaction(Cluster, Sex))) +
  facet_wrap(~ BrainRegion, nrow = 1) +  # Adjust facets for better layout
  geom_boxplot(aes(fill = Sex), position = position_dodge(width = 0.8), outlier.shape = NA) +  # Remove outliers
  # Adjust point position
  scale_fill_manual(values = c("lightpink", "lightblue")) +
  ggtitle("Comparison of Clusters Across Brain Regions in Control Conditions") +
  labs(fill = "Sex") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

control_braincomparison

#cold male vs. female in the different brain regions, without dots
cp_log_3 %>%
  filter(Treatment == "COLD", Sex %in% c("M", "F")) %>%
  ggplot(aes(x = Cluster, y = percentage, group = interaction(Cluster, Sex))) +
  facet_wrap(~ BrainRegion, nrow = 1) +  # Adjust facets for better layout
  geom_boxplot(aes(fill = Sex), position = position_dodge(width = 0.8), outlier.shape = NA) +  # Remove outliers
  # Adjust point position
  scale_fill_manual(values = c("lightpink", "lightblue")) +
  ggtitle("Mouse Dataset: Cold Males vs Females Across Brain Regions") +
  labs(fill = "Sex") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

# First plot for SlideNumber s01-s06
raw_data_section %>%
  filter(BrainRegion == "HYPO", SlideNumber %in% c("s01", "s02", "s03", "s04", "s05", "s06")) %>%
  ggplot(aes(x = Cluster, y = percentage, group = interaction(Cluster, Treatment))) +
  facet_grid(Sex ~ SlideNumber) +
  geom_boxplot(aes(group = interaction(Cluster, Treatment), fill = Treatment), outlier.shape = NA) +
  scale_fill_manual(values = c("white", "lightgrey")) +
  ggtitle("Comparison of Cold vs. Control in the Hypothalamus, Separated by Sex and Section (s01-s06)") +
  labs(y = "Percentage", fill = "Treatment") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

# Second plot for SlideNumber s07-s12
cp_log_3_withslides %>%
  filter(BrainRegion == "HYPO", SlideNumber %in% c("s07", "s08", "s09", "s10", "s11", "s12")) %>%
  ggplot(aes(x = Cluster, y = percentage, group = interaction(Cluster, Treatment))) +
  facet_grid(Sex ~ SlideNumber) +
  geom_boxplot(aes(group = interaction(Cluster, Treatment), fill = Treatment), outlier.shape = NA) +
  scale_fill_manual(values = c("white", "lightgrey")) +
  ggtitle("Comparison of Cold vs. Control in the Hypothalamus, Separated by Sex and Section (s07-s12)") +
  labs(y = "Percentage", fill = "Treatment") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

#mouseIDlabels
cp_log_3_withslides %>%
  filter(BrainRegion == "HYPO", SlideNumber %in% c("s01", "s02", "s03", "s04", "s05", "s06")) %>%
  ggplot(aes(x = Cluster, y = percentage, group = interaction(Cluster, Treatment))) +
  facet_grid(Sex ~ SlideNumber) +
  geom_boxplot(aes(group = interaction(Cluster, Treatment), fill = Treatment)) +  # Keep outliers
  geom_point(aes(fill = Treatment), 
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
             size = 0.75, alpha = 0.6) +  # Add individual points
  geom_text(aes(label = MouseID, fill = Treatment), 
            position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
            vjust = -0.5, size = 3, color = "black") +  # Label all points
  scale_fill_manual(values = c("yellow", "lavender")) +
  ggtitle("Mouse Dataset: Cold vs. Control in the Hypothalamus, Separated by Sex and Section (s01-s06)") +
  labs(y = "Percentage", fill = "Treatment") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

#mouseIDlabels
cp_log_3_withslides %>%
  filter(BrainRegion == "HYPO", SlideNumber %in% c("s07", "s08", "s09", "s10", "s11", "s12")) %>%
  ggplot(aes(x = Cluster, y = percentage, group = interaction(Cluster, Treatment))) +
  facet_grid(Sex ~ SlideNumber) +
  geom_boxplot(aes(group = interaction(Cluster, Treatment), fill = Treatment)) +  # Keep outliers
  geom_point(aes(fill = Treatment), 
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
             size = 0.75, alpha = 0.6) +  # Add individual points
  geom_text(aes(label = MouseID, fill = Treatment), 
            position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8), 
            vjust = -0.5, size = 3, color = "black") +  # Label all points
  scale_fill_manual(values = c("yellow", "lavender")) +
  ggtitle("Mouse Dataset: Cold vs. Control in the Hypothalamus, Separated by Sex and Section (s01-s06)") +
  labs(y = "Percentage", fill = "Treatment") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16))


# Create the combined plot
raw_data_section %>%
  filter(BrainRegion == "HYPO", Treatment == "CON", Sex %in% c("M", "F")) %>%
  ggplot(aes(x = Cluster, y = percentage, group = interaction(Cluster, Sex))) +
  facet_wrap(~ SlideNumber, nrow = 2, scales = "free_x") + 
  geom_boxplot(aes(fill = Sex), outlier.shape = NA) +
  scale_fill_manual(values = c("lightpink", "lightblue")) +  # Different colors for Sex
  ggtitle("Male vs Female Controls in Hypothalamus, Separated by Section") +
  labs(fill = "Sex") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

# Create the combined plot
raw_data_section %>%
  filter(BrainRegion == "HYPO", Treatment == "COLD", Sex %in% c("M", "F")) %>%
  ggplot(aes(x = Cluster, y = percentage, group = interaction(Cluster, Sex))) +
  facet_wrap(~ SlideNumber, nrow = 2, scales = "free_x") + 
  geom_boxplot(aes(fill = Sex), outlier.shape = NA) +
  scale_fill_manual(values = c("lightpink", "lightblue")) +  # Different colors for Sex
  ggtitle("Mouse Dataset: Males vs Females in Cold in Hypothalamus, Separated by Section") +
  labs(fill = "Sex") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

cp_log_3 %>% 
  ggplot(aes(x = Cluster, y = percentage, fill = BrainRegion)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2, outlier.color = "black") +  # Include outliers
  scale_fill_manual(values = c("lightblue", "lightgreen", "lightpink")) +
  ggtitle("Comparison of Cluster Percentages Across Brain Regions") +
  labs(x = "Cluster", y = "Percentage", fill = "Brain Region") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

cp_log_3 %>% 
  ggplot(aes(x = Cluster, y = percentage, fill = BrainRegion)) +
  geom_boxplot(outlier.shape = 21, outlier.size = 2, outlier.color = "black") +  # Include outliers
  scale_fill_manual(values = c("lightblue", "lightgreen", "lightpink")) +
  facet_grid(Sex ~ Treatment) +  # Split by Sex and Treatment
  ggtitle("Comparison of Cluster Percentages Across Brain Regions by Sex and Treatment") +
  labs(x = "Cluster", y = "Percentage", fill = "Brain Region") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))



brainRegionPlot <- cp_log_3 %>% 
  ggplot(aes(x = interaction(Cluster, BrainRegion, lex.order = TRUE), y = percentage)) +
  geom_boxplot(aes(fill = BrainRegion), outlier.shape = NA) +  # Boxplots split by Cluster and BrainRegion
  geom_point(aes(color = Sex), 
             position = position_jitter(width = 0.2, height = 0), size = 2, alpha = 0.7) +  # Points jittered within each group
  scale_fill_manual(values = c("ARC" = "white", "HYPO" = "lightgrey", "PVN" = "darkgrey")) +  # Colors for BrainRegion
  scale_color_manual(values = c("M" = "blue", "F" = "pink")) +  # Colors for Sex
  ggtitle("Cluster Percentages Across Brain Regions") +
  labs(x = "Cluster and Brain Region", y = "Percentage", fill = "Brain Region", color = "Sex") +
  scale_x_discrete(labels = function(x) gsub("\\.", " - ", x)) +  # Replace "." in interaction labels with " - " for readability
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

brainRegionPlot_wfacet <- cp_log_3 %>% 
  ggplot(aes(x = BrainRegion, y = percentage)) +
  geom_boxplot(aes(fill = BrainRegion), outlier.shape = NA) +  # Boxplots split by BrainRegion
  geom_point(aes(color = Sex), 
             position = position_jitter(width = 0.2, height = 0), size = 2, alpha = 0.7) +  # Points jittered within each group
  scale_fill_manual(values = c("ARC" = "white", "HYPO" = "lightgrey", "PVN" = "darkgrey")) +  # Colors for BrainRegion
  scale_color_manual(values = c("M" = "blue", "F" = "pink")) +  # Colors for Sex
  ggtitle("Cluster Percentages Across Brain Regions") +
  labs(x = "Brain Region", y = "Percentage", fill = "Brain Region", color = "Sex") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold")) +
  facet_wrap(~ Cluster, scales = "free_x")  # Facet by Cluster

brainRegionPlot_wfacet

treatmentPlot_wfacet <- cp_log_3 %>% 
  ggplot(aes(x = Cluster, y = percentage)) +
  geom_boxplot(aes(fill = Treatment), outlier.shape = NA) +  # Boxplots split by Treatment
  geom_point(aes(color = Sex, group = interaction(Cluster, Treatment)), 
             position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), 
             size = 2, alpha = 0.5) +  # Points jittered and aligned with boxplots
  scale_fill_manual(values = c("COLD" = "white", "CON" = "lightgrey")) +  # Colors for Treatment
  scale_color_manual(values = c("M" = "blue", "F" = "pink")) +  # Colors for Sex
  ggtitle("Cluster Percentages Across Brain Regions") +
  labs(x = "Brain Region", y = "Percentage", fill = "Treatment", color = "Sex") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold")) +
  facet_wrap(~ BrainRegion, scales = "free_x")  # Facet by Brain Region


treatmentPlot_wfacet

brainRegionPlot_treatment <- cp_log_3 %>% 
  ggplot(aes(x = interaction(Cluster, BrainRegion, lex.order = TRUE), y = percentage)) +
  geom_boxplot(aes(fill = BrainRegion), outlier.shape = NA) +  # Boxplots split by Cluster and BrainRegion
  geom_point(aes(color = Sex), 
             position = position_jitter(width = 0.2, height = 0), size = 2, alpha = 0.7) +  # Points jittered within each group
  scale_fill_manual(values = c("ARC" = "white", "HYPO" = "lightgrey", "PVN" = "darkgrey")) +  # Colors for BrainRegion
  scale_color_manual(values = c("M" = "blue", "F" = "pink")) +  # Colors for Sex
  ggtitle("Cluster Percentages Across Brain Regions") +
  labs(x = "Cluster and Brain Region", y = "Percentage of Cluster", fill = "Brain Region", color = "Sex") +
  scale_x_discrete(labels = function(x) gsub("\\.", " - ", x)) +  # Replace "." in interaction labels with " - " for readability
  facet_wrap(~ Treatment, ncol = 2) +  # Create separate panels for COLD and CON treatments
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

brainRegionPlot_treatment

brainRegionPlot_treatment <- cp_log_3 %>% 
  ggplot(aes(x = BrainRegion, y = percentage)) +  # Use BrainRegion on x-axis
  geom_boxplot(aes(fill = BrainRegion), outlier.shape = NA) +  # Boxplots split by BrainRegion
  geom_point(aes(color = Sex), 
             position = position_jitter(width = 0.2, height = 0), size = 2, alpha = 0.7) +  # Points jittered within each group
  scale_fill_manual(values = c("ARC" = "white", "HYPO" = "lightgrey", "PVN" = "darkgrey")) +  # Colors for BrainRegion
  scale_color_manual(values = c("M" = "blue", "F" = "pink")) +  # Colors for Sex
  ggtitle("Cluster Percentages Across Brain Regions") +
  labs(x = "Brain Region", y = "Percentage of Cluster", fill = "Brain Region", color = "Sex") +
  facet_grid(Treatment ~ Cluster, scales = "free_x") +  # Facet by Cluster (rows) and Treatment (columns)
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

brainRegionPlot_treatment_males <- cp_log_3 %>%
  filter(Sex == "M") %>%  # Filter for males only
  ggplot(aes(x = BrainRegion, y = percentage)) +  # Use BrainRegion on x-axis
  geom_boxplot(aes(fill = BrainRegion), outlier.shape = NA) +  # Boxplots split by BrainRegion
  geom_point(aes(color = Sex), 
             position = position_jitter(width = 0.2, height = 0), size = 2, alpha = 0.7) +  # Points jittered within each group
  scale_fill_manual(values = c("ARC" = "white", "HYPO" = "lightgrey", "PVN" = "darkgrey")) +  # Colors for BrainRegion
  scale_color_manual(values = c("M" = "blue")) +  # Only male color is needed
  ggtitle("Cluster Percentages Across Brain Regions (Males Only)") +
  labs(x = "Brain Region", y = "Percentage of Cluster", fill = "Brain Region", color = "Sex") +
  facet_grid(Cluster ~ Treatment, scales = "free_x") +  # Facet by Cluster (rows) and Treatment (columns)
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

brainRegionPlot_treatment_males

brainRegionPlot_treatment_females <- cp_log_3 %>%
  filter(Sex == "F") %>%  # Filter for males only
  ggplot(aes(x = BrainRegion, y = percentage)) +  # Use BrainRegion on x-axis
  geom_boxplot(aes(fill = BrainRegion), outlier.shape = NA) +  # Boxplots split by BrainRegion
  geom_point(aes(color = Sex), 
             position = position_jitter(width = 0.2, height = 0), size = 2, alpha = 0.7) +  # Points jittered within each group
  scale_fill_manual(values = c("ARC" = "white", "HYPO" = "lightgrey", "PVN" = "darkgrey")) +  # Colors for BrainRegion
  scale_color_manual(values = c("F" = "pink")) +  # Only male color is needed
  ggtitle("Cluster Percentages Across Brain Regions (Females Only)") +
  labs(x = "Brain Region", y = "Percentage of Cluster", fill = "Brain Region", color = "Sex") +
  facet_grid(Cluster ~ Treatment, scales = "free_x") +  # Facet by Cluster (rows) and Treatment (columns)
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

cp_log_3 %>%
  ggplot(aes(x = interaction(Cluster, BrainRegion), y = percentage)) +
  geom_boxplot(aes(fill = BrainRegion), outlier.shape = NA) +
  geom_point(aes(color = Sex), position = position_jitter(width = 0.2, height = 0), size = 2, alpha = 0.7) +
  geom_signif(comparisons = list(c("Ramified.ARC", "Ramified.HYPO")), 
              annotations = c("**"),  # Custom annotations for significance
              y_position = c(1),    # Position of the significance bars
              tip_length = 0.03) +   # Length of the lines
  scale_fill_manual(values = c("ARC" = "white", "HYPO" = "lightgrey", "PVN" = "darkgrey")) +
  scale_color_manual(values = c("M" = "blue", "F" = "pink")) +
  labs(x = "Cluster and Brain Region", y = "Percentage", fill = "Brain Region", color = "Sex") +
  theme_bw(base_size = 14)

ggsave("/Users/alexlawson/Desktop/Masters-Work/boxplots/brainRegionPlotwtreatment-faceted-males.png", plot = brainRegionPlot_treatment_males, dpi = 300, width = 8, height = 6, units = "in")







# ---- Stats ----   
stats_input <- cp_log_3 
stats_input$MouseID <- factor(stats_input$MouseID)
stats_input$Cluster <- factor(stats_input$Cluster)
stats_input$Treatment <- factor(stats_input$Treatment)
stats_input$Sex <- factor(stats_input$Sex)
View(stats_input)
stats_testing_all <- stats_cluster.animal(stats_input,
                                          model = "percentage ~ Cluster*Treatment*BrainRegion*Sex + (1|MouseID)", 
                                          posthoc1 = "~Treatment|Cluster|BrainRegion", 
                                          posthoc2 = "~Treatment|Cluster|BrainRegion|Sex", adjust = "sidak")
stats_testing_all[[1]]
stats_testing_all[[2]]
stats_testing_all[[3]]
stats_testing_all[[4]]
stats_testing_all[[5]]

stats_testing_control <- stats_cluster.animal(data = stats_input %>% filter(Treatment=="CON"),
                                              model = "percentage ~ Cluster*BrainRegion*Sex + (1|MouseID)", 
                                              posthoc1 = "~Sex|Cluster", 
                                              posthoc2 = "~Sex|Cluster|BrainRegion", adjust = "sidak")
stats_testing_control[[1]]
stats_testing_control[[2]]
stats_testing_control[[3]]

stats_testing_treatment <- stats_cluster.animal(data = stats_input %>% filter(Treatment=="COLD"),
                                              model = "percentage ~ Cluster*BrainRegion*Sex + (1|MouseID)", 
                                              posthoc1 = "~Sex|Cluster", 
                                              posthoc2 = "~Sex|Cluster|BrainRegion", adjust = "sidak")
stats_testing_treatment[[1]]
stats_testing_treatment[[2]]
stats_testing_treatment[[3]]

write.csv(stats_testing_all[[3]], "/Users/alexlawson/Desktop/Masters-Work/statistical-results/all-results-separated-by-brain-region-and-sex.csv", row.names = FALSE)


# Calculate median and standard error grouped by Sex and BrainRegion
result <- stats_input %>%
  filter(Treatment == "COLD") %>% # Filter for rows where Treatment is "COLD"
  group_by(Sex, BrainRegion, Cluster) %>%  # Group by Sex and BrainRegion
  summarise(
    Mean = mean(percentage, na.rm = TRUE),                   # Replace `value` with the relevant column
    StdError = sd(percentage, na.rm = TRUE) / sqrt(n()),        # Standard Error = SD / sqrt(sample size)
    .groups = "drop"
  )

# View the result
print(result)
stats_testing_M <- stats_cluster.animal(data = stats_input %>% filter(Sex=="M"&BrainRegion=="PVN"),
                                          model = "percentage ~ Cluster*Treatment + (1|MouseID)", 
                                          posthoc1 = "~Treatment|Cluster", 
                                          posthoc2 = "~Treatment|Cluster", adjust = "sidak")
stats_testing_M[[2]]
stats_testing_M[[3]]

write.csv(stats_testing_M[[3]], "/Users/alexlawson/Desktop/Masters-Work/statistical-results/males-only-separated-by-brain-region.csv", row.names = FALSE)


stats_testing_F <- stats_cluster.animal(data = stats_input %>% filter(Sex=="F"),
                                        model = "percentage ~ Cluster*Treatment*BrainRegion + (1|MouseID)", 
                                        posthoc1 = "~Treatment|Cluster", 
                                        posthoc2 = "~Treatment|Cluster|BrainRegion", adjust = "sidak")
stats_testing_F[[2]]
stats_testing_F[[3]]

write.csv(stats_testing_F[[2]], "/Users/alexlawson/Desktop/Masters-Work/statistical-results/females-only.csv", row.names = FALSE)


stats_testing_hypo <- stats_cluster.animal(data = stats_input %>% filter(BrainRegion=="HYPO"),
                                        model = "percentage ~ Cluster*Treatment*Sex + (1|MouseID)", 
                                        posthoc1 = "~Treatment|Cluster", 
                                        posthoc2 = "~Treatment|Cluster|Sex", adjust = "sidak")
stats_testing_hypo[[2]]
stats_testing_hypo[[3]]

write.csv(stats_testing_hypo[[3]], "/Users/alexlawson/Desktop/Masters-Work/statistical-results/hypothalamus-only-separated-by-sex.csv", row.names = FALSE)

stats_testing_pvn <- stats_cluster.animal(data = stats_input %>% filter(BrainRegion=="PVN"),
                                           model = "percentage ~ Cluster*Treatment*Sex + (1|MouseID)", 
                                           posthoc1 = "~Treatment|Cluster", 
                                           posthoc2 = "~Treatment|Cluster|Sex", adjust = "sidak")
stats_testing_pvn[[2]]
stats_testing_pvn[[3]]

write.csv(stats_testing_pvn[[2]], "/Users/alexlawson/Desktop/Masters-Work/statistical-results/pvn-only.csv", row.names = FALSE)


stats_testing_arc <- stats_cluster.animal(data = stats_input %>% filter(BrainRegion=="ARC"),
                                          model = "percentage ~ Cluster*Treatment*Sex + (1|MouseID)", 
                                          posthoc1 = "~Treatment|Cluster", 
                                          posthoc2 = "~Treatment|Cluster|Sex", adjust = "sidak")
stats_testing_arc[[2]]
stats_testing_arc[[3]]

write.csv(stats_testing_arc[[3]], "/Users/alexlawson/Desktop/Masters-Work/statistical-results/arc-only-separated-by-sex.csv", row.names = FALSE)



stats_input_pvn_M <- stats_input %>% filter(BrainRegion=="PVN"&Sex=="M")
stats_testing_pvn_M <- stats_cluster.animal(data = stats_input_pvn_M,
                                          model = "percentage ~ Cluster*Treatment + (1|MouseID)", 
                                          posthoc1 = "~Treatment|Cluster", 
                                          posthoc2 = "~Treatment|Cluster", adjust = "sidak")

stats_testing_pvn_M[[1]]
stats_testing_pvn_M[[3]]

stats_input_pvn_F <- stats_input %>% filter(BrainRegion=="PVN"&Sex=="F")
stats_testing_pvn_F <- stats_cluster.animal(data = stats_input_pvn_F,
                                            model = "percentage ~ Cluster*Treatment + (1|MouseID)", 
                                            posthoc1 = "~Treatment|Cluster", 
                                            posthoc2 = "~Treatment|Cluster", adjust = "sidak")

stats_testing_pvn_F[[1]]
stats_testing_pvn_F[[3]]



stats_input_slides <- raw_data_section 
stats_input_slides$MouseID <- factor(stats_input_slides$MouseID)
stats_input_slides$Cluster <- factor(stats_input_slides$Cluster)
stats_input_slides$Treatment <- factor(stats_input_slides$Treatment)
stats_input_slides$SlideNumber <- factor(stats_input_slides$SlideNumber)
View(stats_input_slides)

stats_input_slides_hypothalamus <- stats_input_slides %>% filter(BrainRegion=="HYPO")
stats_input_slides_hypo_con <- stats_input_slides_hypothalamus %>% filter(Treatment=="CON")

stats_testing_hypothalamus_con <- stats_cluster.animal(data = stats_input_slides_hypo_con,
                                                   model = "percentage ~ Cluster*SlideNumber*Sex + (1|MouseID)", 
                                                   posthoc1 = "~Sex|Cluster|SlideNumber", 
                                                   posthoc2 = "~Sex|Cluster|SlideNumber", adjust="sidak")
stats_testing_hypothalamus_con[[2]]

stats_testing_hypothalamus <- stats_cluster.animal(data = stats_input_slides_hypo_con %>% filter(SlideNumber=="s01"),
                                                   model = "percentage ~ Cluster*Sex + (1|MouseID)", 
                                                   posthoc1 = "~Sex|Cluster", 
                                                   posthoc2 = "~Sex|Cluster", adjust="sidak")


# Vector of all slide numbers
slide_numbers <- sprintf("s%02d", 1:12)

# Function to run the model for each slide
run_stats <- function(slide) {
  result <- stats_cluster.animal(
    data = stats_input_slides_hypo_con %>% filter(SlideNumber == slide),
    model = "percentage ~ Cluster*Sex + (1|MouseID)",
    posthoc1 = "~Sex|Cluster",
    posthoc2 = "~Sex|Cluster",
    adjust = "sidak"
  )
  
  posthoc <- result[[2]]
  posthoc$SlideNumber <- slide  # Add slide info for tracking
  return(posthoc)
}

# Run the function for each slide and bind results
all_stats <- map_dfr(slide_numbers, run_stats)

write.csv(all_stats, "hypothalamus_posthoc_stats.csv", row.names = FALSE)


stats_input_slides_hypothalamus

# Vector of all slide numbers
slide_numbers <- sprintf("s%02d", 1:12)

# Function to run the model for each section
run_stats <- function(slide) {
  result <- stats_cluster.animal(
    data = stats_input_slides_hypothalamus %>% filter(SlideNumber == slide),
    model = "percentage ~ Treatment*Cluster*Sex + (1|MouseID)",
    posthoc1 = "~Treatment|Cluster",
    posthoc2 = "~Treatment|Cluster|Sex",
    adjust = "sidak"
  )
  
  posthoc <- result[[2]]
  posthoc$SlideNumber <- slide  # Add slide info for tracking
  return(posthoc)
}

# Function to run the model for each slide
run_stats_sex <- function(slide) {
  result <- stats_cluster.animal(
    data = stats_input_slides_hypothalamus %>% filter(SlideNumber == slide),
    model = "percentage ~ Treatment*Cluster*Sex + (1|MouseID)",
    posthoc1 = "~Treatment|Cluster",
    posthoc2 = "~Treatment|Cluster|Sex",
    adjust = "sidak"
  )
  
  posthoc <- result[[3]]
  posthoc$SlideNumber <- slide  # Add slide info for tracking
  return(posthoc)
}

# Run the function for each slide and bind results
all_stats <- map_dfr(slide_numbers, run_stats)
all_stats_sex <- map_dfr(slide_numbers, run_stats_sex)


write.csv(all_stats_sex, "hypothalamus_posthoc_stats_treatment_comparison_sex.csv", row.names = FALSE)


stats_testing_hypothalamus <- stats_cluster.animal(data = stats_input_slides_hypothalamus,
                                                   model = "percentage ~ Cluster*Treatment*SlideNumber*Sex + (1|MouseID)", 
                                                   posthoc1 = "~Treatment|Cluster|SlideNumber", 
                                                   posthoc2 = "~Treatment|Cluster|SlideNumber|Sex", adjust="sidak")

stats_testing_hypothalamus[[1]]
stats_testing_hypothalamus[[2]]
stats_testing_hypothalamus[[3]]

write.csv(stats_testing_hypothalamus[[3]], "/Users/alexlawson/Desktop/Masters-Work/statistical-results/hypothalamus-separated-by-section-and-sex.csv", row.names = FALSE)


stats_input_hypothalamus_s01 <-stats_input_slides_hypothalamus %>% filter(SlideNumber=="s01")
stats_input_hypothalamus_s02 <-stats_input_slides_hypothalamus %>% filter(SlideNumber=="s02")
stats_input_hypothalamus_s03 <-stats_input_slides_hypothalamus %>% filter(SlideNumber=="s03")
stats_input_hypothalamus_s04 <-stats_input_slides_hypothalamus %>% filter(SlideNumber=="s04")
stats_input_hypothalamus_s05 <-stats_input_slides_hypothalamus %>% filter(SlideNumber=="s05")
stats_input_hypothalamus_s06 <-stats_input_slides_hypothalamus %>% filter(SlideNumber=="s06")
stats_input_hypothalamus_s07 <-stats_input_slides_hypothalamus %>% filter(SlideNumber=="s07")
stats_input_hypothalamus_s08 <-stats_input_slides_hypothalamus %>% filter(SlideNumber=="s08")
stats_input_hypothalamus_s09 <-stats_input_slides_hypothalamus %>% filter(SlideNumber=="s09")
stats_input_hypothalamus_s10 <-stats_input_slides_hypothalamus %>% filter(SlideNumber=="s10")
stats_input_hypothalamus_s11 <-stats_input_slides_hypothalamus %>% filter(SlideNumber=="s11")
stats_input_hypothalamus_s12 <-stats_input_slides_hypothalamus %>% filter(SlideNumber=="s12")



stats_testing_hypothalamus_s01 <- stats_cluster.animal(data = stats_input_hypothalamus_s01,
                                                   model = "percentage ~ Cluster*Treatment*Sex+ (1|MouseID)", 
                                                   posthoc1 = "~Treatment|Cluster", 
                                                   posthoc2 = "~Treatment|Cluster|Sex", adjust="sidak")
stats_testing_hypothalamus_s01[[2]]
stats_testing_hypothalamus_s01[[3]]

stats_testing_hypothalamus_s02 <- stats_cluster.animal(data = stats_input_hypothalamus_s02,
                                                       model = "percentage ~ Cluster*Treatment*Sex+ (1|MouseID)", 
                                                       posthoc1 = "~Treatment|Cluster", 
                                                       posthoc2 = "~Treatment|Cluster|Sex", adjust="sidak")

stats_testing_hypothalamus_s02[[2]]
stats_testing_hypothalamus_s02[[3]] 

stats_testing_hypothalamus_s03 <- stats_cluster.animal(data = stats_input_hypothalamus_s03,
                                                       model = "percentage ~ Cluster*Treatment*Sex+ (1|MouseID)", 
                                                       posthoc1 = "~Treatment|Cluster", 
                                                       posthoc2 = "~Treatment|Cluster|Sex", adjust="sidak")
stats_testing_hypothalamus_s03[[2]]
stats_testing_hypothalamus_s03[[3]]


stats_testing_hypothalamus_s04 <- stats_cluster.animal(data = stats_input_hypothalamus_s04,
                                                       model = "percentage ~ Cluster*Treatment*Sex+ (1|MouseID)", 
                                                       posthoc1 = "~Treatment|Cluster", 
                                                       posthoc2 = "~Treatment|Cluster|Sex", adjust="sidak")
stats_testing_hypothalamus_s04[[2]]
stats_testing_hypothalamus_s04[[3]]


stats_testing_hypothalamus_s05 <- stats_cluster.animal(data = stats_input_hypothalamus_s05,
                                                       model = "percentage ~ Cluster*Treatment*Sex+ (1|MouseID)", 
                                                       posthoc1 = "~Treatment|Cluster", 
                                                       posthoc2 = "~Treatment|Cluster|Sex", adjust="sidak")

stats_testing_hypothalamus_s05[[2]]
stats_testing_hypothalamus_s05[[3]]

stats_testing_hypothalamus_s06 <- stats_cluster.animal(data = stats_input_hypothalamus_s06,
                                                       model = "percentage ~ Cluster*Treatment*Sex+ (1|MouseID)", 
                                                       posthoc1 = "~Treatment|Cluster", 
                                                       posthoc2 = "~Treatment|Cluster|Sex", adjust="sidak")

stats_testing_hypothalamus_s06[[2]]
stats_testing_hypothalamus_s06[[3]]

stats_testing_hypothalamus_s07 <- stats_cluster.animal(data = stats_input_hypothalamus_s07,
                                                       model = "percentage ~ Cluster*Treatment*Sex+ (1|MouseID)", 
                                                       posthoc1 = "~Treatment|Cluster", 
                                                       posthoc2 = "~Treatment|Cluster|Sex", adjust="sidak")
stats_testing_hypothalamus_s07[[2]]
stats_testing_hypothalamus_s07[[3]]

stats_testing_hypothalamus_s08 <- stats_cluster.animal(data = stats_input_hypothalamus_s08,
                                                       model = "percentage ~ Cluster*Treatment*Sex+ (1|MouseID)", 
                                                       posthoc1 = "~Treatment|Cluster", 
                                                       posthoc2 = "~Treatment|Cluster|Sex", adjust="sidak")

stats_testing_hypothalamus_s08[[2]]
stats_testing_hypothalamus_s08[[3]]

stats_testing_hypothalamus_s09 <- stats_cluster.animal(data = stats_input_hypothalamus_s09,
                                                       model = "percentage ~ Cluster*Treatment*Sex+ (1|MouseID)", 
                                                       posthoc1 = "~Treatment|Cluster", 
                                                       posthoc2 = "~Treatment|Cluster|Sex", adjust="sidak")

stats_testing_hypothalamus_s09[[2]]
stats_testing_hypothalamus_s09[[3]]

stats_testing_hypothalamus_s10 <- stats_cluster.animal(data = stats_input_hypothalamus_s10,
                                                       model = "percentage ~ Cluster*Treatment*Sex+ (1|MouseID)", 
                                                       posthoc1 = "~Treatment|Cluster", 
                                                       posthoc2 = "~Treatment|Cluster|Sex", adjust="sidak")

stats_testing_hypothalamus_s10[[2]]
stats_testing_hypothalamus_s10[[3]]

stats_testing_hypothalamus_s11 <- stats_cluster.animal(data = stats_input_hypothalamus_s11,
                                                       model = "percentage ~ Cluster*Treatment*Sex+ (1|MouseID)", 
                                                       posthoc1 = "~Treatment|Cluster", 
                                                       posthoc2 = "~Treatment|Cluster|Sex", adjust="sidak")

stats_testing_hypothalamus_s11[[2]]
stats_testing_hypothalamus_s11[[3]]
stats_testing_hypothalamus_s12 <- stats_cluster.animal(data = stats_input_hypothalamus_s12,
                                                       model = "percentage ~ Cluster*Treatment*Sex+ (1|MouseID)", 
                                                       posthoc1 = "~Treatment|Cluster", 
                                                       posthoc2 = "~Treatment|Cluster|Sex", adjust="sidak")

stats_testing_hypothalamus_s12[[2]]
stats_testing_hypothalamus_s12[[3]]


View(stats_testing_hypothalamus_s12[[3]])
stats_testing_hypothalamus_s02[[3]]
stats_testing_hypothalamus_s03[[3]]
stats_testing_hypothalamus_s04[[3]]
stats_testing_hypothalamus_s05[[3]]
stats_testing_hypothalamus_s06[[3]]
stats_testing_hypothalamus_s07[[3]]
stats_testing_hypothalamus_s08[[3]]
stats_testing_hypothalamus_s09[[3]]
stats_testing_hypothalamus_s10[[3]]
stats_testing_hypothalamus_s11[[3]]
stats_testing_hypothalamus_s12[[3]]


stats_testing_brain_comparison<- stats_cluster.animal(data = stats_input,
                                             model = "percentage ~ Cluster*BrainRegion*Sex + (1|MouseID)", 
                                             posthoc1 = "~BrainRegion|Cluster", 
                                             posthoc2 = "~BrainRegion|Cluster|Sex", adjust="sidak")


stats_testing_brain_comparison[[2]]
stats_testing_brain_comparison[[3]]

write.csv(stats_testing_brain_comparison[[2]], "/Users/alexlawson/Desktop/Masters-Work/statistical-results/brainregion-comparison.csv", row.names = FALSE)



stats_input_hypoarc <- stats_input %>% filter(BrainRegion=="ARC"|BrainRegion=="HYPO")
stats_testing_hypoarc<- stats_cluster.animal(data = stats_input_hypoarc,
                                                       model = "percentage ~ Cluster*BrainRegion*Sex + + (1|MouseID)", 
                                                       posthoc1 = "~BrainRegion|Cluster", 
                                                       posthoc2 = "~BrainRegion|Cluster|Sex", adjust="sidak")


stats_testing_hypoarc[[2]]
stats_testing_hypoarc[[3]]
stats_testing_hypoarc[[3]]

write.csv(stats_testing_hypoarc[[3]], "/Users/alexlawson/Desktop/Masters-Work/statistical-results/hypothalamus-arc-separated-by-sex.csv", row.names = FALSE)


stats_input_hypopvn <- stats_input %>% filter(BrainRegion=="PVN"|BrainRegion=="HYPO")
stats_testing_hypopvn<- stats_cluster.animal(data = stats_input_hypopvn,
                                             model = "percentage ~ Cluster*BrainRegion*Sex", 
                                             posthoc1 = "~BrainRegion|Cluster", 
                                             posthoc2 = "~BrainRegion|Cluster|Sex", adjust="sidak")

stats_testing_hypopvn[[2]]
stats_testing_hypopvn[[3]]
write.csv(stats_testing_hypopvn[[2]], "/Users/alexlawson/Desktop/Masters-Work/statistical-results/hypothalamus-pvn.csv", row.names = FALSE)



# Assuming your dataset is named 'data' and the brain region column is named 'Region'
filtered_data <- data %>%
  filter(Region %in% c("ARC", "HYPO"))


# ---- Color By Cluster ----   
#females s03/s04

colorbycluster_log_3clusters_arc_2 <- pca_kmeans03_log %>% 
  filter(BrainRegion=="ARC",MouseID=="MR175", SlideNumber=="s04") %>% select(c(Cluster, ID))
write.csv(colorbycluster_log_3clusters_arc_2, "/Users/alexlawson/Desktop/GliaData/final-data/cbc-arc/ARC_MR175_s04.csv", row.names = FALSE)

library(dplyr)

# Define unique combinations of MouseID and SlideNumber for BrainRegion "PVN"
combinations <- expand.grid(
  MouseID = c("MR138", "MR139", "MR140", "MR141", "MR156", "MR157", "MR158", "MR159", 
              "MR161", "MR169", "MR172", "MR175", "MR196", "MR197", "MR198", "MR199", 
              "MR273", "MR274", "MR277", "MR278", "MR284", "MR287", "MR290", "MR294"),
  SlideNumber = c("s03", "s04", "s05", "s06"),
  BrainRegion = "PVN"
)

# Output directory
output_directory <- "/Users/alexlawson/Desktop/pvn-fem/cbc-data-pvn"

# Iterate through each combination and save CSVs
for (i in 1:nrow(combinations)) {
  mouse_id <- combinations$MouseID[i]
  slide_number <- combinations$SlideNumber[i]
  
  # Filter the data
  filtered_data <- pca_kmeans03_log %>%
    filter(BrainRegion == "PVN", MouseID == mouse_id, SlideNumber == slide_number)
  
  # Check if the filtered data has any rows
  if (nrow(filtered_data) > 0) {
    filtered_data <- filtered_data %>% dplyr::select(Cluster, ID)
    
    # Create the file name
    file_name <- file.path(output_directory, paste0("PVN_", mouse_id, "_", slide_number, ".csv"))
    
    # Write to CSV
    write.csv(filtered_data, file_name, row.names = FALSE)
    
    message("Saved: ", file_name)
  } else {
    message("No data found for: PVN, ", mouse_id, ", ", slide_number)
  }
}

colnames(pca_kmeans03_log)


stats_testing_control <- stats_cluster.animal(data = stats_input %>% filter(Treatment=="CON"),
                                               model = "percentage ~ Cluster*Sex*BrainRegion + (1|MouseID)", 
                                               posthoc1 = "~Sex|Cluster", 
                                               posthoc2 = "~Sex|Cluster|BrainRegion", adjust = "sidak")
stats_testing_control[[1]]
stats_testing_control[[2]]
stats_testing_control[[3]]


stats_testing_cold <- stats_cluster.animal(data = stats_input %>% filter(Treatment=="COLD"),
                                              model = "percentage ~ Cluster*Sex*BrainRegion + (1|MouseID)", 
                                              posthoc1 = "~Sex|Cluster", 
                                              posthoc2 = "~Sex|Cluster|BrainRegion", adjust = "sidak")
stats_testing_cold[[1]]
stats_testing_cold[[2]]
stats_testing_cold[[3]]
# write.csv(stats_testing_control[[3]], "/Users/alexlawson/Desktop/GliaData/final-data/stats-results/sex-cluster-brainregion-control", row.names = FALSE)
# 
# stats_testing_cold <- stats_cluster.animal(data = stats_input %>% filter(Treatment=="COLD"),
#                                            model = "percentage ~ Cluster*BrainRegion*Sex + (1|MouseID)", 
#                                            posthoc1 = "~Sex|Cluster", 
#                                            posthoc2 = "~Sex|Cluster|BrainRegion", adjust="bonferroni")
# stats_testing_cold[[1]]
# stats_testing_cold[[2]]
# stats_testing_cold[[3]]
# 
# stats_testing_all <- stats_cluster.animal(data = stats_input,
#                                           model = "percentage ~ Cluster*Treatment*Sex*BrainRegion + (1|MouseID)", 
#                                           posthoc1 = "~Treatment|Cluster|BrainRegion", 
#                                           posthoc2 = "~Treatment|Cluster|Sex|BrainRegion", adjust = "bonferroni")
# 
# stats_testing_all[[1]]
# stats_testing_all[[2]]
# stats_testing_all[[3]]
# write.csv(stats_testing_all[[2]], "/Users/alexlawson/Desktop/GliaData/final-data/stats-results/treatment-cluster-brainregion-notsexsep.csv", row.names = FALSE)
# 
# stats_testing_males <- stats_cluster.animal(data = stats_input %>% filter(Sex=="M"), 
#                                             model = "percentage ~ Cluster*Treatment*BrainRegion + (1|MouseID)", 
#                                             posthoc1 = "~Treatment|Cluster", 
#                                             posthoc2 = "~Treatment|Cluster|BrainRegion", adjust = "bonferroni")
# 
# stats_testing_males[[1]]
# stats_testing_males[[2]]
# stats_testing_males[[3]]
# 
# 
# 
# stats_testing_fem <- stats_cluster.animal(data = stats_input %>% filter(Sex=="F"), 
#                                           model = "percentage ~ Cluster*Treatment*BrainRegion + (1|MouseID)", 
#                                           posthoc1 = "~Treatment|Cluster", 
#                                           posthoc2 = "~Treatment|Cluster|BrainRegion", adjust = "bonferroni")
# 
# stats_testing_fem[[1]]
# stats_testing_fem[[2]]
# stats_testing_fem[[3]]
# 
# stats_input_pvn <- stats_input %>% filter(BrainRegion=="PVN")
# stats_testing_pvn <- stats_cluster.animal(data = stats_input_pvn,
#                                                model = "percentage ~ Cluster*Treatment*Sex + (1|MouseID)", 
#                                                posthoc1 = "~Treatment|Cluster", 
#                                                posthoc2 = "~Treatment|Cluster|Sex", adjust = "bonferroni")
# 
# stats_testing_pvn[[1]]
# stats_testing_pvn[[2]]
# stats_testing_pvn[[3]]
# 
# 
# 
# write.csv(stats_testing_hypothalamus[[2]], "/Users/alexlawson/Desktop/GliaData/final-data/stats-results/treatment-cluster-sex-hypothalamus.csv", row.names = FALSE)
# 
# 
# stats_testing_hypothalamus_control <- stats_cluster.animal(data = stats_input_slides_hypothalamus %>%filter(Treatment=="CON"),  
#                                                    model = "percentage ~ Cluster*SlideNumber*Sex", 
#                                                    posthoc1 = "~Sex|Cluster|SlideNumber", 
#                                                    posthoc2 = "~Cluster|SlideNumber")
# 
# stats_testing_hypothalamus_control[[1]]
# stats_testing_hypothalamus_control[[2]]
# stats_testing_hypothalamus_control[[3]]
# 
# stats_testing_hypothalamus_cold <- stats_cluster.animal(data = stats_input_slides_hypothalamus %>%filter(Treatment=="COLD"),  
#                                                            model = "percentage ~ Cluster*SlideNumber*Sex", 
#                                                            posthoc1 = "~Sex|Cluster|SlideNumber", 
#                                                            posthoc2 = "~Sex|Cluster")
# 
# stats_testing_hypothalamus_cold[[1]]
# stats_testing_hypothalamus_cold[[2]]
# stats_testing_hypothalamus_cold[[3]]
# 
# s03_hypo <-stats_input_slides_hypothalamus %>% filter(SlideNumber=="s03")
# stats_testing_hypothalamus_s03 <- stats_cluster.animal(data = s03_hypo,  
#                                                    model = "percentage ~ Cluster*Treatment*Sex", 
#                                                    posthoc1 = "~Treatment|Cluster", 
#                                                    posthoc2 = "~Treatment|Cluster|Sex")
# 
# stats_testing_hypothalamus_s03[[2]]
# stats_testing_hypothalamus_s03[[3]]
# 
# stats_input_features <- raw_data01 %>% 
#   group_by(MouseID, Sex, Treatment, BrainRegion) %>% 
#   summarise(across("Foreground.pixels":"Maximum.branch.length", ~mean(.x))) %>% 
#   gather(Measure, Value, "Foreground.pixels":"Maximum.branch.length")
# 
# 
# stats_input_features$Sex <- factor(stats_input_features$Sex)
# stats_input_features$BrainRegion <- factor(stats_input_features$BrainRegion)
# stats_testing_features <- stats_morphologymeasures.animal(data = stats_input_features, 
#                                                  model = "Value ~ Treatment*BrainRegion*Sex", type="lm",
#                                                  posthoc1 = "~Treatment|BrainRegion", 
#                                                  posthoc2 = "~Treatment|BrainRegion|Sex", "bonferroni")
# stats_testing_features[[2]]
# stats_testing_features[[3]]
# 
# stats_input_features_PVN <- stats_input_features %>% filter(BrainRegion == "PVN" & Sex == "M")
# 
# stats_testing_features_PVN <- stats_morphologymeasures.animal(data =stats_input_features_PVN, 
#                                                           model = "Value ~ Treatment", type="lm",
#                                                           posthoc1 = "~Treatment", 
#                                                           posthoc2 = "~Treatment", "bonferroni")
# View(stats_testing_features_PVN[[3]])
# 
# 
# stats_input_features_filtered <- stats_input_features %>%
#   filter(Measure %in% c("Density.of.foreground.pixels.in.hull.area",
#                          "Perimeter",
#                          "Width.of.bounding.rectangle",
#                          "Maximum.radius.from.hull.s.center.of.mass", 
#                          "Mean.radius.from.circle.s.center.of.mass", 
#                          "X..of.junction.voxels", 
#                          "X..of.quadruple.points"))
# stats_testing_features_filtered <- stats_morphologymeasures.animal(data = stats_testing_features_filtered, 
#                                                           model = "Value ~ Treatment*BrainRegion*Sex", type="lm",
#                                                           posthoc1 = "~Treatment|BrainRegion", 
#                                                           posthoc2 = "~Treatment|BrainRegion|Sex", adjust = "fdr")
# 
# View(stats_testing_features_filtered[[3]])
# write.csv(stats_testing_features[[3]], "/Users/alexlawson/Desktop/GliaData/final-data/stats-results/individualfeatures-all.csv", row.names = FALSE)
# 
# #plot for individual features in the pvn in males 
# stats_filtered <- stats_input_features %>%
#   filter(BrainRegion == "PVN" & Sex == "M")
# View(stats_filtered)
# # Create the plot without outliers
# ggplot(stats_filtered, aes(x = Measure, y = Value, fill = Treatment)) +
#   geom_boxplot(outlier.shape = NA) +  # Exclude outliers
#   theme_minimal() +
#   labs(
#     title = "Comparison of Cold vs. Control Individual Features - Males in PVN",
#     x = "Feature",
#     y = "Value"
#   ) +
#   facet_wrap(~ Measure, scales = "free") +
#   scale_fill_manual(values = c("COLD" = "lightblue", "CON" = "white"))
# 



