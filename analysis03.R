library(MicrogliaMorphologyR)
library(factoextra)
library(ppclust)
library(dplyr)
library(ggsignif)
library(ggrepel)


cluster_data <- read.csv("/Users/alexlawson/Desktop/Masters-Work/data-frames/cluster-percentage-data.csv", stringsAsFactors = FALSE)
cluster_data$Treatment <- factor(cluster_data$Treatment, levels = c("CON", "COLD"))
stats_input_final <- cluster_data
View(cluster_data)
stats_input_final$MouseID <- factor(stats_input_final$MouseID)
stats_input_final$Treatment <- factor(stats_input_final$Treatment)
View(stats_input_final)
# ---- Stats ----   

stats_testing_combined <- stats_cluster.animal(data = stats_input_final %>% filter(Treatment=="CON"&BrainRegion=="HYPO"),
                                          model = "percentage ~ Cluster*Sex + (1|MouseID)", 
                                          posthoc1 = "~Sex|Cluster", 
                                          posthoc2 = "~Sex|Cluster|Sex", adjust = "sidak")

stats_testing_combined[[2]]
stats_testing_combined[[3]]

stats_testing_males_pvn <- stats_cluster.animal(data = stats_input_final %>% filter(Sex=="M"&BrainRegion==""),
                                        model = "percentage ~ Cluster*Treatment + (1|MouseID)", 
                                        posthoc1 = "~Treatment|Cluster", 
                                        posthoc2 = "~Treatment|Cluster", adjust = "sidak")
stats_testing_males_pvn[[2]]
stats_testing_males_pvn[[3]]


stats_testing_females_pvn <- stats_cluster.animal(data = stats_input_final %>% filter(Sex=="F"&BrainRegion=="PVN"),
                                                model = "percentage ~ Cluster*Treatment + (1|MouseID)", 
                                                posthoc1 = "~Treatment|Cluster", 
                                                posthoc2 = "~Treatment|Cluster", adjust = "sidak")
stats_testing_females_pvn[[2]]
stats_testing_females_pvn[[3]]

stats_testing_brain_comparison_all<- stats_cluster.animal(data = stats_input_final,
                                                      model = "percentage ~ Cluster*BrainRegion*Sex + (1|MouseID)", 
                                                      posthoc1 = "~BrainRegion|Cluster", 
                                                      posthoc2 = "~BrainRegion|Cluster|Sex", adjust="sidak")


stats_testing_brain_comparison_all[[2]]
stats_testing_brain_comparison_all[[3]]

stats_testing_brain_comparison_control<- stats_cluster.animal(data = stats_input_final %>% filter(Treatment=="CON"),
                                                          model = "percentage ~ Cluster*BrainRegion + (1|MouseID)", 
                                                          posthoc1 = "~BrainRegion|Cluster", 
                                                          posthoc2 = "~BrainRegion|Cluster", adjust="sidak")


stats_testing_brain_comparison_control[[2]]
stats_testing_brain_comparison_control[[3]]

# ---- Plots ----   
#Cold vs Control, dots colored by sex
coldvcontrol_all <- cluster_data %>% 
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
#Cold vs Control faceted by sex
coldvcontrol <- cluster_data %>% 
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

#control m vs. female, with dots 
cluster_data %>%
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
control_braincomparison <- cluster_data %>%
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

#plot split by brain region comparing clusers - faceted 
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


#plot split by brain region and treatment - faceted 
treatmentPlot_wfacet <- cluster_data %>% 
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

#plot split by brain region and treatment - faceted 
brainRegionPlot_treatment <- cluster_data %>% 
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

#plot split by brain region and treatment - males only 
brainRegionPlot_treatment_males <- cluster_data %>%
  filter(Sex == "M")%>%  # Filter for males only
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

#plot split by brain region and treatment - females only  
brainRegionPlot_treatment_females <- cluster_data %>%
  filter(Sex == "F") %>%  # Filter for females only
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


