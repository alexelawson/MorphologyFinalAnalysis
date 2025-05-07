library(MicrogliaMorphologyR)
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

cluster_data <- read.csv("/Users/alexlawson/Desktop/Masters-Work/data-frames/cluster-percentage-data.csv", stringsAsFactors = FALSE)
data_frame_final <- read.csv("/Users/alexlawson/Desktop/Masters-Work/data-frames/processed-dataframe.csv")
stats_input_final <- cluster_data
plot <- clusterplots(data_frame_final, "PC1", "PC2")
plot
View(data_frame_final)

# make sure Cluster is a factor
data_frame_final$Cluster <- factor(data_frame_final$Cluster)

# create the plot
plot <- clusterplots(data_frame_final, "PC1", "PC2")

# add your custom hex‐colours
plot + 
  scale_color_manual(
    name   = "Cluster",
    values = c(
      "1" = "#FF8AB8",
      "2" = "#87CEEB",
      "3" = "#FFEB3B"
    )
  )

  
View(data_frame_final)
View(cluster_data)

View(cluster_data)
View(cluster_data)

cluster_data_PVN <- cluster_data %>% filter(BrainRegion=="PVN"&Cluster=="Ramified")
write_csv(cluster_data_PVN,"/Users/alexlawson/Desktop/Masters-Work/data-frames/cluster-percentage-data-PVN.csv")
stats_input_final$MouseID <- factor(stats_input_final$MouseID)
stats_input_final$Cluster <- factor(stats_input_final$Cluster)
stats_input_final$Treatment <- factor(stats_input_final$Treatment)
stats_input_final$Sex <- factor(stats_input_final$Sex)
View(stats_input_final)
# ---- Summary Data ----
counts_treatment <- cluster_data %>%
  group_by(Treatment) %>%
  summarise(Total_Count = sum(n))

View(counts_treatment)

counts_sex <- cluster_data %>%
  group_by(Sex) %>%
  summarise(Total_Count = sum(n))

counts_brainregion <- cluster_data %>%
  group_by(BrainRegion) %>%
  summarise(Total_Count = sum(n))

counts_treatment_sex <- cluster_data %>%
  group_by(Treatment, Sex) %>%
  summarise(Total_Count = sum(n))

counts_brainregion_sex <- cluster_data %>%
  group_by(BrainRegion, Sex) %>%
  summarise(Total_Count = sum(n))

counts_treatment_brainregion <- cluster_data %>%
  group_by(Treatment, BrainRegion) %>%
  summarise(Total_Count = sum(n))

counts_treatment_brainregion_sex <- cluster_data %>%
  group_by(Treatment, BrainRegion, Sex) %>%
  summarise(Total_Count = sum(n))

View(counts_treatment_brainregion_sex)
# ---- Stats ----   

stats_testing_combined <- stats_cluster.animal(stats_input_final,
                                          model = "percentage ~ Cluster*Treatment*BrainRegion*Sex + (1|MouseID)", 
                                          posthoc1 = "~Treatment|Cluster|BrainRegion", 
                                          posthoc2 = "~Treatment|Cluster|BrainRegion|Sex", adjust = "sidak")

stats_testing_combined[[2]]
stats_testing_combined[[3]]

stats_testing_males_pvn <- stats_cluster.animal(data = stats_input_final %>% filter(Sex=="M"&BrainRegion=="PVN"),
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

stats_testing_males_hypo <- stats_cluster.animal(data = stats_input_final %>% filter(Sex=="M"&BrainRegion=="HYPO"),
                                                   model = "percentage ~ Cluster*Treatment + (1|MouseID)", 
                                                   posthoc1 = "~Treatment|Cluster", 
                                                   posthoc2 = "~Treatment|Cluster", adjust = "sidak")
stats_testing_males_hypo[[2]]
stats_testing_males_hypo[[3]]


stats_testing_females_hypo <- stats_cluster.animal(data = stats_input_final %>% filter(Sex=="F"&BrainRegion=="HYPO"),
                                                  model = "percentage ~ Cluster*Treatment + (1|MouseID)", 
                                                  posthoc1 = "~Treatment|Cluster", 
                                                  posthoc2 = "~Treatment|Cluster", adjust = "sidak")
stats_testing_females_hypo[[2]]
stats_testing_females_hypo[[3]]

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

treatmentPlot_wfacet

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


#plot split by brain region - only controls 
brainRegionPlot_control <- cluster_data %>% 
  filter(Treatment=="CON") %>% #filter for control data only 
  ggplot(aes(x = BrainRegion, y = percentage)) +  # Use BrainRegion on x-axis
  geom_boxplot(aes(fill = BrainRegion), outlier.shape = NA) +  # Boxplots split by BrainRegion
  geom_point(aes(color = Sex), 
             position = position_jitter(width = 0.2, height = 0), size = 2, alpha = 0.6) +  # Points jittered within each group
  scale_fill_manual(values = c("ARC" = "white", "HYPO" = "lightgrey", "PVN" = "darkgrey")) +  # Colors for BrainRegion
  scale_color_manual(values = c("M" = "blue", "F" = "pink")) +  # Colors for Sex
  ggtitle("Comparison of Cluster Percentages in Control Mice") +
  labs(x = "Brain Region", y = "Percentage of Cluster", fill = "Brain Region", color = "Sex") +
  facet_grid(~Cluster, scales = "free_x") +  #Facet by Cluster (rows) and Treatment (columns)
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

brainRegionPlot_control




# ---- CBC Many ----
#procedure for creating multiple df's at once 
combinations <- data.frame(
  BrainRegion = c("HYPO"),
  Sex= c("M"),
  MouseID = c("MR138", "MR139", "MR140", "MR141", "MR156", "MR157", "MR158", "MR159", "MR161", "MR169", "MR172", "MR175", "MR196", "MR197", "MR198", "MR200", "MR273", "MR274", "MR277", "MR278", "MR284", "MR287", "MR290", "MR294"),
  SlideNumber = c("s04", "s05", "s06")
)

# Unique combinations of MouseID and SlideNumber for BrainRegion "ARC"
combinations <- data_frame_final %>%
  filter(BrainRegion == "HYPO"&Sex=="F") %>%
  distinct(MouseID, SlideNumber) %>%
  arrange(MouseID, SlideNumber)

# Unique combinations of MouseID and SlideNumber for BrainRegion "ARC"
combinations <- data_frame_final %>%
  filter(BrainRegion == "HYPO") %>%
  distinct(MouseID, SlideNumber) %>%
  arrange(MouseID, SlideNumber)

colorbycluster <- data_frame_final %>% 
  filter(BrainRegion=="HYPO", MouseID=="MR156", SlideNumber=="s02") %>% dplyr::select(c(Cluster, ID))
write.csv(colorbycluster, "/Users/alexlawson/Masters-Data-Final/Representative-images/Hypothalamus/02/csv/HYPO_MR156_CON_F_s02.csv")
head(colorbycluster)

View(data_frame_final)


View(data_frame_final)

# Output directory
output_directory <- "/Users/alexlawson/Desktop/GliaData/final-data/cbc-pvn"

library(dplyr)
library(stringr)

# Step 1: Get unique combinations (PVN as example, but can be changed)
combinations <- data_frame_final %>%
  filter(BrainRegion == "HYPO") %>%
  distinct(MouseID, SlideNumber, Treatment, Sex, BrainRegion)

# Step 2: Loop through each row
for (i in 1:nrow(combinations)) {
  mouse_id <- combinations$MouseID[i]
  slide_number <- combinations$SlideNumber[i]
  treatment <- combinations$Treatment[i]
  sex <- combinations$Sex[i]
  region <- combinations$BrainRegion[i]
  
  # Step 3: Filter original dataframe for this combo
  colorbycluster <- data_frame_final %>%
    filter(BrainRegion == region,
           MouseID == mouse_id,
           SlideNumber == slide_number) %>%
    dplyr::select(Cluster, ID)
  
  # Step 4: Construct file name
  filename <- paste0(region, "_", mouse_id, "_", treatment, "_", sex, "_", slide_number, ".csv")
  filepath <- file.path("/Users/alexlawson/Masters-Data-Final/Representative-images/Hypothalamus/dapionly/v5/csv", filename)
  
  # Step 5: Write to CSV
  write.csv(colorbycluster, filepath, row.names = FALSE)
}

# ---- Neuron Microglia Counting Stats SEE ANALYSIS 04 ----
interactions <- read.csv("/Users/alexlawson/Desktop/Masters-Work/data-frames/touch-percentages/Microglia-Neuron_Interaction_Percentages_by_Mouse_and_Treatment.csv")
View(interactions)
# Function to calculate cluster percentages within specified groups
calculate_cluster_percentage <- function(data, cluster_var, first_var = NULL, second_var = NULL, third_var=NULL) {
  
  # Dynamically determine grouping variables based on user input
  grouping_vars <- c()
  if (!is.null(first_var)) grouping_vars <- c(grouping_vars, first_var)
  if (!is.null(second_var)) grouping_vars <- c(grouping_vars, second_var)
  if (!is.null(third_var)) grouping_vars <- c(grouping_vars, third_var)
  
  if (length(grouping_vars) == 0) {
    stop("You must specify at least one grouping variable (brain_region or mouse_id).")
  }
  
  # Group by the selected variables
  data %>%
    group_by(across(all_of(grouping_vars))) %>%  # Dynamic grouping
    mutate(total_within_group = sum(n, na.rm = TRUE)) %>%  # Sum only within each group
    group_by(across(all_of(grouping_vars)), !!sym(cluster_var)) %>%
    summarise(
      cluster_count = sum(n, na.rm = TRUE),  # Sum counts for the cluster
      total_within_group = first(total_within_group),  # Keep the total count per group
      output = (cluster_count / total_within_group),  
      .groups = "drop"
    ) %>%
    arrange(across(all_of(grouping_vars)), desc(output))
}

# Function to adjust 0 and 1 values in a specified column
adjust_values <- function(df, column_name) {
  df[[column_name]] <- ifelse(df[[column_name]] == 0, df[[column_name]] + 0.0001, 
                              ifelse(df[[column_name]] == 1, df[[column_name]] - 0.0001, 
                                     df[[column_name]]))
  return(df)
}

cp_treatment_touchdegree_mouseid <- calculate_cluster_percentage(interactions, "Cluster", "Treatment", "TouchDegree", "MouseID")
cp_treatment_mouseid <- calculate_cluster_percentage(interactions, "Cluster", "Treatment", "MouseID")
cp_treatment_nocluster <- calculate_cluster_percentage(interactions, "Treatment", "MouseID")
cp_treatment_nocluster_touchdegree <- calculate_cluster_percentage(interactions, "Treatment", "MouseID", "TouchDegree")
stats_input_neurons<- adjust_values(cp_treatment_touchdegree_mouseid, "output")
write.csv(cp_treatment_mouseid, "/Users/alexlawson/Desktop/Masters-Work/microglia-neurons/stats/clustercounts.csv", row.names = FALSE)

View(cp_treatment_mouseid)
#Overall-Counts
t_test_all <- t.test(cluster_count ~ Treatment, 
                          data = cp_treatment_nocluster,
                          var.equal = FALSE)  # Welch's t-test (recommended)
t_test_all


#By Touch Degree
#One
t_test_td1 <- t.test(cluster_count ~ Treatment, 
                     data = cp_treatment_nocluster_touchdegree   %>% filter(TouchDegree=="1"),
                     var.equal = FALSE)  # Welch's t-test (recommended)
print(t_test_td1)
#Two
t_test_td2 <- t.test(cluster_count ~ Treatment, 
                     data = cp_treatment_nocluster_touchdegree   %>% filter(TouchDegree=="2"),
                     var.equal = FALSE)  # Welch's t-test (recommended)
print(t_test_td2)

anova_model <- aov(cluster_count ~ Treatment * TouchDegree, data = cp_treatment_nocluster_touchdegree)
summary(anova_model)
tukey_results <- TukeyHSD(anova_model)
print(tukey_results)



#By Cluster
#Ramified
t_test_ramified <- t.test(cluster_count ~ Treatment, 
                          data = cp_treatment_mouseid %>% filter(Cluster == "Ramified"),
                          var.equal = FALSE)  # Welch's t-test (recommended)

t_test_ameboid<- t.test(cluster_count ~ Treatment, 
                          data = cp_treatment_mouseid %>% filter(Cluster == "Ameboid"),
                          var.equal = FALSE)  # Welch's t-test (recommended)
print(t_test_ameboid)

t_test_rodlike<- t.test(cluster_count ~ Treatment, 
                        data = cp_treatment_mouseid %>% filter(Cluster == "Rod-Like"),
                        var.equal = FALSE)  # Welch's t-test (recommended)
print(t_test_rodlike)

#ANOVA by cluster

cp_treatment_mouseid$Treatment <- as.factor(cp_treatment_mouseid$Treatment)
cp_treatment_mouseid$Cluster <- as.factor(cp_treatment_mouseid$Cluster)
anova_model <- aov(cluster_count ~ Treatment * Cluster, data = cp_treatment_mouseid)
tukey_results <- TukeyHSD(anova_model)
print(tukey_results)
# Compute estimated marginal means (least-square means)
em_results <- emmeans(anova_model_m, pairwise ~ Treatment * Cluster, adjust = "sidak")
# Print results
em_results

# Convert Tukey results to a data frame
tukey_df <- do.call(rbind, lapply(tukey_results, as.data.frame))
# Add a column to indicate which comparison group the row belongs to
tukey_df$Comparison <- rownames(tukey_df)
rownames(tukey_df) <- NULL
# Print the dataframe
print(tukey_df)
# Keep only comparisons where the same cluster is being compared across treatments
tukey_filtered <- tukey_df[grep("Treatment:Cluster", tukey_df$Comparison), ]
# Remove comparisons that mix clusters (keep only comparisons where the cluster remains the same)# Print filtered results
print(tukey_filtered)
write.csv(tukey_filtered, "/Users/alexlawson/Desktop/Masters-Work/microglia-neurons/stats/2wayanova-clusters.csv", row.names = FALSE)


# Define clusters to analyze
touchdegree <- c("1", "2")

# Create an empty data frame for storing results
t_test_results <- data.frame(TouchDegree = character(),
                             Mean_COLD = numeric(),
                             Mean_CON = numeric(),
                             Mean_Diff = numeric(),
                             t_stat = numeric(),
                             df = numeric(),
                             p_value = numeric(),
                             CI_Lower = numeric(),
                             CI_Upper = numeric(),
                             stringsAsFactors = FALSE)

# Loop through each cluster and perform a t-test
for (cl in touchdegree) {
  subset_data <- cp_treatment_nocluster_touchdegree %>% filter(TouchDegree == cl)
  
  # Perform t-test
  t_test <- t.test(cluster_count ~ Treatment, data = subset_data, var.equal = FALSE)
  
  # Extract means for each group
  means <- t_test$estimate  # Named vector: mean in each treatment group
  mean_cold <- means["mean in group COLD"]
  mean_con <- means["mean in group CON"]
  
  # Store results in the data frame
  t_test_results <- rbind(t_test_results, data.frame(
    Cluster = cl,
    Mean_COLD = mean_cold,
    Mean_CON = mean_con,
    Mean_Diff = mean_con - mean_cold,  # Difference between CON & COLD
    t_stat = t_test$statistic,
    df = t_test$parameter,
    p_value = t_test$p.value,
    CI_Lower = t_test$conf.int[1],
    CI_Upper = t_test$conf.int[2]
  ))
}
# View the results
View(t_test_results)
write.csv(t_test_results, "/Users/alexlawson/Desktop/Masters-Work/microglia-neurons/stats/t_test_results_overall_touchdegree.csv", row.names = FALSE)

t_test_overall <- t.test(cluster_count ~ Treatment, data = cp_treatment_nocluster, var.equal = FALSE)
# Extract means for each group
means <- t_test_overall$estimate  # Named vector: mean in each treatment group
mean_cold <- means["mean in group COLD"]
mean_con <- means["mean in group CON"]

# Store the results in a data frame
t_test_results_overall <- data.frame(
  Mean_COLD = mean_cold,
  Mean_CON = mean_con,
  Mean_Diff = mean_con - mean_cold,  # Difference between CON & COLD
  t_stat = t_test_overall$statistic,
  df = t_test_overall$parameter,
  p_value = t_test_overall$p.value,
  CI_Lower = t_test_overall$conf.int[1],
  CI_Upper = t_test_overall$conf.int[2]
)

# View the results
print(t_test_results_overall)
write.csv(t_test_results_overall, "/Users/alexlawson/Desktop/Masters-Work/microglia-neurons/stats/t_test_results_overall.csv", row.names = FALSE)

#By Cluster + Touch Degree
#Ramified
t_test_td1_ramified <- t.test(cluster_count ~ Treatment, 
                     data = cp_treatment_touchdegree_mouseid %>% filter(Cluster=="Ramified"&TouchDegree=="1"),
                     var.equal = FALSE)  # Welch's t-test (recommended)
print(t_test_td1_ramified)

t_test_td2_ramified <- t.test(cluster_count ~ Treatment, 
                              data = cp_treatment_touchdegree_mouseid %>% filter(Cluster=="Ramified"&TouchDegree=="2"),
                              var.equal = FALSE)  # Welch's t-test (recommended)
print(t_test_td2_ramified)

#Rod-Like
t_test_td1_rodlike <- t.test(cluster_count ~ Treatment, 
                              data = cp_treatment_touchdegree_mouseid %>% filter(Cluster=="Rod-Like"&TouchDegree=="1"),
                              var.equal = FALSE)  # Welch's t-test (recommended)
print(t_test_td1_rodlike)

t_test_td2_rodlike <- t.test(cluster_count ~ Treatment, 
                             data = cp_treatment_touchdegree_mouseid %>% filter(Cluster=="Rod-Like"&TouchDegree=="2"),
                             var.equal = FALSE)  # Welch's t-test (recommended)
print(t_test_td2_rodlike)

#Ameboid
t_test_td1_ameboid <- t.test(cluster_count ~ Treatment, 
                             data = cp_treatment_touchdegree_mouseid %>% filter(Cluster=="Ameboid"&TouchDegree=="1"),
                             var.equal = FALSE)  # Welch's t-test (recommended)
print(t_test_td1_ameboid)

t_test_td2_ameboid <- t.test(cluster_count ~ Treatment, 
                             data = cp_treatment_touchdegree_mouseid %>% filter(Cluster=="Ameboid"&TouchDegree=="2"),
                             var.equal = FALSE)  # Welch's t-test (recommended)
print(t_test_td2_ameboid)


#ANOVA by cluster
cp_treatment_touchdegree_mouseid$Treatment <- as.factor(cp_treatment_touchdegree_mouseid$Treatment)
cp_treatment_touchdegree_mouseid$Cluster <- as.factor(cp_treatment_touchdegree_mouseid$Cluster)
cp_treatment_touchdegree_mouseid$TouchDegree <- as.factor(cp_treatment_touchdegree_mouseid$TouchDegree)
anova_model <- aov(cluster_count ~ Treatment * Cluster * TouchDegree, data = cp_treatment_touchdegree_mouseid)
summary(anova_model)
tukey_results <- TukeyHSD(anova_model)
print(tukey_results)
# Convert Tukey results into a dataframe
tukey_df <- do.call(rbind, lapply(tukey_results[7], as.data.frame))
print(tukey_df)
write.csv(tukey_df, "/Users/alexlawson/Desktop/Masters-Work/microglia-neurons/stats/2wayanova-clusters-td.csv")


# ---- Neuron Microglia Counting Plots----
#count with degree colored
touch_numbers_count <- cp_treatment_touchdegree_mouseid %>% 
  ggplot(aes(x = Cluster, y = cluster_count, group = interaction(Cluster, Treatment))) +
  geom_boxplot(aes(group = interaction(Cluster, Treatment), fill = Treatment), 
               outlier.shape = NA) +  # Remove outliers
  geom_jitter(aes(color = as.factor(TouchDegree)), 
              position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), 
              size = 2, alpha = 0.6) +  # Add jittered dots
  scale_fill_manual(values = c("COLD" = "white", "CON" = "lightgray")) +  # Custom box fill colors
  scale_color_manual(values = c("1" = "lightblue", "2" = "darkblue")) +  # Custom dot colors
  ggtitle("Number of Interactions in Cold and Control Groups in the PVN of Males: Colored by Degree") +
  labs(x = "Cluster", y = "Number of Interactions", fill = "Treatment", color = "Touch Degree") +
  theme_bw(base_size = 14)
touch_numbers_count

#count with degree faceted 
touch_numbers_faceted_count <- cp_treatment_touchdegree_mouseid %>% 
  ggplot(aes(x = Cluster, y = cluster_count, group = interaction(Cluster, Treatment, TouchDegree))) +
  geom_boxplot(aes(group = interaction(Cluster, Treatment), fill = Treatment), 
               outlier.shape = NA) +  # Remove outliers
  geom_jitter(aes(color = as.factor(TouchDegree)), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), 
              size = 2, alpha = 0.6) +  # Add jittered dots
  scale_fill_manual(values = c("COLD" = "white", "CON" = "lightgray")) +  # Custom box fill colors
  scale_color_manual(values = c("1" = "lightblue", "2" = "darkblue")) +  # Custom dot colors
  ggtitle("Number of Interactions in Cold and Control Groups in the PVN of Males: Separated by Degree") +
  labs(x = "Cluster", y = "Number of Interactions", fill = "Treatment", color = "Touch Degree") +
  theme_bw(base_size = 14) +
  facet_wrap(~TouchDegree)  # Separate by TouchDegree
touch_numbers_faceted_count

#count not separated by degree
touch_number_nodegree_count <- cp_treatment_mouseid %>% 
  ggplot(aes(x = Cluster, y = cluster_count, group = interaction(Cluster, Treatment))) +
  geom_boxplot(aes(group = interaction(Cluster, Treatment), fill = Treatment), 
               outlier.shape = NA) +  # Remove outliers
  geom_jitter(aes(color = Treatment),  # Remove TouchDegree, color by Treatment instead
              position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75), 
              size = 2, alpha = 0.8) +  # Add jittered dots
  scale_fill_manual(values = c("COLD" = "white", "CON" = "lightgray")) +  # Custom box fill colors
  scale_color_manual(values = c("COLD" = "lightblue", "CON" = "lightblue")) +  # Jitter colors match Treatment
  ggtitle("Number of Interactions in Cold and Control Groups in the PVN of Males") +
  labs(x = "Cluster", y = "Number of Interactions", fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 14)
touch_number_nodegree_count

#percentage with degree colored
touch_numbers <- cp_treatment_touchdegree_mouseid %>% 
  ggplot(aes(x = Cluster, y = output, group = interaction(Cluster, Treatment))) +
  geom_boxplot(aes(group = interaction(Cluster, Treatment), fill = Treatment), 
               outlier.shape = NA) +  # Remove outliers
  geom_jitter(aes(color = as.factor(TouchDegree)), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), 
              size = 2, alpha = 0.8) +  # Add jittered dots
  scale_fill_manual(values = c("COLD" = "white", "CON" = "lightgray")) +  # Custom box fill colors
  scale_color_manual(values = c("1" = "lightblue", "2" = "darkblue")) +  # Custom dot colors
  ggtitle("Comparison of Cluster Prevalence: Cold v. Control in the PVN of Males") +
  labs(x = "Cluster", y = "Proportion", fill = "Treatment", color = "Touch Degree") +
  theme_bw(base_size = 14)
touch_numbers

#percentage with degree faceted 
touch_numbers_faceted <- cp_treatment_touchdegree_mouseid %>% 
  ggplot(aes(x = Cluster, y = output, group = interaction(Cluster, Treatment, TouchDegree))) +
  geom_boxplot(aes(group = interaction(Cluster, Treatment), fill = Treatment), 
               outlier.shape = NA) +  # Remove outliers
  geom_jitter(aes(color = as.factor(TouchDegree)), 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75), 
              size = 2, alpha = 0.8) +  # Add jittered dots
  scale_fill_manual(values = c("COLD" = "white", "CON" = "lightgray")) +  # Custom box fill colors
  scale_color_manual(values = c("1" = "lightblue", "2" = "darkblue")) +  # Custom dot colors
  ggtitle("Comparison of Cluster Prevalence: Cold v. Control in the PVN of Males") +
  labs(x = "Cluster", y = "Percentage", fill = "Treatment", color = "Touch Degree") +
  theme_bw(base_size = 14) +
  facet_wrap(~TouchDegree)  # Separate by TouchDegree
touch_numbers_faceted

#percentage not separated by degree 
touch_number_nodegree <- cp_treatment_mouseid %>% 
  ggplot(aes(x = Cluster, y = output, group = interaction(Cluster, Treatment))) +
  geom_boxplot(aes(group = interaction(Cluster, Treatment), fill = Treatment), 
               outlier.shape = NA) +  # Remove outliers
  geom_jitter(aes(color = Treatment),  # Remove TouchDegree, color by Treatment instead
              position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75), 
              size = 2, alpha = 0.8) +  # Add jittered dots
  scale_fill_manual(values = c("COLD" = "white", "CON" = "lightgray")) +  # Custom box fill colors
  scale_color_manual(values = c("COLD" = "lightblue", "CON" = "lightblue")) +  # Jitter colors match Treatment
  ggtitle("Comparison of Cluster Prevalence: Cold v. Control in the PVN of Males") +
  labs(x = "Cluster", y = "percentage", fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 14)
touch_number_nodegree


# Create a boxplot instead of a bar plot
touch_number_plot <- ggplot(cp_treatment_nocluster, aes(x = Treatment, y = cluster_count, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +  # Transparent boxplot with no outliers shown
  geom_jitter(aes(color = Treatment), width = 0.15, size = 2, alpha = 0.8) +  # Overlay jittered data points
  scale_fill_manual(values = c("COLD" = "white", "CON" = "lightgray")) +  # Boxplot fill colors
  scale_color_manual(values = c("COLD" = "lightblue", "CON" = "lightblue")) +  # Jitter point colors
  ggtitle("Comparison of Total Number of Interactions in Cold and Control Groups in the PVN of Males") +
  labs(x = "Treatment", y = "Count", fill = "Treatment", color = "Treatment") +
  theme_bw(base_size = 14)

# Display the plot
touch_number_plot

# ---- Old Stats ----

stats_input_neurons$Cluster <- factor(stats_input_neurons$Cluster)
stats_input_neurons$Treatment <- factor(stats_input_neurons$Treatment)
#stats_input_neurons$TouchDegree <- factor(stats_input_neurons$TouchDegree)

View(cp_treatment_mouseid)
stats_testing_neurons <- stats_cluster.animal(stats_input_neurons,
                                              model = "output ~ Cluster*Treatment*TouchDegree + (1|MouseID)", 
                                              posthoc1 = "~Treatment|Cluster", 
                                              posthoc2 = "~Treatment|Cluster|TouchDegree", adjust = "sidak")

stats_testing_neurons[[2]]
stats_testing_neurons[[3]]

View(cp_treatment_nocluster)
# Model for total interactions (not separated by touch degree)
model_lm <- lm(total_within_group ~ Treatment, data = cp_treatment_nocluster)
summary(model_lm)

#model_poisson <- glm(total_within_group ~ Treatment, 
#                     family = poisson(link = "log"), 
#                     data = cp_treatment_nocluster)
#summary(model_poisson)
#dispersion_ratio <- sum(residuals(model_poisson, type = "pearson")^2) / df.residual(model_poisson)
#dispersion_ratio
model_nb <- glm.nb(total_within_group ~ Treatment, data = cp_treatment_nocluster)
#summary(model_nb)
#AIC(model_poisson, model_nb)
#summary(model_nb)
# Get estimated marginal means (EMMs) for Treatment
emmeans_model <- emmeans(model_nb, pairwise ~ Treatment, adjust = "sidak")
# Print results
summary(emmeans_model)
# Fit the Negative Binomial model with interaction
model_nb_degree <- glm.nb(total_within_group ~ Treatment * TouchDegree, data = cp_treatment_nocluster_touchdegree)
# View model summary
summary(model_nb_degree)
emmeans_model_degree <- emmeans(model_nb_degree, pairwise ~ Treatment | TouchDegree, adjust = "sidak")
summary(emmeans_model_degree)

model_nb_cluster <- glm.nb(cluster_count ~ Treatment * Cluster, data = cp_treatment_mouseid)
summary(model_nb_cluster)
emmeans_model_cluster <- emmeans(model_nb_cluster, pairwise ~ Treatment | Cluster, adjust = "sidak")
summary(emmeans_model_cluster)

model_nb_cluster_td <- glm.nb(cluster_count ~ Treatment * Cluster, data = cp_treatment_touchdegree_mouseid %>% filter(TouchDegree=="2"))
emmeans_model_cluster_td <- emmeans(model_nb_cluster_td, pairwise ~ Treatment | Cluster, adjust = "sidak")
summary(emmeans_model_cluster_td)

model_nb_cluster_td <- glm.nb(cluster_count ~ Treatment * Cluster, data = cp_treatment_touchdegree_mouseid %>% filter(TouchDegree=="1"))
emmeans_model_cluster_td <- emmeans(model_nb_cluster_td, pairwise ~ Treatment | Cluster, adjust = "sidak")
summary(emmeans_model_cluster_td)

View(cp_treatment_nocluster_touchdegree)






# ---- FEMALES ----
# ---- Stats ----
interactions_females <- read_excel("/Users/alexlawson/Desktop/Masters-Work/microglia-neurons/females-stats.xlsx")
cp_treatment_touchdegree_mouseid_f <- calculate_cluster_percentage(interactions_females, "Cluster", "Treatment", "TouchDegree", "MouseID")
cp_treatment_mouseid_f <- calculate_cluster_percentage(interactions_females, "Cluster", "Treatment", "MouseID")
cp_treatment_nocluster_f <- calculate_cluster_percentage(interactions_females, "Treatment", "MouseID")
cp_treatment_nocluster_touchdegree_f <- calculate_cluster_percentage(interactions_females, "Treatment", "MouseID", "TouchDegree")
stats_input_neurons_f <- adjust_values(cp_treatment_touchdegree_mouseid_f, "output")

t_test_all_f <- t.test(cluster_count ~ Treatment, 
                     data = cp_treatment_nocluster_f,
                     var.equal = FALSE)  # Welch's t-test (recommended)
t_test_all_f

#By Touch Degree
#One
t_test_td1_f <- t.test(cluster_count ~ Treatment, 
                     data = cp_treatment_nocluster_touchdegree_f   %>% filter(TouchDegree=="1"),
                     var.equal = FALSE)  # Welch's t-test (recommended)
print(t_test_td1_f)
#Two
t_test_td2_f <- t.test(cluster_count ~ Treatment, 
                     data = cp_treatment_nocluster_touchdegree_f   %>% filter(TouchDegree=="2"),
                     var.equal = FALSE)  # Welch's t-test (recommended)
print(t_test_td2_f)

#anova by cluster 
cp_treatment_mouseid_f$Treatment <- as.factor(cp_treatment_mouseid_f$Treatment)
cp_treatment_mouseid_f$Cluster <- as.factor(cp_treatment_mouseid_f$Cluster)
anova_model_f <- aov(cluster_count ~ Treatment * Cluster, data = cp_treatment_mouseid_f)
summary(anova_model_f)
tukey_results_f <- TukeyHSD(anova_model_f)
print(tukey_results_f)

#anova by cluster and touchdegree
#ANOVA by cluster
cp_treatment_touchdegree_mouseid_f$Treatment <- as.factor(cp_treatment_touchdegree_mouseid_f$Treatment)
cp_treatment_touchdegree_mouseid_f$Cluster <- as.factor(cp_treatment_touchdegree_mouseid_f$Cluster)
cp_treatment_touchdegree_mouseid_f$TouchDegree <- as.factor(cp_treatment_touchdegree_mouseid_f$TouchDegree)
anova_model_td_f <- aov(cluster_count ~ Treatment * Cluster * TouchDegree, data = cp_treatment_touchdegree_mouseid_f)
summary(anova_model_td_f)
tukey_results_td_f <- TukeyHSD(anova_model_td_f)
print(tukey_results_td_f)






