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
View(raw_data01)

raw_data01_log <- transform_log(raw_data01,1, start=8, end=34) 
pca_data01_log <- pcadata(raw_data01_log, featurestart=8, featureend=34,
                          pc.start=1, pc.end=10)
pca_data01_log_transformed <- transform_scale(pca_data01_log, start=1, end=3) # scale pca data as input for k-means clustering
kmeans_input01_log <- pca_data01_log_transformed[1:3]
data_kmeans03_log <- kmeans(kmeans_input01_log, centers=3)
pca_kmeans03_log <- cbind(pca_data01_log[1:3], raw_data01, as.data.frame(data_kmeans03_log$cluster)) %>%
  rename(Cluster=`data_kmeans03_log$cluster`) 
clusterfeatures(pca_kmeans03_log, featurestart=11, featureend=37)
#1 - Ameboid
#2 - Rod-Like
#3 - Ramified

plot_log_3 <- clusterplots(pca_kmeans03_log, "PC1", "PC2")
cp_log_3 <- clusterpercentage(pca_kmeans03_log, "Cluster", Treatment, Sex, BrainRegion, MouseID)
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

# ---- PLOTS ----     
#cold vs. control cluster percentages: split by brainregion and sex
cp_log_3 %>% 
  ggplot(aes(x = Cluster, y = percentage, group = interaction(Cluster, Treatment))) +
  facet_grid(Sex ~ BrainRegion) + 
  geom_boxplot(aes(group = interaction(Cluster, Treatment), fill = Treatment), outlier.shape = NA) +  # Remove outliers
  scale_fill_manual(values = c("lightblue", "white")) +
  ggtitle("Cold vs. Control Cluster Percentages Split by Brain Region and Sex") +
  labs(x = "Cluster", y = "Percentage", fill = "Treatment") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

#control male vs. female in the different brain regions, with dots
cp_log_3 %>%
  filter(Treatment == "CON", Sex %in% c("M", "F")) %>%  # Filter for Control Males and Females
  ggplot(aes(x = Cluster, y = percentage, group = interaction(Cluster, Sex))) +
  facet_grid(. ~ BrainRegion) +  # Columns = Brain Region
  geom_boxplot(aes(fill = Sex), position = position_dodge(width = 0.8)) +  # Compare sexes
  scale_fill_manual(values = c("lightpink", "lightblue")) +
  geom_point(position = position_dodge(width = 0.8), size = 0.75) +  # Add points for clarity
  ggtitle("Mouse Dataset: Control Males vs Females Across Brain Regions") +
  labs(fill = "Sex") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

#control male vs. female in the different brain regions, without dots
cp_log_3 %>%
  filter(Treatment == "CON", Sex %in% c("M", "F")) %>%
  ggplot(aes(x = Cluster, y = percentage, group = interaction(Cluster, Sex))) +
  facet_wrap(~ BrainRegion, nrow = 2) +  # Adjust facets for better layout
  geom_boxplot(aes(fill = Sex), position = position_dodge(width = 0.8), outlier.shape = NA) +  # Remove outliers
  # Adjust point position
  scale_fill_manual(values = c("lightpink", "lightblue")) +
  ggtitle("Mouse Dataset: Control Males vs Females Across Brain Regions") +
  labs(fill = "Sex") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

#cold male vs. female in the different brain regions, without dots
cp_log_3 %>%
  filter(Treatment == "COLD", Sex %in% c("M", "F")) %>%
  ggplot(aes(x = Cluster, y = percentage, group = interaction(Cluster, Sex))) +
  facet_wrap(~ BrainRegion, nrow = 2) +  # Adjust facets for better layout
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
cp_log_3_withslides %>%
  filter(BrainRegion == "HYPO", SlideNumber %in% c("s01", "s02", "s03", "s04", "s05", "s06")) %>%
  ggplot(aes(x = Cluster, y = percentage, group = interaction(Cluster, Treatment))) +
  facet_grid(Sex ~ SlideNumber) +
  geom_boxplot(aes(group = interaction(Cluster, Treatment), fill = Treatment), outlier.shape = NA) +
  scale_fill_manual(values = c("skyblue", "lightgreen")) +
  ggtitle("Mouse Dataset: Cold vs. Control in the Hypothalamus, Separated by Sex and Section (s01-s06)") +
  labs(fill = "Treatment") +
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
  scale_fill_manual(values = c("skyblue", "lightgreen")) +
  ggtitle("Mouse Dataset: Cold vs. Control in the Hypothalamus, Separated by Sex and Section (s07-s12)") +
  labs(fill = "Treatment") +
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
cp_log_3_withslides %>%
  filter(BrainRegion == "HYPO", Treatment == "CON", Sex %in% c("M", "F")) %>%
  ggplot(aes(x = Cluster, y = percentage, group = interaction(Cluster, Sex))) +
  facet_wrap(~ SlideNumber, nrow = 2, scales = "free_x") + 
  geom_boxplot(aes(fill = Sex), outlier.shape = NA) +
  scale_fill_manual(values = c("lightpink", "lightblue")) +  # Different colors for Sex
  ggtitle("Mouse Dataset: Males vs Females in Controls in Hypothalamus, Separated by Section") +
  labs(fill = "Sex") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

# Create the combined plot
cp_log_3_withslides %>%
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



# ---- Stats ----   
stats_input <- cp_log_3 
stats_input$MouseID <- factor(stats_input$MouseID)
stats_input$Cluster <- factor(stats_input$Cluster)
stats_input$Treatment <- factor(stats_input$Treatment)
View(stats_input)

stats_testing_control <- stats_cluster.animal(data = stats_input %>% filter(Treatment=="CON"),
                                              model = "percentage ~ Cluster*Sex*BrainRegion + (1|MouseID)", 
                                              posthoc1 = "~Sex|Cluster", 
                                              posthoc2 = "~Sex|Cluster|BrainRegion", adjust = "fdr")
stats_testing_control[[1]]
stats_testing_control[[2]]
stats_testing_control[[3]]
write.csv(stats_testing_control[[3]], "/Users/alexlawson/Desktop/GliaData/final-data/stats-results/sex-cluster-brainregion-control", row.names = FALSE)

stats_testing_cold <- stats_cluster.animal(data = stats_input %>% filter(Treatment=="COLD"),
                                           model = "percentage ~ Cluster*BrainRegion*Sex + (1|MouseID)", 
                                           posthoc1 = "~Sex|Cluster", 
                                           posthoc2 = "~Sex|Cluster|BrainRegion", adjust="fdr")
stats_testing_cold[[1]]
stats_testing_cold[[2]]
stats_testing_cold[[3]]

stats_testing_all <- stats_cluster.animal(data = stats_input,
                                          model = "percentage ~ Cluster*Treatment*Sex*BrainRegion + (1|MouseID)", 
                                          posthoc1 = "~Treatment|Cluster|BrainRegion", 
                                          posthoc2 = "~Treatment|Cluster|Sex|BrainRegion", adjust = "bonferroni")

stats_testing_all[[1]]
stats_testing_all[[2]]
stats_testing_all[[3]]
write.csv(stats_testing_all[[2]], "/Users/alexlawson/Desktop/GliaData/final-data/stats-results/treatment-cluster-brainregion-notsexsep.csv", row.names = FALSE)

stats_testing_males <- stats_cluster.animal(data = stats_input %>% filter(Sex=="M"), 
                                            model = "percentage ~ Cluster*Treatment*BrainRegion + (1|MouseID)", 
                                            posthoc1 = "~Treatment|Cluster", 
                                            posthoc2 = "~Treatment|Cluster|BrainRegion", adjust = "bonferroni")

stats_testing_males[[1]]
stats_testing_males[[2]]
stats_testing_males[[3]]



stats_testing_fem <- stats_cluster.animal(data = stats_input %>% filter(Sex=="F"), 
                                          model = "percentage ~ Cluster*Treatment*BrainRegion + (1|MouseID)", 
                                          posthoc1 = "~Treatment|Cluster", 
                                          posthoc2 = "~Treatment|Cluster|BrainRegion", adjust = "fdr")

stats_testing_fem[[1]]
stats_testing_fem[[2]]
stats_testing_fem[[3]]

stats_input_pvn <- stats_input %>% filter(BrainRegion=="PVN")
stats_testing_males_pvn <- stats_cluster.animal(data = stats_input_pvn,
                                               model = "percentage ~ Cluster*Treatment + (1|MouseID)", 
                                               posthoc1 = "~Treatment|Cluster", 
                                               posthoc2 = "~Treatment|Cluster", adjust = "holm")

stats_testing_males_pvn[[1]]
stats_testing_males_pvn[[2]]
stats_testing_males_pvn[[3]]



stats_input_slides <- cp_log_3_withslides 
stats_input_slides$MouseID <- factor(stats_input_slides$MouseID)
stats_input_slides$Cluster <- factor(stats_input_slides$Cluster)
stats_input_slides$Treatment <- factor(stats_input_slides$Treatment)
stats_input_slides$SlideNumber <- factor(stats_input_slides$SlideNumber)
View(stats_input_slides)

stats_input_slides_hypothalamus <- stats_input_slides %>% filter(BrainRegion=="HYPO")
stats_testing_hypothalamus <- stats_cluster.animal(data = stats_input_slides_hypothalamus,
                                                   model = "percentage ~ Cluster*Treatment*SlideNumber*Sex", 
                                                   posthoc1 = "~Sex|Cluster|Treatment|SlideNumber", 
                                                   posthoc2 = "~Treatment|Cluster|SlideNumber|Sex")

stats_testing_hypothalamus[[1]]
stats_testing_hypothalamus[[2]]
stats_testing_hypothalamus[[3]]

write.csv(stats_testing_hypothalamus[[2]], "/Users/alexlawson/Desktop/GliaData/final-data/stats-results/treatment-cluster-sex-hypothalamus.csv", row.names = FALSE)


stats_testing_hypothalamus_control <- stats_cluster.animal(data = stats_input_slides_hypothalamus %>%filter(Treatment=="CON"),  
                                                   model = "percentage ~ Cluster*SlideNumber*Sex", 
                                                   posthoc1 = "~Sex|Cluster|SlideNumber", 
                                                   posthoc2 = "~Cluster|SlideNumber")

stats_testing_hypothalamus_control[[1]]
stats_testing_hypothalamus_control[[2]]
stats_testing_hypothalamus_control[[3]]

stats_testing_hypothalamus_cold <- stats_cluster.animal(data = stats_input_slides_hypothalamus %>%filter(Treatment=="COLD"),  
                                                           model = "percentage ~ Cluster*SlideNumber*Sex", 
                                                           posthoc1 = "~Sex|Cluster|SlideNumber", 
                                                           posthoc2 = "~Sex|Cluster")

stats_testing_hypothalamus_cold[[1]]
stats_testing_hypothalamus_cold[[2]]
stats_testing_hypothalamus_cold[[3]]

s03_hypo <-stats_input_slides_hypothalamus %>% filter(SlideNumber=="s03")
stats_testing_hypothalamus_s03 <- stats_cluster.animal(data = s03_hypo,  
                                                   model = "percentage ~ Cluster*Treatment*Sex", 
                                                   posthoc1 = "~Treatment|Cluster", 
                                                   posthoc2 = "~Treatment|Cluster|Sex")

stats_testing_hypothalamus_s03[[2]]
stats_testing_hypothalamus_s03[[3]]

stats_input_features <- raw_data01 %>% 
  group_by(MouseID, Sex, Treatment, BrainRegion) %>% 
  summarise(across("Foreground.pixels":"Maximum.branch.length", ~mean(.x))) %>% 
  gather(Measure, Value, "Foreground.pixels":"Maximum.branch.length")


stats_input_features$Sex <- factor(stats_input_features$Sex)
stats_input_features$BrainRegion <- factor(stats_input_features$BrainRegion)
stats_testing_features <- stats_morphologymeasures.animal(data = stats_input_features, 
                                                 model = "Value ~ Treatment*BrainRegion*Sex", type="lm",
                                                 posthoc1 = "~Treatment|BrainRegion", 
                                                 posthoc2 = "~Treatment|BrainRegion|Sex", "bonferroni")
stats_testing_features[[2]]
stats_testing_features[[3]]


stats_input_features_filtered <- stats_input_features %>%
  filter(Measure %in% c("Density.of.foreground.pixels.in.hull.area",
                         "Perimeter",
                         "Width.of.bounding.rectangle",
                         "Maximum.radius.from.hull.s.center.of.mass", 
                         "Mean.radius.from.circle.s.center.of.mass", 
                         "X..of.junction.voxels", 
                         "X..of.quadruple.points"))
stats_testing_features_filtered <- stats_morphologymeasures.animal(data = stats_testing_features_filtered, 
                                                          model = "Value ~ Treatment*BrainRegion*Sex", type="lm",
                                                          posthoc1 = "~Treatment|BrainRegion", 
                                                          posthoc2 = "~Treatment|BrainRegion|Sex", adjust = "fdr")

View(stats_testing_features_filtered[[3]])
write.csv(stats_testing_features[[3]], "/Users/alexlawson/Desktop/GliaData/final-data/stats-results/individualfeatures-all.csv", row.names = FALSE)

#plot for individual features in the pvn in males 
stats_filtered <- stats_input_features %>%
  filter(BrainRegion == "PVN" & Sex == "M")
View(stats_filtered)
# Create the plot without outliers
ggplot(stats_filtered, aes(x = Measure, y = Value, fill = Treatment)) +
  geom_boxplot(outlier.shape = NA) +  # Exclude outliers
  theme_minimal() +
  labs(
    title = "Comparison of Cold vs. Control Individual Features - Males in PVN",
    x = "Feature",
    y = "Value"
  ) +
  facet_wrap(~ Measure, scales = "free") +
  scale_fill_manual(values = c("COLD" = "lightblue", "CON" = "white"))



# ---- Color By Cluster ----   
#females s03/s04

colorbycluster_log_3clusters_arc_2 <- pca_kmeans03_log %>% 
  filter(BrainRegion=="HYPO",MouseID=="MR175", SlideNumber=="s04") %>% select(c(Cluster, ID))
write.csv(colorbycluster_log_3clusters_arc_2, "/Users/alexlawson/Desktop/GliaData/final-data/cbc-hypo-fem-s03s04/cbc-csv/conmr175s04.csv", row.names = FALSE)


