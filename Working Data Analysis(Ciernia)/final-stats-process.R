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


#Loading in the data
#Cluster Percentage Data (not separated by section)
cluster_data <- read.csv("/Users/alexlawson/Desktop/Masters-Work/data-frames/cluster-percentage-data.csv", stringsAsFactors = FALSE)
#Cluster Percentage Data (separated by section - sections only used for the hypothalamus stats)
cluster_data_sections <- read.csv("/Users/alexlawson/Desktop/Masters-Work/data-frames/cluster-percentage-data-w-sections.csv", stringsAsFactors = FALSE)
#Loading in the raw data frame (contains all individual morphological features as well as cluster data and first 2 PCAs)
data_frame_final <- read.csv("/Users/alexlawson/Desktop/Masters-Work/data-frames/final-dataframe.csv")


#preparing the data for downstream processing
stats_input_final <- cluster_data
stats_input_final$MouseID <- factor(stats_input_final$MouseID)
stats_input_final$Cluster <- factor(stats_input_final$Cluster)
stats_input_final$Treatment <- factor(stats_input_final$Treatment)
stats_input_final$Sex <- factor(stats_input_final$Sex)


#All data combined
stats_testing_combined <- stats_cluster.animal(stats_input_final,
                                               model = "percentage ~ Cluster*Treatment*BrainRegion*Sex + (1|MouseID)", 
                                               posthoc1 = "~Treatment|Cluster|BrainRegion", 
                                               posthoc2 = "~Treatment|Cluster|BrainRegion|Sex", adjust = "sidak")

stats_testing_combined[[2]]
stats_testing_combined[[3]]


#Just males in the PVN (because all data showed significance here)
stats_testing_males_pvn <- stats_cluster.animal(data = stats_input_final %>% filter(Sex=="M"&BrainRegion=="PVN"),
                                                model = "percentage ~ Cluster*Treatment + (1|MouseID)", 
                                                posthoc1 = "~Treatment|Cluster", 
                                                posthoc2 = "~Treatment|Cluster", adjust = "sidak")
stats_testing_males_pvn[[2]]
stats_testing_males_pvn[[3]]


#Just females in the PVN
stats_testing_females_pvn <- stats_cluster.animal(data = stats_input_final %>% filter(Sex=="F"&BrainRegion=="PVN"),
                                                  model = "percentage ~ Cluster*Treatment + (1|MouseID)", 
                                                  posthoc1 = "~Treatment|Cluster", 
                                                  posthoc2 = "~Treatment|Cluster", adjust = "sidak")
stats_testing_females_pvn[[2]]
stats_testing_females_pvn[[3]]

#Just males in hypothalamus
stats_testing_males_hypo <- stats_cluster.animal(data = stats_input_final %>% filter(Sex=="M"&BrainRegion=="HYPO"),
                                                 model = "percentage ~ Cluster*Treatment + (1|MouseID)", 
                                                 posthoc1 = "~Treatment|Cluster", 
                                                 posthoc2 = "~Treatment|Cluster", adjust = "sidak")
stats_testing_males_hypo[[2]]
stats_testing_males_hypo[[3]]

#Just females in hypothalamus
stats_testing_females_hypo <- stats_cluster.animal(data = stats_input_final %>% filter(Sex=="F"&BrainRegion=="HYPO"),
                                                   model = "percentage ~ Cluster*Treatment + (1|MouseID)", 
                                                   posthoc1 = "~Treatment|Cluster", 
                                                   posthoc2 = "~Treatment|Cluster", adjust = "sidak")
stats_testing_females_hypo[[2]]
stats_testing_females_hypo[[3]]


#Just males in the ARC
stats_testing_males_arc <- stats_cluster.animal(data = stats_input_final %>% filter(Sex=="M"&BrainRegion=="ARC"),
                                                   model = "percentage ~ Cluster*Treatment + (1|MouseID)", 
                                                   posthoc1 = "~Treatment|Cluster", 
                                                   posthoc2 = "~Treatment|Cluster", adjust = "sidak")
stats_testing_males_arc[[2]]
stats_testing_males_arc[[3]]

#Just females in the ARC
stats_testing_females_arc <- stats_cluster.animal(data = stats_input_final %>% filter(Sex=="F"&BrainRegion=="ARC"),
                                                model = "percentage ~ Cluster*Treatment + (1|MouseID)", 
                                                posthoc1 = "~Treatment|Cluster", 
                                                posthoc2 = "~Treatment|Cluster", adjust = "sidak")
stats_testing_females_arc[[2]]
stats_testing_females_arc[[3]]


#Comparison of brain regions 
stats_testing_brain_comparison_all<- stats_cluster.animal(data = stats_input_final,
                                                          model = "percentage ~ Cluster*BrainRegion*Sex + (1|MouseID)", 
                                                          posthoc1 = "~BrainRegion|Cluster", 
                                                          posthoc2 = "~BrainRegion|Cluster|Sex", adjust="sidak")


stats_testing_brain_comparison_all[[2]]
stats_testing_brain_comparison_all[[3]]

#Comparison of control and cold brain regions
stats_testing_brain_comparison_control<- stats_cluster.animal(data = stats_input_final %>% filter(Treatment=="CON"),
                                                              model = "percentage ~ Cluster*BrainRegion + (1|MouseID)", 
                                                              posthoc1 = "~BrainRegion|Cluster", 
                                                              posthoc2 = "~BrainRegion|Cluster", adjust="sidak")
stats_testing_brain_comparison_control[[2]]
stats_testing_brain_comparison_control[[3]]


#Preparing data separated by section for downstream processing, filtering out the hypothalamus
stats_input_with_sections <- cluster_data_sections %>% filter(BrainRegion=="HYPO")
stats_input_with_sections$MouseID <- factor(stats_input_with_sections$MouseID)
stats_input_with_sections$Cluster <- factor(stats_input_with_sections$Cluster)
stats_input_with_sections$Treatment <- factor(stats_input_with_sections$Treatment)
stats_input_with_sections$Sex <- factor(stats_input_with_sections$Sex)

#Overall comparison
stats_testing_hypothalamus <- stats_cluster.animal(data = stats_input_with_sections,
                                                   model = "percentage ~ Cluster*Treatment*SlideNumber*Sex + (1|MouseID)", 
                                                   posthoc1 = "~Treatment|Cluster|SlideNumber", 
                                                   posthoc2 = "~Treatment|Cluster|SlideNumber|Sex", adjust="sidak")

stats_testing_hypothalamus[[2]]
stats_testing_hypothalamus[[3]]

#Stats section by section to assess any potential differences
stats_testing_hypothalamus_s01 <- stats_cluster.animal(data = stats_input_with_sections%>% filter(SlideNumber=="s01"),
                                                   model = "percentage ~ Cluster*Treatment*Sex + (1|MouseID)", 
                                                   posthoc1 = "~Treatment|Cluster", 
                                                   posthoc2 = "~Treatment|Cluster|Sex", adjust="sidak")

stats_testing_hypothalamus_s01[[2]]
stats_testing_hypothalamus_s01[[3]]

stats_testing_hypothalamus_s02 <- stats_cluster.animal(data = stats_input_with_sections%>% filter(SlideNumber=="s02"),
                                                       model = "percentage ~ Cluster*Treatment*Sex + (1|MouseID)", 
                                                       posthoc1 = "~Treatment|Cluster", 
                                                       posthoc2 = "~Treatment|Cluster|Sex", adjust="sidak")

stats_testing_hypothalamus_s02[[2]]
stats_testing_hypothalamus_s02[[3]]

stats_testing_hypothalamus_s03 <- stats_cluster.animal(data = stats_input_with_sections%>% filter(SlideNumber=="s03"),
                                                       model = "percentage ~ Cluster*Treatment*Sex + (1|MouseID)", 
                                                       posthoc1 = "~Treatment|Cluster", 
                                                       posthoc2 = "~Treatment|Cluster|Sex", adjust="sidak")

stats_testing_hypothalamus_s03[[2]]
stats_testing_hypothalamus_s03[[3]]


stats_testing_hypothalamus_s04 <- stats_cluster.animal(data = stats_input_with_sections%>% filter(SlideNumber=="s04"),
                                                       model = "percentage ~ Cluster*Treatment*Sex + (1|MouseID)", 
                                                       posthoc1 = "~Treatment|Cluster", 
                                                       posthoc2 = "~Treatment|Cluster|Sex", adjust="sidak")

stats_testing_hypothalamus_s04[[2]]
stats_testing_hypothalamus_s04[[3]]

stats_testing_hypothalamus_s05 <- stats_cluster.animal(data = stats_input_with_sections%>% filter(SlideNumber=="s05"),
                                                       model = "percentage ~ Cluster*Treatment*Sex + (1|MouseID)", 
                                                       posthoc1 = "~Treatment|Cluster", 
                                                       posthoc2 = "~Treatment|Cluster|Sex", adjust="sidak")

stats_testing_hypothalamus_s05[[2]]
stats_testing_hypothalamus_s05[[3]]


stats_testing_hypothalamus_s06 <- stats_cluster.animal(data = stats_input_with_sections%>% filter(SlideNumber=="s06"),
                                                       model = "percentage ~ Cluster*Treatment*Sex + (1|MouseID)", 
                                                       posthoc1 = "~Treatment|Cluster", 
                                                       posthoc2 = "~Treatment|Cluster|Sex", adjust="sidak")

stats_testing_hypothalamus_s06[[2]]
stats_testing_hypothalamus_s06[[3]]

stats_testing_hypothalamus_s07 <- stats_cluster.animal(data = stats_input_with_sections%>% filter(SlideNumber=="s07"),
                                                       model = "percentage ~ Cluster*Treatment*Sex + (1|MouseID)", 
                                                       posthoc1 = "~Treatment|Cluster", 
                                                       posthoc2 = "~Treatment|Cluster|Sex", adjust="sidak")

stats_testing_hypothalamus_s07[[2]]
stats_testing_hypothalamus_s07[[3]]


stats_testing_hypothalamus_s08 <- stats_cluster.animal(data = stats_input_with_sections%>% filter(SlideNumber=="s08"),
                                                       model = "percentage ~ Cluster*Treatment*Sex + (1|MouseID)", 
                                                       posthoc1 = "~Treatment|Cluster", 
                                                       posthoc2 = "~Treatment|Cluster|Sex", adjust="sidak")

stats_testing_hypothalamus_s08[[2]]
stats_testing_hypothalamus_s08[[3]]

stats_testing_hypothalamus_s09 <- stats_cluster.animal(data = stats_input_with_sections%>% filter(SlideNumber=="s09"),
                                                       model = "percentage ~ Cluster*Treatment*Sex + (1|MouseID)", 
                                                       posthoc1 = "~Treatment|Cluster", 
                                                       posthoc2 = "~Treatment|Cluster|Sex", adjust="sidak")

stats_testing_hypothalamus_s09[[2]]
stats_testing_hypothalamus_s09[[3]]

stats_testing_hypothalamus_s10 <- stats_cluster.animal(data = stats_input_with_sections%>% filter(SlideNumber=="s10"),
                                                       model = "percentage ~ Cluster*Treatment*Sex + (1|MouseID)", 
                                                       posthoc1 = "~Treatment|Cluster", 
                                                       posthoc2 = "~Treatment|Cluster|Sex", adjust="sidak")

stats_testing_hypothalamus_s10[[2]]
stats_testing_hypothalamus_s10[[3]]

stats_testing_hypothalamus_s11 <- stats_cluster.animal(data = stats_input_with_sections%>% filter(SlideNumber=="s11"),
                                                       model = "percentage ~ Cluster*Treatment*Sex + (1|MouseID)", 
                                                       posthoc1 = "~Treatment|Cluster", 
                                                       posthoc2 = "~Treatment|Cluster|Sex", adjust="sidak")

stats_testing_hypothalamus_s11[[2]]
stats_testing_hypothalamus_s11[[3]]

stats_testing_hypothalamus_s12 <- stats_cluster.animal(data = stats_input_with_sections%>% filter(SlideNumber=="s12"),
                                                       model = "percentage ~ Cluster*Treatment*Sex + (1|MouseID)", 
                                                       posthoc1 = "~Treatment|Cluster", 
                                                       posthoc2 = "~Treatment|Cluster|Sex", adjust="sidak")

stats_testing_hypothalamus_s12[[2]]
stats_testing_hypothalamus_s12[[3]]

#Checking Individual Features 
data <- data_frame_final %>% 
  group_by(MouseID, Sex, Treatment, BrainRegion) %>% 
  summarise(across("Foreground.pixels":"Maximum.branch.length", ~mean(.x))) %>% 
  gather(Measure, Value, "Foreground.pixels":"Maximum.branch.length")

# filter out data you want to run stats on and make sure to make any variables included in model into factors
stats_input <- data %>% filter(BrainRegion=="PVN")
stats_input$Treatment <- factor(stats_input$Treatment)
stats_testing <- stats_morphologymeasures.animal(data = stats_input, 
                                                 model = "Value ~ Treatment*Sex", type="lm",
                                                 posthoc1 = "~Treatment", 
                                                 posthoc2 = "~Treatment|Sex", adjust = "sidak")

stats_testing[[1]]
stats_testing[[2]]
stats_testing[[3]]


stats_input <- data %>% filter(BrainRegion=="HYPO")
stats_input$Treatment <- factor(stats_input$Treatment)
stats_testing <- stats_morphologymeasures.animal(data = stats_input, 
                                                 model = "Value ~ Treatment*Sex", type="lm",
                                                 posthoc1 = "~Treatment", 
                                                 posthoc2 = "~Treatment|Sex", adjust = "sidak")

stats_testing[[1]]
stats_testing[[2]]
stats_testing[[3]]


stats_input <- data %>% filter(BrainRegion=="ARC")
stats_input$Treatment <- factor(stats_input$Treatment)
stats_testing <- stats_morphologymeasures.animal(data = stats_input, 
                                                 model = "Value ~ Treatment*Sex", type="lm",
                                                 posthoc1 = "~Treatment", 
                                                 posthoc2 = "~Treatment|Sex", adjust = "sidak")

stats_testing[[1]]
stats_testing[[2]]
stats_testing[[3]]



