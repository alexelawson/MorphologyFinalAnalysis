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


data_frame_final <- read.csv("/Users/alexlawson/Desktop/Masters-Work/data-frames/processed-dataframe.csv")
View(data_frame_final)
# prepare data for downstream analysis
data_frame_final_pvn <- data_frame_final %>% filter(BrainRegion=="PVN")
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
stats_testing[[5]]

stats_testing <- stats_morphologymeasures.animal(data = stats_input %>% filter(Sex=="F"), 
                                                 model = "Value ~ Treatment", type="lm",
                                                 posthoc1 = "~Treatment", 
                                                 posthoc2 = "~Treatment", adjust = "sidak")
stats_testing[[1]]
stats_testing[[2]]
stats_testing[[3]]
