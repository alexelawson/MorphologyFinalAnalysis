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


data_frame_final <- read.csv("/Users/alexlawson/GitHub/MorphologyFinalAnalysis/Data Files/final-dataframe.csv")

data_fram_PVN <- data_frame_final%>%filter(BrainRegion=="PVN")


data_PVN_forstats <- data_fram_PVN  %>% 
  group_by(MouseID, Sex, Treatment) %>% 
  summarise(across("Foreground.pixels":"Maximum.branch.length", ~mean(.x))) %>% 
  gather(Measure, Value, "Foreground.pixels":"Maximum.branch.length")


data_PVN_forstats$Treatment <- factor(data_PVN_forstats$Treatment)
data_PVN_forstats$Sex <- factor(data_PVN_forstats$Sex)
data_PVN_forstats_output <- stats_morphologymeasures.animal(data =data_PVN_forstats, 
                                                 model = "Value ~ Treatment*Sex", type="lm",
                                                 posthoc1 = "~Treatment|Sex", 
                                                 posthoc2 = "~Treatment", adjust = "sidak")
data_PVN_forstats_output[[1]]
data_PVN_forstats_output[[2]]
data_PVN_forstats_output[[3]]
