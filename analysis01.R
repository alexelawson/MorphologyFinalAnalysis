#THIS IS OLD - DO NOT USE FOR DATA. ONLY FOR REFERENCE (think of it as a lab notebook)
#Code notes: PVN-MR169-s09: redone and merged new data as original slide was off
#            HYPO-MR273-s07: redone and merged new data as there was a massive cell that should have been split 
#            HYPO-MR199-s08: redone and merged new data as there was a massive cell that should have been split 
# Labels for PVN and ARC were swapped. Fixed this in the original dataframe, as well as for the statistical testing. Interim using workingdata01 remains the same.
# Labels on original slides have not been fixed either, so remeber this when using cbc etc. 
# Removed ID's from MR139 s03 because looking back at binary it appeared that the non-microglia had not been removed 

# 
# # raw_data01 <- read.csv("/Users/alexlawson/Desktop/GliaData/workingdata-final-04.csv", stringsAsFactors = FALSE)
# # View(raw_data01)
# # 
# # 
# # 
# # #redone slide mr169_s09  
# # newcellsfinalskeletondir <- "/Users/alexlawson/Desktop/redo-2/skeletondata"
# # newcellsfinalfraclacdir <- "/Users/alexlawson/Desktop/redo-2/fraclac/20241219092849"
# # fraclac_newcellsfinal <- fraclac_tidying(newcellsfinalfraclacdir) 
# # skeleton_newcellsfinal<- skeleton_tidying(newcellsfinalskeletondir)
# # newcellsfinal_merged <- merge_data(fraclac_newcellsfinal, skeleton_newcellsfinal)
# # final_newcellsfinal_merged <- metadata_columns(newcellsfinal_merged, c("BrainRegion","MouseID","Treatment","Sex","SlideNumber"), sep="_")
# # View(final_newcellsfinal_merged)
# # write.csv(final_newcellsfinal_merged, "/Users/alexlawson/Desktop/GliaData/mr278s02mr139s01redone", row.names = FALSE)
# # 
# # final_newcellsfinal_merged<-read.csv("/Users/alexlawson/Desktop/GliaData/mr278s02mr139s01redone")
# # final_newcellsfinal_merged$Sex <- "F"
# # 
# # 
# # # Filter the raw_data01 dataframe
# # filtered_data <- raw_data01[!(
# #   (raw_data01$BrainRegion == "HYPO" & raw_data01$SlideNumber == "s01" & raw_data01$MouseID == "MR139") |
# #     (raw_data01$BrainRegion == "HYPO" & raw_data01$SlideNumber == "s02" & raw_data01$MouseID == "MR278")
# # ), ]
# # 
# # View(filtered_data)
# # 
# # # Merge the filtered dataframe with final_newcellsfinal_merged
# # final_data_testing88 <- rbind(filtered_data, final_newcellsfinal_merged)
# # 
# # # Print the final merged dataframe
# # View(final_data_testing88)
# write.csv(final_data_testing88, "/Users/alexlawson/Desktop/GliaData/workingdata-final-05.csv", row.names = FALSE)
# # 
# # 
# # 
#  raw_data01<-final_data_testing88

#LOADING LIBRARIES AND SETTING WD
library(MicrogliaMorphologyR)
library(factoextra)
library(ppclust)
library(dplyr)
library(ggsignif)
set.seed(1234)
setwd("~/Documents/GitHub/MorphologyAnalysis01")

# ---- Loading and Cleaning the data ----

#LOADING AND CLEANING DATA
#loading in all the skeleton data
skeleton_data_all_dir <- "/Users/alexlawson/Desktop/GliaData/skeleton-data-merg"
skeleton_data_all <- skeleton_tidying(skeleton_data_all_dir)

#loading in all the fraclacdata
hypo_fraclac_dir <- "/Users/alexlawson/Library/CloudStorage/OneDrive-UBC/AlexDataAnalysis/EditedBinary/GliaData/hypo-fraclac/20241210051054"
pvn_fraclac_dir <- "/Users/alexlawson/Library/CloudStorage/OneDrive-UBC/AlexDataAnalysis/EditedBinary/GliaData/pvn-fraclac/20241210074205"
avp_fraclac_dir <- "/Users/alexlawson/Library/CloudStorage/OneDrive-UBC/AlexDataAnalysis/EditedBinary/GliaData/avp-frac-lac/20241210081935"
fraclac_hypo <- fraclac_tidying(hypo_fraclac_dir) 
fraclac_pvn <- fraclac_tidying(pvn_fraclac_dir)
fraclac_avp <- fraclac_tidying(avp_fraclac_dir)
fraclac_data_all <- bind_rows(fraclac_hypo, fraclac_pvn, fraclac_avp)

#merging skeleton and fraclac
data_merged <- merge_data(fraclac_data_all, skeleton_data_all) 
data_merged_cleaned <- data_merged[data_merged$ID!= "0065-3849", ]
View(data_merged_cleaned)

#cleaning column names 
finaldata_merged <- metadata_columns(data_merged_cleaned, c("BrainRegion","MouseID","Treatment","Sex","SlideNumber"), sep="_")
write.csv(finaldata_merged, "/Users/alexlawson/Desktop/GliaData/frac-lac-merged/data_merged_cleaned.csv", row.names = FALSE)
View(finaldata_merged)

#redone slide mr169_s09  
pvn_mr169_s09_skeleton_dir <- "/Users/alexlawson/Desktop/GliaData/pnv-s09169/skeleton-data"
pvn_mr169_s09_fraclac_dir <- "/Users/alexlawson/Desktop/GliaData/pnv-s09169/frac-lac/20241211085426"
fraclac_pvn_mr169_s09 <- fraclac_tidying(pvn_mr169_s09_fraclac_dir) 
skeleton_pvn_mr169_s09 <- skeleton_tidying(pvn_mr169_s09_skeleton_dir)
pvn_mr169_s09_merged <- merge_data(fraclac_pvn_mr169_s09, skeleton_pvn_mr169_s09)
final_pvn_mr169_s09_merged <- metadata_columns(pvn_mr169_s09_merged, c("BrainRegion","MouseID","Treatment","Sex","SlideNumber"), sep="_")
View(final_pvn_mr169_s09_merged)

#redone slide mr273_s07  
hypo_mr273_s07_skeleton_dir <- "/Users/alexlawson/Desktop/GliaData/hypo-s07-273/skeleton-data"
hypo_mr273_s07_fraclac_dir <- "/Users/alexlawson/Desktop/GliaData/hypo-s07-273/fraclac/20241211094727"
fraclac_hypo_mr273_s07 <- fraclac_tidying(hypo_mr273_s07_fraclac_dir) 
skeleton_hypo_mr273_s07 <- skeleton_tidying(hypo_mr273_s07_skeleton_dir)
hypo_mr273_s07_merged <- merge_data(fraclac_hypo_mr273_s07, skeleton_hypo_mr273_s07)
final_hypo_mr273_s07_merged <- metadata_columns(hypo_mr273_s07_merged, c("BrainRegion","MouseID","Treatment","Sex","SlideNumber"), sep="_")
View(final_hypo_mr273_s07_merged)

#redone slide mr199_s08
hypo_mr199_s08_skeleton_dir <- "/Users/alexlawson/Desktop/GliaData/hypo-s08-mr199/skeleton-data"
hypo_mr199_s08_fraclac_dir <- "/Users/alexlawson/Desktop/GliaData/hypo-s08-mr199/fraclac/20241211101301"
fraclac_hypo_mr199_s08 <- fraclac_tidying(hypo_mr199_s08_fraclac_dir) 
skeleton_hypo_mr199_s08 <- skeleton_tidying(hypo_mr199_s08_skeleton_dir)
hypo_mr199_s08_merged <- merge_data(fraclac_hypo_mr199_s08, skeleton_hypo_mr199_s08)
final_hypo_mr199_s08_merged <- metadata_columns(hypo_mr199_s08_merged, c("BrainRegion","MouseID","Treatment","Sex","SlideNumber"), sep="_")
View(final_hypo_mr199_s08_merged)

#removing mr199_s08, mr273_s07, mr169_s09 from the data frame 
#MR169-s09
filtered_finaldata_merged_remove01 <- finaldata_merged %>%
  filter(!(BrainRegion == "PVN" & MouseID == "MR169" & SlideNumber == "s09"))
View(filtered_finaldata_merged_remove01)
#MR273-s07
filtered_finaldata_merged_remove02 <- filtered_finaldata_merged_remove01 %>%
  filter(!(BrainRegion == "HYPO" & MouseID == "MR273" & SlideNumber == "s07"))
View(filtered_finaldata_merged_remove02)
#MR199-s08
filtered_finaldata_merged_remove03 <- filtered_finaldata_merged_remove02 %>%
  filter(!(BrainRegion == "HYPO" & MouseID == "MR199" & SlideNumber == "s08"))
View(filtered_finaldata_merged_remove03)

#merging the cleaned data with the new slides info
cleaned_data01 <- bind_rows(filtered_finaldata_merged_remove03, final_hypo_mr199_s08_merged, final_hypo_mr273_s07_merged, final_pvn_mr169_s09_merged)
View(cleaned_data01)
#saving this file for later use
write.csv(cleaned_data01, "/Users/alexlawson/Desktop/GliaData/workingdata-final-01.csv", row.names = FALSE)
cleaned_data01_labelscorrected_interim <- cleaned_data01 %>%
  mutate(BrainRegion = case_when(
    BrainRegion == "ARC" ~ "TEMP",   # Temporarily rename ARC to TEMP
    BrainRegion == "PVN" ~ "TEMP2",   # Rename PVN to POMC
    TRUE ~ BrainRegion               # Leave other regions unchanged
  ))
View(cleaned_data01_labelscorrected_interim)
cleaned_data01_labelscorrected <- cleaned_data01_labelscorrected_interim %>%
  mutate(BrainRegion = case_when(
    BrainRegion == "TEMP" ~ "PVN",   # Temporarily rename ARC to TEMP
    BrainRegion == "TEMP2" ~ "ARC",   # Rename PVN to POMC
    TRUE ~ BrainRegion               # Leave other regions unchanged
  ))

ids_to_remove <- c(
  "0001-0053", "0002-0084", "0003-0079", "0007-0368", "0008-0407",
  "0010-0496", "0014-0693", "0015-0715", "0016-0744", "0017-0818",
  "0018-0817", "0019-0953", "0022-1127", "0023-1156", "0024-1245",
  "0028-1410", "0032-1605", "0041-1988", "0043-2087", "0046-2264",
  "0047-2297", "0049-2413", "0054-2481", "0055-2523", "0056-2545",
  "0066-2838", "0069-2922", "0072-2986", "0074-3005", "0075-3080",
  "0076-3158", "0079-3206", "0083-3322", "0092-3643", "0094-3672",
  "0095-3700", "0100-3753", "0108-3967", "0116-4180", "0118-4222",
  "0123-4353", "0124-4357", "0126-4389", "0135-4634", "0136-4663",
  "0142-4891", "0151-5106", "0153-5147", "0154-5182", "0157-5212",
  "0158-5224", "0162-5331", "0163-5352", "0169-5988", "0171-6035",
  "0173-6075", "0175-6180", "0179-6302", "0180-6342", "0182-6432",
  "0183-6720", "0184-6739"
)
View(cleaned_data01_labelscorrected)
cleaned_data01_labelscorrected_filtered <- cleaned_data01_labelscorrected[!(cleaned_data01_labelscorrected$MouseID == "MR139" & cleaned_data01_labelscorrected$SlideNumber == "s03" & cleaned_data01_labelscorrected$ID %in% ids_to_remove), ]
View(cleaned_data01_labelscorrected_filtered)
write.csv(cleaned_data01_labelscorrected_filtered, "/Users/alexlawson/Desktop/GliaData/workingdata-final-03-labelscorrected-mr139filtered.csv", row.names = FALSE)

# ---- Investigating features of data ----
#storing in a new variable to work with 
working_data01 <- cleaned_data01_labelscorrected_filtered
View(working_data01)

#INVESTIGATING FEATURES / NORMALIZATION OPTIONS

#feature correlation heamap
featurecorrelations(working_data01,
                    featurestart=8, featureend=34,
                    rthresh=0.8, pthresh=0.05,
                    title="Correlations across features")

#gathering data and investigating outliers
working_data01_gathered <- working_data01 %>% gather(measure, value, 8:ncol(working_data01))
outliers_boxplots(working_data01_gathered)
outliers_distributions(working_data01_gathered)
normalize_logplots(working_data01_gathered,1)
normalize_minmax(working_data01_gathered)
normalize_scaled(working_data01_gathered)

# ---- Scale Data, labels backwards ----
#ScaleOptimization - normalization option 1 
working_data01_scale <- transform_scale(working_data01, start=8, end=34) 
pcadata_elbow(working_data01_scale, featurestart=8, featureend=34)
pca_data01_scale <- pcadata(working_data01_scale, featurestart=8, featureend=34,
                    pc.start=1, pc.end=10)
pca_data01_scale_transformed <- transform_scale(pca_data01_scale, start=1, end=3) # scale pca data as input for k-means clustering
kmeans_input01_scale <- pca_data01_scale_transformed[1:3]
# check for optimal number of clusters using wss and silhouette methods
sampling01_scale <- kmeans_input01_scale[sample(nrow(kmeans_input01_scale), 5000),] #sample 5000 random rows for cluster optimization
fviz_nbclust(sampling01_scale, kmeans, method = 'wss', nstart=25, iter.max=15) # 4 clusters
fviz_nbclust(sampling01_scale, kmeans, method = 'silhouette', nstart=25, iter.max=15)

#4clusters
data_kmeans01_scale <- kmeans(kmeans_input01_scale, centers=4)
pca_kmeans01_scale <- cbind(pca_data01_scale[1:2], working_data01, as.data.frame(data_kmeans01_scale$cluster)) %>%
  rename(Cluster=`data_kmeans01_scale$cluster`) 
clusterfeatures(pca_kmeans01_scale, featurestart=10, featureend=36)

#3clusters 
data_kmeans02_scale <- kmeans(kmeans_input01_scale, centers=3)
pca_kmeans02_scale <- cbind(pca_data01_scale[1:2], working_data01, as.data.frame(data_kmeans02_scale$cluster)) %>%
  rename(Cluster=`data_kmeans02_scale$cluster`) 
clusterfeatures(pca_kmeans02_scale, featurestart=10, featureend=36)

#saving a set of the data to look at on the slide
colorbycluster_scale_3clusters <- pca_kmeans02_scale %>% 
  filter(BrainRegion=="HYPO",MouseID=="MR169", SlideNumber=="s03") %>% select(c(Cluster, ID))
write.csv(colorbycluster_scale_3clusters, "/Users/alexlawson/Desktop/for-jessica/hypo-mr169s03-scale-3clusters.csv", row.names = FALSE)

#saving a set of the data to look at on the slide
colorbycluster_scale_3clusters_02 <- pca_kmeans02_scale %>% 
  filter(BrainRegion=="HYPO",MouseID=="MR278", SlideNumber=="s06") %>% select(c(Cluster, ID))
write.csv(colorbycluster_scale_3clusters_02, "/Users/alexlawson/Desktop/for-jessica/hypo-mr278s06-scale-3clusters.csv", row.names = FALSE)

#saving a set of the data to look at on the slide
colorbycluster_scale_3clusters_03 <- pca_kmeans02_scale %>% 
  filter(BrainRegion=="HYPO",MouseID=="MR290", SlideNumber=="s03") %>% select(c(Cluster, ID))
write.csv(colorbycluster_scale_3clusters_03, "/Users/alexlawson/Desktop/for-jessica/hypo-mr290s03-scale-3clusters.csv", row.names = FALSE)

#saving a set of the data to look at on the slide
colorbycluster_scale_3clusters_04 <- pca_kmeans02_scale %>% 
  filter(BrainRegion=="PVN",MouseID=="MR287", SlideNumber=="s07") %>% select(c(Cluster, ID))
write.csv(colorbycluster_scale_3clusters_04, "/Users/alexlawson/Desktop/for-jessica/pvn-mr287s07-scale-3clusters.csv", row.names = FALSE)

#saving a set of the data to look at on the slide
colorbycluster_scale <- pca_kmeans01_scale %>% 
  filter(BrainRegion=="HYPO",MouseID=="MR169", SlideNumber=="s03") %>% select(c(Cluster, ID))
write.csv(colorbycluster_scale, "/Users/alexlawson/Desktop/for-jessica/hypo-mr169s03-scale.csv", row.names = FALSE)

#saving a set of the data to look at on the slide
colorbycluster_scale2 <- pca_kmeans01_scale %>% 
  filter(BrainRegion=="HYPO",MouseID=="MR290", SlideNumber=="s03") %>% select(c(Cluster, ID))
write.csv(colorbycluster_scale2, "/Users/alexlawson/Desktop/for-jessica/hypo-mr290s03-scale.csv", row.names = FALSE)

#saving a set of the data to look at on the slide
colorbycluster_scale3 <- pca_kmeans01_scale %>% 
  filter(BrainRegion=="HYPO",MouseID=="MR278", SlideNumber=="s06") %>% select(c(Cluster, ID))
write.csv(colorbycluster_scale3, "/Users/alexlawson/Desktop/for-jessica/hypo-mr278s06-scale.csv", row.names = FALSE)

#saving a set of the data to look at on the slide
colorbycluster_scale4 <- pca_kmeans01_scale %>% 
  filter(BrainRegion=="HYPO",MouseID=="MR157", SlideNumber=="s07") %>% select(c(Cluster, ID))
write.csv(colorbycluster_scale4, "/Users/alexlawson/Desktop/for-jessica/hypo-mr157s07-scale.csv", row.names = FALSE)

#saving a set of the data to look at on the slide
colorbycluster_scale5 <- pca_kmeans01_scale %>% 
  filter(BrainRegion=="HYPO",MouseID=="MR139", SlideNumber=="s05") %>% select(c(Cluster, ID))
write.csv(colorbycluster_scale5, "/Users/alexlawson/Desktop/for-jessica/hypo-mr139s05-scale.csv", row.names = FALSE)

# ---- Log optimization ----

#LogOptimization - normalization option 2
working_data01_log <- transform_log(working_data01,1, start=8, end=34) 
pcadata_elbow(working_data01_log, featurestart=8, featureend=34)
pca_data01_log <- pcadata(working_data01_log, featurestart=8, featureend=34,
                       pc.start=1, pc.end=10)

pcfeaturecorrelations(pca_data01_log, pc.start=1, pc.end=3, 
                      feature.start=18, feature.end=44, 
                      rthresh=0.75, pthresh=0.05, 
                      title="Correlation between PCs and features")
pca_data01_log_transformed <- transform_scale(pca_data01_log, start=1, end=3) # scale pca data as input for k-means clustering
kmeans_input01_log <- pca_data01_log_transformed[1:3]
# check for optimal number of clusters using wss and silhouette methods
sampling01_log <- kmeans_input01_log[sample(nrow(kmeans_input01_log), 5000),] #sample 5000 random rows for cluster optimization
fviz_nbclust(sampling01_log, kmeans, method = 'wss', nstart=25, iter.max=15) # 4 clusters
fviz_nbclust(sampling01_log, kmeans, method = 'silhouette', nstart=25, iter.max=15)
# ---- 4 clusters Log Data, labels backwards ----
#4clusters
data_kmeans01_log <- kmeans(kmeans_input01_log, centers=4)
pca_kmeans01_log <- cbind(pca_data01_log[1:2], working_data01, as.data.frame(data_kmeans01_log$cluster)) %>%
  rename(Cluster=`data_kmeans01_log$cluster`) 
clusterfeatures(pca_kmeans01_log, featurestart=10, featureend=36)

# ---- 3 clusters Log Data ----
#3clusters
data_kmeans02_log <- kmeans(kmeans_input01_log, centers=3)
pca_kmeans02_log <- cbind(pca_data01_log[1:3], working_data01, as.data.frame(data_kmeans02_log$cluster)) %>%
  rename(Cluster=`data_kmeans02_log$cluster`) 
View(pca_kmeans02_log)
clusterfeatures(pca_kmeans02_log, featurestart=11, featureend=37)
pca_kmeans02_log_labelscorrected <- pca_kmeans02_log

# ---- Saving 3 cluster log data ----
#saving a set of the data to look at on the slide
colorbycluster_log_3clust <- pca_kmeans02_log %>% 
  filter(BrainRegion=="HYPO",MouseID=="MR169", SlideNumber=="s03") %>% select(c(Cluster, ID))
write.csv(colorbycluster_log_3clust, "/Users/alexlawson/Desktop/for-jessica/hypo-mr169s03-log-3clusters.csv", row.names = FALSE)

#saving a set of the data to look at on the slide
colorbycluster_log_3clusters_02 <- pca_kmeans02_log %>% 
  filter(BrainRegion=="HYPO",MouseID=="MR278", SlideNumber=="s06") %>% select(c(Cluster, ID))
write.csv(colorbycluster_log_3clusters_02, "/Users/alexlawson/Desktop/for-jessica/hypo-mr278s06-log-3clusters.csv", row.names = FALSE)

#saving a set of the data to look at on the slide
colorbycluster_log_3clusters_03 <- pca_kmeans02_log %>% 
  filter(BrainRegion=="HYPO",MouseID=="MR290", SlideNumber=="s03") %>% select(c(Cluster, ID))
write.csv(colorbycluster_log_3clusters_03, "/Users/alexlawson/Desktop/for-jessica/hypo-mr290s03-log-3clusters.csv", row.names = FALSE)

#saving a set of the data to look at on the slide
colorbycluster_log_3clusters_04 <- pca_kmeans02_log %>% 
  filter(BrainRegion=="HYPO",MouseID=="MR277", SlideNumber=="s06") %>% select(c(Cluster, ID))
write.csv(colorbycluster_log_3clusters_04, "/Users/alexlawson/Desktop/for-jessica/hypo-mr277s06-log-3clusters.csv", row.names = FALSE)

#saving a set of the data to look at on the slide
colorbycluster_log_3clusters_05 <- pca_kmeans02_log %>% 
  filter(BrainRegion=="HYPO",MouseID=="MR198", SlideNumber=="s09") %>% select(c(Cluster, ID))
write.csv(colorbycluster_log_3clusters_05, "/Users/alexlawson/Desktop/for-jessica/hypo-mr198s09-log-3clusters.csv", row.names = FALSE)

#saving a set of the data to look at on the slide
colorbycluster_log_3clusters_06 <- pca_kmeans02_log %>% 
  filter(BrainRegion=="PVN",MouseID=="MR287", SlideNumber=="s07") %>% select(c(Cluster, ID))
write.csv(colorbycluster_log_3clusters_06, "/Users/alexlawson/Desktop/for-jessica/pvn-mr287s07-log-3clusters.csv", row.names = FALSE)

#saving a set of the data to look at on the slide
colorbycluster_log_3clusters_07 <- pca_kmeans02_log %>% 
  filter(BrainRegion=="ARC",MouseID=="MR199", SlideNumber=="s05") %>% select(c(Cluster, ID))
write.csv(colorbycluster_log_3clusters_07, "/Users/alexlawson/Desktop/for-jessica/arc-mr199s05-log-3clusters.csv", row.names = FALSE)

colorbycluster_log_3clusters_arc_1 <- pca_kmeans02_log %>% 
  filter(BrainRegion=="ARC",MouseID=="MR140", SlideNumber=="s04") %>% select(c(Cluster, ID))
write.csv(colorbycluster_log_3clusters_arc_1, "/Users/alexlawson/Desktop/GliaData/output-images/arc-cbc/arc-mr140s04-log-3clusters.csv", row.names = FALSE)

colorbycluster_log_3clusters_arc_2 <- pca_kmeans02_log %>% 
  filter(BrainRegion=="ARC",MouseID=="MR199", SlideNumber=="s06") %>% select(c(Cluster, ID))
write.csv(colorbycluster_log_3clusters_arc_2, "/Users/alexlawson/Desktop/GliaData/output-images/arc-cbc/arc-mr199s06-log-3clusters.csv", row.names = FALSE)

colorbycluster_log_3clusters_arc_2 <- pca_kmeans02_log %>% 
  filter(BrainRegion=="ARC",MouseID=="MR287", SlideNumber=="s06") %>% select(c(Cluster, ID))
write.csv(colorbycluster_log_3clusters_arc_2, "/Users/alexlawson/Desktop/GliaData/output-images/arc-cbc/arc-mr169s05-log-3clusters.csv", row.names = FALSE)


#saving a set of the data to look at on the slide
colorbycluster_log <- pca_kmeans01_log %>% 
  filter(BrainRegion=="HYPO",MouseID=="MR169", SlideNumber=="s03") %>% select(c(Cluster, ID))
write.csv(colorbycluster_log, "/Users/alexlawson/Desktop/for-jessica/hypo-mr169s03-log.csv", row.names = FALSE)

# ---- Minmax, labels backwards ----
#Minmax-Optimization - normalization option 3 
working_data01_minmax <- transform_minmax(working_data01, start=8, end=34) 
pcadata_elbow(working_data01_minmax, featurestart=8, featureend=34)
pca_data01_minmax <- pcadata(working_data01_minmax, featurestart=8, featureend=34,
                       pc.start=1, pc.end=10)
pca_data01_minmax_transformed <- transform_scale(pca_data01_minmax, start=1, end=3) # scale pca data as input for k-means clustering
kmeans_input01_minmax <- pca_data01_minmax_transformed[1:3]
# check for optimal number of clusters using wss and silhouette methods
sampling01_minmax <- kmeans_input01_minmax[sample(nrow(kmeans_input01_minmax), 5000),] #sample 5000 random rows for cluster optimization
fviz_nbclust(sampling01_minmax, kmeans, method = 'wss', nstart=25, iter.max=15) # 4 clusters
fviz_nbclust(sampling01_minmax, kmeans, method = 'silhouette', nstart=25, iter.max=15)
data_kmeans01_minmax <- kmeans(kmeans_input01_minmax, centers=4)
pca_kmeans01_minmax <- cbind(pca_data01_minmax[1:2], working_data01, as.data.frame(data_kmeans01_minmax$cluster)) %>%
  rename(Cluster=`data_kmeans01_minmax$cluster`) 
clusterfeatures(pca_kmeans01_minmax, featurestart=10, featureend=36)
#saving a set of the data to look at on the slide
colorbycluster_minmax <- pca_kmeans01_minmax %>% 
  filter(BrainRegion=="HYPO",MouseID=="MR169", SlideNumber=="s03") %>% select(c(Cluster, ID))
  write.csv(colorbycluster_minmax, "/Users/alexlawson/Desktop/for-jessica/hypo-mr169s03-minmax.csv", row.names = FALSE)
# ---- Formatting data for stats ----
#Working through the log data
plot_log_3 <- clusterplots(pca_kmeans02_log_labelscorrected, "PC1", "PC2")
plot_log_3
cp_log_3 <- clusterpercentage(pca_kmeans02_log_labelscorrected, "Cluster", Treatment, Sex, BrainRegion, MouseID)

cp_log_testing <- clusterpercentage(pca_kmeans02_log_labelscorrected, "Cluster", Treatment, Sex)
cp_log_3 <- cp_log_3 %>% mutate(Cluster = 
                      case_when(Cluster=="1" ~ "Ameboid",
                                Cluster=="2" ~ "Rod-Like",
                                Cluster=="3" ~ "Ramified"))
# ---- Plots and such 3 clusters log data, labels correct ----

#graph split by brain region and treatment, mouse id labeled comparing cold vs. control 
cp_log_3 %>% 
  ggplot(aes(x = Cluster, y = percentage, group = interaction(Cluster, Treatment))) +
  facet_grid(Sex ~ BrainRegion) + # Facets by Sex (rows) and BrainRegion (columns)
  geom_boxplot(aes(group = interaction(Cluster, Treatment), fill = Treatment)) +
  scale_fill_manual(values = c("lightblue", "orange")) +
  ggtitle("Mouse dataset: K-means clusters") +
  labs(fill = "Treatment") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

cp_log_3 %>% 
  ggplot(aes(x = Cluster, y = percentage, group = interaction(Cluster, Treatment))) +
  facet_grid(Sex ~ BrainRegion) + 
  geom_boxplot(aes(group = interaction(Cluster, Treatment), fill = Treatment), outlier.shape = NA) +  # Remove outliers
  scale_fill_manual(values = c("lightblue", "white")) +
  ggtitle("Cold VS. Control Cluster Percentages Split by Brain Region and Sex") +
  labs(x = "Cluster", y = "Percentage", fill = "Treatment") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))




#control male vs. female in the different brain regions
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

#cold male vs. female in the different brain regions 
cp_log_3 %>%
  filter(Treatment == "COLD", Sex %in% c("M", "F")) %>%  # Filter for Control Males and Females
  ggplot(aes(x = Cluster, y = percentage, group = interaction(Cluster, Sex))) +
  facet_grid(. ~ BrainRegion) +  # Columns = Brain Region
  geom_boxplot(aes(fill = Sex), position = position_dodge(width = 0.8)) +  # Compare sexes
  scale_fill_manual(values = c("lightpink", "lightblue")) +
  geom_point(position = position_dodge(width = 0.8), size = 0.75) +  # Add points for clarity
  ggtitle("Mouse Dataset: Cold Males vs Females Across Brain Regions") +
  labs(fill = "Sex") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

#graph splitting up slides. looking at cold vs. control in males and females  
cp_log_3_withslides %>%
  filter(BrainRegion == "HYPO") %>%  # Filter to include only HYPO data
  ggplot(aes(x = Cluster, y = percentage, group = interaction(Cluster, Treatment))) +
  facet_grid(Sex ~ SlideNumber) +  # Facet by SlideNumber instead of BrainRegion
  geom_boxplot(aes(group = interaction(Cluster, Treatment), fill = Treatment)) +
  scale_fill_manual(values = c("lightblue", "orange")) +
  geom_point(position = position_dodge(width = 0.8), size = 0.75) +
  ggtitle("Mouse Dataset: Clusters in the Hypothalamus Separated by Section") +
  labs(fill = "Treatment") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

# First plot for SlideNumber s01-s06
plot_s01_s06 <- cp_log_3_withslides %>%
  filter(BrainRegion == "HYPO", SlideNumber %in% c("s01", "s02", "s03", "s04", "s05", "s06")) %>%
  ggplot(aes(x = Cluster, y = percentage, group = interaction(Cluster, Treatment))) +
  facet_grid(Sex ~ SlideNumber) +
  geom_boxplot(aes(group = interaction(Cluster, Treatment), fill = Treatment), outlier.shape = NA) +
  scale_fill_manual(values = c("yellow", "lavender")) +
  ggtitle("Mouse Dataset: Cold vs. Control in the Hypothalamus, Separated by Sex and Section (s01-s06)") +
  labs(fill = "Treatment") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

plot_s01_s06 <- cp_log_3_withslides %>%
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


plot_s01_s06







# Second plot for SlideNumber s07-s12
plot_s07_s12 <- cp_log_3_withslides %>%
  filter(BrainRegion == "HYPO", SlideNumber %in% c("s07", "s08", "s09", "s10", "s11", "s12")) %>%
  ggplot(aes(x = Cluster, y = percentage, group = interaction(Cluster, Treatment))) +
  facet_grid(Sex ~ SlideNumber) +
  geom_boxplot(aes(group = interaction(Cluster, Treatment), fill = Treatment), outlier.shape = NA) +
  scale_fill_manual(values = c("yellow", "lavender")) +
  ggtitle("Mouse Dataset: Cold vs. Control in the Hypothalamus, Separated by Sex and Section (s07-s12)") +
  labs(fill = "Treatment") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

# Print the plots
plot_s01_s06
plot_s07_s12



# Create the combined plot
combined_plot <- cp_log_3_withslides %>%
  filter(BrainRegion == "HYPO", Sex=="F", SlideNumber %in% c("s01", "s02", "s03", "s04", "s05", "s06", 
                                                   "s07", "s08", "s09", "s10", "s11", "s12")) %>%
  ggplot(aes(x = Cluster, y = percentage, group = interaction(Cluster, Treatment))) +
  facet_wrap(~ SlideNumber, nrow = 2) +  # Facet by SlideGroup with two rows
  geom_boxplot(aes(group = interaction(Cluster, Treatment), fill = Treatment), outlier.shape = NA) +
  scale_fill_manual(values = c("yellow", "lavender")) +
  ggtitle("Mouse Dataset: Cold vs. Control in the Hypothalamus, Separated by Section") +
  labs(fill = "Treatment") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

# Print the combined plot
combined_plot



# First plot for SlideNumber s01-s06
plot_s01_s06_con <- cp_log_3_withslides %>%
  filter(BrainRegion == "HYPO", Treatment == "CON", Sex %in% c("M", "F"), SlideNumber %in% c("s01", "s02", "s03", "s04", "s05", "s06")) %>%
  ggplot(aes(x = Cluster, y = percentage, group = interaction(Cluster, Sex))) +
  facet_grid(. ~ SlideNumber) +  # Facet by SlideNumber
  geom_boxplot(aes(fill = Sex), outlier.shape = NA) +
  scale_fill_manual(values = c("lightpink", "lightblue")) +  # Different colors for Sex
  ggtitle("Mouse Dataset: Males vs Females in Controls in Hypothalamus (s01-s06)") +
  labs(fill = "Sex") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

# Second plot for SlideNumber s07-s12
plot_s07_s12_con <- cp_log_3_withslides %>%
  filter(BrainRegion == "HYPO", Treatment == "CON", Sex %in% c("M", "F"), SlideNumber %in% c("s07", "s08", "s09", "s10", "s11", "s12")) %>%
  ggplot(aes(x = Cluster, y = percentage, group = interaction(Cluster, Sex))) +
  facet_grid(. ~ SlideNumber) +  # Facet by SlideNumber
  geom_boxplot(aes(fill = Sex), outlier.shape = NA) +
  scale_fill_manual(values = c("lightpink", "lightblue")) +  # Different colors for Sex
  ggtitle("Mouse Dataset: Males vs Females in Controls in Hypothalamus (s07-s12)") +
  labs(fill = "Sex") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16))

# Print the plots
plot_s01_s06_con
plot_s07_s12_con


# Create the combined plot
combined_plot_controls <- cp_log_3_withslides %>%
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

# Print the combined plot
combined_plot_controls


#graph splitting up slides. looking at female vs. male in cold and control 
cp_log_3_withslides %>%
  filter(BrainRegion == "HYPO") %>%  # Filter to include only HYPO data
  ggplot(aes(x = Cluster, y = percentage, group = interaction(Cluster, Sex))) +
  facet_grid(Treatment ~ SlideNumber) +  # Rows = Treatment, Columns = SlideNumber
  geom_boxplot(aes(fill = Sex), position = position_dodge(width = 0.8)) +  # Compare sexes
  scale_fill_manual(values = c("lightpink", "lightblue")) +
  geom_point(position = position_dodge(width = 0.8), size = 0.75) +  # Add points for clarity
  ggtitle("Mouse Dataset: Difference in clusters between sexes in cold and control conditions in the hypothalamus, split by section number") +
  labs(fill = "Sex") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
        strip.text = element_text(face = "bold", size = 12),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16))


# ---- Stats 3 clusters log data, labels correct ----
# prepare percentages dataset for downstream analysis
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
                                            model = "percentage ~ Cluster*Treatment*BrainRegion*Sex", 
                                            posthoc1 = "~Treatment|Cluster|BrainRegion", 
                                            posthoc2 = "~Treatment|Cluster|Sex|BrainRegion", adjust = "fdr")

stats_testing_all[[1]]
stats_testing_all[[2]]
stats_testing_all[[3]]
stats_testing_all[[4]]

write.csv(stats_testing_all[[3]], "/Users/alexlawson/Desktop/GliaData/final-data/stats-results/treatment-cluster-brainregion.csv", row.names = FALSE)


stats_testing_males <- stats_cluster.animal(data = stats_input %>% filter(Sex=="M"), 
                                           model = "percentage ~ Cluster*Treatment*BrainRegion + (1|MouseID)", 
                                           posthoc1 = "~Treatment|Cluster", 
                                           posthoc2 = "~Treatment|Cluster|BrainRegion", adjust = "fdr")


stats_testing_males[[1]]
stats_testing_males[[2]]
stats_testing_males[[3]]
stats_testing_males[[4]]
stats_testing_males[[5]]

saveRDS(stats_testing_males, file = "/Users/alexlawson/Desktop/GliaData/final-data/stats-results/males.rds")
write.csv(stats_testing_males[[3]], "/Users/alexlawson/Desktop/GliaData/final-data/stats-results/treatment-cluster-brainregion-males.csv", row.names = FALSE)
saveRDS(stats_testing_fem, file = "/Users/alexlawson/Desktop/GliaData/final-data/stats-results/females.rds")
write.csv(stats_testing_fem[[3]], "/Users/alexlawson/Desktop/GliaData/final-data/stats-results/treatment-cluster-brainregion-females.csv", row.names = FALSE)

View(stats_input)
stats_testing_fem <- stats_cluster.animal(data = stats_input %>% filter(Sex=="F"), 
                                            model = "percentage ~ Cluster*Treatment*BrainRegion + (1|MouseID)", 
                                            posthoc1 = "~Treatment|Cluster", 
                                            posthoc2 = "~Treatment|Cluster|BrainRegion", adjust = "fdr")
stats_testing_fem[[1]]
stats_testing_fem[[2]]
stats_testing_fem[[3]]

stats_input_pvn <- stats_input %>% filter(BrainRegion=="PVN")
stats_testing_pvn_male <- stats_cluster.animal(data = stats_input_pvn %>% filter(Sex=="M"), 
                                          model = "percentage ~ Cluster*Treatment", 
                                          posthoc1 = "~Treatment|Cluster", 
                                          posthoc2 = "~Treatment|Cluster", adjust = "fdr")

stats_testing_pvn_male[[2]]


cp_log_3_withslides <- clusterpercentage(pca_kmeans02_log_labelscorrected, "Cluster", Treatment, Sex, BrainRegion, MouseID, SlideNumber)
cp_log_3_withslides <- cp_log_3_withslides %>%
  mutate(SlideNumber = case_when(
    SlideNumber == "s07-1" ~ "s07",   # Temporarily rename ARC to TEMP
    TRUE ~ SlideNumber               # Leave other regions unchanged
  ))

cp_log_3_withslides <- cp_log_3_withslides %>% mutate(Cluster = 
                                  case_when(Cluster=="1" ~ "Ameboid",
                                            Cluster=="2" ~ "Rod-Like",
                                            Cluster=="3" ~ "Ramified"))
stats_input_slides <- cp_log_3_withslides 
stats_input_slides$MouseID <- factor(stats_input_slides$MouseID)
stats_input_slides$Cluster <- factor(stats_input_slides$Cluster)
stats_input_slides$Treatment <- factor(stats_input_slides$Treatment)
stats_input_slides$SlideNumber <- factor(stats_input_slides$SlideNumber)
View(stats_input_slides)


# ---- Stats 3 HYPOTHALAMUS by slide log data, labels correct ----
stats_input_slides_hypothalamus <- stats_input_slides %>% filter(BrainRegion=="HYPO")
View(stats_input_slides_hypothalamus)

stats_testing_hypothalamus <- stats_cluster.animal(data = stats_input_slides_hypothalamus, 
                                                        model = "percentage ~ Cluster*Treatment*SlideNumber*Sex", 
                                                        posthoc1 = "~Treatment|Cluster|SlideNumber|Sex", 
                                                        posthoc2 = "~Treatment|Cluster|Sex|SlideNumber")

stats_testing_hypothalamus[[1]]
stats_testing_hypothalamus[[2]]
stats_testing_hypothalamus[[3]]

write.csv(stats_testing_hypothalamus[[2]], "/Users/alexlawson/Desktop/GliaData/final-data/stats-results/sex-treatment-cluster-slidenumber-nopvalueadjust", row.names = FALSE)


stats_testing_hypothalamus_male <- stats_cluster.animal(data = stats_input_slides_hypothalamus %>% filter(Sex=="M"), 
                                                         model = "percentage ~ Cluster*Treatment*SlideNumber + (1|MouseID)", 
                                                         posthoc1 = "~Treatment|Cluster", 
                                                         posthoc2 = "~Treatment|Cluster|SlideNumber")

stats_testing_hypothalamus_male[[1]]
stats_testing_hypothalamus_male[[2]]
stats_testing_hypothalamus_male[[3]]

write.csv(stats_testing_hypothalamus_male[[3]], "/Users/alexlawson/Desktop/GliaData/final-data/stats-results/treatment-cluster-slidenumber-males-nopvalueadjust.csv", row.names = FALSE)


stats_testing_hypothalamus_con <- stats_cluster.animal(data = stats_input_slides_hypothalamus %>% filter(Treatment=="CON"), 
                                                        model = "percentage ~ Cluster*Sex*SlideNumber + (1|MouseID)", 
                                                        posthoc1 = "~Sex|Cluster", 
                                                        posthoc2 = "~Sex|Cluster|SlideNumber", adjust = "fdr")
stats_testing_hypothalamus_con[[1]]
stats_testing_hypothalamus_con[[2]]
stats_testing_hypothalamus_con[[3]]

stats_testing_hypothalamus_cold <- stats_cluster.animal(data = stats_input_slides_hypothalamus %>% filter(Treatment=="COLD"), 
                                                       model = "percentage ~ Cluster*Sex*SlideNumber + (1|MouseID)", 
                                                       posthoc1 = "~Cluster|Sex", 
                                                       posthoc2 = "~Sex|Cluster|SlideNumber", adjust = "fdr")
stats_testing_hypothalamus_cold[[1]]
stats_testing_hypothalamus_cold[[2]]
stats_testing_hypothalamus_cold[[3]]



stats_testing_hypothalamus_female <- stats_cluster.animal(data = stats_input_slides_hypothalamus %>% filter(Sex=="F"), 
                                                        model = "percentage ~ Cluster*Treatment*SlideNumber + (1|MouseID)", 
                                                        posthoc1 = "~Treatment|Cluster", 
                                                        posthoc2 = "~Treatment|Cluster|SlideNumber")

stats_testing_hypothalamus_female[[1]]
stats_testing_hypothalamus_female[[2]]
stats_testing_hypothalamus_female[[3]]

stats_testing_slides_hypothalamus01 <- stats_cluster.animal(data = stats_input_slides_hypothalamus %>% filter(SlideNumber=="s01"), 
                                          model = "percentage ~ Cluster*Treatment*Sex + (1|MouseID)", 
                                          posthoc1 = "~Sex|Cluster", 
                                          posthoc2 = "~Treatment|Cluster|Sex")

stats_testing_slides_hypothalamus01[[1]]
stats_testing_slides_hypothalamus01[[2]]
stats_testing_slides_hypothalamus01[[3]]
write.csv(stats_testing_hypothalamus_female[[3]], "/Users/alexlawson/Desktop/GliaData/final-data/stats-results/treatment-cluster-slidenumber-females-nopvalueadjust", row.names = FALSE)


stats_testing_slides_hypothalamus02 <- stats_cluster.animal(data = stats_input_slides_hypothalamus %>% filter(SlideNumber=="s02"), 
                                                            model = "percentage ~ Cluster*Treatment*Sex + (1|MouseID)", 
                                                            posthoc1 = "~Treatment|Cluster", 
                                                            posthoc2 = "~Treatment|Cluster|Sex", adjust = "fdr")


stats_testing_slides_hypothalamus02[[1]]
stats_testing_slides_hypothalamus02[[2]]
stats_testing_slides_hypothalamus02[[3]]

stats_testing_slides_hypothalamus03 <- stats_cluster.animal(data = stats_input_slides_hypothalamus %>% filter(SlideNumber=="s03"), 
                                                            model = "percentage ~ Cluster*Treatment*Sex", 
                                                            posthoc1 = "~Treatment|Cluster", 
                                                            posthoc2 = "~Treatment|Cluster|Sex")


stats_testing_slides_hypothalamus03[[1]]
stats_testing_slides_hypothalamus03[[2]]
stats_testing_slides_hypothalamus03[[3]]

stats_testing_slides_hypothalamus04 <- stats_cluster.animal(data = stats_input_slides_hypothalamus %>% filter(SlideNumber=="s04"), 
                                                            model = "percentage ~ Cluster*Treatment*Sex + (1|MouseID)", 
                                                            posthoc1 = "~Treatment|Cluster", 
                                                            posthoc2 = "~Treatment|Cluster|Sex", adjust = "fdr")


stats_testing_slides_hypothalamus04[[1]]
stats_testing_slides_hypothalamus04[[2]]
stats_testing_slides_hypothalamus04[[3]]

stats_testing_slides_hypothalamus05 <- stats_cluster.animal(data = stats_input_slides_hypothalamus %>% filter(SlideNumber=="s05"), 
                                                            model = "percentage ~ Cluster*Treatment*Sex + (1|MouseID)", 
                                                            posthoc1 = "~Treatment|Cluster", 
                                                            posthoc2 = "~Treatment|Cluster|Sex", adjust = "fdr")


stats_testing_slides_hypothalamus05[[1]]
stats_testing_slides_hypothalamus05[[2]]
stats_testing_slides_hypothalamus05[[3]]

stats_testing_slides_hypothalamus06 <- stats_cluster.animal(data = stats_input_slides_hypothalamus %>% filter(SlideNumber=="s06"), 
                                                            model = "percentage ~ Cluster*Treatment*Sex + (1|MouseID)", 
                                                            posthoc1 = "~Treatment|Cluster", 
                                                            posthoc2 = "~Treatment|Cluster|Sex", adjust = "fdr")


stats_testing_slides_hypothalamus06[[1]]
stats_testing_slides_hypothalamus06[[2]]
stats_testing_slides_hypothalamus06[[3]]

stats_testing_slides_hypothalamus07 <- stats_cluster.animal(data = stats_input_slides_hypothalamus %>% filter(SlideNumber=="s07"), 
                                                            model = "percentage ~ Cluster*Treatment*Sex + (1|MouseID)", 
                                                            posthoc1 = "~Treatment|Cluster", 
                                                            posthoc2 = "~Treatment|Cluster|Sex", adjust = "fdr")


stats_testing_slides_hypothalamus07[[1]]
stats_testing_slides_hypothalamus07[[2]]
stats_testing_slides_hypothalamus07[[3]]

stats_testing_slides_hypothalamus08 <- stats_cluster.animal(data = stats_input_slides_hypothalamus %>% filter(SlideNumber=="s08"), 
                                                            model = "percentage ~ Cluster*Treatment*Sex + (1|MouseID)", 
                                                            posthoc1 = "~Treatment|Cluster", 
                                                            posthoc2 = "~Treatment|Cluster|Sex", adjust = "fdr")


stats_testing_slides_hypothalamus08[[1]]
stats_testing_slides_hypothalamus08[[2]]
stats_testing_slides_hypothalamus08[[3]]

stats_testing_slides_hypothalamus09 <- stats_cluster.animal(data = stats_input_slides_hypothalamus %>% filter(SlideNumber=="s09"), 
                                                            model = "percentage ~ Cluster*Treatment*Sex + (1|MouseID)", 
                                                            posthoc1 = "~Treatment|Cluster", 
                                                            posthoc2 = "~Treatment|Cluster|Sex", adjust = "fdr")


stats_testing_slides_hypothalamus09[[1]]
stats_testing_slides_hypothalamus09[[2]]
stats_testing_slides_hypothalamus09[[3]]

stats_testing_slides_hypothalamus10 <- stats_cluster.animal(data = stats_input_slides_hypothalamus %>% filter(SlideNumber=="s10"), 
                                                            model = "percentage ~ Cluster*Treatment*Sex + (1|MouseID)", 
                                                            posthoc1 = "~Treatment|Cluster", 
                                                            posthoc2 = "~Treatment|Cluster|Sex", adjust = "fdr")


stats_testing_slides_hypothalamus10[[1]]
stats_testing_slides_hypothalamus10[[2]]
stats_testing_slides_hypothalamus10[[3]]

stats_testing_slides_hypothalamus11 <- stats_cluster.animal(data = stats_input_slides_hypothalamus %>% filter(SlideNumber=="s11"), 
                                                            model = "percentage ~ Cluster*Treatment*Sex + (1|MouseID)", 
                                                            posthoc1 = "~Treatment|Cluster", 
                                                            posthoc2 = "~Treatment|Cluster|Sex", adjust = "fdr")


stats_testing_slides_hypothalamus11[[1]]
stats_testing_slides_hypothalamus11[[2]]
stats_testing_slides_hypothalamus11[[3]]

stats_testing_slides_hypothalamus12 <- stats_cluster.animal(data = stats_input_slides_hypothalamus %>% filter(SlideNumber=="s12"), 
                                                            model = "percentage ~ Cluster*Treatment*Sex + (1|MouseID)", 
                                                            posthoc1 = "~Treatment|Cluster", 
                                                            posthoc2 = "~Treatment|Cluster|Sex", adjust = "fdr")


stats_testing_slides_hypothalamus12[[1]]
stats_testing_slides_hypothalamus12[[2]]
stats_testing_slides_hypothalamus12[[3]]


# ---- Fuzzyclustering  ----
#fuzzy-clustering
fuzzy_clustering <- fcm(kmeans_input01_log, centers=3, nstart=25)

View(fuzzy_clustering)
summary(fuzzy_clustering)

saveRDS(fuzzy_clustering, file = "/Users/alexlawson/Library/CloudStorage/OneDrive-UBC/AlexDataAnalysis/final-data/fuzzy_clustering.rds")

fuzzyclustering_dataframe <- cbind(workingdat01_wcluster, fuzzy_clustering$u)
View(fuzzyclustering_dataframe)
clusterfeatures(fuzzyclustering_dataframe, featurestart=8, featureend=34)
fuzzyclustering_dataframe_highfit <- fuzzyclustering_dataframe %>% 
  filter(`Cluster 1` > 0.70|
           `Cluster 2` > 0.70|
           `Cluster 3` > 0.70)

cp_highfit <- clusterpercentage(fuzzyclustering_dataframe_highfit, "Cluster", MouseID, Treatment, Sex, BrainRegion)
cp_highfit$Treatment <- factor(cp_highfit$Treatment, levels=c("CON","COLD"))

cp_highfit %>% 
  filter(BrainRegion=="HYPO") %>%
  ggplot(aes(x=Cluster, y=percentage, group=interaction(Cluster, Treatment))) +
  facet_wrap(~Sex) +
  geom_boxplot(aes(group=interaction(Cluster, Treatment), fill=Treatment)) +
  scale_fill_manual(values=c("#fde725","#482878")) +
  geom_point(position=position_dodge(width=0.8), size=0.75, aes(group=interaction(Cluster,Treatment), color=BrainRegion)) +
  ggtitle("2xLPS mouse dataset: K-means clusters") +
  labs(fill="Treatment") +
  theme_bw(base_size=14) +
  theme(axis.text.x=element_text(angle=45, vjust=1, hjust=1))
# ---- OLDCODE ----
# #ramified group
# View(workingdata01_wcluster)
# 
# ramified_workingdata <- workingdat01_wcluster %>% filter(Cluster=="1")
# working_data01_log_ramified <- transform_log(ramified_workingdata,1, start=8, end=34) 
# pcadata_elbow(working_data01_log_ramified, featurestart=8, featureend=34)
# pca_data01_log_ramified <- pcadata(working_data01_log_ramified, featurestart=8, featureend=34,
#                           pc.start=1, pc.end=10)
# pca_data01_log_transformed_ramified <- transform_scale(pca_data01_log_ramified, start=1, end=3) # scale pca data as input for k-means clustering
# kmeans_input01_log_ramified <- pca_data01_log_transformed_ramified[1:3]
# # check for optimal number of clusters using wss and silhouette methods
# sampling01_log_ramified <- kmeans_input01_log_ramified[sample(nrow(kmeans_input01_log_ramified), 5000),] #sample 5000 random rows for cluster optimization
# fviz_nbclust(sampling01_log_ramified, kmeans, method = 'wss', nstart=25, iter.max=15) # 4 clusters
# fviz_nbclust(sampling01_log_ramified, kmeans, method = 'silhouette', nstart=25, iter.max=15)
# 
# #4clusters
# View(ramified_workingdata)
# View(data_kmeans01_log_ramified)
# data_kmeans01_log_ramified <- kmeans(kmeans_input01_log_ramified, centers=4)
# pca_kmeans01_log_ramified <- cbind(ramified_workingdata[1:34], as.data.frame(data_kmeans01_log_ramified$cluster)) %>%
#   rename(Cluster=`data_kmeans01_log_ramified$cluster`) 
# clusterfeatures(pca_kmeans01_log_ramified, featurestart=10, featureend=34)

# # prepare data for downstream analysis
# workingdata01_wcluster<- cbind(working_data01, as.data.frame(data_kmeans02_log$cluster)) %>%
#   rename(Cluster=`data_kmeans02_log$cluster`) 
# View(workingdata01_wcluster)
# workingdata01_ramified <- workingdata01_wcluster %>% filter(Cluster=="1")
# View(workingdata01_ramified)
# 
# data_forstats <- workingdata01_ramified %>% 
#   group_by(MouseID, Sex, Treatment, BrainRegion) %>% 
#   summarise(across("Foreground pixels":"Maximum branch length", ~mean(.x))) %>% 
#   gather(Measure, Value, "Foreground pixels":"Maximum branch length")
# 
# View(data_forstats)
# 
# # filter out data you want to run stats on and make sure to make any variables included in model into factors
# stats_input_individual <- data_forstats 
# stats_input_individual$Treatment <- factor(stats_input_individual $Treatment)
# stats_input_individual$Sex <- factor(stats_input_individual$Sex)
# stats_input_individual$BrainRegion <- factor(stats_input_individual$BrainRegion)
# 
# # run stats analysis for changes in individual morphology measures
# # you can specify up to two posthoc comparisons (posthoc1 and posthoc2 arguments) - if you only have one set of posthocs to run, specify the same comparison twice for both arguments. you will just get the same results in output[[2]] and output[[3]].
# stats_testing_individual <- stats_morphologymeasures.animal(data = stats_input_individual,
#                                                  model = "Value ~ Treatment*Sex", type="lm",
#                                                  posthoc1 = "~Treatment|Sex", 
#                                                  posthoc2 = "~Treatment*Sex")
# stats_testing_individual[[1]]
# stats_testing_individual[[2]]
# stats_testing_individual[[3]]


# #Log/Scale-Optimization
# working_data02_log <- transform_log(working_data01,1, start=8, end=34) 
# working_data01_log_scale <- transform_scale(working_data02_log, start=8, end=34) 
# pcadata_elbow(working_data01_log_scale, featurestart=8, featureend=34)
# pca_data01_logscale <- pcadata(working_data01_log_scale, featurestart=8, featureend=34,
#                                pc.start=1, pc.end=10)
# pca_data01_logscale_transformed <- transform_scale(pca_data01_logscale, start=1, end=3) # scale pca data as input for k-means clustering
# kmeans_input01_logscale <- pca_data01_logscale_transformed[1:3]
# sampling01_logscale <- kmeans_input01_logscale[sample(nrow(kmeans_input01_logscale), 5000),] #sample 5000 random rows for cluster optimization
# fviz_nbclust(sampling01_logscale, kmeans, method = 'wss', nstart=25, iter.max=15) # 4 clusters
# # fviz_nbclust(sampling_03, kmeans, method = 'silhouette', nstart=25, iter.max=15)
# data_kmeans01_logscale <- kmeans(kmeans_input01_logscale, centers=4)
# pca_kmeans01_logscale <- cbind(pca_data01_logscale[1:3], working_data01, as.data.frame(data_kmeans01_logscale$cluster)) %>%
#   rename(Cluster=`data_kmeans01_logscale$cluster`) 
# clusterfeatures(pca_kmeans01_logscale, featurestart=10, featureend=36)

# arc_test_dir <- "/Users/alexlawson/Library/CloudStorage/OneDrive-UBC/AlexDataAnalysis/EditedBinary/GliaData/arc-skeleton-data"
# arc_skeleton_test <- skeleton_tidying(arc_test_dir)
# View(arc_skeleton_test)
# arc_frac_dir <- "/Users/alexlawson/Library/CloudStorage/OneDrive-UBC/AlexDataAnalysis/EditedBinary/GliaData/avp-frac-lac/20241210081935"
# fraclac_arc <- fraclac_tidying(arc_frac_dir) 
# 
# pvn_test_dir <- "/Users/alexlawson/Library/CloudStorage/OneDrive-UBC/AlexDataAnalysis/EditedBinary/GliaData/pvn-skeleton-data"
# pvn_skeleton_test <- skeleton_tidying(pvn_test_dir)
# View(pvn_skeleton_test)
# pvn_frac_dir <- "/Users/alexlawson/Library/CloudStorage/OneDrive-UBC/AlexDataAnalysis/EditedBinary/GliaData/pvn-fraclac/20241210074205"
# fraclac_pvn <- fraclac_tidying(pvn_frac_dir) 
# 
# hypo_dir <- "/Users/alexlawson/Library/CloudStorage/OneDrive-UBC/AlexDataAnalysis/EditedBinary/GliaData/hypo-skeleton-data"
# hypo_skeleton <- skeleton_tidying(hypo_dir)
# View(hypo_skeleton)
# hypo_frac_dir <- "/Users/alexlawson/Library/CloudStorage/OneDrive-UBC/AlexDataAnalysis/EditedBinary/GliaData/hypo-fraclac/20241210051054"
# fraclac_hypo <- fraclac_tidying(hypo_frac_dir) 
# fraclac_hypo_filtered<- fraclac_hypo[fraclac_hypo$ID %in% hypo_skeleton_test$ID, ]
# View(fraclac_hypo)
# View(fraclac_hypo_filtered)


#  /Users/alexlawson/Library/CloudStorage/OneDrive-UBC/AlexDataAnalysis/EditedBinary/GliaData/hypo-skeleton-data
#/Users/alexlawson/Library/CloudStorage/OneDrive-UBC/AlexDataAnalysis/EditedBinary/GliaData/pvn-skeleton-data
# fraclac_hypothalamus_dir <- "/Users/alexlawson/Library/CloudStorage/OneDrive-UBC/AlexDataAnalysis/EditedBinary/GliaData/fraclac/hypo-frac-lac"
# skeleton_hypothalamus_dir <- "/Users/alexlawson/Library/CloudStorage/OneDrive-UBC/AlexDataAnalysis/EditedBinary/GliaData/hypothalamus-skelton-data"
# 
# fraclac_pvn_dir <- "/Users/alexlawson/Library/CloudStorage/OneDrive-UBC/AlexDataAnalysis/EditedBinary/GliaData/fraclac/pvn-frac-lac"
# skeleton_pvn_dir <- "/Users/alexlawson/Library/CloudStorage/OneDrive-UBC/AlexDataAnalysis/EditedBinary/GliaData/pvn-skeleton-data"
# 
# fraclac_arc_dir <- "/Users/alexlawson/Library/CloudStorage/OneDrive-UBC/AlexDataAnalysis/EditedBinary/GliaData/fraclac/arc-frac-lac"
# skeleton_arc_dir <- "/Users/alexlawson/Library/CloudStorage/OneDrive-UBC/AlexDataAnalysis/EditedBinary/GliaData/arc-skeleton-data"
# 
# # these steps may be very time-intensive, depending on how many cells you are analyzing (i.e., on the order of 1000s of cells). 
# fraclac_hypothalamus <- fraclac_tidying(fraclac_hypothalamus_dir) 
# skeleton_hypothalamus <- skeleton_tidying(skeleton_hypothalamus_dir)
# 
# fraclac_pvn <- fraclac_tidying(fraclac_pvn_dir) 
# skeleton_pvn <- skeleton_tidying(skeleton_pvn_dir)
# 
# fraclac_arc <- fraclac_tidying(fraclac_arc_dir)
# skeleton_arc <- skeleton_tidying(skeleton_arc_dir)
# 
# data_hypothalamus <- merge_data(fraclac_hypothalamus, skeleton_hypothalamus)
# data_pvn <- merge_data(fraclac_pvn, skeleton_pvn)
# data_arc <- merge_data(fraclac_arc, skeleton_arc)
# 
# finaldata_hypothalamus <- metadata_columns(data_hypothalamus, c("BrainRegion","MouseID","Treatment","Sex","SlideNumber"), sep="_")
# finaldata_hypothalamus <- finaldata_hypothalamus[finaldata_hypothalamus$Area >= 10, ]
# finaldata_pvn <- metadata_columns(data_pvn, c("BrainRegion","MouseID","Treatment","Sex","SlideNumber"), sep="_")
# finaldata_pvn <- finaldata_pvn[finaldata_pvn$Area >= 10, ]
# finaldata_arc <- metadata_columns(data_arc, c("BrainRegion","MouseID","Treatment","Sex","SlideNumber"), sep="_")
# finaldata_arc <- finaldata_arc[finaldata_arc$Area >= 10, ]
# 
# combined_data <- rbind(finaldata_hypothalamus, finaldata_pvn, finaldata_arc)
# 
# write.csv(finaldata_hypothalamus, "/Users/alexlawson/Documents/GitHub/MorphologyAnalysis01/data-files/hypothalamus-data.csv", row.names = FALSE)
# write.csv(finaldata_pvn, "/Users/alexlawson/Documents/GitHub/MorphologyAnalysis01/data-files/pvn-data.csv", row.names = FALSE)
# write.csv(finaldata_arc, "/Users/alexlawson/Documents/GitHub/MorphologyAnalysis01/data-files/arc-data.csv", row.names = FALSE)
# write.csv(combined_data, "/Users/alexlawson/Desktop/GliaData/frac-lac-merged", row.names = FALSE)
# 
# hypo_data <- read.csv("/Users/alexlawson/Documents/GitHub/MorphologyAnalysis01/data-files/hypothalamus-data.csv")
# pvn_data <- read.csv("/Users/alexlawson/Documents/GitHub/MorphologyAnalysis01/data-files/pvn-data.csv")
# arc_data <- read.csv("/Users/alexlawson/Documents/GitHub/MorphologyAnalysis01/data-files/arc-data.csv")
# combined_data <- rbind(hypo_data, pvn_data, arc_data)

# gather your numerical morphology data into one column ('measure') which contains the feature name, and another column ('value') which contains measured values
# combined_data_gathered <- combined_data %>% gather(measure, value, 8:ncol(combined_data))
# featurecorrelations(combined_data,
#                     featurestart=8, featureend=34,
#                     rthresh=0.8, pthresh=0.05,
#                     title="Correlations across features")
# correlationstats <- featurecorrelations_stats(combined_data,
#                                               featurestart=8, featureend=34,
#                                               rthresh=0.8, pthresh=0.05)
# View(correlationstats)
# #check for outliers
# outliers_boxplots(combined_data_gathered)
# outliers_distributions(combined_data_gathered)
# normalize_logplots(combined_data_gathered,1)
# normalize_minmax(combined_data_gathered)
# normalize_scaled(combined_data_gathered)
# 
# combined_data_log <- transform_log(combined_data, start=8, end=34, 1)
# #pcadata_elbow(combined_data_log, featurestart=8, featureend=34)
# #pcfeaturecorrelations(pca_data_combined, pc.start=1, pc.end=5, 
# #                      feature.start=18, feature.end=44, 
# #                      rthresh=0.75, pthresh=0.05, 
# #                      title="Correlation between PCs and features")
# pca_data_combined <- pcadata(combined_data_log, featurestart=8, featureend=34,
#                              pc.start=1, pc.end=10)
# pca_data_scale_combined <- transform_scale(pca_data_combined, start=1, end=10)
# kmeans_input_combined <- pca_data_scale_combined[1:3]
# sampling <- kmeans_input_combined[sample(nrow(kmeans_input_combined), 5000),]
# 
# fviz_nbclust(sampling, kmeans, method = 'silhouette', nstart=50, iter.max=20)
# fviz_nbclust(sampling, kmeans, method = 'wss', nstart=50, iter.max=20)
# data_kmeans_combined <- kmeans(kmeans_input_combined, centers=4)
# pca_kmeans_combined <- cbind(pca_data_combined[1:2], combined_data, as.data.frame(data_kmeans_combined$cluster)) %>%
#   rename(Cluster=`data_kmeans_combined$cluster`) 
# 
# # gather your data by experimental variables (e.g., Treatment, Sex, MouseID, etc.)
# gathered_expvariables <- pca_data_combined %>% gather(variable, value, 11:16) 
# 
# plots_expvariable(gathered_expvariables, "PC1", "PC2")
# 
# 
# plot <- clusterplots(pca_kmeans_combined, "PC1", "PC2")
# plot
# clusterfeatures(pca_kmeans_combined, featurestart=10, featureend=36)
# cp_combined <- clusterpercentage(pca_kmeans_combined, "Cluster", MouseID, BrainRegion, Treatment, Sex)
# 
# 
# colorbycluster <- pca_kmeans_combined %>% 
#   filter(BrainRegion=="HYPO",MouseID=="MR169", SlideNumber=="s03") %>% select(c(Cluster, ID))
# write.csv(colorbycluster, "/Users/alexlawson/Library/CloudStorage/OneDrive-UBC/AlexDataAnalysis/color-by-cluster/HYPO_MR169_s3_4pca.csv", row.names = FALSE)
# 
# 
# colorbycluster <- pca_kmeans_combined %>% 
#   filter(BrainRegion=="HYPO",MouseID=="MR290", SlideNumber=="s02") %>% select(c(Cluster, ID))
# write.csv(colorbycluster, "/Users/alexlawson/Library/CloudStorage/OneDrive-UBC/AlexDataAnalysis/color-by-cluster/HYPO_MR290_s02.csv", row.names = FALSE)
# 
# cp_combined <- clusterpercentage(pca_kmeans_combined, "Cluster", MouseID, BrainRegion, Treatment, Sex)
# # prepare percentages dataset for downstream analysis
# stats_input <- cp_combined_male
# stats_input$MouseID <- factor(stats_input$MouseID)
# stats_input$Cluster <- factor(stats_input$Cluster)
# stats_input$Treatment <- factor(stats_input$Treatment)
# # # 
#  # run stats analysis for changes in cluster percentages, at the animal level
# # you can specify up to two posthoc comparisons (posthoc1 and posthoc2 arguments) - if you only have one set of posthocs to run, specify the same comparison twice for both arguments. you will just get the same results in output[[2]] and output[[3]].
# stats_testing_new <- stats_cluster.animal(data = stats_input, 
#                                       model = "percentage ~ Cluster*BrainRegion*Treatment + (1|MouseID)", 
#                                       posthoc1 = "~Treatment|Cluster|BrainRegion", 
#                                       posthoc2 = "~Cluster|Treatment", adjust = "holm")
# # # 
# stats_testing_new[[1]]
# stats_testing_new[[2]]
# stats_testing_new[[3]]
# 
# View(combined_data)
# 
# morphology_data_grouped <- combined_data %>% 
#   group_by(MouseID, Sex, Treatment, BrainRegion) %>% 
#   summarise(across("Foreground.pixels":"Maximum.branch.length", ~mean(.x))) %>% 
#   gather(Measure, Value, "Foreground.pixels":"Maximum.branch.length")
# 
# morphology_stats_input <- morphology_data_grouped 
# morphology_stats_input$Treatment <- factor(morphology_stats_input$Treatment)
# 
# 
# colorbycluster_MR290_HYPO_S09 <- pca_kmeans_combined %>% 
#   filter(MouseID=="MR290", BrainRegion=="HYPO", SlideNumber=="s09") %>% select(c(Cluster, ID))
# head(colorbycluster_MR290_HYPO_S09)
# 
# write.csv(colorbycluster_MR290_HYPO_S09, "/Users/alexlawson/Library/CloudStorage/OneDrive-UBC/AlexDataAnalysis/color-by-cluster/MR290_HYPO_S09.csv")
# 
# colorbycluster_MR290_HYPO_S08 <- pca_kmeans_combined %>% 
#   filter(MouseID=="MR290", BrainRegion=="HYPO", SlideNumber=="s08") %>% select(c(Cluster, ID))
# 
# View(colorbycluster_MR290_HYPO_S08)
# 
# 
# write.csv(colorbycluster_MR290_HYPO_S08, "/Users/alexlawson/Library/CloudStorage/OneDrive-UBC/AlexDataAnalysis/color-by-cluster/MR290_HYPO_S08.csv")
# 
# arc_data_gathered <- arc_data %>% gather(measure, value, 8:ncol(arc_data))
# arc_data_gathered_logtransformed <- transform_log(arc_data, 1, start=8, end=34)
# pcadata_elbow(arc_data_gathered_logtransformed, featurestart=8, featureend=34)
# pca_data_arc <- pcadata(arc_data_gathered_logtransformed, featurestart=8, featureend=34,
#                              pc.start=1, pc.end=10)
# 
# 
# pcfeaturecorrelations(pca_data_arc, pc.start=1, pc.end=4, 
#                       feature.start=18, feature.end=44, 
#                       rthresh=0.75, pthresh=0.05, 
#                       title="Correlation between PCs and features")
# 
# pca_data_scale_arc <- transform_scale(pca_data_arc, start=1, end=4)
# kmeans_input_arc <- pca_data_scale_arc[1:4]
# sampling_arc <- kmeans_input_arc[sample(nrow(kmeans_input_arc), 1000),]
# 
# fviz_nbclust(sampling_arc, kmeans, method = 'silhouette', nstart=25, iter.max=50)
# fviz_nbclust(sampling, kmeans, method = 'wss', nstart=25, iter.max=50)
# data_kmeans_arc <- kmeans(kmeans_input_arc, centers=5)
# pca_kmeans_arc <- cbind(pca_data_arc[1:2], arc_data, as.data.frame(data_kmeans_arc$cluster)) %>%
#   rename(Cluster=`data_kmeans_arc$cluster`) 
# plot <- clusterplots(pca_kmeans_arc, "PC1", "PC2")
# plot
# clusterfeatures(pca_kmeans_arc, featurestart=10, featureend=36)
# 
# cp_arc <- clusterpercentage(pca_kmeans_arc, "Cluster", MouseID, BrainRegion, Treatment, Sex)
# # prepare percentages dataset for downstream analysis
# stats_input_arc <- cp_arc
# stats_input_arc$MouseID <- factor(stats_input_arc$MouseID)
# stats_input_arc$Cluster <- factor(stats_input_arc$Cluster)
# stats_input_arc$Treatment <- factor(stats_input_arc$Treatment)
# 
# # run stats analysis for changes in cluster percentages, at the animal level
# # you can specify up to two posthoc comparisons (posthoc1 and posthoc2 arguments) - if you only have one set of posthocs to run, specify the same comparison twice for both arguments. you will just get the same results in output[[2]] and output[[3]].
# stats_testing_arc <- stats_cluster.animal(data = stats_input_arc, 
#                                       model = "percentage ~ Cluster*Treatment + Sex + (1|MouseID)", 
#                                       posthoc1 = "~Treatment|Cluster", 
#                                       posthoc2 = "~Treatment|Cluster|Sex")
# stats_testing_arc[[1]]
# stats_testing_arc[[2]]
# stats_testing_arc[[3]]
# 
# 
# filtered_combined_data <- combined_data %>%
#   select(-Width.of.bounding.rectangle, -Height.of.bounding.rectangle, -Perimeter, -Mean.radius, -X..of.slab.voxels)
# filtered_combined_data_gathered <- filtered_combined_data %>% gather(measure, value, 8:ncol(filtered_combined_data))
# normalize_logplots(filtered_combined_data_gathered,1)
# normalize_minmax(filtered_combined_data_gathered)
# normalize_scaled(filtered_combined_data_gathered)
# filtered_combined_data_minmax <- transform_minmax(filtered_combined_data, start=8, end=29)
# filtered_pca_data_combined <- pcadata(filtered_combined_data_minmax, featurestart=8, featureend=29,
#                                       pc.start=1, pc.end=10)
# pcadata_elbow(filtered_combined_data_scale, featurestart=8, featureend=29)
# filtered_pca_data_scale_combined <- transform_scale(filtered_pca_data_combined, start=1, end=10)
# filtered_kmeans_input_combined <- filtered_pca_data_scale_combined[1:4]
# filtered_sampling <- filtered_kmeans_input_combined[sample(nrow(filtered_kmeans_input_combined), 5000),]
# fviz_nbclust(filtered_sampling, kmeans, method = 'silhouette', nstart=50, iter.max=20)
# fviz_nbclust(filtered_sampling, kmeans, method = 'wss', nstart=50, iter.max=20)
# filtered_kmeans_input_combined <- filtered_pca_data_scale_combined[1:4]
# filtered_data_kmeans_combined <- kmeans(filtered_kmeans_input_combined, centers=5)
# filtered_pca_kmeans <- cbind(combined_data, as.data.frame(filtered_data_kmeans_combined$cluster)) %>%
#   rename(Cluster=`filtered_data_kmeans_combined$cluster`)
# clusterfeatures(filtered_pca_kmeans, featurestart=8, featureend=29)
