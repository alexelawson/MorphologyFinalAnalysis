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
stats_input_final <- cluster_data
View(data_frame_final)

data_frame_pvn <- data_frame_final%>% filter (data_frame_final$BrainRegion=="PVN")

data_frame_pvn <- data_frame_pvn[, !names(data_frame_pvn) %in% c("PC1", "PC2", "PC3", "Cluster")]

View(data_frame_pvn)
data_frame_pvn_logtransformed <- transform_log(data_frame_pvn, 1, start=8, end=34) 
pcadata_elbow(data_frame_pvn_logtransformed, featurestart=8, featureend=34)
pca_data <- pcadata(data_frame_pvn_logtransformed, featurestart=8, featureend=34,
                    pc.start=1, pc.end=10)
pca_data_scale <- transform_scale(pca_data, start=1, end=3) # scale pca data as input for k-means clustering
kmeans_input <- pca_data_scale[1:3]
sampling <- kmeans_input[sample(nrow(kmeans_input), 200),] #sample 5000 random rows for cluster optimization
fviz_nbclust(sampling, kmeans, method = 'silhouette', nstart=25, iter.max=50) # 4 clusters
data_kmeans <- kmeans(kmeans_input, centers=4)
pca_kmeans <- cbind(pca_data[1:2], data_frame_pvn, as.data.frame(data_kmeans$cluster)) %>%
  rename(Cluster=`data_kmeans$cluster`) 
clusterfeatures(pca_kmeans, featurestart=10, featureend=36)


