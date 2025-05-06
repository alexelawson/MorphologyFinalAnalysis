# List of necessary libraries
libraries <- c("gplots", "ggplot2", "dplyr", "tidyr")
# Install any libraries that are not already installed
for (lib in libraries) {
  if (!requireNamespace(lib, quietly = TRUE)) {
    install.packages(lib)
  }
}

# Install corrplot package if not already installed
if (!requireNamespace("corrplot", quietly = TRUE)) {
  install.packages("corrplot")
}

# Load corrplot library
library(corrplot)

# Load the libraries
lapply(libraries, library, character.only = TRUE)

View(data)
# Load data and prepare it
data <- read.csv("/Users/alexlawson/Masters-Data-Final/Morphoglia/PVN_all/Morphoglia/MorphoGlia/Data/Morphology_HDBSCAN_30_0.1_150_5_reformatted.csv")
data <- data[data$Clusters != -1, ]  # Filter out noise

data_M <- data %>% filter(tissue=="M")
data_F <- data %>% filter(tissue=="F")


data_M$Cluster_Labels <- factor(data_M$Clusters)
data_M$categories <- factor(gsub("_", " ", data_M$treatment)) 

data_F$Cluster_Labels <- factor(data_F$Clusters)
data_F$categories <- factor(gsub("_", " ", data_F$treatment))


# Create and print the contingency table
contingency_table_M <- table(data_M$treatment, data_M$Clusters)
print(contingency_table_M)

contingency_table_F <- table(data_F$treatment, data_F$Clusters)
print(contingency_table_F)


# Chi-square test
chisq_results_M <- chisq.test(contingency_table_M)
print(chisq_results_M)

# Chi-square test
chisq_results_F <- chisq.test(contingency_table_F)
print(chisq_results_F)

Plot_path = "/Users/alexlawson/Masters-Data-Final/Morphoglia/PVN_all/Morphoglia/MorphoGlia/Data/stats"
# Visualization: Balloon Plot
png(filename = file.path(Plot_path, "Cluster_Frequencies-M-reformatted.png"), width = 800, height = 600)
balloonplot(t(contingency_table_M), main = "Frequencies", xlab = "Clusters", ylab = "Experimental Group",
            label = TRUE, show.margins = FALSE, dotcolor = "lightblue", text.size = 1.5, font = 2)
dev.off()

library(corrplot)
library(RColorBrewer)

custom_palette <- colorRampPalette(c(
  "#2166AC", "#92C5DE",  # negative residuals (pale blue)
  "#F0F0F0", "#F7F7F7", "#FFFFFF", "#F7F7F7",  # widened white zone from ~–0.5 to 0.5
  "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"  # positive residuals (strong red)
))(100)

# Save high-res PNG with outlined boxes
png(filename = file.path(Plot_path, "Chi_Square_Correlogram-F-reformatted.png"),
    width = 1600, height = 1400, res = 300)

corrplot(
  chisq_results_F$residuals,
  is.cor = FALSE,
  method = "color",
  col = custom_palette,
  tl.col = "black",
  tl.cex = 1.4,
  font = 2,
  number.cex = 1,
  addCoef.col = "black",
  cl.pos = "r",
  cl.cex = 0.8,
  cl.ratio = 0.2,
  tl.srt = 45,
  zlim = c(-1.67, 1.67), 
  na.label = " ",
  mar = c(1, 1, 1, 5),
  addgrid.col = "black"  # ← outlines each cell
)

dev.off()


# Contribution Plot
contrib_M <- 100 * contingency_table_F^2 / sum(contingency_table_F)
contrib_M <- round(contrib_M, 3)
png(filename = file.path(Plot_path, "Contribution_Plot-F-v2.png"), width = 800, height = 600)
corrplot(contrib_M, is.cor = FALSE, tl.col = "black", tl.cex = 1.2, font = 2,
         addCoef.col = "black", col = colorRampPalette(c("yellow", "red"))(10), cl.pos = "n", cl.cex = 1.8, tl.srt = 0)
dev.off()
chisq.residuals <- chisq.test(contingency_table)$residuals
print(chisq.residuals)
library(vcd)
assocstats(contingency_table)






