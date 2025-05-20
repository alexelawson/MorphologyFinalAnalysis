library(dplyr)
library(Seurat)
library(patchwork)
library(clustree)

#use markers from paper to cluster the microglia 
#go over representation for pathway 
#everything with xist is female and everything else is male 

cold_data_dir <- "/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/COLD"
con_data_dir <- "/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/CON"

# Load Control
control_counts <- Read10X(data.dir = con_data_dir)
#Load Cold Stress
cold_counts <- Read10X(data.dir = cold_data_dir)

control <- CreateSeuratObject(counts = control_counts, project = "Control", min.cells = 3, min.features = 200)
cold <- CreateSeuratObject(counts = cold_counts, project = "COLD", min.cells = 3, min.features = 200)

control[["percent.mt"]] <- PercentageFeatureSet(control, pattern = "^mt-")
cold[["percent.mt"]]    <- PercentageFeatureSet(cold,    pattern = "^mt-")
#
control <- subset(control, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)	
cold <- subset(cold, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)	

# Tag each object with a condition label
control$condition <- "CON"
cold$condition    <- "COLD"

# 2. Put them in a list and preprocess each
objs <- list(control, cold)

for (i in seq_along(objs)) {
 objs[[i]] <- NormalizeData(objs[[i]], verbose = FALSE)
 objs[[i]] <- FindVariableFeatures(objs[[i]], selection.method = "vst")
}


# 3. Find the anchors (matching cells across conditions)
anchors <- FindIntegrationAnchors(object.list = objs,
                                  dims = 1:20)        # match along first 20 PCs

# 4. Integrate the data
combined <- IntegrateData(anchorset = anchors,
                          dims = 1:20)           # again, use 20 dims


# 5. Switch to the new “integrated” assay, then run your standard workflow
DefaultAssay(combined) <- "integrated"
# # Merge into one Seurat object
# combined <- merge(control,
#                   y = cold,
#                   add.cell.ids = c("CON","COLD"),
#                   project = "Combined")

# # Then run the standard workflow on combined:
# combined <- NormalizeData(combined)
# combined <- FindVariableFeatures(combined, selection.method="vst")

all_genes_combined <- rownames(combined)
combined <- ScaleData(combined, npcs = 20, features = all_genes_combined)
combined <- RunPCA(combined, features = VariableFeatures(object = combined))
ElbowPlot(combined, ndims=40, reduction = "pca")
combined <- FindNeighbors(combined, dims=1:20)
combined <- FindClusters(combined, resolution=0.5)
combined <- RunUMAP(combined, dims=1:20)
umap_combined <- DimPlot(combined, reduction="umap", group.by="seurat_clusters")
umap_combined
ggsave(filename = '/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/PLOTS/UMAP_combined.png', umap_combined, width = 10, height = 5)


#-------macrophage cell markers-------
plot_macrophage_feature <- FeaturePlot(combined, features = c("Ifitm3",
                                                "S100a6",
                                                "Lgals3",
                                                "Ifi27l2a",
                                                "S100a11",
                                                "Plac8",
                                                "S100a10",
                                                "Ccr2",
                                                "Ifitm6",
                                                "Cybb"))

ggsave(filename = '/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/PLOTS/macrofage_features.png', plot_macrophage_feature , width = 10, height = 5)

plot_macrophage_vln <- VlnPlot(combined, features = c("Ifitm3",
                                                                   "S100a6",
                                                                   "Lgals3",
                                                                   "Ifi27l2a",
                                                                   "S100a11",
                                                                   "Plac8",
                                                                   "S100a10",
                                                                   "Ccr2",
                                                                   "Ifitm6",
                                                                   "Cybb"))

ggsave(filename = '/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/PLOTS/macrofage_vln.png', plot_macrophage_vln, width = 10, height = 5)

#Macrophage = Cluster 10

#-------microglia 1 cell markers-------

plot_cluster1_feature <- FeaturePlot(combined, features = c("Pf4",
                                                "F13a1",
                                                "Mrc1",
                                                "Lyve1",
                                                "Dab2",
                                                "Ccl24",
                                                "Ms4a7",
                                                "Ms4a6c",
                                                "Blvrb",
                                                "Fcgrt"))

ggsave(filename = '/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/PLOTS/clust01_features.png', plot_cluster1_feature, width = 10, height = 5)


plot_cluster1_vln <- VlnPlot(combined, features = c("Pf4",
                                            "F13a1",
                                            "Mrc1",
                                            "Lyve1",
                                            "Dab2",
                                            "Ccl24",
                                            "Ms4a7",
                                            "Ms4a6c",
                                            "Blvrb",
                                            "Fcgrt"))

ggsave(filename = '/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/PLOTS/clust01_vln.png', plot_cluster1_vln, width = 10, height = 5)

#Microglia 1 - Cluster 8

#-------microglia 2 cell markers----
plot_cluster2_features <- FeaturePlot(combined, features = c("Spp1",
                                                "Fabp5",
                                                "Ccl4",
                                                "Ccl3",
                                                "Csf1",
                                                "Lpl",
                                                "Ctsb",
                                                "Fam20c",
                                                "Mif",
                                                "Fabp3"))

ggsave(filename = '/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/PLOTS/clust02_features.png', plot_cluster2_features, width = 10, height = 5)


plot_cluster2_vln <- VlnPlot(combined, features = c("Spp1",
                                            "Fabp5",
                                            "Ccl4",
                                            "Ccl3",
                                            "Csf1",
                                            "Lpl",
                                            "Ctsb",
                                            "Fam20c",
                                            "Mif",
                                            "Fabp3"))

ggsave(filename = '/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/PLOTS/clust02_vln.png', plot_cluster2_vln, width = 10, height = 5)
plot_cluster2_vln
#microglia 2 is cluster 1

#-------microglia 3 cell markers----
plot_cluster3_features <- FeaturePlot(combined, features = c("Ube2c",
                                                "Hist1h2ap",
                                                "Birc5",
                                                "Hmgb2",
                                                "H2afx",
                                                "Cenpa",
                                                "Mki67",
                                                "Pclaf",
                                                "Top2a",
                                                "Pttg1"))

ggsave(filename = '/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/PLOTS/clust03_features.png', plot_cluster3_features, width = 10, height = 5)


plot_cluster_3_vln <- VlnPlot(combined, features = c("Ube2c",
                                            "Hist1h2ap",
                                            "Birc5",
                                            "Hmgb2",
                                            "H2afx",
                                            "Cenpa",
                                            "Mki67",
                                            "Pclaf",
                                            "Top2a",
                                            "Pttg1"))

ggsave(filename = '/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/PLOTS/clust03_vln.png', plot_cluster_3_vln, width = 10, height = 5)

plot_cluster_3_vln
plot_cluster_3_vln
#cluster 3

#-------microglia 4 cell markers-----
plot_cluster4_feature <- FeaturePlot(combined, features = c("Crybb1",
                                                "P2ry12",
                                                "Sparc",
                                                "Ccr5",
                                                "Tgfbr1",
                                                "Lpcat2",
                                                "Siglech",
                                                "C1qa",
                                                "Tmem119",
                                                "Sall1"))

ggsave(filename = '/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/PLOTS/clust04_features.png', plot_cluster4_feature, width = 10, height = 5)


plot_cluster4_vln <- VlnPlot(combined, features = c("Crybb1",
                                            "P2ry12",
                                            "Sparc",
                                            "Ccr5",
                                            "Tgfbr1",
                                            "Lpcat2",
                                            "Siglech",
                                            "C1qa",
                                            "Tmem119",
                                            "Sall1"))

plot_cluster4_vln

ggsave(filename = '/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/PLOTS/clust04_vln.png', plot_cluster4_vln, width = 10, height = 5)

#cluster 5, 2, 0, 7

# prolif cell markers ----------------------------------------------

plot_prolif_feature <- FeaturePlot(combined, features = c("Pif1",
                                                "Sox2"))
plot_prolif_vln <- VlnPlot(combined, features = c("Pif1",
                                            "Sox2"))
plot_prolif_vln

#OUTPUT CLUSTERS
#Macrophage: Cluster 10
#Microglia 1: Cluster 8 
#Microglia 2: +Cluster 1
#Microglia 3: +Cluster 3
#Microglia 4: Cluster 5, +Cluster 2, +Cluster 0, Cluster 7 
#Proliferative: Cluster 4, Cluster 6


new.cluster.ids <- c("Microglia 4", "Microglia 2", "Microglia 4", "Microglia 3", "Proliferative", "Microglia 4", "Proliferative", "Microglia 4", "Microglia 1", "Pericyte", "Macrophage")
names(new.cluster.ids) <- levels(combined)
combined <- RenameIdents(combined, new.cluster.ids)
DimPlot(combined, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
combined$renamed_clusters <- Idents(combined)
table(Idents(combined)) 
Idents(combined) <- "renamed_clusters"

View(combined)

keep <- paste0("Microglia ", 1:4)
combined_clean <- subset(combined, idents = keep)

control <- subset(combined_clean, subset = condition == "CON")
cold <- subset(combined_clean, subset = condition == "COLD")

# Define sex marker gene for female classification
female_genes <- c("Xist")
expr_threshold <- 0  # Any non-zero expression will count

# Re-pull your expression matrices (so you haven’t overwritten them)
DefaultAssay(control) <- "RNA"
expr_control <- GetAssayData(control, assay="RNA", slot="data")

DefaultAssay(cold) <- "RNA"
expr_cold    <- GetAssayData(cold,   assay="RNA", slot="data")

# Identify female cells (any cell with non-zero expression of Xist)
female_cells_control <- colnames(expr_control)[colSums(expr_control[female_genes, , drop = FALSE] > expr_threshold) > 0]

# Identify male cells (remaining cells that do not meet the female criteria)
male_cells_control <- setdiff(colnames(expr_control), female_cells_control)

# Identify female cells (any cell with non-zero expression of Xist)
female_cells_cold <- colnames(expr_cold)[colSums(expr_cold[female_genes, , drop = FALSE] > expr_threshold) > 0]

# Identify male cells (remaining cells that do not meet the female criteria)
male_cells_cold <- setdiff(colnames(expr_cold), female_cells_cold)

# Optional: confirm number of cells
length(female_cells_control) # 269
length(male_cells_control)   # 557

# Optional: confirm number of cells
length(female_cells_cold) # 1113
length(male_cells_cold)   # 90


# Subset Seurat objects for female and male cells
female_cells_control_subset <- subset(combined_clean, cells = female_cells_control)
male_cells_control_subset  <- subset(combined_clean, cells = male_cells_control)

# Subset Seurat objects for female and male cells
female_cells_cold_subset <- subset(combined_clean, cells = female_cells_cold)
male_cells_cold_subset  <- subset(combined_clean, cells = male_cells_cold)

# Merge the male and female subsets (keeping original names intact)
combined_control_subset <- merge(female_cells_control_subset, y = male_cells_control_subset, add.cell.ids = c("Female", "Male"), project = "Microglia_CON")
# Add metadata for sex
combined_control_subset$sex <- ifelse(grepl("^Female", colnames(combined_control_subset)), "Female", "Male")

# Merge the male and female subsets (keeping original names intact)
combined_cold_subset <- merge(female_cells_cold_subset, y = male_cells_cold_subset, add.cell.ids = c("Female", "Male"), project = "Microglia_CON")
# Add metadata for sex
combined_cold_subset$sex <- ifelse(grepl("^Female", colnames(combined_cold_subset)), "Female", "Male")

combined_all_microglia <- merge(combined_cold_subset, y = combined_control_subset, add.cell.ids = c("Female", "Male"), project = "Microglia_CON") 

View(combined_all_microglia)

# Normalize and scale (post-merge)
DefaultAssay(combined_all_microglia) <- "RNA"
combined_all_microglia <- NormalizeData(combined_all_microglia)
combined_all_microglia <- FindVariableFeatures(combined_all_microglia)
combined_all_microglia <- ScaleData(combined_all_microglia)

combined_female <- subset(
  combined_all_microglia,
  subset = sex == "Female"
)

combined_male <- subset(
  combined_all_microglia,
  subset = sex == "Male"
)


# Set the identities for DE analysis
Idents(combined_female) <- "condition"
# JoinLayers to properly link for analysis 
combined_female <- JoinLayers(combined_female)
# Differential expression analysis
deg_combined_female <- FindMarkers(
 combined_female,
  ident.1 = "CON",
  ident.2 = "COLD",
  assay = "RNA",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
# View and save results
head(deg_combined_female)
write.csv(deg_combined_female, "/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/Data Tables/female_cold_vs_con.csv")

# Set the identities for DE analysis
Idents(combined_male) <- "condition"
# JoinLayers to properly link for analysis 
combined_male <- JoinLayers(combined_male)
# Differential expression analysis
deg_combined_male <- FindMarkers(
  combined_male,
  ident.1 = "CON",
  ident.2 = "COLD",
  assay = "RNA",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
# View and save results
head(deg_combined_male)
write.csv(deg_combined_male, "/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/Data Tables/male_cold_vs_con.csv")




