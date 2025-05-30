library(dplyr)
library(Seurat)
library(patchwork)
library(clustree)
library(clusterProfiler)
library(org.Mm.eg.db)
library(GOSemSim)
library(pheatmap)
library(enrichplot)

#use markers from paper to cluster the microglia 
#go over representation for pathway 
#everything with xist is female and everything else is male 

set.seed(89)

#setting directories
cold_data_dir <- "/Users/alexlawson/MorphologyFinalAnalysis/rnasec/COLD"
con_data_dir <- "/Users/alexlawson/MorphologyFinalAnalysis/rnasec/CON"

# Load Control
control_counts <- Read10X(data.dir = con_data_dir)
#Load Cold Stress
cold_counts <- Read10X(data.dir = cold_data_dir)

#creating seurat objects
control <- CreateSeuratObject(counts = control_counts, project = "Control", min.cells = 3, min.features = 200)
cold <- CreateSeuratObject(counts = cold_counts, project = "COLD", min.cells = 3, min.features = 200)

#filtering cells
control[["percent.mt"]] <- PercentageFeatureSet(control, pattern = "^mt-")
cold[["percent.mt"]]    <- PercentageFeatureSet(cold,    pattern = "^mt-")
control <- subset(control, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)	
cold <- subset(cold, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)	

#Tagging each object with a condition label 
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
# 5a) Cell‐cycle scoring
s.features   <- cc.genes.updated.2019$s.genes
g2m.features <- cc.genes.updated.2019$g2m.genes

combined <- CellCycleScoring(
  combined,
  s.features   = s.features,
  g2m.features = g2m.features,
  set.ident    = FALSE
)

#Scale while regressing out cell‐cycle scores
all_genes <- VariableFeatures(combined)
combined <- ScaleData(
  combined,
  features       = all_genes,
  verbose        = FALSE
)

combined <- RunPCA(combined, features = VariableFeatures(object = combined))
combined <- FindNeighbors(combined, dims=1:20)
combined <- FindClusters(combined, resolution=0.5)
combined <- RunUMAP(combined, dims=1:20)
umap_combined <- DimPlot(combined, reduction="umap", group.by="seurat_clusters")

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

#ggsave(filename = '/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/PLOTS/macrofage_features.png', plot_macrophage_feature , width = 10, height = 5)

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

#ggsave(filename = '/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/PLOTS/macrofage_vln.png', plot_macrophage_vln, width = 10, height = 5)

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

#ggsave(filename = '/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/PLOTS/clust01_features.png', plot_cluster1_feature, width = 10, height = 5)


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

#ggsave(filename = '/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/PLOTS/clust01_vln.png', plot_cluster1_vln, width = 10, height = 5)

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

#ggsave(filename = '/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/PLOTS/clust02_features.png', plot_cluster2_features, width = 10, height = 5)

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

#ggsave(filename = '/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/PLOTS/clust02_vln.png', plot_cluster2_vln, width = 10, height = 5)
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

#ggsave(filename = '/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/PLOTS/clust03_features.png', plot_cluster3_features, width = 10, height = 5)


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

#ggsave(filename = '/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/PLOTS/clust03_vln.png', plot_cluster_3_vln, width = 10, height = 5)

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

#ggsave(filename = '/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/PLOTS/clust04_features.png', plot_cluster4_feature, width = 10, height = 5)


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


#ggsave(filename = '/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/PLOTS/clust04_vln.png', plot_cluster4_vln, width = 10, height = 5)

#cluster 5, 2, 0, 7

# prolif cell markers ----------------------------------------------

plot_prolif_feature <- FeaturePlot(combined, features = c("Pif1",
                                                          "Sox2"))
plot_prolif_vln <- VlnPlot(combined, features = c("Pif1",
                                                  "Sox2"))
#OUTPUT CLUSTERS
#Macrophage: Cluster 10
#Microglia 1: Cluster 8 
#Microglia 2: +Cluster 1
#Microglia 3: +Cluster 3
#Microglia 4: Cluster 5, +Cluster 2, +Cluster 0, Cluster 7 
#Proliferative: Cluster 4, Cluster 6

# rest of workflow----------------------------------------------
new.cluster.ids <- c("Microglia 4", "Microglia 2", "Microglia 4", "Microglia 3", "Proliferative", "Microglia 4", "Proliferative", "Microglia 4", "Microglia 1", "Pericyte", "Macrophage")
names(new.cluster.ids) <- levels(combined)
combined <- RenameIdents(combined, new.cluster.ids)
DimPlot(combined, reduction = "umap", label = TRUE, pt.size = 0.5)
combined$renamed_clusters <- Idents(combined)
table(Idents(combined)) 
Idents(combined) <- "renamed_clusters"
DimPlot(combined, reduction="umap", group.by="renamed_clusters")

keep <- paste0("Microglia ", 1:4)
combined_clean <- subset(combined, idents = keep)
DimPlot(combined_clean, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
DefaultAssay(combined_clean) <- "RNA"


control <- subset(combined_clean, subset = condition == "CON")
cold <- subset(combined_clean, subset = condition == "COLD")
# Define sex marker gene for female classification
female_genes <- c("Xist")
expr_threshold <- 0  # Any non-zero expression will count

# Re-pull your expression matrices (so you haven’t overwritten them)
DefaultAssay(control) <- "RNA"
expr_control <- GetAssayData(control, assay="RNA", layer="data")

DefaultAssay(cold) <- "RNA"
expr_cold <- GetAssayData(cold,   assay="RNA", layer="data")

# Identify female cells (any cell with non-zero expression of Xist)
female_cells_control <- colnames(expr_control)[colSums(expr_control[female_genes, , drop = FALSE] > expr_threshold) > 0]


# Identify male cells (remaining cells that do not meet the female criteria)
male_cells_control <- setdiff(colnames(expr_control), female_cells_control)

# Identify female cells (any cell with non-zero expression of Xist)
female_cells_cold <- colnames(expr_cold)[colSums(expr_cold[female_genes, , drop = FALSE] > expr_threshold) > 0]

# Identify male cells (remaining cells that do not meet the female criteria)
male_cells_cold <- setdiff(colnames(expr_cold), female_cells_cold)

#confirm number of cells
length(female_cells_control) # 269
length(male_cells_control)   # 557

#confirm number of cells
length(female_cells_cold) # 1113
length(male_cells_cold)   # 90


# Subset Seurat objects for female and male cells
female_cells_control_subset <- subset(combined_clean, cells = female_cells_control)
male_cells_control_subset  <- subset(combined_clean, cells = male_cells_control)

# Subset Seurat objects for female and male cells
female_cells_cold_subset <- subset(combined_clean, cells = female_cells_cold)
male_cells_cold_subset  <- subset(combined_clean, cells = male_cells_cold)

# Control merge
combined_control_subset <- merge(
  female_cells_control_subset, 
  y = male_cells_control_subset, 
  add.cell.ids = c("F_CON","M_CON"), 
  project    = "Microglia_CON"
)

combined_control_subset$sex <- ifelse(
  grepl("^F_CON", colnames(combined_control_subset)), 
  "Female", 
  "Male"
)
combined_control_subset$condition <- "CON"
# Cold merge
combined_cold_subset <- merge(
  female_cells_cold_subset, 
  y = male_cells_cold_subset, 
  add.cell.ids = c("F_COLD","M_COLD"), 
  project    = "Microglia_COLD"
)
combined_cold_subset$sex <- ifelse(
  grepl("^F_COLD", colnames(combined_cold_subset)), 
  "Female", 
  "Male"
)
combined_cold_subset$condition <- "COLD"

# Final merge
combined_all_microglia <- merge(
  combined_cold_subset,
  y = combined_control_subset,
  project    = "Microglia_ALL"
)

# Normalize and scale (post-merge)
DefaultAssay(combined_all_microglia) <- "RNA"
combined_all_microglia <- NormalizeData(combined_all_microglia)
combined_all_microglia <- FindVariableFeatures(combined_all_microglia)
combined_all_microglia <- ScaleData(combined_all_microglia)
View(combined_all_microglia)

combined_all_microglia$group <- paste(combined_all_microglia$condition,
                                      combined_all_microglia$sex, sep = "_")
combined_all_microglia <- JoinLayers(combined_all_microglia)

Idents(combined_all_microglia) <- "group"

View(combined_all_microglia)

# Differential expression analysis
deg_combined_female <- FindMarkers(
   combined_all_microglia,
   ident.1 = "COLD_Female",
   ident.2 = "CON_Female",
   assay = "RNA",
   logfc.threshold = 0.25,
   min.pct = 0.1
 )

# Filter upregulated genes (avg_log2FC > 0)
up_female <- deg_combined_female[deg_combined_female$avg_log2FC > 0,]
# Filter downregulated genes (avg_log2FC < 0)
down_female <- deg_combined_female[deg_combined_female$avg_log2FC < 0,]

cat("Female — Upregulated:", nrow(up_female), "\n")
cat("Female — Downregulated:", nrow(down_female), "\n")


# Differential expression analysis
deg_combined_male <- FindMarkers(
  combined_all_microglia,
  ident.1 = "COLD_Male",
  ident.2 = "CON_Male",
  assay = "RNA",
  logfc.threshold = 0.25,
  min.pct = 0.1
)

# DEGs
# Overepresentation
# how many degs you use impacts

# Filter upregulated genes (avg_log2FC > 0)
up_male <- deg_combined_male[deg_combined_male$avg_log2FC > 0,]

# Filter downregulated genes (avg_log2FC < 0)
down_male <- deg_combined_male[deg_combined_male$avg_log2FC < 0,]


cat("Male — Upregulated:", nrow(up_male), "\n")
cat("Male — Downregulated:", nrow(down_male), "\n")

#write.csv(up_male, "/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/Data Tables/upregulated-male.csv")
#write.csv(down_male, "/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/Data Tables/downregulated-male.csv")
#write.csv(up_female, "/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/Data Tables/upregulated-female.csv")
#write.csv(down_female, "/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/Data Tables/downregulated-female.csv")

# Save top 500 upregulated and downregulated DEGs for males and females (positive log2FC)
up_male_top250 <- up_male %>%
  arrange(desc(avg_log2FC)) %>%
  filter(p_val < 0.05) %>%
  head(250)

up_female_top250 <- up_female %>%
  arrange(desc(avg_log2FC)) %>%
  filter(p_val < 0.05) %>%
  head(250)


down_female_top250 <- down_female %>%
  arrange(avg_log2FC) %>%
  filter(p_val < 0.05) %>%
  head(250)

down_male_top250 <- down_male %>%
  arrange(avg_log2FC) %>%
  filter(p_val < 0.05) %>%
  head(250)

#write.csv(up_male_top500, "/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/Data Tables/upregulated-male-500.csv")
#write.csv(up_female_top500, "/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/Data Tables/upregulated-female-500.csv")
#write.csv(down_female_top500, "/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/Data Tables/downregulated-female-500.csv")
#write.csv(down_male_top500, "/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/Data Tables/downregulated-male-500.csv")

write.csv(up_male_top250, "/Users/alexlawson/MorphologyFinalAnalysis/rnasec/Data Tables/upregulated-male-250.csv")
write.csv(up_female_top250, "/Users/alexlawson/MorphologyFinalAnalysis/rnasec/Data Tables/upregulated-female-250.csv")
write.csv(down_female_top250, "/Users/alexlawson/MorphologyFinalAnalysis/rnasec/Data Tables/downregulated-female-250.csv")
write.csv(down_male_top250, "/Users/alexlawson/MorphologyFinalAnalysis/rnasec/Data Tables/downregulated-male-250.csv")


up_genes_male <- rownames(up_male_top250)
down_genes_male <- rownames(down_male_top250)
up_genes_female <- rownames(up_female_top250)
down_genes_female <- rownames(down_female_top250)

# Subset final merged object into female and male Seurat objects
combined_female <- subset(combined_all_microglia, subset = sex == "Female")
combined_male   <- subset(combined_all_microglia, subset = sex == "Male")

#extract RNA data from objects
rna_data_female <- GetAssayData(combined_female, assay = "RNA", layer = "data")
rna_data_male <- GetAssayData(combined_male, assay = "RNA", layer = "data")


#Female filtering for GO analysis
pct_expr_female <- rowSums(rna_data_female > 0) / ncol(rna_data_female) # Calculate percent of cells expressing each gene
gene_universe_female <- names(pct_expr_female[pct_expr_female > 0.1]) 
all_genes_to_convert_female <- unique(c(up_genes_female, down_genes_female, gene_universe_female)) # Combine all genes for conversion
gene_conversion_female <- bitr(all_genes_to_convert_female,
                        fromType = "SYMBOL",
                        toType = "ENSEMBL",
                        OrgDb = org.Mm.eg.db) # Convert to ENSEMBL IDs


#overrepresentation or gene set enrighment analysis 
#different in terms of calculations -> different things 
#Male filtering for GO analysis
pct_expr_male <- rowSums(rna_data_male > 0) / ncol(rna_data_male) # Calculate percent of cells expressing each gene
gene_universe_male <- names(pct_expr_male[pct_expr_male > 0.1]) # Filter genes expressed in >10% of cells (min.pct = 0.1)
all_genes_to_convert_male <- unique(c(up_genes_male, down_genes_male, gene_universe_male)) # Combine all genes for conversion
gene_conversion_male <- bitr(all_genes_to_convert_male,
                               fromType = "SYMBOL",
                               toType = "ENSEMBL",
                               OrgDb = org.Mm.eg.db) # Convert to ENSEMBL IDs

# Match converted gene lists for females
up_ensembl_female <- gene_conversion_female %>%
  filter(SYMBOL %in% up_genes_female) %>%
  pull(ENSEMBL)

down_ensembl_female <- gene_conversion_female %>%
  filter(SYMBOL %in% down_genes_female) %>%
  pull(ENSEMBL)

universe_ensembl_female <- gene_conversion_female %>%
  filter(SYMBOL %in% gene_universe_female) %>%
  pull(ENSEMBL)

# Match converted gene lists for males
up_ensembl_male <- gene_conversion_male %>%
  filter(SYMBOL %in% up_genes_male) %>%
  pull(ENSEMBL)

down_ensembl_male <- gene_conversion_male %>%
  filter(SYMBOL %in% down_genes_male) %>%
  pull(ENSEMBL)

universe_ensembl_male <- gene_conversion_male %>% #important bc then you are looking at the genes in your experiment 
  filter(SYMBOL %in% gene_universe_male) %>%
  pull(ENSEMBL)

View(ego_up_female)
ego_up_female <- enrichGO(gene          = up_ensembl_female,
                   universe      = universe_ensembl_female,
                   OrgDb         = org.Mm.eg.db,
                   keyType       = "ENSEMBL",
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05,
                   readable      = TRUE)


ego_down_female <- enrichGO(gene          = down_ensembl_female,
                     universe      = universe_ensembl_female,
                     OrgDb         = org.Mm.eg.db,
                     keyType       = "ENSEMBL",
                     ont           = "BP",
                     pAdjustMethod = "BH",
                     qvalueCutoff  = 0.05,   
                     readable      = TRUE)

ego_up_male <- enrichGO(gene          = up_ensembl_male,
                          universe      = universe_ensembl_male,
                          OrgDb         = org.Mm.eg.db,
                          keyType       = "ENSEMBL",
                          ont           = "BP",
                          pAdjustMethod = "BH",
                          qvalueCutoff  = 0.05,  
                          readable      = TRUE)

ego_down_male <- enrichGO(gene          = down_ensembl_male,
                            universe      = universe_ensembl_male,
                            OrgDb         = org.Mm.eg.db,
                            keyType       = "ENSEMBL",
                            ont           = "BP",
                            qvalueCutoff  = 0.05,   
                            pAdjustMethod = "BH",
                            readable      = TRUE)


write.csv(ego_down_female, "/Users/alexlawson/MorphologyFinalAnalysis/rnasec/Data Tables/GOdownregulated-female.csv")
write.csv(ego_up_female, "/Users/alexlawson/MorphologyFinalAnalysis/rnasec/Data Tables/GOupregulated-female.csv")
write.csv(ego_down_male, "/Users/alexlawson/MorphologyFinalAnalysis/rnasec/Data Tables/GOdownregulated-male.csv")
write.csv(ego_up_male, "/Users/alexlawson/MorphologyFinalAnalysis/rnasec/Data Tables/GOupregulated-male.csv")

# Barplot
#select genes and put cutoff. often good to run both separately and together. Run together first and then apart. universe = genes that you had some counts for 
#
# Barplot
barplot(ego_up_female, showCategory = 20, title = "Top GO BP Terms - Upregulated Female")
barplot(ego_down_female, showCategory = 20, title = "Top GO BP Terms - Downregulated Female")
barplot(ego_up_male, showCategory = 20, title = "Top GO BP Terms - Upregulated Male")
barplot(ego_down_male, showCategory = 20, title = "Top GO BP Terms - Downregulated Male")


ego_up_female_MF <- enrichGO(gene          = up_ensembl_female,
                          universe      = universe_ensembl_female,
                          OrgDb         = org.Mm.eg.db,
                          keyType       = "ENSEMBL",
                          ont           = "MF",
                          pAdjustMethod = "BH",
                          qvalueCutoff  = 0.05, 
                          readable      = TRUE)


ego_down_female_MF <- enrichGO(gene          = down_ensembl_female,
                            universe      = universe_ensembl_female,
                            OrgDb         = org.Mm.eg.db,
                            keyType       = "ENSEMBL",
                            ont           = "MF",
                            pAdjustMethod = "BH",
                            qvalueCutoff  = 0.05,  
                            readable      = TRUE)

ego_up_male_MF <- enrichGO(gene          = up_ensembl_male,
                        universe      = universe_ensembl_male,
                        OrgDb         = org.Mm.eg.db,
                        keyType       = "ENSEMBL",
                        ont           = "MF",
                        pAdjustMethod = "BH",
                        qvalueCutoff  = 0.05,  
                        readable      = TRUE)

#terms that have a given number of genes 
ego_down_male_MF <- enrichGO(gene          = down_ensembl_male,
                          universe      = universe_ensembl_male,
                          OrgDb         = org.Mm.eg.db,
                          keyType       = "ENSEMBL",
                          ont           = "MF",
                          qvalueCutoff  = 0.05,   
                          pAdjustMethod = "BH",
                          readable      = TRUE)

barplot(ego_up_female_MF, showCategory = 20, title = "Top GO MF Terms - Upregulated Female")
barplot(ego_down_female_MF, showCategory = 20, title = "Top GO MF Terms - Downregulated Female")
barplot(ego_up_male_MF, showCategory = 20, title = "Top GO MF Terms - Upregulated Male")
barplot(ego_down_male_MF, showCategory = 20, title = "Top GO MF Terms - Downregulated Male")


# Set the identity class to renamed_clusters
Idents(combined_all_microglia) <- "renamed_clusters"

# Subset only cells from cluster 2
cluster2_microglia <- subset(combined_all_microglia, idents = "Microglia 2")
Idents(cluster2_microglia) <- "group"

# Differential expression analysis
deg_combined_female_2 <- FindMarkers(
  cluster2_microglia,
  ident.1 = "COLD_Female",
  ident.2 = "CON_Female",
  assay = "RNA",
  logfc.threshold = 0.25,
  min.pct = 0.2
)

# Filter upregulated genes (avg_log2FC > 0)
up_female_2 <- deg_combined_female_2[deg_combined_female_2$avg_log2FC > 0.2,]
# Filter downregulated genes (avg_log2FC < 0)
down_female_2 <- deg_combined_female_2[deg_combined_female_2$avg_log2FC < 0.2,]

cat("Female — Upregulated:", nrow(up_female_2), "\n")
cat("Female — Downregulated:", nrow(down_female_2), "\n")


# Differential expression analysis
deg_combined_male_2 <- FindMarkers(
  cluster2_microglia,
  ident.1 = "COLD_Male",
  ident.2 = "CON_Male",
  assay = "RNA",
  logfc.threshold = 0.25,
  min.pct = 0.2
)

# DEGs
# Overepresentation
# how many degs you use impacts

# Filter upregulated genes (avg_log2FC > 0)
up_male_2 <- deg_combined_male_2[deg_combined_male_2$avg_log2FC > 0.2,]

# Filter downregulated genes (avg_log2FC < 0)
down_male_2 <- deg_combined_male_2[deg_combined_male_2$avg_log2FC < 0.2,]


cat("Male — Upregulated:", nrow(up_male_2), "\n")
cat("Male — Downregulated:", nrow(down_male_2), "\n")

# Save top 500 up-/downregulated DEGs for males and females in cluster 2
up_male_top500_2 <- up_male_2 %>%
  arrange(desc(avg_log2FC)) %>%
  head(500)

up_female_top500_2 <- up_female_2 %>%
  arrange(desc(avg_log2FC)) %>%
  head(500)

down_female_top500_2 <- down_female_2 %>%
  arrange(avg_log2FC) %>%
  head(500)

down_male_top500_2 <- down_male_2 %>%
  arrange(avg_log2FC) %>%
  head(500)

# (optional) write out
# write.csv(up_male_top500_2, "/…/upregulated-male-500.csv")
# write.csv(up_female_top500_2, "/…/upregulated-female-500.csv")
# write.csv(down_female_top500_2, "/…/downregulated-female-500.csv")
# write.csv(down_male_top500_2, "/…/downregulated-male-500.csv")

# Gene name vectors
up_genes_male_2   <- rownames(up_male_top500_2)
down_genes_male_2 <- rownames(down_male_top500_2)
up_genes_female_2   <- rownames(up_female_top500_2)
down_genes_female_2 <- rownames(down_female_top500_2)

# Subset the cluster-2 Seurat object by sex
combined_female_2 <- subset(combined_all_microglia, subset = sex == "Female")
combined_male_2   <- subset(combined_all_microglia, subset = sex == "Male")

# Extract RNA assay data
rna_data_female_2 <- GetAssayData(combined_female_2, assay = "RNA", layer = "data")
rna_data_male_2   <- GetAssayData(combined_male_2,   assay = "RNA", layer = "data")

# Female filtering for GO
pct_expr_female_2        <- rowSums(rna_data_female_2 > 0) / ncol(rna_data_female_2)
gene_universe_female_2   <- names(pct_expr_female_2[pct_expr_female_2 > 0.1])
all_genes_to_convert_female_2 <- unique(c(up_genes_female_2, down_genes_female_2, gene_universe_female_2))
gene_conversion_female_2 <- bitr(all_genes_to_convert_female_2,
                                 fromType = "SYMBOL",
                                 toType   = "ENSEMBL",
                                 OrgDb    = org.Mm.eg.db)

# Male filtering for GO
pct_expr_male_2        <- rowSums(rna_data_male_2 > 0) / ncol(rna_data_male_2)
gene_universe_male_2   <- names(pct_expr_male_2[pct_expr_male_2 > 0.1])
all_genes_to_convert_male_2 <- unique(c(up_genes_male_2, down_genes_male_2, gene_universe_male_2))
gene_conversion_male_2 <- bitr(all_genes_to_convert_male_2,
                               fromType = "SYMBOL",
                               toType   = "ENSEMBL",
                               OrgDb    = org.Mm.eg.db)

# Match to ENSEMBL IDs
up_ensembl_female_2   <- gene_conversion_female_2 %>% filter(SYMBOL %in% up_genes_female_2)   %>% pull(ENSEMBL)
down_ensembl_female_2 <- gene_conversion_female_2 %>% filter(SYMBOL %in% down_genes_female_2) %>% pull(ENSEMBL)
universe_ensembl_female_2 <- gene_conversion_female_2 %>% filter(SYMBOL %in% gene_universe_female_2) %>% pull(ENSEMBL)

up_ensembl_male_2     <- gene_conversion_male_2 %>% filter(SYMBOL %in% up_genes_male_2)     %>% pull(ENSEMBL)
down_ensembl_male_2   <- gene_conversion_male_2 %>% filter(SYMBOL %in% down_genes_male_2)   %>% pull(ENSEMBL)
universe_ensembl_male_2   <- gene_conversion_male_2 %>% filter(SYMBOL %in% gene_universe_male_2)   %>% pull(ENSEMBL)

# GO enrichment
ego_up_female_2 <- enrichGO(gene          = up_ensembl_female_2,
                            universe      = universe_ensembl_female_2,
                            OrgDb         = org.Mm.eg.db,
                            keyType       = "ENSEMBL",
                            ont           = "BP",
                            pAdjustMethod = "BH",
                            qvalueCutoff  = 0.05,
                            readable      = TRUE)

ego_down_female_2 <- enrichGO(gene          = down_ensembl_female_2,
                              universe      = universe_ensembl_female_2,
                              OrgDb         = org.Mm.eg.db,
                              keyType       = "ENSEMBL",
                              ont           = "BP",
                              pAdjustMethod = "BH",
                              qvalueCutoff  = 0.05,
                              readable      = TRUE)

ego_up_male_2 <- enrichGO(gene          = up_ensembl_male_2,
                          universe      = universe_ensembl_male_2,
                          OrgDb         = org.Mm.eg.db,
                          keyType       = "ENSEMBL",
                          ont           = "BP",
                          pAdjustMethod = "BH",
                          qvalueCutoff  = 0.05,
                          readable      = TRUE)

ego_down_male_2 <- enrichGO(gene          = down_ensembl_male_2,
                            universe      = universe_ensembl_male_2,
                            OrgDb         = org.Mm.eg.db,
                            keyType       = "ENSEMBL",
                            ont           = "BP",
                            pAdjustMethod = "BH",
                            qvalueCutoff  = 0.05,
                            readable      = TRUE)

# Barplots
barplot(ego_up_female_2,   showCategory = 20, title = "Top GO BP Terms - Upregulated Female (Cluster 2)")
barplot(ego_down_female_2, showCategory = 20, title = "Top GO BP Terms - Downregulated Female (Cluster 2)")
barplot(ego_up_male_2,     showCategory = 20, title = "Top GO BP Terms - Upregulated Male (Cluster 2)")
barplot(ego_down_male_2,   showCategory = 20, title = "Top GO BP Terms - Downregulated Male (Cluster 2)")









