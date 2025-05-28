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

set.seed(90)

#setting directories
cold_data_dir <- "/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/COLD"
con_data_dir <- "/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/CON"

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


#### Step 5: switch to integrated assay, then… ####
DefaultAssay(combined) <- "integrated"

# ───── 5a) Cell-cycle scoring ─────
# Seurat ships with two vectors of S- and G2/M-phase markers called `cc.genes`
# (or use the updated 2019 list: cc.genes.updated.2019)
s.features  <- cc.genes.updated.2019$s.genes
g2m.features <- cc.genes.updated.2019$g2m.genes

combined <- CellCycleScoring(
  object      = combined,
  s.features  = s.features,
  g2m.features= g2m.features,
  set.ident   = FALSE        # don’t overwrite your existing identities
)

# now each cell has metadata columns: S.Score, G2M.Score, and Phase

# ───── 5b) ScaleData while regressing out cell-cycle ─────
all_genes_combined <- rownames(combined)
combined <- ScaleData(
  object       = combined,
  features     = all_genes_combined,
  vars.to.regress = c("S.Score", "G2M.Score")
)

# ───── 5c) then continue with PCA/clustering/UMAP as before ─────
combined <- RunPCA(combined, features = VariableFeatures(combined))
features = VariableFeatures(combined)
#combined <- ScaleData(combined, features = all_genes_combined)
#combined <- RunPCA(combined, features = VariableFeatures(object = combined))
ElbowPlot(combined, ndims=40, reduction = "pca")
combined <- FindNeighbors(combined, dims=1:20)
combined <- FindClusters(combined, resolution=0.5)
combined <- RunUMAP(combined, dims=1:20)
umap_combined <- DimPlot(combined, reduction="umap", group.by="seurat_clusters")
umap_combined


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

plot_macrophage_vln

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

plot_cluster1_vln

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

plot_cluster2_vln

#ggsave(filename = '/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/PLOTS/clust02_vln.png', plot_cluster2_vln, width = 10, height = 5)
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

plot_cluster4_vln

#ggsave(filename = '/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/PLOTS/clust04_vln.png', plot_cluster4_vln, width = 10, height = 5)

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


new.cluster.ids <- c("Microglia 4", "Microglia 4", "Microglia 4", "Microglia 2", "Proliferative", "Microglia 3", "Proliferative", "Microglia 1", "Macrophage")
names(new.cluster.ids) <- levels(combined)
combined <- RenameIdents(combined, new.cluster.ids)
DimPlot(combined, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(combined, reduction="umap", group.by="condition")
combined$renamed_clusters <- Idents(combined)
table(Idents(combined)) 
Idents(combined) <- "renamed_clusters"

keep <- paste0("Microglia ", 1:4)
keep_wmacro <- c(paste0("Microglia ", 1:4), "Macrophage")
combined_clean <- subset(combined, idents = keep)
combined_clean_wmacro <- subset(combined, idents = keep_wmacro)
updated_umap <- DimPlot(combined_clean, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
#ggsave(filename = '/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/PLOTS/umap-with-clusters.png', updated_umap, width = 10, height = 5)


control <- subset(combined_clean, subset = condition == "CON")
cold <- subset(combined_clean, subset = condition == "COLD")

# Define sex marker gene for female classification
female_genes <- c("Xist")
expr_threshold <- 0  # Any non-zero expression will count

# Re-pull your expression matrices (so you haven’t overwritten them)
DefaultAssay(control) <- "RNA"
expr_control <- GetAssayData(control, assay="RNA", slot="data")

DefaultAssay(cold) <- "RNA"
expr_cold <- GetAssayData(cold,   assay="RNA", slot="data")

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

View(combined_all_microglia)
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


combined_female <- subset(
  combined_all_microglia,
  subset = sex == "Female"
)

combined_male <- subset(
  combined_all_microglia,
  subset = sex == "Male"
)

# Differential expression analysis
deg_combined_female <- FindMarkers(
  combined_all_microglia,
  ident.1 = "COLD_Female",
  ident.2 = "CON_Female",
  assay = "RNA",
  logfc.threshold = 0.25,
  min.pct = 0.1
)

# Filter female DEGs
up_female <- deg_combined_female[deg_combined_female$avg_log2FC > 0 & deg_combined_female$p_val_adj < 0.1,]
down_female <- deg_combined_female[deg_combined_female$avg_log2FC < 0 & deg_combined_female$p_val_adj < 0.1,]
View(up_female)
View(down_female)
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
# Filter upregulated genes (avg_log2FC > 0)
upregulated_genes_male_convcold <- deg_combined_male[deg_combined_male$avg_log2FC > 0 & deg_combined_male$p_val_adj < 0.1,]

# Filter downregulated genes (avg_log2FC < 0)
downregulated_genes_male_convcold <- deg_combined_male[deg_combined_male$avg_log2FC < 0 & deg_combined_male$p_val_adj < 0.1,]

cat("Male — Upregulated:", nrow(upregulated_genes_male_convcold), "\n")
cat("Male — Downregulated:", nrow(downregulated_genes_male_convcold), "\n")

View(combined_all_microglia)

# Convert gene symbols to Entrez IDs
up_male_entrez <- bitr(rownames(up_male), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID
up_female_entrez <- bitr(rownames(up_female), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID
down_male_entrez <- bitr(rownames(down_male), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID
down_female_entrez <- bitr(rownames(down_female), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)$ENTREZID

go_bp <- godata('org.Mm.eg.db', ont = "BP")  # You can change to "MF" or "CC" if needed

# Upregulated comparison
sim_up <- mgeneSim(up_male_entrez, up_female_entrez, semData = go_bp, measure = "Wang", combine = "BMA")

# Downregulated comparison
sim_down <- mgeneSim(down_male_entrez, down_female_entrez, semData = go_bp, measure = "Wang", combine = "BMA")

# Summary scores
mean(sim_up, na.rm = TRUE)
mean(sim_down, na.rm = TRUE)

pheatmap(sim_up, cluster_rows = TRUE, cluster_cols = TRUE, main = "Semantic Similarity: Upregulated (Male vs Female)")
pheatmap(sim_down, cluster_rows = TRUE, cluster_cols = TRUE, main = "Semantic Similarity: Downregulated (Male vs Female)")

ego_up_male <- enrichGO(gene = up_male_entrez, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05, readable = TRUE)
ego_down_male <- enrichGO(gene = down_male_entrez,   OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05, readable = TRUE)
ego_up_female <- enrichGO(gene = up_female_entrez,   OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05, readable = TRUE)
ego_down_female <- enrichGO(gene = down_female_entrez, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH", pvalueCutoff  = 0.05, qvalueCutoff  = 0.05, readable = TRUE)

dotplot(ego_up_male, showCategory = 20) + 
  ggtitle("GO-BP: Upregulated Genes in Male (COLD vs CON)")

dotplot(ego_up_female, showCategory = 20) + 
  ggtitle("GO-BP: Upregulated Genes in Female (COLD vs CON)")

dotplot(ego_down_male, showCategory = 20) + 
  ggtitle("GO-BP: Downregulated Genes in Male (COLD vs CON)")

dotplot(ego_down_female, showCategory = 20) + 
  ggtitle("GO-BP: Downregulated Genes in Female (COLD vs CON)")


res_up_male <- as.data.frame(ego_up_male)
res_down_male <- as.data.frame(ego_down_male)

gene_lists <- list(
  Male_Up   = up_male_entrez,
  Female_Up = up_female_entrez,
  Male_Down = down_male_entrez,
  Female_Down = down_female_entrez
)

cc <- compareCluster(
  geneCluster   = gene_lists,
  fun           = "enrichGO",
  OrgDb         = org.Mm.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

dotplot(cc, showCategory = 10) +
  facet_wrap(~ Cluster, ncol = 4) +
  ggtitle("GO-BP Enrichment: Male vs Female (Up/Down)")


# Differential expression analysis
deg_combined_control <- FindMarkers(
  combined_all_microglia,
  ident.1 = "COLD_Male",
  ident.2 = "CON_Male",
  assay = "RNA",
  logfc.threshold = 0.25,
  min.pct = 0.1
)

# DEGs
# Filter upregulated genes (avg_log2FC > 0)
upregulated_genes_control <- deg_combined_control[deg_combined_control$avg_log2FC > 0 & deg_combined_control$p_val_adj < 0.1,]

# Filter downregulated genes (avg_log2FC < 0)
downregulated_genes_control <- deg_combined_control[deg_combined_control$avg_log2FC < 0 & deg_combined_control$p_val_adj < 0.1,]

cat("Female — Upregulated:", nrow(upregulated_genes_control), "\n")
cat("Female — Downregulated:", nrow(downregulated_genes_control), "\n")



# cluster_2 <- subset(
#   combined_all_microglia,
#   subset = renamed_clusters == "Microglia 2"
# )
# 
# combined_female_cluster2 <- subset(
#   cluster_2,
#   subset = sex == "Female"
# )
# 
# combined_male_cluster2 <- subset(
#   cluster_2,
#   subset = sex == "Male"
# )
# 
# Idents(combined_female_cluster2) <- "condition"
# table(Idents(combined_female_cluster2))
# # JoinLayers to properly link for analysis 
# combined_female_cluster2 <- JoinLayers(combined_female_cluster2)
# # Differential expression analysis
# deg_combined_female_cluster2 <- FindMarkers(
#   combined_female_cluster2,
#   ident.1 = "COLD",
#   ident.2 = "CON",
#   assay = "RNA",
#   test.use = "wilcox",
#   logfc.threshold = 0.25,
#   min.pct = 0.1
# )
# # View and save results
# View(deg_combined_female_cluster2)
# write.csv(deg_combined_female_cluster2, "/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/Data Tables/female_cold_vs_con-cluster2.csv")
# 
# # Set the identities for DE analysis
# Idents(combined_male_cluster2) <- "condition"
# # JoinLayers to properly link for analysis 
# combined_male_cluster2 <- JoinLayers(combined_male_cluster2)
# # Differential expression analysis
# deg_combined_male_cluster2 <- FindMarkers(
#   combined_male_cluster2,
#   ident.1 = "COLD",
#   ident.2 = "CON",
#   assay = "RNA",
#   test.use = "wilcox",
#   logfc.threshold = 0.25,
#   min.pct = 0.1
# )
# # View and save results
# head(deg_combined_male_cluster2)
# write.csv(deg_combined_male_cluster2, "/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/Data Tables/male_cold_vs_con-cluster2.csv")
# 
# # DEGs
# # Filter upregulated genes (avg_log2FC > 0)
# upregulated_genes_male_convcold_cluster2 <- deg_combined_male_cluster2[deg_combined_male_cluster2$avg_log2FC > 0 & deg_combined_male_cluster2$p_val_adj < 0.1, ]
# 
# # Filter downregulated genes (avg_log2FC < 0)
# downregulated_genes_male_convcold_cluster2 <- deg_combined_male_cluster2[deg_combined_male_cluster2$avg_log2FC < 0 & deg_combined_male_cluster2$p_val_adj < 0.1, ]
# 
# # Save male DEGs
# write.csv(deg_combined_male_cluster2, "male_cold_vs_con-cluster2.csv")
# 
# # Filter male DEGs
# up_male_cluster2 <- deg_combined_male_cluster2[deg_combined_male_cluster2$avg_log2FC > 0 & deg_combined_male_cluster2$p_val_adj < 0.1, ]
# down_male_cluster2 <- deg_combined_male_cluster2[deg_combined_male_cluster2$avg_log2FC < 0 & deg_combined_male_cluster2$p_val_adj < 0.1, ]
# 
# cat("Male — Upregulated:", nrow(up_male_cluster2), "\n")
# cat("Male — Downregulated:", nrow(down_male_cluster2), "\n")
# 
# # Save female DEGs
# write.csv(deg_combined_female, "female_cold_vs_con.csv")
# 
# # Filter female DEGs
# up_female <- deg_combined_female[deg_combined_female$avg_log2FC > 0 & deg_combined_female$p_val_adj < 0.1, ]
# down_female <- deg_combined_female[deg_combined_female$avg_log2FC < 0 & deg_combined_female$p_val_adj < 0.1, ]
# View(up_female)
# View(down_female)
# cat("Female — Upregulated:", nrow(up_female), "\n")
# cat("Female — Downregulated:", nrow(down_female), "\n")
# 

