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

set.seed(99)

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
View(combined)
# 5a) Cell‐cycle scoring
s.features   <- cc.genes.updated.2019$s.genes
g2m.features <- cc.genes.updated.2019$g2m.genes

combined <- CellCycleScoring(
  combined,
  s.features   = s.features,
  g2m.features = g2m.features,
  set.ident    = FALSE
)

# 5b) Scale while regressing out cell‐cycle scores
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
plot_macrophage_feature
plot_macrophage_vln

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

plot_cluster2_vln
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

# rest of workflow----------------------------------------------
new.cluster.ids <- c("Microglia 4", "Microglia 2", "Microglia 4", "Microglia 3", "Proliferative", "Microglia 4", "Proliferative", "Microglia 4", "Microglia 1", "Pericyte", "Macrophage")
names(new.cluster.ids) <- levels(combined)
combined <- RenameIdents(combined, new.cluster.ids)
DimPlot(combined, reduction = "umap", label = TRUE, pt.size = 0.5)
DimPlot(combined, reduction="umap", group.by="condition")
combined$renamed_clusters <- Idents(combined)
table(Idents(combined)) 
Idents(combined) <- "renamed_clusters"

keep <- paste0("Microglia ", 1:4)
combined_clean <- subset(combined, idents = keep)
DimPlot(combined_clean, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

combined_clean <- JoinLayers(combined_clean)
View(combined_clean)
# RNA assay -> active and pull raw normalized data
DefaultAssay(combined_clean) <- "RNA"
expr <- GetAssayData(combined_clean, assay = "RNA", slot = "data")

#Annotate sex based on Xist expression
female_cells <- colnames(expr)[expr["Xist", ] > 0]
combined_clean$sex <- ifelse(colnames(combined_clean) %in% female_cells, "Female", "Male")

#Annotate condition (if you already have combined_clean$condition, skip this)
if (!"condition" %in% colnames(combined_clean@meta.data)) {
  combined_clean$condition <- ifelse(grepl("COLD$", colnames(combined_clean)), "COLD", "CON")
}

#Define a combined factor for easy subsetting and DE
combined_clean$group <- paste(combined_clean$condition, combined_clean$sex, sep = "_")
Idents(combined_clean) <- "group"

#(Re)normalize, find variable features, and scale—only if you need to rerun downstream analyses
combined_clean <- NormalizeData(combined_clean)
combined_clean <- FindVariableFeatures(combined_clean)
combined_clean <- ScaleData(combined_clean)

#Female COLD vs CON
deg_female <- FindMarkers(
  combined_clean,
  ident.1 = "COLD_Female",
  ident.2 = "CON_Female",
  assay = "RNA",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
up_female   <- subset(deg_female, avg_log2FC > 0.2)
down_female <- subset(deg_female, avg_log2FC < -0.2)
cat("Female — Upregulated:",   nrow(up_female),   "\n")
cat("Female — Downregulated:", nrow(down_female), "\n")

#Male COLD vs CON
deg_male <- FindMarkers(
  combined_clean,
  ident.1 = "COLD_Male",
  ident.2 = "CON_Male",
  assay = "RNA",
  logfc.threshold = 0.25,
  min.pct = 0.1
)
up_male   <- subset(deg_male, avg_log2FC > 0.2)
down_male <- subset(deg_male, avg_log2FC < -0.2)
cat("Male — Upregulated:",   nrow(up_male),   "\n")
cat("Male — Downregulated:", nrow(down_male), "\n")




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

write.csv(ego_up_male, "/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/Data Tables/GO-up-males.csv")
write.csv(ego_up_female, "/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/Data Tables/GO-up-females.csv")
write.csv(ego_down_male, "/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/Data Tables/GO-down-males.csv")
write.csv(ego_down_female, "/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/Data Tables/GO-down-females.csv")


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



# #Calculate percent mitochondrial reads
# combined_all_microglia_2 <- combined_all_microglia
# combined_all_microglia_2[["percent.mt"]] <- PercentageFeatureSet(
#   object = combined_all_microglia_2,
#   pattern = "^mt-"
# )
# 
# View(combined_all_microglia_2)
# # Violin plots for quick QC
# VlnPlot(combined_all_microglia_2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), group.by = "renamed_clusters")
# 
# # Extract metadata
# md <- combined_all_microglia_2@meta.data
# 
# # Pairwise t-tests with Bonferroni correction
# 
# # 1. Gene count (number of detected genes per cell)
# gene_pairwise <- pairwise.t.test(
#   x = md$nFeature_RNA,
#   g = md$renamed_clusters,
#   p.adjust.method = "bonferroni"
# )
# 
# # 2. UMI count (total transcript count per cell)
# umi_pairwise <- pairwise.t.test(
#   x = md$nCount_RNA,
#   g = md$renamed_clusters,
#   p.adjust.method = "bonferroni"
# )
# 
# # 3. Mitochondrial content (already done, but keeping for reference)
# mt_pairwise <- pairwise.t.test(
#   x = md$percent.mt,
#   g = md$renamed_clusters,
#   p.adjust.method = "bonferroni"
# )
# 
# # View p-values for each
# gene_pairwise$p.value      # gene count differences
# umi_pairwise$p.value       # UMI count differences
# mt_pairwise$p.value        # mitochondrial differences
# 

