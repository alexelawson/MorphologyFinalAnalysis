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
cold_data_dir <- "/Users/alexlawson/GitHub/MorphologyFinalAnalysis/rnasec/COLD"
con_data_dir <- "/Users/alexlawson/GitHub/MorphologyFinalAnalysis/rnasec/CON"

# Load Control
control_counts <- Read10X(data.dir = con_data_dir)
#Load Cold Stress
cold_counts <- Read10X(data.dir = cold_data_dir)

#creating seurat objects
control <- CreateSeuratObject(counts = control_counts, project = "Control", min.cells = 3, min.features = 200)
cold <- CreateSeuratObject(counts = cold_counts, project = "Cold", min.cells = 3, min.features = 200)

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
up_female <- deg_combined_female[deg_combined_female$avg_log2FC > 0 & deg_combined_female$p_val_adj < 0.05,]

# Filter downregulated genes (avg_log2FC < 0)
down_female <- deg_combined_female[deg_combined_female$avg_log2FC < 0 & deg_combined_female$p_val_adj < 0.05,]
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
up_male <- deg_combined_male[deg_combined_male$avg_log2FC > 0 & deg_combined_male$p_val_adj < 0.05,]

# Filter downregulated genes (avg_log2FC < 0)
down_male <- deg_combined_male[deg_combined_male$avg_log2FC < 0 & deg_combined_male$p_val_adj < 0.05,]


cat("Male — Upregulated:", nrow(up_male), "\n")
cat("Male — Downregulated:", nrow(down_male), "\n")

#write.csv(up_male, "/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/Data Tables/upregulated-male.csv")
#write.csv(down_male, "/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/Data Tables/downregulated-male.csv")
#write.csv(up_female, "/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/Data Tables/upregulated-female.csv")
#write.csv(down_female, "/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/Data Tables/downregulated-female.csv")

# Save top upregulated and downregulated DEGs for males and females (positive log2FC)
up_male_all <- up_male %>%
  arrange(desc(avg_log2FC))

up_female_all <- up_female %>%
  arrange(desc(avg_log2FC))

down_female_all <- down_female %>%
  arrange(avg_log2FC)

down_male_all<- down_male %>%
  arrange(avg_log2FC)

View(down_male_all)

#write.csv(up_male_all, "/Users/alexlawson/GitHub/MorphologyFinalAnalysis/rnasec/Data Tables/upregulated-male-pvaladjust.csv")
#write.csv(up_female_all, "/Users/alexlawson/GitHub/MorphologyFinalAnalysis/rnasec/Data Tables/upregulated-female--pvaladjust.csv")
#write.csv(down_female_all, "/Users/alexlawson/GitHub/MorphologyFinalAnalysis/rnasec/Data Tables/downregulated-female--pvaladjust.csv")
#write.csv(down_male_all, "/Users/alexlawson/GitHub/MorphologyFinalAnalysis/rnasec/Data Tables/downregulated-male--pvaladjust.csv")

up_genes_male <- rownames(up_male_all)
down_genes_male <- rownames(down_male_all)
up_genes_female <- rownames(up_female_all)
down_genes_female <- rownames(down_female_all)

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



#write.csv(ego_down_female, "/Users/alexlawson/GitHub/MorphologyFinalAnalysis/rnasec/Data Tables/GOdownregulated-female.csv")
#write.csv(ego_up_female, "/Users/alexlawson/GitHub/MorphologyFinalAnalysis/rnasec/Data Tables/GOupregulated-female.csv")
#write.csv(ego_down_male, "/Users/alexlawson/GitHub/MorphologyFinalAnalysis/rnasec/Data Tables/GOdownregulated-male.csv")
#write.csv(ego_up_male, "/Users/alexlawson/GitHub/MorphologyFinalAnalysis/rnasec/Data Tables/GOupregulated-male.csv")

# Barplot
#select genes and put cutoff. often good to run both separately and together. Run together first and then apart. universe = genes that you had some counts for 


#enrichment analysis 

library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(tibble)

pct_expr_female <- rowSums(rna_data_female > 0) / ncol(rna_data_female)
gene_universe_female <- names(pct_expr_female)[ pct_expr_female > 0.1 ] 
#changing the row column to actually have a header (symbol)
deg_combined_female_sym <- deg_combined_female %>%
  rownames_to_column(var = "SYMBOL")
female_DE_filtered <- deg_combined_female_sym %>%
  filter(SYMBOL %in% gene_universe_female) 

gene_conversion_female <- bitr(
  female_DE_filtered$SYMBOL,
  fromType = "SYMBOL",
  toType   = "ENSEMBL",
  OrgDb    = org.Mm.eg.db
)

female_DE_annotated <- female_DE_filtered %>%
  inner_join(gene_conversion_female, by = "SYMBOL") %>%
  select(ENSEMBL, avg_log2FC)

geneList_female <- female_DE_annotated$avg_log2FC
names(geneList_female) <- female_DE_annotated$ENSEMBL
geneList_female <- geneList_female[ !is.na(geneList_female) ]
geneList_female <- sort(geneList_female, decreasing = TRUE)

gse_results_female <- gseGO(
  geneList     = geneList_female,
  OrgDb        = org.Mm.eg.db,
  ont          = "BP",
  keyType      = "ENSEMBL",
  pvalueCutoff = 0.05,
  verbose      = TRUE
)

#males
pct_expr_male <- rowSums(rna_data_male > 0) / ncol(rna_data_male)
gene_universe_male <- names(pct_expr_male)[ pct_expr_male > 0.1 ] 
deg_combined_male_sym <- deg_combined_male %>%
  rownames_to_column(var = "SYMBOL")
male_DE_filtered <- deg_combined_male_sym %>%
  filter(SYMBOL %in% gene_universe_male) 

gene_conversion_male <- bitr(
  male_DE_filtered$SYMBOL,
  fromType = "SYMBOL",
  toType   = "ENSEMBL",
  OrgDb    = org.Mm.eg.db
)
male_DE_annotated <- male_DE_filtered %>%
  inner_join(gene_conversion_male, by = "SYMBOL") %>%
  select(ENSEMBL, avg_log2FC)

geneList_male <- male_DE_annotated$avg_log2FC
names(geneList_male) <- male_DE_annotated$ENSEMBL
geneList_male <- geneList_male[ !is.na(geneList_male) ]
geneList_male <- sort(geneList_male, decreasing = TRUE)

gse_results_male <- gseGO(
  geneList     = geneList_male,
  OrgDb        = org.Mm.eg.db,
  ont          = "BP",
  keyType      = "ENSEMBL",
  pvalueCutoff = 0.05,
  verbose      = TRUE
)


# Ridgeplot of top 10 – shows distribution of enrichment scores
ridgeplot(
  gse_results_male,
  showCategory = 20,
)



goplot(ego_down_female, GOANCESTOR = GOBPANCESTOR)



