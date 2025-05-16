library(dplyr)
library(Seurat)
library(patchwork)
library(clustree)
library(sctransform)

# use markers from paper to cluster the microglia 
# go over representation for pathway 
# everything with xist is female and everything else is male 

cold_data_dir <- "/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/COLD"
con_data_dir  <- "/Users/alexlawson/Documents/GitHub/MorphologyFinalAnalysis/rnasec/CON"

# Load Control
control_counts <- Read10X(data.dir = con_data_dir)
# Load Cold Stress
cold_counts    <- Read10X(data.dir = cold_data_dir)

control <- CreateSeuratObject(counts = control_counts, project = "Control", min.cells = 3, min.features = 200)
cold    <- CreateSeuratObject(counts = cold_counts,    project = "COLD",   min.cells = 3, min.features = 200)

control[["percent.mt"]] <- PercentageFeatureSet(control, pattern = "^mt-")
cold[["percent.mt"]]    <- PercentageFeatureSet(cold,    pattern = "^mt-")

control <- subset(control, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)
cold    <- subset(cold,    subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

control$condition <- "CON"
cold$condition    <- "COLD"

# 2. Put them in a list and preprocess each
objs <- list(control, cold)

for (i in seq_along(objs)) {
  objs[[i]] <- SCTransform(objs[[i]], verbose = FALSE)
}

anchors <- FindIntegrationAnchors(
  object.list   = objs,
  reduction     = "rpca",
  dims          = 1:20
)
combined <- IntegrateData(anchorset = anchors, dims = 1:20)


# 5. Switch to the new “integrated” assay
DefaultAssay(combined) <- "integrated"

# ── Cell‐cycle scoring and regression ────────────────────────────────────────

# rename gene lists so you don't overwrite anything
cc_s_genes   <- cc.genes$s.genes
cc_g2m_genes <- cc.genes$g2m.genes

# score cells for S vs. G2/M phase
combined_cc <- CellCycleScoring(
  object       = combined,
  s.features   = cc_s_genes,
  g2m.features = cc_g2m_genes,
  set.ident    = FALSE
)

# regress out cell‐cycle variation during scaling
all_genes_cc <- rownames(combined_cc)
combined_cc  <- ScaleData(
  combined_cc,
  features       = all_genes_cc,
  vars.to.regress = c("S.Score", "G2M.Score"),
  verbose        = FALSE
)

# ── Continue with PCA, clustering, UMAP ─────────────────────────────────────
combined_cc <- RunPCA(combined_cc, features = VariableFeatures(combined_cc))
ElbowPlot(combined_cc, ndims = 40, reduction = "pca")

combined_cc <- FindNeighbors(combined_cc, dims = 1:20)
combined_cc <- FindClusters(combined_cc, resolution = 0.5)
combined_cc <- RunUMAP(combined_cc, dims = 1:20)

# plot
umap_cc <- DimPlot(combined_cc, reduction = "umap", group.by = "seurat_clusters")
umap_cc
