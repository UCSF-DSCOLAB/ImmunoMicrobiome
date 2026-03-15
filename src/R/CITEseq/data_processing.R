library(tidyverse)
library(Seurat)
library(dsb)
library(harmony)

proc_dir = "data/single_cell_GEX/processed"
out_dir = "data/scRNA_CITE/DB1-7_tmp/"
samples = list.files(proc_dir, pattern = "XHLT2") %>% grep(x=., "tmp", value = T, invert = T)


sobjs = list()
stat = list()
for( s in samples ) {
  print(paste0("Processing data for ", s))
  stat_local = list()

  sobj_local = readRDS( file.path(proc_dir, s, paste0("/qc_and_addModalities/",s,"_seurat_object_qcAndAddModalities.rds")) )

  # Perform dsb normalization for ADT data.
  raw_data <- Read10X_h5( file.path(proc_dir, s, "/cellranger/raw_feature_bc_matrix.h5") )
  empty_drops <- colnames(raw_data$`Gene Expression`)[colSums(raw_data$`Gene Expression`) < 100 & colSums(raw_data$`Antibody Capture`) > 10 ]
  empty_drops <- empty_drops[! empty_drops %in% Cells(sobj_local)]
  ADT_background <- raw_data$`Antibody Capture`[, empty_drops]
  rownames(ADT_background) = gsub("_", "-", rownames(ADT_background))
  ADT_counts <- sobj_local[['ADT']]@counts
  #Preserving the raw ADT data
  sobj_local[['raw.ADT']] = sobj_local[['ADT']]

  stat_local[[ 'empty_drops_count' ]] <- length(empty_drops)

  adt_norm = DSBNormalizeProtein(
    cell_protein_matrix = ADT_counts,
    empty_drop_matrix = ADT_background,
    denoise.counts = TRUE,
    use.isotype.control = TRUE,
    isotype.control.name.vec = grep("isotype", rownames(ADT_counts), value = T, ignore.case = T)
  )

  adt_norm = apply(adt_norm,
                   2,
                   function(x){ ifelse(test = x < -10, yes = 0, no = x)})

  sobj_local[['ADT']] = CreateAssayObject(data = adt_norm)


  p = VlnPlot(sobj_local, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"), ncol = 4, pt.size = 0)
  print(p)

  print(dim(sobj_local))
  stat_local[[ "dim_raw" ]] <- dim(sobj_local)
  stat_local[[ "droplet.types" ]] = table(sobj_local$DROPLET.TYPE.FINAL)

  stat[[ s ]] = stat_local
  sobjs[[ s ]] = sobj_local
}

sobjs = readRDS(file.path(out_dir,paste0("sobjs_dsbNormedADT_20220303.rds")))

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes


## Merge data.
sobjs_all = merge( sobjs[[1]], y = sobjs[-1], project = "XHLT2", add.cell.ids = names(sobjs) )


## I see that the pool-specific effect comes only from the ADT data (pools 4 and 5 (preped by one technician of the three) are very different from other 5 pools). Hence, integrating the ADT data. Then correct the pool-effect in RNA using harmony, and use the PCA from integrated ADT data and harmony from RNA data for neighbors and UMAP/clustering.

##### Integrate ADT data across libraries (orid.ident)
DefaultAssay(sobjs_all) <- "ADT"
sobjs.list <- SplitObject(sobjs_all, split.by = "orig.ident")

# Select features that are repeatedly variable across datasets for integration
#features <- SelectIntegrationFeatures(object.list = sobjs.list)
# Above command fails since the normalized data is imported from DSB and not calculated inside Seurat. Hence just using all non-isotype markers for the following analyses.
features <- grep("isotype", rownames(sobjs_all[["ADT"]]), ignore.case = T, value = T, invert = T)

# Remove the isotype controls
features <- features[ features %in% grep("isotype", rownames(sobjs_all[["ADT"]]), ignore.case = T, value = T, invert = T) ]


# Using first library from each pool as reference to speed up the process. Also using "rpca" speeds up the process. https://github.com/satijalab/seurat/discussions/3999
# Testing out the original FindIntegrationAnchors call using parallel processing (disabled in my RStudio) using R at data/scRNA_CITE/DB1-7/adt_integration/

## For "rpca", Data needs to be scaled and PCA needs to be calculated for each individual object.
sobjs.list <- lapply(X = sobjs.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = TRUE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})


immune.anchors <- FindIntegrationAnchors(object.list = sobjs.list, anchor.features = features, reference = grep("SCG1", names(sobjs.list)), reduction = "rpca" )

# this command creates an 'integrated' data assay
sobjs_all_adt_corr <- IntegrateData(anchorset = immune.anchors, new.assay.name = "integrated.ADT")

#####saveRDS(sobjs_all_adt_corr, file = file.path(out_dir,paste0("sobjs_all_adt_corr_tmp_20220303.rds")))

#sobjs_all_adt_corr = readRDS(file.path(out_dir,paste0("sobjs_all_adt_corr_tmp_20220127.rds")))


# Calculate PCA for the integrated.ADT
DefaultAssay(sobjs_all_adt_corr) <- "integrated.ADT"
sobjs_all_adt_corr <- ScaleData(sobjs_all_adt_corr, vars.to.regress = c("nCount_ADT","nFeature_ADT","S.Score", "G2M.Score") ) %>%
                RunPCA(npcs = 30, verbose = FALSE, reduction.name = "int.apca")

# Perform batch correction in RNA using harmony
DefaultAssay(sobjs_all_adt_corr) <- "RNA"
sobjs_all_both_corr <- NormalizeData(sobjs_all_adt_corr) %>%
                FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
                CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE) %>%
                ScaleData(vars.to.regress = c("percent.mt","percent.ribo","nCount_RNA","nFeature_RNA","S.Score", "G2M.Score") ) %>%
                RunPCA(features = VariableFeatures(object = .), npcs=50) %>%
                RunHarmony("orig.ident")


# Calculate neighbors using harmony and int.apca
sobjs_all_both_corr <- FindMultiModalNeighbors(
  sobjs_all_both_corr, reduction.list = list("harmony", "int.apca"),
  dims.list = list(1:30, 1:18), modality.weight.name = "RNA.weight", weighted.nn.name = "weighted.nn", prune.SNN = 1/20)

sobjs_all_both_corr <- RunUMAP(sobjs_all_both_corr, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
sobjs_all_both_corr <- RunUMAP(sobjs_all_both_corr, reduction = "harmony", reduction.name = "harmony.umap", reduction.key = "harmony.UMAP_", dims = 1:30)
sobjs_all_both_corr <- RunUMAP(sobjs_all_both_corr, reduction = "int.apca", reduction.name = "int.aumap", reduction.key = "int.aUMAP_", dims = 1:18)


sobjs_all_both_corr = FindClusters(sobjs_all_both_corr, graph.name = "wsnn", algorithm = 3, resolution = 2, verbose = TRUE, n.start = 10)

saveRDS(sobjs_all_both_corr, file = file.path(out_dir,paste0("sobjs_all_both_corr_clust.rds")))


