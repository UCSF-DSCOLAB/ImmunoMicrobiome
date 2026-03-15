library(tidyverse)
library(Azimuth)
library(Seurat)
# The seurat object of baseline data was used to make Azimuth reference dataset, and the longitudinal scRNA-seq data was mapped this reference using Azimuth
out_dir = "202008_.._XHLT2/data/scRNA_CITE/DB8-14/azimuth_ref_DB1-7/"
sobjs = readRDS("202008_.._XHLT2/data/integration/final_obj_files_20230105/sobjs_all_both_corr_clust.rds")
sobjs = RunUMAP(sobjs, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_", return.model=T)
sobjs = SCTransform(sobjs, assay='RNA', new.assay.name="SCT")
#saveRDS(object = sobjs, file = file.path("fullref.Rds"))
ref = AzimuthReference(sobjs, refUMAP="wnn.umap", refDR="pca", refAssay="SCT", metadata="manual.celltype", dims=1:50, k.param=31, reference.version="1.0.0", plotref = "wnn.umap")
SaveAnnoyIndex(object = ref[["refdr.annoy.neighbors"]], file = file.path(out_dir, "idx.annoy"))
saveRDS(object = ref, file = file.path(out_dir, "ref.Rds"))

