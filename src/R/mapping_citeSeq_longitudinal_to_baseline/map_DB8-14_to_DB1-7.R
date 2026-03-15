library(tidyverse)
library(Azimuth)
library(Seurat)
fs = list.files("data/longitudinal/single_cell_GEX/", pattern="seurat_object_qcAndAddModalities.rds", recursive=T, full.names=T)
mapping = read.csv("fmlCluster_to_sample_mapping.csv")
mapping$key = paste0(mapping$BEST.GUESS, mapping$orig.ident)
sobjs = list()
for (f in fs) {
  id = str_extract(f, "XHLT2-POOL-DB\\d+-SCG\\d+")
  print(id)
  obj = readRDS(f)
  obj = RunAzimuth(obj, reference = "202008_.._XHLT2/data/scRNA_CITE/DB8-14/azimuth_ref_DB1-7/")
  obj$key = paste0(obj$BEST.GUESS, obj$orig.ident)
  obj$donor_id = mapping$donor_id[ match(obj$key, mapping$key) ]
  obj$colabs_id = mapping$colabs_id[ match(obj$key, mapping$key) ]
  sobjs[[id]] = obj
}
sobjs_all = merge( sobjs[[1]], y = sobjs[-1], project = "XHLT2", add.cell.ids = names(sobjs) )
saveRDS(sobjs_all, "202008_.._XHLT2/data/scRNA_CITE/DB8-14/azimuth_ref_DB1-7/sobjs_all_mapped.rds") # The mapped annotation is available on GEO

