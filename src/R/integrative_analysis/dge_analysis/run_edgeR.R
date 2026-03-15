library(Seurat)
library(tidyverse)
library(edgeR)
sobj = readRDS( file.path("sobjs_all_both_corr_clust.rds"))   

obj = readRDS("mofa_model_final_object_v1.rds")


# Get f1-top, f1-bottom, f3-top, and f3-bottom donors
f1_donor_bin1 = names(head(sort(na.omit(obj$factors.noOutl[,1])), round(length(na.omit(obj$factors.noOutl[,1]))*0.2) ))
f1_donor_bin5 = names(tail(sort(na.omit(obj$factors.noOutl[,1])), round(length(na.omit(obj$factors.noOutl[,1]))*0.2) ))
f3_donor_bin1 = names(head(sort(na.omit(obj$factors.noOutl[,3])), round(length(na.omit(obj$factors.noOutl[,3]))*0.2) ))
f3_donor_bin5 = names(tail(sort(na.omit(obj$factors.noOutl[,3])), round(length(na.omit(obj$factors.noOutl[,3]))*0.2) ))

# Get f3-top vs f1-top donors
n_donor=20
f1_pos_ids = obj$factors.noOutl %>% as.data.frame() %>% rownames_to_column("donor_id") %>% arrange(Factor1) %>% pull(donor_id) %>% tail(n_donor)
f3_pos_ids = obj$factors.noOutl %>% as.data.frame() %>% rownames_to_column("donor_id") %>% arrange(Factor3) %>% pull(donor_id) %>% tail(n_donor)
# Remove the overlapping ones.
f1_pos_ids_new = f1_pos_ids[!f1_pos_ids %in% f3_pos_ids]
f3_pos_ids_new = f3_pos_ids[!f3_pos_ids %in% f1_pos_ids]

# Get scRNA-seq pool ids
scR_sample_pool_map <- read.table("202008_.._XHLT2/data/integration/final_obj_files_20230105/sample_to_pool_mapping.tsv", header = T) %>% mutate(sample = paste0("XHLT2-",sample) ) %>% pull(pool, name = "sample")
scR_sample_pool_map_df <- as.data.frame(scR_sample_pool_map)

ctypes = names(obj$dataset_list)[7:42]

comparisons = list(
  "f1top_f1btm" = list(f1_donor_bin5, f1_donor_bin1),
  "f3top_f3btm" = list(f3_donor_bin5, f3_donor_bin1),
  "f3top_f1top" = list(f3_pos_ids_new, f1_pos_ids_new)
  )

dge_res = list()
for(ctype in ctypes) {
  print(ctype)
  dge_res[[ctype]] = list()
  sub = subset(sobj, subset = manual.celltype == ctype)

  # Use only those samples that have more than 20 cells.
  min_cell_per_donor = 20
  allowed_donors = names(which(table(sub$donor_id) > min_cell_per_donor))
  sub = subset(sub, cells = Cells(sub)[sub$donor_id %in% allowed_donors] )

  pb_counts <- t(rowsum(t(as.matrix(sub[['RNA']]@counts)), sub$donor_id))
  pb_counts <- pb_counts[, ! colnames(pb_counts) %in% c("XHLT2-HS33","XHLT2-HS99") ]
  for(comp in names(comparisons)) {
    print(comp)
    num_ids = intersect( comparisons[[comp]][[1]] , colnames(pb_counts))
    den_ids = intersect( comparisons[[comp]][[2]] , colnames(pb_counts))
    # Process only those comparisons where there are at least 2 samples per group
    if(length(num_ids) <= 1 | length(den_ids) <= 1) {
      print(paste0("Skipping, because do not have enough samples for ", comp, " for ", ctype))
      next
    }

    pb_sample_anno <- data.frame(
        sample = c(num_ids, den_ids),
        group = c(rep("num",length(num_ids)), rep("den",length(den_ids))),
        batch = scR_sample_pool_map[c(num_ids, den_ids)]
    )

    pb_dge <- DGEList(
        counts = pb_counts[,c(num_ids, den_ids)],
        samples = pb_sample_anno,
        group = pb_sample_anno$cell_type
    )

    pb_dge <- calcNormFactors(pb_dge)

    idx = rowSums(cpm(pb_dge) > 0) > 0.30*ncol(pb_dge)
    pb_dge = pb_dge[idx,]

    design <- model.matrix(~ batch + group, data = pb_dge$samples)

    pb_dge <- estimateDisp(pb_dge, design)

    pb_fit <- glmQLFit(pb_dge, design)
    pb_lrt <- glmQLFTest(pb_fit, coef = "groupnum")
    dge_res[[ctype]][[comp]] = topTags(pb_lrt, Inf)
  }
}
saveRDS(dge_res, "202008_.._XHLT2/data/integration/dge_analysis/dge_res.rds")

