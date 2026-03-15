library(Seurat)
library(tidyverse)
library(sva)


out_dir = "202008_.._XHLT2/data/scRNA_CITE/DB8-14/mapping_to_DB1-7/"
sobjs_all = readRDS(file.path(out_dir, "sobjs_all_mapped.rds"))



#' @param mat: feature (row) x sample (column) data matrix
#' @param batch: named vector indicating batch labels where names match the sample names
#' @param meta: sample metadata
run_pvca <- function(mat, batch, meta, plot_title) {
  # step 1: load required packages
  library(package = "pvca") # devtools::install_github("dleelab/pvca")
  VarCorr <- lme4::VarCorr

  # Match samples
  cmn = colnames(mat) #intersect(colnames(mat), rownames(meta))
  mat = mat[,cmn]
  meta = meta[cmn,]
  batch = batch[cmn]

  # step 2: specify which assay is going to be used as input
  pvca_input_matrix <- mat
  if (sum(is.na(pvca_input_matrix)) > 0) {
    stop("pvca input matrix should not have missing values")
  }

  pvca_input_phenodata = data.frame(batch = batch,
                         age = meta$demographicSurvey_age,
                         sex = ifelse( meta$demographicSurvey_genderFemale > meta$demographicSurvey_genderMale, "F", "M"),
                         row.names = colnames(mat))

  # step 5: run PVCA
  fit <- PVCA(counts    = pvca_input_matrix,
              meta      = pvca_input_phenodata,
              inter     = FALSE,
              threshold = 0.6)

  # step 6: plot barplot with PVCA results
  pvca_barplot_data <- data.frame(explained = as.vector(fit),
                                  effect    = names(fit)) %>%
    arrange(explained) %>%
    mutate(effect = factor(effect, levels = effect))

  pvca_barplot <- ggplot(data    = pvca_barplot_data,
                         mapping = aes(x = effect, y = explained)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = signif(explained, digits = 3)),
              nudge_y   = 0.01,
              size      = 3) +
    labs(x = NULL, y = "Proportion of the variance explained") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45,
                                     vjust = 1,
                                     hjust = 1)) +
    ggtitle(plot_title)
  print(pvca_barplot)
}


## count_mat is a matrix/data.frame of genes (rows) x samples (columns)
calc_tmm_data <- function(count_mat) {
  d = edgeR::DGEList(counts = count_mat)
  d = edgeR::calcNormFactors(d, method = "TMM")
  tmm = edgeR::cpm(count_mat, lib.size=d$samples$lib.size * d$samples$norm.factors)
  return(tmm)
}


# Keep only BLE1, AAE5, and AAE6 timepoints
sobjs_all = subset(sobjs_all, cells=Cells(sobjs_all)[grepl("BLE1|AAE6", sobjs_all$colabs_id)]  ) 
gc()

#### Batch correction before preparing matrices.
# Load batch and demographic information
scR_sample_pool_map = unique(
    sobjs_all@meta.data[, c("colabs_id","orig.ident")] %>% 
    mutate(pool = str_extract(orig.ident, "DB\\d+")) %>% 
    dplyr::select(c("colabs_id","pool"))
    ) %>% 
  pull(pool, name="colabs_id")
scR_sample_pool_map = scR_sample_pool_map[ ! is.na(names(scR_sample_pool_map)) ]


demographic = read.csv("202008_.._XHLT2/data/integration/demographics.csv", row.names = 1)
demographic = demographic[ 
                str_extract(names(scR_sample_pool_map), "XHLT2-HS\\d+"), 
                c("demographicSurvey_age","demographicSurvey_genderFemale","demographicSurvey_genderMale")
                ] %>% 
                `rownames<-`(names(scR_sample_pool_map))


# Either use sum of raw counts or use average of normalized counts across single-cells to calculate pseudobulk. In case of the raw counts, the Combat_seq is run on the raw counts and the corrected counts are cpm'ed.
useRawSum = FALSE
useNormAverage = ! useRawSum
old = FALSE ## If useNormAverage and old are true, the Combat_seq is used on AverageExpression output and the data is not log-transformed. If useNormAverage is true but old is false, Combat is used on AverageExpression output and the data is log-transformed. If both are false, the summed raw read counts will be batch-corrected using Combat_seq and the data is log-transformed.



pdf(file.path(out_dir,"scRNA_batch_correction_2tp.pdf"))
rna_matrices = list(preScale = list())
adt_matrices = list(preScale = list())
for(n in unique(sobjs_all$predicted.manual.celltype) ) {
  print(n)
  subs = subset(sobjs_all,
                cells = Cells(sobjs_all)[sobjs_all$predicted.manual.celltype == n ]
                )

  # Process only those cell types that have more than 20 cells in more than 10 samples.
  min_cell_per_donor = 20
  allowed_donors = names(which(table(subs$colabs_id) > min_cell_per_donor))
  if(length(allowed_donors) <= 10) {
    print(paste0("Skipping '", n, "' because of to few samples (<=10)"))
    next
  }
  subs = subset(subs, cells = Cells(subs)[subs$colabs_id %in% allowed_donors] )
  #######


  # Calculate pseudobulk expression
  if(useRawSum) {
    ## Aggregating the counts in sparse matrix format (adapted the code from Seurat's PseudobulkExpression function, this is extremely fast).
    data <- FetchData(object = subs, vars = rev(x = "colabs_id"))
    category.matrix <- Matrix::sparse.model.matrix(object = as.formula(object = paste0("~0+",paste0("data[,", 1, "]", collapse = ":"))))
    rna_mat = as.matrix((GetAssayData(subs, "count", "RNA")) %*% category.matrix)
    colnames(rna_mat) = gsub("data\\[, 1\\]","",colnames(rna_mat))
  } else if(useNormAverage) {
    rna = AverageExpression(subs, assays = "RNA", features = rownames(subs[['RNA']]), slot = "data", group.by = "colabs_id")
    rna_mat = rna$RNA
    if(old)
      warning("When useNormAverage and old are true, the averaged expression data isn't log-transformed. It should be log-transformed because the output of AverageExpression is in linear space.")
  }

 subs[['ADT']]@counts = subs[['ADT']]@data
  adt = AverageExpression(subs, assays = "ADT", features = rownames(subs[['ADT']]), slot = "counts", group.by = "colabs_id")
  adt_mat = adt$ADT

  # Get batch labels
  scR_sample_pool_map_n = scR_sample_pool_map[colnames(rna_mat)]
  print(table(scR_sample_pool_map_n))

  # Test if the number of genes captured are different across batches.
  if(useRawSum) {
    t = colSums(rna_mat > 0)
    pval = kruskal.test(t, scR_sample_pool_map_n)$p.value
    boxplot(t ~ scR_sample_pool_map_n, ylab="#genes", main=paste("ks.pval", formatC(pval, format = "e", digits =  2)))
  }

  # Remove pools with single sample
  one_sample_pool = names(which(table(scR_sample_pool_map_n) == 1))
  if( length(one_sample_pool) > 0 ) {
    singleton_sample = names(scR_sample_pool_map_n)[ scR_sample_pool_map_n %in% one_sample_pool ]
    scR_sample_pool_map_n = scR_sample_pool_map_n[ ! scR_sample_pool_map_n %in% one_sample_pool ]
    rna_mat = rna_mat[, !colnames(rna_mat) %in% singleton_sample ]
    adt_mat = adt_mat[, !colnames(adt_mat) %in% singleton_sample ]
  }


  ## Remove lowly expressed features: keep features expressed in at least 30% samples
  rna_mat <- rna_mat[ rowSums(rna_mat > 0) > 0.3*ncol(rna_mat) , ]
  adt_mat <- adt_mat[ rowSums(adt_mat > 0) > 0.3*ncol(adt_mat) , ]


  # Skip the cell type if the demographics information do not have enough levels for the final set of samples. This happens only in cases where the final set of samples includes samples from a very small number of participants.
  skip_ctype_due_to_lack_of_diversity = FALSE
  for(i in 1:ncol(demographic)) {
    if(length(table(demographic[colnames(rna_mat), i])) == 1)
      skip_ctype_due_to_lack_of_diversity = TRUE
  }
  if(skip_ctype_due_to_lack_of_diversity) {
    print(paste0("Skipping '", n, "' because of lack of demographic diversity"))
    next
  }


  #### Perform Batch-correction
  tmp_mat = rna_mat
  if(useRawSum)
    tmp_mat = calc_tmm_data(tmp_mat)
  run_pvca(tmp_mat, scR_sample_pool_map, demographic, paste0(n, " - RNA.before corr"))

  if(useRawSum | old) {
    rna_corrected = ComBat_seq(rna_mat, scR_sample_pool_map_n)
  } else if(useNormAverage & !old) {
    rna_corrected = ComBat(rna_mat, scR_sample_pool_map_n)
  }

  tmp_mat = rna_corrected
  if(useRawSum)
    tmp_mat = calc_tmm_data(tmp_mat)
  run_pvca(tmp_mat, scR_sample_pool_map, demographic, paste0(n, " - RNA.after corr"))

  if(useRawSum) {
    rna_corrected = calc_tmm_data(rna_corrected)
  }

  # Remove scRNA features correlating with scRNA pool
  tmp_ds <- rna_corrected
  pool <- scR_sample_pool_map_n
  pvals <- apply(tmp_ds, 1, function(x) kruskal.test(x, pool )$p.value )
  pvals_adj <- p.adjust(pvals, method = "BH")
  keep_feature <- pvals_adj > 0.01
  keep_feature[is.na(keep_feature)] <- FALSE
  rna_corrected <- tmp_ds[keep_feature,]
  cat("Removed", nrow(tmp_ds) - nrow(rna_corrected), "out of", nrow(tmp_ds), "features for", n, "\n")


  run_pvca(adt_mat, scR_sample_pool_map, demographic, paste0(n, " - ADT.before corr"))
  adt_corrected = ComBat(adt_mat, scR_sample_pool_map[colnames(adt_mat)])
  run_pvca(adt_corrected, scR_sample_pool_map, demographic, paste0(n, " - ADT.after corr"))

  # Select top ngenes variable genes across the cohort.
  ngenes = 100000 # Take all genes
  madv = apply(rna_corrected, 1, mad)
  avgv = apply(rna_corrected, 1, mean)
  var_genes = names(tail(sort(madv), ngenes))
  mad_cutoff = min(tail(sort(madv), ngenes))

  p = data.frame(avgv, madv) %>% dplyr::filter( ! is.na(madv) & madv != 0 ) %>% ggplot(aes(log1p(avgv), log1p(madv))) + geom_point(size=0.1, alpha=0.5) + theme_classic() + ggtitle(n) + geom_hline(yintercept = log1p(mad_cutoff))
  print(p)

  # Log transformation after MADs are calculated (log transformation flips the mean vs MAD correlations, and we ended up picking only the lowly expressed features)
  if(useRawSum) {
    rna_corrected = log(rna_corrected+1, 10)
  }
  if(!old) {
    min_val = min(rna_corrected)
    if(min_val < 0)
      rna_corrected = rna_corrected + abs(min_val) # If there are negative values, add the non-negative of that value to the whole matrix to make the smallest value to be 0.
    rna_corrected = log(rna_corrected+1, 10)
  }

  # Subset expression data and scale the data.
  #rna_matrices[[n]] = rna$RNA[var_genes,] %>% as.matrix()
  #adt_matrices[[n]] = adt$ADT %>% as.matrix()

  rna_matrices[[n]] = apply(rna_corrected[var_genes,], 1, scale) %>% t() %>% data.frame() %>% setNames(colnames(rna_corrected)) %>% as.matrix()
  rna_matrices[["preScale"]][[n]] = rna_corrected
  adt_matrices[[n]] = apply(adt_corrected, 1, scale) %>% t() %>% data.frame() %>% setNames(colnames(adt_corrected)) %>% as.matrix()
  adt_matrices[["preScale"]][[n]] = adt_corrected

  run_pvca(rna_matrices[[n]], scR_sample_pool_map, demographic, paste0(n, " - RNA final var. feat."))

  cat("RNA dim", dim(rna_matrices[[n]]), "\n", "ADT dim", dim(adt_matrices[[n]]), "\n" )


  saveRDS(rna_matrices, file.path(out_dir, "allgenes_pseudobulk_rna_followup.rds" ) )
  saveRDS(adt_matrices, file.path(out_dir, "allgenes_pseudobulk_adt_followup.rds" ) )
}
dev.off()



