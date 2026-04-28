#```{r functions}

run_pca <- function(data, return_loadings=T, title="title", draw_plots = T) {

  pca=prcomp(data)
  pca_out=as.data.frame(pca$x)
  
  percentVar = round((pca$sdev)^2 / sum(pca$sdev^2)*100) %>% setNames(colnames(pca_out))
  p1 = ggplot(pca_out, aes(PC1, PC2)) + geom_point(size=2) +
    xlab(paste0("PC1: ",percentVar[1],"% variance")) +
    ylab(paste0("PC2: ",percentVar[2],"% variance")) +
    #coord_fixed(ratio = 1) +
    #ggtitle("Biogenic amine PCA") +
    guides(colour = guide_legend(override.aes = list(shape = 15))) +
    theme_minimal() + 
    ggtitle(title)
  
  if(draw_plots) {
    print(p1)
    barplot(percentVar)
  }
  
  if(return_loadings)
    return( list( "pca_out"=pca_out, "percentVar"=percentVar, plot = p1) )
  
  if(! return_loadings)
    return(NULL)

}



##### PVCA analysis
#' @param mat: feature (row) x sample (column) data matrix
#' @param batch: named vector or data frame indicating batch labels where names or rownames match the sample names
#' @param meta: sample metadata, where the rownames match the sample names
run_pvca <- function(mat, batch, meta, plot_title = "") {  
  # step 1: load required packages
  library(package = "pvca") # devtools::install_github("dleelab/pvca")
  VarCorr <- lme4::VarCorr
  
  if(! is.data.frame(batch)) {
    batch = as.data.frame(batch) %>% setNames("batch")
  }
  
  # Match samples
  cmn = intersect(intersect(colnames(mat), rownames(meta)), rownames(batch))
  mat = mat[,cmn]
  meta = meta[cmn,]
  batch = batch[cmn,]
  
  # step 2: specify which assay is going to be used as input
  pvca_input_matrix <- mat
  if (sum(is.na(pvca_input_matrix)) > 0) {
    stop("pvca input matrix should not have missing values")
  }
  
  #meta = meta[colnames(mat),]
  #meta = apply(meta, 2, function(x){ if(is.numeric(x)){ x[is.na(x)] = 999; x } else { x[is.na(x)] = "NA"; x }  }) %>% as.data.frame()
  pvca_input_phenodata = data.frame(age = meta$demographicSurvey_age,
                         sex = ifelse( meta$demographicSurvey_genderFemale > meta$demographicSurvey_genderMale, "F", "M"),
                         row.names = colnames(mat))
  pvca_input_phenodata = cbind(batch, pvca_input_phenodata)
  
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






# Same as plot_data_heatmap function from MOFA2, except that this function doesn't make a plot but just returns the matrix used for pheatmap
get_heatmap_data = function (object, factor, view = 1, groups = "all", features = 50, 
    annotation_features = NULL, annotation_samples = NULL, transpose = FALSE, 
    imputed = FALSE, denoise = FALSE, max.value = NULL, min.value = NULL, 
    ...) 
{
    if (!is(object, "MOFA")) 
        stop("'object' has to be an instance of MOFA")
    stopifnot(length(factor) == 1)
    stopifnot(length(view) == 1)
    groups <- MOFA2:::.check_and_get_groups(object, groups)
    factor <- MOFA2:::.check_and_get_factors(object, factor)
    view <- MOFA2:::.check_and_get_views(object, view)
    W <- do.call(rbind, get_weights(object, views = view, factors = factor, 
        as.data.frame = FALSE))
    Z <- lapply(get_factors(object)[groups], function(z) as.matrix(z[, 
        factor]))
    Z <- do.call(rbind, Z)[, 1]
    Z <- Z[!is.na(Z)]
    if (isTRUE(denoise)) {
        data <- predict(object, views = view, groups = groups)[[1]]
    }
    else {
        if (isTRUE(imputed)) {
            data <- get_imputed_data(object, view, groups)[[1]]
        }
        else {
            data <- get_data(object, views = view, groups = groups)[[1]]
        }
    }
    if (is(data, "list")) {
        data <- do.call(cbind, data)
    }
    if (is(features, "numeric")) {
        if (length(features) == 1) {
            features <- rownames(W)[tail(order(abs(W)), n = features)]
        }
        else {
            features <- rownames(W)[order(-abs(W))[features]]
        }
        features <- names(W[features, ])[order(W[features, ])]
    }
    else if (is(features, "character")) {
        stopifnot(all(features %in% features_names(object)[[view]]))
    }
    else {
        stop("Features need to be either a numeric or character vector")
    }
    data <- data[features, ]
    data <- data[, names(Z)]
    data <- data[, apply(data, 2, function(x) !all(is.na(x)))]
    order_samples <- names(sort(Z, decreasing = TRUE))
    order_samples <- order_samples[order_samples %in% colnames(data)]
    data <- data[, order_samples]
    if (!is.null(annotation_samples)) {
        if (is.data.frame(annotation_samples)) {
            message("'annotation_samples' provided as a data.frame, please make sure that the rownames match the sample names")
            if (any(!colnames(data) %in% rownames(annotation_samples))) {
                stop("There are rownames in annotation_samples that do not correspond to sample names in the model")
            }
            annotation_samples <- annotation_samples[colnames(data), 
                , drop = FALSE]
        }
        else if (is.character(annotation_samples)) {
            stopifnot(annotation_samples %in% colnames(object@samples_metadata))
            tmp <- object@samples_metadata
            rownames(tmp) <- tmp$sample
            tmp$sample <- NULL
            tmp <- tmp[order_samples, , drop = FALSE]
            annotation_samples <- tmp[, annotation_samples, drop = FALSE]
            rownames(annotation_samples) <- rownames(tmp)
        }
        else {
            stop("Input format for 'annotation_samples' not recognised ")
        }
        foo <- sapply(annotation_samples, function(x) is.logical(x) || 
            is.character(x))
        if (any(foo)) 
            annotation_samples[, which(foo)] <- lapply(annotation_samples[, 
                which(foo), drop = FALSE], as.factor)
    }
    if (!is.null(annotation_features)) {
        stop("'annotation_features' is currently not implemented")
    }
    if (transpose) {
        data <- t(data)
        if (!is.null(annotation_samples)) {
            annotation_features <- annotation_samples
            annotation_samples <- NULL
        }
        if (!is.null(annotation_features)) {
            annotation_samples <- annotation_features
            annotation_features <- NULL
        }
    }
    if (!is.null(max.value)) 
        data[data >= max.value] <- max.value
    if (!is.null(min.value)) 
        data[data <= min.value] <- min.value
    return(data)
    #pheatmap(data, annotation_row = annotation_features, annotation_col = annotation_samples,...)
}


flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
    )
}

enrichment_dotplot <- function(comb_res, pval.adj.cutoff, max_pathway = 50) {
  all_stat = data.frame()
  all_pval = data.frame()
  all_sig_path = c()
  for(ds in names(comb_res)) {
    if(class(comb_res[[ds]]) == "gseaResult") {
      stat_tmp = comb_res[[ds]]@result[,c("Description","NES"),drop=F] %>% setNames(c("rn","NES"))
      pval_tmp = comb_res[[ds]]@result[,c("Description","p.adjust"),drop=F] %>% setNames(c("rn","padj"))
    } else {
      stat_tmp = comb_res[[ds]]$set.statistics %>% as.data.frame() %>% rownames_to_column("rn")
      pval_tmp = comb_res[[ds]]$pval.adj %>% as.data.frame() %>% setNames("padj") %>% rownames_to_column("rn")
    }
    if(nrow(all_stat) == 0) {
      all_stat = stat_tmp
      all_pval = pval_tmp
    } else {
      all_stat = full_join(all_stat, stat_tmp, by = "rn")
      all_pval = full_join(all_pval, pval_tmp, by = "rn")
    }
    all_sig_path = unique(c(all_sig_path, pval_tmp$rn[ pval_tmp$padj <= pval.adj.cutoff ] ))
  }
  all_stat = all_stat %>% column_to_rownames("rn") %>% setNames(names(comb_res))
  all_pval = all_pval %>% column_to_rownames("rn") %>% setNames(names(comb_res))
  
  all_stat = all_stat[all_sig_path,]
  all_pval = all_pval[all_sig_path,]
  
  
  
  driving_assays = rep(1, length(all_pval))
  all_pval$joint = apply(all_pval,1,function(z){
    ii = !is.na(z)
    t = qnorm(z[ii])
    t = sum(t*driving_assays[ii])/sqrt(sum(driving_assays[ii]^2))
    pnorm(t, lower.tail = T)
  })
  all_stat$joint = 0
  
  
  pathway_levels = row.names(all_pval %>% arrange(joint))
  

  all_stat = pivot_longer(all_stat %>% rownames_to_column("rn"), cols = -"rn") %>% setNames(c("pathway","dataset","stat"))
  all_pval = pivot_longer(all_pval %>% rownames_to_column("rn"), cols = -"rn") %>% setNames(c("pathway","dataset","pval"))
  all_stat$pval.adj = all_pval$pval
  
  
  all_stat$pathway = factor(all_stat$pathway, levels = rev(pathway_levels))
  all_stat = all_stat %>% dplyr::filter(pathway %in% head(pathway_levels, max_pathway) )
  p = ggplot(all_stat, aes(dataset, pathway, size=-log10(pval.adj), color=as.character(sign(stat)))) + 
    geom_point(shape=1) + scale_color_manual(values = c("1"="red","-1"="blue")) + 
    geom_point(data = all_stat %>% dplyr::filter(pval.adj < pval.adj.cutoff)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle=90))
  print(p)
}



run_mHG <- function (signif_feat_names, mofa_weights_vec, Term2Gene, backround = "im", p.adjust.methods = "BH", minGSSize=2)
{
    mHG_result = list()
    # Positive-end
    signif_feat_names_pos = names(which(mofa_weights_vec[ signif_feat_names ] > 0))
    geneList = mofa_weights_vec[signif_feat_names_pos]
    geneList = abs(geneList)
    geneList = sort(geneList, decreasing = T)
    universe = NULL
    if (backround == "im") {
        universe = names(mofa_weights_vec)
    } else {
        universe = unique(Term2Gene[, 2])
    }
    if (length(geneList) > 0) {
        mHG_result[["pos"]] = mgh_test_wrapper(geneList = geneList, 
          minGSSize = minGSSize, maxGSSize = 500, pvalueCutoff = 2, 
          universe = universe, TERM2GENE = Term2Gene, 
          p.adjust.methods = p.adjust.methods)
    }

    # Negative-end
    signif_feat_names_neg = names(which(mofa_weights_vec[ signif_feat_names ] < 0))
    geneList = mofa_weights_vec[signif_feat_names_neg]
    geneList = abs(geneList)
    geneList = sort(geneList, decreasing = T)
    universe = NULL
    if (backround == "im") {
        universe = names(mofa_weights_vec)
    } else {
        universe = unique(Term2Gene[, 2])
    }
    if (length(geneList) > 0) {
        mHG_result[["neg"]] = mgh_test_wrapper(geneList = geneList, 
          minGSSize = minGSSize, maxGSSize = 500, pvalueCutoff = 2, 
          universe = universe, TERM2GENE = Term2Gene, 
          p.adjust.methods = "BH")
    }
    
    if(! is.null(mHG_result$pos))
      mHG_result$pos = mHG_result$pos %>% 
        dplyr::select("Description","pvalue","setSize","rank","p.adj","core_enrichment") %>% 
        mutate(NES = 1, qvalues = p.adj)
    if(! is.null(mHG_result$neg))
      mHG_result$neg = mHG_result$neg %>% 
        dplyr::select("Description","pvalue","setSize","rank","p.adj","core_enrichment") %>% 
        mutate(NES = -1, qvalues = p.adj)
    
    if(length(mHG_result) == 2) {
      mHG_result$combined = rbind(mHG_result$pos, mHG_result$neg) %>% 
        group_by(Description) %>% 
        slice_min(pvalue, n=1, with_ties = F) %>%
        ungroup() %>% 
        mutate(p.adj = p.adjust(pvalue, p.adjust.methods), qvalues = p.adj) %>%
        arrange(p.adj)
    } else if(length(mHG_result) == 1) {
      mHG_result$combined = mHG_result[[1]]
    }
    
    if(! is.null(mHG_result))
      return(mHG_result)
}


#' @description Wrapper for running minimum hypergeometric (mHG) tests. Selected features (e.g., significantly DE genes, features significantly associated with a factor, etc) are ordered by some feature statistics, and enrichment of genesets is calculated. Only positive values can be used as feature statistics. If using log2FC as feature statistics, make sure to run the analysis once for the positive log2FC and second time for negative log2FC (after abs(log2FC)).  The results can be combined by using the direction with strongest pvalue for each geneset.
#' @param geneList: named vector of feature statistics (e.g., log2FC) for selected (e.g., significantly DE) features. names = feature names. Vector must be ordered in the decreasing order.
#' @param universe: vector of background features, e.g., all features that were measured in the analysis
#' @param TERM2GENE: data.frame with two columns, first = geneset name, second = features in that geneset. Columns names can be anything.
#' @param minGSSize: minimum size of the geneset to be used
#' @param maxGSSize: maximum size of the geneset to be used
#' @param pvalueCutoff: Ignore; currently not in use
#' @param p.adjust.methods: method for pvalue adjustment for multiple testing. See p.adjust for accepted methods.
mgh_test_wrapper <- function(geneList,
                            universe,
                            TERM2GENE,
                            minGSSize = 2,
                            maxGSSize = 500,
                            pvalueCutoff = 2,
                            p.adjust.methods="BH"){
  TERM2GENE=TERM2GENE[TERM2GENE[,2]%in%universe,,drop = F]
  genes = names(geneList)
  genes = intersect(genes, universe)
  ##keep only terms satisfying the size restriction
  size_counts = table(TERM2GENE[,1])
  terms_keep = names(size_counts[size_counts>=minGSSize & size_counts <=maxGSSize])
  TERM2GENE = TERM2GENE[TERM2GENE[,1]%in%terms_keep,,drop = F]
  pathway_names = unique(TERM2GENE[,1])
  genes_all = c(genes, setdiff(universe,genes))
  result_table = data.frame(ID = pathway_names)
  result_table[["Description"]] = pathway_names
  result_table[["setSize"]] = NA
  result_table[["pvalue"]] = NA
  result_table[["p.adj"]] = NA
  result_table[["rank"]] = NA
  result_table[["coreSetSiz"]] = NA
  result_table[["core_enrichment"]] = NA
  m = length(genes)
  for(pathway_id in 1:length(pathway_names)){
    pathway_name0 = pathway_names[pathway_id]
    term_genes = TERM2GENE[TERM2GENE[,1] ==pathway_name0,2]
    result_table[["setSize"]][pathway_id] = length(term_genes)
    mgh_vec = rep(0,length(genes_all))
    names(mgh_vec) = genes_all
    mgh_vec[term_genes]=1
    mHGtest_res = mHG::mHG.test(lambdas = mgh_vec, n_max = m)
    result_table[["pvalue"]][pathway_id] = mHGtest_res$p.value
    result_table[["rank"]][pathway_id] = mHGtest_res$n
    result_table[["coreSetSiz"]][pathway_id] = mHGtest_res$b
    result_table[["core_enrichment"]][pathway_id] = paste0(names(mgh_vec)[mgh_vec==1][1:mHGtest_res$b],collapse = "#")
  }
  result_table=result_table[order(result_table[["pvalue"]]),]
  result_table[["p.adj"]] = stats::p.adjust(p =result_table[["pvalue"]] , method = p.adjust.methods)
  #return(list(vec=mgh_vec, res=mHGtest_res ))
  ##turn to GSEA class object?
  return(result_table)
}






feature_selectionI = function(factor_loadings_one,  factor_loading_pval_one, FEWR = 0.05, edge_thr = 0.3, percent_thr_low = c(0.1, 0.1, 0.1, 0.05, 0.01, 0.01), use_bonferroni = FALSE, use_fdr = FALSE, top_feat = NULL, top_feat_by_direction = NULL, return_pvals = F){
  if(!is.null(top_feat) & !is.null(top_feat_by_direction)) {
    stop("Only one of top_feat and top_feat_by_direction can be defined. Both are not allowed.")
  }
  message("Note: It is assumed that the features in factor_loadings_one and factor_loading_pval_one are in the same order.")
  selected_list = list()
  for(d in names(factor_loadings_one)){
    tmp = factor_loadings_one[[d]][,1]
    names(tmp) = rownames(factor_loadings_one[[d]])
    tmp1 = factor_loading_pval_one[[d]][,1]
    names(tmp1) = rownames(factor_loading_pval_one[[d]])
    #tmp1 = tmp1[names(tmp)]
    if(use_bonferroni)
      tmp1 = tmp1 * length(tmp1)
    if(use_fdr)
      tmp1 = p.adjust(tmp1, method = "BH")
    idx1 = which(tmp1<=FEWR)
    #if(use_raw_pvals)
    #  idx1 = which(tmp1<=FEWR)
    #edge_thr1 = min(edge_thr,quantile(abs(tmp),1-percent_thr_low[d]))
    edge_thr1 = edge_thr
    idx2 = which(abs(tmp)>edge_thr1)
    if(! is.null(top_feat)) {
      tmp_feat = names(head(sort(tmp1), top_feat))
      idx1 = which(names(tmp1) %in% tmp_feat)
      idx2 = 1:length(tmp1)
      #prob = min(top_feat/length(tmp1), 1)
      #top_thr=quantile(tmp1, prob)
      #idx1 = which(tmp1 <= top_thr)
      #idx2 = 1:length(tmp1)
    }
    if(! is.null(top_feat_by_direction)) {
      tmp_feat = c( names(head(sort(tmp), top_feat_by_direction)),
                    names(tail(sort(tmp), top_feat_by_direction)))
      idx1 = which(names(tmp1) %in% tmp_feat)
      idx2 = 1:length(tmp1)
    }
    idx = intersect(idx1, idx2)
    selected_list[[d]] = tmp[idx]
    if(return_pvals) {
      selected_list[[d]] = data.frame(weight=tmp[idx], pval=tmp1[idx], is.sig = tmp1[idx] <= FEWR & abs(tmp[idx])>edge_thr1)
    }
  }
  return(selected_list)
}



enrichment_dotplot_v2 <- function(comb_res, pval.adj.cutoff = NULL, max_pathway = NULL, topN_per_dataset = NULL, title = "") {
  all_stat = data.frame()
  all_pval = data.frame()
  all_sig_path = c()
  for(ds in names(comb_res)) {
    if(class(comb_res[[ds]]) == "gseaResult") {
      stat_tmp = comb_res[[ds]]@result[,c("Description","NES"),drop=F] %>% setNames(c("rn","NES"))
      pval_tmp = comb_res[[ds]]@result[,c("Description","p.adjust"),drop=F] %>% setNames(c("rn","padj"))
    } else if(class(comb_res[[ds]]) == "list") {
      #comb_res[[ds]]$combined = comb_res[[ds]]$combined[ comb_res[[ds]]$combined$Description != "\n", ]
      if(nrow(comb_res[[ds]]$combined) == 0) {
        stat_tmp = data.frame(rn = "tmp", NES = NA)
        pval_tmp = data.frame(rn = "tmp", padj = NA)
      } else {
        stat_tmp = comb_res[[ds]]$combined[,c("Description","NES"),drop=F] %>% setNames(c("rn","NES"))
        pval_tmp = comb_res[[ds]]$combined[,c("Description","p.adj"),drop=F] %>% setNames(c("rn","padj"))
      }
    } else {
      #stat_tmp = comb_res[[ds]]$set.statistics %>% as.data.frame() %>% rownames_to_column("rn")
      #pval_tmp = comb_res[[ds]]$pval.adj %>% as.data.frame() %>% setNames("padj") %>% rownames_to_column("rn")
      stat_tmp = comb_res[[ds]][,c("Description","NES"),drop=F] %>% setNames(c("rn","NES"))
      pval_tmp = comb_res[[ds]][,c("Description","p.adj"),drop=F] %>% setNames(c("rn","padj"))
    }
    if(nrow(all_stat) == 0) {
      all_stat = stat_tmp
      all_pval = pval_tmp
    } else {
      all_stat = full_join(all_stat, stat_tmp, by = "rn")
      all_pval = full_join(all_pval, pval_tmp, by = "rn")
    }
    if(!is.null(pval.adj.cutoff)) {
      all_sig_path = unique(c(all_sig_path, pval_tmp$rn[ pval_tmp$padj <= pval.adj.cutoff ] ))
    }
    if(!is.null(topN_per_dataset)) {
      all_sig_path = unique(c(all_sig_path, pval_tmp %>% dplyr::filter(padj <= pval.adj.cutoff) %>% arrange(padj) %>% head(topN_per_dataset) %>% pull(rn) ))
    }
  }
  
  # Remove "NA" values for path that are introduced due to this condition: (if(nrow(comb_res[[ds]]$combined) == 0)) above.
  all_sig_path = all_sig_path[ ! is.na(all_sig_path) ]
  if(length(all_sig_path) == 0) {
    return(ggplot() + theme_void())
  }
  
  dup = all_stat$rn[ duplicated(all_stat$rn) ]
  cat("Found", length(dup), "duplicate pathway ids.\n")
  print(head(dup))
  all_stat = all_stat[!duplicated(all_stat$rn),]
  all_pval = all_pval[!duplicated(all_pval$rn),]


  all_stat = all_stat %>% remove_rownames() %>% column_to_rownames("rn") %>% setNames(names(comb_res))
  all_pval = all_pval %>% remove_rownames() %>% column_to_rownames("rn") %>% setNames(names(comb_res))
  
  all_stat = all_stat[all_sig_path,]
  all_pval = all_pval[all_sig_path,]
  
  
  
  driving_assays = rep(1, length(all_pval))
  all_pval$joint = apply(all_pval,1,function(z){
    ii = !is.na(z)
    t = qnorm(z[ii])
    t = sum(t*driving_assays[ii])/sqrt(sum(driving_assays[ii]^2))
    pnorm(t, lower.tail = T)
  })
  all_stat$joint = 0
  
  
  if(!is.null(max_pathway))
    pathway_levels = all_pval %>% rownames_to_column("rn") %>% arrange(joint) %>% pull(rn)
  
  if(!is.null(topN_per_dataset))
    pathway_levels = all_sig_path


  all_stat = pivot_longer(all_stat %>% rownames_to_column("rn"), cols = -"rn") %>% setNames(c("pathway","dataset","stat"))
  all_pval = pivot_longer(all_pval %>% rownames_to_column("rn"), cols = -"rn") %>% setNames(c("pathway","dataset","pval"))
  all_stat$pval.adj = all_pval$pval
  
  
  all_stat$pathway = factor(as.vector(all_stat$pathway), levels = rev(pathway_levels))
  all_stat$dataset = factor(as.vector(all_stat$dataset), levels = c(names(comb_res), "joint"))
  if(!is.null(max_pathway))
    all_stat = all_stat %>% dplyr::filter(pathway %in% head(pathway_levels, max_pathway) )
  p = ggplot(all_stat, aes(dataset, pathway, size=-log10(pval.adj), color=as.character(sign(stat)))) + 
    geom_point(shape=1) + scale_color_manual(values = c("1"="red","-1"="blue")) + 
    geom_point(data = all_stat %>% dplyr::filter(pval.adj < pval.adj.cutoff)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
    ggtitle(title)
  return(p)
}



# Using Cauchy method for joint p-val calculation. Cauchy method is used on pvalues to calculate joint pvalues, which are then converted to padj using "BH"
enrichment_dotplot_v3 <- function(comb_res, pval.adj.cutoff = NULL, max_pathway = NULL, topN_per_dataset = NULL, title = "", cauchy_thr = 0.8) {
  all_stat = data.frame()
  all_pval = data.frame()
  all_sig_path = c()
  all_sig_top_path = c()
  for(ds in names(comb_res)) {
    #if(class(comb_res[[ds]]) == "gseaResult") {
    #  stat_tmp = comb_res[[ds]]@result[,c("Description","NES"),drop=F] %>% setNames(c("rn","NES"))
    #  pval_tmp = comb_res[[ds]]@result[,c("Description","p.adjust"),drop=F] %>% setNames(c("rn","padj"))
    #} else if(class(comb_res[[ds]]) == "list") {
      #comb_res[[ds]]$combined = comb_res[[ds]]$combined[ comb_res[[ds]]$combined$Description != "\n", ]
      if(nrow(comb_res[[ds]]$combined) == 0) {
        stat_tmp = data.frame(rn = "tmp", NES = NA)
        pval_tmp = data.frame(rn = "tmp", pval = NA)
        padj_tmp = data.frame(rn = "tmp", padj = NA)
      } else {
        stat_tmp = comb_res[[ds]]$combined[,c("Description","NES"),drop=F] %>% setNames(c("rn","NES"))
        pval_tmp = comb_res[[ds]]$combined[,c("Description","pvalue"),drop=F] %>% setNames(c("rn","pval"))
        padj_tmp = comb_res[[ds]]$combined[,c("Description","p.adj"),drop=F] %>% setNames(c("rn","padj"))
      }
    #} else {
      #stat_tmp = comb_res[[ds]]$set.statistics %>% as.data.frame() %>% rownames_to_column("rn")
      #pval_tmp = comb_res[[ds]]$pval.adj %>% as.data.frame() %>% setNames("padj") %>% rownames_to_column("rn")
    #  stat_tmp = comb_res[[ds]][,c("Description","NES"),drop=F] %>% setNames(c("rn","NES"))
    #  pval_tmp = comb_res[[ds]][,c("Description","p.adj"),drop=F] %>% setNames(c("rn","padj"))
    #}
    if(nrow(all_stat) == 0) {
      all_stat = stat_tmp
      all_pval = pval_tmp
      all_padj = padj_tmp
    } else {
      all_stat = full_join(all_stat, stat_tmp, by = "rn")
      all_pval = full_join(all_pval, pval_tmp, by = "rn")
      all_padj = full_join(all_padj, padj_tmp, by = "rn")
    }
    if(!is.null(pval.adj.cutoff)) {
      all_sig_path = unique(c(all_sig_path, padj_tmp$rn[ padj_tmp$padj <= pval.adj.cutoff ] ))
    }
    if(!is.null(topN_per_dataset)) {
      all_sig_top_path = unique(c(all_sig_top_path, padj_tmp %>% dplyr::filter(padj <= pval.adj.cutoff) %>% arrange(padj) %>% head(topN_per_dataset) %>% pull(rn) ))
    }
  }
  
  if(!is.null(topN_per_dataset)) {
    all_sig_path = all_sig_top_path
  }
  
  # Remove "NA" values for path that are introduced due to this condition: (if(nrow(comb_res[[ds]]$combined) == 0)) above.
  all_sig_path = all_sig_path[ ! is.na(all_sig_path) ]
  if(length(all_sig_path) == 0) {
    return(ggplot() + theme_void())
  }
  
  dup = all_stat$rn[ duplicated(all_stat$rn) ]
  cat("Found", length(dup), "duplicate pathway ids.\n")
  print(head(dup))
  all_stat = all_stat[!duplicated(all_stat$rn),]
  all_pval = all_pval[!duplicated(all_pval$rn),]
  all_padj = all_padj[!duplicated(all_padj$rn),]

  all_stat = all_stat %>% remove_rownames() %>% column_to_rownames("rn") %>% setNames(names(comb_res))
  all_pval = all_pval %>% remove_rownames() %>% column_to_rownames("rn") %>% setNames(names(comb_res))
  all_padj = all_padj %>% remove_rownames() %>% column_to_rownames("rn") %>% setNames(names(comb_res))
  
  pval_mat_rd = all_pval
  pval_mat_rd[is.na(pval_mat_rd)] = -1
  pval_mat_rd[pval_mat_rd > cauchy_thr ] = (cauchy_thr+1.0)/2
  pval_mat_rd[pval_mat_rd == -1] = NA
  pval_tan = tan((0.5- pval_mat_rd) * pi)
  pval_combine = apply(pval_tan,1,mean, na.rm = T)
  pval_combine = pcauchy(pval_combine, location = 0, scale = 1, lower.tail = F, log.p = F)
  pval_combine[is.na(pval_combine)] = 1.0
  all_pval$joint = pval_combine
  all_padj$joint = pval_combine %>% p.adjust("BH")
  all_stat$joint = 0

  all_stat = all_stat[all_sig_path,]
  all_pval = all_pval[all_sig_path,]
  all_padj = all_padj[all_sig_path,]
  
  if(!is.null(max_pathway))
    pathway_levels = all_padj %>% rownames_to_column("rn") %>% arrange(joint) %>% pull(rn)
  
  if(!is.null(topN_per_dataset))
    pathway_levels = all_sig_path


  all_stat = pivot_longer(all_stat %>% rownames_to_column("rn"), cols = -"rn") %>% setNames(c("pathway","dataset","stat"))
  all_pval = pivot_longer(all_pval %>% rownames_to_column("rn"), cols = -"rn") %>% setNames(c("pathway","dataset","pval"))
  all_padj = pivot_longer(all_padj %>% rownames_to_column("rn"), cols = -"rn") %>% setNames(c("pathway","dataset","padj"))
  all_stat$pval.adj = all_padj$padj
  
  
  all_stat$pathway = factor(as.vector(all_stat$pathway), levels = rev(pathway_levels))
  all_stat$dataset = factor(as.vector(all_stat$dataset), levels = c(names(comb_res), "joint"))
  if(!is.null(max_pathway))
    all_stat = all_stat %>% dplyr::filter(pathway %in% head(pathway_levels, max_pathway) )
  p = ggplot(all_stat, aes(dataset, pathway, size=-log10(pval.adj), color=as.character(sign(stat)))) + 
    geom_point(shape=1) + scale_color_manual(values = c("1"="red","-1"="blue")) + 
    geom_point(data = all_stat %>% dplyr::filter(pval.adj < pval.adj.cutoff)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.5)) +
    ggtitle(title)
  return(p)
}





# Calculate composite score - term score. Adapted from Leying's code.
# ds_DF must have samples as rows and features as columns
calc_term_score <- function(ds_DF, 
                                        projection_coef_DF,
                                        leading_edges,
                                        sep = "/",
                                        is_vector = FALSE,
                                        secondary_weights = NULL # The secondary weights should be a numerical vector with same order as the leading_edges vector. Used only when is_vector is TRUE.
                            ){
  #leading_edges can be either "/" (sep) seperated or as a vector of features
  if(!is_vector){
    leading_edges = strsplit(leading_edges, sep)[[1]]
    leading_edges = intersect(rownames(projection_coef_DF),leading_edges)
  }
  leading_edges = leading_edges[ leading_edges %in% colnames(ds_DF) ]
  if(length(leading_edges)>0){
    w =  abs(projection_coef_DF[leading_edges,1])
    if(is_vector & ! is.null(secondary_weights) & length(w) == length(secondary_weights)) {
      warning("Using the secondary weights.")
      w = w * secondary_weights
    }
    x = as.matrix(ds_DF[,leading_edges])
    z = x %*% matrix(w, nrow = length(w), ncol = 1)
    z = as.data.frame(z) %>% `rownames<-`(rownames(ds_DF))
  }else{
    z = NULL
  }

  return(z)
}


# factor: factor name, e.g. "Factor1"
# database: allows selecting a term explicitly from a selected database.
composite_feature_creation <- function(obj, enrichment_mod_list, factor, term, res_tbl_type = "combined", database = NULL, qval_cutoff = 0.1) {
    enrichment_results_mod_factor = enrichment_mod_list[[factor]]
    datasets = names(enrichment_results_mod_factor)
    res_return = list()
    for(ds in datasets) {
        res_tbl = enrichment_results_mod_factor[[ds]][[res_tbl_type]]
        idx = res_tbl$Description == term & res_tbl$p.adj < qval_cutoff
        if(!is.null(database)) {
            idx = idx & res_tbl$db == database
        }
        if(sum(idx) > 1) {
            print("Warning::: Found the selected term in multiple pathway databases. The first one is selected by default. To avoid this behavior, specify the pathway database of interest.")
        }
        leading_edges = res_tbl$core_enrichment[idx]
        if(length(leading_edges)!=0) {
            ds_DF = obj$dataset_list[[ds]]$group1
            projection_coef_DF = obj$weights[[ds]][,factor,drop=F]
            res_return[[ds]] = calc_term_score(t(ds_DF), projection_coef_DF, leading_edges, sep = "#")
        }
    }
    return(res_return)
}


# "features" can accept number or a vector of feature names
# row_annot_lbls must be in the same order as the feature names in "features"
# clm_annot_lbls must be a named vector with names corresponding to the sample ids.
make_topFeat_heatmap <- function(model, view, factor, features = 40,
                                 cluster_rows=F, cluster_cols=F,
                                 addJoelsGrpRug=F, default_colors = F,
                                 show_rownames = T, show_colnames = T, title = "", 
                                 draw_plot = T, use_raster = F, row_annot_lbls = NULL, 
                                 clm_annot_lbls = NULL, ...) {
  tmp_mat = get_heatmap_data(model,
    view = view,         # view of interest
    factor = factor,             # factor of interest
    features = features,          # number of features to plot (they are selected by weight)
    ...
  )
  
  if(default_colors) {
    mycol =  circlize::colorRamp2(c(min(tmp_mat), quantile(as.vector(tmp_mat), 0.05), 0, quantile(as.vector(tmp_mat), 0.95), max(tmp_mat)), rev(brewer.pal(11,"RdBu")[c(1,3,6,9,11)]))
  } else {
    mycol = circlize::colorRamp2(c(-2, -1, 1, 2, 8), rev(brewer.pal(11,"RdYlBu")[c(1,3,6,9,11)]))
  }
  
  my_row_annot = NA
  if(! is.null(row_annot_lbls)) {
    my_row_annot = data.frame(
      row_annot = row_annot_lbls,
      row.names = rownames(tmp_mat)
    )
    my_row_annot_color = list(
      row_annot = kelly(length(unique(row_annot_lbls))+1)[-1] %>% setNames(sort(unique(my_row_annot$row_annot)))
    )
  }
  
  my_clm_annot = NA
  if(! is.null(clm_annot_lbls)) {
    my_clm_annot = data.frame(
      clm_annot = clm_annot_lbls[colnames(tmp_mat)],
      row.names = colnames(tmp_mat)
    )
    my_col_annot_color = list(
      clm_annot = kelly(length(unique(clm_annot_lbls))+1)[-1] %>% setNames(sort(unique(my_clm_annot$clm_annot)))
    )
  }
  
  
  if(addJoelsGrpRug) {
    my_clm_annot = data.frame(
      joelsGrp = case_when(
        colnames(tmp_mat) %in% lachno_donor ~ "Lachnospiraceae_rich",
        colnames(tmp_mat) %in% eggr_donor ~ "Eggertherllaceae_rich",
        TRUE ~ "other"
      ),
      row.names = colnames(tmp_mat)
    )
    my_col_annot_color = list(
      joelsGrp = kelly(6)[c(3,6,4)] %>% setNames(sort(unique(my_clm_annot$joelsGrp)))
    )
  }
  
  if(exists("my_row_annot_color") & exists("my_col_annot_color")) {
    my_annot_color = append(my_row_annot_color, my_col_annot_color)
  } else if(exists("my_row_annot_color")) {
    my_annot_color = my_row_annot_color
  } else if(exists("my_col_annot_color")) {
    my_annot_color = my_col_annot_color
  } else {
    my_annot_color = NA
  }

  #pdf("stoolSpecies_11.pdf", width=14, height=8)
  ComplexHeatmap::pheatmap(tmp_mat, color = mycol, cluster_cols = cluster_cols, cluster_rows = cluster_rows, annotation_col = my_clm_annot, annotation_row = my_row_annot, annotation_colors = my_annot_color, show_rownames = show_rownames, show_colnames = show_colnames, main = title, run_draw = draw_plot, use_raster = use_raster)
  #dev.off()
}




colMedians = function(mat, na.rm = F) {
    apply(mat, 2, median, na.rm=na.rm)
}




sc_celltypes <- c(
"B memory non-switched",
"B memory switched",
"B naive",
"Basophil",
"CD4 A4B7 memory",
"CD4 A6 memory",
"CD4 cytotoxic",
"CD4 Naive A",
"CD4 Naive B",
"CD4 Naive C",
"CD4 Naive D",
"CD4 Naive E",
#"CD4 TEMRA",
"CD4 Th1-like/Tem",
"CD4 Th2-like/Tcm",
"CD8 A1B1/A6B1 Tem",
"CD8 AEB7 memory",
"CD8 cytotoxic",
"CD8 Naive A",
"CD8 Naive B",
"CD8 Naive C",
"CD8 Naive D",
"CD8 Tem",
"cDC2",
"cMo A",
"cMo B",
"cMo C",
"cMo_T",
"gdT",
"iNK",
"MAIT DNAM-1low",
"mNK A",
"mNK B",
"NK adaptive memory",
"ncMo",
"pDC",
"Treg CD25high"
#"Unlabelled"
)

# The left_mat and right_mat must have matching rows.
my_cor_func <- function(left_mat, right_mat, method="spearman", padj_method="BH") {
  colnames(left_mat) = paste0("left_", colnames(left_mat))
  colnames(right_mat) = paste0("right_", colnames(right_mat))
  concat_cors = rcorr(as.matrix(left_mat), as.matrix(right_mat), type = method)
  concat_cors = flattenCorrMatrix(concat_cors$r, concat_cors$P)
  concat_cors = concat_cors[ grepl("left_", concat_cors$row) != grepl("left_", concat_cors$column) | grepl("right_", concat_cors$row) != grepl("right_", concat_cors$column), ]
  concat_cors$row = str_remove(concat_cors$row, "left_|right_")
  concat_cors$column = str_remove(concat_cors$column, "left_|right_")
  concat_cors$p_adj = p.adjust(concat_cors$p, method = padj_method)
  return(concat_cors)
}




# Needs the "obj" loaded in memory
plot_gex_weights <- function(features, factor, title) {
  my_weights = data.frame()
  my_signif = data.frame()
  for(ctype in sc_celltypes) {
      my_weights = rbind(my_weights, as.data.frame(t(obj$weights.noOutl[[ctype]][match(features, rownames(obj$weights.noOutl[[ctype]])), factor])) %>% setNames(features) )
      signif_lbls = rep("", length(features))
      signif_lbls[ features %in% names(obj$signif_features[[factor]][[ctype]]) ] = "*"
      signif_lbls[ features %in% names(obj$signif_features_afterCorr[[factor]][[ctype]]) ] = "**"
      my_signif = rbind(my_signif, data.frame(
          data.frame(t(signif_lbls))
      ) )
  }
  my_weights[is.na(my_weights)] = 0
  rownames(my_weights) = sc_celltypes
  rownames(my_signif) = sc_celltypes
  colnames(my_signif) = features
  ComplexHeatmap::Heatmap(my_weights, clustering_method_rows = "ward.D2", clustering_method_columns = "ward.D2", cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sprintf("%s", my_signif[i, j]), x, y, gp = gpar(fontsize = 10))
  }, column_title = title)
}




compare_LE_between_factors <- function(enrichment_results_mod, term, dedup=FALSE) {
    tmp_df = data.frame()
    all_f1s = c()
    all_f1n3s = c()
    all_f3s = c()
    for(ctype in sc_celltypes) {
        res1 = enrichment_results_mod$Factor1[[ctype]]$combined[enrichment_results_mod$Factor1[[ctype]]$combined$Description == term,]
        res3 = enrichment_results_mod$Factor3[[ctype]]$combined[enrichment_results_mod$Factor3[[ctype]]$combined$Description == term,]
        if(nrow(res1) == 0) {
            res1 = data.frame(qvalues = 1, core_enrichment = "", NES=0)
        }
        if(nrow(res3) == 0) {
            res3 = data.frame(qvalues = 1, core_enrichment = "", NES=0)
        }
        if(res1$qvalues < 0.1 || res3$qvalues < 0.1) {
            if(res1$qvalues >= 0.1) {
                res1$core_enrichment = ""
                res1$NES = 0
            }
            if(res3$qvalues >= 0.1) {
                res3$core_enrichment = ""
                res3$NES = 0
            }
            v1 = unlist(str_split(res1$core_enrichment, "#"))
            v3 = unlist(str_split(res3$core_enrichment, "#"))
            all_f1s = unique(c(all_f1s, setdiff(v1,v3)))
            all_f1n3s = unique(c(all_f1n3s, intersect(v1,v3)))
            all_f3s = unique(c(all_f3s, setdiff(v3,v1)))
            tmp_df = rbind(
                tmp_df,
                data.frame(
                    f1 = paste0(setdiff(v1,v3), collapse = ", "),
                    f1n3 = paste0(intersect(v1,v3), collapse = ", "),
                    f3 = paste0(setdiff(v3,v1), collapse=", "),
                    dir_f1 = res1$NES,
                    dir_f3 = res3$NES,
                    row.names = ctype
                )
            )
        }
    }
    if(dedup) {
        print("dedup=TRUE is now deprecated.")
        if(FALSE) {
            for(i in 1:nrow(tmp_df)) {
                f1 = unlist(str_split(tmp_df[i, "f1"], ", "))
                f1n3 = unlist(str_split(tmp_df[i, "f1n3"], ", "))
                f3 = unlist(str_split(tmp_df[i, "f3"], ", "))
                tmp_df[i, c("f1","f1n3","f3")] = c(
                    paste0(f1[! f1 %in% c(all_f3s)], collapse = ", "),
                    paste0(f1n3, collapse = ", "),
                    paste0(f3[! f3 %in% c(all_f1s)], collapse = ", ")
                )
            }
        }
    }
    f1_clm_union = grep("^$",sort(unique(str_split(paste(tmp_df$f1, collapse = ", "), ",\\s*")[[1]])), value = T, invert = T)
    f1n3_clm_union = grep("^$",sort(unique(str_split(paste(tmp_df$f1n3, collapse = ", "), ",\\s*")[[1]])), value = T, invert = T)
    f3_clm_union = grep("^$",sort(unique(str_split(paste(tmp_df$f3, collapse = ", "), ",\\s*")[[1]])), value = T, invert = T)
    tmp_df = rbind(
        tmp_df,
        data.frame(
            f1 = paste(f1_clm_union, collapse=", "),
            f1n3 = paste(f1n3_clm_union, collapse=", "),
            f3 = paste(f3_clm_union, collapse=", "),
            dir_f1 = NA,
            dir_f3 = NA,
            row.names = "All celltypes combined - union of above rows"
        )
    )
    tmp_df = rbind(
        tmp_df,
        data.frame(
            f1 = paste(f1_clm_union[! f1_clm_union %in% c(f1n3_clm_union, f3_clm_union)], collapse=", "),
            f1n3 = paste(intersect( unique(c(f1_clm_union, f1n3_clm_union)), unique(c(f3_clm_union, f1n3_clm_union)) ), collapse=", "),
            f3 = paste(f3_clm_union[! f3_clm_union %in% c(f1n3_clm_union, f1_clm_union)], collapse=", "),
            dir_f1 = NA,
            dir_f3 = NA,
            row.names = "Across all celltypes"
        )
    )
    return(tmp_df) 
}


sc_parents1 = c("B","B","B","Basophil","CD4","CD4","CD4","CD4","CD4","CD4","CD4","CD4","CD4","CD4","CD8","CD8","CD8","CD8","CD8","CD8","CD8","CD8","DC","Mo","Mo","Mo","Mo","gdT","NK","MAIT","NK","NK","NK","Mo","DC","Treg") %>% setNames(sc_celltypes)
sc_parents2 = c("B","B","B","Basophil","CD4 mem","CD4 mem","CD4 cyto","CD4 naive","CD4 naive","CD4 naive","CD4 naive","CD4 naive","CD4 mem","CD4 mem","CD8 mem","CD8 mem","CD8 cyto","CD8 naive","CD8 naive","CD8 naive","CD8 naive","CD8 mem","cDC","cMo","cMo","cMo","cMo","gdT","NK","MAIT","NK","NK","NK","ncMo","pDC","Treg") %>% setNames(sc_celltypes)
# Broad populations for calculating 'effective expression' across the sub-celltypes
sc_parents3 = c("B non-naive", "B non-naive", "B naive", "Basophil", "CD4 non-naive", "CD4 non-naive", "CD4 non-naive", "CD4 naive", "CD4 naive", "CD4 naive", "CD4 naive", "CD4 naive", "CD4 non-naive", "CD4 non-naive", "CD8 non-naive", "CD8 non-naive", "CD8 non-naive", "CD8 naive", "CD8 naive", "CD8 naive", "CD8 naive", "CD8 non-naive", "myeloid", "myeloid", "myeloid", "myeloid", "cMo_T", "gdT", "NK", "MAIT DNAM-1low", "NK", "NK", "NK", "myeloid", "pDC", "Treg CD25high") %>% setNames(sc_celltypes)



