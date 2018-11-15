#' Filter MT and Ribsomal protein genes
#'
#' Removes all mitochondrial and ribosomal protein genes in the count matrix and returns the reduces matrix.
#'
#' @return remove_ribosomal_proteins_and_mitochondrial_genes_from_matrix returns a reduced matrix with all mitochondrial and ribosomal protein genes removed
#'

#' @export
remove_ribosomal_proteins_and_mitochondrial_genes_from_matrix <- function(cell_by_ensembl_mat = matrix(), species = 'hsapiens', rm_cell_cycle_genes = FALSE)  {
  library(biomaRt)
  if(!all(grepl("^ENS", colnames(cell_by_ensembl_mat))))  {
    stop("Column names should be ensembl ids and should start with 'ENS'")
  }
  if(species == 'hsapiens')  {
    ensembl <- biomaRt::useMart("ensembl", dataset="hsapiens_gene_ensembl")
    useDataset("hsapiens_gene_ensembl", mart=ensembl)
    ensembl_gene_metadata_df <- biomaRt::select(ensembl, keys = colnames(cell_by_ensembl_mat), keytype = 'ensembl_gene_id', columns=c('ensembl_gene_id','chromosome_name', "hgnc_symbol", "description")) %>% dplyr::rename(ENSEMBL = ensembl_gene_id)
    ensembl_gene_metadata_MT_and_ribosomal_df <- dplyr::filter(ensembl_gene_metadata_df, grepl("^M?RP[SL]", hgnc_symbol) | (chromosome_name == 'MT'))
  }
  else if(species == 'mmusculus') {
    ensembl <- biomaRt::useMart("ensembl", dataset="mmusculus_gene_ensembl")
    useDataset("mmusculus_gene_ensembl", mart=ensembl)
    ensembl_gene_metadata_df <- biomaRt::select(ensembl, keys = colnames(cell_by_ensembl_mat), keytype = 'ensembl_gene_id', columns=c('ensembl_gene_id','chromosome_name', "mgi_symbol", "description")) %>% dplyr::rename(ENSEMBL = ensembl_gene_id)
    ensembl_gene_metadata_MT_and_ribosomal_df <- dplyr::filter(ensembl_gene_metadata_df, grepl("^M?RP[SL]", mgi_symbol, ignore.case = TRUE) | (chromosome_name == 'MT'))
  }
  else {
    warning("This is not a valid species option")
    stop()
  }
  #At this point we know we have either hsapiens or mmusculus species, so we go ahead and determine whether to remove the cell cycle genes
  cell_cycle_genes <- character()
  if(rm_cell_cycle_genes)  {
    cell_cycle_df <- read.table(dir(file.path(system.file(package = "rnaseqUtils"), 'extdata'), paste0(species, "_cell_cycle_genes.tsv"), full.names = TRUE), sep = '\t', header = TRUE, stringsAsFactors = FALSE)
    #current_gene_count <- ncol(ensembl_gene_metadata_MT_and_ribosomal_df)
    cell_cycle_genes <- cell_cycle_df$ENSEMBL
  }
  print(head(ensembl_gene_metadata_MT_and_ribosomal_df))
  ensembl_overlaps <- colnames(cell_by_ensembl_mat)[which(colnames(cell_by_ensembl_mat) %in% unique(c(ensembl_gene_metadata_MT_and_ribosomal_df$ENSEMBL, cell_cycle_genes)))]
  message(paste(length(ensembl_overlaps), 'of the', ncol(cell_by_ensembl_mat), 'ensembl ids are being removed'))
  return(cell_by_ensembl_mat[,!colnames(cell_by_ensembl_mat) %in% ensembl_overlaps])
}


#' Filter random genes from matrix
#'
#' Remove genes from raw count matrix that are unlikely contributing to the biology based on their expression patterns
#'
#' This function performs the following filtering tasks. It initially removes all columns (genes) that have less than two non-zero entries, ie. that are expressed in zero or one cell.
#' This c

#' @export
get_random_genes_from_matrix <-function(count_matrix = matrix(), recompute_dist_removal_threshold = 500, min_cell_count_threshold = 4, distance_matrix = NULL, distance_function = dist, normalize_matrix = FALSE, heuristic_threshold = 200, heuristic_step_size = 1000, n_sample_draws = 1000, tail_area = 0.001)  {
  #remove zero counts
  sample_names <- rownames(count_matrix)
  gene_counts <- apply(count_matrix, 2, function(.col)  {sum(.col > 0)})
  message(paste("Removing", sum(gene_counts < min_cell_count_threshold), ' genes expressed in less than', min_cell_count_threshold, 'cells from initial matrix of', nrow(count_matrix), 'by', ncol(count_matrix), 'matrix, leaving',  sum(gene_counts >= min_cell_count_threshold), 'genes'))
  min_thresh_mat <- count_matrix[,gene_counts >= min_cell_count_threshold]
  min_thresh_counts <- apply(min_thresh_mat, 2, function(.col)  {sum(.col > 0)})
  min_thresh_mat_keep_logical <- rep(TRUE, ncol(min_thresh_mat))
  names(min_thresh_mat_keep_logical) <- colnames(min_thresh_mat)
  
  if(normalize_matrix)  {
    min_thresh_mat <- t(apply(min_thresh_mat, 1, function(.x) {log2(.x/sum(.x) * 1e5 + 1)}))
  }
  if(is.null(distance_matrix))  {
    dist_mat <- distance_function(min_thresh_mat)
  }
  else {
    dist_mat <- distance_matrix
  }
  for(i in 1:1000)  {
    restart = FALSE
    print(paste("Iteration", i))
    #Create a list of vectors, with each vector containing cell counts corersponding to genes that will share the same sample of randomly selected cell distances
    valid_counts_less_than_threshold <- sort(unique(min_thresh_counts[min_thresh_mat_keep_logical][(min_thresh_counts[min_thresh_mat_keep_logical] < heuristic_threshold) & (min_thresh_counts[min_thresh_mat_keep_logical] >= min_cell_count_threshold)]))
    #print(valid_counts_less_than_threshold)
    list_of_cell_counts_to_test <- lapply(c(valid_counts_less_than_threshold, seq(heuristic_threshold, max(min_thresh_counts[min_thresh_mat_keep_logical]), by = heuristic_step_size)), function(.n)  {
      if(.n < heuristic_threshold)  {
        return(unique((min_thresh_counts[min_thresh_mat_keep_logical])[min_thresh_counts[min_thresh_mat_keep_logical] == .n]))
      }
      else {
        return(unique(min_thresh_counts[min_thresh_mat_keep_logical][(min_thresh_counts[min_thresh_mat_keep_logical]  >= .n) & (min_thresh_counts[min_thresh_mat_keep_logical] < (.n + heuristic_step_size))]))
      }
    })
    n_removed_since_last_dist_update <- 0
    pvals_and_dif_mats_list <- lapply(1:length(list_of_cell_counts_to_test), function(.x)  {return(NULL)})
    for(i in 1:length(list_of_cell_counts_to_test))  {
      counts <- list_of_cell_counts_to_test[[i]]
      min_count <- min(counts)
      sampled_mean_dists <- rnaseqUtils:::get_sampled_distances(dist_mat, n_samples = n_sample_draws, n_cells_per_sample = min_count)
      #sampled_mean_mean_dist <- mean(sampled_mean_dists)
      pval_and_dif_by_gene_count_mat <- sapply(names(min_thresh_counts[min_thresh_mat_keep_logical][min_thresh_counts[min_thresh_mat_keep_logical] %in% counts]), function(.name)  {
        expr_values_for_gene <- min_thresh_mat[,.name]
        logical_of_expressed_genes <- expr_values_for_gene > 0
        n_expr_genes <- sum(logical_of_expressed_genes)
        if(n_expr_genes == min_count)  {
          gene_mean_dist <- rnaseqUtils:::get_mean_dist_from_matrix(dist_mat[logical_of_expressed_genes, logical_of_expressed_genes])
        }
        else {
          gene_mean_dist <- mean(sapply(1:5, function(.iter)  {
            subsampled_inds <- sort(sample(which(logical_of_expressed_genes), min_count, replace = FALSE))
            rnaseqUtils:::get_mean_dist_from_matrix(dist_mat[subsampled_inds, subsampled_inds])
          }))
        }
        return(c(mean(sampled_mean_dists < gene_mean_dist), log10(gene_mean_dist/min(sampled_mean_dists))))
      })
      logical_failing_genes <- pval_and_dif_by_gene_count_mat[1,] > tail_area
      updated_min_thresh_mat_keep_logical <- min_thresh_mat_keep_logical
      updated_min_thresh_mat_keep_logical[names(which(logical_failing_genes))] <- FALSE
      min_thresh_mat_keep_logical <- updated_min_thresh_mat_keep_logical
      
      message(paste("Removing", sum(logical_failing_genes), "of", ncol(pval_and_dif_by_gene_count_mat), "genes that are expressed in", min_count, 'to', max(counts), 'cells'))
      updated_n_removed_since_last_dist_update <- n_removed_since_last_dist_update + sum(logical_failing_genes)
      n_removed_since_last_dist_update <- updated_n_removed_since_last_dist_update
      #print(paste("Nowa at", n_removed_since_last_dist_update,'removed of', recompute_dist_removal_threshold,'before recomputing'))
      if(n_removed_since_last_dist_update >= recompute_dist_removal_threshold)  {
        n_removed_since_last_dist_update <- 0
        message(paste("Over", recompute_dist_removal_threshold, "genes_have been removed since last distance matrix update, triggering new distance matrix computation of", sum(min_thresh_mat_keep_logical), "genes" ))
        min_thresh_mat_keep_logical[colnames(min_thresh_mat) %in% colnames(pval_and_dif_by_gene_count_mat[,logical_failing_genes])] <- FALSE
        dist_mat <- distance_function(min_thresh_mat[,min_thresh_mat_keep_logical])
        restart = TRUE
        break
      }
      pvals_and_dif_mats_list[[i]] <- pval_and_dif_by_gene_count_mat
    }
    if(restart)  {
      next
    }
    else{
      message(paste("Final filtered_count matrix has", ncol(min_thresh_mat), 'after removal of', ncol(count_matrix) - ncol(min_thresh_mat), ' genes'))
    }
    return(list(min_thresh_mat, dist_mat, do.call(cbind, pvals_and_dif_mats_list)))
  }
}

#` Calculate mean value from matrix
#'
#' The function calculates and returns the mean of the upper-triangular portion of a matrix
#'
#' @return get_mean_dist_from_matrix returns the upper-triangular mean of temp_mat, exluding the zeros from the diagonal or the lower-triangular portion of the matrix
#'
get_mean_dist_from_matrix <- function(temp_mat)  {
  temp_mat[!upper.tri(temp_mat)] <- 0
  mean(temp_mat[temp_mat != 0])
}

get_sampled_distances <- function(dist_mat, n_samples, n_cells_per_sample)  {
  mean_pw_cor <- sapply(1:n_samples, function(.ind)  {
    sub_sample_inds <- sort(sample(1:ncol(dist_mat), n_cells_per_sample, replace = FALSE))
    sampled_dist_mat <- dist_mat[sub_sample_inds,sub_sample_inds]
    get_mean_dist_from_matrix(sampled_dist_mat)
  })
}
