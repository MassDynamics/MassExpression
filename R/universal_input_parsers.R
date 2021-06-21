# ui <- list(intensity_colnames, colname_runs, colname_condition, colname_protein_id)

create_se <- function(assay,
                      sample_info,
                      metadata,
                      colname_runs,
                      feature_id){

  # Check that names of runs in sample_info are available as columns of essays
  return_overlap_assay_sample_info <-

  se <- SummarizedExperiment(assays = assay,
                             colData = sample_info,
                             metadata = metadata)
}


return_overlap_assay_sample_info <- function(assay, sample_info, colname_runs){
  # assay should have as many columns as sample_info rows
  # There should be 1:1 between colnames in assay and rows in sample_info

  run_names_from_essay <- colnames(assay)
  run_names_from_sample_info <- unique(sample_info[, colname_runs])
  ncol_assay <- ncol(assay)
  nrow_sample_info <- nrow(sample_info)

  intersect_names <- intersect(run_names_from_essay, run_names_from_sample_info)

  if(length(intersect_names) == 0){
    warning("Runs in assay cannot be matched to runs in the experimental design.
         Check that each column name in assay is present in the experimental
         design.")
    assay <- NULL
    sample_info <- NULL
  }

  if(ncol_assay != nrow_sample_info){
    warning("The number of runs available in essay do not match the number
            of runs available in the experimental design.
            Only overlapping runs are going to be considered.")
  }

  assay <- assay[, intersect_names]
  sample_info <- sample_info %>% filter(get(colname_runs) %in% intersect_names)

  assay <- assay[,sort(colnames(assay))]
  sample_info <- sample_info %>% arrange(colname_runs)

  list(assay=assay, sample_info=sample_info)

}
