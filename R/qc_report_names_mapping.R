#' QC plots names
#' @description QC report names of all the plots generated with the MassExpression workflow.

get_names_qc <- function(){
  return(c(
    "DE_summary_proteins",
    "PCA_proteins",
    "PCA_screeplot_proteins",
    "PCA_DE_proteins",
    "CV_proteins",
    "CV_table_proteins",
    "samples_correlations_proteins",
    "samples_correlations_DE_proteins",
    "missing_by_proteins",
    "missing_by_sample_proteins",
    "intensities_boxplots_proteins",
    "RLE_raw_proteins",
    "RLE_normalised_proteins",
    "imputed_proteins"
  ))
}
