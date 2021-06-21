
#' This function does MNAR imputation
#'
#' @param myQuantDT quantification data
#' @param id_type the id column of myQuantDT that indicates
#' @param int_type the column of myQuantDT to be imputed
#' @param f_imputStDev The Standard Deviation parameter for MNAR Imputation
#' @param f_imputePosition The Position parameter for MNAR Imputation
#' @return quantification data with missing values imputed
#' @examples
#'  ## dt_int <- imputeLFQ(myQuantDT = dt_int,
#'  ## id_type = "id",
#'  ## int_type = "log2NInt", #log2Intensity
#'  ## 0.3,
#'  ## 1.8,
#'  ## )
#' @export
#' @importFrom data.table setnames

imputeLFQ <- function(myQuantDT,
                       id_type,
                       int_type,
                       f_imputeStDev,
                       f_imputePosition) {
  
  myQuantDT[, int_impute := get(int_type)]
  myQuantDT <- myQuantDT[, colnames(myQuantDT) != int_type, with = F]
  impute.StdDev <-
    f_imputeStDev * myQuantDT[Imputed == 0, sd(int_impute)]
  impute.position <-
    myQuantDT[Imputed == 0, mean(int_impute)] - f_imputePosition * myQuantDT[Imputed == 0, sd(int_impute)]
  set.seed(255)
  myQuantDT[Imputed == 1, int_impute := rnorm(.N, (impute.position), (impute.StdDev))]
  
  # Add columns for centered and centered + scaled intensities
  myQuantDT[, nRLE := scale(int_impute, center = TRUE, scale = FALSE), by = .(get(id_type))]
  myQuantDT[, z_norm := scale(int_impute, center = TRUE, scale = TRUE), by = .(get(id_type))]
  
  setnames(myQuantDT, "int_impute", int_type)
  
  return(myQuantDT)
}