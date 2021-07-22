#' Filter protein IDs that are not present in at least one experiment

#' @param funDT data.table
#' @param condition_col_name name of column indicating the condition of interest
#' @param run_id_col_name name of column indicating the protein id
#' 
#' @import data.table
#' @return `list`: filtered data.table and list of proteins 
#' @export prepareForFiltering

prepareForFiltering <- function(funDT, condition_col_name, run_id_col_name){
  filterDT <- unique(funDT[, .(get(condition_col_name), get(run_id_col_name))])[, .(max_count = .N), by = .(V1)]
  filterDT <- merge(funDT, filterDT, by.x = condition_col_name, by.y = "V1", all.x = T)
  filterDT <- filterDT[Imputed == 0, .(count_rep = .N, 
                                       max_count = max(max_count, na.rm = T)), 
                       by = c(condition_col_name, "ID")][, repPC := count_rep/max_count]
  
  return(filterDT)
}