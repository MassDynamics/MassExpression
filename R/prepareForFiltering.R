#' Filter protein IDs that are not present in at least one experiment

#' @param longIntensityDT data.table
#' @param conditionColname name of column indicating the condition of interest
#' @param runIdColname name of column indicating the protein id
#' 
#' @import data.table
#' @return `list`: filtered data.table and list of proteins 
#' @export prepareForFiltering

#TODO: when there will be another condition I can stratify by the other condition too

prepareForFiltering <- function(longIntensityDT, conditionColname, runIdColname){
  # how many replicates in each conditions
  filterDT <- unique(longIntensityDT[, .(get(conditionColname), get(runIdColname))])[, .(max_count = .N), by = .(V1)]
  filterDT <- merge(longIntensityDT, filterDT, by.x = conditionColname, by.y = "V1", all.x = T)
  filterDT <- filterDT[Imputed == 0, .(count_rep = .N, 
                                       max_count = max(max_count, na.rm = T)), 
                       by = c(conditionColname, "ID")][, repPC := count_rep/max_count]
  
  return(filterDT)
}
