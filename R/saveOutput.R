
#' Save RData with SummarizedExperiment object and ProteinViz json
#' 
#' @param IntensityExperiment XX.
#' @param CompleteIntensityExperiment XX
#' @param longIntensityDT XX
#' @param output_folder XX
#' @import data.table
#' @export

saveOutput <- function(IntensityExperiment, CompleteIntensityExperiment, 
                       longIntensityDT, 
                       output_folder){
  dir.create(file.path(output_folder), recursive = TRUE, showWarnings = FALSE)
  
  comparisonExperiments <- 
    listComparisonExperiments(CompleteIntensityExperiment)
  
  save(IntensityExperiment, 
       CompleteIntensityExperiment, 
       comparisonExperiments,
       file = file.path(output_folder, "DiscoveryQuant.RData"))
  
  writeProteinViz(CompleteIntensityExperiment, output_folder)
  
  # Write Replicates tab data
  limmaStats <- SummarizedExperiment::rowData(CompleteIntensityExperiment)[, c("ProteinId", "GeneName", "Description")]
  longDTProt <- merge(longIntensityDT, limmaStats, by="ProteinId", all.x=TRUE)
  longDTProt <- as.data.table(longDTProt)
  
  writeReplicateData(longDTProt, output_folder)
}


#' Write limma statistics for users to read
#' @param CompleteIntensityExperiment produced by runLimmaPipeline
#' @param outputFolder A path to a folder to write the statistics to.
#' @export writeLimmaStatisticsTable

writeLimmaStatisticsTable <- function(CompleteIntensityExperiment, outputFolder){
  write.table(rowData(CompleteIntensityExperiment), sep = "\t",
              file = file.path(outputFolder, "protein_limma_statistics.tsv"))
}


#' Write data needed for the Replicates tab
#' This is a wrapper of `oneProteinReplData`, looping over all proteins.
#' 
#' @param longDTProt data.table protein information stored in long format. 
#' Columns required: `ProteinId`, `GeneName`, `Description`, `log2NInt`, `Condition`,
#'  `Replicate`, `Imputed`. 
#' @param outputFolder A path to a folder to write the statistics to.
#' @importFrom jsonlite write_json
#' @export writeReplicateData

writeReplicateData <- function(longDTProt, outputFolder){
  proteinSet <- unique(longDTProt$ProteinId)
  protList <- lapply(proteinSet, function(prot) oneProteinReplData(longDTProt[ProteinId %in% prot,]))
  protDF <- do.call(rbind, protList)
  
  dir.create(outputFolder, showWarnings = FALSE)
  outPath = file.path(outputFolder,"protein_counts_and_intensity.json")
  write_json(protDF, outPath, digits = NA, na = "null")
}


#' Transform data for one protein from a long format to the nested data structure needed for the Replicates tab.
#' @param oneProt data.table single protein information stored in long format. 
#' Columns required: `ProteinId`, `GeneName`, `Description`, `log2NInt`, `Condition`,
#'  `Replicate`, `Imputed`. 
#' @export oneProteinReplData
#' 
oneProteinReplData <- function(oneProt){
  infoProt <- unique(oneProt[,c("ProteinId", "GeneName", "Description")])
  
  infoConds <- oneProt[, numberOfReplicateCount:= length(Replicate), by = Condition ]
  infoConds <- infoConds[, precentageOfReplicates:= sum(Imputed==0)/length(Replicate), by = Condition ]
  infoConds <- unique(infoConds[, c("Condition", "numberOfReplicateCount", "precentageOfReplicates")])
  setnames(infoConds, old = "Condition", new = "name")
  
  conditions <- data.frame(matrix(NA, nrow = length(infoConds$name), ncol = 4))
  for(cond_idx in 1:length(infoConds$name)){
    cond <- infoConds$name[cond_idx]
    infoOneCond <- infoConds[name %in% cond, ]
    
    oneCondRepl <- data.table(oneProt)[Condition %in% cond, c("log2NInt", "Imputed")] 
    oneCondRepl$replicateNum <- 1:nrow(oneCondRepl)
    
    entryCond <- dplyr::tibble(infoOneCond, intensityValues=list(oneCondRepl))
    conditions[cond_idx, ] <- entryCond
  }
  
  colnames(conditions) <- c("name", "numberOfReplicateCount", "precentageOfReplicates", "intensityValues")
  # Combine with protein Infos
  oneProtNested <- dplyr::tibble(infoProt, conditions=list(conditions))
}