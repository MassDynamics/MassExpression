
#' Save RData with SummarizedExperiment object and ProteinViz json
#' 
#' @param IntensityExperiment Summarised experiment object returned by the `runGenericDiscovery` function. 
#' This contains the initial data (row intensities and expeirment design) provided in input.
#' @param CompleteIntensityExperiment Summarised experiment object returned by the `runGenericDiscovery` function. 
#' This contains the results of the differential expression analysis, the intensities used for the analysis and a masking logical matrix 
#' showing which intensities have been imputed.   
#' @param longIntensityDT Tabular data structure resturned by the `runGenericDiscovery` function. 
#' This contains, for each protein and sample, the raw, log2-transformed, normalised (when required) intensities, 
#' and a column indicating whether imputation has occurred for each intensity.   
#' @param outputFolder str. Path to folder where `data` should be saved
#' @import data.table
#' @export

saveOutput <- function(IntensityExperiment, CompleteIntensityExperiment, 
                       longIntensityDT, 
                       outputFolder){
  dir.create(file.path(outputFolder), recursive = TRUE, showWarnings = FALSE)
  
  comparisonExperiments <- 
    listComparisonExperiments(CompleteIntensityExperiment)
  
  save(IntensityExperiment, 
       CompleteIntensityExperiment, 
       comparisonExperiments,
       file = file.path(outputFolder, "DiscoveryQuant.RData"))
  
  # Write protein Viz
  writeProteinViz(CompleteIntensityExperiment, outputFolder)
  
  # Write Replicates tab data
  limmaStats <- SummarizedExperiment::rowData(CompleteIntensityExperiment)[, c("ProteinId", "GeneName", "Description")]
  longDTProt <- merge(longIntensityDT, limmaStats, by="ProteinId", all.x=TRUE)
  longDTProt <- as.data.table(longDTProt)
  
  writeReplicateData(longDTProt, outputFolder)
  
  # Write results files
  writeProteinQuantIntensities(CompleteIntensityExperiment, outputFolder)
  writeOutputFile(longIntensityDT, outputFolder, "protein_quant_intensities.txt")
  
  
}


#' Generic function to write `data` to `outputFolder` with `fileName`
#' 
#' The data are saved in tab delimited format.
#' 
#' @param data tabular data
#' @param outputFolder str. Path to folder where `data` should be saved.
#' @param fileName str. Name of file to write in output.
#' @import data.table
#' @import parallel
#' @export
writeOutputFile <- function(data, outputFolder, fileName){
  dir.create(outputFolder, showWarnings = FALSE)
  fwrite(data, file.path(outputFolder, fileName), row.names = FALSE,
         sep = "\t", quote = FALSE,
         showProgress = T, verbose = T, nThread = parallel::detectCores() - 1)
  
}

#' Save txt file with differential expression analysis results and the intensities used for the analysis.   
#' 
#' @param CompleteIntensityExperiment Summarised experiment object returned by the `runGenericDiscovery` function. 
#' This contains the results of the differential expression analysis, the intensities used for the analysis and a masking logical matrix 
#' showing which intensities have been imputed.   
#' @param outputFolder str. Path to folder where `data` should be saved.
#' @import data.table
#' @export
writeProteinQuantIntensities <- function(CompleteIntensityExperiment, outputFolder){
  limmaStatsS4 <- SummarizedExperiment::rowData(CompleteIntensityExperiment)
  limmaStats <- as.data.table(limmaStatsS4)
  colnames(limmaStats) <- colnames(limmaStatsS4)
  
  descriptionColumnsDT <- limmaStats[, c("ProteinId", "GeneName", "Description")]
  anovaCol <- ifelse("F" %in% colnames(limmaStats), "F", "t")
  resultsColumns <- c("AveExpr", anovaCol, "adj.P.Val", 
                      colnames(limmaStats)[grep("logFC |CI.L |CI.R |P.Value |adj.P.Val |NImputed|NReplicates",
                      colnames(limmaStats))])
  resultsDT <- limmaStats[, ..resultsColumns]
  
  # Intensities and imputed T/F
  intensities <- data.table(assay(CompleteIntensityExperiment))
  colnames(intensities) <- paste0("intensities_", colnames(intensities))
  
  imputed <- data.table(assays(CompleteIntensityExperiment)$imputedMask)
  colnames(imputed) <- paste0("imputed_", colnames(imputed))
  
  tot_samples <- ncol(intensities)
  order_columns <- do.call(c, lapply(1:tot_samples, function(x) c(x, x+tot_samples)))
  
  inten_imputed <- cbind(intensities, imputed)
  inten_imputed <- inten_imputed[,..order_columns]
  
  complete_results <- cbind(descriptionColumnsDT, inten_imputed, resultsDT)
  
  writeOutputFile(complete_results, outputFolder, "protein_quant.txt")
}


#' Write limma statistics for users to read
#' @param CompleteIntensityExperiment Summarised experiment object returned by the `runGenericDiscovery` function. 
#' This contains the results of the differential expression analysis, the intensities used for the analysis and a masking logical matrix 
#' showing which intensities have been imputed.   
#' @param outputFolder str. Path to folder where `data` should be saved.
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
#' @param outputFolder str. Path to folder where `data` should be saved.
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
  setnames(infoProt, old = "Description", new = "ProteinDescription")
  infoProt <- infoProt[, ProteinGroupId:=ProteinId]
  infoProt <- infoProt[,c("ProteinGroupId", "ProteinId", "GeneName", "ProteinDescription")]
  
  infoConds <- oneProt[, numberOfReplicateCount:= sum(Imputed==0), by = Condition ]
  infoConds <- infoConds[, precentageOfReplicates:= sum(Imputed==0)/length(Replicate), by = Condition ]
  infoConds <- unique(infoConds[, c("Condition", "precentageOfReplicates","numberOfReplicateCount")])
  setnames(infoConds, old = "Condition", new = "name")
  
  conditions <- data.frame(matrix(NA, nrow = length(infoConds$name), ncol = 4))
  for(cond_idx in 1:length(infoConds$name)){
    cond <- infoConds$name[cond_idx]
    infoOneCond <- infoConds[name %in% cond, ]
    
    oneCondRepl <- data.table(oneProt)[Condition %in% cond, c("SampleName","log2NIntNorm", "Imputed")] 
    setnames(oneCondRepl, old = "log2NIntNorm", new = "log2NInt_ProteinGroupId")
    setnames(oneCondRepl, old = "SampleName", new = "Replicate")
    oneCondRepl <- oneCondRepl[,c("Replicate", "log2NInt_ProteinGroupId", "Imputed")]
    
    entryCond <- dplyr::tibble(infoOneCond, intensityValues=list(oneCondRepl))
    conditions[cond_idx, ] <- entryCond
  }
  
  colnames(conditions) <- c("name", "precentageOfReplicates", "numberOfReplicateCount", "intensityValues")
  # Combine with protein Infos
  oneProtNested <- dplyr::tibble(infoProt, conditions=list(conditions))
}