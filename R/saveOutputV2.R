#' MassExpression output to MassDynamics V2 Protein Counts JSON
#'
#' @description Save output from MassExpression and MDFlexiComparisons in the V2 formats required for protein counts intensities json.
#'
#' @export MassExpressionTOProtCountV2
#'
#' @importFrom SummarizedExperiment rowData
#' @importFrom S4Vectors metadata
#' @import data.table

MassExpressionTOProtCountV2 <- function(CompleteIntensityExperiment,
                                        longIntensityDT,
                                        outputFolder){
  dir.create(file.path(outputFolder), recursive = TRUE, showWarnings = FALSE)
  
  comparisonExperiments <-
    MassExpression::listComparisonExperiments(CompleteIntensityExperiment)
  
  # Write Replicates tab data
  SummarizedExperiment::rowData(CompleteIntensityExperiment)$GroupId <- 1:nrow(CompleteIntensityExperiment)
  GroupLabelType <- "Protein"
  GroupColumnName <- "ProteinId"
  
  limmaStats <- SummarizedExperiment::rowData(CompleteIntensityExperiment)[, c("GroupId", "ProteinId",
                                                                               "GeneName", "Description")]
  longDTProt <- merge(longIntensityDT, limmaStats, by=GroupColumnName, all.x=TRUE)
  longDTProt <- data.table::as.data.table(longDTProt)
  
  #three_proteins <- longDTProt[longDTProt$GroupId %in% c(1,2,3),]
  writeReplicateDataV2(longDTProt, outputFolder,
                            GroupLabelType = GroupLabelType,
                            GroupColumnName = GroupColumnName)
  
}


#' Write protein_counts_and_intensity.json
#' This is a wrapper of `oneProteinReplDataGeneric`, looping over all proteins.
#'
#' @param longDTProt data.table protein information stored in long format.
#' Columns required: `ProteinId`, `GeneName`, `Description`, `log2NInt`, `Condition`,
#'  `Replicate`, `Imputed`.
#' @param outputFolder str. Path to folder where `data` should be saved.
#'
#' @importFrom jsonlite write_json
#' @export writeReplicateDataV2

writeReplicateDataV2 <- function(longDTProt, outputFolder, GroupColumnName, GroupLabelType){
  proteinSet <- unique(longDTProt[, GroupId])
  # prot <- proteinSet[1]
  protList <- lapply(proteinSet, function(prot){
    oneGroupReplDataV2(longDTProt[GroupId %in% prot,], GroupColumnName, GroupLabelType)
  })
  
  
  protDF <- do.call(rbind, protList)
  
  dir.create(outputFolder, showWarnings = FALSE)
  outPath = file.path(outputFolder,"protein_counts_and_intensity.json")
  write_json(protDF, outPath, digits = NA, na = "null")
}



#' Transform data for one protein from a long format to the nested data structure needed for the Replicates tab.
#' @param oneProt data.table single protein information stored in long format.
#' Columns required: `ProteinId`, `GeneName`, `Description`, `log2NInt`, `Condition`,
#'  `Replicate`, `Imputed`.
#' @export oneGroupReplDataV2
#'
oneGroupReplDataV2 <- function(oneProt, GroupColumnName, GroupLabelType){
  infoProt <- unique(oneProt[,c("GroupId","ProteinId", "GeneName", "Description")])
  infoProt <- infoProt[, GroupLabel:=get(GroupColumnName)] #TODO clarify qhat's this field
  infoProt <- infoProt[, GroupLabelType:=GroupLabelType]
  setnames(infoProt, "ProteinId", "ProteinIds")
  setnames(infoProt, "GeneName", "GeneNames")
  infoProt <- infoProt[,c("GroupId", "GroupLabel","GroupLabelType","ProteinIds", "GeneNames", "Description")]
  
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
    setnames(oneCondRepl, old = "SampleName", new = "replicate")
    oneCondRepl <- oneCondRepl[,c("replicate", "log2NInt_ProteinGroupId", "Imputed")]
    
    entryCond <- dplyr::tibble(infoOneCond, intensityValues=list(oneCondRepl))
    conditions[cond_idx, ] <- entryCond
  }
  
  colnames(conditions) <- c("name", "precentageOfReplicates", "numberOfReplicateCount", "intensityValues")
  # Combine with protein Infos
  oneProtNested <- dplyr::tibble(infoProt, conditions=list(conditions))
}

