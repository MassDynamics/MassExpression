#' This function creates a summarized experiment object from a protein intensity table and experiment design. 
#' 
#' @param experimentDesign data.frame. Experiment design provided in input by the user. Required columsn are: `SampleName` and `Condition`.
#' @param proteinIntensities data.frame. Wide matrix of intensities. Rows are proteins and columns are SampleNames. Required column: `ProteinId`. 
#' @param listMetadata list of metadata: `Species`, `LabellingMethod`, `NormalisationAppliedToAssay`.
#' @return A SummarizedExperiment object
#' @export
#' @importFrom SummarizedExperiment SummarizedExperiment rowData colData
#' @importFrom S4Vectors SimpleList
#' @importFrom dplyr group_by mutate row_number

createSummarizedExperiment <- function(experimentDesign, 
                                       proteinIntensities, 
                                       listMetadata, 
                                       subjectCol = "Subject",
                                       techReplCol = "TechRepl"){
  
  print("Sanity check experiment design and detect experiment type.")
  # check and update the design
  if(!("Condition" %in% colnames(experimentDesign))){
    stop("'Condition' column is not available in the experiment design.")
  }
  experimentDesign$Condition <- as.character(experimentDesign$Condition)
  
  if(!("SampleName" %in% colnames(experimentDesign))){
    stop("'SampleName' column is not available in the experiment design.")
  }
  experimentDesign$SampleName <- as.character(experimentDesign$SampleName)
  
  freqSampleName <- table(experimentDesign$SampleName)
  if(any(freqSampleName > 1)){
    stop("SampleName in experiment design has to be unique.")
  }
  
  if(!(all(experimentDesign$SampleName %in% colnames(proteinIntensities)))){
    missing_samples <- experimentDesign$SampleName[!(experimentDesign$SampleName %in% colnames(proteinIntensities))] 
    stop(paste0("The following SampleNames in the experiment design are missing from the protein intensities table: ", 
         paste(missing_samples, collapse = ", ")))
  }
  
  if(!("ProteinId" %in% colnames(proteinIntensities))){
    stop("'ProteinId' column is not available in the protein intensities table.")
  }
  proteinIntensities$ProteinId <- as.character(proteinIntensities$ProteinId)
  
  # Update experiment design
  # Add subject columns
  if(!(subjectCol %in% colnames(experimentDesign))){
    experimentDesign[, subjectCol] <- 1:nrow(experimentDesign)
    print("Subject information is added automatically assuming samples derive from different subjects.")
  }
  
  designTypeInfo <- detectDesignType(experimentDesign = experimentDesign, 
                                     techReplCol = "TechRepl", subjectCol = "Subject")
  experimentDesign <- designTypeInfo$experimentDesign
  expType <- designTypeInfo$expType
  
  # Add Replicate column
  if(is.null(expType$condition2Name)){
    experimentDesign = as.data.table(experimentDesign)[,Replicate:=1:.N, by = c("Condition")]
  }else{
    cond2Name <- expType$condition2Name
    experimentDesign = as.data.table(experimentDesign)[,Replicate:=1:.N, by = c("Condition", cond2Name)]
  }
  
  # metadata 
  if(!(listMetadata$NormalisationAppliedToAssay %in% c("None", "Median"))){
    stop("normalisationMethod should be one one of: 'Median' or `None`.")
  }
  listMetadata <- list(inputMetadata = listMetadata, experimentType = expType)
  
  print("Prepare rowData")
  rowDataPossible <-  c("ProteinId","GeneName","Description")
  rowDataPresent <- intersect(rowDataPossible, colnames(proteinIntensities))
  rowDataAbsent <- rowDataPossible[!(rowDataPossible %in% rowDataPresent)]
  
  assayData <- as.matrix(proteinIntensities[,experimentDesign$SampleName])
  colnames(assayData) <- experimentDesign$SampleName

  rownames(assayData) <- proteinIntensities$ProteinId
  
  rowFeatures <- proteinIntensities[,rowDataPresent,drop=FALSE] 
  rowFeatures[,rowDataAbsent] <- NA
  
  # construct summarized experiment object
  IntensityExperiment <- SummarizedExperiment(rowData = rowFeatures,
                                              assays= SimpleList(raw=assayData),
                                              colData = experimentDesign, 
                                              metadata = listMetadata)
  
  colnames(IntensityExperiment) <- experimentDesign$SampleName
  rownames(IntensityExperiment) <- rowData(IntensityExperiment)$ProteinId
  stopifnot(colnames(IntensityExperiment) == 
              SummarizedExperiment::colData(IntensityExperiment)$SampleName)
  
  return(IntensityExperiment)
}

#' Detect the design type 
#' 
#' @param experimentDesign data.frame. Experiment design provided in input by the user.
#' @param techReplCol Name of Technical Replicate column name. Deafult `TechRepl`
#' @param subjectCol Name of Subject column name. Deafult `Subject`

#' @export

detectDesignType <- function(experimentDesign, techReplCol, subjectCol){
  
  expType <- NULL
  cond2Name <- NULL
  cond2Levels <- NULL
  conditionTimeOrDose = FALSE
  
  if(hasConditionOnly(experimentDesign)){
    cond1Name <- getConditionName(experimentDesign)
    cond1Levels <- unique(experimentDesign[, cond1Name, drop=TRUE])
    repCondName <- cond1Name
    
  } else if (hasConditionWithTimeDose(experimentDesign)) {
    condNames <- getConditionName(experimentDesign)
    cond1Name <- "Condition"
    cond2Name <- condNames[!(condNames %in% "Condition")]
    repCondName <- cond2Name
    
    conditionTimeOrDose = TRUE
    cond1Levels <- unique(experimentDesign[, cond1Name, drop=TRUE])
    cond2Levels <- unique(experimentDesign[, cond2Name, drop=TRUE])
  } else {
    stop("Experiment design not supported. Check details for allowed column names.")
  }
  
  # check tech repl
  experimentDesign <- addTechRepl(experimentDesign, 
                                  techReplCol = techReplCol, 
                                  subjectCol = subjectCol, 
                                  cond1Name = cond1Name, 
                                  cond2Name = cond2Name)
  print(experimentDesign)
  techRepl <- sum(table(experimentDesign[,techReplCol]) > 1) > 0
  
  # check subject for repeated measurements
  repMeas <-  hasRepMeas(experimentDesign, 
                         repCondName = repCondName, 
                         subjectCol = subjectCol)
  
  # Number of levels in condition
  if(!is.null(cond2Name)){
    enoughLevels <- hasEnoughLevelsInConditions(experimentDesign, cond2Name)
    if(!enoughLevels){
      stop(paste0("Condition: ", cond2Name, " in the experiment design should have at least 2 levels."))
    }
    enoughLevels <- hasEnoughLevelsInConditions(experimentDesign, cond1Name)
    if(!enoughLevels){
      stop(paste0("Condition: ", cond1Name, " in the experiment design should have at least 2 levels."))
    }
  }else{
    enoughLevels <- hasEnoughLevelsInConditions(experimentDesign, cond1Name)
    if(!enoughLevels){
      stop(paste0("Condition: ", cond1Name, " in the experiment design should have at least 2 levels."))
    }
  }
  
  # check enough replicates within condition
  enoughRepl <- hasEnoughNumberReplPerGroup(experimentDesign, 
                                         cond1Name = cond1Name, 
                                         cond2Name = cond2Name)
  if(!enoughRepl){
    warning("The experiment might not have enough biological replicates within the given conditions 
            to fit stable linear models.")
  }
  
  expType <- list(conditionOnly = hasConditionOnly(experimentDesign), 
                  conditionTimeOrDose = conditionTimeOrDose,
                  condition1Name = cond1Name,
                  condition1Levels = cond1Levels,
                  condition2Name = cond2Name,
                  condition2Levels = cond2Levels,
                  hasEnoughBioRepl = enoughRepl, 
                  hasTechRepl = techRepl, 
                  hasRepMeas = repMeas)
  
  list(expType = expType, experimentDesign = experimentDesign)
}

#' Add technical replicate information to design
#' @keywords internal
#' @noRd

addTechRepl <- function(experimentDesign, 
                        techReplCol, 
                        subjectCol, 
                        cond1Name, 
                        cond2Name = NULL){
  
  if(!(techReplCol %in% colnames(experimentDesign))){
    # Technical replicates
    if(is.null(cond2Name)){
      experimentDesign[,techReplCol] <- paste0(experimentDesign[,subjectCol, drop=TRUE],"-", 
                                               experimentDesign[,cond1Name, drop=TRUE])
    } else {
      experimentDesign[,techReplCol] <- paste0(experimentDesign[,subjectCol, drop=TRUE],"-", 
                                               experimentDesign[,cond1Name, drop=TRUE], "-", 
                                               experimentDesign[,cond2Name, drop=TRUE])
    }
  }
  experimentDesign[,techReplCol] <- as.numeric(as.factor(experimentDesign[,techReplCol, drop=TRUE]))
  
  return(experimentDesign)
}

#' Check for repeated measurements in design
#' @keywords internal
#' @noRd

hasRepMeas <- function(experimentDesign,
                       subjectCol = "Subject",
                       repCondName){
  
  hasReps = FALSE
  hasSubjCol <- sum(subjectCol %in% colnames(experimentDesign)) == 1
  if(hasSubjCol){
    # Repeated measurements have different conditions but same subject 
    repMeas <- experimentDesign %>% 
      group_by(get(subjectCol)) %>%
      summarise(reps = length(unique(get(repCondName))))
    hasReps <- sum(repMeas$reps > 1) > 0
  } 
  return(hasReps)
}


#' Check for enough replicates given the design
#' @keywords internal
#' @noRd

hasEnoughNumberReplPerGroup <- function(experimentDesign, cond1Name, cond2Name = NULL){
  if(is.null(cond2Name)){
    freqRepl <- table(experimentDesign[,cond1Name])
    enoughRepl <- as.logical(1 - sum(freqRepl < 3) )
  } else {
    freqRepl <- table(experimentDesign[,cond1Name], experimentDesign[,cond2Name])
    enoughRepl <- as.logical(1 - (sum(freqRepl < 3) > 0) )
  }
  return(enoughRepl)
}

#' Check for enough levels in conditions in the experiment design
#' @keywords internal
#' @noRd

hasEnoughLevelsInConditions <- function(experimentDesign, condName){
  len_levels_condition <- length(names(table(experimentDesign[,condName])))
  return(len_levels_condition > 1)
}

#' Check for only one condition in experiment design
#' @keywords internal
#' @noRd

hasConditionOnly <- function(inputDesign){
  colnamesInputDesign <- colnames(inputDesign)
  return(sum(c("Condition","Time", "Dose") %in% colnamesInputDesign) == 1)
}

#' Check for only one condition + time/dose in experiment design
#' @keywords internal
#' @noRd

hasConditionWithTimeDose <- function(inputDesign){
  colnamesInputDesign <- colnames(inputDesign)
  opt1 <- sum(c("Condition","Time") %in% colnamesInputDesign)
  opt2 <- sum(c("Condition", "Dose") %in% colnamesInputDesign)
  return(opt1 == 2 | opt2 == 2)
}

#' Extract condition names from experiment design
#' @keywords internal
#' @noRd

getConditionName <- function(inputDesign){
  colnamesInputDesign <- colnames(inputDesign)
  which_condition <- c("Condition","Time", "Dose") %in% colnamesInputDesign
  return(c("Condition","Time", "Dose")[which_condition])
}