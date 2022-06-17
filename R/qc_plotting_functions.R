#' Extract intensity and design from a SummarizedExperiment object
#' 
#' @param Experiment SummarizedExperiment object
#' @param log logical. Whether data should be logged before plotting
#' 
#' @export preparePlottingData
#' 
#' @import ggplot2 
#' @importFrom dplyr as_tibble
#' @importFrom SummarizedExperiment colData assay

preparePlottingData <- function(Experiment, log=FALSE){
  if(length(assay(IntensityExperiment)) > 1){
    warnings("More than one essay is present in the SummarizedExperiment. 
             Only the first one will be used.")
  }
  intensities <- assay(Experiment)
  design <- as_tibble(colData(Experiment))
  
  if(log){
    intensities <- log2(intensities+0.5)
  }
  
  list(intensities=intensities, design=design)
} 


#' Plot first 2 dimensions of PCA
#' 
#' @param Experiment SummarizedExperiment object
#' @param format 'pdf' or 'html'. Prepare image to be rendered for pdf or html Rmd output
#' @param log logical. Whether data should be logged before plotting
#' @param onlyDEProteins logical. TRUE if only DE proteins should be considered. 
#' @format "pdf" or "html". If html an interactive plot is returned using plotly. 
#' 
#' @export plot_pca_experiment
#' @details #' A protein is defined DE is the adjusted PValue is less than 0.05. 
#' When multiple comparisons are available, the adjusted PValues of the F-test is used. 
#' @import ggplot2 
#' @importFrom uuid UUIDgenerate

plot_pca_experiment <- function(Experiment, 
                                format="pdf", 
                                log=FALSE, 
                                onlyDEProteins=FALSE, 
                                title = "PCA plot"){
  # prepare data for plotting
  toPlot <- preparePlottingData(Experiment, log)
  intensities <- toPlot$intensities
  design <- toPlot$design
  
  if(onlyDEProteins){
    limmaStats <- rowData(Experiment)
    if(!("adj.P.Val" %in% colnames(limmaStats))){
      stop("No adjusted PVlaues in rowData of the Experiment provided.")
    }
    protDE <- limmaStats$ProteinId[limmaStats$adj.P.Val < 0.05]
    intensities <- intensities[rownames(intensities) %in% protDE, ]
    nDE <- nrow(intensities)
    title <- paste0("PCA plot using only N=",nDE," DE proteins.")
  }
  
  if(nrow(intensities)<=4){
    warning("Not enough DE proteins to produce a PCA plot.")
    return(NULL)
  }else{
  
    rownames(intensities) <- UUIDgenerate(use.time = NA, n = nrow(intensities))
    res.pca <- FactoMineR::PCA(t(intensities), graph = FALSE, ncp = 2)
    
    eig.val <- factoextra::get_eigenvalue(res.pca)
    eig.val <- data.table(dims = rownames(eig.val), eig.val)
    
    samples.pca <- factoextra::get_pca_ind(res.pca)
    samples.coord <- as_tibble(samples.pca$coord)
    samples.coord$SampleName = design$SampleName
    
    samples.coord <- merge(samples.coord, design)
    samples.coord <- samples.coord %>% tidyr::unite(plotSampleName, Condition, Replicate, 
                                                    sep="_", remove=FALSE)
    
    p <- ggplot(as_tibble(samples.coord), aes(x = Dim.1, y=Dim.2, colour=Condition, fill=Condition, label=Replicate)) +
      stat_ellipse(geom = "polygon", alpha=0.1) +
      geom_point(size = 3, alpha = 0.7) +
      theme_minimal() +
      scale_x_continuous(str_c("PCA 1 - ", eig.val[dims == "Dim.1", round(variance.percent,1)], "%")) +
      scale_y_continuous(str_c("PCA 2 - ", eig.val[dims == "Dim.2", round(variance.percent,1)], "%"))+
      ggtitle(title)
    
    
    if(format == "pdf"){
      p <- p + ggrepel::geom_label_repel(aes(label = plotSampleName, fill = NULL),
                                box.padding   = 0.35, 
                                point.padding = 0.5,
                                segment.color = 'grey50',
                                show.legend = FALSE)
    }else if(format == "html"){
      p <- plotly::ggplotly(p) %>% plotly::config(displayModeBar = T, 
                             modeBarButtons = list(list('toImage')),
                             displaylogo = F)
    }else{
      stop(paste0("Format ", format," not available."))
    }
    
    # Scree plot
    num_dimensions = sum(eig.val$eigenvalue>0.01)-2
    scree_plot <- ggplot(eig.val[1:num_dimensions], aes(x=reorder(dims, -`variance.percent`), y=`variance.percent`)) +
      geom_bar(stat="identity", fill = "skyblue2") +
      theme_minimal() +
      # geom_text_repel(aes(label=(round(`variance.percent`,1))), direction = 'y') +
      scale_x_discrete("PCA components") +
      scale_y_continuous("% Variance") +
      ggtitle("Scree plot")
    
    if (format == "html"){
      scree_plot <- plotly::ggplotly(scree_plot) %>%  plotly::config(displayModeBar = T, 
                                  modeBarButtons = list(list('toImage')),
                                  displaylogo = F)
    }
    
    list(p, scree_plot)
  }
  
}


#' Plot first 2 dimensions of Multi-Dimensional Scaling plot
#' 
#' @param Experiment SummarizedExperiment object
#' @param log logical. Whether data should be logged before plotting
#' 
#' @export plot_mds_experiment

plot_mds_experiment <- function(Experiment, log=FALSE){
  
  # prepare data for plotting
  toPlot <- preparePlottingData(Experiment, log)
  intensities <- toPlot$intensities
  design <- toPlot$design
  
  pi <- limma::plotMDS(intensities, plot=FALSE)
  x <- pi$x
  y <- pi$y
  ve <- round(pi$var.explained*100)
  varExpl <- data.frame(VarianceExplained = ve, 
                        Components = 1:length(ve))
  
  expNames <- colnames(intensities)
  data_plot <- data.frame(x=x, y=y, SampleName=expNames) %>% left_join(design)
  
  # MDS plot
  p <- ggplot(data_plot, aes(x = x, y=y, colour=Condition, fill=Condition, label=Replicate)) +
    stat_ellipse(geom = "polygon", alpha=0.1) +
    geom_point(size = 3, alpha = 0.7) +
    theme_minimal() +
    scale_x_continuous(str_c("leading logFC 1 - ", varExpl[varExpl$Components == 1, "VarianceExplained"], "%")) +
    scale_y_continuous(str_c("leading logFC 2 - ", varExpl[varExpl$Components == 2, "VarianceExplained"], "%"))
  
  
  num_dimensions = sum(varExpl$VarianceExplained>1)-2
  scree_plot <- ggplot(varExpl[1:num_dimensions,], aes(x=reorder(Components, rev(VarianceExplained)), y=VarianceExplained)) +
    geom_bar(stat="identity", fill = "skyblue2") +
    theme_minimal() +
    # geom_text_repel(aes(label=(round(`variance.percent`,1))), direction = 'y') +
    scale_x_discrete("Components") +
    scale_y_continuous("% Variance") +
    ggtitle("Scree plot")
  
  list(p, scree_plot)
  
} 




#' Lollipop plot of missingness by intensity column
#' @description Lollipop plot of missingness by intensity column where missingness is defined as a value equal to 0 .
#' 
#' @param Experiment SummarizedExperiment object
#' 
#' @export plot_replicate_missingness

plot_replicate_missingness <- function(Experiment){
  # prepare data for plotting
  prep_data <- preparePlottingData(Experiment)
  intensities <- prep_data[['intensities']]
  design <- prep_data[['design']]

  num.proteins <- dim(intensities)[1]
  missing.vector <- round(colSums(0 == intensities)/num.proteins * 100, 1)
  missing.table <- data.frame(SampleName = as.character(names(missing.vector)),
                              MissingValues = as.numeric(missing.vector))
  missing.table <- missing.table %>% left_join(design) %>%
    tidyr::unite(plotSampleName, Condition, Replicate, 
                 sep="_", remove=FALSE) %>%
    mutate(plotSampleName = reorder(plotSampleName, as.numeric(as.factor(Condition))))
  
  
  p <- ggplot(missing.table, aes(x = plotSampleName, y = as.numeric(MissingValues), 
                                 color = Condition, label=as.factor(plotSampleName))) +
    geom_segment( aes(x= plotSampleName, xend= plotSampleName, y=0, yend=as.numeric(MissingValues)), 
                  color="grey") +
    geom_point(size=2, alpha=0.9) +
    coord_flip() +
    theme_minimal() +
    scale_x_discrete("Sample names") +
    scale_y_continuous("% Missing Values") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.2),
          panel.grid.major.y = element_blank(),
          panel.border = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position="bottom"
    ) +
    guides(fill=guide_legend(nrow=2,byrow=TRUE))
  p
  
}


#' Lollipop plot of the number of proteins identified by replicate
#' 
#' @param Experiment SummarizedExperiment object
#' @export plot_protein_counts_by_replicate


plot_protein_counts_by_replicate <- function(Experiment){
  
  longIntensityDF <- as_tibble(SEToLongDT(Experiment))
  if(sum(str_detect(colnames(longIntensityDF), "NImputed")) != 0){
    imputedColNames <- str_detect(colnames(longIntensityDF), "NImputed")
    longIntensityDF$Imputed  <- rowSums(longIntensityDF[,imputedColNames]) > 0 
    longIntensityDF <- longIntensityDF[!longIntensityDF$Imputed, ]
  }else{
    longIntensityDF <- longIntensityDF[longIntensityDF$Intensity > 0, ]
  }
  
  dt <- longIntensityDF %>% 
    group_by(Condition, Replicate) %>%
    summarise(N = n()) %>%
    mutate(Replicate = forcats::fct_reorder(Replicate, as.numeric(as.factor(Condition))))
  
  p <- ggplot(dt, aes(x = as.factor(Replicate), y = N, color = Condition, label=Replicate)) +
    geom_segment( aes(x=as.factor(Replicate), xend=as.factor(Replicate), y=0, yend=N), 
                  color="light grey") +
    geom_point(size=2, alpha=0.9) +
    coord_flip() +
    theme_minimal() +
    scale_x_discrete("Condition - Replicate") +
    scale_y_continuous("# Identified proteins") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.2),
          panel.grid.major.y = element_blank(),
          panel.border = element_blank(),
          axis.ticks.y = element_blank()
    )
  
  p
}




#' Histogram of the distribution of missingness by protein 
#' @description Histogram of the distribution of missingness by protein  where missingness is defined as 0 values.
#'
#' @param Experiment SummarizedExperiment object
#' 
#' @export plot_protein_missingness

plot_protein_missingness <- function(Experiment){
  # prepare data for plotting
  intensities <- preparePlottingData(Experiment)[['intensities']]
  
  num.samples <- dim(intensities)[2]
  missing.vector <- rowSums(0 == intensities)
  missing.vector <- (missing.vector/num.samples)
  tot.features <- nrow(intensities) 
  
  # number of proteins with the max number of missing values
  ymax = max(table(missing.vector)) #/length(missing.vector)
  
  dt = as.data.frame(list(missing.vector = missing.vector))
  prot30perc <- sum(missing.vector <= 0.3)
  p <- ggplot(dt, aes(x = 1-missing.vector)) +
    annotate('rect', xmin = 0.7, xmax = 1.05, ymin = 0, ymax = ymax, alpha=0.2)  +
    geom_histogram(binwidth = max(0.1, round(1/max(num.samples), 2)), fill="skyblue2") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.2),
          #panel.grid.major.y = element_blank(),
          #panel.border = element_blank(),
          #axis.ticks.y = element_blank()
    ) +
    scale_x_continuous("% of measurements per protein", 
                       labels = scales::percent, 
                       limits = c(-0.1, 1.15), 
                       breaks = seq(0, 1, 0.1)) +
    labs(y = "Number of proteins") +
    annotate('text', x = 0.3, y = 0.8*ymax, 
             label=str_c("Number of proteins with\n at least 30% available values: ", 
                         prot30perc, "\n",
                         round(prot30perc/tot.features,2)*100, "% of the total proteins") )
  
  p
  
}




#' Boxplot of the relative log expression (RLE) values across samples
#' @description RLE values are computed using the data from the CompleteIntensityExperiment which contains 
#' the intensity values used for the DE analysis (log-transformed, normalised when required). 
#'
#' @param IntensityExperiment SummarizedExperiment object of the raw data, i.e. including missing values
#' @param CompleteIntensityExperiment  SummarizedExperiment object of the data used for the DE analysis 
#' (log-transformed, normalised when required). 
#' @param includeImputed logical. Whether imputed values should be considered when computing the RLE values.
#' @param plotRawRLE logical. If TRUE RLE using the raw log2 initial data is used. 
#' Otherwise, normalised data are used.
#' 
#' @export plot_rle_boxplot
#' @importFrom tidyr pivot_longer
#' @importFrom data.table data.table

plot_rle_boxplot <- function(IntensityExperiment, CompleteIntensityExperiment, 
                             includeImputed = FALSE, plotRawRLE = FALSE){
  
  
  longRawDF <- SEToLongDT(IntensityExperiment)
  longRawDF$Imputed <- longRawDF$Intensity == 0
  longRawDF <- longRawDF[,c("ProteinId","Imputed","SampleName","Intensity")]
  setnames(longRawDF, "Intensity","rawIntensity")
  longRawDF <- longRawDF[Imputed == FALSE, rawlog2Int:= log2(rawIntensity)]
  
  if(plotRawRLE){
    longIntensityDF <- longRawDF[Imputed == FALSE, 
                                       RLE := rawlog2Int - median(rawlog2Int), 
                                       by = ProteinId] 
    longIntensityDF <- longIntensityDF[!longIntensityDF$Imputed]
    
  }else{
    # the intensity plotted are the ones present in the assay experiment
    longIntensityDF <- SEToLongDT(CompleteIntensityExperiment)
    longIntensityDF <- as_tibble(longIntensityDF) %>% 
      left_join(as_tibble(longRawDF))
    longIntensityDF <- data.table(longIntensityDF[!longIntensityDF$Imputed,
                                                  c("ProteinId","Imputed", "Intensity","SampleName")])
    
    if(includeImputed){
      longIntensityDF <- longIntensityDF[,RLE := Intensity - median(Intensity), 
                                         by = ProteinId]
    }else{
      longIntensityDF <- longIntensityDF[!longIntensityDF$Imputed, 
                                         RLE := Intensity - median(Intensity), 
                                         by = ProteinId] 
      longIntensityDF <- longIntensityDF[!longIntensityDF$Imputed]
    }
  }
  
  
  # Merge with design
  design <- colData(IntensityExperiment)
  longIntensityDF <- merge(longIntensityDF, design, by = "SampleName", all.x=TRUE)
  longIntensityDF <- as_tibble(longIntensityDF) %>% tidyr::unite(plotSampleName, Condition, Replicate, 
                                                      sep="_", remove=FALSE)
  
  p = ggplot(longIntensityDF, aes(x = plotSampleName, y = RLE)) + 
    geom_boxplot(aes(fill = Condition), 
                 outlier.alpha=0.3,
                 notch = TRUE, 
                 notchwidth = 0.8) +
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.2),
          panel.grid.major.y = element_blank(),
          panel.border = element_blank(),
          axis.ticks.y = element_blank()) +
    scale_x_discrete("Sample names") +
    scale_y_continuous("RLE") +
    geom_hline(yintercept = 0, linetype="dotted")

  p

}

#' Convert output of `plot_rle_boxplot` to interactive
#' @param fig non interactive ggplot object
#' @export makeRLEBoxplotInteractive

makeRLEBoxplotInteractive <- function(fig){
  fig <- plotly::ggplotly(fig, tooltip = c("y")) %>% 
    plotly::config(displayModeBar = T, 
                   modeBarButtons = list(list('toImage')),
                   displaylogo = F)
  # this code will be needed if you want to make an interactive version
  fig$x$data <- lapply(fig$x$data, FUN = function(x){
    # When creating plot p with ggplot if you specify 
    #fill = cut use x$fill$color instead of $line$color
    x$marker$outliercolor = x$line$color 
    # When creating plot p with ggplot if you specify 
    # fill = cut use x$fill$color instead $line$color
    x$marker$color = x$line$color
    # When creating plot p with ggplot if you specify fill = 
    #cut use x$fill$color instead $line$color
    x$marker$line = x$line$color 
    return(x)
  })
  
  return(fig)
}

#' Boxplot of log expression values across samples. 
#'
#' @param Experiment SummarizedExperiment object of raw intensities (pre log-transformation)
#' @export plot_log_measurement_boxplot
#' @importFrom stringr str_detect
#' 
#' @details Zero values are treated as missing values and excluded. 

plot_log_measurement_boxplot <- function(Experiment){ 
  
  longIntensityDF <- as_tibble(SEToLongDT(Experiment))
  longIntensityDF <- longIntensityDF[longIntensityDF$Intensity > 0, ]
  longIntensityDF$Intensity <- log2(longIntensityDF$Intensity)
  
  dataToPlot <- as_tibble(longIntensityDF) %>% tidyr::unite(plotSampleName, 
                                                            Condition, Replicate, 
                                                            sep="_", remove=FALSE)
  
  p <- ggplot(dataToPlot , aes(x=plotSampleName, y=Intensity, 
                                    fill=Condition)) +
    geom_boxplot(outlier.alpha=0.3, notch = TRUE, notchwidth = 0.8) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.2),
          panel.grid.major.y = element_blank(),
          panel.border = element_blank(),
          axis.ticks.y = element_blank()
    ) +
    geom_hline(yintercept = 0, linetype="dotted", colour="grey") +
    scale_x_discrete("Sample names") +
    scale_y_continuous("Log2 Intensity")
  
  p
}


#' Density distribution of intensity values
#'
#' @param Experiment SummarizedExperiment object
#' @param log logical. If TRUE log2 transformation is applied to the intensity 
#' before plotting
#' 
#' @export plot_density_distr
plot_density_distr <- function(Experiment, log=FALSE){
  # prepare data for plotting
  toPlot <- preparePlottingData(Experiment, log)
  intensities <- toPlot$intensities
  design <- toPlot$design
  
  long_int <- as_tibble(intensities) %>% pivot_longer(cols = colnames(intensities), 
                                           names_to = "SampleName", 
                                           values_to = "Intensity")
  long_int <- long_int %>% left_join(design)
  p = ggplot(long_int, aes(x = Intensity, fill = Condition, colour=Condition)) + 
    geom_density(alpha=0.6) +
    theme_minimal() +
    labs(x = "log2 Intensity", y = "Density")
  p
}


#' Distribution of imputed versus non imputed values. 
#'
#' @param CompleteIntensityExperiment SummarizedExperiment object of the data used for the DE analysis with limma
#' @param byCondition logical. TRUE to produce separate density distributions by Condition.
#' @export plot_imputed_vs_not

plot_imputed_vs_not <- function(CompleteIntensityExperiment, byCondition=FALSE,
                                title = "Intensity (Imputed vs Actual)"){
  
  assay1 <- as_tibble(assay(CompleteIntensityExperiment))
  assay1$ProteinId <- rownames(assay(CompleteIntensityExperiment))
  long1 <- assay1 %>% 
    pivot_longer(all_of(colnames(assay(CompleteIntensityExperiment))),
                 names_to = "SampleName",values_to = "Intensity")
  
  assay2 <- as_tibble(assays(CompleteIntensityExperiment)[[2]])
  assay2$ProteinId <- rownames(assays(CompleteIntensityExperiment)[[2]])
  long2 <- assay2 %>% 
    pivot_longer(all_of(colnames(assays(CompleteIntensityExperiment)[[2]])),
                 names_to = "SampleName",values_to = "Imputed")
  
  long3 <- long1 %>% left_join(long2)
  design <- colData(CompleteIntensityExperiment)
  long3 <- long3 %>% left_join(as_tibble(design))
  long3$Imputed <- as.factor(long3$Imputed)
  long3$Imputed <- long3$Imputed == 1
  
  p <- ggplot(long3, aes(x=Intensity, fill=Imputed, 
                                   colour = Imputed)) +
    geom_density(alpha=0.4) +
    theme_minimal() +
    ggtitle(title) +
    labs(x = "Log2 Intensity")
  
  if(byCondition){
    p <- p + facet_wrap(~Condition)
  }
  
  p
  
}



#' CV distributuions of a protein summarized by condition of interest
#' @param Experiment SummarizedExperiment object
#' @export plot_condition_cv_distribution

plot_condition_cv_distribution <- function(Experiment, title = "Protein Intensity CV"){
  longIntensityDF <- SEToLongDT(Experiment)
  longIntensityDF$Imputed <- longIntensityDF$Intensity == 0
  
  cvdt <- longIntensityDF[Imputed == 0][, countRep := .N, by = .(ProteinId, Condition)]
  cvdt[, countRepMax := max(countRep), by = .(ProteinId, Condition)]
  cvdt[, ReplicatePC := countRep/countRepMax]
  cvdt[, intensity := as.double(Intensity)]
  cvdt <- cvdt[ReplicatePC >= 0.5]
  
  cvdt <- cvdt[ReplicatePC >= 0.5, .(cv = sd(intensity)/mean(intensity)), by = .(ProteinId, Condition)]
  
  p <- ggplot(cvdt, aes(x=cv, fill=Condition, colour=Condition)) +
    geom_density(alpha=0.4) +
    theme_minimal() +
    scale_x_continuous("% CV", labels = scales::percent) +
    ggtitle(title)
  
  list(p, as_tibble(cvdt))
  
}


#' Volcano plot of results
#' @param comparison.statistics data.frame containing 
#' `FoldChange`, `PValue`, `AdjustedPValue` and `ProteinId` columns 
#' @export plot_volcano

plot_volcano <- function(comparison.statistics){
  p <- ggplot(comparison.statistics, 
              aes(FoldChange, -log10(PValue), fdr = AdjustedPValue, ProteinId = ProteinId)) + 
    geom_point() +
    ggtitle(comparison.statistics$comparison[1]) +
    geom_hline(yintercept=-log10(0.05)) 
  p
}



#' MA plot of results
#' @param comparison.statistics data.frame containing 
#' `FoldChange`, `AveExpr`, and `ProteinId` columns.
#' @export plot_ma

plot_ma <- function(comparison.statistics){
  p <- ggplot(comparison.statistics, 
              aes(x = AveExpr, y = FoldChange, ProteinId = ProteinId)) + 
    geom_point() +
    ggtitle(comparison.statistics$comparison[1]) +
    geom_hline(yintercept=0)
  p
}


#' Heatmap of the binary intensity matrix of missing values 
#' @param Experiment SummarizedExperiment object of raw data. 
#' Intensities values equal to 0 are considered as missing values.
#' @param complete logical. TRUE to use all the matrix of intensities. 
#' FALSE to include only features with at least one missing value across samples.   
#' @export plot_heatmap_missingness
#' @import pheatmap

plot_heatmap_missingness <- function(Experiment, complete=FALSE){
  data <- assay(Experiment)
  if(!complete){
    data <- data[rowSums(data==0) > 0,]
  }
  mt_binary <- ifelse(data != 0,1,0)
  design <- colData(Experiment)[c("SampleName", "Condition")]
  column_annotation <- data.frame(SampleName = colnames(mt_binary)) %>% 
    left_join(as_tibble(design)) 
  rownames(column_annotation) <- column_annotation$SampleName
  column_annotation$SampleName <- NULL
  colour_cells <- grDevices::colorRampPalette(c("#e0f3db", "#43a2ca"))(2)
  rownames(mt_binary) <- NULL
  return(pheatmap(mt_binary, color = colour_cells, 
                  annotation_col = column_annotation, fontsize = 8, 
                  show_colnames=TRUE, legend_labels = c("0","1"), legend_breaks = 0:1))
  
}




#' Correlation plot of samples using DE proteins
#' @param Experiment SummarizedExperiment object of imputed data. 
#' @param onlyDEProteins logical. TRUE if only DE proteins should be considered. 
#' @export plot_samples_correlation_matrix
#' 
#' @details #' A protein is defined DE is the adjusted PValue is less than 0.05. 
#' When multiple comparisons are available, the adjusted PValues of the F-test is used. 
#' At least 5 DE proteins are required to produce the correlation plot.  
#' 

plot_samples_correlation_matrix <- function(Experiment, onlyDEProteins=FALSE, title = "All proteins"){
  
  # prepare data for plotting
  toPlot <- preparePlottingData(Experiment, log=FALSE)
  intensities <- toPlot$intensities
  design <- as_tibble(toPlot$design) %>% tidyr::unite(plotSampleName, 
                                                      Condition, Replicate, 
                                                      sep="_", remove=FALSE)
  
  
  column_name_order <- sapply(design$SampleName, function(z) which(colnames(intensities) %in% z))
  colnames(intensities)[column_name_order] <- design$plotSampleName
  
  if(onlyDEProteins){
    limmaStats <- rowData(Experiment)
    if(!("adj.P.Val" %in% colnames(limmaStats))){
      stop("No adjusted PVlaues in rowData of the Experiment provided.")
    }
    protDE <- limmaStats$ProteinId[limmaStats$adj.P.Val < 0.05]
    intensities <- intensities[rownames(intensities) %in% protDE, ]
    nDE <- nrow(intensities)
    title <- paste0("Using only N=", nDE, " proteins")
  }
  
  if(nrow(intensities) > 4){
    DT_corMatrix <- Hmisc::rcorr(intensities)
    DT_corMatrix <- DT_corMatrix$r
    
    DT_corMatrix[DT_corMatrix <= -1] = -1
    DT_corMatrix[DT_corMatrix >= 1] = 1
  
    
    corrplot::corrplot(DT_corMatrix, type = "upper", method = "square",
             title = title,
             tl.cex = 0.5, mar = c(0,0,1.5,0),
             tl.col = "black", order = "hclust")
  }else{
    warning("Not enough DE proteins to produce a correlation plot.")
    return(NULL)
  }
  
}
