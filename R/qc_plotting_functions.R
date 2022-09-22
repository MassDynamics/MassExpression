#' Extract intensity and design from a SummarizedExperiment object
#' 
#' @param Experiment SummarizedExperiment object
#' @param assayName name of assay to use
#' @param log logical. Whether data should be logged before plotting
#' 
#' @export preparePlottingData
#' 
#' @import ggplot2 
#' @importFrom dplyr as_tibble
#' @importFrom SummarizedExperiment colData assay

preparePlottingData <- function(Experiment,  assayName="intensities", log=FALSE){

  intensities <- assays(Experiment)[[assayName]]
  design <- as_tibble(colData(Experiment))
  
  if(log){
    intensities <- log2(intensities+0.5)
  }
  
  list(intensities=intensities, design=design)
} 


#' Subset features to be used for PCA
#' 
#' @param Experiment SummarizedExperiment object.
#' @param auto_select_features str. One of 'de' or empty string.

#' @importFrom SummarizedExperiment rowData

select_features_for_pca <- function(Experiment,
                                auto_select_features=NULL){
  
  if(is.null(auto_select_features)){
    return(Experiment)
  }
  
  if(!(auto_select_features %in% c("de"))){
    stop(paste0("Invalid `auto_select_features` argument:", auto_select_features))
  }
  
  if(auto_select_features == "de"){
    limmaStats <- rowData(Experiment)
    if(!("adj.P.Val" %in% colnames(limmaStats))){
      stop("No adjusted PVlaues in rowData of the Experiment provided.")
    }
    protDE <- limmaStats$ProteinId[limmaStats$adj.P.Val < 0.05]
    Experiment <- Experiment[rownames(Experiment) %in% protDE, ]
  }
  
  # TODO
  # add select highly variable genes
  # add provide vector of proteins
  
  return(Experiment)
}

#' Compute PCAs
#' 
#' @param Experiment SummarizedExperiment object
#' @param assayName name of assay to use
#' @param ndim number of dimensions kept in result 

#' @importFrom uuid UUIDgenerate
#' @import FactoMineR
#' @import factoextra

compute_pcas <- function(Experiment, assayName, log, ndim=2){
  
  # prepare data to compute PC
  toCompute <- preparePlottingData(Experiment, assayName, log)
  intensities <- toCompute$intensities
  design <- toCompute$design
  
  rownames(intensities) <- UUIDgenerate(use.time = NA, n = nrow(intensities))
  res.pca <- FactoMineR::PCA(t(intensities), graph = FALSE, ncp = ndim)
  
  eig.val <- factoextra::get_eigenvalue(res.pca)
  eig.val <- data.table(dims = rownames(eig.val), eig.val)
  
  samples.pca <- factoextra::get_pca_ind(res.pca)
  samples.coord <- as_tibble(samples.pca$coord)
  samples.coord$SampleName = design$SampleName
  
  samples.coord <- merge(samples.coord, design)
  
  list(pcas = samples.coord, eigenval = eig.val)
}  


#' Plot principal components
#' 
#' @param Experiment SummarizedExperiment object
#' @param assayName name of assay to use
#' @param dimPlot vector of size 2 specifying the dimensions to plot.  
#' @param log logical. Whether data should be logged before plotting
#' @param auto_select_features str. One of 'de' (differentially expressed features see Details), 
#' 'hvf' (highly variable features). If not provided all features are kept.
#' @param format 'pdf' or 'html'. Prepare image to be rendered for pdf or html Rmd output
#' @param title_pca str. title on PCA plot
#' @param title_screeplot str. title on screeplot

#' @export plot_chosen_pca_experiment
#' @details #' A protein is defined DE if the adjusted PValue of the t-test or ANOVA (with multiple groups)
#' is less than 0.05. 
#' 
#' @import ggplot2 
#' @importFrom uuid UUIDgenerate
#' @importFrom stringr str_c

plot_chosen_pca_experiment <- function(Experiment, 
                                assayName="intensities",
                                dimPlot = c(1,2), 
                                log=FALSE, 
                                auto_select_features=NULL, 
                                title_pca = "PCA plot",
                                title_screeplot = "Scree plot",
                                format="html"){
  
  # Subset experiment with features required
  Experiment <- select_features_for_pca(Experiment, auto_select_features = auto_select_features)
  
  if(nrow(Experiment)<=4){
    warning("Not enough features to produce a PCA plot (less than 5).")
    return(NULL)
  }
  
  if(!is.null(auto_select_features)){
    if(auto_select_features == "de"){
      nDE <- nrow(Experiment)
      title_pca <- paste0("PCA plot using only N=",nDE," DE proteins.")
    }
  }
  
  # Calculate PCs
  pcaToPlot <- compute_pcas(Experiment, assayName, log, ndim=max(dimPlot))
  dim1 <- paste0("Dim.",dimPlot[1])
  dim2 <- paste0("Dim.",dimPlot[2])
  
  samples.coord <- pcaToPlot$pcas[,c(dim1, dim2, "Condition", "Replicate")]
  samples.coord <- samples.coord %>% tidyr::unite(plotSampleName, Condition, Replicate, 
                                                  sep="_", remove=FALSE)
  eig.val <- pcaToPlot$eigenval
    
  # Prepare plot
  p <- ggplot(samples.coord, aes(x = get(dim1), y=get(dim2), colour=Condition, 
                                 fill=Condition, 
                                 label=Replicate)) +
    geom_point(size = 3, alpha = 0.7) +
    stat_ellipse(geom = "polygon", alpha=0.1) +
    theme_minimal() +
    scale_x_continuous(str_c("PCA",  dimPlot[1], " - ", eig.val[dims == dim1, round(variance.percent,1)], "%")) +
    scale_y_continuous(str_c("PCA", dimPlot[2], " - ", eig.val[dims == dim2, round(variance.percent,1)], "%")) +
    ggtitle(title_pca) +
    theme(plot.title = element_text(hjust = 0.5))
  
  
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
    ggtitle(title_screeplot) + 
    theme(plot.title = element_text(hjust = 0.5))
  
  if (format == "html"){
    scree_plot <- plotly::ggplotly(scree_plot) %>%  plotly::config(displayModeBar = T, 
                                                                   modeBarButtons = list(list('toImage')),
                                                                   displaylogo = F)
  }
  
  list(p, scree_plot)
  
}



#' Plot first 2 dimensions of PCA
#' 
#' @param Experiment SummarizedExperiment object
#' @param assayName name of assay to use
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
                                assayName = "intensities",
                                format="pdf", 
                                log=FALSE, 
                                onlyDEProteins=FALSE, 
                                title = "PCA plot"){
  # prepare data for plotting
  toPlot <- preparePlottingData(Experiment, assayName=assayName, log=log)
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
      scale_y_continuous("% Variance")
    
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
#' @param assayName name of assay to use
#' @param log logical. Whether data should be logged before plotting
#' 
#' @export plot_mds_experiment

plot_mds_experiment <- function(Experiment, assayName="intensities", log=FALSE){
  
  # prepare data for plotting
  toPlot <- preparePlottingData(Experiment, assayName = assayName, log=log)
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
#' @param assayName name of assay to use
#' @param title str. Plot title 
#' 
#' @export plot_replicate_missingness

plot_replicate_missingness <- function(Experiment, assayName="raw", title = "Missingness by samples using protein measurements"){
  # prepare data for plotting
  prep_data <- preparePlottingData(Experiment, assayName = assayName)
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
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5)) + 
    guides(fill=guide_legend(nrow=2,byrow=TRUE))
  p
  
}


#' Lollipop plot of the number of proteins identified by replicate
#' 
#' @param Experiment SummarizedExperiment object
#' @param assayName name of assay to use
#' @param title str. Plot title 
#' 
#' @import dplyr
#' @import tidyr
#' @import tibble
#' @import ggplot2

#' @export plot_n_identified_proteins_by_replicate


plot_n_identified_proteins_by_replicate <- function(Experiment, 
                                             assayName = "raw", 
                                             title = "# Identified proteins by sample"){
  
  longIntensityDF <- as_tibble(SEToLongDT(Experiment, assayName = assayName))
  if(sum(str_detect(colnames(longIntensityDF), "NImputed")) != 0){
    imputedColNames <- str_detect(colnames(longIntensityDF), "NImputed")
    longIntensityDF$Imputed  <- rowSums(longIntensityDF[,imputedColNames]) > 0 
    longIntensityDF <- longIntensityDF[!longIntensityDF$Imputed, ]
  }else{
    longIntensityDF <- longIntensityDF[longIntensityDF$Intensity > 0, ]
  }
  
  longIntensityDF$plotSampleName <- as.factor(paste(longIntensityDF$Condition, longIntensityDF$Replicate, sep = "_"))
  dt <- longIntensityDF %>% 
    group_by(plotSampleName) %>%
    summarise(N = n())
  
  dt <- dt %>% left_join(unique(longIntensityDF[,c("plotSampleName", "Condition")])) %>%
    mutate(plotSampleName = reorder(plotSampleName, as.numeric(as.factor(Condition))))
    
  
  p <- ggplot(dt, aes(x = plotSampleName, y = N, color = Condition, label=plotSampleName)) +
    geom_segment( aes(x=plotSampleName, xend=plotSampleName, y=0, yend=N), 
                  color="light grey") +
    geom_point(size=2, alpha=0.9) +
    coord_flip() +
    theme_minimal() +
    scale_x_discrete("Condition - Replicate") +
    scale_y_continuous("# Identified proteins") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.2),
          panel.grid.major.y = element_blank(),
          panel.border = element_blank(),
          axis.ticks.y = element_blank(), 
          legend.position="bottom"
    ) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5)) + 
    guides(fill=guide_legend(nrow=2,byrow=TRUE))
  
  p
}


#' Lollipop plot of the number of proteins identified across all replicates of a condition
#' 
#' @param Experiment SummarizedExperiment object
#' @param assayName name of assay to use
#' @param condition_colname str. Name of grouping condition.
#' @param title str. Plot title 
#' 
#' @export plot_consistent_proteins_by_replicate

plot_consistent_proteins_by_replicate <- function(Experiment, 
                                                    assayName = "raw", 
                                                    condition_colname = "Condition", 
                                                    title = "# Consistently identified proteins by condition"){
  
  longIntensityDF <- as_tibble(SEToLongDT(Experiment, assayName = assayName))
  if(sum(str_detect(colnames(longIntensityDF), "NImputed")) != 0){
    imputedColNames <- str_detect(colnames(longIntensityDF), "NImputed")
    longIntensityDF$Imputed  <- rowSums(longIntensityDF[,imputedColNames]) > 0 
    longIntensityDF <- longIntensityDF[!longIntensityDF$Imputed, ]
  }else{
    longIntensityDF <- longIntensityDF[longIntensityDF$Intensity > 0, ]
  }
  
  # Replicates in each condition
  replicateCond <- longIntensityDF %>% 
    group_by(get(condition_colname)) %>%
    summarise(NRepl = length(unique(Replicate)))
  colnames(replicateCond)[1] <- condition_colname
  
  # Number of avail proteins in each conditions
  NAvailCond <- longIntensityDF %>% 
    group_by(get(condition_colname), ProteinId) %>%
    summarise(NAvail = n())
  colnames(NAvailCond)[1] <- condition_colname
  
  # Complete proteins in condition
  NAvailCond <- NAvailCond %>% 
    left_join(replicateCond)
  NAvailCond <- NAvailCond[NAvailCond$NAvail == NAvailCond$NRepl,]
  
  dt <- NAvailCond %>% 
    group_by(get(condition_colname)) %>%
    summarise(N = n())
  colnames(dt)[1] <- condition_colname
  
  p <- ggplot(dt, aes(x = .data[[condition_colname]], 
                      y = N, 
                      color = .data[[condition_colname]], 
                      label = .data[[condition_colname]])) +
    geom_segment( aes(x=.data[[condition_colname]], xend=.data[[condition_colname]], y=0, yend=N), 
                  color="light grey") +
    geom_point(size=2, alpha=0.9) +
    coord_flip() +
    theme_minimal() +
    scale_x_discrete(condition_colname) +
    scale_y_continuous("# Consistently identified proteins in condition") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.2),
          panel.grid.major.y = element_blank(),
          panel.border = element_blank(),
          axis.ticks.y = element_blank(), 
          legend.position="bottom"
    ) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5)) + 
    guides(fill=guide_legend(nrow=2,byrow=TRUE))
  
  p
}



#' Heatmap showing pattern of missingness
#' 
#' @param Experiment SummarizedExperiment object
#' @param assayName name of assay to use
#' @param condition_colname str. Name of grouping condition.
#' @param title str. Plot title 
#' 
#' @description Columns (samples) have been clustered based on pattern of missingness. 
#' 
#' @import ComplexHeatmap
#' @import SummarizedExperiment
#' 
#' @export plot_missingness_heatmap

plot_missingness_heatmap <- function(Experiment, 
                                     assayName = "raw", 
                                     condition_colname = "Condition", 
                                     title = "Missingness pattern"){
  
  condition <- colData(Experiment)[, condition_colname]
  
  if (dim(Experiment)[1] > 10000) {
    set.seed(255)
    random_sample <- sample(1:dim(Experiment)[1], 10000)
    Experiment <- Experiment[random_sample, ]
  }
  
  y <- assays(Experiment)[[assayName]]
  y_missing = t(apply(y, 1, function(x) ifelse(x == 0, 1, 0)))
  
  ha_column <- HeatmapAnnotation(Condition = condition)
  
  hm <- Heatmap(y_missing,
                column_title = title,
                name = "Intensity",
                col = c("#8FBC8F", "#FFEFDB"),
                show_row_names = FALSE,
                show_column_names = FALSE,
                cluster_rows = TRUE,
                cluster_columns = TRUE,
                show_column_dend = FALSE,
                show_row_dend = FALSE,
                top_annotation = ha_column,
                row_names_gp =  grid::gpar(fontsize = 7),
                column_names_gp = grid::gpar(fontsize = 8),
                heatmap_legend_param = list(#direction = "horizontal",
                  heatmap_legend_side = "bottom",
                  labels = c("missing","observed"),
                  legend_width = unit(6, "cm")),
  )
  hm <- draw(hm, heatmap_legend_side = "left")
}

#' Histogram of the distribution of missingness by protein 
#' @description Histogram of the distribution of missingness by protein  where missingness is defined as 0 values.
#'
#' @param Experiment SummarizedExperiment object
#' @param assayName name of assay to use
#' @param title str. Plot title
#' 
#' @export plot_protein_missingness

plot_protein_missingness <- function(Experiment, assayName="raw", 
                                     title = "Protein missingness"){
  # prepare data for plotting
  intensities <- preparePlottingData(Experiment, assayName = assayName)[['intensities']]
  
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
                         round(prot30perc/tot.features,2)*100, "% of the total proteins") ) +
    ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
  
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
#' @param title str. Plot title
#' @param format str. 'pdf' or 'html'
#' 
#' @export plot_rle_boxplot
#' @importFrom tidyr pivot_longer
#' @importFrom data.table data.table

plot_rle_boxplot <- function(IntensityExperiment, CompleteIntensityExperiment, 
                             includeImputed = FALSE, 
                             plotRawRLE = FALSE, 
                             title="RLE plot", 
                             format = "html"){
  
  
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
    longIntensityDF <- SEToLongDT(CompleteIntensityExperiment, assayName = "intensities")
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
    geom_hline(yintercept = 0, linetype="dotted") +
    ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
  
  
  if(format =="html"){
    plotly::ggplotly(p) %>% plotly::config(displayModeBar = T, 
                                             modeBarButtons = list(list('toImage')),
                                             displaylogo = F)
  } else{
    p
  }

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
#' @param title str. Plot title
#' @param format str. 'pdf' or 'html'
#' 
#' @export plot_log_measurement_boxplot
#' @importFrom stringr str_detect
#' @details Zero values are treated as missing values and excluded. 

plot_log_measurement_boxplot <- function(Experiment, 
                                         title = "Intensities distribution", 
                                         format = "html"){ 
  
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
    scale_y_continuous("Log2 Intensity") +
    ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
  
  if(format =="html"){
    plotly::ggplotly(p) %>% plotly::config(displayModeBar = T, 
                                           modeBarButtons = list(list('toImage')),
                                           displaylogo = F)
  } else{
    p
  }
}


#' Density distribution of intensity values
#'
#' @param Experiment SummarizedExperiment object
#' @param assayName name of assay to use
#' @param log logical. If TRUE log2 transformation is applied to the intensity 
#' before plotting
#' 
#' @export plot_density_distr
plot_density_distr <- function(Experiment, assayname, log=FALSE){
  # prepare data for plotting
  toPlot <- preparePlottingData(Experiment, assayName, log)
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
#' @param title str. Plot title
#' @param format str. 'pdf' or html'
#' 
#' @export plot_imputed_vs_not

plot_imputed_vs_not <- function(CompleteIntensityExperiment, byCondition=FALSE,
                                title = "Intensity (Imputed vs Actual)",
                                format = "html"){
  
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
    theme(plot.title = element_text(hjust = 0.5)) +
    labs(x = "Log2 Intensity") +
    ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
  
  if(byCondition){
    p <- p + facet_wrap(~Condition)
  }
  
  if(format =="html"){
    plotly::ggplotly(p) %>% plotly::config(displayModeBar = T, 
                                           modeBarButtons = list(list('toImage')),
                                           displaylogo = F)
  } else{
    p
  }
  
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
    ggtitle(title) + 
    theme(plot.title = element_text(hjust = 0.5))
    
  
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




#' Correlation plot of samples using DE proteins
#' @param Experiment SummarizedExperiment object of imputed data. 
#' @param assayName name of assay to use
#' @param onlyDEProteins logical. TRUE if only DE proteins should be considered. 
#' @export plot_samples_correlation_matrix
#' 
#' @details #' A protein is defined DE is the adjusted PValue is less than 0.05. 
#' When multiple comparisons are available, the adjusted PValues of the F-test is used. 
#' At least 5 DE proteins are required to produce the correlation plot.  
#' 

plot_samples_correlation_matrix <- function(Experiment, assayName="intensities",
                                            onlyDEProteins=FALSE, title = "All proteins"){
  
  # prepare data for plotting
  toPlot <- preparePlottingData(Experiment, assayName, log=FALSE)
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
    title <- paste0("Using only N=", nDE, " DE proteins")
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
    warning("Not enough DE proteins to produce a correlation plot. At least 4 are needed.")
    return(NULL)
  }
  
}
