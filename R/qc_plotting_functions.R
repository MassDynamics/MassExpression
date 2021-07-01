#' Extract the essay and colData (design data frame) from a SummarizedExperiment object
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
#' @param log logical. Whether data should be logged before plotting
#' 
#' @export pca_plot_experiment
#' 
#' @import ggplot2 

pca_plot_experiment <- function(Experiment, log=FALSE){
  # prepare data for plotting
  toPlot <- preparePlottingData(Experiment, log)
  intensities <- toPlot$intensities
  design <- toPlot$design
  
  res.pca <- FactoMineR::PCA(t(intensities), graph = FALSE, ncp = 2)
  
  eig.val <- factoextra::get_eigenvalue(res.pca)
  eig.val <- data.table(dims = rownames(eig.val), eig.val)
  
  samples.pca <- factoextra::get_pca_ind(res.pca)
  samples.coord <- as_tibble(samples.pca$coord)
  samples.coord$IntensityColumn = design$IntensityColumn
  
  samples.coord <- merge(samples.coord, design)
  
  p <- ggplot(as_tibble(samples.coord), aes(x = Dim.1, y=Dim.2, colour=Condition, fill=Condition, label=Replicate)) +
    stat_ellipse(geom = "polygon", alpha=0.1) +
    geom_point(size = 3, alpha = 0.7) +
    theme_minimal() +
    scale_x_continuous(str_c("PCA 1 - ", eig.val[dims == "Dim.1", round(variance.percent,1)], "%")) +
    scale_y_continuous(str_c("PCA 2 - ", eig.val[dims == "Dim.2", round(variance.percent,1)], "%"))+
    ggtitle("PCA plot")
  
  
  num_dimensions = sum(eig.val$eigenvalue>0.01)-2
  scree_plot <- ggplot(eig.val[1:num_dimensions], aes(x=reorder(dims, -`variance.percent`), y=`variance.percent`)) +
    geom_bar(stat="identity", fill = "skyblue2") +
    theme_minimal() +
    # geom_text_repel(aes(label=(round(`variance.percent`,1))), direction = 'y') +
    scale_x_discrete("PCA components") +
    scale_y_continuous("% Variance") +
    ggtitle("Scree plot")
  
  list(p, scree_plot)
  
}


#' Plot first 2 dimensions of Multi-Dimensional Scaling plot
#' 
#' @param Experiment SummarizedExperiment object
#' @param log logical. Whether data should be logged before plotting
#' 
#' @export mds_plot_experiment

mds_plot_experiment <- function(Experiment, log=FALSE){
  
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
  data_plot <- data.frame(x=x, y=y, IntensityColumn=expNames) %>% left_join(design)
  
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




#' this function creates a lollipop plot of missingness by intensity column where missingness is defined as a value equal to 0 .
#' 
#' @param Experiment SummarizedExperiment object
#' 
#' @export replicate_missingness_experiment

replicate_missingness_experiment <- function(Experiment){
  # prepare data for plotting
  prep_data <- preparePlottingData(Experiment)
  intensities <- prep_data[['intensities']]
  design <- prep_data[['design']]

  num.proteins <- dim(intensities)[1]
  missing.vector <- round(colSums(0 == intensities)/num.proteins * 100, 1)
  missing.table <- data.frame(IntensityColumn = as.character(names(missing.vector)),
                              MissingValues = as.numeric(missing.vector))
  missing.table <- missing.table %>% left_join(design) %>%
    mutate(Replicate = reorder(Replicate, as.numeric(as.factor(Condition))))
  
  p <- ggplot(missing.table, aes(x = Replicate, y = as.numeric(MissingValues), 
                                 color = Condition, label=as.factor(Replicate))) +
    geom_segment( aes(x= Replicate, xend= Replicate, y=0, yend=as.numeric(MissingValues)), 
                  color="grey") +
    geom_point(size=2, alpha=0.9) +
    coord_flip() +
    theme_minimal() +
    scale_x_discrete("Replicate") +
    scale_y_continuous("% Missing Values") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.2),
          panel.grid.major.y = element_blank(),
          panel.border = element_blank(),
          axis.ticks.y = element_blank()
    )
  p
  
}


#' Lollipop plot of the number of proteins identified by replicate
#' 
#'
#' @param Experiment SummarizedExperiment object
#' @export protein_counts_by_replicate


protein_counts_by_replicate <- function(Experiment){
  
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




#' Histogram of the distribution of missingness by protein where missingness is defined as 0 values.
#'
#' @param Experiment SummarizedExperiment object
#' 
#' @export protein_missingness_experiment

protein_missingness_experiment <- function(Experiment){
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




#' Boxplot of the distribution of the relative log expression (RLE) values across samples. 
#' RLE values are computed using the data from the CompleteIntensityExperiment which contains 
#' the intensity values used for the DE analysis (log-transformed, normalised when required).   
#'
#' @param IntensityExperiment SummarizedExperiment object of the raw data, i.e. including missing values
#' @param CompleteIntensityExperiment  SummarizedExperiment object of the data used for the DE analysis 
#' (log-transformed, normalised when required). 
#' @param includeImputed logical. Whether imputed values should be considered when computing the RLE values.
#'  
#' @export plot_rle_boxplot
#' @importFrom tidyr pivot_longer
#' @importFrom data.table data.table

plot_rle_boxplot <- function(IntensityExperiment, CompleteIntensityExperiment, 
                             includeImputed = FALSE){
  longRawDF <- SEToLongDT(IntensityExperiment)
  longRawDF$Imputed <- longRawDF$Intensity == 0
  longRawDF <- longRawDF[,c("ProteinId","Imputed","IntensityColumn")]
  
  # the intensity plotted are the ones present in the assay experiment
  longIntensityDF <- SEToLongDT(CompleteIntensityExperiment)
  longIntensityDF <- as_tibble(longIntensityDF) %>% 
    left_join(as_tibble(longRawDF))
  longIntensityDF <- data.table(longIntensityDF[!longIntensityDF$Imputed,
                                                c("ProteinId","Imputed", "Intensity","IntensityColumn")])
  
  if(includeImputed){
    longIntensityDF <- longIntensityDF[,RLE := Intensity - median(Intensity), 
                                       by = ProteinId]
  }else{
    longIntensityDF <- longIntensityDF[!longIntensityDF$Imputed, 
                                       RLE := Intensity - median(Intensity), 
                                       by = ProteinId] 
  }
  
  # Merge with design
  design <- colData(CompleteIntensityExperiment)
  longIntensityDF <- merge(longIntensityDF, design, by = "IntensityColumn", all.x=TRUE)
  
  p = ggplot(longIntensityDF, aes(x = IntensityColumn, y = RLE)) + 
    geom_boxplot(aes(fill = Condition)) + 
    theme_minimal() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.2),
          panel.grid.major.y = element_blank(),
          panel.border = element_blank(),
          axis.ticks.y = element_blank()) +
    scale_x_discrete("Replicate") +
    scale_y_continuous("RLE") +
    ggtitle("Intensities are centered on protein's medians") + 
    geom_hline(yintercept = 0, linetype="dotted")

  p

}

#' Boxplot of the distribution of the log expression values across samples. 
#'
#' @param Experiment SummarizedExperiment object
#' @param use_imputed logical
#' @export plot_measurement_boxplot
#' @importFrom stringr str_detect

plot_measurement_boxplot <- function(Experiment, log=FALSE){ 
  
  longIntensityDF <- as_tibble(SEToLongDT(Experiment))
  if(log){
    longIntensityDF$Intensity <- log2(longIntensityDF$Intensity+0.5)
  }
  
  p <- ggplot(longIntensityDF , aes(x=Replicate, y=Intensity, 
                                    fill=Condition)) +
    geom_boxplot() +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.2),
          panel.grid.major.y = element_blank(),
          panel.border = element_blank(),
          axis.ticks.y = element_blank()
    ) +
    geom_hline(yintercept = 0, linetype="dotted", colour="grey") +
    scale_x_discrete("Replicate") +
    scale_y_continuous("Log2 Intensity")
  
  p
}


#' Density distribution of intensity values
#'
#' @param Experiment SummarizedExperiment object
#' 
#' @export plot_density_distr
plot_density_distr <- function(Experiment, log=FALSE){
  # prepare data for plotting
  toPlot <- preparePlottingData(Experiment, log)
  intensities <- toPlot$intensities
  design <- toPlot$design
  
  long_int <- as_tibble(intensities) %>% pivot_longer(cols = colnames(intensities), 
                                           names_to = "IntensityColumn", 
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
#' @param IntensityExperiment SummarizedExperiment object of raw data
#' @param CompleteIntensityExperiment SummarizedExperiment object of the data used for the DE analysis with limma
#' @export plot_imputed_vs_not

plot_imputed_vs_not <- function(IntensityExperiment, CompleteIntensityExperiment){
  longRawDF <- SEToLongDT(IntensityExperiment)
  longRawDF$Imputed <- longRawDF$Intensity == 0
  longRawDF <- longRawDF[,c("ProteinId","Imputed","IntensityColumn")]
  
  # the intensity plotted are the ones present in the assay experiment
  longIntensityDF <- SEToLongDT(CompleteIntensityExperiment)
  longIntensityDF <- as_tibble(longIntensityDF) %>% left_join(as_tibble(longRawDF))
  
  p <- ggplot(longIntensityDF, aes(x=Intensity, fill=Imputed, 
                                   colour = Imputed)) +
   facet_wrap(~Condition) +
    geom_density(alpha=0.4) +
    theme_minimal() +
    ggtitle("Intensity (Imputed vs Not)") +
    labs(x = "Log2 Intensity")
  
  p
  
}

#' CV distributuions of a protein summarized by condition of interest
#' @param Experiment SummarizedExperiment object
#' @export plot_condition_cv_distribution

plot_condition_cv_distribution <- function(Experiment){
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
    ggtitle("Protein Intensity CV")
  
  p
  
}


#' Volcano plot of results
#' @param Experiment SummarizedExperiment object
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
#' @param Experiment SummarizedExperiment object
#' @export plot_ma

plot_ma <- function(comparison.statistics){
  p <- ggplot(comparison.statistics, 
              aes(x = AveExpr, y = FoldChange, ProteinId = ProteinId)) + 
    geom_point() +
    ggtitle(comparison.statistics$comparison[1]) +
    geom_hline(yintercept=0)
  p
}