
#' Plot first 2 dimensions of PCA
#' 
#' @param intensities Matrix of intensities (rows are features, columns are samples)
#' @param design Experimental design
#' @param log logical
#' 
#' @export pca_plot_experiment
#' 
#' @import ggplot2 
#' @importFrom dplyr as_data_frame 

pca_plot_experiment <- function(intensities, design, log=FALSE){
  
  if(log){
    intensities <- log2(intensities+0.5)
  }
  
  res.pca <- FactoMineR::PCA(t(intensities), graph = FALSE, ncp = 2)
  
  eig.val <- factoextra::get_eigenvalue(res.pca)
  eig.val <- data.table(dims = rownames(eig.val), eig.val)
  
  samples.pca <- factoextra::get_pca_ind(res.pca)
  samples.coord <- as_data_frame(samples.pca$coord)
  samples.coord$IntensityColumn = design$IntensityColumn
  
  samples.coord <- merge(samples.coord, design)
  print(colnames(samples.coord))
  
  
  p <- ggplot(as_data_frame(samples.coord), aes(x = Dim.1, y=Dim.2, colour=Condition, fill=Condition, label=Replicate)) +
    stat_ellipse(geom = "polygon", alpha=0.1) +
    geom_point(size = 3, alpha = 0.7) +
    theme_minimal() +
    scale_x_continuous(str_c("PCA 1 - ", eig.val[dims == "Dim.1", round(variance.percent,1)], "%")) +
    scale_y_continuous(str_c("PCA 2 - ", eig.val[dims == "Dim.2", round(variance.percent,1)], "%"))
  
  
  num_dimensions = sum(eig.val>0.01)-2
  scree_plot <- ggplot(eig.val[1:num_dimensions], aes(x=reorder(dims, -`variance.percent`), y=`variance.percent`)) +
    geom_bar(stat="identity", fill = "skyblue2") +
    theme_minimal() +
    # geom_text_repel(aes(label=(round(`variance.percent`,1))), direction = 'y') +
    scale_x_discrete("PCA components") +
    scale_y_continuous("% Variance") +
    ggtitle("Scree plot")
  
  list(plotly::ggplotly(p), scree_plot)
  
}


#' Plot first 2 dimensions of Multi-Dimensional Scaling plot
#' 
#' @param intensities Matrix of intensities (rows are features, columns are samples)
#' @param design Experimental design
#' @param log logical
#' 
#' @export mds_plot_experiment

mds_plot_experiment <- function(intensities, design, log=FALSE){
  
  if(log){
    intensities <- log2(intensities+0.5)
  }
  
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
  
  list(plotly::ggplotly(p), scree_plot)
  
} 




#' this function creates a lollipop plot of missingness by intensity column where missingness is defined as a value equal to 0 .
#' 
#' @param intensities Matrix of intensities (rows are features, columns are samples)
#' @param design Experimental design
#' 
#' @export replicate_missingness_experiment

replicate_missingness_experiment <- function(intensities, design){
  
  missing.vector <- colSums(0 == intensities)
  missing.table <- as_data_frame(cbind(names(missing.vector), missing.vector))
  colnames(missing.table) <- c("IntensityColumn", "MissingValues")
  missing.table <- merge(missing.table, design)
  missing.table <- as_data_frame(missing.table)
  
  
  p <- ggplot(missing.table, aes(x = Replicate, y = as.numeric(MissingValues), color = Condition, label=as.factor(Replicate))) +
    geom_segment( aes(x= Replicate, xend= Replicate, y=0, yend=as.numeric(MissingValues)), color="grey") +
    geom_point(size=2, alpha=0.9) +
    coord_flip() +
    theme_minimal() +
    scale_x_discrete("Replicate") +
    scale_y_continuous("# Missing Values") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.2),
          panel.grid.major.y = element_blank(),
          panel.border = element_blank(),
          axis.ticks.y = element_blank()
    )
  p
  
}

#' Histogram of the distribution of missingness by protein where missingness is defined as 0 values.
#' 
#' @param intensities Matrix of intensities (rows are features, columns are samples)
#' 
#' @export protein_missingness_experiment


protein_missingness_experiment <- function(intensities){
  
  num.samples <- dim(intensities)[2]
  missing.vector <- rowSums(0 == intensities)
  missing.vector <- 1-missing.vector/num.samples
  tot.features <- nrow(intensities) 
  
  # number of proteins with the max number of missing values
  ymax = max(table(missing.vector)) #/length(missing.vector)
  
  dt = as.data.frame(list(missing.vector = missing.vector))
  prot30perc <- sum(missing.vector <= 0.3)
  p <- ggplot(dt, aes(x = missing.vector)) +
    annotate('rect', xmin = -0.05, xmax = 0.3, ymin = 0, ymax = ymax, alpha=0.2)  +
    geom_histogram(binwidth = max(0.1, round(1/max(num.samples), 2)), fill="skyblue2") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.2),
          panel.grid.major.y = element_blank(),
          panel.border = element_blank(),
          axis.ticks.y = element_blank()
    ) +
    scale_x_continuous("% of missing values per protein", 
                       labels = scales::percent, 
                       limits = c(-0.05, 1.15), 
                       breaks = seq(0, 1, 0.1)) +
    labs(y = "Number of proteins") +
    annotate('text', x = 0.5, y = 0.8*ymax, 
             label=str_c("Number of proteins with\n <= 30% missing values: ", prot30perc, "\n",
                         round(prot30perc/tot.features,2)*100, "% of total proteins") )
  
  p
  
}

