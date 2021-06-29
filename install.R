#!/usr/bin/Rscript
install.packages("devtools", repos = "http://cran.rstudio.com/")
install.packages("data.table", repos = "http://cran.rstudio.com/")
install.packages("foreach", repos = "http://cran.rstudio.com/")
install.packages("stringr", repos = "http://cran.rstudio.com/")
install.packages("dplyr", repos = "http://cran.rstudio.com/")
install.packages("S4Vectors", repos = "http://cran.rstudio.com/")
install.packages("jsonlite", repos = "http://cran.rstudio.com/")
install.packages("tidyr", repos = "http://cran.rstudio.com/")
install.packages("plotly", repos = "http://cran.rstudio.com/")
install.packages("ggplot2", repos = "http://cran.rstudio.com/")
install.packages("forcats", repos = "http://cran.rstudio.com/")
install.packages("knitr", repos = "http://cran.rstudio.com/")
install.packages("rmarkdown", repos = "http://cran.rstudio.com/")


#### Install packages Bioconductor

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("Biobase", "limma", "SummarizedExperiment"))
