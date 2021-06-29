#!/usr/bin/Rscript

ensure_package_installed <- function (package, repos = repos) {
  install.packages(package, repos = repos)
  library(package, character.only=TRUE)
}

ensure_package_installed("devtools", repos = "http://cran.rstudio.com/")
ensure_package_installed("data.table", repos = "http://cran.rstudio.com/")
ensure_package_installed("foreach", repos = "http://cran.rstudio.com/")
ensure_package_installed("stringr", repos = "http://cran.rstudio.com/")
ensure_package_installed("dplyr", repos = "http://cran.rstudio.com/")
ensure_package_installed("S4Vectors", repos = "http://cran.rstudio.com/")
ensure_package_installed("jsonlite", repos = "http://cran.rstudio.com/")
ensure_package_installed("tidyr", repos = "http://cran.rstudio.com/")
ensure_package_installed("plotly", repos = "http://cran.rstudio.com/")
ensure_package_installed("ggplot2", repos = "http://cran.rstudio.com/")
ensure_package_installed("forcats", repos = "http://cran.rstudio.com/")
ensure_package_installed("knitr", repos = "http://cran.rstudio.com/")
ensure_package_installed("rmarkdown", repos = "http://cran.rstudio.com/")


#### Install packages Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("Biobase", "limma", "GenomicRanges", "GenomeInfoDb", "SummarizedExperiment"))
