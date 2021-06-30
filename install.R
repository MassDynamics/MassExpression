#!/usr/bin/Rscript

ensure_package_installed <- function (package, repos = repos) {
  if(!require(package, character.only=TRUE)) {
    install.packages(package, repos = repos)
    library(package, character.only=TRUE)
  }
}

ensure_package_installed("gert", repos = "http://cran.rstudio.com/")
ensure_package_installed("usethis", repos = "http://cran.rstudio.com/")
ensure_package_installed("devtools", repos = "http://cran.rstudio.com/")
ensure_package_installed("data.table", repos = "http://cran.rstudio.com/")
ensure_package_installed("foreach", repos = "http://cran.rstudio.com/")
ensure_package_installed("stringr", repos = "http://cran.rstudio.com/")
ensure_package_installed("dplyr", repos = "http://cran.rstudio.com/")
ensure_package_installed("jsonlite", repos = "http://cran.rstudio.com/")
ensure_package_installed("tidyr", repos = "http://cran.rstudio.com/")
ensure_package_installed("plotly", repos = "http://cran.rstudio.com/")
ensure_package_installed("ggplot2", repos = "http://cran.rstudio.com/")
ensure_package_installed("forcats", repos = "http://cran.rstudio.com/")
ensure_package_installed("knitr", repos = "http://cran.rstudio.com/")
ensure_package_installed("rmarkdown", repos = "http://cran.rstudio.com/")
ensure_package_installed("statmod", repos = "http://cran.rstudio.com/")

#### Install packages Bioconductor
ensure_package_installed("BiocManager", repos = "http://cran.rstudio.com/")
BiocManager::install(c("Biobase", "limma", "S4Vectors", "GenomicRanges", "GenomeInfoDb", "SummarizedExperiment"))

if (!BiocManager::valid()){
  stop("BiocManager is not valid")
}
