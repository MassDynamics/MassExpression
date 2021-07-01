#!/usr/bin/Rscript

ensure_package_installed <- function (package, repos = repos) {
  if(!require(package, character.only=TRUE)) {
    install.packages(package, repos = repos)
    library(package, character.only=TRUE)
  }
}

ensure_package_installed_with_version <- function (package, version, repos = repos) {
  if(!require(package, character.only=TRUE)) {
    install_version(package, version = version, repos = repos)
    library(package, character.only=TRUE)
  } else if (packageVersion(package) != version) {
    install_version(package, version = version, repos = repos)
    library(package, character.only=TRUE)
  }

  if(packageVersion(package) != version) {
    stop("issue with: ", package, " ", version)
  }

}


ensure_package_installed("devtools", repos = "http://cran.rstudio.com/")


ensure_package_installed("gert", repos = "http://cran.rstudio.com/")
ensure_package_installed("usethis", repos = "http://cran.rstudio.com/")
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
ensure_package_installed_with_version("FactoMineR", '2.4', repos = "http://cran.rstudio.com/")
ensure_package_installed_with_version("factoextra", '1.0.7', repos = "http://cran.rstudio.com/")

#### Install packages Bioconductor
ensure_package_installed("BiocManager", repos = "http://cran.rstudio.com/")
BiocManager::install(c("Biobase", "limma", "S4Vectors", "GenomicRanges", "GenomeInfoDb", "SummarizedExperiment"))

if (BiocManager::valid() != TRUE){
  BiocManager::valid()$out_of_date
  sessionInfo()
  stop("BiocManager is not valid")
}

sessionInfo()