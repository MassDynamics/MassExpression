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

ensure_package_installed("roxygen2", repos = list("http://cran.rstudio.com/", "https://cran.ms.unimelb.edu.au/"))
ensure_package_installed("devtools", repos = list("http://cran.rstudio.com/", "https://cran.ms.unimelb.edu.au/"))
ensure_package_installed("grDevices", repos = list("http://cran.rstudio.com/", "https://cran.ms.unimelb.edu.au/"))


ensure_package_installed_with_version("gert", "1.3.1",repos = "http://cran.rstudio.com/")
ensure_package_installed_with_version("usethis","2.0.1", repos = "http://cran.rstudio.com/")
ensure_package_installed_with_version("data.table", "1.14.0",repos = "http://cran.rstudio.com/")
ensure_package_installed_with_version("foreach", "1.5.1", repos = "http://cran.rstudio.com/")
ensure_package_installed_with_version("stringr", "1.4.0", repos = "http://cran.rstudio.com/")
ensure_package_installed_with_version("dplyr", "1.0.7", repos = "http://cran.rstudio.com/")
ensure_package_installed_with_version("jsonlite", "1.7.2", repos = "http://cran.rstudio.com/")
ensure_package_installed_with_version("tidyr", "1.1.3", repos = "http://cran.rstudio.com/")
ensure_package_installed_with_version("plotly", "4.9.4.1", repos = "http://cran.rstudio.com/")
ensure_package_installed_with_version("ggplot2", "3.3.5", repos = "http://cran.rstudio.com/")
ensure_package_installed_with_version("forcats", "0.5.1", repos = "http://cran.rstudio.com/")
ensure_package_installed_with_version("knitr", "1.33", repos = "http://cran.rstudio.com/")
ensure_package_installed_with_version("rmarkdown","2.9", repos = "http://cran.rstudio.com/")
ensure_package_installed_with_version("statmod","1.4.36", repos = "http://cran.rstudio.com/")
ensure_package_installed_with_version("FactoMineR", '2.4', repos = "http://cran.rstudio.com/")
ensure_package_installed_with_version("factoextra", '1.0.7', repos = "http://cran.rstudio.com/")
ensure_package_installed_with_version("testthat", '3.0.4', repos = "http://cran.rstudio.com/")
ensure_package_installed_with_version("ggrepel", '0.9.1', repos = "http://cran.rstudio.com/")
ensure_package_installed_with_version("pheatmap", '1.0.12', repos = "http://cran.rstudio.com/")
ensure_package_installed_with_version("Hmisc", '4.5-0', repos = "http://cran.rstudio.com/")
ensure_package_installed_with_version("corrplot", '0.90', repos = "http://cran.rstudio.com/")
ensure_package_installed_with_version("lme4", '1.1-28', repos = "http://cran.rstudio.com/")
ensure_package_installed_with_version("pbkrtest", '0.5.1', repos = "http://cran.rstudio.com/")
ensure_package_installed_with_version("car", '3.0-12', repos = "http://cran.rstudio.com/")
ensure_package_installed_with_version("DT", '0.18', repos = "http://cran.rstudio.com/")
ensure_package_installed_with_version("Rcpp", '1.0.7', repos = "http://cran.rstudio.com/")
ensure_package_installed_with_version("RcppArmadillo", '0.10.6.0.0', repos = "http://cran.rstudio.com/")
ensure_package_installed_with_version("readr", '2.0.0', repos = "http://cran.rstudio.com/")
ensure_package_installed_with_version("credentials", '1.3.1', repos = "http://cran.rstudio.com/")
ensure_package_installed_with_version("matrixStats", '0.60.0', repos = "http://cran.rstudio.com/")
ensure_package_installed_with_version("broom", '0.7.9', repos = "http://cran.rstudio.com/")
ensure_package_installed_with_version("uuid", '1.0-3', repos = "http://cran.rstudio.com/")
ensure_package_installed_with_version("parallel", '4.1.0', repos = "http://cran.rstudio.com/")


#### Install packages Bioconductor
ensure_package_installed_with_version("BiocManager", '1.30.16', repos = "http://cran.rstudio.com/")
BiocManager::install(c("Biobase", "limma", "S4Vectors", "GenomicRanges", "GenomeInfoDb", "SummarizedExperiment"), ask=FALSE, type="source")

biocManager_valid <- BiocManager::valid()
if (typeof(biocManager_valid) == 'list'){
  BiocManager::valid()$out_of_date
  print(BiocManager::valid()$out_of_date)
  sessionInfo()
  stop("BiocManager is not valid")
}

sessionInfo()
