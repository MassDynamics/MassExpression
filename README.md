MassExpression
================

-   [Workflow example](#workflow-example)
    -   [Setup](#setup)
    -   [Internal data available](#internal-data-available)
    -   [Run workflow with sample data](#run-workflow-with-sample-data)
-   [Generate vignette](#generate-vignette)
-   [Save output to display in the app and create QC
    reports](#save-output-to-display-in-the-app-and-create-qc-reports)
    -   [Render QC](#render-qc)
-   [Session Information](#session-information)

[![DOI](https://zenodo.org/badge/377716166.svg)](https://zenodo.org/badge/latestdoi/377716166)

Universal Imports + High Quality QC + Differential Expression Analysis =
Awesomeness.

# Workflow example

## Setup

``` r
library(devtools)
library(tibble)
library(plotly)
library(stringr)
```

``` r
devtools::install_github("MassDynamics/MassExpression")
```

``` r
library(MassExpression)
utils::packageVersion("MassExpression")
#> [1] '0.0.78'
```

## Internal data available

``` r
data(package="MassExpression")
```

## Run workflow with sample data

``` r
intensities <- mq_lfq_data$intensities
design <- mq_lfq_data$design
parameters <- mq_lfq_data$parameters
normalisation_method <- parameters[parameters[,1] == "UseNormalisationMethod",2]
species <- parameters[parameters[,1] == "Species",2]
labellingMethod <- parameters[parameters[,1] == "LabellingMethod",2]


results <- runGenericDiscovery(experimentDesign = design, 
                               proteinIntensities = intensities, 
                               normalisationMethod = normalisation_method, 
                               species = species, 
                               labellingMethod = labellingMethod)

IntensityExperiment <- results$IntensityExperiment
CompleteIntensityExperiment <-  results$CompleteIntensityExperiment
longIntensityDT <- results$longIntensityDT
design <- colData(CompleteIntensityExperiment)
```

`results` is a list containing two `SummarizedExperiment` objects:

-   `IntensityExperiment`: contains the raw data (including missing
    values)

-   `CompleteIntensityExperiment`: contains the imputed data and summary
    statistics about the number of replicates and imputed proteins in
    each group of the conditions of interest.

# Generate vignette

``` r
tools::buildVignettes(dir = ".", tangle=TRUE)
```

# Save output to display in the app and create QC reports

``` r
output_folder <- "path/to/output/folder"

saveOutput(IntensityExperiment = IntensityExperiment, 
CompleteIntensityExperiment = CompleteIntensityExperiment,
longIntensityDT = longIntensityDT, 
outputFolder =  output_folder)
```

## Render QC

``` r
# Render and save QC report 
qc_report <- system.file("rmd","QC_report.Rmd", package = "MassExpression")

rmarkdown::render(qc_report,
                  params = list(listInt = results,
                                experiment = "Mass Dynamics QC report",
                                output_figure = file.path(output_folder, "figure_html/"),
                                format = "html"),
                  output_file = file.path(output_folder, "QC_Report.html"),
                  output_format=rmarkdown::html_document(
                            self_contained=FALSE,
                            lib_dir=file.path(output_folder,"qc_report_files"),
                            code_folding= "hide",
                            theme="united",
                            toc = TRUE,
                            toc_float = TRUE,
                            fig_caption= TRUE,
                            df_print="paged"))
# Render PDF
rmarkdown::render(qc_report,
                  params = list(listInt = results,
                                experiment = "Mass Dynamics QC report",
                                output_figure = file.path(output_folder_pdf, "figure_pdf/"),
                                format = "pdf"),
                  output_file = file.path(output_folder_pdf, "QC_Report.pdf"),
                  output_format=rmarkdown::pdf_document(
                    toc = TRUE,
                    fig_caption= TRUE))
```

# Session Information

``` r
sessionInfo()
#> R version 4.1.0 (2021-05-18)
#> Platform: x86_64-apple-darwin17.0 (64-bit)
#> Running under: macOS Big Sur 10.16
#> 
#> Matrix products: default
#> BLAS:   /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRblas.dylib
#> LAPACK: /Library/Frameworks/R.framework/Versions/4.1/Resources/lib/libRlapack.dylib
#> 
#> locale:
#> [1] en_AU.UTF-8/en_AU.UTF-8/en_AU.UTF-8/C/en_AU.UTF-8/en_AU.UTF-8
#> 
#> attached base packages:
#> [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
#> [8] methods   base     
#> 
#> other attached packages:
#>  [1] MassExpression_0.0.78       SummarizedExperiment_1.22.0
#>  [3] GenomicRanges_1.44.0        GenomeInfoDb_1.28.4        
#>  [5] IRanges_2.26.0              S4Vectors_0.30.2           
#>  [7] MatrixGenerics_1.4.3        matrixStats_0.61.0         
#>  [9] Biobase_2.52.0              BiocGenerics_0.38.0        
#> 
#> loaded via a namespace (and not attached):
#>  [1] tidyselect_1.1.2       xfun_0.29              purrr_0.3.4           
#>  [4] lattice_0.20-45        colorspace_2.0-2       vctrs_0.3.8           
#>  [7] generics_0.1.2         htmltools_0.5.2        yaml_2.2.2            
#> [10] utf8_1.2.2             rlang_1.0.1            pillar_1.7.0          
#> [13] glue_1.6.2             DBI_1.1.2              RColorBrewer_1.1-2    
#> [16] uuid_1.0-3             GenomeInfoDbData_1.2.6 foreach_1.5.2         
#> [19] lifecycle_1.0.1        stringr_1.4.0          zlibbioc_1.38.0       
#> [22] munsell_0.5.0          gtable_0.3.0           codetools_0.2-18      
#> [25] evaluate_0.14          knitr_1.37             fastmap_1.1.0         
#> [28] fansi_1.0.2            scales_1.1.1           limma_3.48.3          
#> [31] DelayedArray_0.18.0    jsonlite_1.8.0         XVector_0.32.0        
#> [34] ggplot2_3.3.5          digest_0.6.29          stringi_1.7.6         
#> [37] dplyr_1.0.8            grid_4.1.0             cli_3.2.0             
#> [40] tools_4.1.0            bitops_1.0-7           magrittr_2.0.2        
#> [43] RCurl_1.98-1.6         tibble_3.1.6           tidyr_1.2.0           
#> [46] crayon_1.5.0           pkgconfig_2.0.3        pheatmap_1.0.12       
#> [49] ellipsis_0.3.2         Matrix_1.4-0           data.table_1.14.2     
#> [52] assertthat_0.2.1       rmarkdown_2.11         rstudioapi_0.13       
#> [55] iterators_1.0.14       R6_2.5.1               compiler_4.1.0
```
