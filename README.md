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
devtools::install_github("MassExpression")
```

``` r
library(MassExpression)
utils::packageVersion("MassExpression")
#> [1] '0.0.68'
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
CompleteIntensityExperiment = CompleteExperiment,
longIntensityDT = longIntensityDT, 
output_folder =  output_folder)
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
#>  [1] MassExpression_0.0.68       SummarizedExperiment_1.22.0
#>  [3] GenomicRanges_1.44.0        GenomeInfoDb_1.28.1        
#>  [5] IRanges_2.26.0              S4Vectors_0.30.0           
#>  [7] MatrixGenerics_1.4.2        matrixStats_0.60.1         
#>  [9] Biobase_2.52.0              BiocGenerics_0.38.0        
#> 
#> loaded via a namespace (and not attached):
#>  [1] tidyselect_1.1.1       xfun_0.25              purrr_0.3.4           
#>  [4] lattice_0.20-44        colorspace_2.0-2       vctrs_0.3.8           
#>  [7] generics_0.1.0         htmltools_0.5.1.1      yaml_2.2.1            
#> [10] utf8_1.2.2             rlang_0.4.11           pillar_1.6.2          
#> [13] glue_1.4.2             DBI_1.1.1              RColorBrewer_1.1-2    
#> [16] GenomeInfoDbData_1.2.6 foreach_1.5.1          lifecycle_1.0.0       
#> [19] stringr_1.4.0          zlibbioc_1.38.0        munsell_0.5.0         
#> [22] gtable_0.3.0           codetools_0.2-18       evaluate_0.14         
#> [25] knitr_1.33             fansi_0.5.0            scales_1.1.1          
#> [28] limma_3.48.3           DelayedArray_0.18.0    jsonlite_1.7.2        
#> [31] XVector_0.32.0         ggplot2_3.3.5          digest_0.6.27         
#> [34] stringi_1.7.3          dplyr_1.0.7            grid_4.1.0            
#> [37] tools_4.1.0            bitops_1.0-7           magrittr_2.0.1        
#> [40] RCurl_1.98-1.4         tibble_3.1.3           tidyr_1.1.3           
#> [43] crayon_1.4.1           pkgconfig_2.0.3        ellipsis_0.3.2        
#> [46] pheatmap_1.0.12        Matrix_1.3-4           data.table_1.14.0     
#> [49] assertthat_0.2.1       rmarkdown_2.10         iterators_1.0.13      
#> [52] R6_2.5.1               compiler_4.1.0
```
