MassExpression
================

The LFQProcessing R package enables users to load and process a simple generic format for 
proteomic LC-MS/MS shotgun experiment data. This package contains utilities for performing differential expression DE statistics and Quality Control reports.

A sister R package LFQProcessing, and precedessor, provides the same utilities for output from the [MaxQuant Computational Platform](https://www.maxquant.org/).

Both LFQProcessing and MassExpression were developed by MassDynamics to enable greater accessibility, reproducibility and insight into proteomics datasets and insights. All feedback is welcome.

For more details please see our biorXiv paper: [Mass Dynamics 1.0: A streamlined, web-based environment for analyzing, sharing and integrating Label-Free Data](https://doi.org/10.1101/2021.03.03.433806).

-   [RData available in the package for
    testing](#rdata-available-in-the-package-for-testing)
-   [Run example end-to-end](#run-example-end-to-end)
    -   [Render QC](#render-qc)
    -   [Special Artifacts Export](#special-artifacts-export)

## RData available in the package for testing

``` r
data(package="MassExpression")
```

# Run example end-to-end

Load data and run workflow with the runner `runGenericDiscovery`.

``` r
library(MassExpression)

output_folder <- "./test_output"

design <- fragpipe_data$design
intensities <- fragpipe_data$intensities

intensities <- fragpipe_data$intensities
design <- fragpipe_data$design
parameters <- fragpipe_data$parameters
normalisation_method <- parameters[parameters[,1] == "UseNormalisationMethod",2]
species <- parameters[parameters[,1] == "Species",2]
labellingMethod <- parameters[parameters[,1] == "LabellingMethod",2]


results <- runGenericDiscovery(experimentDesign = design, 
                               proteinIntensities = intensities, 
                               normalisationMethod = normalisation_method, 
                               species = species, 
                               labellingMethod = labellingMethod)

CompleteExperiment <- results$CompleteIntensityExperiment
IntensityExperiment <- results$IntensityExperiment

comparisonExperiments <- 
    listComparisonExperiments(CompleteExperiment)
  
saveOutput(IntensityExperiment = IntensityExperiment, 
CompleteIntensityExperiment = CompleteExperiment, output_folder =  output_folder)
```

`results` is a list containing two `SummarizedExperiment` objects:

-   `IntensityExperiment`: contains the raw data (including missing
    values)

-   `CompleteIntensityExperiment`: contains the imputed data and summary
    statistics about the number of replicates and imputed proteins in
    each group of the conditions of interest.
    
## Render QC

``` r
# Render and save QC report 
qc_report <- system.file("rmd","QC_report.Rmd", package = "MassExpression")
listIntensityExperiments = list(IntensityExperiment=IntensityExperiment,
                                CompleteIntensityExperiment=CompleteExperiment)

qc_report_folder = file.path(getwd(), output_folder)
rmarkdown::render(qc_report,
                  params = list(listInt = listIntensityExperiments,
                                experiment = "Mass Dynamics QC report",
                                output_figure = file.path(qc_report_folder, "figure_html/"),
                                format = "html"),
                  output_dir = qc_report_folder,
                  output_format=rmarkdown::html_document(
                    self_contained=FALSE,
                    lib_dir=file.path(qc_report_folder,"qc_report_files"),
                    code_folding= "hide",
                    theme="united",
                    toc_float = TRUE,
                    fig_caption= TRUE,
                    df_print="paged"))
# Render PDF
rmarkdown::render(qc_report,
                  params = list(listInt = listIntensityExperiments,
                                experiment = "Mass Dynamics QC report",
                                output_figure = file.path(qc_report_folder, "figure_pdf/"),
                                format = "pdf"),
                  output_dir = qc_report_folder,
                  output_format=rmarkdown::pdf_document(
                    toc = TRUE,
                    fig_caption= TRUE))
```


## Special Artifacts Export

The DE results from `IntensityExperiment` can be parsed into a `json` output (using
the `writeProteinViz` function).

``` r
CompleteIntensityExperiment <- listIntensityExperiments$CompleteIntensityExperiment
results <- listIntensityExperiments$IntensityExperiment
writeProteinViz(outputFolder = outputFolder, IntensityExperiment=results$IntensityExperiment)
```
