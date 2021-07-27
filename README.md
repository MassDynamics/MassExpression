MassExpression
================

-   [Download sample data](#download-sample-data)
-   [RData available in the package for
    testing](#rdata-available-in-the-package-for-testing)
-   [Run example end-to-end](#run-example-end-to-end)
    -   [Render QC](#render-qc)
    -   [Save artefacts](#save-artefacts)

Universal Imports + High Quality QC + Differential Expression Analysis =
Awesomeness.

## Download sample data

aws s3 sync s3://md-test-experiments/universal-input/data-formats/ .

## RData available in the package for testing

``` r
data(package="MassExpression")
```

# Run example end-to-end

Load data and run workflow with the runner `runGenericDiscovery`.

``` r
library(MassExpression)

output_folder <- "path/to/output"

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
    listComparisonExperiments(completeExperiment)
  
saveOutput(IntensityExperiment = IntensityExperiment, 
CompleteIntensityExperiment = CompleteExperiment, output_folder =  output_folder)
```

## Render QC

``` r
# Render and save QC report 
qc_report <- system.file("rmd","QC_report.Rmd", package = "MassExpression")

rmarkdown::render(qc_report,
                  params = list(listInt = listIntensityExperiments,
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
                  params = list(listInt = listIntensityExperiments,
                                experiment = "Mass Dynamics QC report",
                                output_figure = file.path(output_folder_pdf, "figure_pdf/"),
                                format = "pdf"),
                  output_file = file.path(output_folder_pdf, "QC_Report.pdf"),
                  output_format=rmarkdown::pdf_document(
                    toc = TRUE,
                    fig_caption= TRUE))
```

`results` is a list containing two `SummarizedExperiment` objects:

-   `IntensityExperiment`: contains the raw data (including missing
    values)

-   `CompleteIntensityExperiment`: contains the imputed data and summary
    statistics about the number of replicates and imputed proteins in
    each group of the conditions of interest.

## Save artefacts

The DE results from `IntensityExperiment` are going to be displayed for
a user and therefore they need to be parsed into a `json` output (using
the `writeProteinViz` function).

``` r
CompleteIntensityExperiment <- listIntensityExperiments$CompleteIntensityExperiment
results <- listIntensityExperiments$IntensityExperiment
writeProteinViz(outputFolder = outputFolder, IntensityExperiment=results$IntensityExperiment)
```
