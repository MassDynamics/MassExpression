# MassExpression

Universal Imports + High Quality QC + Differential Expression Analysis = Awesomeness. 

## Download sample data

aws s3 sync s3://md-test-experiments/universal-input/data-formats/ . 

## RData available in the package for testing 

```{r}
data(package="MassExpression")
```

# Run example end-to-end

Load data and run workflow with the runner `runGenericDiscovery`.

```{r fragpipe}
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

Render QC

```{r render-qc}
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
  - `IntensityExperiment`: contains the raw data (including missing values)
  - `CompleteIntensityExperiment`: contains the imputed data and summary statistics about the number of replicates and imputed proteins in each group of the conditions of interest. 


# Outputs

The MassExpression package should produce outputs in several different formats, each serving a different purpose:

1. Results Tables (these are standard human readable tables with parsed/transformed data, processed data and statistics)
2. JSON Objects. (these are designed as inputs to a RAILS database for internal MD use)
    a. protein_viz.json (this is the data required for the volcano plot analysis feature. It is the LIMMA statistics for each binary comparison).
    b. protein_counts_and_intensities.json (this is the protein intensities in long form broken down by each protein/experimental condition with associate statistics and imputation flag.)
3. .RData object (this is a catchall for experiment data we want to )

**Save artefacts**

The DE results from `IntensityExperiment` are going to be displayed for a user and therefore they need to be parsed into a `json` output (using the `writeProteinViz` function).  

```{r}
CompleteIntensityExperiment <- listIntensityExperiments$CompleteIntensityExperiment
results <- listIntensityExperiments$IntensityExperiment
writeProteinViz(outputFolder = outputFolder, IntensityExperiment=results$IntensityExperiment)
```

**Schema of the protein_counts_and_intensities.json output**

```
[
	{
		"ProteinGroupId": "100",
		"ProteinId": "O35593",
		"GeneName": "Psmd14",
		"ProteinDescription": "26S proteasome non-ATPase regulatory subunit 14",
		"FastaHeaders": "sp|O35593|PSDE_MOUSE 26S proteasome non-ATPase regulatory subunit 14 OS=Mus musculus OX=10090 GN=Psmd14 PE=1 SV=2",
		"ProteinQValue": 0,
		"ProteinScore": 10.73,
		"conditions": [
			{
				"name": "Cerebellum",
				"precentageOfReplicates": 0.5,
				"numberOfReplicateCount": 4,
				"intensityValues": [
					{
						"replicateNum": 1,
						"centeredIntensity": 0.749918461968428,
						"z_norm": 0.894384217678135,
						"log2NInt_ProteinGroupId": -2.96634794751,
						"Imputed": 0
					},
					{
						"replicateNum": 2,
						"centeredIntensity": 0.647910901219001,
						"z_norm": 0.772725721394879,
						"log2NInt_ProteinGroupId": -3.07935607412021,
						"Imputed": 0
					},
			}
	}
]

```