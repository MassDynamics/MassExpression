# MassExpression

Universal Imports + High Quality QC + Differential Expression Analysis = Awesomeness. 

## Download sample data

aws s3 sync s3://md-test-experiments/universal-input/data-formats/ . 

## RData available in the package for testing 

```{r}
data(package="MassExpression")
```

# Run example end-to-end

```{r fragpipe}
library(MassExpression)

output_folder <- "path/to/output"

design <- fragpipe_data$design
intensities <- fragpipe_data$intensities

# The flag should eventually come from the metadata 
normalise_flag <- FALSE

# Workflow runner
# If save = FALSE the output is printed into memory 
listIntensityExperiments = runGenericDiscovery(experimentDesign = design, 
                    proteinIntensities = intensities,
                    normalise = normalise_flag)
                    
# Save output to folder
Intensity <- listIntensityExperiments$IntensityExperiment
CompleteIntensity <-  listIntensityExperiments$CompleteIntensityExperiment
saveOutput(Intensity, CompleteIntensity, output_folder)

# Render and save QC report 
qc_report <- system.file("rmd","QC_report.Rmd", package = "MassExpression")
rmarkdown::render(qc_report, 
                  params = list(listInt = listIntensityExperiments,
                                experiment = "Mass Dynamics QC report",
                                output_figure = output_folder),
                                output_dir = output_folder,
                  output_file = "test-data-qc.html",
                  output_format=rmarkdown::html_document(
                            self_contained=FALSE,
                            lib_dir=file.path(output_folder,"qc_report_files")))
```




`listIntensityExperiments` is a list containing two `SummarizedExperiment` objects:
  - `IntensityExperiment`: contains the raw data (including missing values) and the results of the limma statistics (which can be accessed with `rowData(IntensityExperiment)`) 
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