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
outputFolder <- "../../output"
design <- fragpipe_data$design
intensities <- fragpipe_data$intensities

listIntensityExperiments <- runGenericDiscovery(experimentDesign = design, proteinIntensities = intensities)
```

`listIntensityExperiments` is a list containing two `SummarizedExperiment` objects:
  - `IntensityExperiment`: contains the raw data (including missing values) and the results of the limma statistics (which can be accessed with `rowData(IntensityExperiment)`) 
  - `CompleteIntensityExperiment`: contains the imputed data and summary statistics about the number of replicates and imputed proteins in each group of the conditions of interest. 

**Save artefacts**

The DE results from `IntensityExperiment` are going to be displayed for a user and therefore they need to be parsed into a `json` output (using the `writeProteinViz` function).  

```{r}
CompleteIntensityExperiment <- listIntensityExperiments$CompleteIntensityExperiment
IntensityExperiment <- listIntensityExperiments$IntensityExperiment
writeProteinViz(outputFolder = outputFolder, IntensityExperiment=IntensityExperiment)
```

**Schema of the json output**

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