# MassExpression

Universal Imports + High Quality QC + Differential Expression Analysis = Awesomeness. 

## Download sample data

aws s3 sync s3://md-test-experiments/universal-input/data-formats/ . 

## RData available in the package for testing 

```{r}
data(package="MassExpression")
```

## Run example end-to-end

```{r }
outputFolder <- "output"
intensities <- fragpipe_data$intensities
design <- fragpipe_data$design

IntensityExperiment <- runGenericDiscovery(experimentDesign = design, proteinIntensities = intensities)
writeProteinViz(outputFolder = outputFolder, IntensityExperiment=IntensityExperiment)
```