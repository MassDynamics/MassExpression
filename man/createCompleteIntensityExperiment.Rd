% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createCompleteIntensityExperiment.R
\name{createCompleteIntensityExperiment}
\alias{createCompleteIntensityExperiment}
\title{Adds `longIntensityDT` information to `IntensityExperiment` with statistics about 
the number of replicates and number of imputed values in each condition of interest}
\usage{
createCompleteIntensityExperiment(
  IntensityExperiment,
  limmaStats,
  normalisationAppliedToAssay,
  longIntensityDT,
  conditionComparisonMapping
)
}
\arguments{
\item{IntensityExperiment}{output from `createSummarizedExperiment`}

\item{limmaStats}{Statistic table from limma output}

\item{normalisationAppliedToAssay}{one of 'None' or 'Median'}

\item{longIntensityDT}{Long format table with raw and imputed intensities. Each row is a feature (protein/peptide).}

\item{conditionComparisonMapping}{a dataframe matching different limma statistics comparison strings to up and down conditions}
}
\description{
Adds `longIntensityDT` information to `IntensityExperiment` with statistics about 
the number of replicates and number of imputed values in each condition of interest
}
