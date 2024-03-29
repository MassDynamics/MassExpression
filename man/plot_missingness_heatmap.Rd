% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qc_plotting_functions.R
\name{plot_missingness_heatmap}
\alias{plot_missingness_heatmap}
\title{Heatmap showing pattern of missingness}
\usage{
plot_missingness_heatmap(
  Experiment,
  assayName = "raw",
  condition_colname = "Condition",
  title = "Missingness pattern"
)
}
\arguments{
\item{Experiment}{SummarizedExperiment object}

\item{assayName}{name of assay to use}

\item{condition_colname}{str. Name of grouping condition.}

\item{title}{str. Plot title}
}
\description{
Columns (samples) have been clustered based on pattern of missingness.
}
