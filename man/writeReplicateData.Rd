% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/saveOutput.R
\name{writeReplicateData}
\alias{writeReplicateData}
\title{Write data needed for the Replicates tab
This is a wrapper of `oneProteinReplData`, looping over all proteins.}
\usage{
writeReplicateData(longDTProt, outputFolder)
}
\arguments{
\item{longDTProt}{data.table protein information stored in long format. 
Columns required: `ProteinId`, `GeneName`, `Description`, `log2NInt`, `Condition`,
 `Replicate`, `Imputed`.}

\item{outputFolder}{str. Path to folder where `data` should be saved.}
}
\description{
Write data needed for the Replicates tab
This is a wrapper of `oneProteinReplData`, looping over all proteins.
}
