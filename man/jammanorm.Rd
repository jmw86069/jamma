% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/jamma-norm.R
\name{jammanorm}
\alias{jammanorm}
\title{Normalize data using MA-plot logic}
\usage{
jammanorm(
  x,
  controlGenes = NULL,
  minimum_mean = 0,
  controlSamples = NULL,
  centerGroups = NULL,
  useMean = TRUE,
  noise_floor = -Inf,
  noise_floor_value = noise_floor,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{x}{numeric matrix with expression data suitable for use
by \code{jammaplot()}. Gene expression data is typically transformed
using \code{log2(1+x)} to roughly normally distributed data.}
}
\value{
numeric matrix whose columns have been normalized using
the y-axis mean offset following MA-plot calculations. The output
contains attributes: \code{"nf"} numeric vector of normalization factors;
\code{"hk"} list of controlGenes used for each sample; \code{"hk_count"} the
number of controlGenes used for each sample.
}
\description{
Normalize data using MA-plot logic
}
\details{
This function normalizes data using the MA-plot y-axis
offset to define the normalization factor. It can optionally
use only a subset of \code{controlGenes}, such as known
housekeeper genes. It can optionally apply a \code{minimum_mean}
threshold, to ensure the mean abundance of each row
is at least this value.

The approach is very similar to other housekeeper normalization
methods, which typically use the geometric mean (or the mean
of the log-transformed abundances.) However, the strategy here
is to assert that the mean difference from average for the
selected \code{controlGenes} is equal to zero.

In order for this approach to be valid:
\itemize{
\item Data on the MA-plots should be horizontal for all samples.
\item The majority of controlGenes should not be changing, which
is seen as the controlGenes representing the mean expression.
\item For data with a small number of genes, known controlGenes
should be selected and verified to be unchanging across samples.
}
}
