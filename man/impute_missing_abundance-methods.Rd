% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R, R/methods_SE.R
\docType{methods}
\name{impute_missing_abundance}
\alias{impute_missing_abundance}
\alias{impute_missing_abundance,spec_tbl_df-method}
\alias{impute_missing_abundance,tbl_df-method}
\alias{impute_missing_abundance,tidybulk-method}
\alias{impute_missing_abundance,SummarizedExperiment-method}
\alias{impute_missing_abundance,RangedSummarizedExperiment-method}
\title{impute transcript abundance if missing from sample-transcript pairs}
\usage{
impute_missing_abundance(
  .data,
  .formula,
  .sample = NULL,
  .transcript = NULL,
  .abundance = NULL,
  suffix = "",
  force_scaling = FALSE
)

\S4method{impute_missing_abundance}{spec_tbl_df}(
  .data,
  .formula,
  .sample = NULL,
  .transcript = NULL,
  .abundance = NULL,
  suffix = "",
  force_scaling = FALSE
)

\S4method{impute_missing_abundance}{tbl_df}(
  .data,
  .formula,
  .sample = NULL,
  .transcript = NULL,
  .abundance = NULL,
  suffix = "",
  force_scaling = FALSE
)

\S4method{impute_missing_abundance}{tidybulk}(
  .data,
  .formula,
  .sample = NULL,
  .transcript = NULL,
  .abundance = NULL,
  suffix = "",
  force_scaling = FALSE
)

\S4method{impute_missing_abundance}{SummarizedExperiment}(
  .data,
  .formula,
  .sample = NULL,
  .transcript = NULL,
  .abundance = NULL,
  suffix = "",
  force_scaling = FALSE
)

\S4method{impute_missing_abundance}{RangedSummarizedExperiment}(
  .data,
  .formula,
  .sample = NULL,
  .transcript = NULL,
  .abundance = NULL,
  suffix = "",
  force_scaling = FALSE
)
}
\arguments{
\item{.data}{A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))}

\item{.formula}{A formula with no response variable, representing the desired linear model where the first covariate is the factor of interest and the second covariate is the unwanted variation (of the kind ~ factor_of_interest + batch)}

\item{.sample}{The name of the sample column}

\item{.transcript}{The name of the transcript/gene column}

\item{.abundance}{The name of the transcript/gene abundance column}

\item{suffix}{A character string. This is added to the imputed count column names. If empty the count column are overwritten}

\item{force_scaling}{A boolean. In case a abundance-containing column is not scaled (columns with _scale suffix), setting force_scaling = TRUE will result in a scaling by library size, to compensating for a possible difference in sequencing depth.}
}
\value{
A consistent object (to the input) non-sparse abundance

A consistent object (to the input) with imputed abundance

A consistent object (to the input) with imputed abundance

A consistent object (to the input) with imputed abundance

A `SummarizedExperiment` object

A `SummarizedExperiment` object
}
\description{
impute_missing_abundance() takes as input A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment)) and returns a consistent object (to the input) with additional sample-transcript pairs with imputed transcript abundance.
}
\details{
`r lifecycle::badge("maturing")`

This function imputes the abundance of missing sample-transcript pair using the median of the sample group defined by the formula
}
\examples{


res =
	impute_missing_abundance(
		tidybulk::se_mini,
	~ condition
)


}
