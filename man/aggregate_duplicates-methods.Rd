% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R, R/methods_SE.R
\docType{methods}
\name{aggregate_duplicates}
\alias{aggregate_duplicates}
\alias{aggregate_duplicates,spec_tbl_df-method}
\alias{aggregate_duplicates,tbl_df-method}
\alias{aggregate_duplicates,tidybulk-method}
\alias{aggregate_duplicates,SummarizedExperiment-method}
\alias{aggregate_duplicates,RangedSummarizedExperiment-method}
\title{Aggregates multiple counts from the same samples (e.g., from isoforms), concatenates other character columns, and averages other numeric columns}
\usage{
aggregate_duplicates(
  .data,
  .sample = NULL,
  .transcript = NULL,
  .abundance = NULL,
  aggregation_function = sum,
  keep_integer = TRUE
)

\S4method{aggregate_duplicates}{spec_tbl_df}(
  .data,
  .sample = NULL,
  .transcript = NULL,
  .abundance = NULL,
  aggregation_function = sum,
  keep_integer = TRUE
)

\S4method{aggregate_duplicates}{tbl_df}(
  .data,
  .sample = NULL,
  .transcript = NULL,
  .abundance = NULL,
  aggregation_function = sum,
  keep_integer = TRUE
)

\S4method{aggregate_duplicates}{tidybulk}(
  .data,
  .sample = NULL,
  .transcript = NULL,
  .abundance = NULL,
  aggregation_function = sum,
  keep_integer = TRUE
)

\S4method{aggregate_duplicates}{SummarizedExperiment}(
  .data,
  .sample = NULL,
  .transcript = NULL,
  .abundance = NULL,
  aggregation_function = sum,
  keep_integer = TRUE
)

\S4method{aggregate_duplicates}{RangedSummarizedExperiment}(
  .data,
  .sample = NULL,
  .transcript = NULL,
  .abundance = NULL,
  aggregation_function = sum,
  keep_integer = TRUE
)
}
\arguments{
\item{.data}{A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment))}

\item{.sample}{The name of the sample column}

\item{.transcript}{The name of the transcript/gene column}

\item{.abundance}{The name of the transcript/gene abundance column}

\item{aggregation_function}{A function for counts aggregation (e.g., sum,  median, or mean)}

\item{keep_integer}{A boolean. Whether to force the aggregated counts to integer}
}
\value{
A consistent object (to the input) with aggregated transcript abundance and annotation

A consistent object (to the input) with aggregated transcript abundance and annotation

A consistent object (to the input) with aggregated transcript abundance and annotation

A consistent object (to the input) with aggregated transcript abundance and annotation

A `SummarizedExperiment` object

A `SummarizedExperiment` object
}
\description{
aggregate_duplicates() takes as input A `tbl` (with at least three columns for sample, feature and transcript abundance) or `SummarizedExperiment` (more convenient if abstracted to tibble with library(tidySummarizedExperiment)) and returns a consistent object (to the input) with aggregated transcripts that were duplicated.
}
\details{
`r lifecycle::badge("maturing")`

This function aggregates duplicated transcripts (e.g., isoforms, ensembl).
For example, we often have to convert ensembl symbols to gene/transcript symbol,
 but in doing so we have to deal with duplicates. `aggregate_duplicates` takes a tibble
 and column names (as symbols; for `sample`, `transcript` and `count`) as arguments and
 returns a tibble with aggregate transcript with the same name. All the rest of the column
 are appended, and factors and boolean are appended as characters.

 Underlying custom method:
 data |>
		filter(n_aggr > 1) |>
		group_by(!!.sample,!!.transcript) |>
		dplyr::mutate(!!.abundance := !!.abundance |> aggregation_function())
}
\examples{

# Create a aggregation column
se_mini = tidybulk::se_mini
SummarizedExperiment::rowData(se_mini )$gene_name = rownames(se_mini )

   aggregate_duplicates(
     se_mini,
   .transcript = gene_name
   )


}
