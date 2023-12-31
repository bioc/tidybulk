% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tidyr_methods.R
\name{unnest}
\alias{unnest}
\alias{nest}
\title{unnest}
\arguments{
\item{data}{A tbl. (See tidyr)}

\item{cols}{<[`tidy-select`][tidyr_tidy_select]> Columns to unnest.
If you `unnest()` multiple columns, parallel entries must be of
compatibble sizes, i.e. they're either equal or length 1 (following the
standard tidyverse recycling rules).}

\item{names_sep}{If `NULL`, the default, the names will be left
  as is. In `nest()`, inner names will come from the former outer names;
  in `unnest()`, the new outer names will come from the inner names.

  If a string, the inner and outer names will be used together. In `nest()`,
  the names of the new outer columns will be formed by pasting together the
  outer and the inner column names, separated by `names_sep`. In `unnest()`,
  the new inner names will have the outer names (+ `names_sep`) automatically
  stripped. This makes `names_sep` roughly symmetric between nesting and unnesting.}

\item{keep_empty}{See tidyr::unnest}

\item{names_repair}{See tidyr::unnest}

\item{ptype}{See tidyr::unnest}

\item{.drop}{See tidyr::unnest}

\item{.id}{tidyr::unnest}

\item{.sep}{tidyr::unnest}

\item{.preserve}{See tidyr::unnest}

\item{.data}{A tbl. (See tidyr)}

\item{...}{Name-variable pairs of the form new_col = c(col1, col2, col3) (See tidyr)}
}
\value{
A tidySummarizedExperiment objector a tibble depending on input

A tt object
}
\description{
unnest

nest
}
\examples{


tidybulk::se_mini |> tidybulk() |> nest( data = -.feature) |> unnest(data)


tidybulk::se_mini \%>\% tidybulk() \%>\% nest( data = -.feature)

}
