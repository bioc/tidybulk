% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dplyr_methods.R
\name{rowwise}
\alias{rowwise}
\title{Group input by rows}
\arguments{
\item{data}{Input data frame.}

\item{...}{Variables to be preserved when calling summarise(). This is typically a set of variables whose combination uniquely identify each row. NB: unlike group_by() you can not create new variables here but instead you can select multiple variables with (e.g.) everything().}
}
\value{
A consistent object (to the input)

  A `tbl`
}
\description{
See [this repository](https://github.com/jennybc/row-oriented-workflows)
for alternative ways to perform row-wise operations.
}
\details{
`rowwise()` is used for the results of [do()] when you
create list-variables. It is also useful to support arbitrary
complex operations that need to be applied to each row.

Currently, rowwise grouping only works with data frames. Its
main impact is to allow you to work with list-variables in
[summarise()] and [mutate()] without having to
use \code{[[1]]}. This makes `summarise()` on a rowwise tbl
effectively equivalent to [plyr::ldply()].
}
\examples{

df <- expand.grid(x = 1:3, y = 3:1)
df_done <- df |> rowwise() 

}
