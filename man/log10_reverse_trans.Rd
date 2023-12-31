% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ggplot.R
\name{log10_reverse_trans}
\alias{log10_reverse_trans}
\title{log10_reverse_trans}
\usage{
log10_reverse_trans()
}
\value{
A scales object
}
\description{
it perform log scaling and reverse the axis. Useful to plot negative log probabilities. To not be used directly but with ggplot (e.g. scale_y_continuous(trans = "log10_reverse") )
}
\details{
`r lifecycle::badge("maturing")`
}
\examples{

library(ggplot2)
library(tibble)

tibble(pvalue = c(0.001, 0.05, 0.1), fold_change = 1:3) \%>\%
 ggplot(aes(fold_change , pvalue)) +
 geom_point() +
 scale_y_continuous(trans = "log10_reverse")

}
