% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_model.R
\name{bounds}
\alias{bounds}
\title{function that finds maximum widow size to searxh for a cutoff}
\usage{
bounds(x, ns = 25)
}
\arguments{
\item{x}{\itemize{
\item score data
}}

\item{ns}{\itemize{
\item minimum number of data points on each side of the cutoff to which cubic
parts are fitted
}}
}
\value{
list with ubl - minimum value of the window's left boundary point, ubr - maximum value of the window's right boundary point
}
\description{
function that finds maximum widow size to searxh for a cutoff
}
