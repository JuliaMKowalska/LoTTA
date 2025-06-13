% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_model.R
\name{Initial_treatment_DIS}
\alias{Initial_treatment_DIS}
\title{function that samples initial values for the treatment model with a discrete prior}
\usage{
Initial_treatment_DIS(x, t, Ct_start, cstart, grid, lb, ubr, ubl, s, jlb = 0.2)
}
\arguments{
\item{x}{\itemize{
\item score data
}}

\item{t}{\itemize{
\item treatment data
}}

\item{Ct_start}{\itemize{
\item posterior samples of cutoff location (categorized by natural numbers)
obtained through "cutoff_initial_dis.txt"c
}}

\item{cstart}{\itemize{
\item the first point with a positive prior mass
}}

\item{grid}{\itemize{
\item distance between two consecutive points with nonzero prior mass
}}

\item{lb}{\itemize{
\item minimum window size (grid size in case of discrete score)
}}

\item{ubr}{\itemize{
\item maximum value of the window's right boundary point
}}

\item{ubl}{\itemize{
\item minimum value of the window's left boundary point
}}

\item{s}{\itemize{
\item seed
}}

\item{jlb}{\itemize{
\item minimum jump size
}}
}
\value{
list with initial parameters values for treatment model with discrete score
and .RNG.seed value
}
\description{
function that samples initial values for the treatment model with a discrete prior
}
