% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/error_prop.R
\name{error_prop}
\alias{error_prop}
\title{Estimate Errors of Regression Coefficients}
\usage{
error_prop(
  data_obj,
  pairscan_obj,
  perm = FALSE,
  verbose = FALSE,
  run_parallel = FALSE,
  n_cores = 4,
  just_m = FALSE
)
}
\arguments{
\item{data_obj}{a \code{\link{Cape}} object}

\item{pairscan_obj}{a pairscan object from \code{\link{pairscan}}}

\item{perm}{A logical value to indicate whether error propagation should be 
performed on the test statistics (FALSE) or the permuted test statistics (TRUE).}

\item{verbose}{A logical value to indicate whether the progress of the function 
should be printed to the screen.}

\item{run_parallel}{boolean, default = FALSE}

\item{n_cores}{The number of cores to use if run_parallel is TRUE, default = 4}

\item{just_m}{If TRUE only the m12 and m21 values are calculated. If FALSE, the
default, the standard deviations are also calculated.}
}
\value{
This function returns the data object with a new list element: var_to_var_influences 
if perm is set to FALSE and var_to_var_influences_perm if perm is set to TRUE. These tables 
include the errors calculated for the marker1 to marker2 (m21) influences as well as the
marker2 to marker1 (m12) influences. These results are used by \code{\link{calc_p}} to
calculate empirical p values.
}
\description{
This function uses error propagation formulas for 
quantities computed from regression coefficients to 
estimate the error for all regression coefficients.
}
