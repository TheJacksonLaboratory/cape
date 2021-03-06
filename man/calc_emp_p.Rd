% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_emp_p.R
\name{calc_emp_p}
\alias{calc_emp_p}
\title{Calculate empirical p-values}
\usage{
calc_emp_p(obs_dist, null_dist)
}
\arguments{
\item{obs_dist}{The observed distribution}

\item{null_dist}{The null distribution}
}
\value{
An empirical p value for each observed value
}
\description{
This function uses ecdf to calculate empirical p values
given a null distribution and an observed distribution
}
