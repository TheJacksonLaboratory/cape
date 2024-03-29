% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/run_cape.R
\name{run_cape}
\alias{run_cape}
\title{Runs CAPE}
\usage{
run_cape(
  pheno_obj,
  geno_obj,
  results_file = "cross.RDS",
  p_or_q = 0.05,
  n_cores = 4,
  initialize_only = FALSE,
  verbose = TRUE,
  run_parallel = FALSE,
  param_file = NULL,
  yaml_params = NULL,
  results_path = NULL,
  plot_pdf = TRUE
)
}
\arguments{
\item{pheno_obj}{the cape object holding the phenotype data returned by}

\item{geno_obj}{the genotype object}

\item{results_file}{the name of the saved data_obj RDS file. The base 
name is used as the base name for all saved RDS files.}

\item{p_or_q}{A threshold indicating the maximum adjusted p value 
considered significant. Or, if FDR p value correction is used, the
the maximum q value considered significant.}

\item{n_cores}{The number of CPUs to use if run_parallel is set to TRUE}

\item{initialize_only, }{If TRUE, cape will not be run. Instead an initialized
data object will be returned. This data object will contain normalized and mean-centered
trait values, and eigentraits, and will have covariates specified. However, the 
singlescan, pairscan, etc. will not be run.}

\item{verbose}{Whether progress should be printed to the screen}

\item{run_parallel}{boolean, if TRUE runs certain parts of the code as parallel blocks}

\item{param_file}{yaml full path to the parameter file}

\item{yaml_params}{yaml string containing the parameters. Either the param_file or 
yaml_params can be null.}

\item{results_path}{path that results should be written to.}

\item{plot_pdf}{boolean, TRUE by default. If FALSE no pdf will be generated by the analysis.}
}
\value{
This function invisibly returns the data object with all final 
data included. In addition, data saved to the data_obj$results_path directory
}
\description{
This function takes in a data object and genotype object that
have been formatted for cape, as well as a string identifying
a parameter file. It runs cape on the data using the parameters
specified in the file.
}
\details{
This function assumes you already have all required libraries 
and functions loaded.
}
\examples{
\dontrun{
final_data_obj <- run_cape(pheno_obj, geno_obj)
}

}
