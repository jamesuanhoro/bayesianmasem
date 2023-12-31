% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/00_class_defs.R
\docType{class}
\name{bmasem-class}
\alias{bmasem-class}
\title{A class for models fitted with the package.}
\description{
A class for models fitted with the package.
}
\section{Slots}{

\describe{
\item{\code{major_parameters}}{Summary statistics for structural parameters.}

\item{\code{minor_factor_matrix}}{Summary statistics for standardized
residual covariances.}

\item{\code{data_list}}{Data used to fit the model.}

\item{\code{priors}}{Priors used to fit the model.}

\item{\code{stan_fit}}{Fitted CmdStan model.}

\item{\code{version}}{Package version used to fit model.}
}}

