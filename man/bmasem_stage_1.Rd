% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bmasem_stage_1.R
\name{bmasem_stage_1}
\alias{bmasem_stage_1}
\title{Fit random-effects Bayesian meta-analytic CFAs with minor factors assumed.}
\usage{
bmasem_stage_1(
  sample_cov = NULL,
  sample_nobs = NULL,
  type = "re",
  seed = 12345,
  warmup = 1000,
  sampling = 1000,
  refresh = (warmup + sampling)/10,
  adapt_delta = 0.9,
  max_treedepth = 10,
  chains = 3,
  ncores = max(parallel::detectCores() - 2, 1),
  priors = new_bmasempriors(),
  show = TRUE,
  show_messages = TRUE,
  cluster = NULL,
  conditional_re = TRUE
)
}
\arguments{
\item{sample_cov}{(list of matrices) sample covariance or correlation
matrices.
The rownames and/or colnames must contain the observed variable names.
For now, assumes there are no missing elements in the covariance matrices.}

\item{sample_nobs}{(vector of positive integer) Number of observations
for each study.}

\item{type}{(character) One of "fe", "re", or "dep" for fixed-effects,
random-effects, and dependent-samples MASEM respectively.
The "dep" argument is experimental, see details below.}

\item{seed}{(positive integer) seed, set to obtain replicable results.}

\item{warmup}{(positive integer) The number of warmup iterations to run per
chain.}

\item{sampling}{(positive integer) The number of post-warmup iterations to
run per chain, retained for inference.}

\item{refresh}{(positive integer) How often to print the status
of the sampler.}

\item{adapt_delta}{(real in (0, 1)) Increase to resolve divergent
transitions.}

\item{max_treedepth}{(positive integer) Increase to resolve problems with
maximum tree depth.}

\item{chains}{(positive integer) The number of Markov chains to run.}

\item{ncores}{(positive integer) The number of chains to run in parallel.}

\item{priors}{An object of \code{\link{bmasempriors-class}}.
See \code{\link{new_bmasempriors}} for more information.}

\item{show}{(Logical) If TRUE, show table of results, if FALSE, do not
show table of results. As an example, use FALSE for simulation studies.}

\item{show_messages}{(Logical) If TRUE, show messages from Stan sampler,
if FALSE, hide messages.}

\item{cluster}{An optional integer vector identifying the cluster each group
belongs to.
Asssume there are five groups, the first three belong to cluster 1
and the last two belong to cluster 2,
then the argument would be: \code{cluster = c(1, 1, 1, 2, 2)}.
This feature is experimental, see details below.}

\item{conditional_re}{(LOGICAL)
Only relevant for analysis of correlation structures.
If TRUE, sample levels of the study-level random effect;
If FALSE, don't.}
}
\value{
A list containing fit indices, pooled correlation matrix and
its asymptotic covariance matrix, pooled partial correlation matrix and
its asymptotic covariance matrix, Stan object and data_list used to fit
Stan object.
}
\description{
A function to pool correlation matrices
permitting fixed-, random-effects, and clustered-samples pooling.
Correlation matrices must be complete. This will change in the near future.
}
\details{
When \code{type = "dep"}, the user must supply the cluster IDs, see cluster
parameter documentation above. However, this feature is experimental.
Additionally, the cluster inputs are not validated.
}
\examples{
\dontrun{
pool_fit <- bmasem_stage_1(
  sample_cov = issp89$data, sample_nobs = issp89$n
)
# Pooled correlation matrix
pool_fit$r_mat
pool_fit$r_mat_cov
Pooled partial correlation matrix
pool_fit$p_mat
pool_fit$p_mat_cov
}
}
\references{
\insertAllCited{}
}
