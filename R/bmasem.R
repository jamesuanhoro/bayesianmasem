#' Fit random-effects Bayesian meta-analytic CFAs with minor factors assumed.
#'
#' @description A function to fit fixed-, random-effects, and clustered
#' Bayesian meta-analytic CFAs with minor factors assumed.
#' Does not yet accomodate moderators,
#' and correlation matrices must be complete.
#' This will change in the near future.
#' @param model A description of the user-specified model, lavaan syntax.
#' @param data An optional data frame containing the observed variables used in
#' the model.
#' @param group An optional string identifying the grouping variable in
#' the data object.
#' @param sample_cov (list of matrices) sample covariance or correlation
#' matrices.
#' The rownames and/or colnames must contain the observed variable names.
#' For now, assumes there are no missing elements in the covariance matrices.
#' @param sample_nobs (vector of positive integer) Number of observations
#' for each study.
#' @param correlation (LOGICAL)
#' If TRUE (default): assume simple structure;
#' If FALSE: estimate all cross-loadings using generalized
#' @param method (character) One of "normal", "lasso", "logistic",
#' "GDP", or "none". See details below.
#' @param type (character) One of "fe", "re", or "dep" for fixed-effects,
#' random-effects, and dependent-samples MASEM respectively.
#' The "dep" argument is experimental, see details below.
#' @param orthogonal (LOGICAL)
#' If TRUE: constrain all factors orthogonal (overrides model syntax);
#' If FALSE: according to model syntax.
#' @param simple_struc (LOGICAL) Only relevant for CFAs.
#' If TRUE: analyze correlation matrices;
#' If FALSE: analyze covariance matrices.
#' @param cluster An optional integer vector identifying the cluster each group
#' belongs to.
#' Asssume there are five groups, the first three belong to cluster 1
#' and the last two belong to cluster 2,
#' then the argument would be: \code{cluster = c(1, 1, 1, 2, 2)}.
#' This feature is experimental, see details below.
#' @param seed (positive integer) seed, set to obtain replicable results.
#' @param warmup (positive integer) The number of warmup iterations to run per
#' chain.
#' @param sampling (positive integer) The number of post-warmup iterations to
#' run per chain, retained for inference.
#' @param refresh (positive integer) How often to print the status
#' of the sampler.
#' @param adapt_delta (real in (0, 1)) Increase to resolve divergent
#' transitions.
#' @param max_treedepth (positive integer) Increase to resolve problems with
#' maximum tree depth.
#' @param chains (positive integer) The number of Markov chains to run.
#' @param ncores (positive integer) The number of chains to run in parallel.
#' @param priors An object of \code{\link{bmasempriors-class}}.
#' See \code{\link{new_bmasempriors}} for more information.
#' @param show (Logical) If TRUE, show table of results, if FALSE, do not
#' show table of results. As an example, use FALSE for simulation studies.
#' @param show_messages (Logical) If TRUE, show messages from Stan sampler,
#' if FALSE, hide messages.
#' @param target (character) One of "rstan" or "cmdstan". If "cmdstan",
#' CmdStan and CmdStanR need to be installed on the device.
#' @returns An object of \code{\link{bmasem-class}}
#' @details
#' CFAs assume standardized factors.
#' Latent variable regression models are not yet implemented.
#'
#' There are different methods for estimating models in this package:
#'
#' - \code{normal}: under belief that minor factor influences are
#' on average zero with continuous deviations away from zero.
#' - \code{lasso}: under belief that minor factor influences are largely
#' zero with a small number of non-zero residual covariances.
#' - \code{logistic}: for similar belief as normal but more readily
#' accomodates extreme outliers.
#' - \code{GDP}: to mimic a global-local approach, i.e.
#' attempt to shrink near 0 residual covariances to 0
#' with minimal shrinking for larger residual covariances
#' \insertCite{armagan_generalized_2013}{bayesianmasem}.
#' - \code{none}: if intending to ignore the influence of minor factors.
#'
#' When \code{type = "dep"}, the user must supply the cluster IDs, see cluster
#' parameter documentation above. However, this feature is experimental and only
#' available when \code{target = "cmdstan"}.
#' Additionally, the cluster inputs are not validated.
#' @examples
#' \dontrun{
#' model_syntax <- "# latent variable definitions
#' F1 =~ JP1 + JP2 + JP3
#' F2 =~ JN1 + JN2 + JN4 + JN4
#' F3 =~ TD1 + TD2"
#' bmasem(
#'   model_syntax,
#'   sample_cov = issp89$data, sample_nobs = issp89$n
#' )
#' model_syntax <- paste0(
#'   "distress =~ ", paste0("x", 1:14, collapse = " + "), "\n",
#'   "anxiety =~ ", paste0("x", seq(1, 14, 2), collapse = " + "), "\n",
#'   "depression =~ ", paste0("x", seq(2, 14, 2), collapse = " + ")
#' )
#' bmasem(
#'   model_syntax,
#'   sample_cov = Norton13$data, sample_nobs = Norton13$n, orthogonal = TRUE
#' )
#' }
#' @references \insertAllCited{}
#' @export
bmasem <- function(
    model = NULL,
    data = NULL,
    group = NULL,
    sample_cov = NULL,
    sample_nobs = NULL,
    correlation = TRUE,
    method = "normal",
    type = "re",
    orthogonal = FALSE,
    simple_struc = TRUE,
    seed = 12345,
    warmup = 1000,
    sampling = 1000,
    refresh = (warmup + sampling) / 10,
    adapt_delta = .9,
    max_treedepth = 10,
    chains = 3,
    ncores = max(parallel::detectCores() - 2, 1),
    priors = new_bmasempriors(),
    show = TRUE,
    show_messages = TRUE,
    cluster = NULL,
    target = "rstan") {
  message("Processing user input ...")

  # Model cannot be NULL
  .user_input_check("model", model)

  # Priors must be class bmasempriors
  .user_input_check("priors", priors)

  # method must be valid
  .user_input_check("method", method)

  # type must be valid
  # Future: when moderators are added, warn user that fixed-effects
  # ignores moderators.
  .user_input_check("type", type)

  # Must provide either data and group or sample_cov and sample_nobs
  .user_input_check("data", data, group, sample_cov, sample_nobs)

  # target must be valid
  .user_input_check("target", target)

  # check for cluster when type = "dep"
  .user_input_check("cluster", type, cluster)

  # Run lavaan fit
  if (!is.null(data)) {
    lav_fit <- lavaan::cfa(
      model,
      data = data, group = group, std.lv = TRUE,
      likelihood = "wishart", ceq.simple = TRUE,
      do.fit = FALSE, orthogonal = orthogonal
    )
  } else {
    lav_fit <- lavaan::cfa(
      model,
      sample.cov = sample_cov, sample.nobs = sample_nobs, std.lv = TRUE,
      likelihood = "wishart", ceq.simple = TRUE,
      do.fit = FALSE, orthogonal = orthogonal
    )
  }

  par_table <- lavaan::lavaanify(
    model,
    ceq.simple = TRUE, std.lv = TRUE
  )

  # Obtain data list for Stan
  data_list <- .create_data_list_meta(
    lavaan_object = lav_fit,
    method = method,
    type = type,
    simple_struc = simple_struc,
    priors = priors,
    cluster = cluster,
    correlation = correlation,
    partab = par_table
  )

  message("User input fully processed :)\n Now to modeling.")

  stan_fit <- .target_fitter(
    target, data_list, seed, warmup, sampling, refresh,
    adapt_delta, max_treedepth, chains, ncores, show_messages
  )

  bma_results <- .clean_up_stan_fit(
    stan_fit = stan_fit, data_list = data_list, priors = priors
  )
  if (show) {
    show(bma_results)
  }

  return(bma_results)
}
