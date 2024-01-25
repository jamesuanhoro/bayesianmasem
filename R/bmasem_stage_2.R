#' Fit random-effects Bayesian meta-analytic CFAs with minor factors assumed.
#'
#' @description A function that takes a pooled correlation matrix object
#' and returns a fitted SEM.
#' @param pool_fit The list returned by the \code{bmasem_stage_1} function.
#' @inheritParams bmasem
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
#' @examples
#' \dontrun{
#' pool_fit_issp <- bmasem_stage_1(
#'   sample_cov = issp89$data, sample_nobs = issp89$n
#' )
#' model_syntax_issp <- "# latent variable definitions
#' F1 =~ JP1 + JP2 + JP3
#' F2 =~ JN1 + JN2 + JN4 + JN4
#' F3 =~ TD1 + TD2"
#' bmasem_stage_2(model_syntax_issp, pool_fit_issp)
#' pool_fit_norton <- bmasem_stage_1(
#'   sample_cov = Norton13$data, sample_nobs = Norton13$n
#' )
#' model_syntax_norton <- paste0(
#'   "distress =~ ", paste0("x", 1:14, collapse = " + "), "\n",
#'   "anxiety =~ ", paste0("x", seq(1, 14, 2), collapse = " + "), "\n",
#'   "depression =~ ", paste0("x", seq(2, 14, 2), collapse = " + ")
#' )
#' bmasem_stage_2(model_syntax_norton, pool_fit_norton, orthogonal = TRUE)
#' }
#' @references \insertAllCited{}
#' @export
bmasem_stage_2 <- function(
    model = NULL,
    pool_fit = NULL,
    method = "normal",
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
    show_messages = TRUE) {
  message("Processing user input ...")

  # Model cannot be NULL
  .user_input_check("model", model)

  # Priors must be class bmasempriors
  .user_input_check("priors", priors)

  # method must be valid
  .user_input_check("method", method)

  # Run lavaan fit
  lav_fit <- lavaan::cfa(
    model,
    sample.cov = pool_fit$r_mat,
    sample.nobs = sum(pool_fit$data_list$Np),
    std.lv = TRUE,
    likelihood = "wishart", ceq.simple = TRUE,
    do.fit = FALSE, orthogonal = orthogonal
  )

  par_table <- lavaan::lavaanify(
    model,
    ceq.simple = TRUE, std.lv = TRUE
  )

  # Obtain data list for Stan
  data_list <- .create_data_list_pooled(
    lavaan_object = lav_fit,
    method = method,
    simple_struc = simple_struc,
    priors = priors,
    partab = par_table,
    acov_mat = pool_fit$r_mat_cov,
    old_names = rownames(pool_fit$r_mat)
  )

  message("User input fully processed :)\n Now to modeling.")

  cfa_model <- instantiate::stan_package_model(
    name = "cfa_cor", package = "bayesianmasem"
  )

  message("Fitting Stan model ...")

  stan_fit <- cfa_model$sample(
    data = data_list,
    seed = seed,
    iter_warmup = warmup,
    iter_sampling = sampling,
    refresh = refresh,
    init = .init_fx(data_list),
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    chains = chains,
    parallel_chains = ncores,
    show_messages = show_messages
  )

  bma_results <- .clean_up_stan_fit(
    stan_fit = stan_fit, data_list = data_list, priors = priors
  )
  if (show) {
    show(bma_results)
  }

  return(bma_results)
}
