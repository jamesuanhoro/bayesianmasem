#' Fit random-effects Bayesian meta-analytic CFAs with minor factors assumed.
#'
#' @description A function to pool correlation matrices
#' permitting fixed-, random-effects, and clustered-samples pooling.
#' Correlation matrices must be complete. This will change in the near future.
#' @inheritParams bmasem
#' @returns A list containing fit indices, pooled correlation matrix,
#' asymptotic covariance matrix, Stan object and
#' data_list used to fit Stan object
#' @details
#' When \code{type = "dep"}, the user must supply the cluster IDs, see cluster
#' parameter documentation above. However, this feature is experimental.
#' Additionally, the cluster inputs are not validated.
#' @examples
#' \dontrun{
#' pool_fit <- bmasem_stage_1(
#'   sample_cov = issp89$data, sample_nobs = issp89$n
#' )
#' pool_fit$r_mat
#' pool_fit$r_mat_cov
#' }
#' @references \insertAllCited{}
#' @export
bmasem_stage_1 <- function(
    sample_cov = NULL,
    sample_nobs = NULL,
    type = "re",
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
    cluster = NULL) {
  message("Processing user input ...")

  # Priors must be class bmasempriors
  .user_input_check("priors", priors)

  # type must be valid
  # Future: when moderators are added, warn user that fixed-effects
  # ignores moderators.
  .user_input_check("type", type)

  # Must provide either data and group or sample_cov and sample_nobs
  .user_input_check("data", NULL, NULL, sample_cov, sample_nobs)

  # check for cluster when type = "dep"
  .user_input_check("cluster", type, cluster)

  # create fake unidimensional model syntax
  fac_name <- paste0(letters[sample(26, 5)], collapse = "")
  var_names <- colnames(sample_cov[[1]])
  while (fac_name %in% var_names) {
    fac_name <- paste0(letters[sample(26, 5)], collapse = "")
  }
  model <- paste0(
    fac_name, " =~ ", paste0(var_names, collapse = " + ")
  )

  # Run lavaan fit
  lav_fit <- lavaan::cfa(
    model,
    sample.cov = sample_cov, sample.nobs = sample_nobs, std.lv = TRUE,
    likelihood = "wishart", ceq.simple = TRUE,
    do.fit = FALSE
  )

  par_table <- lavaan::lavaanify(
    model,
    ceq.simple = TRUE, std.lv = TRUE
  )

  # Obtain data list for Stan
  data_list <- .create_data_list_meta(
    lavaan_object = lav_fit,
    method = "none",
    type = type,
    simple_struc = TRUE,
    priors = priors,
    cluster = cluster,
    correlation = TRUE,
    partab = par_table
  )
  var_names <- rownames(data_list$loading_pattern)
  data_list <- data_list[c(
    "Ng", "Np", "Ni", "r_obs_vec", "r_obs_vec_cov", "rm_i_l_par", "rm_i_s_par",
    "Nc", "C_ID", "type"
  )]

  message("User input fully processed :)\n Now to modeling.")

  pool_model <- instantiate::stan_package_model(
    name = "pool_cor", package = "bayesianmasem"
  )

  message("Fitting Stan model ...")

  pool_fit <- pool_model$sample(
    data = data_list,
    seed = seed,
    iter_warmup = warmup,
    iter_sampling = sampling,
    refresh = refresh,
    init = function() {
      list(
        ln_v_int_wi = array(data_list$rm_i_l_par, (data_list$type >= 2) * 1),
        ln_v_int_be = array(data_list$rm_i_l_par, (data_list$type == 3) * 1)
      )
    },
    adapt_delta = adapt_delta,
    max_treedepth = max_treedepth,
    chains = chains,
    parallel_chains = ncores,
    show_messages = show_messages
  )

  fit_indices <- posterior::summarise_draws(
    pool_fit$draws(c("ppp", "rmsea_mn", "rmsea_be", "rmsea_wi", "prop_be"))
  )
  fit_indices$variable <- c(
    "PPP", "Overall RMSEA", "Between RMSEA", "Within RMSEA", "% between"
  )
  r_mat_draws <- posterior::as_draws_matrix(pool_fit$draws("r_mat"))
  r_mat <- matrix(colMeans(r_mat_draws), nrow = data_list$Ni)
  r_mat_cov <- stats::cov(r_mat_draws[, which(lower.tri(r_mat))])
  colnames(r_mat) <- rownames(r_mat) <- var_names
  result <- list(
    fit_indices = fit_indices,
    r_mat = r_mat,
    r_mat_cov = r_mat_cov,
    data_list = data_list,
    stan_fit = pool_fit
  )

  return(result)
}
