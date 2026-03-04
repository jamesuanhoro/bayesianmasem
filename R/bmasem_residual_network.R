#' Compute residual network
#'
#' @description Interpret the error correlations as a residual network model
#' \insertCite{epskamp_generalized_2017}{bayesianmasem}.
#' @param object (bmasem) An object of \code{\link{bmasem-class}}
#' returned by \code{\link{bmasem}}.
#' @returns A data.frame containing posterior samples of the partial
#' correlation matrix.
#' @examples
#' \dontrun{
#' model_syntax <- paste("distress =~", paste0("x", 1:14, collapse = " + "))
#' fit <- bmasem(
#'   model_syntax,
#'   sample_cov = Norton13$data, sample_nobs = Norton13$n
#' )
#' res_net <- bmasem_residual_network(fit)
#' p_corr_df <- posterior::summarise_draws(res_net)
#' n_items <- sqrt(nrow(p_corr_df))
#' p_corr_mat <- matrix(p_corr_df$mean, n_items)
#' p_corr_mat
#' qgraph::qgraph(p_corr_mat)
#' }
#' @references \insertAllCited{}
#' @export
bmasem_residual_network <- function(object) {
  stopifnot(inherits(object, "bmasem"))

  if (object@data_list$method >= 90) {
    stop(paste0(
      "There are no residuals to plot when method == \"none\"."
    ))
  }

  resids <- posterior::as_draws_df(object@stan_fit$draws("Resid"))
  len <- ncol(resids) - 3 # exclude special variables for draws
  p <- sqrt(len)
  result <- t(apply(resids, 1, function(x) {
    mat <- matrix(x[seq_len(len)], p, p)
    diag(mat) <- 1
    inv_mat <- -1 * stats::cov2cor(pracma::pinv(mat))
    diag(inv_mat) <- 1
    inv_vec <- c(as.numeric(inv_mat), x[(len + 1):(len + 3)])
    inv_vec
  }))

  colnames(result) <- gsub("Resid", "p_corr", colnames(resids))
  result <- posterior::as_draws_df(as.data.frame(result))

  result
}
