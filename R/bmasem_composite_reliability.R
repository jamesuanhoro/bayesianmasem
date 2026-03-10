#' Obtain composite reliability from a fitted model
#'
#' @description Get composite reliability
#' @param object (bmasem) An object of \code{\link{bmasem-class}}
#' returned by \code{\link{bmasem}}.
#' @param interval (real in (0, 1)) Credible interval to return.
#' @param return_draws (LOGICAL) If `TRUE`, returns the realised distribution
#' of reliability.
#' If `FALSE` (default), return the distribution summary.
#' @returns
#' - If \code{return_draws = FALSE},
#'   returns the summary of posterior distribution of composite reliability
#'   for each factor.
#' - If \code{return_draws = TRUE},
#'   returns draws from the realised distribution of composite reliability
#'   for each factor.
#' @details This metric is the ratio of common variance (CV) to
#' total variance (TV) for each factor.
#'
#' - **Common Variance**: For factor \eqn{f}, the common variance is:
#'   \deqn{
#'   \mathrm{CV}_f = (\sum_{i \in I_f} \lambda_{if})^2 \phi_{ff},
#'   }
#'   where \eqn{I_f} indexes the indicators of factor \eqn{f}.
#' - **Total Variance**: For factor \eqn{f}, the total variance is:
#'   \deqn{
#'   \mathrm{TV}_f = \sum_{i \in I_f} \sum_{j \in I_f} \omega_{ij},
#'   }
#'   where \eqn{\Omega = = (\omega_{ij})} is the model-implied covariance
#'   matrix of the indicators.
#' @examples
#' \dontrun{
#' model_syntax <- paste0(
#'   "distress =~ ", paste0("x", 1:14, collapse = " + "), "\n",
#'   "anxiety =~ ", paste0("x", seq(1, 14, 2), collapse = " + "), "\n",
#'   "depression =~ ", paste0("x", seq(2, 14, 2), collapse = " + ")
#' )
#' fit <- bmasem(
#'   model_syntax,
#'   sample_cov = Norton13$data, sample_nobs = Norton13$n, orthogonal = TRUE
#' )
#' bmasem_composite_reliability(fit)
#' }
#' @export
bmasem_composite_reliability <- function(
    object, interval = .9, return_draws = FALSE) {
  stopifnot(inherits(object, "bmasem"))

  if (interval <= 0 || interval >= 1 || !is.numeric(interval)) {
    stop("Interval must be a number between 0 and 1")
  }

  data_list <- object@data_list
  if (data_list$pa_indicator == 1) {
    stop("Model must be a CFA not a path analysis.")
  }

  stan_draws <- posterior::as_draws(object@stan_fit)

  index_mat <- 1 * (
    data_list$loading_pattern >= ifelse(data_list$complex_struc == 1, -999, 1) |
      data_list$loading_fixed != -999
  )
  n_factors <- ncol(data_list$loading_pattern)
  n_items <- nrow(data_list$loading_pattern)
  load_draws <- posterior::subset_draws(stan_draws, variable = "Load_mat")
  phi_draws <- posterior::subset_draws(stan_draws, variable = "phi_mat")
  omega_draws <- posterior::subset_draws(stan_draws, variable = "Omega")

  draw_dims <- dim(load_draws)

  reliability_array <- array(
    dim = c(draw_dims[1:2], n_factors),
    dimnames = list(NULL, NULL, colnames(data_list$loading_pattern))
  )
  for (i in seq_len(draw_dims[1])) {
    for (j in seq_len(draw_dims[2])) {
      load_mat <- matrix(load_draws[i, j, ], n_items)
      phi_mat <- matrix(phi_draws[i, j, ], n_factors)
      omega_mat <- matrix(omega_draws[i, j, ], n_items)
      for (k in seq_len(n_factors)) {
        which_items <- which(index_mat[, k] == 1)
        load_mat_k <- load_mat[which_items, k, drop = FALSE]
        reliability_array[i, j, k] <- sum(
          tcrossprod(
            load_mat_k %*% phi_mat[k, k, drop = FALSE], load_mat_k
          )
        ) / sum(omega_mat[which_items, which_items, drop = FALSE])
      }
    }
  }

  reliability_array <- posterior::as_draws_array(reliability_array)
  if (isTRUE(return_draws)) {
    return(reliability_array)
  }

  result <- .bmasem_post_sum_draws(reliability_array, interval)
  return(result)
}
