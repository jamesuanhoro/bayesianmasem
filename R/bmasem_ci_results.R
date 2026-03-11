#' Conditional indepepence estimates for path analysis
#'
#' @description Obtain conditional independence estimates for path analysis.
#' Estimates are standardized.
#' @param object (bmasem) An object of \code{\link{bmasem-class}}
#' returned by \code{\link{bayesianmasem}}.
#' @param interval Confidence interval to select
#' @param summarize (LOGICAL)
#' If TRUE (default): Return posterior summary as data.frame;
#' If FALSE: Return posterior draws as data.frame.
#' @returns Table of parameters related to conditional independence.
#' @examples
#' \dontrun{
#' model_syntax <- paste(
#'   "BI ~ ATT + SN + PBC",
#'   "BEH ~ PBC + BI",
#'   sep = "\n"
#' )
#' fit <- bmasem(
#'   model_syntax,
#'   sample_cov = cooke16$data, sample_nobs = cooke16$n
#' )
#' ci_table <- bmasem_ci_results(fit)
#' }
#' @export
bmasem_ci_results <- function(object, interval = .9, summarize = TRUE) {
  stopifnot(inherits(object, "bmasem"))

  if (object@data_list$method >= 90) {
    stop(paste0(
      "There are no residuals to plot when ",
      "method == \"none\", \"WB\", \"WB-cond\", \"WW\"."
    ))
  }

  if (object@data_list$pa_indicator != 1) {
    stop(paste0(
      "This method only works for objects produced by ",
      "minorbpa()."
    ))
  }

  if (sum(object@data_list$cond_ind_mat) == 0) {
    message(
      "All possible associations are modelled."
    )
    return(NULL)
  }

  ind_names <- rownames(object@data_list$loading_pattern)
  ci_details <- object@data_list$cond_ind_details
  x_var <- ind_names[ci_details[, 1]]
  y_var <- ind_names[ci_details[, 2]]
  z_var <- vector("character", nrow(ci_details))
  for (i in seq_len(nrow(ci_details))) {
    z_var[i] <- paste0(
      ind_names[as.logical(ci_details[i, 3:ncol(ci_details)])],
      collapse = ", "
    )
  }
  var_s <- paste0(x_var, " _||_ ", y_var, " | ", z_var)

  resid_idxs <- paste0("Resid[", ci_details[, 1], ",", ci_details[, 2], "]")

  if (isFALSE(summarize)) {
    result <- posterior::as_draws_df(object@stan_fit$draws(resid_idxs))
    result <- result[
      ,
      c(resid_idxs, ".chain", ".iteration", ".draw")
    ]
    colnames(result)[seq_along(var_s)] <- var_s
  } else {
    resid_idxs_i <- as.integer(factor(resid_idxs, unique(resid_idxs)))
    post_sum <- .bmasem_post_sum(
      object@stan_fit, resid_idxs, interval = interval
    )
    result <- cbind(relation = var_s, post_sum[resid_idxs_i, -1])
    rownames(result) <- NULL
  }

  return(result)
}
