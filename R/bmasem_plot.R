#' Pretty print model results
#'
#' @description Plot model parameters
#' @param object (bmasem) An object of \code{\link{bmasem-class}}
#' returned by \code{\link{bmasem}}.
#' @param type (character) Type of plot:
#' "trace" for traceplots and "hist" for histograms.
#' @param subset (character) Subset of parameters:
#' NULL (Default) showing all estimated parameters;
#' Any other response will be used as regular expressions to subset
#' the parameters. It can be loading names or types of parameters.
#' @param ... additional arguments to relevant bayesplot function
#' @returns bayesplot object
#' @examples
#' \dontrun{
#' model_syntax <- paste0(
#'   "distress =~ ", paste0("x", 1:14, collapse = " + "), "\n",
#'   "anxiety =~ ", paste0("x", seq(1, 14, 2), collapse = " + "), "\n",
#'   "depression =~ ", paste0("x", seq(2, 14, 2), collapse = " + ")
#' )
#' bmasem(
#'   model_syntax,
#'   sample_cov = Norton13$data, sample_nobs = Norton13$n, orthogonal = TRUE
#' )
#' bmasem_plot(fit)
#' }
#' @export
bmasem_plot <- function(object, type = "trace", subset = NULL, ...) {
  stopifnot(inherits(object, "bmasem"))
  stopifnot(type %in% c("hist", "trace"))

  param_list <- .get_param_list(object@data_list)
  if (!is.null(subset)) {
    param_names <- names(param_list)
    param_list <- param_list[grep(subset, param_names, ignore.case = TRUE)]
  }
  pd_mat <- posterior::as_draws_matrix(object@stan_fit$draws(param_list))
  colnames(pd_mat)[seq_along(param_list)] <- names(param_list)

  if (type == "trace") {
    result <- bayesplot::mcmc_trace(pd_mat, ...)
  } else if (type == "hist") {
    result <- bayesplot::mcmc_hist(pd_mat, ...)
  }

  return(result)
}
