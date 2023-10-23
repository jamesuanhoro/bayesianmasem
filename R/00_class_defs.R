#' A class for setting up priors.
#'
#' @slot lkj_shape (positive real) The shape parameter of the LKJ-prior on the
#' interfactor correlation matrix in confirmatory factor models.
#' @slot ml_par (real) The location parameter of the normal prior on loadings.
#' @slot sl_par (positive real) The scale parameter
#' of the normal prior on loadings.
#' @slot rc_par (positive real) The shape parameter of the Beta(rc_par, rc_par)
#' prior on the residual error correlations.
#' @slot mr_par (real) The location parameter of the normal prior
#' on the log-RMSEA.
#' @slot sr_par (positive real) The scale parameter of the normal prior
#' on the log-RMSEA.
#' @slot br_par (positive real) The scale parameter of the normal prior
#' on the regression coefficients for the log-RMSEA.
#' @slot rm_par (positive real) The scale parameter of the normal prior
#' on the tau / CRMR parameter.
#'
#' @name bmasempriors-class
#' @rdname bmasempriors-class
#' @export
methods::setClass(
  "bmasempriors",
  methods::representation(
    lkj_shape = "numeric",
    ml_par = "numeric",
    sl_par = "numeric",
    rc_par = "numeric",
    mr_par = "numeric",
    sr_par = "numeric",
    br_par = "numeric",
    rm_par = "numeric"
  ),
  prototype = list(
    lkj_shape = 2.0,
    ml_par = 0.0,
    sl_par = 0.5,
    rc_par = 2.0,
    mr_par = log(0.8),
    sr_par = 0.7,
    br_par = 0.5,
    rm_par = 0.15
  )
)

#' Set priors in package
#'
#' @description Modify default priors in package.
#' @param lkj_shape (positive real) The shape parameter of the LKJ-prior on the
#' interfactor correlation matrix in confirmatory factor models.
#' @param ml_par (real) The location parameter of the normal prior on loadings.
#' @param sl_par (positive real) The scale parameter
#' of the normal prior on loadings.
#' @param rc_par (positive real) The shape parameter of the Beta(rc_par, rc_par)
#' prior on the residual error correlations.
#' @param mr_par (real) The location parameter of the normal prior
#' on the log-RMSEA.
#' @param sr_par (positive real) The scale parameter of the normal prior
#' on the log-RMSEA.
#' @param br_par (positive real) The scale parameter of the normal prior
#' on the regression coefficients for the log-RMSEA.
#' @param rm_par (positive real) The scale parameter of the normal prior
#' on the tau / CRMR parameter.
#' @returns An object of \code{\link{bmasempriors-class}}
#' @examples
#' \dontrun{
#' # Change LKJ shape parameter only
#' custom_priors <- new_bmasempriors(lkj_shape = 1.0)
#' model_syntax <- paste0(
#'   "distress =~ ", paste0("x", 1:14, collapse = " + "), "\n",
#'   "anxiety =~ ", paste0("x", seq(1, 14, 2), collapse = " + "), "\n",
#'   "depression =~ ", paste0("x", seq(2, 14, 2), collapse = " + ")
#' )
#' bmasem(
#'   model_syntax,
#'   sample_cov = Norton13$data, sample_nobs = Norton13$n, orthogonal = TRUE,
#'   priors = custom_priors
#' )
#' }
#' @export
new_bmasempriors <- function(
    lkj_shape = 2.0,
    ml_par = 0.0,
    sl_par = 0.5,
    rc_par = 2.0,
    mr_par = log(0.8),
    sr_par = 0.7,
    br_par = 0.5,
    rm_par = 0.15) {
  bm_priors_object <- methods::new("bmasempriors")
  bm_priors_object <- methods::initialize(
    bm_priors_object,
    lkj_shape = lkj_shape,
    ml_par = ml_par,
    sl_par = sl_par,
    rc_par = rc_par,
    mr_par = mr_par,
    sr_par = sr_par,
    br_par = br_par,
    rm_par = rm_par
  )
  return(bm_priors_object)
}
