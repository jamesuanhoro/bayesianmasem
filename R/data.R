#' Studies on the Hospital Anxiety and Depression Scale.
#'
#' The data set includes 28 studies on 14 items measuring the Hospital Anxiety
#' and Depression Scale (HADS) Reported by
#' \insertCite{norton_hospital_2013;textual}{bayesianmasem}, with data and
#' documentation copied over from the metaSEM package
#' \insertCite{Cheung-metaSEM}{bayesianmasem}.
#'
#' @format ## `Norton13`
#' A list of data with the following structure:
#' \describe{
#'   \item{data}{
#'      A list of 28 studies of correlation matrices.
#'      The variables are 14 items (x1 to x14) measuring HADS}
#'   \item{n}{A vector of sample sizes}
#'   \item{population}{A vector of the population of the data}
#'   \item{group}{
#'      A vector of classification into patients vs. non-patients
#'      based on population}
#' }
#' @references \insertAllCited{}
#' @source <https://cran.r-project.org/web/packages/metaSEM/index.html>
"Norton13"

#' Another dataset from the metaSEM package.
#'
#' The data set includes 11 studies on 9 items measuring work-related attitudes
#' \insertCite{noauthor_international_1992}{bayesianmasem}.
#' Data and documentation copied over from the metaSEM package
#' \insertCite{Cheung-metaSEM}{bayesianmasem}.
#'
#' @format ## `issp89`
#' A list of data with the following structure:
#' \describe{
#'   \item{data}{A list of 11 studies of covariance matrices}
#'   \item{n}{A vector of sample sizes}
#' }
#' @references \insertAllCited{}
#' @source <https://cran.r-project.org/web/packages/metaSEM/index.html>
"issp89"
