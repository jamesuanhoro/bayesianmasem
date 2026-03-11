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

#' Another dataset from the metaSEM package.
#'
#' The data set includes 131 correlation matrices on the
#' mediating role of leader–member exchange in transformational
#' leadership
#' \insertCite{boer_revisiting_2016}{bayesianmasem}.
#' Data and documentation copied over from the metaSEM package
#' \insertCite{Cheung-metaSEM}{bayesianmasem}.
#'
#' @format ## `boer16`
#' A list of data with the following structure:
#' \describe{
#'   \item{data}{
#'     A list of correlation matrices. The variables are LMX
#'     (leader-member exchange), TFL (transformational leadership),
#'     JS (job satisfaction), OC (organizational commitment), and
#'     LE (leader effectiveness).
#'   }
#'   \item{n}{A vector of sample sizes}
#'   \item{RelLMX}{The reliability of LMX}
#'   \item{RelTFL}{The reliability of TFL}
#' }
#' @examples
#' \dontrun{
#' model_syntax <- paste(
#'   "LMX ~ TFL",
#'   "JS + OC + LE ~ TFL + LMX",
#'   sep = "\n"
#' )
#' bmasem(model_syntax, sample_cov = boer16$data, sample_nobs = boer16$n)
#' }
#' @references \insertAllCited{}
#' @source <https://cran.r-project.org/web/packages/metaSEM/index.html>
"boer16"

#' Another dataset from the metaSEM package.
#'
#' The data set includes 34 correlation matrices on the
#' how well the theory of planned behaviour predicts alcohol consumption
#' \insertCite{cooke_how_2016}{bayesianmasem}.
#' We adjusted the third correlation matrix to the nearest
#' valid correlation matrix.
#' Data and documentation copied over from the metaSEM package
#' \insertCite{Cheung-metaSEM}{bayesianmasem}.
#'
#' @format ## `cooke16`
#' A list of data with the following structure:
#' \describe{
#'   \item{data}{
#'     A list of correlation matrices.
#'     The variables are SN (subjective norm), ATT (attitude),
#'     PBC (perceived behavior control), BI (behavioral intention),
#'     and BEH (behavior).
#'   }
#'   \item{n}{A vector of sample sizes}
#'   \item{MeanAge}{
#'     Mean age of the participants except for Ajzen and Sheikh (2013),
#'     which is the median age, and Glassman, et al. (2010a) to
#'     Glassman, et al. (2010d), which are based on the range of 18 to 24.
#'   }
#'   \item{Female}{Percentage of female participants.}
#' }
#' @examples
#' \dontrun{
#' model_syntax <- paste(
#'   "BI ~ ATT + SN + PBC",
#'   "BEH ~ PBC + BI",
#'   sep = "\n"
#' )
#' bmasem(model_syntax, sample_cov = cooke16$data, sample_nobs = cooke16$n)
#' }
#' @references \insertAllCited{}
#' @source <https://cran.r-project.org/web/packages/metaSEM/index.html>
"cooke16"
