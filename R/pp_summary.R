#' Pretty print model results
#'
#' @description Nice printing of model results,
#' optionally produces HTML document
#' @param object (bmasem) An object of \code{\link{bmasem-class}}
#' returned by \code{\link{bmasem}}.
#' @param interval (real in (0, 1)) Credible interval to return.
#' @param digits (positive integer) Number of decimal places to print in table
#' @param simple (Logical) TRUE to produce table with less information
#' about parameters;
#' FALSE: produces table with more information
#' @param save_html (string) Optional file name to save table as HTML
#' @returns NULL
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
#' pp_summary(fit)
#' pp_summary(fit, simple = FALSE)
#' }
#' @export
pp_summary <- function(
    object, interval = .9, digits = 3, simple = TRUE,
    save_html = NULL) {
  stopifnot(inherits(object, "bmasem"))

  if (interval <= 0 || interval >= 1 || !is.numeric(interval)) {
    stop("Interval must be a number between 0 and 1")
  }

  table_to_print <- .create_major_params(
    stan_fit = object@stan_fit,
    data_list = object@data_list,
    interval = interval
  )
  method_str <- .method_hash(object@data_list$method)
  n_obs <- paste0(object@data_list$Np, collapse = ", ")

  if (isTRUE(simple)) {
    if (object@data_list$method == 4) {
      table_to_print[2, "mean"] <- table_to_print[2, "median"]
      table_to_print[2, "sd"] <- table_to_print[2, "mad"]
    }
    # Use index because of quantile name changes
    table_to_print <- table_to_print[, c(1:5, 7, 9:10, 11:12)]
  }

  if (object@data_list$method == 100) {
    # Remove RMSE index
    table_to_print[1, 5:ncol(table_to_print)] <- NA_real_
  }

  if (object@data_list$type == 1) {
    # Remove all RMSEAs
    table_to_print[2:5, 5:ncol(table_to_print)] <- NA_real_
  } else if (object@data_list$type == 2) {
    table_to_print[3:5, 5:ncol(table_to_print)] <- NA_real_
  }

  type_str <- paste0(
    .type_hash(object@data_list$type, elaborate = TRUE),
    ", "
  )

  caption_str <- paste0(
    "Parameter estimates (",
    type_str,
    "method = ", method_str,
    ", sample size(s) = ", n_obs, ")"
  )

  result <- huxtable::huxtable(table_to_print) |>
    huxtable::column_to_header("group") |>
    huxtable::set_caption(caption_str)

  gof_row <- which(result$from == "Goodness of fit")
  disp_row <- which(result$from == "Dispersion between and within clusters")
  gof_disp_rows <- c(gof_row + 1, disp_row + 1:4)

  for (row_id in gof_disp_rows) {
    result <- huxtable::merge_cells(result, row_id, 1:3)
  }

  footnote_str <- ""

  if (isTRUE(simple) && object@data_list$method == 4) {
    footnote_str <- paste0(
      footnote_str, "\n",
      "Mean and SD RMSE are median and mad respectively ",
      "because RMSE is usually extremely right-skewed."
    )
  }

  result <- huxtable::add_footnote(result, text = footnote_str)

  if (isTRUE(simple)) {
    huxtable::number_format(result)[, 4:8] <- list(
      function(x) trimws(format(round(x, digits), nsmall = digits))
    )
    result <- huxtable::set_number_format(result, col = 9, value = 0)
  } else if (isFALSE(simple)) {
    huxtable::number_format(result)[, 4:10] <- list(
      function(x) trimws(format(round(x, digits), nsmall = digits))
    )
    result <- huxtable::set_number_format(result, col = 11:12, value = 0) |>
      huxtable::insert_row(
        "Relation", rep("", 2),
        "Location", rep("", 1),
        "Dispersion", rep("", 3),
        "Parameter convergence", rep("", 2),
        after = 0
      ) |>
      huxtable::merge_cells(1, 1:3) |>
      huxtable::merge_cells(1, 4:5) |>
      huxtable::merge_cells(1, 6:9) |>
      huxtable::merge_cells(1, 10:12) |>
      huxtable::set_align(1, huxtable::everywhere, "center") |>
      huxtable::set_valign(1, huxtable::everywhere, "bottom") |>
      huxtable::set_tb_padding(1, huxtable::everywhere, 10) |>
      huxtable::set_header_rows(1, TRUE) |>
      huxtable::set_bottom_border(
        1,
        col = huxtable::everywhere
      ) |>
      huxtable::set_right_border(
        row = huxtable::everywhere, col = c(3, 5, 9)
      )
  }

  header_rows <- rev(which(
    rownames(result) == "" | substr(rownames(result), 1, 1) == "."
  ))[-1]
  header_rows <- c(header_rows, header_rows - 1)
  result <- huxtable::set_bottom_border(
    result, header_rows, huxtable::everywhere
  ) |>
    huxtable::style_headers(bold = TRUE)

  if (!is.null(save_html) && is.character(save_html)) {
    huxtable::quick_html(result, file = save_html)
  }

  print(result)

  return()
}
