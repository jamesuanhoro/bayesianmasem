bmasem_test_hist <- function(fit) {
  testthat::expect_error(
    gg <- bmasem_plot(fit, type = "hist"),
    NA
  )
  testthat::expect_true(inherits(gg, "ggplot"))
}

bmasem_test_trace <- function(fit) {
  testthat::expect_error(
    gg <- bmasem_plot(fit, type = "trace"),
    NA
  )
  testthat::expect_true(inherits(gg, "ggplot"))
}

bmasem_test_pp_shared <- function(
    print_out, method, simple = TRUE, error = FALSE, corr = TRUE) {
  testthat::expect_true(grepl(method, print_out, ignore.case = TRUE))
  testthat::expect_true(grepl("Goodness of fit", print_out, ignore.case = TRUE))
  testthat::expect_true(grepl("PPP", print_out, ignore.case = TRUE))
  testthat::expect_true(grepl("RMSE", print_out, ignore.case = TRUE))
  testthat::expect_true(grepl("RMSEA", print_out, ignore.case = TRUE))
  testthat::expect_true(
    grepl("Inter-factor correlations", print_out, ignore.case = TRUE)
  )
  testthat::expect_true(grepl("Factor loadings", print_out, ignore.case = TRUE))
  if (isFALSE(simple)) {
    testthat::expect_true(grepl("Location", print_out, ignore.case = TRUE))
    testthat::expect_true(grepl("Dispersion", print_out, ignore.case = TRUE))
    testthat::expect_true(grepl("convergence", print_out, ignore.case = TRUE))
    testthat::expect_true(grepl("median", print_out, ignore.case = TRUE))
    testthat::expect_true(grepl("mad", print_out, ignore.case = TRUE))
  }
  if (isTRUE(error)) {
    testthat::expect_true(
      grepl("Error correlations", print_out, ignore.case = TRUE)
    )
  }
  if (isFALSE(corr)) {
    testthat::expect_true(
      grepl("Residual variances", print_out, ignore.case = TRUE)
    )
  }
}
