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

bmasem_test_comp_rel <- function(fit, ret_draws = FALSE) {
  testthat::expect_error(
    cr <- bmasem_composite_reliability(fit, return_draws = ret_draws),
    NA
  )
  if (isTRUE(ret_draws)) {
    testthat::expect_true(inherits(cr, "draws_array"))
  } else {
    testthat::expect_true(inherits(cr, "data.frame"))
  }
}

bmasem_test_res_net <- function(fit, method) {
  if (.method_hash(method) >= 90) {
    testthat::expect_error(
      bmasem_residual_network(fit),
      "There are no residuals to plot when method == \"none\"."
    )
  } else {
    testthat::expect_error(
      b_res_net <- bmasem_residual_network(fit),
      NA
    )
    testthat::expect_true(inherits(b_res_net, "draws_df"))
  }
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

bmasem_test_pp_pa_shared <- function(
    print_out, method, simple = TRUE) {
  testthat::expect_true(grepl(method, print_out, ignore.case = TRUE))
  testthat::expect_true(grepl("Goodness of fit", print_out, ignore.case = TRUE))
  testthat::expect_true(grepl("PPP", print_out, ignore.case = TRUE))
  testthat::expect_true(grepl("RMSE", print_out, ignore.case = TRUE))
  testthat::expect_true(grepl("RMSEA", print_out, ignore.case = TRUE))
  testthat::expect_true(
    grepl("R square", print_out, ignore.case = TRUE)
  )
  testthat::expect_true(
    grepl("Regression coefficients", print_out, ignore.case = TRUE)
  )
  if (isFALSE(simple)) {
    testthat::expect_true(grepl("Location", print_out, ignore.case = TRUE))
    testthat::expect_true(grepl("Dispersion", print_out, ignore.case = TRUE))
    testthat::expect_true(grepl("convergence", print_out, ignore.case = TRUE))
    testthat::expect_true(grepl("median", print_out, ignore.case = TRUE))
    testthat::expect_true(grepl("mad", print_out, ignore.case = TRUE))
  }
  testthat::expect_true(
    grepl("Error correlations", print_out, ignore.case = TRUE)
  )
}

bmasem_test_pa_ci <- function(fit, summarize) {
  if (fit@data_list$method < 90 && sum(fit@data_list$cond_ind_mat) > 0) {
    if (isFALSE(summarize)) {
      testthat::expect_true(
        "draws_df" %in% class(bmasem_ci_results(fit, summarize = FALSE))
      )
    } else {
      testthat::expect_true(
        class(bmasem_ci_results(fit, summarize = TRUE)) == "data.frame"
      )
    }
  } else {
    if (fit@data_list$method >= 90) {
      err_msg <- paste0(
        "There are no residuals to plot when ",
        "method == \"none\", \"WB\", \"WB-cond\", \"WW\"."
      )
      testthat::expect_error(bmasem_ci_results(fit), err_msg)
    } else {
      if (sum(fit@data_list$cond_ind_mat) == 0) {
        msg <- paste0(
          "All possible associations are modelled."
        )
        testthat::expect_message(bmasem_ci_results(fit), msg)
      }
    }
  }
}

test_warm <- 300
test_samp <- 300
test_chns <- 1
