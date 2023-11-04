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

test_warm <- 300
test_samp <- 300
test_chns <- 1

method <- sample(.method_hash(), 1)
model_syntax <- paste0(
  "distress =~ ", paste0("x", 1:14, collapse = " + "), "\n",
  "anxiety =~ ", paste0("x", seq(1, 14, 2), collapse = " + "), "\n",
  "depression =~ ", paste0("x", seq(2, 14, 2), collapse = " + ")
)
cluster_names <- gsub(
  ",", "", lapply(strsplit(names(Norton13$data), " "), "[[", 1)
)
subset <- c(1, 19:21)
cluster_ids <- as.integer(as.factor(cluster_names[subset]))
fit <- bmasem(
  model_syntax,
  sample_cov = Norton13$data[subset], sample_nobs = Norton13$n[subset],
  cluster = cluster_ids,
  simple_struc = TRUE, orthogonal = TRUE,
  warmup = test_warm, sampling = test_samp, chains = test_chns,
  method = method, type = "dep",
  refresh = 0, show_messages = TRUE
)
