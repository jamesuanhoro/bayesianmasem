test_that("Dependent samples works for meta-CFA on issp89", {
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
  expect_error(fit <- bmasem(
    model_syntax,
    sample_cov = Norton13$data[subset], sample_nobs = Norton13$n[subset],
    cluster = cluster_ids,
    simple_struc = TRUE, orthogonal = TRUE,
    warmup = 500, sampling = 500, chains = 1,
    method = method, type = "dep",
    refresh = 0, show_messages = FALSE
  ), NA)
  expect_true(all(slotNames(fit) %in% c(
    "major_parameters", "minor_factor_matrix", "data_list",
    "priors", "stan_fit", "version"
  )))
  bmasem_test_hist(fit)
  bmasem_test_trace(fit)
  expect_error(
    print_out <- capture_output(pp_summary(fit), width = 300),
    NA
  )
  bmasem_test_pp_shared(print_out, method)
  expect_error(
    print_out <- capture_output(pp_summary(fit, simple = FALSE), width = 300),
    NA
  )
  bmasem_test_pp_shared(print_out, method, FALSE)
})
