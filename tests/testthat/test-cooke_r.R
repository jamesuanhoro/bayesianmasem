test_that("Corr: Dependent samples works for PA on Cooke16", {
  method <- sample(.method_hash(), 1)
  model_syntax <- paste(
    "BI ~ ATT + SN + PBC",
    "BEH ~ PBC + BI",
    sep = "\n"
  )
  cluster_names <- gsub(
    ",", "", lapply(strsplit(names(cooke16$data), " "), "[[", 1)
  )
  subset <- 1:5
  cluster_ids <- as.integer(as.factor(cluster_names[subset]))
  expect_error(fit <- bmasem(
    model_syntax,
    sample_cov = cooke16$data[subset], sample_nobs = cooke16$n[subset],
    cluster = cluster_ids,
    warmup = test_warm, sampling = test_samp, chains = test_chns,
    method = method, type = "dep",
    refresh = 0, show_messages = FALSE
  ), NA)
  expect_true(all(slotNames(fit) %in% c(
    "major_parameters", "minor_factor_matrix", "data_list",
    "priors", "stan_fit", "version"
  )))
  bmasem_test_hist(fit)
  bmasem_test_trace(fit)
  bmasem_test_pa_ci(fit, TRUE)
  bmasem_test_pa_ci(fit, FALSE)
  expect_error(
    print_out <- capture_output(pp_summary(fit), width = 300),
    NA
  )
  bmasem_test_pp_pa_shared(print_out, method)
  expect_error(
    print_out <- capture_output(pp_summary(fit, simple = FALSE), width = 300),
    NA
  )
  bmasem_test_pp_pa_shared(print_out, method, FALSE)
})
