test_that("Corr: Dependent samples works for meta-CFA on Norton13", {
  # expect_error(fit <- bmasem(
  #   model_syntax,
  #   sample_cov = Norton13$data[subset], sample_nobs = Norton13$n[subset],
  #   cluster = cluster_ids,
  #   simple_struc = TRUE, orthogonal = TRUE,
  #   warmup = test_warm, sampling = test_samp, chains = test_chns,
  #   method = method, type = "dep",
  #   refresh = 0, show_messages = TRUE
  # ), NA)
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
