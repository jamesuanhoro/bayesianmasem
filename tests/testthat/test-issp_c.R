test_that("Cov: Random method works for meta-CFA on issp89", {
  print(type <- sample(c("fe", "re"), 1))
  method <- sample(.method_hash(), 1)
  model_syntax <- "# latent variable definitions
    F1 =~ JP1 + JP2 + JP3
    F2 =~ JN1 + JN2 + JN3 + JN4
    F3 =~ TD1 + TD2
    JN4 ~~ JP2"

  expect_error(fit <- bmasem(
    model_syntax,
    sample_cov = issp89$data[1:3], sample_nobs = issp89$n[1:3],
    orthogonal = sample(c(TRUE, FALSE), 1),
    simple_struc = sample(c(TRUE, FALSE), 1),
    warmup = test_warm, sampling = test_samp, chains = test_chns,
    method = method, type = type,
    refresh = 0, show_messages = FALSE, correlation = FALSE
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
  bmasem_test_pp_shared(print_out, method, error = TRUE, corr = FALSE)
  expect_error(
    print_out <- capture_output(pp_summary(fit, simple = FALSE), width = 300),
    NA
  )
  bmasem_test_pp_shared(print_out, method, FALSE, TRUE, corr = FALSE)
})
