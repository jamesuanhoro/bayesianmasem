test_that("Corr: Random method works for PA on boer16", {
  print(type <- sample(c("fe", "re"), 1))
  method <- sample(.method_hash(), 1)
  model_syntax <- paste(
    "LMX ~ TFL",
    "JS + OC + LE ~ TFL + LMX",
    sep = "\n"
  )
  subset <- c(23, 48, 66, 77)
  expect_error(fit <- bmasem(
    model_syntax,
    sample_cov = boer16$data[subset], sample_nobs = boer16$n[subset],
    warmup = test_warm, sampling = test_samp, chains = test_chns,
    method = method, type = type,
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
  bmasem_test_pp_pa_shared(print_out, method)
  expect_error(
    print_out <- capture_output(pp_summary(fit, simple = FALSE), width = 300),
    NA
  )
  bmasem_test_pp_pa_shared(print_out, method, FALSE)
})
