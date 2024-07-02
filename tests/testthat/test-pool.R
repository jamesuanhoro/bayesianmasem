test_that("Random type works for pooling on issp89", {
  print(type <- sample(c("fe", "re"), 1))

  expect_error(stage_1 <- bmasem_stage_1(
    sample_cov = issp89$data[1:3], sample_nobs = issp89$n[1:3],
    warmup = test_warm, sampling = test_samp, chains = test_chns,
    type = type, refresh = 0, show_messages = FALSE
  ), NA)
  expect_true(all(names(stage_1) %in% c(
    "fit_indices", "r_mat", "r_mat_cov", "p_mat", "p_mat_cov",
    "data_list", "stan_fit"
  )))

  method <- sample(.method_hash(), 1)
  model_syntax <- "# latent variable definitions
    F1 =~ JP1 + JP2 + JP3
    F2 =~ JN1 + JN2 + JN3 + JN4
    F3 =~ TD1 + TD2
    JN4 ~~ JP2"
  expect_error(stage_2 <- bmasem_stage_2(
    model_syntax, stage_1,
    method = method,
    orthogonal = sample(c(TRUE, FALSE), 1),
    simple_struc = sample(c(TRUE, FALSE), 1),
    warmup = test_warm, sampling = test_samp, chains = test_chns,
    refresh = 0, show_messages = FALSE
  ), NA)
  expect_true(all(slotNames(stage_2) %in% c(
    "major_parameters", "minor_factor_matrix", "data_list",
    "priors", "stan_fit", "version"
  )))
  bmasem_test_hist(stage_2)
  bmasem_test_trace(stage_2)
  expect_error(
    print_out <- capture_output(pp_summary(stage_2), width = 300),
    NA
  )
  bmasem_test_pp_shared(print_out, method, error = TRUE)
  expect_error(
    print_out <- capture_output(
      pp_summary(stage_2, simple = FALSE), width = 300
    ),
    NA
  )
  bmasem_test_pp_shared(print_out, method, FALSE, TRUE)
})

test_that("Dependent-samples works for pooling on Norton13", {
  cluster_names <- gsub(
    ",", "", lapply(strsplit(names(Norton13$data), " "), "[[", 1)
  )
  subset <- c(1, 19:21)
  cluster_ids <- as.integer(as.factor(cluster_names[subset]))

  expect_error(stage_1 <- bmasem_stage_1(
    sample_cov = Norton13$data[subset], sample_nobs = Norton13$n[subset],
    cluster = cluster_ids,
    warmup = test_warm, sampling = test_samp, chains = test_chns,
    type = "dep", refresh = 0, show_messages = FALSE
  ), NA)
  expect_true(all(names(stage_1) %in% c(
    "fit_indices", "r_mat", "r_mat_cov", "p_mat", "p_mat_cov",
    "data_list", "stan_fit"
  )))

  method <- sample(.method_hash(), 1)
  model_syntax <- paste0(
    "distress =~ ", paste0("x", 1:14, collapse = " + "), "\n",
    "anxiety =~ ", paste0("x", seq(1, 14, 2), collapse = " + "), "\n",
    "depression =~ ", paste0("x", seq(2, 14, 2), collapse = " + ")
  )
  expect_error(stage_2 <- bmasem_stage_2(
    model_syntax, stage_1,
    method = method,
    orthogonal = TRUE,
    simple_struc = TRUE,
    warmup = test_warm, sampling = test_samp, chains = test_chns,
    refresh = 0, show_messages = FALSE
  ), NA)
  expect_true(all(slotNames(stage_2) %in% c(
    "major_parameters", "minor_factor_matrix", "data_list",
    "priors", "stan_fit", "version"
  )))
  bmasem_test_hist(stage_2)
  bmasem_test_trace(stage_2)
  expect_error(
    print_out <- capture_output(pp_summary(stage_2), width = 300),
    NA
  )
  bmasem_test_pp_shared(print_out, method)
  expect_error(
    print_out <- capture_output(
      pp_summary(stage_2, simple = FALSE),
      width = 300
    ),
    NA
  )
  bmasem_test_pp_shared(print_out, method, FALSE)
})
