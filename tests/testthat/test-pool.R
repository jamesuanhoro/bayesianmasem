test_that("Random type works for pooling on issp89", {
  print(type <- sample(c("fe", "re"), 1))

  expect_error(fit <- bmasem_pool(
    sample_cov = issp89$data[1:3], sample_nobs = issp89$n[1:3],
    warmup = test_warm, sampling = test_samp, chains = test_chns,
    type = type, refresh = 0, show_messages = FALSE
  ), NA)
  expect_true(all(names(fit) %in% c(
    "fit_indices", "r_mat", "r_mat_cov",
    "data_list", "stan_fit"
  )))
})

test_that("Dependent-samples works for pooling on Norton13", {
  cluster_names <- gsub(
    ",", "", lapply(strsplit(names(Norton13$data), " "), "[[", 1)
  )
  subset <- c(1, 19:21)
  cluster_ids <- as.integer(as.factor(cluster_names[subset]))

  expect_error(fit <- bmasem_pool(
    sample_cov = Norton13$data[subset], sample_nobs = Norton13$n[subset],
    cluster = cluster_ids,
    warmup = test_warm, sampling = test_samp, chains = test_chns,
    type = "dep", refresh = 0, show_messages = FALSE
  ), NA)
  expect_true(all(names(fit) %in% c(
    "fit_indices", "r_mat", "r_mat_cov",
    "data_list", "stan_fit"
  )))
})
