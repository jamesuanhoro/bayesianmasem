# Helper functions in package

#' Package on attach message
#' @returns R version minimum as string
#' @keywords internal
.onAttach <- function(libname, pkgname) {
  version <- .bayesianmasem_version()
  packageStartupMessage(" ")
  packageStartupMessage(strrep("#", 79))
  packageStartupMessage("This is ", paste(pkgname, version))
  packageStartupMessage(
    "\nAll users of R (or SEM) are invited to report bugs, ",
    "submit functions or ideas\nfor functions. ",
    "An efficient way to do this is to open an issue on GitHub\n",
    "https://github.com/jamesuanhoro/bayesianmasem/issues/."
  )
  packageStartupMessage(strrep("#", 79))
}

#' Get package version function
#' @returns Package version as string
#' @keywords internal
.bayesianmasem_version <- function() {
  version <- read.dcf(
    system.file("DESCRIPTION", package = "bayesianmasem"),
    fields = "Version"
  )[1]
  return(version)
}

#' Log matrix function using eigendecomposition
#' @param r_mat (matrix) Correlation matrix
#' @returns Logarithm of the correlation matrix
#' @export
log_m <- function(r_mat) {
  eig <- eigen(r_mat, symmetric = TRUE)
  return(eig$vectors %*% diag(log(eig$values)) %*% t(eig$vectors))
}

#' Log matrix function using eigendecomposition
#' @param r_vec (vector) Strict lower half vector of correlation matrix
#' @returns Strict lower half vector of log(correlation matrix)
#' @export
g_map <- function(r_vec) {
  d <- 0.5 * (1 + sqrt(1 + 8 * length(r_vec)))
  r_mat <- matrix(0, d, d)
  r_mat[lower.tri(r_mat, diag = FALSE)] <- r_vec
  r_mat <- r_mat + t(r_mat)
  diag(r_mat) <- 1
  r_log_mat <- log_m(r_mat)
  r_log_vec <- r_log_mat[lower.tri(r_log_mat, diag = FALSE)]
  return(r_log_vec)
}

#' Asymptotic variance matrix of lower half vector of log(correlation matrix)
#' @param n (positive integer) Sample size
#' @param acov_mat (matrix) Asymptotic variance matrix of lower half vector
#' of correlation matrix.
#' @inheritParams log_m
#' @returns Asymptotic variance matrix of strict lower half vector
#' of log(correlation matrix).
#' @export
get_avar_mat <- function(r_mat, n, acov_mat = NULL) {
  r_mat <- cov2cor(r_mat)
  if (is.null(acov_mat)) {
    omega_r <- get_asy_cov(r_mat) / n
  } else {
    omega_r <- acov_mat
  }
  r_vec <- r_mat[lower.tri(r_mat, diag = FALSE)]
  a_mat_inv <- pracma::jacobian(g_map, r_vec)
  omega_gamma <- a_mat_inv %*% omega_r %*% a_mat_inv
  return(omega_gamma)
}

#' Transform unbounded vector to correlation matrix
#' @param gamma_in (vector) Strict lower half vector of log(correlation matrix)
#' @param tol_val (positive real) Tolerance for calculation
#' @param ret_vec (Logical) If TRUE, return strict lower half vector,
#' ELSE: return correlation matrix.
#' @returns Transformed correlation matrix
#' @export
g_inv_map <- function(gamma_in, tol_val = 1e-13, ret_vec = FALSE) {
  mat_c <- matrix(nrow = 0, ncol = 0)
  iter <- -1

  d <- 0.5 * (1 + sqrt(1 + 8 * length(gamma_in)))
  if (!(is.vector(gamma_in)) && d %% 1 == 0) {
    stop("Dimension of 'gamma' is incorrect")
  } else if (!(tol_val >= 1e-14 && tol_val <= 1e-4)) {
    stop("Incorrect tolerance value")
  } else {
    mat_a <- matrix(0, d, d)
    mat_a[lower.tri(mat_a, diag = FALSE)] <- gamma_in
    mat_a <- mat_a + t(mat_a)

    diag_a <- diag(mat_a)

    dist <- sqrt(d)
    while (dist > sqrt(d) * tol_val) {
      diag_delta <- log(diag(expm::expm(mat_a)))
      diag_a <- diag_a - diag_delta
      diag(mat_a) <- diag_a
      dist <- norm(diag_delta, type = "2")
      iter <- iter + 1
    }

    mat_c <- expm::expm(mat_a)
    diag(mat_c) <- rep(1, d)
  }

  ret <- list(C = mat_c, iter = iter)
  if (isTRUE(ret_vec)) {
    ret <- mat_c[lower.tri(mat_c, diag = FALSE)]
  }
  return(ret)
}

#' Asymptotic variance matrix of lower half vector of correlation matrix
#' @inheritParams log_m
#' @returns Asymptotic variance matrix of strict lower half vector
#' of correlation matrix.
#' @export
get_asy_cov <- function(r_mat) {
  ltri_idxs <- which(lower.tri(r_mat), arr.ind = TRUE)
  p_ast <- nrow(ltri_idxs)
  omega <- .omega_computer(r_mat, p_ast, ltri_idxs)
  return(omega)
}

#' Asymptotic variance matrix of lower half vector of correlation matrix
#' @inheritParams log_m
#' @param p_ast (positive integer) Dimension of strict lower half vector
#' @param ltri_idxs (array) Indices of strict lower half vector
#' @returns Asymptotic variance matrix of strict lower half vector
#' of correlation matrix.
#' @keywords internal
.omega_computer <- function(r_mat, p_ast, ltri_idxs) {
  omega <- matrix(0, nrow = p_ast, ncol = p_ast)
  for (col in 1:p_ast) {
    sub_idx <- col:p_ast
    idx_id <- matrix(nrow = length(sub_idx), ncol = 4)
    idx_id[, 1] <- ltri_idxs[col, 1]
    idx_id[, 2] <- ltri_idxs[col, 2]
    idx_id[, 3] <- ltri_idxs[sub_idx, 1]
    idx_id[, 4] <- ltri_idxs[sub_idx, 2]
    res_vec <- .four_idx_solver(r_mat, idx_id)
    omega[col, sub_idx] <- res_vec
    omega[sub_idx, col] <- res_vec
  }
  return(omega)
}

#' An internal solver for omega
#' @inheritParams log_m
#' @param idx_four (array) Contains four-way indices for computing omega.
#' @returns asymptotic covariances
#' @keywords internal
.four_idx_solver <- function(r_mat, idx_four) {
  i <- idx_four[, 1]
  j <- idx_four[, 2]
  k <- idx_four[, 3]
  l <- idx_four[, 4]
  r_s <- matrix(nrow = nrow(idx_four), ncol = 6)
  p <- nrow(r_mat)
  r_s[, 1] <- r_mat[(j - 1) * p + i]
  r_s[, 2] <- r_mat[(k - 1) * p + i]
  r_s[, 3] <- r_mat[(l - 1) * p + i]
  r_s[, 4] <- r_mat[(k - 1) * p + j]
  r_s[, 5] <- r_mat[(l - 1) * p + j]
  r_s[, 6] <- r_mat[(l - 1) * p + k]
  ret <- .5 * r_s[, 1] * r_s[, 6] *
    (r_s[, 2]^2 + r_s[, 3]^2 + r_s[, 4]^2 + r_s[, 5]^2) +
    r_s[, 2] * r_s[, 5] + r_s[, 3] * r_s[, 4] - (
      r_s[, 1] * r_s[, 2] * r_s[, 3] + r_s[, 1] * r_s[, 4] * r_s[, 5] +
        r_s[, 2] * r_s[, 4] * r_s[, 6] + r_s[, 3] * r_s[, 5] * r_s[, 6]
    )
  return(ret)
}

#' Method hash function
#' @description A function that swaps method from string to integer
#' and vice-versa
#' @inheritParams .converter_helper
#' @returns If search_term is integer, returns string and vice-versa
#' @keywords internal
.method_hash <- function(search_term = NULL) {
  # Reserving 90+ for methods with no residuals
  list_methods <- c(
    "normal" = 1,
    "lasso" = 2,
    "logistic" = 3,
    "GDP" = 4,
    "WB" = 90,
    "WB-cond" = 91,
    "WW" = 92,
    "none" = 100
  )
  converted_value <- .converter_helper(search_term, list_methods)
  return(converted_value)
}

#' converter function
#' @description A function that swaps type from string to integer
#' and vice-versa
#' @param search_term Integer or string or NULL
#' @param str_list Character vector of accepted values
#' @returns If search_term is integer, returns string and vice-versa
#' @keywords internal
.converter_helper <- function(search_term, str_list) {
  if (is.null(search_term)) {
    converted_value <- names(str_list)
  } else if (is.integer(search_term) || is.numeric(search_term)) {
    search_term <- as.integer(search_term)
    converted_value <- names(str_list)[which(search_term == str_list)]
  } else if (is.character(search_term)) {
    idx <- which(tolower(search_term) == tolower(names(str_list)))
    converted_value <- as.integer(str_list[idx])
  }

  return(converted_value)
}

#' Initial values for Stan fitter
#' @description A function that provides initial values for Stan fitter
#' @param data_list A data list object
#' @returns A function with initial values
#' @keywords internal
.init_fx <- function(dl) {
  function() {
    list(
      resids = rep(0, (dl$method < 90) * ncol(dl$r_obs_vec)),
      loadings_complex = rep(
        0,
        dl$complex_struc * (
          prod(dim(dl$loading_pattern)) - max(dl$loading_pattern) -
            sum(dl$loading_fixed != -999)
        )
      )
    )
  }
}

#' target fitter function
#' @description A function that takes user input and fits the
#' Stan model.
#' @param data_list Data list object passed to Stan
#' @inheritParams bmasem
#' @returns Fitted Stan model
#' @keywords internal
.target_fitter <- function(
    target,
    data_list,
    seed,
    warmup,
    sampling,
    refresh,
    adapt_delta,
    max_treedepth,
    chains,
    ncores,
    show_messages) {
  init_resid <- .init_fx(data_list)

  if (target == "rstan") {
    suppressWarnings(stan_fit <- rstan::sampling(
      stanmodels$mcfar,
      data = data_list,
      chains = chains,
      cores = ncores,
      seed = seed,
      warmup = warmup,
      iter = warmup + sampling,
      refresh = refresh,
      init = init_resid,
      control = list(
        adapt_delta = adapt_delta,
        max_treedepth = max_treedepth
      ),
      show_messages = show_messages
    ))
    message(rstan::check_hmc_diagnostics(stan_fit))
  } else if (target == "cmdstan") {
    message(paste0(
      "Compiling Stan code ...\n",
      "This takes a while the first time you run a CFA ",
      "and the first time you run an SEM"
    ))

    mcfar <- cmdstanr::cmdstan_model(
      system.file("cmdstan/mcfar.stan", package = "bayesianmasem"),
      stanc_options = list("O1")
    )

    message("Fitting Stan model ...")

    stan_fit <- mcfar$sample(
      data = data_list,
      seed = seed,
      iter_warmup = warmup,
      iter_sampling = sampling,
      refresh = refresh,
      init = init_resid,
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth,
      chains = chains,
      parallel_chains = ncores,
      show_messages = show_messages
    )
  }

  return(stan_fit)
}
