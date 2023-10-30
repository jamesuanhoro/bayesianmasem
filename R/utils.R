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
  r_mat <- stats::cov2cor(r_mat)
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

#' Type hash function
#' @description A function that swaps type from string to integer
#' and vice-versa
#' @param elaborate (LOGICAL) If TRUE, print full names, otherwise
#' print abbreviations
#' @inheritParams .converter_helper
#' @returns If search_term is integer, returns string and vice-versa
#' @keywords internal
.type_hash <- function(search_term = NULL, elaborate = FALSE) {
  list_types <- c(
    "fe" = 1,
    "re" = 2,
    "dep" = 3
  )
  if (isTRUE(elaborate)) {
    # Only use when feeding integers
    list_types <- c(
      "Fixed-effects" = 1,
      "Random-effects" = 2,
      "Dependent-samples" = 3
    )
  }
  converted_value <- .converter_helper(search_term, list_types)
  return(converted_value)
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
.init_fx <- function(data_list) {
  function() {
    list(
      resids = rep(0, (data_list$method < 90) * ncol(data_list$r_obs_vec)),
      loadings_complex = rep(
        0,
        data_list$complex_struc * sum(
          data_list$loading_pattern == 0 & data_list$loading_fixed == -999
        )
      )
    )
  }
}

#' mini-check for estimation method
#' @description mini-check for estimation method
#' @inheritParams bmasem
#' @returns NULL
#' @keywords internal
.check_bm_method <- function(method) {
  accepted_methods <- .method_hash()
  err_msg <- paste0(
    "method must be one of the following: ",
    paste0("\"", accepted_methods, "\"", collapse = ", ")
  )
  if (!tolower(method) %in% tolower(accepted_methods)) stop(err_msg)
}

#' mini-check for data inputs
#' @description mini-check for data inputs
#' @inheritParams bmasem
#' @returns NULL
#' @keywords internal
.check_bm_data <- function(data, group, sample_cov, sample_nobs) {
  if (
    (is.null(data) || is.null(group)) &&
      (is.null(sample_cov) || is.null(sample_nobs))
  ) {
    stop(paste0(
      "User must provide either:\n\t",
      "(i) a dataset and group variable or\n\t",
      "(ii) sample covariance matrices and sample sizes"
    ))
  }
}

#' mini-check for type
#' @description mini-check for type
#' @inheritParams bmasem
#' @returns NULL
#' @keywords internal
.check_bm_type <- function(type) {
  accepted_types <- .type_hash()
  err_msg <- paste0(
    "type must be one of the following: ",
    paste0(
      "\"", accepted_types, "\"",
      collapse = ", "
    )
  )
  if (is.null(type)) stop(err_msg)
  if (!tolower(type) %in% accepted_types) stop(err_msg)
}

#' Check user input function
#' @description A function that checks user input for adequacy
#' and fails on inadequate input.
#' @param type string for type of check
#' @param object_1 Object to check
#' @param object_2 Object to check
#' @param object_3 Object to check
#' @param object_4 Object to check
#' @returns NULL
#' @keywords internal
.user_input_check <- function(
    type,
    object_1 = NULL,
    object_2 = NULL,
    object_3 = NULL,
    object_4 = NULL) {
  if (type == "model") {
    stopifnot("Model cannot be null" = !is.null(object_1))
  }

  if (type == "priors") {
    stopifnot(
      "See ?new_bmasempriors for how to set up priors." =
        inherits(object_1, "bmasempriors")
    )
  }

  if (type == "method") .check_bm_method(object_1)

  if (type == "data") {
    .check_bm_data(object_1, object_2, object_3, object_4)
  }

  if (type == "type") .check_bm_type(object_1)

  if (type == "cluster" && object_1 == "dep") {
    stopifnot(
      "supply cluster information when type = \"dep\"" = !is.null(object_2)
    )
  }

  return(NULL)
}

#' Posterior summary helper function
#' @description A function that slightly modifies the default summary function
#' in posterior package.
#' @param stan_fit Fitted Stan object
#' @param variable Variable(s) to search for in Stan
#' @param interval Confidence interval to select
#' @param major If TRUE, add some preamble for printing the major
#' parameters table.
#' @returns Summary of posterior draws
#' @keywords internal
.bmasem_post_sum <- function(stan_fit, variable, interval = .9, major = FALSE) {
  draws <- posterior::subset_draws(
    posterior::as_draws(stan_fit),
    variable = variable
  )

  lo_lim <- (1.0 - interval) / 2.0
  up_lim <- 1.0 - lo_lim # nolint
  result <- as.data.frame(posterior::summarise_draws(
    draws, "mean", "median", "sd", "mad",
    ~ quantile(.x, probs = c(lo_lim, up_lim), na.rm = TRUE),
    posterior::default_convergence_measures()
  ))

  if (isTRUE(major)) {
    result <- cbind(
      variable = result$variable,
      group = "", from = "", op = "", to = "",
      result[, -1]
    )
  }

  return(result)
}

#' Create parameter list for plotting
#' @param data_list Data list object passed to Stan
#' @returns Name of varying model parameters
#' @keywords internal
.get_param_list <- function(data_list) {
  param_structure <- data_list$loading_pattern
  fac_names <- colnames(param_structure)
  ind_names <- rownames(param_structure)

  rms_params <- c()
  if (data_list$method < 90) {
    rms_params <- c("RMSE" = "rms_src")
  }

  rmsea_params <- c()
  if (data_list$type >= 2) {
    rmsea_params <- c("overall RMSEA" = "rmsea_mn")
  }
  if (data_list$type == 3) {
    rmsea_params <- c(
      rmsea_params,
      "between RMSEA" = "rmsea_be",
      "within RMSEA" = "rmsea_wi", "% between" = "prop_be"
    )
  }
  params <- c(rms_params, rmsea_params)

  phi_idxs <- which(
    lower.tri(data_list$corr_mask) & data_list$corr_mask == 1,
    arr.ind = TRUE
  )
  if (nrow(phi_idxs) > 0) {
    phi_params <- paste0("phi_mat[", apply(
      phi_idxs, 1, paste0,
      collapse = ","
    ), "]")
    names(phi_params) <- apply(phi_idxs, 1, function(x) {
      paste0(fac_names[x[1]], "~~", fac_names[x[2]])
    })
    params <- c(params, phi_params)
  }

  load_idxs <- which(
    data_list$loading_pattern >= ifelse(data_list$complex_struc == 1, -999, 1) &
      data_list$loading_fixed == -999,
    arr.ind = TRUE
  )
  if (nrow(load_idxs) > 0) {
    load_params <- paste0(
      "Load_mat[", apply(load_idxs, 1, paste0, collapse = ","), "]"
    )
    names(load_params) <- apply(load_idxs, 1, function(x) {
      paste0(fac_names[x[2]], "=~", ind_names[x[1]])
    })
    params <- c(params, load_params)
  }

  if (data_list$Nce > 0) {
    rc_idxs <- data_list$error_mat
    rc_params <- paste0("res_cor[", seq_len(data_list$Nce), "]")
    names(rc_params) <- apply(rc_idxs, 1, function(x) {
      paste0(ind_names[x[1]], "~~", ind_names[x[2]])
    })
    params <- c(params, rc_params)
  }

  return(params)
}

#' Create major parameters helper function
#' @description A function that creates the table of major parameters
#' @param stan_fit Fitted Stan object
#' @param data_list Data list object passed to Stan
#' @param interval Confidence interval to select
#' @returns Summary of posterior draws
#' @keywords internal
.create_major_params <- function(stan_fit, data_list, interval = .9) {
  indicator_labels <- rownames(data_list$loading_pattern)
  factor_labels <- colnames(data_list$loading_pattern)

  rmsea_params <- c("rmsea_mn", "rmsea_be", "rmsea_wi", "prop_be")
  rmsea_names <- c(
    paste0("RMSEA (", c("Overall", "Between", "Within"), ")"),
    "% dispersion between"
  )
  params <- c("ppp", "rms_src", rmsea_params)
  from_list <- c("PPP", "RMSE", rmsea_names)

  phi_idxs <- which(lower.tri(data_list$corr_mask), arr.ind = TRUE)
  if (nrow(phi_idxs) > 0) {
    phi_idxs <- paste0("phi_mat[", apply(
      phi_idxs, 1, paste0,
      collapse = ","
    ), "]")
    params <- c(params, phi_idxs)
  }

  load_idxs <- paste0("Load_mat[", apply(which(
    data_list$loading_pattern >= ifelse(data_list$complex_struc == 1, -999, 1) |
      data_list$loading_fixed != -999,
    arr.ind = TRUE
  ), 1, paste0, collapse = ","), "]")
  params <- c(params, load_idxs)

  if (data_list$Nce > 0) {
    params <- c(params, "res_cor")
  }

  major_parameters <- .bmasem_post_sum(
    stan_fit = stan_fit, variable = params, interval = interval, major = TRUE
  )

  major_parameters <- .modify_major_params(
    major_parameters,
    which(major_parameters$variable %in% c("ppp", "rms_src")),
    group = "Goodness of fit",
    from = from_list[1:2]
  )

  idxs <- which(major_parameters$variable %in% rmsea_params)
  major_parameters <- .modify_major_params(
    major_parameters, idxs,
    group = "Dispersion between and within clusters", op = "",
    from = rmsea_names
  )

  idxs <- grep("Load\\_mat", major_parameters$variable)
  major_parameters <- .modify_major_params(
    major_parameters, idxs,
    group = "Factor loadings", op = "=~",
    from = factor_labels[as.integer(
      gsub("Load_mat\\[\\d+,|\\]", "", major_parameters[idxs, ]$variable)
    )],
    to = indicator_labels[as.integer(
      gsub("Load_mat\\[|,\\d+\\]", "", major_parameters[idxs, ]$variable)
    )]
  )

  idxs <- grep("res\\_cor", major_parameters$variable)
  major_parameters <- .modify_major_params(
    major_parameters, idxs,
    group = "Error correlations", op = "~~",
    from = indicator_labels[data_list$error_mat[, 1]],
    to = indicator_labels[data_list$error_mat[, 2]]
  )

  idxs <- grep("phi\\_mat", major_parameters$variable)
  major_parameters <- .modify_major_params(
    major_parameters, idxs,
    group = "Inter-factor correlations", op = "~~",
    from = factor_labels[as.integer(
      gsub("phi_mat\\[\\d+,|\\]", "", major_parameters[idxs, ]$variable)
    )],
    to = factor_labels[as.integer(
      gsub("phi_mat\\[|,\\d+\\]", "", major_parameters[idxs, ]$variable)
    )]
  )

  major_parameters$ess_bulk <- round(major_parameters$ess_bulk, 1)
  major_parameters$ess_tail <- round(major_parameters$ess_tail, 1)

  major_parameters <- major_parameters[, -1]

  return(major_parameters)
}

#' Modify major parameters table helper function
#' @description A function that adds user friendly descriptions to the
#' major parameters table
#' @param major_parameters Major paramters table
#' @param idxs Relevant rows indexes
#' @param group Parameter group
#' @param op lavaan style operator
#' @param from Variable from
#' @param to Variable to
#' @returns Updated major parameters table
#' @keywords internal
.modify_major_params <- function(
    major_parameters,
    idxs,
    group = "",
    op = "",
    from = "",
    to = "") {
  result <- major_parameters

  if (length(idxs) > 0) {
    result[idxs, ]$group <- group
    result[idxs, ]$op <- op
    result[idxs, ]$from <- from
    result[idxs, ]$to <- to
  }

  return(result)
}
