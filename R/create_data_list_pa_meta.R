#' Stan data helper function
#' @description A function that creates data list object passed to Stan
#' @param lavaan_object lavaan fit object of corresponding model
#' @param partab lavaan parameter table output with ceq.simple & std.lv = TRUE
#' @param old_data original list of matrices, needed for PA
#' @inheritParams bmasem
#' @returns Data list object used in fitting Stan model
#' @keywords internal
.create_data_list_pa_meta <- function(
    lavaan_object = NULL,
    method = "normal",
    type = "re",
    priors = NULL,
    cluster = NULL,
    partab = NULL,
    x_mat = NULL,
    conditional_re = TRUE,
    old_data = NULL) {
  data_list <- list()

  # Get number of groups
  data_list$Ng <- lavaan::lavInspect(lavaan_object, "ngroups")

  if (data_list$Ng < 2) {
    stop(paste0(
      "Only one group found. ",
      "Meta-analysis requires multiple samples."
    ))
  }

  # Retrieve parameter structure from lavaan
  param_structure <- lavaan::lavInspect(lavaan_object)[[1]]

  # Set method
  data_list$method <- .method_hash(method)

  # Set type
  data_list$type <- 2 # assume random-effects by default
  data_list$Nc <- 0 # number of clusters
  data_list$C_ID <- vector("integer", data_list$Ng) # cluster ID
  if (type == "fe") {
    data_list$type <- 1
  } else if (type == "dep") {
    data_list$type <- 3
    data_list$C_ID <- cluster
    data_list$Nc <- max(data_list$C_ID)
  }

  # Set simple structure to 0 by default, change within CFA section
  data_list$complex_struc <- 0

  methods::validObject(priors) # validate priors
  # Shape parameter for LKJ of interfactor corr
  data_list$shape_phi_c <- priors@lkj_shape
  data_list$rm_par <- priors@rm_par # sigma(tau) parameter
  data_list$rs_par <- priors@rs_par # sigma(res-sd) parameter
  data_list$rc_par <- priors@rc_par # residual corr parameter
  data_list$rm_i_l_par <- priors@mr_par # meta-reg location(intercept)
  data_list$rm_i_s_par <- priors@sr_par # meta-reg scale(intercept)
  data_list$rm_b_s_par <- priors@br_par # meta-reg scale(beta)

  # Sample size
  data_list$Np <- lavaan::lavInspect(lavaan_object, "nobs")
  # Sample cov
  data_list$S_new <- lapply(lavaan::lavInspect(
    lavaan_object, "SampStat"
  ), "[[", "cov")
  new_var_names <- rownames(data_list$S_new[[1]])
  data_list$S <- lapply(old_data, \(x) {
    x[new_var_names, new_var_names]
  })
  # Number of items
  data_list$Ni <- nrow(data_list$S[[1]])

  # Loading pattern, 0s and 1s
  data_list$loading_pattern <- (param_structure$lambda > 0) * 1

  # Assume CFA by default
  data_list$pa_indicator <- 1
  data_list$sem_indicator <- 0

  # Loading pattern
  data_list$loading_pattern <- param_structure$lambda
  # Number of factors
  data_list$Nf <- ncol(data_list$loading_pattern)
  # location(loading) parameter
  data_list$load_est <- matrix(
    priors@ml_par, data_list$Ni, data_list$Nf,
    dimnames = dimnames(data_list$loading_pattern)
  )
  # sigma(loading) parameter
  data_list$load_se <- matrix(
    priors@sl_par, data_list$Ni, data_list$Nf,
    dimnames = dimnames(data_list$loading_pattern)
  )
  # get fixed loadings
  data_list$loading_fixed <- matrix(
    -999, data_list$Ni, data_list$Nf,
    dimnames = dimnames(data_list$loading_pattern)
  )
  fix_load <- partab[partab$op == "=~" & partab$free == 0, ]
  fix_col_ids <- unname(sapply(
    fix_load$lhs, function(x) which(x == colnames(param_structure$lambda))
  ))
  fix_row_ids <- unname(sapply(
    fix_load$rhs, function(x) which(x == rownames(param_structure$lambda))
  ))
  for (i in seq_len(length(fix_col_ids))) {
    data_list$loading_fixed[
      fix_row_ids[i], fix_col_ids[i]
    ] <- fix_load$ustart[i]
  }

  all_zero_loadings <- all(data_list$loading_pattern == 0)
  all_zero_loadings_fix <- all(data_list$loading_fixed == -999)
  if (isFALSE(all_zero_loadings) || isFALSE(all_zero_loadings_fix)) {
    err_msg <- paste0(
      "This model has latent variables, ",
      "this function is only for path analysis"
    )
    stop(err_msg)
  }

  # res-var
  theta_var_diag <- diag(param_structure$psi)
  data_list$res_var_pattern <- theta_var_diag
  theta_zeroes <- theta_var_diag != 0
  if (sum(theta_zeroes) > 0) {
    theta_var_diag[theta_zeroes] <-
      theta_var_diag[theta_zeroes] - min(theta_var_diag[theta_zeroes]) + 1
    data_list$res_var_pattern <- theta_var_diag
  }
  data_list$res_var_fixed <- array(999, data_list$Ni)
  ind_names <- rownames(param_structure$psi)
  fix_rv <- partab[
    partab$op == "~~" & partab$free == 0 &
      partab$lhs %in% ind_names & partab$lhs == partab$rhs,
  ]
  fix_ind_ids <- unname(sapply(
    fix_rv$lhs, function(x) which(x == ind_names)
  ))
  for (i in seq_len(length(fix_ind_ids))) {
    data_list$res_var_fixed[fix_ind_ids[i]] <- fix_rv$ustart[i]
  }

  # Set up coefficient table
  data_list$coef_pattern <- matrix(0, data_list$Ni, data_list$Ni)
  data_list$coef_est <- matrix(0, data_list$Ni, data_list$Ni)
  data_list$coef_se <- matrix(priors@sc_par, data_list$Ni, data_list$Ni)
  data_list$coef_fixed <- matrix(-999, data_list$Ni, data_list$Ni)
  if (!is.null(param_structure$beta)) {
    # This is an SEM
    data_list$sem_indicator <- 1
    # Factor coefficient matrix
    data_list$coef_pattern <- param_structure$beta
    dimnames(data_list$coef_est) <- dimnames(data_list$coef_pattern)
    dimnames(data_list$coef_se) <- dimnames(data_list$coef_pattern)
    dimnames(data_list$coef_fixed) <- dimnames(data_list$coef_pattern)
    beta_zeroes <- data_list$coef_pattern != 0
    data_list$coef_pattern[beta_zeroes] <-
      data_list$coef_pattern[beta_zeroes] -
      min(data_list$coef_pattern[beta_zeroes]) + 1
    # get fixed coefs
    fix_coef <- partab[partab$op == "~" & partab$free == 0, ]
    fix_col_ids <- unname(sapply(
      fix_coef$lhs, function(x) which(x == colnames(param_structure$lambda))
    ))
    fix_row_ids <- unname(sapply(
      fix_coef$rhs, function(x) which(x == colnames(param_structure$lambda))
    ))
    for (i in seq_len(length(fix_col_ids))) {
      data_list$coef_fixed[
        fix_row_ids[i], fix_col_ids[i]
      ] <- fix_coef$ustart[i]
    }
  } else {
    err_msg <- "Model does not specify any coefficients. Not a path analysis"
    stop(err_msg)
  }

  # Check for correlated error terms
  # Number of correlated errors
  data_list$error_mat <- matrix(ncol = 2, nrow = 0)
  theta_corr_mat <- param_structure$psi
  diag(theta_corr_mat) <- 0
  theta_zeroes <- theta_corr_mat != 0
  if (sum(theta_zeroes) > 0) {
    theta_corr_mat[theta_zeroes] <-
      theta_corr_mat[theta_zeroes] - min(theta_corr_mat[theta_zeroes]) + 1
    data_list$error_mat <- which(theta_zeroes, arr.ind = TRUE)
    data_list$error_mat <- data_list$error_mat[
      data_list$error_mat[, 1] > data_list$error_mat[, 2], ,
      drop = FALSE
    ]
  }
  data_list$error_pattern <- theta_corr_mat
  data_list$Nce <- nrow(data_list$error_mat)

  # Check for conditional independence
  cond_ind_list <- dagitty::impliedConditionalIndependencies(
    dagitty::lavaanToGraph(partab)
  )
  ind_names <- rownames(param_structure$psi)
  data_list$cond_ind_mat <- matrix(
    0,
    nrow = nrow(param_structure$psi), ncol = ncol(param_structure$psi),
    dimnames = dimnames(param_structure$psi)
  )
  data_list$cond_ind_details <- matrix(nrow = 0, ncol = 2 + data_list$Ni)
  if (length(cond_ind_list) > 0) {
    for (i in seq_along(cond_ind_list)) {
      ci <- cond_ind_list[[i]]
      x_ind <- which(ind_names == ci$X)
      y_ind <- which(ind_names == ci$Y)
      if (data_list$error_pattern[x_ind, y_ind] == 0) {
        ci_row <- vector("integer", 2 + data_list$Ni)
        ci_row[1] <- x_ind
        ci_row[2] <- y_ind
        ci_row[3:length(ci_row)] <- ifelse(ind_names %in% ci$Z, 1, 0)
        data_list$cond_ind_details <- rbind(data_list$cond_ind_details, ci_row)
        data_list$cond_ind_mat[x_ind, y_ind] <- 1
        data_list$cond_ind_mat[y_ind, x_ind] <- 1
        # keep using cond_ind_mat to know number of cis,
        # ci_row repeats unique relations
      }
    }
  }

  data_list$correlation <- 1

  data_list$S <- lapply(data_list$S, \(x) {
    suppressWarnings(x <- stats::cov2cor(x))
    x
  })
  data_list$r_mat_s <- lapply(data_list$S, \(x) {
    x[is.na(x)] <- 999
    x
  })

  # For now, no moderators
  data_list$X <- .process_x_mat(
    x_mat, data_list$Ng, data_list$type
  )
  data_list$p <- ncol(data_list$X)
  data_list$X_c <- matrix(nrow = 0, ncol = 0)
  data_list$p_c <- ncol(data_list$X_c)

  data_list$conditional_re <- as.integer(isTRUE(conditional_re))

  data_list$compute_r <- 1
  ni_sq <- (data_list$Ni * (data_list$Ni - 1)) %/% 2
  data_list$r_obs_vec_inp <- array(0, dim = c(data_list$Ng, ni_sq))
  data_list$r_obs_vec_cov_inp <- array(0, dim = c(data_list$Ng, ni_sq, ni_sq))

  data_list <- data_list[!names(data_list) %in% c("S", "S_new")]

  return(data_list)
}
