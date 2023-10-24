#' Stan data helper function
#' @description A function that creates data list object passed to Stan
#' @param lavaan_object lavaan fit object of corresponding model
#' @param partab lavaan parameter table output with ceq.simple & std.lv = TRUE
#' @inheritParams bmasem
#' @returns Data list object used in fitting Stan model
#' @keywords internal
.create_data_list_meta <- function(
    lavaan_object = NULL,
    method = "normal",
    type = "re",
    simple_struc = TRUE,
    priors = NULL,
    cluster = NULL,
    correlation = TRUE,
    partab = NULL) {
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
  data_list$rc_par <- priors@rc_par # residual corr parameter
  data_list$rm_i_l_par <- priors@mr_par # meta-reg location(intercept)
  data_list$rm_i_s_par <- priors@sr_par # meta-reg scale(intercept)
  data_list$rm_b_s_par <- priors@br_par # meta-reg scale(beta)

  # Sample cov
  data_list$S <- lapply(lavaan::lavInspect(
    lavaan_object, "SampStat"
  ), "[[", "cov")
  if (isTRUE(correlation)) {
    data_list$S <- lapply(data_list$S, stats::cov2cor)
  } else {
    stop("Only correlation analysis methods are implemented right now.")
  }
  # Number of items
  data_list$Ni <- nrow(data_list$S[[1]])
  # Sample size
  data_list$Np <- lavaan::lavInspect(lavaan_object, "nobs")

  # Loading pattern
  data_list$loading_pattern <- param_structure$lambda
  # Number of factors
  data_list$Nf <- ncol(data_list$loading_pattern)
  # location(loading) parameter
  data_list$load_est <- matrix(
    priors@ml_par, data_list$Ni, data_list$Nf
  )
  # sigma(loading) parameter
  data_list$load_se <- matrix(
    priors@sl_par, data_list$Ni, data_list$Nf
  )
  # get fixed loadings
  data_list$loading_fixed <- matrix(-999, data_list$Ni, data_list$Nf)
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

  # Is this an SEM or a CFA?
  psi_mat <- param_structure$psi
  sum_off_diag_psi <- sum(psi_mat[lower.tri(psi_mat)])
  if (is.null(param_structure$beta)) {
    # This is a CFA
    data_list$sem_indicator <- 0
    data_list$complex_struc <- as.integer(ifelse(isFALSE(simple_struc), 1, 0))
    data_list$corr_mask <- diag(data_list$Nf)
    data_list$corr_mask[lower.tri(data_list$corr_mask)] <-
      (psi_mat[lower.tri(psi_mat)] != 0) + 0
    data_list$corr_mask[upper.tri(data_list$corr_mask)] <-
      t(data_list$corr_mask)[upper.tri(data_list$corr_mask)]
  } else {
    stop("Only CFAs are implemented right now.")
    # This is an SEM
    data_list$sem_indicator <- 1
    # Factor correlation matrix
    data_list$F_corr_mat <- matrix(byrow = TRUE, ncol = 2, nrow = 0)
    if (sum_off_diag_psi > 0) {
      # Get which elements in theta are non-zero
      data_list$F_corr_mat <- which(psi_mat != 0, arr.ind = TRUE)
      # Eliminate diagonal elements
      data_list$F_corr_mat <- data_list$F_corr_mat[
        data_list$F_corr_mat[, 1] != data_list$F_corr_mat[, 2],
      ]
      # Eliminate duplicate rows
      data_list$F_corr_mat <- unique(t(apply(data_list$F_corr_mat, 1, sort)))
    }
    data_list$Nf_corr <- nrow(data_list$F_corr_mat)
    data_list$coef_pattern <- (param_structure$beta > 0) * 1
  }

  # Marker variables per factor
  # Each factor should have one unique indicator or stop!
  data_list$markers <- array(dim = data_list$Nf)
  for (j in seq_len(ncol(data_list$loading_pattern))) {
    data_list$markers[j] <- which(
      data_list$loading_pattern[, j] != 0
    )[1]
  }

  # Check for correlated error terms
  # Number of correlated errors
  data_list$error_mat <- matrix(byrow = TRUE, ncol = 2, nrow = 0)
  theta_mat <- param_structure$theta
  sum_off_diag_theta <- sum(theta_mat[lower.tri(theta_mat)])
  if (sum_off_diag_theta > 0) {
    # Get which elements in theta are non-zero
    data_list$error_mat <- which(theta_mat != 0, arr.ind = TRUE)
    # Eliminate diagonal elements
    data_list$error_mat <- data_list$error_mat[
      data_list$error_mat[, 1] != data_list$error_mat[, 2],
    ]
    # Eliminate duplicate rows
    data_list$error_mat <- unique(t(apply(data_list$error_mat, 1, sort)))
  }
  data_list$Nce <- nrow(data_list$error_mat)

  # For now, no moderators
  data_list$p <- 0
  data_list$X <- matrix(nrow = data_list$Ng, ncol = data_list$p)
  data_list$p_c <- 0
  data_list$X_c <- matrix(nrow = 0, ncol = 0)

  # Fail on missing data
  if (any(unlist(lapply(data_list$S, function(x) sum(is.na(x)) > 0)))) {
    stop("All correlation matrices must have complete data")
  } else {
    data_list$r_obs_vec <- do.call("rbind", lapply(data_list$S, function(s) {
      vec <- stats::cov2cor(s)[lower.tri(s, diag = FALSE)]
      g_map(vec)
    }))
    ni_sq <- ncol(data_list$r_obs_vec)
    data_list$r_obs_vec_cov <- array(dim = c(data_list$Ng, ni_sq, ni_sq))
    for (i in seq_len(data_list$Ng)) {
      data_list$r_obs_vec_cov[i, , ] <- get_avar_mat(
        data_list$S[[i]], data_list$Np[i]
      )
    }
  }

  return(data_list)
}
