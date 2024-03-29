#' Stan data helper function
#' @description A function that creates data list object passed to Stan
#' @param lavaan_object lavaan fit object of corresponding model
#' @param partab lavaan parameter table output with ceq.simple & std.lv = TRUE
#' @param acov_mat (Optional) The asymptotic variance matrix of
#' lower triangular half (column-order) of the correlation matrix
#' to be used for correlation structure analysis.
#' @param old_names (Optional) Variable name order of original
#' correlation matrix, used to reorder acov_mat.
#' @inheritParams bmasem
#' @returns Data list object used in fitting Stan model
#' @keywords internal
.create_data_list_pooled <- function(
    lavaan_object = NULL,
    method = "normal",
    simple_struc = TRUE,
    priors = NULL,
    partab = NULL,
    acov_mat = NULL,
    old_names = NULL) {
  data_list <- list()

  # Get number of groups
  data_list$Ng <- 1

  # Retrieve parameter structure from lavaan
  param_structure <- lavaan::lavInspect(lavaan_object)

  # Set method
  data_list$method <- .method_hash(method)

  # Set type
  data_list$type <- 1 # use fixed-effects
  data_list$Nc <- 0 # no clusters
  data_list$C_ID <- vector("integer", 1) # cluster ID

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
  data_list$Np <- array(lavaan_object@SampleStats@nobs[[1]], 1)
  # Sample cov
  data_list$S <- lavaan_object@SampleStats@cov[[1]]
  # Number of items
  data_list$Ni <- nrow(data_list$S)

  data_list$correlation <- 1
  data_list$S <- stats::cov2cor(data_list$S)
  r_vec <- data_list$S[lower.tri(data_list$S, diag = FALSE)]
  ni_sq <- length(r_vec)
  data_list$r_obs_vec <- array(g_map(r_vec), dim = c(1, ni_sq))
  data_list$r_obs_vec_cov <- array(dim = c(1, ni_sq, ni_sq))
  tmp_acov <- acov_mat
  if (!isTRUE(all.equal(old_names, rownames(param_structure$lambda)))) {
    tmp_acov <- .fix_acov(
      acov_mat, old_names, rownames(param_structure$lambda)
    )
  }
  data_list$r_obs_vec_cov[1, , ] <- get_avar_mat(
    data_list$S, data_list$Np[1], tmp_acov
  )
  data_list$S <- NULL

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
  data_list$error_mat <- matrix(nrow = 0, ncol = 2)
  theta_corr_mat <- param_structure$theta
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

  # For now, no moderators
  data_list$p <- 0
  data_list$X <- matrix(nrow = 1, ncol = data_list$p)
  data_list$p_c <- 0
  data_list$X_c <- matrix(nrow = 0, ncol = 0)

  data_list$conditional_re <- 0

  return(data_list)
}
