// Each S is dist according to some E
// And each E follows some structured matrix
// And that structured matrix is IW (see Wu & Browne)
functions {
  array[] int find_row_zero(matrix sq_mat, array[] int is_row_fixed) {
    int ni = rows(sq_mat);
    int des_zero = ni - sum(is_row_fixed);
    array[ni] int tmp_res = rep_array(0, ni);
    int idx = 0;
    for (i in 1:ni) {
      if (is_row_fixed[i] == 0) {
        int temp = 0;
        for (j in 1:ni) {
          if (is_row_fixed[j] == 0) {
            if (sq_mat[i, j] == 0) temp += 1;
          }
        }
        if (temp == des_zero) {
          idx += 1;
          tmp_res[idx] = i;
        }
      }
    }
    array[idx] int result = tmp_res[1:idx];
    return(result);
  }
  array[,] int find_recursive_set(matrix coef_mat) {
    int ni = rows(coef_mat);
    array[ni, ni] int result = rep_array(0, ni, ni);
    array[ni] int fix_variable = rep_array(0, ni);
    int ni_sofar = 0;
    int i = 0;
    while (ni_sofar < ni) {
      array[size(find_row_zero(coef_mat, fix_variable))] int temp = find_row_zero(coef_mat, fix_variable);
      fix_variable[temp] = rep_array(1, size(temp));
      i += 1;
      result[i, 1:size(temp)] = temp;
      ni_sofar += size(temp);
    }
    return(result);
  }
  vector find_factor_res_var(matrix coef_mat, matrix cor_psi) {
    int ni = rows(coef_mat);
    array[ni, ni] int full_set = find_recursive_set(coef_mat);
    vector[ni] error_var = rep_vector(1.0, ni);
    vector[ni] total_var_psi = rep_vector(1.0, ni);
    array[ni] int iv = rep_array(0, ni);

    int n_set_1 = 0;
    for (i in 1:ni) {
      if (full_set[1, i] != 0) n_set_1 += 1;
    }
    array[n_set_1] int idx_set_1 = rep_array(0, n_set_1);
    {
      int counter = 0;
      for (i in 1:ni) {
        if (full_set[1, i] != 0) {
          counter += 1;
          idx_set_1[counter] = full_set[1, i];
        }
      }
    }
    error_var[idx_set_1] = total_var_psi[idx_set_1];

    matrix[n_set_1, n_set_1] ic_cor = cor_psi[idx_set_1, idx_set_1];
    vector[n_set_1] start_var = total_var_psi[idx_set_1];
    matrix[ni, ni] iv_cov;
    iv_cov[1:n_set_1, 1:n_set_1] = quad_form_diag(ic_cor, sqrt(start_var));

    int n_sets = 0;
    for (i in 1:ni) {
      if (sum(full_set[i, ]) > 0) n_sets += 1;
    }
    array[n_sets, ni] int set = full_set[1:n_sets, ];

    int n_set_start;
    int n_set_end = 0;
    for (i in 1:(n_sets - 1)) {
      n_set_start = n_set_end + 1;
      int n_set_i = 0;
      int n_set_ip1 = 0;
      for (j in 1:ni) {
        if (set[i, j] != 0) n_set_i += 1;
        if (set[i + 1, j] != 0) n_set_ip1 += 1;
      }
      n_set_end += n_set_i;
      array[n_set_i] int idx_set_i = rep_array(0, n_set_i);
      array[n_set_ip1] int idx_set_ip1 = rep_array(0, n_set_ip1);
      {
        int counter_i = 0;
        int counter_ip1 = 0;
        for (j in 1:ni) {
          if (set[i, j] != 0) {
            counter_i += 1;
            idx_set_i[counter_i] = set[i, j];
          }
          if (set[i + 1, j] != 0) {
            counter_ip1 += 1;
            idx_set_ip1[counter_ip1] = set[i + 1, j];
          }
        }
      }
      iv[n_set_start:n_set_end] = idx_set_i;

      matrix[n_set_ip1, n_set_end] tmp_beta = coef_mat[idx_set_ip1, iv[1:n_set_end]];

      matrix[n_set_ip1, n_set_ip1] var_reg = quad_form_sym(
        iv_cov[1:n_set_end, 1:n_set_end], tmp_beta'
      );

      matrix[n_set_ip1, n_set_ip1] temp_psi = cor_psi[idx_set_ip1, idx_set_ip1];
      vector[n_set_ip1] temp_psi_sd = rep_vector(0.0, n_set_ip1);

      for (j in 1:n_set_ip1) {
        error_var[idx_set_ip1[j]] = total_var_psi[idx_set_ip1[j]] - var_reg[j, j];
        temp_psi_sd[j] = fmax(0.0, sqrt(error_var[idx_set_ip1[j]]));
      }

      if (i < (n_sets - 1)) {
        temp_psi = quad_form_diag(temp_psi, temp_psi_sd);
        int n_agg = n_set_end + n_set_ip1;
        matrix[n_agg, n_agg] real_temp_psi = rep_matrix(0, n_agg, n_agg);
        real_temp_psi[1:n_set_end, 1:n_set_end] = iv_cov[1:n_set_end, 1:n_set_end];
        real_temp_psi[(n_set_end + 1):n_agg, (n_set_end + 1):n_agg] = temp_psi;
        array[n_agg] int agg;
        agg[1:n_set_end] = iv[1:n_set_end];
        agg[(n_set_end + 1):n_agg] = idx_set_ip1;
        matrix[n_agg, n_agg] temp_path2 = rep_matrix(0, n_agg, n_agg);
        temp_path2[(n_set_end + 1):n_agg, ] = coef_mat[idx_set_ip1, agg];
        matrix[n_agg, n_agg] id_mat = identity_matrix(n_agg);
        matrix[n_agg, n_agg] id_tmp2_inv = inverse(id_mat - temp_path2);
        iv_cov[1:n_agg, 1:n_agg] = quad_form(real_temp_psi, id_tmp2_inv');
      }
    }

    return(error_var);
  }
  real generalized_double_pareto_lpdf(vector x, real alpha, real scale) {
    // generalized double Pareto
    // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3903426/
    return(sum(
      -log(2) - log(scale) - (alpha + 1.0) * log1p(abs(x) / (scale * alpha))
    ));
  }
  matrix cov2cor (matrix C_mat) {
    int p = dims(C_mat)[1];
    vector[p] S_i = 1 ./ sqrt(diagonal(C_mat));
    matrix[p, p] R_mat = quad_form_diag(C_mat, S_i);
    for (i in 1:p) R_mat[i, i] = 1;
    return R_mat;
  }
  int sign(real x) {
    if (x > 0)
      return 1;
    else
      return -1;
  }
  vector four_idx_solver(
    matrix r_mat, array[] int i, array[] int j, array[] int k, array[] int l
  ) {
    int n = size(i);
    vector[n] ret;

    for (t in 1:n) {
      real rij = r_mat[i[t], j[t]];
      real rik = r_mat[i[t], k[t]];
      real ril = r_mat[i[t], l[t]];
      real rjk = r_mat[j[t], k[t]];
      real rjl = r_mat[j[t], l[t]];
      real rkl = r_mat[k[t], l[t]];
      ret[t] =
        0.5 * rij * rkl *
        (square(rik) + square(ril) + square(rjk) + square(rjl))
        + rik * rjl + ril * rjk - (
          rij * rik * ril + rij * rjk * rjl + rik * rjk * rkl + ril * rjl * rkl
        );
    }

    return(ret);
  }
  matrix omega_computer(matrix r_mat) {
    int p = rows(r_mat);
    int d = (p * (p - 1)) %/% 2;
    matrix[d, d] omega;

    array[d] int ii;
    array[d] int jj;

    int pos = 0;
    for (j in 1:(p - 1)) {
      for (i in (j + 1):p) {
        pos += 1;
        ii[pos] = i;
        jj[pos] = j;
      }
    }

    omega = rep_matrix(0, d, d);

    for (col in 1:d) {
      int m = d - col + 1;

      array[m] int i;
      array[m] int j;
      array[m] int k;
      array[m] int l;

      for (t in 1:m) {
        i[t] = ii[col];
        j[t] = jj[col];
        k[t] = ii[col + t - 1];
        l[t] = jj[col + t - 1];
      }

      vector[m] vals = four_idx_solver(r_mat, i, j, k, l);

      for (t in 1:m) {
        omega[col, col + t - 1] = vals[t];
        omega[col + t - 1, col] = vals[t];
      }
    }

    return omega;
  }
  matrix frechet_log_apply(matrix r_mat, matrix h_mat) {
    int p = rows(r_mat);
    vector[p] lambda;
    matrix[p, p] q_vecs;

    (q_vecs, lambda) = eigendecompose_sym(r_mat);

    vector[p] ln_lambda = log(lambda);
    matrix[p, p] f_mat;

    for (i in 1:p) {
      for (j in 1:p) {
        if (i == j) {
          f_mat[i, j] = 1.0 / lambda[i];
        } else {
          f_mat[i,j] = (ln_lambda[i] - ln_lambda[j]) / (lambda[i] - lambda[j]);
        }
      }
    }

    matrix[p, p] b_mat;

    b_mat = quad_form(h_mat, q_vecs);
    b_mat = f_mat .* b_mat;

    return(quad_form(b_mat, q_vecs'));
  }
  vector vechs(matrix mat) {
    int p = rows(mat);
    int d = (p * (p - 1)) %/% 2;
    vector[d] v;

    int k = 0;
    for (j in 1:(p - 1)) {
      for (i in (j + 1):p) {
        k += 1;
        v[k] = mat[i, j];
      }
    }
    return(v);
  }
  matrix get_jacob(matrix r_mat) {
    int p = rows(r_mat);
    int d = (p * (p - 1)) %/% 2;
    matrix[d, d] jacob_mat;

    for (k in 1:d) {
      matrix[p, p] h_mat = rep_matrix(0, p, p);
      // fill lower triangle basis
      int idx = 1;
      for (j in 1:(p - 1) ) {
        for (i in (j + 1):p) {
          if (idx == k) {
            h_mat[i, j] = 1;
            h_mat[j, i] = 1;
          }
          idx += 1;
        }
      }

      matrix[p, p] l_mat = frechet_log_apply(r_mat, h_mat);

      jacob_mat[, k] = vechs(l_mat);
    }

    return(jacob_mat);
  }
  matrix gamma_avar(matrix omega, matrix jacob) {
    int p = rows(omega);

    matrix[p, p] omega_gamma = quad_form(omega, jacob');
    omega_gamma = (omega_gamma + omega_gamma') / 2;

    return(omega_gamma);
  }
  vector matrix_log_vech(matrix A) {
    int p = rows(A);
    int p_ast = (p * (p - 1)) %/% 2;
    vector[p] lambda_tmp = eigenvalues_sym(A);
    matrix[p, p] q_mat_tmp = eigenvectors_sym(A);
    vector[p] ln_lambda;
    matrix[p, p] q_mat;
    matrix[p, p] ln_A;
    vector[p_ast] ln_A_lower_tri;

    for (i in 1:p) {
      ln_lambda[i] = log(lambda_tmp[p - i + 1]);
      q_mat[, i] = q_mat_tmp[, p - i + 1];
    }

    ln_A = quad_form(diag_matrix(ln_lambda), q_mat');

    {
      int pos = 0;
      for (j in 1:(p - 1)) {
        for (i in (j + 1):p) {
          pos += 1;
          ln_A_lower_tri[pos] = ln_A[i, j];
        }
      }
    }

    return(ln_A_lower_tri);
  }
}
data {
  int<lower = 0> Ng;  // number of groups
  array[Ng] int<lower = 0> Np;  // number persons by matrix
  int<lower = 0> Ni;  // number items
  array[Ng] matrix[Ni, Ni] r_mat_s;  // covariance matrices
  array[Ng] vector[(Ni * (Ni - 1)) %/% 2] r_obs_vec_inp;  // covariance matrices
  array[Ng] matrix[(Ni * (Ni - 1)) %/% 2, (Ni * (Ni - 1)) %/% 2] r_obs_vec_cov_inp;  // covariance matrices
  int p;  // number of moderators
  matrix[Ng, p] X;  // moderator matrix
  real rm_i_l_par;  // meta-reg int loc hyper-parameter -- log(.08)
  real<lower = 0> rm_i_s_par;  // meta-reg int scale hyper-parameter -- .7
  real<lower = 0> rm_b_s_par;  // meta-reg beta scale hyper-parameter -- .5
  int<lower = 0> Nc;  // number of clusters
  array[Ng] int<lower = 0, upper = Nc> C_ID;  // cluster ID
  int<lower = 1, upper = 3> type; // which type
  int<lower = 0, upper = 1> conditional_re;
  int<lower = 0> Nce;  // number correlated errors
  array[Ni, Ni] int error_pattern; // cor error matrix
  array[Ni, Ni] int coef_pattern;  // coef pattern
  array[Ni, Ni] real coef_fixed;  // coef fixed
  real<lower = 0> rm_par;  // rms scale parameter
  real<lower = 1> rc_par;  // residual corr parameter
  int<lower = 1, upper = 100> method; // which method
  array[Ni, Ni] int cond_ind_mat; // conditional independence locations
  matrix[Ni, Ni] coef_est;
  matrix[Ni, Ni] coef_se;
  int<lower = 0, upper = 1> compute_r;
}
transformed data {
  int Nisqd2 = (Ni * (Ni - 1)) %/% 2;
  int N_type_wi = 1;
  int N_type_be = 1;
  array[Ng] vector[Nisqd2] r_obs_vec;  // transformed correlations
  array[Ng] matrix[Nisqd2, Nisqd2] r_obs_vec_cov;  // variance of transformed correlations
  array[Ng] matrix[Nisqd2, Nisqd2] L_vec_cov;
  array[Ng] vector[Nisqd2] r_obs_vec_white;
  array[Ng] matrix[Nisqd2, Nisqd2] L_vec_cov_inv;
  array[Ng] int<lower = 0, upper = 1> missing_var_indicator = rep_array(1, Ng);
  array[Ng] int<lower = 0, upper = 1> missing_element_indicator = rep_array(0, Ng);
  array[Ng] int N_missing_elements_by_study = rep_array(0, Ng);
  int N_missing_elements = 0;
  int N_missing_items = 0;
  // specified with full matrix, but only using sub-matrix indicators!
  array[Ni, Ng] int N_missing_items_array = rep_array(0, Ni, Ng);
  array[Ni, Ng] int presence_indicator = rep_array(0, Ni, Ng);
  array[Ni, Ng] int presence = rep_array(0, Ni, Ng);
  array[Ng] int presence_count = rep_array(0, Ng);

  for (i in 1:Ng) {
    if (max(vechs(r_mat_s[i])) < 900) {
      missing_var_indicator[i] = 0;
      if (compute_r == 1) {
        matrix[Nisqd2, Nisqd2] omega_r = omega_computer(r_mat_s[i]) / (1.0 * Np[i]);
        matrix[Nisqd2, Nisqd2] jacob_mat = get_jacob(r_mat_s[i]);
        r_obs_vec[i] = matrix_log_vech(r_mat_s[i]);
        r_obs_vec_cov[i] = gamma_avar(omega_r, jacob_mat);
      } else {
        r_obs_vec[i] = r_obs_vec_inp[i];
        r_obs_vec_cov[i] = r_obs_vec_cov_inp[i];
      }
      L_vec_cov[i] = cholesky_decompose(r_obs_vec_cov[i]);
      L_vec_cov_inv[i] = inverse(L_vec_cov[i]);
      r_obs_vec_white[i] = L_vec_cov_inv[i] * r_obs_vec[i];
      for (j in 1:Ni) presence_indicator[j, i] = j;
      presence_count[i] = Ni;
    } else {
      for (j in 1:(Ni - 1)) {
        for (k in (j + 1):Ni) {
          if (r_mat_s[i][k, j] < 900) {
            presence[k, i] = 1;
            presence[j, i] = 1;
          }
        }
      }

      for (k in 1:Ni) {
        if (presence[k, i] == 1) {
          presence_count[i] += 1;
          presence_indicator[presence_count[i], i] = k;
        }
      }

      int Ni_i = presence_count[i];
      matrix[Ni_i, Ni_i] sub_mat = r_mat_s[i][
        presence_indicator[1:Ni_i, i],
        presence_indicator[1:Ni_i, i]
      ];

      for (j in 1:(Ni_i - 1)) {
        for (k in (j + 1):Ni_i) {
          if (sub_mat[k, j] > 900) {
            missing_element_indicator[i] = 1;
            N_missing_elements_by_study[i] += 1;
            N_missing_elements += 1;
            N_missing_items_array[k, i] = 1;
            N_missing_items_array[j, i] = 1;
          }
        }
      }

      if (missing_element_indicator[i] == 0) {
        int Nisqd2_i = (Ni_i * (Ni_i - 1)) %/% 2;
        matrix[Nisqd2_i, Nisqd2_i] omega_r = omega_computer(sub_mat) / (1.0 * Np[i]);
        matrix[Nisqd2_i, Nisqd2_i] jacob_mat = get_jacob(sub_mat);
        r_obs_vec[i][1:Nisqd2_i] = matrix_log_vech(sub_mat);
        r_obs_vec_cov[i][1:Nisqd2_i, 1:Nisqd2_i] = gamma_avar(omega_r, jacob_mat);
        L_vec_cov[i][1:Nisqd2_i, 1:Nisqd2_i] = cholesky_decompose(
          r_obs_vec_cov[i][1:Nisqd2_i, 1:Nisqd2_i]
        );
        L_vec_cov_inv[i][1:Nisqd2_i, 1:Nisqd2_i] = inverse(
          L_vec_cov[i][1:Nisqd2_i, 1:Nisqd2_i]
        );
        r_obs_vec_white[i][1:Nisqd2_i] = L_vec_cov_inv[i][1:Nisqd2_i, 1:Nisqd2_i] *
          r_obs_vec[i][1:Nisqd2_i];
      }
    }
    N_missing_items += sum(N_missing_items_array[, i]);
  }

  int N_complete_studies = Ng - sum(missing_var_indicator);

  if (type < 2) N_type_wi = 0;
  if (type < 3) N_type_be = 0;

  real sqrt_two = sqrt(2.0);
  real pi_sqrt_three = pi() / sqrt(3.0);
  int<lower = 0> Nco_uniq = 0;  // N_non-zero coef unique
  int<lower = 0> Nco = 0;  // N_non-zero coef
  int<lower = 0> Nco_fixed = 0;  // N_non-zero coef unique
  int Ncond_ind = 0;
  int N_rms = 1;
  int N_alpha = 0;
  int<lower = 0> Nce_uniq = 0; // N_correlated errors unique

  if (method < 90) {
    for (i in 2:Ni) {
      for (j in 1:(i - 1)) {
        if (cond_ind_mat[i, j] == 1) {
          Ncond_ind += 1;
        }
      }
    }
  }

  if (method == 100) {
    N_rms = 0;
  }

  if (method == 4) N_alpha = 1;

  for (i in 1:Ni) {
    for (j in 1:Ni) {
      if (coef_pattern[i, j] != 0) {
        Nco += 1;
        if (coef_pattern[i, j] > Nco_uniq) Nco_uniq = coef_pattern[i, j];
      } else if (coef_fixed[i, j] > -990) {
        Nco_fixed += 1;
      }
    }
  }

  for (j in 1:(Ni - 1)) {
    for (i in (j + 1):Ni) {
      if (error_pattern[i, j] > Nce_uniq) {
        Nce_uniq = error_pattern[i, j];
      }
    }
  }
}
parameters {
  vector<lower = 0>[N_missing_items] var_shifts;
  vector[N_type_wi] ln_v_int_wi;
  vector[p] ln_v_beta_wi;
  vector[N_type_be] ln_v_int_be;
  array[N_complete_studies * N_type_wi * conditional_re] vector[Nisqd2] c_clus;
  matrix[Nisqd2, Nc] g_clus;
  vector<lower = -1.0, upper = 1.0>[N_missing_elements] cor_missing;
  vector<lower = 0.0, upper = 1.0>[N_rms] rms_src_p;
  vector<lower = 2.0>[N_alpha] gdp_alpha;
  vector[Ncond_ind] resids;  // residual vector
  vector<lower = 0, upper = 1>[Nce_uniq] res_cor_01;  // correlated errors on 01
  vector<lower = -1, upper = 1>[Nco_uniq] coefs;  // may need to be limited to (-1, 1)?
}
transformed parameters {
  real rms_src_tmp = 0.0;

  if (method != 100) rms_src_tmp = rms_src_p[1];
  if (method == 2) {
    rms_src_tmp /= sqrt_two;
  } else if (method == 3) {
    rms_src_tmp /= pi_sqrt_three;
  } else if (method == 4) {
    rms_src_tmp /= sqrt_two * gdp_alpha[1] / sqrt(
      (gdp_alpha[1] - 1.0) * (gdp_alpha[1] - 2.0)
    );
  }
}
model {
  rms_src_p ~ normal(0, rm_par);
  if (method == 1) {
    // normal
    resids ~ std_normal();
  } else if (method == 2) {
    // lasso
    resids ~ double_exponential(0, 1);
  } else if (method == 3) {
    // logistic
    resids ~ logistic(0, 1);
  } else if (method == 4) {
    // https://www.mdpi.com/2075-1680/11/9/462
    gdp_alpha ~ lognormal(1, 1);
    target += generalized_double_pareto_lpdf(
      resids | gdp_alpha[1], 1.0);
  }

  res_cor_01 ~ beta(rc_par, rc_par);

  var_shifts ~ std_normal();

  ln_v_int_wi ~ normal(rm_i_l_par, rm_i_s_par);
  ln_v_int_be ~ normal(rm_i_l_par, rm_i_s_par);
  ln_v_beta_wi ~ normal(0, rm_b_s_par);

  cor_missing ~ normal(0, .5);

  {
    matrix[Ni, Ni] r_mat;
    vector[Nisqd2] r_vec;
    vector[Ni] res_var;
    vector[Nce_uniq] res_cor_u = res_cor_01 * 2 - 1;
    matrix[Ni, Ni] Coef_mat = rep_matrix(0, Ni, Ni);
    matrix[Ni, Ni] One_min_Beta_inv;
    matrix[Ni, Nce] loading_par_exp = rep_matrix(0, Ni, Nce);
    matrix[Ni, Ni] loading_par_exp_2;
    vector[Ni] delta_mat_ast;
    vector[Ni] total_var;
    matrix[Ni, Ni] cor_psi = identity_matrix(Ni);

    for (i in 1:Ni) {
      for (j in 1:Ni) {
        if (coef_pattern[i, j] != 0) {
          coefs[coef_pattern[i, j]] ~ normal(coef_est[i, j], coef_se[i, j]);
          Coef_mat[i, j] = coefs[coef_pattern[i, j]];
        } else if (coef_fixed[i, j] > -990) {
          Coef_mat[i, j] = coef_fixed[i, j];
        }
      }
    }

    One_min_Beta_inv = inverse(identity_matrix(Ni) - Coef_mat);

    {
      int pos_err = 0;
      for (j in 1:(Ni - 1)) {
        for (i in (j + 1):Ni) {
          if (error_pattern[i, j] != 0) {
            pos_err += 1;
            cor_psi[i, j] = res_cor_u[error_pattern[i, j]];
            cor_psi[j, i] = cor_psi[i, j];
          }
        }
      }
    }

    // fix this, error correlation matrix is NaN
    res_var = find_factor_res_var(Coef_mat, cor_psi);

    {
      int pos_err = 0;
      for (j in 1:(Ni - 1)) {
        for (i in (j + 1):Ni) {
          if (error_pattern[i, j] != 0) {
            pos_err += 1;
            loading_par_exp[i, pos_err] = sqrt(
              abs(res_cor_u[error_pattern[i, j]]) * res_var[i]
            );
            loading_par_exp[j, pos_err] =
              sign(res_cor_u[error_pattern[i, j]]) * sqrt(
                abs(res_cor_u[error_pattern[i, j]]) * res_var[j]
              );
          }
        }
      }
    }

    loading_par_exp_2 = tcrossprod(loading_par_exp);
    delta_mat_ast = res_var - diagonal(loading_par_exp_2);

    r_mat = quad_form_sym(
      add_diag(loading_par_exp_2, delta_mat_ast),
      One_min_Beta_inv'
    );

    total_var = diagonal(r_mat);

    if (method < 90) {
      int pos = 0;
      for (i in 2:Ni) {
        for (j in 1:(i - 1)) {
          if (cond_ind_mat[i, j] == 1) {
             pos += 1;
             r_mat[i, j] += resids[pos] * rms_src_tmp * sqrt(total_var[i] * total_var[j]);
             r_mat[j, i] = r_mat[i, j];
          }
        }
      }
    }

    r_vec = matrix_log_vech(r_mat);

    int complete_pos = 0;
    int missing_pos = 0;
    int missing_pos_i = 0;

    if (type == 3) {
      to_vector(g_clus) ~ normal(0, exp(ln_v_int_be[1]));
    }

    for (i in 1:Ng) {
      real m_val;

      if (missing_var_indicator[i] == 0) {
        if (type == 1) {
          target += normal_lupdf(r_obs_vec_white[i] | L_vec_cov_inv[i] * r_vec, 1);
        } else if (type >= 2) {
          m_val = exp(ln_v_int_wi[1] + X[i, ] * ln_v_beta_wi);
          if (conditional_re == 1) {
            complete_pos += 1;
            c_clus[complete_pos] ~ normal(r_vec, m_val);
            if (type == 2) {
              target += normal_lupdf(r_obs_vec_white[i] | L_vec_cov_inv[i] * c_clus[complete_pos], 1);
            } else if (type == 3) {
              target += normal_lupdf(
                r_obs_vec_white[i] | L_vec_cov_inv[i] * (c_clus[complete_pos] + g_clus[, C_ID[i]]), 1
              );
            }
          } else if (conditional_re == 0) {
            if (type == 2) {
              target += multi_normal_cholesky_lupdf(
                r_obs_vec[i] | r_vec,
                cholesky_decompose(add_diag(r_obs_vec_cov[i], square(m_val)))
              );
            } else if (type == 3) {
              target += multi_normal_cholesky_lupdf(
                r_obs_vec[i] | r_vec + g_clus[, C_ID[i]],
                cholesky_decompose(add_diag(r_obs_vec_cov[i], square(m_val)))
              );
            }
          }
        }
      } else {
        int Ni_i = presence_count[i];
        int Nisqd2_i = (Ni_i * (Ni_i - 1)) %/% 2;
        m_val = 0;
        if (type >= 2) m_val = exp(ln_v_int_wi[1] + X[i, ] * ln_v_beta_wi);
        vector[Nisqd2_i] r_obs_vec_i;
        matrix[Nisqd2_i, Nisqd2_i] r_obs_vec_cov_i;
        matrix[Nisqd2_i, Nisqd2_i] L_vec_cov_i;

        // create imputed matrix
        matrix[Ni_i, Ni_i] r_mat_i = r_mat_s[i][
          presence_indicator[1:Ni_i, i],
          presence_indicator[1:Ni_i, i]
        ];

        // get i-specific log-scale vector
        vector[Nisqd2_i] r_vec_i = matrix_log_vech(
          r_mat[
            presence_indicator[1:Ni_i, i],
            presence_indicator[1:Ni_i, i]
          ]
        );

        if (missing_element_indicator[i] == 1) {
          for (j in 1:Ni_i) {
            if (N_missing_items_array[j, i] == 1) {
              missing_pos_i += 1;
              r_mat_i[j, j] += var_shifts[missing_pos_i];
            }
          }

          for (j in 1:(Ni_i - 1)) {
            for (k in (j + 1):Ni_i) {
              if (r_mat_i[k, j] > 900) {
                missing_pos += 1;
                r_mat_i[k, j] = cor_missing[missing_pos] / sqrt(r_mat_i[j, j] * r_mat_i[k, k]);
                r_mat_i[j, k] = r_mat_i[k, j];
              }
            }
          }

          for (j in 1:Ni_i) r_mat_i[j, j] = 1.0;

          r_obs_vec_i = matrix_log_vech(r_mat_i);
          matrix[Nisqd2_i, Nisqd2_i] omega_r = omega_computer(r_mat_i) / (1.0 * Np[i]);
          matrix[Nisqd2_i, Nisqd2_i] jacob_mat = get_jacob(r_mat_i);
          r_obs_vec_cov_i = add_diag(gamma_avar(omega_r, jacob_mat), square(m_val));
        } else {
          r_obs_vec_i = r_obs_vec[i][1:Nisqd2_i];
          r_obs_vec_cov_i = add_diag(r_obs_vec_cov[i][1:Nisqd2_i, 1:Nisqd2_i], square(m_val));
        }

        L_vec_cov_i = cholesky_decompose(r_obs_vec_cov_i);
        // diag has been added already!
        if (type <= 2) {
          target += multi_normal_cholesky_lpdf(r_obs_vec_i | r_vec_i, L_vec_cov_i);
        } else if (type == 3) {
          target += multi_normal_cholesky_lpdf(
            r_obs_vec_i | (r_vec_i + g_clus[1:Nisqd2_i, C_ID[i]]), L_vec_cov_i
          );
        }
      }
    }
  }
}
generated quantities {
  real D_rep = 0.0;
  real D_obs = 0.0;
  real<lower = 0, upper = 1> ppp;
  matrix[Ni, Ni] r_mat;
  vector[Nisqd2] r_vec;
  real v_mn = 0.0;
  real rmsea_mn = sqrt(v_mn);
  real v_wi = 0.0;
  real rmsea_wi = sqrt(v_wi);
  real v_be = 0.0;
  real rmsea_be = sqrt(v_be);
  vector[p] rmsea_beta_wi;
  real prop_be = 0.0;
  vector[Ng] log_lik;
  real<lower = 0> rms_src = 0.0;  // RMSE of residuals
  matrix[Ni, Ni] Coef_mat = rep_matrix(0, Ni, Ni);
  vector[Ni] r_square = rep_vector(0, Ni);
  vector[Ni] res_sds;
  vector[Ni] res_var;
  vector[Nce] res_cor;
  vector[Nce] res_cov;
  matrix[Ni, Ni] Resid = rep_matrix(0.0, Ni, Ni);

  if (method != 100) rms_src = rms_src_p[1];

  if (method < 90) {
    int pos = 0;
    for (i in 2:Ni) {
      for (j in 1:(i - 1)) {
        if (cond_ind_mat[i, j] == 1) {
          pos += 1;
          Resid[i, j] = resids[pos] * rms_src_tmp;
          Resid[j, i] = Resid[i, j];
        }
      }
    }
    if (Ncond_ind == 0) rms_src = 0;
  }

  {
    vector[Nce_uniq] res_cor_u = res_cor_01 * 2 - 1;
    matrix[Ni, Ni] One_min_Beta_inv;
    matrix[Ni, Nce] loading_par_exp = rep_matrix(0, Ni, Nce);
    matrix[Ni, Ni] loading_par_exp_2;
    vector[Ni] delta_mat_ast;
    vector[Ni] total_var;
    matrix[Ni, Ni] cor_psi = identity_matrix(Ni);

    for (i in 1:Ni) {
      for (j in 1:Ni) {
        if (coef_pattern[i, j] != 0) {
          Coef_mat[i, j] = coefs[coef_pattern[i, j]];
        } else if (coef_fixed[i, j] > -990) {
          Coef_mat[i, j] = coef_fixed[i, j];
        }
      }
    }

    One_min_Beta_inv = inverse(identity_matrix(Ni) - Coef_mat);

    {
      int pos_err = 0;
      for (j in 1:(Ni - 1)) {
        for (i in (j + 1):Ni) {
          if (error_pattern[i, j] != 0) {
            pos_err += 1;
            cor_psi[i, j] = res_cor_u[error_pattern[i, j]];
            cor_psi[j, i] = cor_psi[i, j];
          }
        }
      }
    }

    // fix this, error correlation matrix is NaN
    res_var = find_factor_res_var(Coef_mat, cor_psi);
    res_sds = sqrt(res_var);

    {
      int pos_err = 0;
      for (j in 1:(Ni - 1)) {
        for (i in (j + 1):Ni) {
          if (error_pattern[i, j] != 0) {
            pos_err += 1;
            res_cor[pos_err] = res_cor_u[error_pattern[i, j]];
            res_cov[pos_err] = res_cor[pos_err] * res_sds[i] * res_sds[j];
            loading_par_exp[i, pos_err] = sqrt(
              abs(res_cor_u[error_pattern[i, j]]) * res_var[i]
            );
            loading_par_exp[j, pos_err] =
              sign(res_cor_u[error_pattern[i, j]]) * sqrt(
                abs(res_cor_u[error_pattern[i, j]]) * res_var[j]
              );
          }
        }
      }
    }

    loading_par_exp_2 = tcrossprod(loading_par_exp);
    delta_mat_ast = res_var - diagonal(loading_par_exp_2);

    r_mat = quad_form_sym(
      add_diag(loading_par_exp_2, delta_mat_ast),
      One_min_Beta_inv'
    );

    total_var = diagonal(r_mat);
    r_square = 1.0 - res_var ./ total_var;

    if (method < 90) {
      int pos = 0;
      for (i in 2:Ni) {
        for (j in 1:(i - 1)) {
          if (cond_ind_mat[i, j] == 1) {
             pos += 1;
             r_mat[i, j] += resids[pos] * rms_src_tmp * sqrt(total_var[i] * total_var[j]);
             r_mat[j, i] = r_mat[i, j];
          }
        }
      }
    }

    r_vec = matrix_log_vech(r_mat);

    int complete_pos = 0;
    int missing_pos = 0;
    int missing_pos_i = 0;

    for (i in 1:Ng) {
      real m_val;

      if (missing_var_indicator[i] == 0) {
        vector[Nisqd2] r_vec_sim;
        vector[Nisqd2] tmp_loc = r_vec;
        matrix[Nisqd2, Nisqd2] tmp_cov = r_obs_vec_cov[i];
        if (type >= 2) {
          m_val = exp(ln_v_int_wi[1] + X[i, ] * ln_v_beta_wi);
          tmp_cov = add_diag(r_obs_vec_cov[i], square(m_val));
          if (type == 3) {
            tmp_loc = r_vec + g_clus[, C_ID[i]];
          }
        }
        log_lik[i] = multi_normal_lpdf(r_obs_vec[i] | tmp_loc, tmp_cov);
        r_vec_sim = multi_normal_rng(tmp_loc, tmp_cov);
        D_obs += -2.0 * multi_normal_lpdf(r_obs_vec[i] | tmp_loc, tmp_cov);
        D_rep += -2.0 * multi_normal_lpdf(r_vec_sim | tmp_loc, tmp_cov);
      } else {
        int Ni_i = presence_count[i];
        int Nisqd2_i = (Ni_i * (Ni_i - 1)) %/% 2;
        vector[Nisqd2_i] r_vec_sim_i;
        m_val = 0;
        if (type >= 2) m_val = exp(ln_v_int_wi[1] + X[i, ] * ln_v_beta_wi);
        vector[Nisqd2_i] r_obs_vec_i;
        matrix[Nisqd2_i, Nisqd2_i] r_obs_vec_cov_i;
        matrix[Nisqd2_i, Nisqd2_i] L_vec_cov_i;
        vector[Nisqd2_i] tmp_loc;

        // create imputed matrix
        matrix[Ni_i, Ni_i] r_mat_i = r_mat_s[i][
          presence_indicator[1:Ni_i, i],
          presence_indicator[1:Ni_i, i]
        ];

        // get i-specific log-scale vector
        vector[Nisqd2_i] r_vec_i = matrix_log_vech(
          r_mat[
            presence_indicator[1:Ni_i, i],
            presence_indicator[1:Ni_i, i]
          ]
        );

        if (missing_element_indicator[i] == 1) {
          for (j in 1:Ni_i) {
            if (N_missing_items_array[j, i] == 1) {
              missing_pos_i += 1;
              r_mat_i[j, j] += var_shifts[missing_pos_i];
            }
          }

          for (j in 1:(Ni_i - 1)) {
            for (k in (j + 1):Ni_i) {
              if (r_mat_i[k, j] > 900) {
                missing_pos += 1;
                r_mat_i[k, j] = cor_missing[missing_pos] / sqrt(r_mat_i[j, j] * r_mat_i[k, k]);
                r_mat_i[j, k] = r_mat_i[k, j];
              }
            }
          }

          for (j in 1:Ni_i) r_mat_i[j, j] = 1.0;

          r_obs_vec_i = matrix_log_vech(r_mat_i);
          matrix[Nisqd2_i, Nisqd2_i] omega_r = omega_computer(r_mat_i) / (1.0 * Np[i]);
          matrix[Nisqd2_i, Nisqd2_i] jacob_mat = get_jacob(r_mat_i);
          r_obs_vec_cov_i = add_diag(gamma_avar(omega_r, jacob_mat), square(m_val));
        } else {
          r_obs_vec_i = r_obs_vec[i][1:Nisqd2_i];
          r_obs_vec_cov_i = add_diag(r_obs_vec_cov[i][1:Nisqd2_i, 1:Nisqd2_i], square(m_val));
        }

        L_vec_cov_i = cholesky_decompose(r_obs_vec_cov_i);
        // diag has been added already!
        if (type <= 2) {
          tmp_loc = r_vec_i;
        } else if (type == 3) {
          tmp_loc = r_vec_i + g_clus[1:Nisqd2_i, C_ID[i]];
        }
        log_lik[i] = multi_normal_lpdf(r_obs_vec_i | tmp_loc, r_obs_vec_cov_i);
        r_vec_sim_i = multi_normal_rng(tmp_loc, r_obs_vec_cov_i);
        D_obs += -2.0 * multi_normal_lpdf(r_obs_vec_i | tmp_loc, r_obs_vec_cov_i);
        D_rep += -2.0 * multi_normal_lpdf(r_vec_sim_i | tmp_loc, r_obs_vec_cov_i);
      }
    }
    ppp = D_rep > D_obs ? 1.0 : 0.0;
  }

  if (type == 2) {
    vector[Ng] ebx_wi = exp(ln_v_int_wi[1] + X * ln_v_beta_wi);

    rmsea_mn = mean(ebx_wi);
    v_mn = square(rmsea_mn);
    rmsea_beta_wi = ln_v_beta_wi * rmsea_mn;
  } else if (type == 3) {
    vector[Ng] ebx_wi = exp(ln_v_int_wi[1] + X * ln_v_beta_wi);

    rmsea_wi = mean(ebx_wi);
    rmsea_be = exp(ln_v_int_be[1]);
    v_wi = square(rmsea_wi);
    v_be = square(rmsea_be);
    v_mn = v_wi + v_be;
    prop_be = v_be / v_mn;
    rmsea_mn = sqrt(v_mn);
    rmsea_beta_wi = ln_v_beta_wi * rmsea_wi;
  }
}
