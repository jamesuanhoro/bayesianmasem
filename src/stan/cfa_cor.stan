// Each S is dist according to some E
// And each E follows some structured matrix
// And that structured matrix is IW (see Wu & Browne)
functions {
  int sign(real x) {
    if (x > 0)
      return 1;
    else
      return -1;
  }
  real generalized_double_pareto_lpdf(vector x, real alpha, real scale) {
    // generalized double Pareto
    // https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3903426/
    return(sum(
      -log(2) - log(scale) - (alpha + 1.0) * log1p(abs(x) / (scale * alpha))
    ));
  }
  real eff(int p, real x) {
    return(
      2 * lmgamma(p, x / 2) - x * p * log(x / 2) + x * p
    );
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
  array[Ng] vector[(Ni * (Ni - 1)) %/% 2] r_obs_vec;  // covariance matrices
  array[Ng] matrix[(Ni * (Ni - 1)) %/% 2, (Ni * (Ni - 1)) %/% 2] r_obs_vec_cov;  // covariance matrices
  int<lower = 0> Nf; // N_factors
  int<lower = 0> Nce; // N_correlated errors
  array[Ni, Ni] int error_pattern; // cor error matrix
  array[Ni, Nf] int loading_pattern;  // loading pattern
  array[Ni, Nf] real loading_fixed;  // loading pattern
  array[Nf] int markers; // markers
  matrix[Nf, Nf] corr_mask;  // 1 for correlated factors, 0 otherwise
  real<lower = 1> shape_phi_c; // lkj prior shape for phi
  real<lower = 0> rm_par;  // rms scale parameter
  real<lower = 1> rc_par;  // residual corr parameter
  int<lower = 1, upper = 100> method; // which method
  int<lower = 0, upper = 1> complex_struc;
  int p;  // number of moderators
  matrix[Ng, p] X;  // moderator matrix
  real rm_i_l_par;  // meta-reg int loc hyper-parameter -- log(.08)
  real<lower = 0> rm_i_s_par;  // meta-reg int scale hyper-parameter -- .7
  real<lower = 0> rm_b_s_par;  // meta-reg beta scale hyper-parameter -- .5
  int<lower = 0> Nc;  // number of clusters
  array[Ng] int<lower = 0, upper = Nc> C_ID;  // cluster ID
  int<lower = 1, upper = 3> type; // which type
  matrix[Ni, Nf] load_est;
  matrix[Ni, Nf] load_se;
  int<lower = 0, upper = 1> marginal_re;
}
transformed data {
  real sqrt_two = sqrt(2.0);
  real pi_sqrt_three = pi() / sqrt(3.0);
  int<lower = 0> Nl_uniq = 0;  // N_non-zero loadings unique
  int<lower = 0> Nl = 0;  // N_non-zero loadings
  int<lower = 0> Nl_fixed = 0;  // N_non-zero loadings unique
  int Nisqd2 = (Ni * (Ni - 1)) %/% 2;
  int Nisqd2_vec = (Ni * (Ni - 1)) %/% 2;
  int N_rms = 1;
  int N_alpha = 0;
  int N_complex = 0;
  int N_type_wi = 1;
  int N_type_be = 1;
  array[Ng] matrix[Nisqd2_vec, Nisqd2_vec] L_vec_cov;
  int<lower = 0> Nce_uniq = 0; // N_correlated errors unique

  for (i in 1:Ng) {
    L_vec_cov[i] = cholesky_decompose(r_obs_vec_cov[i]);
  }

  if (method >= 90) {
    Nisqd2 = 0;
  }

  if (method == 100) {
    N_rms = 0;
  }

  if (method == 4) N_alpha = 1;

  if (type < 2) N_type_wi = 0;
  if (type < 3) N_type_be = 0;

  for (i in 1:Ni) {
    for (j in 1:Nf) {
      if (loading_pattern[i, j] != 0) {
        Nl += 1;
        if (loading_pattern[i, j] > Nl_uniq) Nl_uniq = loading_pattern[i, j];
      } else if (loading_fixed[i, j] > -990) {
        Nl_fixed += 1;
      }
    }
  }

  if (complex_struc == 1) {
    N_complex = Ni * Nf - Nl - Nl_fixed;
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
  vector<lower = 0.0, upper = 1.0>[N_rms] rms_src_p;
  vector<lower = 2.0>[N_alpha] gdp_alpha;
  vector[Nisqd2] resids;  // residual vector
  vector<lower = -1.0, upper = 1.0>[Nl_uniq] loadings;  // loadings
  cholesky_factor_corr[Nf] phi_mat_chol;
  vector<lower = 0, upper = 1>[Nce_uniq] res_cor_01;  // correlated errors on 01
  vector<lower = -1.0, upper = 1.0>[N_complex] loadings_complex;
  vector<lower = 0, upper = 1.0>[complex_struc] sigma_loadings_complex;
  vector<lower = 2.0>[complex_struc] gdp_loadings_complex;
  vector[N_type_wi] ln_v_int_wi;
  vector[p] ln_v_beta_wi;
  vector[N_type_be] ln_v_int_be;
  array[Ng * N_type_wi * marginal_re] vector[Nisqd2_vec] c_clus;
  matrix[Nisqd2_vec, Nc] g_clus;
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

  if (complex_struc == 1) {
    sigma_loadings_complex ~ normal(0, 0.25);
    gdp_loadings_complex ~ lognormal(1, 1);
    target += generalized_double_pareto_lpdf(
      loadings_complex | gdp_loadings_complex[1], sigma_loadings_complex[1]);
  }

  phi_mat_chol ~ lkj_corr_cholesky(shape_phi_c);
  res_cor_01 ~ beta(rc_par, rc_par);

  ln_v_int_wi ~ normal(rm_i_l_par, rm_i_s_par);
  ln_v_int_be ~ normal(rm_i_l_par, rm_i_s_par);
  ln_v_beta_wi ~ normal(0, rm_b_s_par);

  {
    matrix[Ni, Ni] Omega;
    vector[Nisqd2_vec] r_vec;

    {
      vector[Ni] res_var;
      matrix[Nf, Nf] phi_mat = multiply_lower_tri_self_transpose(phi_mat_chol) .* corr_mask;
      vector[Nce_uniq] res_cor_u = res_cor_01 * 2 - 1;
      matrix[Ni, Nf] Load_mat = rep_matrix(0, Ni, Nf);
      matrix[Ni, Ni] lamb_phi_lamb;
      matrix[Ni, Nce] loading_par_exp = rep_matrix(0, Ni, Nce);
      matrix[Ni, Ni] loading_par_exp_2;
      vector[Ni] delta_mat_ast;

      {
        int pos_complex = 0;
        for (i in 1:Ni) {
          for (j in 1:Nf) {
            if (loading_pattern[i, j] != 0) {
              loadings[loading_pattern[i, j]] ~ normal(load_est[i, j], load_se[i, j]);
              Load_mat[i, j] = loadings[loading_pattern[i, j]];
            } else if (loading_fixed[i, j] > -990) {
              Load_mat[i, j] = loading_fixed[i, j];
            } else if (complex_struc == 1) {
              pos_complex += 1;
              Load_mat[i, j] = loadings_complex[pos_complex];
            }
          }
        }
      }

      lamb_phi_lamb = quad_form_sym(phi_mat, Load_mat');

      res_var = 1.0 - diagonal(lamb_phi_lamb);

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
      Omega = add_diag(lamb_phi_lamb + loading_par_exp_2, delta_mat_ast);

      if (method < 90) {
        int pos = 0;
        for (i in 2:Ni) {
          for (j in 1:(i - 1)) {
            pos += 1;
            Omega[i, j] += resids[pos] * rms_src_tmp;
            Omega[j, i] = Omega[i, j];
          }
        }
      }
    }

    r_vec = matrix_log_vech(Omega);

    if (type == 3) {
      to_vector(g_clus) ~ normal(0, exp(ln_v_int_be[1]));
    }

    for (i in 1:Ng) {
      real m_val;

      if (type == 1) {
        target += multi_normal_cholesky_lupdf(
          r_obs_vec[i] | r_vec, L_vec_cov[i]
        );
      } else if (type >= 2) {
        m_val = exp(ln_v_int_wi[1] + X[i, ] * ln_v_beta_wi);
        if (marginal_re == 1) {
          c_clus[i] ~ normal(r_vec, m_val);
          if (type == 2) {
            target += multi_normal_cholesky_lupdf(
              r_obs_vec[i] | c_clus[i], L_vec_cov[i]
            );
          } else if (type == 3) {
            target += multi_normal_cholesky_lupdf(
              r_obs_vec[i] | c_clus[i] + g_clus[, C_ID[i]], L_vec_cov[i]
            );
          }
        } else (marginal_re == 0) {
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
    }
  }
}
generated quantities {
  real D_rep = 0.0;
  real D_obs = 0.0;
  real<lower = 0, upper = 1> ppp;
  real<lower = 0> rms_src = 0.0;  // RMSE of residuals
  matrix[Ni, Nf] Load_mat = rep_matrix(0, Ni, Nf);
  matrix[Nf, Nf] phi_mat = multiply_lower_tri_self_transpose(phi_mat_chol) .* corr_mask;
  vector[Ni] res_sds;
  vector[Ni] res_var;
  vector[Nce] res_cor;
  matrix[Ni, Ni] Resid = rep_matrix(0.0, Ni, Ni);
  real v_mn = 0.0;
  real rmsea_mn = sqrt(v_mn);
  real v_wi = 0.0;
  real rmsea_wi = sqrt(v_wi);
  real v_be = 0.0;
  real rmsea_be = sqrt(v_be);
  vector[p] rmsea_beta_wi;
  real prop_be = 0.0;
  vector[Ng] log_lik;
  matrix[Ni, Ni] Omega;

  if (method != 100) rms_src = rms_src_p[1];

  {
    vector[Nisqd2_vec] r_vec;
    vector[Nisqd2_vec] r_vec_sim;

    {
      vector[Nce_uniq] res_cor_u = res_cor_01 * 2 - 1;
      matrix[Ni, Ni] lamb_phi_lamb;
      matrix[Ni, Nce] loading_par_exp = rep_matrix(0, Ni, Nce);
      matrix[Ni, Ni] loading_par_exp_2;
      vector[Ni] delta_mat_ast;

      {
        int pos_complex = 0;
        for (i in 1:Ni) {
          for (j in 1:Nf) {
            if (loading_pattern[i, j] != 0) {
              Load_mat[i, j] = loadings[loading_pattern[i, j]];
            } else if (loading_fixed[i, j] > -990) {
              Load_mat[i, j] = loading_fixed[i, j];
            } else if (complex_struc == 1) {
              pos_complex += 1;
              Load_mat[i, j] = loadings_complex[pos_complex];
            }
          }
        }
      }

      lamb_phi_lamb = quad_form_sym(phi_mat, Load_mat');

      res_var = 1.0 - diagonal(lamb_phi_lamb);

      res_sds = sqrt(res_var);

      {
        int pos_err = 0;
        for (j in 1:(Ni - 1)) {
          for (i in (j + 1):Ni) {
            if (error_pattern[i, j] != 0) {
              pos_err += 1;
              res_cor[pos_err] = res_cor_u[error_pattern[i, j]];
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
      Omega = add_diag(lamb_phi_lamb + loading_par_exp_2, delta_mat_ast);

      if (method < 90) {
        int pos = 0;
        for (i in 2:Ni) {
          for (j in 1:(i - 1)) {
            pos += 1;
            Omega[i, j] += resids[pos] * rms_src_tmp;
            Omega[j, i] = Omega[i, j];
          }
        }
      }
    }

    r_vec = matrix_log_vech(Omega);

    for (i in 1:Ng) {
      real m_val;
      vector[Nisqd2_vec] tmp_loc = r_vec;
      matrix[Nisqd2_vec, Nisqd2_vec] tmp_cov = r_obs_vec_cov[i];

      if (type >= 2) {
        m_val = exp(ln_v_int_wi[1] + X[i, ] * ln_v_beta_wi);
        tmp_cov = add_diag(r_obs_vec_cov[i], square(m_val));
        if (type == 3) {
          tmp_cov = add_diag(tmp_cov, square(exp(ln_v_int_be[1])));
        }
      }
      log_lik[i] = multi_normal_lpdf(r_obs_vec[i] | tmp_loc, tmp_cov);
      r_vec_sim = multi_normal_rng(tmp_loc, tmp_cov);
      D_obs += -2.0 * multi_normal_lpdf(r_obs_vec[i] | tmp_loc, tmp_cov);
      D_rep += -2.0 * multi_normal_lpdf(r_vec_sim | tmp_loc, tmp_cov);
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

  if (method < 90) {
    int pos = 0;

    for (i in 2:Ni) {
      for (j in 1:(i - 1)) {
        pos += 1;
        Resid[i, j] = resids[pos] * rms_src_tmp;
        Resid[j, i] = Resid[i, j];
      }
    }
  }

  for (j in 1:Nf) {
    if (Load_mat[markers[j], j] < 0) {
      Load_mat[, j] *= -1.0;
      phi_mat[, j] *= -1.0;
      phi_mat[j, ] *= -1.0;
    }
  }
}
