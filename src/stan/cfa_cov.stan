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
  real gen_matrix_beta_ii_lpdf(matrix S, matrix Omega, real n, real m) {
    int p = rows(S);
    real F_1 = eff(p, m) + eff(p, n) - eff(p, m + n);
    real F_2 = -((n - p - 1) * log_determinant_spd(S)) - (m * log_determinant_spd(Omega)) +
      ((m + n) * log_determinant_spd((m * Omega + n * S) / (m + n)));
    real ll = (F_1 + F_2) / -2.0;
    return(ll);
  }
}
data {
  int<lower = 0> Ng;  // number of groups
  array[Ng] int<lower = 0> Np;  // number persons by matrix
  int<lower = 0> Ni;  // number items
  array[Ng] matrix[Ni, Ni] S;  // covariance matrices
  int<lower = 0> Nf; // N_factors
  int<lower = 0> Nce; // N_correlated errors
  array[Ni, Ni] int error_pattern; // cor error matrix
  array[Ni, Nf] int loading_pattern;  // loading pattern
  array[Ni, Nf] real loading_fixed;  // loading pattern
  array[Ni] int res_var_pattern;  // res_var pattern
  array[Ni] real<lower = 0.0> res_var_fixed;  // res_var pattern
  array[Nf] int markers; // markers
  matrix[Nf, Nf] corr_mask;  // 1 for correlated factors, 0 otherwise
  real<lower = 1> shape_phi_c; // lkj prior shape for phi
  real<lower = 0> rm_par;  // rms scale parameter
  real<lower = 0> rs_par;  // residual sd parameter
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
}
transformed data {
  real sqrt_two = sqrt(2.0);
  real pi_sqrt_three = pi() / sqrt(3.0);
  int<lower = 0> Nl_uniq = 0;  // N_non-zero loadings unique
  int<lower = 0> Nl = 0;  // N_non-zero loadings
  int<lower = 0> Nl_fixed = 0;  // N_non-zero loadings unique
  int<lower = 0> Nrv_uniq = 0;  // N_non-zero res_var unique
  int<lower = 0> Nrv = 0;  // N_non-zero res_var
  int<lower = 0> Nrv_fixed = 0;  // N_non-zero res_var unique
  int Nisqd2 = (Ni * (Ni - 1)) %/% 2;
  int N_rms = 1;
  int N_alpha = 0;
  int N_complex = 0;
  int N_type_wi = 1;
  int N_type_be = 1;
  int<lower = 0> Nce_uniq = 0; // N_correlated errors unique

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

  for (i in 1:Ni) {
    if (res_var_pattern[i] != 0) {
      Nrv += 1;
      if (res_var_pattern[i] > Nrv_uniq) Nrv_uniq = res_var_pattern[i];
    } else if (res_var_fixed[i] < 990) {
      Nrv_fixed += 1;
    }
  }
}
parameters {
  vector<lower = 0.0, upper = 1.0>[N_rms] rms_src_p;
  vector<lower = 2.0>[N_alpha] gdp_alpha;
  vector[Nisqd2] resids;  // residual vector
  vector[Nl_uniq] loadings;  // loadings
  vector<lower = 0>[Nrv_uniq] res_sds_u;  // item residual sds heteroskedastic
  cholesky_factor_corr[Nf] phi_mat_chol;
  vector<lower = 0, upper = 1>[Nce_uniq] res_cor_01;  // correlated errors on 01
  vector[N_complex] loadings_complex;
  vector<lower = 0>[complex_struc] sigma_loadings_complex;
  vector<lower = 2.0>[complex_struc] gdp_loadings_complex;
  vector[N_type_wi] ln_v_int_wi;
  vector[p] ln_v_beta_wi;
  vector[N_type_be] ln_v_int_be;
  array[Nc] cov_matrix[Ni] E_clus;
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

  res_sds_u ~ student_t(3, 0, rs_par);
  phi_mat_chol ~ lkj_corr_cholesky(shape_phi_c);
  res_cor_01 ~ beta(rc_par, rc_par);

  ln_v_int_wi ~ normal(rm_i_l_par, rm_i_s_par);
  ln_v_int_be ~ normal(rm_i_l_par, rm_i_s_par);
  ln_v_beta_wi ~ normal(0, rm_b_s_par);

  {
    matrix[Ni, Ni] Omega;
    vector[Ni] total_var;

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

      for (i in 1:Ni) {
        if (res_var_pattern[i] != 0) {
          res_var[i] = square(res_sds_u[res_var_pattern[i]]);
        } else if (res_var_fixed[i] < 990) {
          res_var[i] = res_var_fixed[i];
        }
      }

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

      total_var = diagonal(Omega);

      if (method < 90) {
        int pos = 0;
        for (i in 2:Ni) {
          for (j in 1:(i - 1)) {
            pos += 1;
            Omega[i, j] += resids[pos] * rms_src_tmp * sqrt(total_var[i] * total_var[j]);
            Omega[j, i] = Omega[i, j];
          }
        }
      }
    }

    if (type == 3) {
      real m_val_c = 1.0 / square(exp(ln_v_int_be[1])) + Ni - 1;
      for (i in 1:Nc) {
        E_clus[i] ~ wishart(m_val_c, Omega / m_val_c);
      }
    }

    for (i in 1:Ng) {
      real m_val;

      if (type == 1) {
        target += wishart_lupdf(S[i] | Np[i] - 1.0, Omega / (Np[i] - 1.0));
      } else if (type >= 2) {
        m_val = 1.0 / square(exp(ln_v_int_wi[1] + X[i, ] * ln_v_beta_wi)) + Ni - 1;
        if (type == 2) {
          target += gen_matrix_beta_ii_lpdf(S[i] | Omega, Np[i] - 1.0, m_val);
        } else if (type == 3) {
          target += gen_matrix_beta_ii_lpdf(S[i] | E_clus[C_ID[i]], Np[i] - 1.0, m_val);
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
  vector[Nce] res_cov;
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
    vector[Ni] total_var;
    matrix[Ni, Ni] S_sim;
    matrix[Ni, Ni] Sigma_sim;

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

      for (i in 1:Ni) {
        if (res_var_pattern[i] != 0) {
          res_var[i] = square(res_sds_u[res_var_pattern[i]]);
        } else if (res_var_fixed[i] < 990) {
          res_var[i] = res_var_fixed[i];
        }
      }

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
      Omega = add_diag(lamb_phi_lamb + loading_par_exp_2, delta_mat_ast);

      total_var = diagonal(Omega);

      if (method < 90) {
        int pos = 0;
        for (i in 2:Ni) {
          for (j in 1:(i - 1)) {
            pos += 1;
            Omega[i, j] += resids[pos] * rms_src_tmp * sqrt(total_var[i] * total_var[j]);
            Omega[j, i] = Omega[i, j];
          }
        }
      }
    }

    for (i in 1:Ng) {
      real m_val;

      if (type == 1) {
        log_lik[i] = wishart_lpdf(S[i] | Np[i] - 1.0, Omega / (Np[i] - 1.0));
        S_sim = wishart_rng(Np[i] - 1.0, Omega / (Np[i] - 1.0));
        D_obs += -2.0 * wishart_lpdf(S[i] | Np[i] - 1.0, Omega / (Np[i] - 1.0));
        D_rep += -2.0 * wishart_lpdf(S_sim | Np[i] - 1.0, Omega / (Np[i] - 1.0));
      } else if (type >= 2) {
        m_val = 1.0 / square(exp(ln_v_int_wi[1] + X[i, ] * ln_v_beta_wi)) + Ni - 1;
        if (type == 2) {
          log_lik[i] = gen_matrix_beta_ii_lpdf(S[i] | Omega, Np[i] - 1.0, m_val);
          Sigma_sim = inv_wishart_rng(m_val, m_val * Omega);
          S_sim = wishart_rng(Np[i] - 1.0, Sigma_sim / (Np[i] - 1.0));
        } else if (type == 3) {
          log_lik[i] = gen_matrix_beta_ii_lpdf(
            S[i] | E_clus[C_ID[i]], Np[i] - 1.0, m_val
          );
          Sigma_sim = inv_wishart_rng(m_val, m_val * E_clus[C_ID[i]]);
          S_sim = wishart_rng(Np[i] - 1.0, Sigma_sim / (Np[i] - 1.0));
        }
        D_obs += -2.0 * gen_matrix_beta_ii_lpdf(S[i] | Omega, Np[i] - 1.0, m_val);
        D_rep += -2.0 * gen_matrix_beta_ii_lpdf(S_sim | Omega, Np[i] - 1.0, m_val);
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
