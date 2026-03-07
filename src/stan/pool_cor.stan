// Each S is dist according to some E
// And each E follows some structured matrix
// And that structured matrix is IW (see Wu & Browne)
functions {
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
  real rm_i_l_par;  // meta-reg int loc hyper-parameter -- log(.08)
  real<lower = 0> rm_i_s_par;  // meta-reg int scale hyper-parameter -- .7
  int<lower = 0> Nc;  // number of clusters
  array[Ng] int<lower = 0, upper = Nc> C_ID;  // cluster ID
  int<lower = 1, upper = 3> type; // which type
  int<lower = 0, upper = 1> conditional_re;
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
      matrix[Nisqd2, Nisqd2] omega_r = omega_computer(r_mat_s[i]) / (1.0 * Np[i]);
      matrix[Nisqd2, Nisqd2] jacob_mat = get_jacob(r_mat_s[i]);
      r_obs_vec[i] = matrix_log_vech(r_mat_s[i]);
      r_obs_vec_cov[i] = quad_form(omega_r, jacob_mat');
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
    }
    N_missing_items += sum(N_missing_items_array[, i]);
  }

  int N_complete_studies = Ng - sum(missing_var_indicator);

  if (type < 2) N_type_wi = 0;
  if (type < 3) N_type_be = 0;
}
parameters {
  vector<lower = 0>[N_missing_items] var_shifts;
  cholesky_factor_corr[Ni] r_chol;
  vector[N_type_wi] ln_v_int_wi;
  vector[N_type_be] ln_v_int_be;
  array[N_complete_studies * N_type_wi * conditional_re] vector[Nisqd2] c_clus;
  matrix[Nisqd2, Nc] g_clus;
  vector<lower = -1.0, upper = 1.0>[N_missing_elements] cor_missing;
}
model {
  var_shifts ~ std_normal();
  r_chol ~ lkj_corr_cholesky(1.0);

  ln_v_int_wi ~ normal(rm_i_l_par, rm_i_s_par);
  ln_v_int_be ~ normal(rm_i_l_par, rm_i_s_par);

  cor_missing ~ normal(0, .5);

  {
    matrix[Ni, Ni] r_mat = multiply_lower_tri_self_transpose(r_chol);
    vector[Nisqd2] r_vec = matrix_log_vech(r_mat);
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
          m_val = exp(ln_v_int_wi[1]);
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

        // create imputed matrix
        matrix[Ni_i, Ni_i] r_mat_i = r_mat_s[i][
          presence_indicator[1:Ni_i, i],
          presence_indicator[1:Ni_i, i]
        ];

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
        }

        // get i-specific log-scale vector
        vector[Nisqd2_i] r_vec_i = matrix_log_vech(
          r_mat[
            presence_indicator[1:Ni_i, i],
            presence_indicator[1:Ni_i, i]
          ]
        );

        // create required vector, matrices
        m_val = 0;
        if (type >= 2) m_val = exp(ln_v_int_wi[1]);
        vector[Nisqd2_i] r_obs_vec_i = matrix_log_vech(r_mat_i);
        matrix[Nisqd2_i, Nisqd2_i] omega_r = omega_computer(r_mat_i) / (1.0 * Np[i]);
        matrix[Nisqd2_i, Nisqd2_i] jacob_mat = get_jacob(r_mat_i);
        matrix[Nisqd2_i, Nisqd2_i] r_obs_vec_cov_i = add_diag(
          quad_form(omega_r, jacob_mat'), square(m_val)
        );
        matrix[Nisqd2_i, Nisqd2_i] L_vec_cov_i = cholesky_decompose(r_obs_vec_cov_i);

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
  matrix[Ni, Ni] r_mat = multiply_lower_tri_self_transpose(r_chol);
  vector[Nisqd2] r_vec = matrix_log_vech(r_mat);
  matrix[Ni, Ni] p_mat = -1 * cov2cor(chol2inv(r_chol));
  real v_mn = 0.0;
  real rmsea_mn = sqrt(v_mn);
  real v_wi = 0.0;
  real rmsea_wi = sqrt(v_wi);
  real v_be = 0.0;
  real rmsea_be = sqrt(v_be);
  real prop_be = 0.0;
  vector[Ng] log_lik;

  for (i in 1:Ni) p_mat[i, i] = 1.0;

  {
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
          m_val = exp(ln_v_int_wi[1]);
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
        vector[Nisqd2_i] tmp_loc;

        // create imputed matrix
        matrix[Ni_i, Ni_i] r_mat_i = r_mat_s[i][
          presence_indicator[1:Ni_i, i],
          presence_indicator[1:Ni_i, i]
        ];

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
        }

        // get i-specific log-scale vector
        vector[Nisqd2_i] r_vec_i = matrix_log_vech(
          r_mat[
            presence_indicator[1:Ni_i, i],
            presence_indicator[1:Ni_i, i]
          ]
        );

        // create required vector, matrices
        m_val = 0;
        if (type >= 2) m_val = exp(ln_v_int_wi[1]);
        vector[Nisqd2_i] r_obs_vec_i = matrix_log_vech(r_mat_i);
        matrix[Nisqd2_i, Nisqd2_i] omega_r = omega_computer(r_mat_i) / (1.0 * Np[i]);
        matrix[Nisqd2_i, Nisqd2_i] jacob_mat = get_jacob(r_mat_i);
        matrix[Nisqd2_i, Nisqd2_i] r_obs_vec_cov_i = add_diag(
          quad_form(omega_r, jacob_mat'), square(m_val)
        );

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
    real ebx_wi = exp(ln_v_int_wi[1]);

    rmsea_mn = ebx_wi;
    v_mn = square(rmsea_mn);
  } else if (type == 3) {
    real ebx_wi = exp(ln_v_int_wi[1]);

    rmsea_wi = ebx_wi;
    rmsea_be = exp(ln_v_int_be[1]);
    v_wi = square(rmsea_wi);
    v_be = square(rmsea_be);
    v_mn = v_wi + v_be;
    prop_be = v_be / v_mn;
    rmsea_mn = sqrt(v_mn);
  }
}
