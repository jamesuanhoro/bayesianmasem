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
  real rm_i_l_par;  // meta-reg int loc hyper-parameter -- log(.08)
  real<lower = 0> rm_i_s_par;  // meta-reg int scale hyper-parameter -- .7
  int<lower = 0> Nc;  // number of clusters
  array[Ng] int<lower = 0, upper = Nc> C_ID;  // cluster ID
  int<lower = 1, upper = 3> type; // which type
}
transformed data {
  int Nisqd2 = (Ni * (Ni - 1)) %/% 2;
  int N_type_wi = 1;
  int N_type_be = 1;
  array[Ng] matrix[Nisqd2, Nisqd2] L_vec_cov;

  for (i in 1:Ng) {
    L_vec_cov[i] = cholesky_decompose(r_obs_vec_cov[i]);
  }

  if (type < 2) N_type_wi = 0;
  if (type < 3) N_type_be = 0;
}
parameters {
  cholesky_factor_corr[Ni] r_chol;
  vector[N_type_wi] ln_v_int_wi;
  vector[N_type_be] ln_v_int_be;
  array[Ng * N_type_wi] vector[Nisqd2] c_clus;
  matrix[Nisqd2, Nc] g_clus;
}
model {
  r_chol ~ lkj_corr_cholesky(1.0);

  ln_v_int_wi ~ normal(rm_i_l_par, rm_i_s_par);
  ln_v_int_be ~ normal(rm_i_l_par, rm_i_s_par);

  {
    vector[Nisqd2] r_vec = matrix_log_vech(multiply_lower_tri_self_transpose(r_chol));

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
        m_val = exp(ln_v_int_wi[1]);
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
      }
    }
  }
}
generated quantities {
  real D_rep = 0.0;
  real D_obs = 0.0;
  real<lower = 0, upper = 1> ppp;
  vector[Nisqd2] r_vec = matrix_log_vech(multiply_lower_tri_self_transpose(r_chol));
  matrix[Ni, Ni] r_mat = multiply_lower_tri_self_transpose(r_chol);
  real v_mn = 0.0;
  real rmsea_mn = sqrt(v_mn);
  real v_wi = 0.0;
  real rmsea_wi = sqrt(v_wi);
  real v_be = 0.0;
  real rmsea_be = sqrt(v_be);
  real prop_be = 0.0;
  vector[Ng] log_lik;

  {
    vector[Nisqd2] r_vec_sim;

    for (i in 1:Ng) {
      real m_val;
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
