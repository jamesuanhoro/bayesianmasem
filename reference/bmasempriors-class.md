# A class for setting up priors.

A class for setting up priors.

## Slots

- `lkj_shape`:

  (positive real) The shape parameter of the LKJ-prior on the
  interfactor correlation matrix in confirmatory factor models.

- `ml_par`:

  (real) The location parameter of the normal prior on loadings.

- `sl_par`:

  (positive real) The scale parameter of the normal prior on loadings.

- `rs_par`:

  (positive real) The scale parameter of the Student-t(3,0,) prior on
  residual standard deviations.

- `rc_par`:

  (positive real) The shape parameter of the Beta(rc_par, rc_par) prior
  on the residual error correlations.

- `mr_par`:

  (real) The location parameter of the normal prior on the log-RMSEA.

- `sr_par`:

  (positive real) The scale parameter of the normal prior on the
  log-RMSEA.

- `br_par`:

  (positive real) The scale parameter of the normal prior on the
  regression coefficients for the log-RMSEA.

- `rm_par`:

  (positive real) The scale parameter of the normal prior on the tau /
  CRMR parameter.
