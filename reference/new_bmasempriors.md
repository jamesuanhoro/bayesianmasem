# Set priors in package

Modify default priors in package.

## Usage

``` r
new_bmasempriors(
  lkj_shape = 2,
  ml_par = 0,
  sl_par = 0.5,
  rs_par = 1,
  rc_par = 2,
  mr_par = log(0.08),
  sr_par = 0.7,
  br_par = 0.5,
  rm_par = 0.15
)
```

## Arguments

- lkj_shape:

  (positive real) The shape parameter of the LKJ-prior on the
  interfactor correlation matrix in confirmatory factor models.

- ml_par:

  (real) The location parameter of the normal prior on loadings.

- sl_par:

  (positive real) The scale parameter of the normal prior on loadings.

- rs_par:

  (positive real) The scale parameter of the Student-t(3,0,) prior on
  residual standard deviations.

- rc_par:

  (positive real) The shape parameter of the Beta(rc_par, rc_par) prior
  on the residual error correlations.

- mr_par:

  (real) The location parameter of the normal prior on the log-RMSEA.

- sr_par:

  (positive real) The scale parameter of the normal prior on the
  log-RMSEA.

- br_par:

  (positive real) The scale parameter of the normal prior on the
  regression coefficients for the log-RMSEA.

- rm_par:

  (positive real) The scale parameter of the normal prior on the tau /
  CRMR parameter.

## Value

An object of
[`bmasempriors-class`](https://jamesuanhoro.github.io/bayesianmasem/reference/bmasempriors-class.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# Change LKJ shape parameter only
custom_priors <- new_bmasempriors(lkj_shape = 1.0)
model_syntax <- paste0(
  "distress =~ ", paste0("x", 1:14, collapse = " + "), "\n",
  "anxiety =~ ", paste0("x", seq(1, 14, 2), collapse = " + "), "\n",
  "depression =~ ", paste0("x", seq(2, 14, 2), collapse = " + ")
)
bmasem(
  model_syntax,
  sample_cov = Norton13$data, sample_nobs = Norton13$n, orthogonal = TRUE,
  priors = custom_priors
)
} # }
```
