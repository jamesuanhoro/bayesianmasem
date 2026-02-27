# Obtain composite reliability from a fitted model

Get composite reliability

## Usage

``` r
bmasem_composite_reliability(object, interval = 0.9, return_draws = FALSE)
```

## Arguments

- object:

  (bmasem) An object of
  [`bmasem-class`](https://jamesuanhoro.github.io/bayesianmasem/reference/bmasem-class.md)
  returned by
  [`bmasem`](https://jamesuanhoro.github.io/bayesianmasem/reference/bmasem.md).

- interval:

  (real in (0, 1)) Credible interval to return.

- return_draws:

  (LOGICAL) If `TRUE`, returns the realised distribution of reliability.
  If `FALSE` (default), return the distribution summary.

## Value

- If `return_draws = FALSE`, returns the summary of posterior
  distribution of composite reliability for each factor.

- If `return_draws = TRUE`, returns draws from the realised distribution
  of composite reliability for each factor.

## Details

This metric is the ratio of common variance (CV) to total variance (TV)
for each factor.

- **Common Variance**: For factor \\f\\, the common variance is: \$\$
  \mathrm{CV}\_f = (\sum\_{i \in I_f} \lambda\_{if})^2 \phi\_{ff}, \$\$
  where \\I_f\\ indexes the indicators of factor \\f\\.

- **Total Variance**: For factor \\f\\, the total variance is: \$\$
  \mathrm{TV}\_f = \sum\_{i \in I_f} \sum\_{j \in I_f} \omega\_{ij},
  \$\$ where \\\Omega = = (\omega\_{ij})\\ is the model-implied
  covariance matrix of the indicators.

## Examples

``` r
if (FALSE) { # \dontrun{
model_syntax <- paste0(
  "distress =~ ", paste0("x", 1:14, collapse = " + "), "\n",
  "anxiety =~ ", paste0("x", seq(1, 14, 2), collapse = " + "), "\n",
  "depression =~ ", paste0("x", seq(2, 14, 2), collapse = " + ")
)
fit <- bmasem(
  model_syntax,
  sample_cov = Norton13$data, sample_nobs = Norton13$n, orthogonal = TRUE
)
bmasem_composite_reliability(fit)
} # }
```
