# Conditional indepepence estimates for path analysis

Obtain conditional independence estimates for path analysis. Estimates
are standardized.

## Usage

``` r
bmasem_ci_results(object, interval = 0.9, summarize = TRUE)
```

## Arguments

- object:

  (bmasem) An object of
  [`bmasem-class`](https://jamesuanhoro.github.io/bayesianmasem/reference/bmasem-class.md)
  returned by
  [`bayesianmasem`](https://jamesuanhoro.github.io/bayesianmasem/reference/bayesianmasem.md).

- interval:

  Confidence interval to select

- summarize:

  (LOGICAL) If TRUE (default): Return posterior summary as data.frame;
  If FALSE: Return posterior draws as data.frame.

## Value

Table of parameters related to conditional independence.

## Examples

``` r
if (FALSE) { # \dontrun{
model_syntax <- paste(
  "BI ~ ATT + SN + PBC",
  "BEH ~ PBC + BI",
  sep = "\n"
)
fit <- bmasem(
  model_syntax,
  sample_cov = cooke16$data, sample_nobs = cooke16$n
)
ci_table <- bmasem_ci_results(fit)
} # }
```
