# Pretty print model results

Nice printing of model results, optionally produces HTML document

## Usage

``` r
pp_summary(object, interval = 0.9, digits = 3, simple = TRUE, save_html = NULL)
```

## Arguments

- object:

  (bmasem) An object of
  [`bmasem-class`](https://jamesuanhoro.github.io/bayesianmasem/reference/bmasem-class.md)
  returned by
  [`bmasem`](https://jamesuanhoro.github.io/bayesianmasem/reference/bmasem.md).

- interval:

  (real in (0, 1)) Credible interval to return.

- digits:

  (positive integer) Number of decimal places to print in table

- simple:

  (Logical) TRUE to produce table with less information about
  parameters; FALSE: produces table with more information

- save_html:

  (string) Optional file name to save table as HTML

## Examples

``` r
if (FALSE) { # \dontrun{
model_syntax <- paste0(
  "distress =~ ", paste0("x", 1:14, collapse = " + "), "\n",
  "anxiety =~ ", paste0("x", seq(1, 14, 2), collapse = " + "), "\n",
  "depression =~ ", paste0("x", seq(2, 14, 2), collapse = " + ")
)
bmasem(
  model_syntax,
  sample_cov = Norton13$data, sample_nobs = Norton13$n, orthogonal = TRUE
)
pp_summary(fit)
pp_summary(fit, simple = FALSE)
} # }
```
