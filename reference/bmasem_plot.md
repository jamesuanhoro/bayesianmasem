# Plot model results

Plot model parameters

## Usage

``` r
bmasem_plot(object, type = "trace", subset = NULL, ...)
```

## Arguments

- object:

  (bmasem) An object of
  [`bmasem-class`](https://jamesuanhoro.github.io/bayesianmasem/reference/bmasem-class.md)
  returned by
  [`bmasem`](https://jamesuanhoro.github.io/bayesianmasem/reference/bmasem.md).

- type:

  (character) Type of plot: "trace" for traceplots and "hist" for
  histograms.

- subset:

  (character) Subset of parameters: NULL (Default) showing all estimated
  parameters; Any other response will be used as regular expressions to
  subset the parameters. It can be loading names or types of parameters.

- ...:

  additional arguments to relevant bayesplot function

## Value

bayesplot object

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
bmasem_plot(fit)
} # }
```
