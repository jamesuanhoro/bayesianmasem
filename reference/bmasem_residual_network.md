# Compute residual network

Interpret the error correlations as a residual network model (Epskamp et
al. 2017) .

## Usage

``` r
bmasem_residual_network(object)
```

## Arguments

- object:

  (bmasem) An object of
  [`bmasem-class`](https://jamesuanhoro.github.io/bayesianmasem/reference/bmasem-class.md)
  returned by
  [`bmasem`](https://jamesuanhoro.github.io/bayesianmasem/reference/bmasem.md).

## Value

A data.frame containing posterior samples of the partial correlation
matrix.

## References

Epskamp S, Rhemtulla M, Borsboom D (2017). “Generalized Network
Psychometrics: Combining Network and Latent Variable Models.”
*Psychometrika*, **82**(4), 904–927. ISSN 1860-0980,
[doi:10.1007/s11336-017-9557-x](https://doi.org/10.1007/s11336-017-9557-x)
,
[2024-07-01](https://jamesuanhoro.github.io/bayesianmasem/reference/2024-07-01).

## Examples

``` r
if (FALSE) { # \dontrun{
model_syntax <- paste("distress =~", paste0("x", 1:14, collapse = " + "))
fit <- bmasem(
  model_syntax,
  sample_cov = Norton13$data, sample_nobs = Norton13$n
)
res_net <- bmasem_residual_network(fit)
p_corr_df <- posterior::summarise_draws(res_net)
n_items <- sqrt(nrow(p_corr_df))
p_corr_mat <- matrix(p_corr_df$mean, n_items)
p_corr_mat
qgraph::qgraph(p_corr_mat)
} # }
```
