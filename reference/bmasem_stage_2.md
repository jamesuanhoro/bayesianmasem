# Fit Bayesian meta-analytic CFAs and path models with minor factors assumed.

A function that takes a pooled correlation matrix object and returns a
fitted CFA or path model.

## Usage

``` r
bmasem_stage_2(
  model = NULL,
  pool_fit = NULL,
  method = "normal",
  orthogonal = FALSE,
  simple_struc = TRUE,
  seed = 12345,
  warmup = 1000,
  sampling = 1000,
  refresh = (warmup + sampling)/10,
  adapt_delta = 0.9,
  max_treedepth = 10,
  chains = 3,
  ncores = max(parallel::detectCores() - 2, 1),
  priors = new_bmasempriors(),
  show = TRUE,
  show_messages = TRUE
)
```

## Arguments

- model:

  A description of the user-specified model, lavaan syntax.

- pool_fit:

  The list returned by the `bmasem_stage_1` function.

- method:

  (character) One of "normal", "lasso", "logistic", "GDP", or "none".
  See details below.

- orthogonal:

  (LOGICAL) If TRUE: constrain all factors orthogonal (overrides model
  syntax); If FALSE (default): according to model syntax.

- simple_struc:

  (LOGICAL) Only relevant for CFAs. If TRUE (default): assume simple
  structure; If FALSE: estimate all cross-loadings using generalized

- seed:

  (positive integer) seed, set to obtain replicable results.

- warmup:

  (positive integer) The number of warmup iterations to run per chain.

- sampling:

  (positive integer) The number of post-warmup iterations to run per
  chain, retained for inference.

- refresh:

  (positive integer) How often to print the status of the sampler.

- adapt_delta:

  (real in (0, 1)) Increase to resolve divergent transitions.

- max_treedepth:

  (positive integer) Increase to resolve problems with maximum tree
  depth.

- chains:

  (positive integer) The number of Markov chains to run.

- ncores:

  (positive integer) The number of chains to run in parallel.

- priors:

  An object of
  [`bmasempriors-class`](https://jamesuanhoro.github.io/bayesianmasem/reference/bmasempriors-class.md).
  See
  [`new_bmasempriors`](https://jamesuanhoro.github.io/bayesianmasem/reference/new_bmasempriors.md)
  for more information.

- show:

  (Logical) If TRUE, show table of results, if FALSE, do not show table
  of results. As an example, use FALSE for simulation studies.

- show_messages:

  (Logical) If TRUE, show messages from Stan sampler, if FALSE, hide
  messages.

## Value

An object of
[`bmasem-class`](https://jamesuanhoro.github.io/bayesianmasem/reference/bmasem-class.md)

## Details

CFAs assume standardized factors. Latent variable regression models are
not yet implemented.

There are different methods for estimating models in this package:

- `normal`: under belief that minor factor influences are on average
  zero with continuous deviations away from zero.

- `lasso`: under belief that minor factor influences are largely zero
  with a small number of non-zero residual covariances.

- `logistic`: for similar belief as normal but more readily accomodates
  extreme outliers.

- `GDP`: to mimic a global-local approach, i.e. attempt to shrink near 0
  residual covariances to 0 with minimal shrinking for larger residual
  covariances (Armagan et al. 2013) .

- `none`: if intending to ignore the influence of minor factors.

## References

Armagan A, Dunson DB, Lee J (2013). “Generalized double Pareto
shrinkage.” *Statistica Sinica*, **23**(1), 119–143. ISSN 1017-0405,
<https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3903426/>.

## Examples

``` r
if (FALSE) { # \dontrun{
pool_fit_issp <- bmasem_stage_1(
  sample_cov = issp89$data, sample_nobs = issp89$n
)
model_syntax_issp <- "# latent variable definitions
F1 =~ JP1 + JP2 + JP3
F2 =~ JN1 + JN2 + JN4 + JN4
F3 =~ TD1 + TD2"
bmasem_stage_2(model_syntax_issp, pool_fit_issp)
pool_fit_norton <- bmasem_stage_1(
  sample_cov = Norton13$data, sample_nobs = Norton13$n
)
model_syntax_norton <- paste0(
  "distress =~ ", paste0("x", 1:14, collapse = " + "), "\n",
  "anxiety =~ ", paste0("x", seq(1, 14, 2), collapse = " + "), "\n",
  "depression =~ ", paste0("x", seq(2, 14, 2), collapse = " + ")
)
bmasem_stage_2(model_syntax_norton, pool_fit_norton, orthogonal = TRUE)
} # }
```
