# Fit Bayesian meta-analytic CFAs and path models with minor factors assumed.

A function to fit fixed-, random-effects, and clustered Bayesian
meta-analytic CFAs for covariance and correlation matrices and path
models for correlation matrices with minor factors assumed. The function
accomodate moderators. When fitting path models, correlation matrices
may be incomplete. Input matrices for CFA must be complete.

## Usage

``` r
bmasem(
  model = NULL,
  sample_cov = NULL,
  sample_nobs = NULL,
  correlation = TRUE,
  method = "normal",
  type = "re",
  orthogonal = FALSE,
  simple_struc = TRUE,
  x_mat = NULL,
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
  show_messages = TRUE,
  cluster = NULL,
  conditional_re = TRUE
)
```

## Arguments

- model:

  A description of the user-specified model, lavaan syntax.

- sample_cov:

  (list of matrices) sample covariance or correlation matrices. The
  rownames and/or colnames must contain the observed variable names.

- sample_nobs:

  (vector of positive integer) Number of observations for each study.

- correlation:

  (LOGICAL) If TRUE (default): analyze correlation matrices based on
  logarithm of a matrix transformation (Archakov and Hansen 2021) ; If
  FALSE: analyze covariance matrices methods in (Uanhoro 2023) .

- method:

  (character) One of "normal", "lasso", "logistic", "GDP", or "none".
  See details below.

- type:

  (character) One of "fe", "re", or "dep" for fixed-effects,
  random-effects, and dependent-samples MASEM respectively. The "dep"
  argument is experimental, see details below.

- orthogonal:

  (LOGICAL) If TRUE: constrain all factors orthogonal (overrides model
  syntax); If FALSE (default): according to model syntax.

- simple_struc:

  (LOGICAL) Only relevant for CFAs. If TRUE (default): assume simple
  structure; If FALSE: estimate all cross-loadings using generalized

- x_mat:

  (data.frame) Meta-analytic predictors. Each row contains the data for
  a group and each column is a variable. Column names should be
  labelled.

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

- cluster:

  An optional integer vector identifying the cluster each group belongs
  to. Asssume there are five groups, the first three belong to cluster 1
  and the last two belong to cluster 2, then the argument would be:
  `cluster = c(1, 1, 1, 2, 2)`. This feature is experimental, see
  details below.

- conditional_re:

  (LOGICAL) Only relevant for analysis of correlation structures. If
  TRUE, sample levels of the study-level random effect (usually faster);
  If FALSE, don't (usually more efficient).

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

When `type = "dep"`, the user must supply the cluster IDs, see cluster
parameter documentation above. However, this feature is experimental.
Additionally, the cluster inputs are not validated.

## References

Archakov I, Hansen PR (2021). “A New Parametrization of Correlation
Matrices.” *Econometrica*, **89**(4), 1699–1715. ISSN 1468-0262,
[doi:10.3982/ECTA16910](https://doi.org/10.3982/ECTA16910) ,
<https://onlinelibrary.wiley.com/doi/abs/10.3982/ECTA16910>.  
  
Armagan A, Dunson DB, Lee J (2013). “Generalized double Pareto
shrinkage.” *Statistica Sinica*, **23**(1), 119–143. ISSN 1017-0405,
<https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3903426/>.  
  
Uanhoro JO (2023). “Hierarchical Covariance Estimation Approach to
Meta-Analytic Structural Equation Modeling.” *Structural Equation
Modeling: A Multidisciplinary Journal*, **30**(4), 532–546. ISSN
1070-5511,
[doi:10.1080/10705511.2022.2142128](https://doi.org/10.1080/10705511.2022.2142128)
,
[2023-06-19](https://jamesuanhoro.github.io/bayesianmasem/reference/2023-06-19).

## Examples

``` r
if (FALSE) { # \dontrun{
model_syntax <- "# latent variable definitions
F1 =~ JP1 + JP2 + JP3
F2 =~ JN1 + JN2 + JN4 + JN4
F3 =~ TD1 + TD2"
bmasem(
  model_syntax,
  sample_cov = issp89$data, sample_nobs = issp89$n
)
model_syntax <- paste0(
  "distress =~ ", paste0("x", 1:14, collapse = " + "), "\n",
  "anxiety =~ ", paste0("x", seq(1, 14, 2), collapse = " + "), "\n",
  "depression =~ ", paste0("x", seq(2, 14, 2), collapse = " + ")
)
bmasem(
  model_syntax,
  sample_cov = Norton13$data, sample_nobs = Norton13$n, orthogonal = TRUE
)
} # }
```
