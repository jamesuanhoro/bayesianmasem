# Fit random-effects Bayesian meta-analytic CFAs with minor factors assumed.

A function to pool correlation matrices permitting fixed-,
random-effects, and clustered-samples pooling. Correlation matrices must
be complete. This will change in the near future.

## Usage

``` r
bmasem_stage_1(
  sample_cov = NULL,
  sample_nobs = NULL,
  type = "re",
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

- sample_cov:

  (list of matrices) sample covariance or correlation matrices. The
  rownames and/or colnames must contain the observed variable names. For
  now, assumes there are no missing elements in the covariance matrices.

- sample_nobs:

  (vector of positive integer) Number of observations for each study.

- type:

  (character) One of "fe", "re", or "dep" for fixed-effects,
  random-effects, and dependent-samples MASEM respectively. The "dep"
  argument is experimental, see details below.

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
  If FALSE, don't.

## Value

A list containing fit indices, pooled correlation matrix and its
asymptotic covariance matrix, pooled partial correlation matrix and its
asymptotic covariance matrix, Stan object and data_list used to fit Stan
object.

## Details

When `type = "dep"`, the user must supply the cluster IDs, see cluster
parameter documentation above. However, this feature is experimental.
Additionally, the cluster inputs are not validated.

## References

There are no references for Rd macro `\insertAllCites` on this help
page.

## Examples
