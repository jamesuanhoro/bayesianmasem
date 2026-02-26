# Stan data helper function

A function that creates data list object passed to Stan

## Usage

``` r
.create_data_list_meta(
  lavaan_object = NULL,
  method = "normal",
  type = "re",
  simple_struc = TRUE,
  priors = NULL,
  cluster = NULL,
  correlation = TRUE,
  partab = NULL,
  x_mat = NULL,
  conditional_re = TRUE
)
```

## Arguments

- lavaan_object:

  lavaan fit object of corresponding model

- method:

  (character) One of "normal", "lasso", "logistic", "GDP", or "none".
  See details below.

- type:

  (character) One of "fe", "re", or "dep" for fixed-effects,
  random-effects, and dependent-samples MASEM respectively. The "dep"
  argument is experimental, see details below.

- simple_struc:

  (LOGICAL) Only relevant for CFAs. If TRUE (default): assume simple
  structure; If FALSE: estimate all cross-loadings using generalized

- priors:

  An object of
  [`bmasempriors-class`](https://jamesuanhoro.github.io/bayesianmasem/reference/bmasempriors-class.md).
  See
  [`new_bmasempriors`](https://jamesuanhoro.github.io/bayesianmasem/reference/new_bmasempriors.md)
  for more information.

- cluster:

  An optional integer vector identifying the cluster each group belongs
  to. Asssume there are five groups, the first three belong to cluster 1
  and the last two belong to cluster 2, then the argument would be:
  `cluster = c(1, 1, 1, 2, 2)`. This feature is experimental, see
  details below.

- correlation:

  (LOGICAL) If TRUE (default): analyze correlation matrices based on
  logarithm of a matrix transformation (Archakov and Hansen 2021) ; If
  FALSE: analyze covariance matrices methods in (Uanhoro 2023) .

- partab:

  lavaan parameter table output with ceq.simple & std.lv = TRUE

- x_mat:

  (data.frame) Meta-analytic predictors. Each row contains the data for
  a group and each column is a variable. Column names should be
  labelled.

- conditional_re:

  (LOGICAL) Only relevant for analysis of correlation structures. If
  TRUE, sample levels of the study-level random effect; If FALSE, don't.

## Value

Data list object used in fitting Stan model
