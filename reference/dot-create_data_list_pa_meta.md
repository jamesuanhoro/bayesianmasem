# Stan data helper function

A function that creates data list object passed to Stan

## Usage

``` r
.create_data_list_pa_meta(
  lavaan_object = NULL,
  method = "normal",
  type = "re",
  priors = NULL,
  cluster = NULL,
  partab = NULL,
  x_mat = NULL,
  conditional_re = TRUE,
  old_data = NULL
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

- partab:

  lavaan parameter table output with ceq.simple & std.lv = TRUE

- x_mat:

  (data.frame) Meta-analytic predictors. Each row contains the data for
  a group and each column is a variable. Column names should be
  labelled.

- conditional_re:

  (LOGICAL) Only relevant for analysis of correlation structures. If
  TRUE, sample levels of the study-level random effect (usually faster);
  If FALSE, don't (usually more efficient).

- old_data:

  original list of matrices, needed for PA

## Value

Data list object used in fitting Stan model
