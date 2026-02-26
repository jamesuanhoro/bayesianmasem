# Stan data helper function

A function that creates data list object passed to Stan

## Usage

``` r
.create_data_list_pooled(
  lavaan_object = NULL,
  method = "normal",
  simple_struc = TRUE,
  priors = NULL,
  partab = NULL,
  acov_mat = NULL,
  old_names = NULL
)
```

## Arguments

- lavaan_object:

  lavaan fit object of corresponding model

- method:

  (character) One of "normal", "lasso", "logistic", "GDP", or "none".
  See details below.

- simple_struc:

  (LOGICAL) Only relevant for CFAs. If TRUE (default): assume simple
  structure; If FALSE: estimate all cross-loadings using generalized

- priors:

  An object of
  [`bmasempriors-class`](https://jamesuanhoro.github.io/bayesianmasem/reference/bmasempriors-class.md).
  See
  [`new_bmasempriors`](https://jamesuanhoro.github.io/bayesianmasem/reference/new_bmasempriors.md)
  for more information.

- partab:

  lavaan parameter table output with ceq.simple & std.lv = TRUE

- acov_mat:

  (Optional) The asymptotic variance matrix of lower triangular half
  (column-order) of the correlation matrix to be used for correlation
  structure analysis.

- old_names:

  (Optional) Variable name order of original correlation matrix, used to
  reorder acov_mat.

## Value

Data list object used in fitting Stan model
