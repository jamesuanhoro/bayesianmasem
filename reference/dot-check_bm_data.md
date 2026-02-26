# mini-check for data inputs

mini-check for data inputs

## Usage

``` r
.check_bm_data(data, group, sample_cov, sample_nobs)
```

## Arguments

- data:

  An optional data frame containing the observed variables used in the
  model.

- group:

  An optional string identifying the grouping variable in the data
  object.

- sample_cov:

  (list of matrices) sample covariance or correlation matrices. The
  rownames and/or colnames must contain the observed variable names. For
  now, assumes there are no missing elements in the covariance matrices.

- sample_nobs:

  (vector of positive integer) Number of observations for each study.
