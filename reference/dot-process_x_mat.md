# A function to process the meta-analytic SEM predictor matrix

A function to process the meta-analytic SEM predictor matrix

## Usage

``` r
.process_x_mat(x_mat = NULL, n_groups, type)
```

## Arguments

- x_mat:

  (data.frame) Meta-analytic predictors. Each row contains the data for
  a group and each column is a variable. Column names should be
  labelled.

- n_groups:

  Number of groups

- type:

  Type of meta-analysis

## Value

Returns numeric matrix of predictor variables
