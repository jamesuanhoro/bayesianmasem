# Transform unbounded vector to correlation matrix

Transform unbounded vector to correlation matrix

## Usage

``` r
g_inv_map(gamma_in, tol_val = 1e-13, ret_vec = FALSE)
```

## Arguments

- gamma_in:

  (vector) Strict lower half vector of log(correlation matrix)

- tol_val:

  (positive real) Tolerance for calculation

- ret_vec:

  (Logical) If TRUE, return strict lower half vector, ELSE: return
  correlation matrix.

## Value

Transformed correlation matrix
