# Clean up Stan fit helper function

A function that cleans up model returned by Stan

## Usage

``` r
.clean_up_stan_fit(stan_fit, data_list, interval = 0.9, priors)
```

## Arguments

- stan_fit:

  Stan fit

- data_list:

  Data list object passed to Stan

- interval:

  (real in (0, 1)) Credible interval to return.

- priors:

  An object of
  [`bmasempriors-class`](https://jamesuanhoro.github.io/bayesianmasem/reference/bmasempriors-class.md).
  See
  [`new_bmasempriors`](https://jamesuanhoro.github.io/bayesianmasem/reference/new_bmasempriors.md)
  for more information.

## Value

An object of
[`bmasem-class`](https://jamesuanhoro.github.io/bayesianmasem/reference/bmasem-class.md)
