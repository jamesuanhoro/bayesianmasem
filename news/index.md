# Changelog

## bayesianmasem 0.2.5

- Added composite reliability calculation for each factor in
  [`bmasem_composite_reliability()`](https://jamesuanhoro.github.io/bayesianmasem/reference/bmasem_composite_reliability.md).
- Added residual network (<https://doi.org/10.1007/s11336-017-9557-x>)
  in
  [`bmasem_residual_network()`](https://jamesuanhoro.github.io/bayesianmasem/reference/bmasem_residual_network.md).
  This simply transforms the random error matrix into a partial
  correlation matrix.
- Whitening vectors to use univariate instead of multi-normal
  likelihoods for correlation matrix inputs resulting in speedup for
  larger matrices.
- Can handle missing elements in correlation matrices when pooling with
  [`bmasem_stage_1()`](https://jamesuanhoro.github.io/bayesianmasem/reference/bmasem_stage_1.md)
- Some code cleanup

## bayesianmasem 0.2.4

- Fixed default prior for location on log-RMSEA

## bayesianmasem 0.2.3

- The stage 1 analysis in
  [`bmasem_stage_1()`](https://jamesuanhoro.github.io/bayesianmasem/reference/bmasem_stage_1.md)
  now additionally returns a pooled partial correlation matrix for
  network analysis.

## bayesianmasem 0.2.2

- Added option to marginalize the random-effects in `type = "RE"` models
  (longer runtime, larger ESS and ultimately more efficient)

## bayesianmasem 0.2.1

- If fixed-effects, then set no moderators
- Removed equality test on real numbers in Stan
- Some instantiate-forced changes

## bayesianmasem 0.2.0

- Added meta-analytic SEM predictor matrix
- Export model-implied matrix

## bayesianmasem 0.1.3

- Ensure asymptotic variance of log-correlation matrix is symmetric

## bayesianmasem 0.1.2

- Fixed bug in
  [`bmasem_stage_2()`](https://jamesuanhoro.github.io/bayesianmasem/reference/bmasem_stage_2.md)
  function. The `acov_mat` is now correctly ordered based on lavaan
  object

## bayesianmasem 0.1.1

- Pooled object can also be analyzed using bayesianmasem package with
  [`bmasem_stage_2()`](https://jamesuanhoro.github.io/bayesianmasem/reference/bmasem_stage_2.md)
  function

## bayesianmasem 0.1.0

- Now allows for pooling correlation matrices
  [`bmasem_stage_1()`](https://jamesuanhoro.github.io/bayesianmasem/reference/bmasem_stage_1.md)
  function, pooled matrix can be analyzed in a different software

## bayesianmasem 0.0.2

- Added covariance matrix analysis methods
- Added parameter equality constraints for residual correlations

## bayesianmasem 0.0.1

- Package works and has basic tests
- Added a `NEWS.md` file to track changes to the package.
