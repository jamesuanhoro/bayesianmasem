# bayesianmasem 0.2.2

* Added option to marginalize the random-effects in `type = "RE"` models (longer runtime, larger ESS and ultimately more efficient)

# bayesianmasem 0.2.1

* If fixed-effects, then set no moderators
* Removed equality test on real numbers in Stan
* Some instantiate-forced changes

# bayesianmasem 0.2.0

* Added meta-analytic SEM predictor matrix
* Export model-implied matrix

# bayesianmasem 0.1.3

* Ensure asymptotic variance of log-correlation matrix is symmetric

# bayesianmasem 0.1.2

* Fixed bug in `bmasem_stage_2()` function. The `acov_mat` is now correctly ordered based on lavaan object

# bayesianmasem 0.1.1

* Pooled object can also be analyzed using bayesianmasem package with `bmasem_stage_2()` function

# bayesianmasem 0.1.0

* Now allows for pooling correlation matrices `bmasem_stage_1()` function, pooled matrix can be analyzed in a different software

# bayesianmasem 0.0.2

* Added covariance matrix analysis methods
* Added parameter equality constraints for residual correlations

# bayesianmasem 0.0.1

* Package works and has basic tests
* Added a `NEWS.md` file to track changes to the package.
