# Another dataset from the metaSEM package.

The data set includes 34 correlation matrices on the how well the theory
of planned behaviour predicts alcohol consumption (Cooke et al. 2016) .
We adjusted the third correlation matrix to the nearest valid
correlation matrix. Data and documentation copied over from the metaSEM
package (Cheung 2015) .

## Usage

``` r
cooke16
```

## Format

### `cooke16`

A list of data with the following structure:

- data:

  A list of correlation matrices. The variables are SN (subjective
  norm), ATT (attitude), PBC (perceived behavior control), BI
  (behavioral intention), and BEH (behavior).

- n:

  A vector of sample sizes

- MeanAge:

  Mean age of the participants except for Ajzen and Sheikh (2013), which
  is the median age, and Glassman, et al. (2010a) to Glassman, et al.
  (2010d), which are based on the range of 18 to 24.

- Female:

  Percentage of female participants.

## Source

<https://cran.r-project.org/web/packages/metaSEM/index.html>

## References

Cheung MW (2015). “metaSEM: An R Package for Meta-Analysis using
Structural Equation Modeling.” *Frontiers in Psychology*, **5**(1521).
[doi:10.3389/fpsyg.2014.01521](https://doi.org/10.3389/fpsyg.2014.01521)
.  
  
Cooke R, Dahdah M, Norman P, French DP (2016). “How well does the theory
of planned behaviour predict alcohol consumption? A systematic review
and meta-analysis.” *Health Psychology Review*, **10**(2), 148–167. ISSN
1743-7199, 1743-7202,
[doi:10.1080/17437199.2014.947547](https://doi.org/10.1080/17437199.2014.947547)
.

## Examples

``` r
if (FALSE) { # \dontrun{
model_syntax <- paste(
  "BI ~ ATT + SN + PBC",
  "BEH ~ PBC + BI",
  sep = "\n"
)
bmasem(model_syntax, sample_cov = cooke16$data, sample_nobs = cooke16$n)
} # }
```
