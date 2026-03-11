# Another dataset from the metaSEM package.

The data set includes 131 correlation matrices on the mediating role of
leader–member exchange in transformational leadership (Boer et al. 2016)
. Data and documentation copied over from the metaSEM package (Cheung
2015) .

## Usage

``` r
boer16
```

## Format

### `boer16`

A list of data with the following structure:

- data:

  A list of correlation matrices. The variables are LMX (leader-member
  exchange), TFL (transformational leadership), JS (job satisfaction),
  OC (organizational commitment), and LE (leader effectiveness).

- n:

  A vector of sample sizes

- RelLMX:

  The reliability of LMX

- RelTFL:

  The reliability of TFL

## Source

<https://cran.r-project.org/web/packages/metaSEM/index.html>

## References

Boer D, Deinert A, Homan AC, Voelpel SC (2016). “Revisiting the
mediating role of leader-member exchange in transformational leadership:
the differential impact model.” *European Journal of Work and
Organizational Psychology*, **25**(6), 883–899. ISSN 1359-432X,
[doi:10.1080/1359432X.2016.1170007](https://doi.org/10.1080/1359432X.2016.1170007)
,
[2022-03-24](https://jamesuanhoro.github.io/bayesianmasem/reference/2022-03-24).  
  
Cheung MW (2015). “metaSEM: An R Package for Meta-Analysis using
Structural Equation Modeling.” *Frontiers in Psychology*, **5**(1521).
[doi:10.3389/fpsyg.2014.01521](https://doi.org/10.3389/fpsyg.2014.01521)
.

## Examples

``` r
if (FALSE) { # \dontrun{
model_syntax <- paste(
  "LMX ~ TFL",
  "JS + OC + LE ~ TFL + LMX",
  sep = "\n"
)
bmasem(model_syntax, sample_cov = boer16$data, sample_nobs = boer16$n)
} # }
```
