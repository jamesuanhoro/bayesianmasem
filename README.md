
# bayesianmasem

[![Project Status: Active The project has reached a stable, usable state
and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
![GitHub](https://img.shields.io/github/license/jamesuanhoro/bayesianmasem)
[![Codecov test
coverage](https://codecov.io/gh/jamesuanhoro/bayesianmasem/branch/main/graph/badge.svg)](https://app.codecov.io/gh/jamesuanhoro/bayesianmasem?branch=main)
[![minimal R
version](https://img.shields.io/badge/R%3E%3D-4.1.0-6666ff.svg)](https://cran.r-project.org/)
![GitHub R package
version](https://img.shields.io/github/r-package/v/jamesuanhoro/bayesianmasem)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/bayesianmasem)](https://cran.r-project.org/package=bayesianmasem)
![GitHub last
commit](https://img.shields.io/github/last-commit/jamesuanhoro/bayesianmasem)
[![R-CMD-check](https://github.com/jamesuanhoro/bayesianmasem/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/jamesuanhoro/bayesianmasem/actions/workflows/check-standard.yaml)
[![bayesianmasem status
badge](https://jamesuanhoro.r-universe.dev/badges/bayesianmasem)](https://jamesuanhoro.r-universe.dev)

#### Table of Contents

- [Package overview](#package-overview)
- [Features](#features)
- [Installation](#installation)
- [An example](#an-example)
- [Citations](#citations)

## Package overview

This is a package for Bayesian meta-analytic structural equation models
(MASEM, Cheung and Chan 2009), specifically meta-analysis of
confirmatory factor models. The package fits fixed-, random- and
dependent-samples MASEMs.

## Features

Features of `bayesianmasem` include:

- Uses lavaan-style syntax for model specification
- Fits fixed-effects, random-effects and **dependent-samples** MASEM
- Allows fixing loadings to known values or constraining loadings equal
- Estimates minor factor influences (Uanhoro 2023), i.e. full residual
  correlation matrix reflecting misspecification in the true shared
  correlation matrix
- Ability to relax simple-structure, i.e. allow for estimation of all
  cross-loadings using shrinkage-to-zero priors leading to more accurate
  estimates of interfactor correlations.
- Built atop [Stan](https://mc-stan.org/) (an efficient Bayesian
  sampler).

## Installation

`bayesianmasem` is hosted on GitHub, so we need the `remotes` package to
install it. We also need to install the `cmdstanr` package and CmdStan
in order to use Stan.

Instructions:

``` r
install.packages("remotes")  # install remotes

# next install cmdstanr and CmdStan:
install.packages(
  "cmdstanr",
  repos = c("https://mc-stan.org/r-packages/", getOption("repos"))
)
cmdstanr::check_cmdstan_toolchain(fix = TRUE)
cmdstanr::install_cmdstan()

# Then finally bayesianmasem:
remotes::install_github("jamesuanhoro/bayesianmasem")
```

## An example

We use the example in the Norton data (Norton et al. 2013). Let’s fit
the Zigmond-Snaith anxiety-depression model to the HADS assuming a
positive factor on select items.

``` r
model_syntax <- "
anxiety =~ x1 + x3 + x5 + x7 + x9 + x11 + x13
depression =~ x2 + x4 + x6 + x8 + x10 + x12 + x14
pos =~ a * x2 + a * x4 + a * x6 + a * x7 + a * x12 + a * x14
pos =~ 0 * x1 + 0 * x3 + 0 * x5 + 0 * x8 + 0 * x9 + 0 * x10 + 0 * x11 + 0 * x13
pos ~~ 0 * anxiety + 0 * depression
"
```

There are 14 items. The model assumes odd-numbered items load onto
anxiety, even-numbered item load on depression. Items 2, 4, 6, 7, 12,
and 14 load on a positive-wording factor, i.e. a method effect. We
assume the loadings on this positive factor to be equal, and this factor
assumed uncorrelated to the substantive factors. We also specify that
other items have a fixed loading of 0 on this positive factor – this
helps when we relax simple structure, we want to ensure the positive
method factor has 0 loading on other items.

We can fit a random-effects model to the input correlation matrices and
provide sample size information. The `Norton13` data come from the
[metaSEM package](https://cran.r-project.org/package=metaSEM) (Cheung
2015):

``` r
fit_ad <- bmasem(
  model_syntax,
  sample_cov = Norton13$data, sample_nobs = Norton13$n,
  method = "normal", simple_struc = FALSE
)
```

Setting `method = "normal"` implies the true correlation matrix
underlying the data differs from the structured correlation matrix due
to minor factors. Setting `simple_struc = FALSE` allows for all
cross-loadings to be estimated using priors that attempt to shrink the
cross-loadings to 0. Parameters for all functions including `bmasem()`
are documented: `?bmasem`.

## Citations

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-Cheung-metaSEM" class="csl-entry">

Cheung, Mike W.-L. 2015. “<span class="nocase">metaSEM</span>: An r
Package for Meta-Analysis Using Structural Equation Modeling.”
*Frontiers in Psychology* 5 (1521).
<https://doi.org/10.3389/fpsyg.2014.01521>.

</div>

<div id="ref-cheung_two-stage_2009" class="csl-entry">

Cheung, Mike W.-L., and Wai Chan. 2009. “A Two-Stage Approach to
Synthesizing Covariance Matrices in Meta-Analytic Structural Equation
Modeling.” *Structural Equation Modeling: A Multidisciplinary Journal*
16 (1): 28–53. <https://doi.org/10.1080/10705510802561295>.

</div>

<div id="ref-norton_hospital_2013" class="csl-entry">

Norton, Sam, Theodore Cosco, Frank Doyle, John Done, and Amanda Sacker.
2013. “The Hospital Anxiety and Depression Scale: A Meta Confirmatory
Factor Analysis.” *Journal of Psychosomatic Research* 74 (1): 74–81.
<https://doi.org/10.1016/j.jpsychores.2012.10.010>.

</div>

<div id="ref-uanhoro_modeling_2023" class="csl-entry">

Uanhoro, James Ohisei. 2023. “Modeling Misspecification as a Parameter
in Bayesian Structural Equation Models.” *Educational and Psychological
Measurement* 0 (0): 00131644231165306.
<https://doi.org/10.1177/00131644231165306>.

</div>

</div>
