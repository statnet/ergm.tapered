# `ergm.tapered`: Tapered Exponential-Family Models for Networks

<img src="man/figures/ergm.tapered_hl.png" align="right" width="250" height="250" alt="RDS network"/>

[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/ergm.tapered?color=2ED968)](https://cranlogs.r-pkg.org/)
[![cran version](https://www.r-pkg.org/badges/version/ergm.tapered)](https://cran.r-project.org/package=ergm.tapered)
[![Coverage status](https://codecov.io/gh/statnet/ergm.tapered/branch/master/graph/badge.svg)](https://codecov.io/github/statnet/ergm.tapered?branch=master)
[![R build status](https://github.com/statnet/ergm.tapered/workflows/R-CMD-check/badge.svg)](https://github.com/statnet/ergm.tapered/actions)

A set of terms and functions implementing Tapered exponential-family random
graph models (ERGMs).  Tapered ERGMs are a modification of ERGMs that reduce the
effects of phase transitions, and with properly chosen hyper-parameters,
provably removes all multiphase behavior.

Each ERGM has a corresponding Tapered ERGM. Indeed, the `ergm.tapered` package fits any `ergm` as it is based on `ergm` itself.

# Installation

<!-- The package is available on CRAN and can be installed using -->

<!--```{r} -->
<!--install.packages("ergm.tapered") -->
<!--``` -->

To install the latest development version from github, you can also use:

```{r}
# If devtools is not installed:
# install.packages("devtools")

devtools::install_github("statnet/ergm.tapered")
```
<!-- devtools::install_github("statnet/ergm", rev="tapered") -->
For now, this will install a variant of the current version of `ergm` (the `tapered` branch) that is needed. In a bit the CRAN version of `ergm` will have it.

# Implementation

Load package and example data

```{r}
library(ergm.tapered)
data(sampson)
```

Sampson (1969) recorded the social interactions among a group of monks while he was a resident as an experimenter at the cloister.
     Of particular interest are the data on positive affect relations
     ("liking," in which each monk was asked if he had positive relations
     to each of the other monks. Each monk ranked only his top three
     choices (or four, in the case of ties) on "liking".  Here, we
     consider a directed edge from monk A to monk B to exist if A
     nominated B among these top choices at any one of three time points during the year.
     For details see:
     
```
help(sampson)
```

We can make a quick visualization of the network

```{r}
plot(sampson)
```

![](.github/figures/samplikeplot.png)<!-- -->

A natural model is one that includes a term measuring the transitivity
of triples in the network, defined as a set of edges {(i,j), (j,k),
(i,k)}.

``` r
    fit <- ergm(samplike ~ edges + ttriple)
```

    Starting maximum pseudolikelihood estimation (MPLE):

    ...

    Optimizing with step length 0.0054.
    The log-likelihood improved by 27.9793.
    Estimating equations are not within tolerance region.
    Iteration 3 of at most 60:
    Error in ergm.MCMLE(init, nw, model, initialfit = (initialfit <- NULL),  : 
      Unconstrained MCMC sampling did not mix at all. Optimization cannot continue.
    Calls: ergm -> ergm.MCMLE
    In addition: Warning message:
    In ergm_MCMC_sample(s, control, theta = mcmc.init, verbose = max(verbose -  :
      Unable to reach target effective size in iterations alotted.
    Calls: ergm -> ergm.MCMLE -> ergm_MCMC_sample
   

This fit fails to converge computationally as the model is near
degenerate. We could try to get it to fit by working on the
computational algorithm. However, `ergm.tapered` considers a variant of
the ERGM that reflects our prior belief that the true generating process
is non-degenerate:

``` r
fit <- ergm.tapered(samplike ~ edges + ttriple)
```

 Starting maximum pseudolikelihood estimation (MPLE):

    ...

    Iteration 2 of at most 60:
    Optimizing with step length 1.0000.
    The log-likelihood improved by 0.0047.
    Precision adequate twice. Stopping.
    Finished MCMLE.
    Evaluating log-likelihood at the estimate. Fitting the dyad-independent submodel...
    Bridging between the dyad-independent submodel and the full model...
    Setting up bridge sampling...
    Using 16 bridges: 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 .
    Bridging finished.
    This model was fit using MCMC.  To examine model diagnostics and check
    for degeneracy, use the mcmc.diagnostics() function.

``` r
summary(fit)
```

     Results:

            Estimate Std. Error MCMC % z value Pr(>|z|)    
    edges   -1.83992    0.27737      0  -6.634  < 1e-04 ***
    ttriple  0.21810    0.05816      0   3.750 0.000177 ***
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    The estimated tapering scaling factor is 2.

         Null Deviance: 424.2  on 306  degrees of freedom
     Residual Deviance: 358.0  on 304  degrees of freedom
     
    AIC: 362  BIC: 369.4  (Smaller is better. MC Std. Err. = 0)

This tapered ERGM fits. The summary indicates that the coefficient on
the transitive triple term is positive (about 0.22) and statistically
above zero. This coefficient has the same interpretation as those in a standard ERGM.

This model fixes the tapering parameter at 2 units. Let’s try to taper
less by increasing the tapering parameter to 3 (`r=3`):

``` r
fit <- ergm.tapered(samplike ~ edges + ttriple, r=3)
```

    Starting maximum pseudolikelihood estimation (MPLE):

    ...

    Iteration 2 of at most 60:
    This model was fit using MCMC.  To examine model diagnostics and check
    for degeneracy, use the mcmc.diagnostics() function.

``` r
summary(fit)
```

     Results:

            Estimate Std. Error MCMC % z value Pr(>|z|)    
    edges   -1.82890    0.27246      0  -6.712   <1e-04 ***
    ttriple  0.21174    0.05255      0   4.029   <1e-04 ***
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    The estimated tapering scaling factor is 3.

         Null Deviance: 424.2  on 306  degrees of freedom
     Residual Deviance: 358.6  on 304  degrees of freedom
     
    AIC: 362.6  BIC: 370.1  (Smaller is better. MC Std. Err. = 0)

It does not effect it much.

The software allows the tapering to be estimated based on the shape of
the distributions of the model statistics. Let’s try that:

``` r
fit <- ergm.tapered(samplike ~ edges + ttriple, fixed=FALSE)
```

    Starting maximum pseudolikelihood estimation (MPLE):

    Iteration 4 of at most 60:
    Optimizing with step length 1.0000.
    The log-likelihood improved by 0.0003.
    Precision adequate twice. Stopping.
    Finished MCMLE.

``` r
summary(fit)
```

     Results:

            Estimate Std. Error MCMC % z value Pr(>|z|)    
    edges   -1.83372    0.28426      0  -6.451   <1e-04 ***
    ttriple  0.21236    0.09938      0   2.137   0.0326 *  
    ---
    Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    The estimated tapering scaling factor is 2.813.

         Null Deviance: 424.2  on 306  degrees of freedom
     Residual Deviance: 359.0  on 304  degrees of freedom
     
    AIC: 365  BIC: 376.2  (Smaller is better. MC Std. Err. = 0)

The estimated tapering parameter is about `2.8` (It is printed under the
coefficient table). This is between the default value and the second
guess. The coefficients of the ERGM terms are about the same as before.

Enjoy trying `ergm.tapering`!

<!-- A more detailed vignette with information on measurement error and diagnostics can be found here: [[link to katie's page]] -->

See the following papers for more information and examples:

#### Statistical Methodology

* Fellows, Ian E. and Handcock, Mark S. (2017) [Removing Phase Transitions from Gibbs Measures](https://proceedings.mlr.press/v54/fellows17a/fellows17a.pdf), *Proceedings of the 20th International Conference on Artificial Intelligence and Statistics (AISTATS)s*, Volume 54.
* Blackburn, Bart and Handcock, Mark S. (2022) [Practical Network Modeling via Tapered Exponential-family Random Graph Models](https://doi.org/10.1080/10618600.2022.2116444), *Journal of Computational and Graphical Statistics*.

## Public and Private repositories

To facilitate open development of the package while giving the core developers an opportunity to publish on their developments before opening them up for general use, this project comprises two repositories:
* A public repository `statnet/ergm.tapered`
* A private repository `statnet/ergm.tapered-private`

The intention is that all developments in `statnet/ergm.tapered-private` will eventually make their way into `statnet/ergm.tapered` and onto CRAN.

Developers and Contributing Users to the Statnet Project should read https://statnet.github.io/private/ for information about the relationship between the public and the private repository and the workflows involved.

## Latest Windows and MacOS binaries

A set of binaries is built after every commit to the repository. We strongly encourage testing against them before filing a bug report, as they may contain fixes that have not yet been sent to CRAN. They can be downloaded through the following links:

* [MacOS binary (a `.tgz` file in a `.zip` file)](https://nightly.link/statnet/ergm.tapered/workflows/R-CMD-check.yaml/master/macOS-rrelease-binaries.zip)
* [Windows binary (a `.zip` file in a `.zip` file)](https://nightly.link/statnet/ergm.tapered/workflows/R-CMD-check.yaml/master/Windows-rrelease-binaries.zip)

You will need to extract the MacOS `.tgz` or the Windows `.zip` file from the outer `.zip` file before installing. These binaries are usually built under the latest version of R and their operating system and may not work under other versions.

You may also want to install the corresponding latest binaries for packages on which `ergm.tapered` depends, in particular [`statnet.common`](https://github.com/statnet/statnet.common), [`network`](https://github.com/statnet/network), and [`ergm`](https://github.com/statnet/ergm).
