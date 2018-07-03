# MiSPU [![License](http://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](http://www.gnu.org/licenses/gpl-2.0.html) [![CRAN](http://www.r-pkg.org/badges/version/MiSPU)](http://cran.rstudio.com/package=MiSPU) 

This package implements a new test called adaptive interaction sum of powered score (aiSPU) test for testing high-dimensional parameters under generalized linear models (GLMs) with high-dimensional nuisance parameters. Some related methods have been implemented as well. *The package will be available on CRAN later.* 

## Installation
To install the stable version from CRAN, simply run the following from an R console (*Not applicable now*):

```r
install.packages("aispu")
```

To install the latest development builds directly from GitHub, run this instead:

```r
if (!require("devtools"))
  install.packages("devtools")
devtools::install_github("ChongWu-Biostat/aispu")
```

## aispu

In this package, we implement the adaptive interaction sum of powered score (aiSPU) test based on the truncated Lasso penalty and an adaptive testing idea, which can maintain high statistical power and correct Type I error rates across a wide range of alternatives. There are two ways to calculate its p-values: a bootstrap-based way and an asymptotics-based way.


```r
library(aispu)

# Generate the data (codes for the simulations in the manuscript)n = 30signal.r = 0nInformative = 3p = 40seed = 1s = 0.01non.zero = floor((p/2) * s)alpha = c(rep(0,p/2 - non.zero), runif(non.zero,-signal.r,signal.r))beta = c(rep(2,nInformative), rep(0,(p/2- 3)), alpha)

dat = sim_data(seed, n = n, p = p, beta = beta)X = dat$XY = dat$Ycov = NULLX.tmp = Xcov2 = X.tmp[,1:(p/2)]X = X.tmp[,(p/2 + 1):p]# run the test, bootstrap-based methodaispu(Y, X,cov = NULL, cov2, pow = c(1:6, Inf), model= "gaussian",penalty = "tlp", n.perm = 1000,resample = "boot")

# run the test, asymptotics-based methodaispu(Y, X,cov = NULL, cov2, pow = c(1:6, Inf), model= "gaussian",penalty = "tlp", n.perm = 1000,resample = "asy")

```


