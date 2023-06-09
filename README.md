
# mosumvar

<!-- badges: start -->
[![Build Status](https://travis-ci.org/Dom-Owens-UoB/VAR_MOSUM.svg?branch=master)](https://travis-ci.org/Dom-Owens-UoB/VAR_MOSUM)
<!-- badges: end -->

A variety of methods for the analysis of multiple change points in time series data. 
All methods use Moving Sum (MOSUM) statistics for segmentation, and dependence properties are extracted with Vector Autoregression (VAR) models. 

## Installation

You can install the released version of mosumvar from [github](https://github.com/) with:

``` r
library(devtools)
devtools::install_github("https://github.com/Dom-Owens-UoB/mosumvar")
```

## Example


``` r
library(mosumvar)

## Simulate VAR data
n <- 1000
p <- 4
A <- matrix(-.1, p, p)
diag(A) <- 0.7
simdata <-  rbind(VAR.sim(n/2, coeffs = A), VAR.sim(n/2, coeffs = -A))
ts.plot(simdata)


## Score 
ts <- mosumvar(simdata, order=1,  method = "Score")
ts

## Wald 
tw <- mosumvar(simdata, order=1,  method = "Wald")
tw

## Multiscale
ms <- mosumvar.ms(simdata, order=1, method = "Score")
ms

## Subsample
ss <- mosumvar.sub(simdata, order=1, method = "Score")
ss

## MOSUM recursive segmentation
bs <- mosumvar.bs(simdata, order=1)
bs

## Dimension Reduction
dr <- mosumvar.uni(simdata, order=1, method = "Score", rm.cross.terms = T, global.resids = T, do.bootstrap = T)
dr
```

