
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
devtools::install_github("https://github.com/Dom-Owens-UoB/VAR_MOSUM", subdir = "mosumvar")
```

## Example


``` r
library(mosumvar)
data(voldata)

## Score test
ts <- mosumvar(voldata[,2:5], p=1, G=250, method = "Score")
ts

## Wald test
tw <- mosumvar(voldata[,2:5], p=1, G=250, method = "Wald")
tw

## Multiple Filtre
mf <- MFA( as.matrix(voldata[,2:5]), p=1, Gset=c(100,200,300), test = "Score" )
mf

## Subsample
ss <- mosum_sub( as.matrix(voldata[,2:5]), p=1, G=250, method = "Score")
ss

## MOSUM Binary Segmentation
bs <- MOSUMBS( as.matrix(voldata[,2:5]), p=1, G=250)
bs

## Dimension Reduction
dr <- mosum_univ( as.matrix(voldata[,2:8]), p=1, G=250, method = "Score", rm_cross_terms = T, global_resids = T)
dr

## Dimension Reduction with Multiplier Bootstrap
dr_mb <- mosum_univ( as.matrix(voldata[,2:8]), p=1, G=250, method = "Score", rm_cross_terms = T, global_resids = T, do_bootstrap = "multiplier")
dr_mb

## Simulate VAR data
n <- 1000
d <- 5
A <- diag(0.5, d) + rnorm(d^2, 0, 0.1)
simdata <-  VAR.sim(n,coeffs = A)
plot.ts(simdata)
```

