
#' @title Default bandwidth
#' @keywords internal
default.G <- function(n, p, order){
  # dim <- order*(order*p + 1 ) + p*(p+1)/2
  # dimScore <- p*(p*order + 1 )/2 + order*(order+1)/2
  # fl <- floor(4/3 * n^(2/3))
  # if(method == "Wald") out <- min(fl, dim)
  # if(method == "Score") out <- min(fl, dimScore)
  c1 <-   -0.77720; c2 <- 0.19567; c3 <- 1.70305
  err <- 0.598#1/log(log(n))
  out <- exp( -(c2*log(log(sqrt(order*p^2)) ) + c3*log(log(n)) - log(err) ) / c1 )
  out <- max(round(out), 10)
  return(out)
}

#' @title Dimension warning
#' @keywords internal
dim_warning <- function(n, G, p, order, method) {
  dim <- p*(p*order + 1 ) + p*(p+1)/2
  dimScore <- p*(p*order + 1 )/2 + p*(p+1)/2
  fl <- floor(n^(2/3))

  W3dim <- paste0("Bandwidth too small relative to model dimensions: set G > p(p*order + 1) * log(p(p*order + 1)) = ", 
                  p*(p*order + 1) * log(p*(p*order + 1)), "\n")
  Wfl <- paste0("Bandwidth small relative to sample size: consider setting G > floor(n^(2/3)) = ", fl, "\n" )
  Wlarge <- "Large dimensions: consider `option = univariate`\n"

  if(G < dim & method == "Wald") warning(paste0("Not enough degrees of freedom for Wald method: set G > p(p*order + 1) + p(p+1)/2 = ", 
                                                p*(p*order + 1)* log(p*(p*order + 1)), "\n"))
  if(G < dimScore & method == "Score")warning(paste0("Not enough degrees of freedom for Score method: set G > p(p*order + 1)/2 + p(p+1)/2 = ", dimScore, "\n"))
  if(G < p*(p*order + 1)* log(p*(p*order + 1)) ) warning(W3dim)
  if(G < fl ) warning(Wfl)
  if(p*(p*order + 1 ) > 30) warning(Wlarge)
  
} 

#' Segment data under a VAR model
#'
#' @param x data matrix
#' @param order integer VAR model order
#' @param G integer MOSUM bandwidth; see reference for default
#' @param method detector, one of \code{"Wald", "Score"}
#' @param estim estimator method, one of \code{"C", "H"}
#' @param var.estim variance estimator method, one of \code{"Local", "Global"}
#' @param alpha Numeric significance level
#' @param criterion location procedure, one of \code{"eps", "eta"}
#' @param nu Numeric location procedure hyperparameter
#' @param do.bootstrap Boolean, determine threshold via bootstrap method
#' @param n.bootstrap Integer; number of bootstrap replicates
#' @param thresh rejection threshold; see reference for default
#' @param do.plot Boolean, return plot
#' @return List containing
#' \itemize{
#'   \item{\code{thresh}}{ input}
#'   \item{\code{mosum}}{ vector of mosum detector values}
#'   \item{\code{cps}}{ estimated change points}
#'   \item{\code{plot}}{ detector plot}
#'   \item{\code{estim}}{ input}
#' }
#' @examples
#' data(voldata)
#' mosumvar(voldata[,2:5], 1)
mosumvar <- function(x, order, G = NULL, method = c("Wald","Score"), estim = c("C","H"), var.estim = c("Local","Global"),  
                     alpha = 0.05, criterion= c("eps","eta"), nu= 0.25, do.bootstrap = FALSE, n.bootstrap = 1000, thresh = NULL, do.plot = TRUE){
  x <- as.matrix(x)
  order <- as.integer(order)
  method <- match.arg(method, c("Wald","Score"))
  estim <- match.arg(estim, c("C","H"))
  var.estim <- match.arg(var.estim, c("Local","Global"))
  alpha <- min(1, max(0, alpha))
  criterion <- match.arg(criterion, c("eps","eta"))
  nu <- min(1, max(0, nu))
  if(is.null(G)) G <- default.G(nrow(x), ncol(x), order)
  if(order + ncol(x) == 1){
    out <- mosumvar.uni(x, order, G, alpha = alpha, estim= estim, method = method,
                        criterion = criterion, nu=nu, do.bootstrap = do.bootstrap, n.bootstrap = n.bootstrap, thresh = thresh, do.plot = do.plot)
  } else{
    if(method== "Wald") out <- wald(x, order, G, alpha = alpha, estim= estim, 
                                    criterion = criterion, nu=nu, do.bootstrap = do.bootstrap, n.bootstrap = n.bootstrap, thresh = thresh, do.plot = do.plot)
    if(method== "Score") out <- score(x, order, G, alpha = alpha, estim= estim, var.estim = var.estim,
                                      criterion = criterion, nu=nu, do.bootstrap = do.bootstrap, n.bootstrap = n.bootstrap, thresh = thresh, do.plot = do.plot)
  }
  attr(out, "class") <- "mosumvar"
  return(out)
}


#' Simulate multiple time series from a VAR model
#'
#' @param n integer data length
#' @param mu Numeric vector of means, defaults to zero 
#' @param Sigma error covariance matrix, defaults to identity
#' @param coeffs list or matrix of VAR coefficients; model dimension and order are inferred from this
#' @param error.dist string for error distribution, one of \code{"normal", "t", "garch"}
#' @param P1 Covariance matrix for BEKK garch(1)
#' @param Q1 Autoregression matrix for BEKK garch(1)
#' @param df Integer degrees of freedom for t-distribution
#' @return data frame of simulated time series
#' @examples
#' A <- diag(0.7,4)
#' data <- VAR.sim(100, coeffs=A)
#' plot.ts(data)
VAR.sim <- function(n, mu = NULL, Sigma = NULL, coeffs, error.dist = c("normal","t","garch"), P1 = NULL, Q1 = NULL, df = 3){
  error.dist <- match.arg(error.dist, c("normal","t","garch"))
  if(is.matrix(coeffs)) coeffs <- list(coeffs)
  p <- ncol(coeffs[[1]])
  if(is.null(Sigma)) Sigma <- diag(1, p)
  if(is.null(mu)) mu <- rep(0, p)
  if(is.null(P1)) P1 <- matrix(1)
  if(is.null(Q1)) Q1 <- matrix(1)
  return( VAR_sim(n, mu, Sigma, coeffs, error.dist, P1, Q1, df) )
} 


#' Fit a piecewise VAR model to data
#'
#' @param x data matrix
#' @param cps change point vector
#' @param order integer VAR model order (optional, uses AIC otherwise)
#' @param pen penalty scalar; defaults to sSIC with exponent `1.01`
#' @return list of model list and cost 
#' @examples
#' data(voldata)
#' run_mosum <- mosumvar(voldata[,2:5], 1, 250)
#' mosumvar.fit(voldata[,2:5],run_mosum$cps)
mosumvar.fit <- function(x,cps, order=NULL, pen = log(nrow(x))^1.01 ){
  n <- nrow(x) 
  starts <- c(0, cps); ends <- c(cps, n)
  
  q <- length(cps)
  RSS <- 0
  out <- as.list(1:(q+1) )
  for (ii in 1:(q+1)) {
    out[[ii]] <- ar.ols(x[(starts[ii]+1):ends[ii],] , order.max = order, aic = is.null(order))
    #out[[ii]]$resid <- na.omit(out[[ii]]$resid) 
    #V <- out[[ii]]$var.pred
    RSS <- RSS + (ends[ii] - starts[ii]+1) *  norm( out[[ii]]$var.pred , type="F")^2 #   sum(diag( t(V) %*% V ))  
  }
  
  sSIC <- pen*q + (n/2) * log(RSS / n)
  return(list(model = out, sSIC = sSIC))
}

##### data documentation
 

#' Volatility data of five technology assets (IBM, AAPL, INTC, MSFT, ORCL), the S&P technology sector (XLK), and the S&P index (SP) 
#'
#' @name voldata
#' @docType data
#' @usage data(voldata)
NULL