
#' Segment data under univariate AR models, optionally using dimension reduction 
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
#' @param rm.cross.terms Boolean, perform dimension reduction
#' @param do.bootstrap Boolean, determine threshold via bootstrap method
#' @param n.bootstrap integer number of simulations for bootstrap
#' @param global.resids Boolean, use residuals from full VAR model
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
#' mosumvar.uni(voldata[,2:5], 1, 250)
mosumvar.uni <- function(x, order, G = NULL,  method = c("Wald","Score"), estim = c("C","H"), var.estim = c("Local","Global"),  alpha = 0.05,  
                         criterion= c("eps","eta"), nu=0.25,
                         rm.cross.terms =FALSE, do.bootstrap = FALSE, n.bootstrap = 1000, global.resids = FALSE, thresh = NULL, do.plot = TRUE){
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  order <- as.integer(order)
  method <- match.arg(method, c("Wald","Score"))
  estim <- match.arg(estim, c("C","H"))
  var.estim <- match.arg(var.estim, c("Local","Global"))
  alpha <- min(1, max(0, alpha))
  criterion <- match.arg(criterion, c("eps","eta"))
  nu <- min(1, max(0, nu))
  if(is.null(G)) G <- default.G(nrow(x), ncol(x), order)
  if(global.resids) eps_global <- ar(x, order.max = order, demean = T, method = "ols", aic = F)$resid
  if(rm.cross.terms) {
    z <- remove_cross_terms(x, order, p)
    Phi <- as.list(1:p)
    for (ii in 1:p) {
      xi <- c()
      for (jj in 1:order) {
        interv <- (order-1+jj):(n-jj)
        xi <- cbind(xi,x[interv,ii])
      }
      Phi[[ii]] <- t(as.matrix(lm(z[-(1:order),ii] ~ xi)$coefficients))
    }
    eps <- (ar(x, order.max = order, demean = T, method = "ols", aic = F)$resid)
  } else {
    ##Test setup----------------------------
    z <- x
    mod <- Phi <- as.list(1:p)
    eps <- matrix(0, nrow = n, ncol = p)
    for(ii in 1:p){ ##fit univariate score models
      mod[[ii]] <- ar(x[,ii], order.max = order, demean = T, method = "ols", aic = F)
      Phi[[ii]] <- mod[[ii]]$x.intercept
      if(global.resids){
        eps[,ii] <- eps_global[,ii]
      } else {
        eps[,ii] <- as.matrix(mod[[ii]]$resid);
      }
      eps[,ii][1:order] <- rnorm(order)*1e-4 ##solve NA
      if(order==1) Phi[[ii]] <- matrix(c(Phi[[ii]],mod[[ii]]$ar[,,1]),1,2 ) #cbind(Phi,  matrix( mod$ar, nrow=d, ncol=d))
      if(order>1){ 
        for (jj in 1:order){ #collect parameters into mat
          Phi[[ii]] <- cbind(Phi[[ii]],  mod[[ii]]$ar[jj,,])
        } 
      } 
    }
  }
  if(is.null(thresh)){
    c_alpha <- -log(log( (1-alpha)^(-1/2))) #critical value
    a <- sqrt(2*log(n/G)) #test transform multipliers
    b <- 2*log(n/G) + (p*order+1)/2 * log(log(n/G)) - log(2/3 * gamma( (p*order+1)/2)) ##CORRECTED
    thresh <- (b+c_alpha)/a #threshold
    thresh <- max(thresh, sqrt(2*log(n)) + c_alpha/sqrt(2*log(n)) )##ASYMPTOTIC
  }
  
  Reject <- FALSE
  ##Run test-----------------------------
  stat <- rep(0, n) #initialise statistic vector
  if(method == "Wald"){
    for (ii in 1:p) {
      stat <- stat +  get_W_RCPP( as.matrix(z[,ii]),order,G,estim)#
    }
  }
  if(method == "Score") {
    for (ii in 1:p) {
      stat <- stat +  get_T_RCPP( as.matrix(z[,ii]),order,G,as.matrix(Phi[[ii]]), as.matrix(eps[,ii]),PhiList = list(), estim,var.estim)#
    }
  }
    
  if(do.bootstrap){
    boot <- bootstrap.sim(x, order = order, G = G, uni = TRUE, method = method,
                          estim = estim, var.estim = var.estim, rm.cross.terms = rm.cross.terms, global.resids = global.resids,
                          n.bootstrap = n.bootstrap)
    thresh <- quantile(boot, 1 - alpha)
  }
  
  cps <- c() #empty changepoint vector
  if(max(stat) > thresh){ #compare test stat with threshold
    Reject <- TRUE
    cps <- get_cps(stat,thresh,G, nu=nu, criterion)
    if( is.null(cps) ) Reject <- FALSE #doesn't pass nu-test
  } 
  ##Multiplier Bootstrap--------------------------
  # if (do.bootstrap){
  #   if( is.null(cps)) cps <- c(0)
  #   mbs <- multiplier_bootstrap(z, z, order, G, Phi, eps, cps, L = floor(n/4), n.bootstrap, estim, var.estim)
  #   thresh <- quantile(mbs, 1-alpha) ##overwrite threshold with bootstrap quantile 
  #   if(max(stat) > thresh){ #compare test stat with new threshold
  #     Reject <- TRUE
  #     cps <- get_cps(stat,thresh,G, nu=nu, criterion)
  #     if( is.null(cps) ) Reject <- FALSE #doesn't pass nu-test
  #   } 
  # }

  ##Plot------------------------------------
  if(do.plot){
    plot.ts(stat, ylab = "Statistic") # plot test statistic
    abline(h = thresh, col = "blue") #add threshold
    if(Reject==TRUE) abline(v = cps, col = "red")  #if rejecting H0, add estimated cps
    pl <- recordPlot()
    #plot( a*Tn - b); abline(h=c_alpha, col="blue") #rescaled plot
  } else pl <- NULL
  
  ##Output------------------------------------
  out <- list(Reject = Reject, Threshold = thresh, mosum = stat, cps = cps, plot = pl, estim = estim)
  attr(out, "class") <- "mosumvar"
  return(out)
}

#' @title Remove cross terms
#' @keywords internal
remove_cross_terms <- function(x,order,p){
  mod <- ar(x, order.max = order, demean = T, method = "ols", aic = F)
  
  xlist <- as.list(1:order) 
  for (jj in 1:order){ #
    aa <- matrix(mod$ar[jj,,], p,p)
    diag(aa) <- rep(0, p)
    xt <-  x %*% t(aa)
    xt <- xt[-(1:jj),]
    xt <- rbind(xt, matrix(0,jj,p))
    xlist[[jj]] <- xt
  } 
  for (jj in 1:order) {
    x <- x - xlist[[jj]]
  }
  return(x)
  
} 
