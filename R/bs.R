#' @title binary segmentation inner loop
#' @keywords internal
MBS_RECUR <- function(x, order, p, s, e, thresh, G, estim = "C", var_estim = "Local", cps, stat = list(), iter =1, criterion="eps", 
                      nu=.25, alpha = 0.05){
  if(e-s>2*G + order && iter < length(stat)){ ## segment is long enough, and recursion is shallow
    iter <- iter +1
    mod <- ar.ols(x, aic=F, order.max = order)
    mod_a <- mod$x.intercept
    eps <- as.matrix(mod$resid); eps[1:order,] <- 1e-4 ##solve NA
    if(order==1) mod_a <- cbind(mod_a,  matrix( mod$ar, nrow=p, ncol=p))
    if(order>1){
      for (jj in 1:order){ #collect parameters into mat
        mod_a <- cbind(mod_a,  mod$ar[jj,,])
      }
    }
    sub_vector <- get_T_RCPP(x[s:e,, drop=FALSE], order,G-order-1,mod_a,eps, PhiList = list(), estim,var_estim) ##calculate statistic on subsample ##index here
    stat[[iter]][(s+G):(e-G)] <- sub_vector[G:(e-s-G),]
    statW <- max(sub_vector) ##collect into output vector ##[which(ind==k)]
    if(is.null(thresh)){
      n <- e-s+1
      c_alpha <- -log(log( (1-alpha)^(-1/2))) #critical value
      a <- sqrt(2*log(n/G)) #test transform multipliers
      b <- 2*log(n/G) + p*(p*order+1)/2 * log(log(n/G)) - log(2/3 * gamma(p*(p*order+1)/2))
      thresh <- (b+c_alpha)/a #threshold
      thresh <- max(thresh, sqrt(2*log(n)) + c_alpha/sqrt(2*log(n)) )
    }
    if(statW> thresh){ #test
      cps_se <- get_cps(sub_vector[,], thresh, G, nu=nu, criterion = criterion)# get_local_maxima(sub_vector, thresh, G, nu=nu) #MOSUM locate change points in interval
      if(length(cps_se)>0){
        cps_se <- c(s, s+cps_se, e) #add start and end points
        for(jj in 2:length(cps_se)){ #BS on each interval
          cps_jj <- MBS_RECUR(x,order, p, s= cps_se[jj-1], e= cps_se[jj], thresh = NULL, G = G, 
                              estim,var_estim, cps_se, stat = stat, iter = iter, nu=nu ) ##RECURSION
          stat <- cps_jj$stat
          cps_se <- union(cps_se, cps_jj$cps) ##add change points
          for(cp in cps_se){
            if(min(abs(cp-cps)) > nu * G) cps <- c(cps, cp) ##add change points to output
          }
        }
      }
      
    }
    #else break
  }
  return(list(cps = sort(cps), stat=stat, thresh = thresh))
}

#' Segment data under a VAR model with MOSUM-Binary Segmentation
#'
#' @param x data matrix
#' @param order integer VAR model order
#' @param G integer MOSUM bandwidth
#' @param estim string estimation method
#' @param var.estim string variance estimation method
#' @param alpha Numeric significance level
#' @param criterion string location procedure
#' @param nu Numeric location procedure hyperparameter
#' @param max.iter integer maximum number of splits
#' @param do.plot Boolean, return plot
#' @return List containing
#' \itemize{ 
#'   \item{\code{mosum}}{ vector of mosum detector values}
#'   \item{\code{cps}}{ estimated change points}
#'   \item{\code{plot}}{ detector plot}
#'   \item{\code{estim}}{ input}
#' }
#' @examples
#' data(voldata)
#' mosumvar.bs(voldata[,2:5], 1, 250)
mosumvar.bs <- function(x, order, G = NULL, estim = c("C","H"), var.estim = c("Local","Global"), alpha = 0.05, 
                        criterion= c("eps","eta"), nu = 0.25, max.iter = 3, do.plot = TRUE){
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  order <- as.integer(order) 
  estim <- match.arg(estim, c("C","H"))
  var.estim <- match.arg(var.estim, c("Local","Global"))
  alpha <- min(1, max(0, alpha)) 
  criterion <- match.arg(criterion, c("eps","eta"))
  nu <- min(1, max(0, nu))
  
  if(is.null(G)) G <- default.G(n, p, order)
  dim_warning(n,G,p,order, "Score")
  ##Test setup----------------------------
  Reject <- FALSE
  ##Run test-----------------------------
  stat <- list(1:max.iter) #initialise statistic vector
  for (ii in 1:max.iter) {
    stat[[ii]] <- rep(0,n)
  }
  ##
  callBS <- MBS_RECUR(x,order,p,s=1,e=n,thresh=NULL,G = G,cps = c(), stat=stat, iter=0, criterion = criterion, nu = nu, alpha = alpha)
  cps <- callBS$cps[-c(1, length(callBS$cps))]
  #par(mfrow = c(max.iter,1))
  if(do.plot){
    plot.ts(callBS$stat[[1]], ylab="mosum") # plot test statistic
    abline(h = callBS$thresh[[1]], col = 2)
    for (ii in 2:max.iter) {
      lines(callBS$stat[[ii]], col = ii+1)
      # abline(h = callBS$thresh[[ii]], col = ii+1, lty = 2)
    }
    legend("topright", legend =1:max.iter,title="Iteration", fill = c(1,2:max.iter +1) ) #add threshold
    if(length(cps) > 0) abline(v = cps, col = "red")  #if rejecting H0, add estimated cps
    pl <- recordPlot()
  } else pl <- NULL
  out <- list(mosum = callBS$stat, cps = cps, plot=pl, estim = estim)
  attr(out, "class") <- "mosumvar"
  return(out)
}