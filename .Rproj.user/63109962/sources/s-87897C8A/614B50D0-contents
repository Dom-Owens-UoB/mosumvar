# mosum_lm

#' @title get H (fixed)
#' @keywords internal
get_H_fixed <- function(mod, X){
  return( as.matrix(-mod$residuals * X))
}

#' @title get A (fixed)
#' @keywords internal
get_A_fixed <- function(H, k,G){
  return(colSums(H[(k+1):(k+G),,drop=F])-colSums(H[(k-G+1):k,,drop=F]))
}

#' @title get Tk (fixed)
#' @keywords internal
get_Tk_fixed <- function(mod, k, G, H, C12 ){ 
  A <- (get_A_fixed(H,k,G)) #difference matrix
  sig <- get_local2(mod,k,G)
  ##Sigma estimator ------
  sgd <- sig^(-.5) * C12
  #------------------------------
  Tkn <-  norm(sgd %*% A,type="2") /sqrt(2*G)#(2*G)^(-.5) 
  return(Tkn)
}

#' @title get Wk (fixed)
#' @keywords internal
get_Wk_fixed <- function(x, k, G, p, C_root ){ 
  mod1 <- lm(x[(k-G+1):k,1] ~ x[(k-G+1):k,2:p] -1 ) ##left model
  mod2 <- lm(x[(k+1):(k+G),1] ~ x[(k+1):(k+G),2:p] - 1) ##right model
  
  mod1$coefficients[is.na(mod1$coefficients)] <- 0 #handle NAs
  mod2$coefficients[is.na(mod2$coefficients)] <- 0
  
  
  A <- mod2$coefficients - mod1$coefficients
  
  sig <- get_local1(mod1,mod2,G)
  ##Sigma estimator ------
  sgd <- sig^(-.5) * C_root  ##not sqrt of C
  #------------------------------
  Wkn <- sqrt(G/2) *  norm(sgd %*% A,type = "2")#as.numeric( sqrt(t(A) %*% sgd %*% A) )#
  return(Wkn)
}

#' @title Variance estimator LOCAL1
#' @keywords internal
get_local1 <- function(mod1,mod2, G){
  out <- (sum(mod1$residuals^2) + sum(mod2$residuals^2))/(2*G) ##local1 variance estim
  return(out)
}

#' @title Variance estimator LOCAL2
#' @keywords internal
get_local2 <- function(mod,k,G){
  eps <- mod$residuals
  upper <- var(eps[(k+1):(k+G)])
  lower <- var(eps[(k-G+1):(k)])
  out <- (upper + lower)* (G-1)/(2*G)
  return(out)
}

#' @title Variance estimator GLOBAL
#' @keywords internal
get_global <- function(mod) var(mod$residuals)



## main

#' Segment data under a linear regression model
#'
#' @param X matrix of covariates
#' @param y vector of responses
#' @param G integer MOSUM bandwidth
#' @param intercept include intercept in regression
#' @param method detector, one of \code{"Wald", "Score"}
#' @param alpha Numeric significance level
#' @param criterion location procedure, one of \code{"eps", "eta"}
#' @param nu Numeric location procedure hyperparameter
#' @param thresh rejection threshold; see reference for default
#' @param do.plot Boolean, return plot
#' @return List containing 
#' \itemize{
#'   \item{\code{thresh}}{ input}
#'   \item{\code{mosum}}{ vector of mosum detector values}
#'   \item{\code{cps}}{ estimated change points}
#'   \item{\code{plot}}{ detector plot}
#' }
#' @examples
#' data(voldata)
#' mosumlm(voldata[,3:5], voldata[,2])
mosumlm <- function(X, y, G = NULL, intercept = TRUE, method = c("Wald", "Score"), alpha = 0.05, criterion= c("eps","eta"), nu=0.25, thresh = NULL, do.plot = TRUE){
  X <- as.matrix(X)
  if(intercept) X <- cbind(X, 1)
  out <- NULL
  n <- dim(X)[1]
  p <- dim(X)[2]+1 
  yX <- as.matrix(cbind(y,X))
  if(is.null(G)) G <- default.G(n, sqrt(p-1), 1)
  
  method <- match.arg(method, c("Wald","Score"))
  criterion <- match.arg(criterion, c("eps","eta"))
  nu <- min(1, max(0, nu))
  alpha <- min(1, max(0, alpha)) 
  
  ##Test setup----------------------------
  if(is.null(thresh)){
    c_alpha <- -log(log( (1-alpha)^(-1/2))) #critical value
    a <- sqrt(2*log(n/G)) #test transform multipliers
    b <- 2*log(n/G) + (p-1)/2 * log(log(n/G)) - log(2/3 * gamma((p-1)/2)) ##CORRECTED
    thresh <- (b+c_alpha)/a #threshold
    thresh <- max(thresh, sqrt(2*log(n)) + c_alpha/sqrt(2*log(n)) )##ASYMPTOTIC
  }

  Reject <- FALSE
  ##Run test-----------------------------
  stat <- rep(0, n) #initialise statistic vector
  C <-t(yX[,2:p]) %*% (yX[,2:p]) / n
  ev <- eigen(C)
  ##t(ev$vectors) %*% diag( (ev$values)^(-.5) ) %*% (ev$vectors)
  
  if(method == "Wald") { 
    Croot <-  (ev$vectors) %*% diag( (ev$values)^(.5) ) %*% t(ev$vectors) 
    for (tt in (G+1):(n-G)) {
      stat[tt] <- get_Wk_fixed(yX, k=tt, G, p = p, Croot ) 
    }
  } 
  if(method == "Score"){  
    C12 <- (ev$vectors) %*% diag( (ev$values)^(-.5) ) %*% t(ev$vectors)
    mod <- lm(yX[,1] ~ yX[,2:p]-1)
    H <- get_H_fixed(mod,yX[,2:p])
    for (tt in (G+1):(n-G)) {
      stat[tt] <- get_Tk_fixed(mod,k=tt,G,H,C12)
    }
  }
  
  cps <- get_cps(stat,thresh,G, nu=nu, criterion) #locate
  if(length(cps)>0) Reject<-TRUE
  
  ##Plot------------------------------------
  if(do.plot){
    plot.ts(stat, ylab="mosum") # plot test statistic
    abline(h = thresh, col = "blue") #add threshold
    abline(v = cps, col = "red")  #if rejecting H0, add estimated cps
    pl <- recordPlot()
  #plot( a*Tn - b); abline(h=c_alpha, col="blue") #rescaled plot
  } else pl <- NULL

  ##Output------------------------------------
  out <- list(thresh = thresh, mosum = stat, cps = cps, plot = pl)
  
  return(out)
}



#' Segment data under a linear regression model on a coarse grid 
#' @param X matrix of covariates
#' @param y vector of responses
#' @param G integer MOSUM bandwidth
#' @param intercept include intercept in regression
#' @param method detector, one of \code{"Wald", "Score"}
#' @param kap Numeric subsampling resolution constant
#' @param alpha Numeric significance level
#' @param criterion location procedure, one of \code{"eps", "eta"}
#' @param nu Numeric location procedure hyperparameter
#' @param thresh rejection threshold; see reference for default
#' @param do.plot Boolean, return plot
#' @return List containing 
#' \itemize{
#'   \item{\code{thresh}}{ input}
#'   \item{\code{mosum}}{ vector of mosum detector values}
#'   \item{\code{cps}}{ estimated change points}
#'   \item{\code{plot}}{ detector plot}
#' }
#' @examples
#' data(voldata)
#' mosumlm.sub(voldata[,3:5], voldata[,2])
mosumlm.sub <- function(X, y, G = NULL, intercept = TRUE, method = c("Wald", "Score"), kap = 0.1,  alpha = 0.05, 
                            criterion= c("eps","eta"), nu=0.25, thresh = NULL, do.plot = TRUE){
  X <- as.matrix(X)
  if(intercept) X <- cbind(X, 1) 
  n <- dim(X)[1]
  p <- dim(X)[2]+1 
  yX <- cbind(y,X)
  if(is.null(G)) G <- default.G(n, sqrt(p-1), 1)
  
  method <- match.arg(method, c("Wald","Score"))
  criterion <- match.arg(criterion, c("eps","eta"))
  nu <- min(1, max(0, nu))
  alpha <- min(1, max(0, alpha)) 
  kap <- min(1, max(0, kap))
  R <- floor(kap*G)
  
  ##Test setup----------------------------
  if(is.null(thresh)){
    c_alpha <- -log(log( (1-alpha)^(-1/2))) #critical value
    a <- sqrt(2*log(n/G)) #test transform multipliers
    b <- 2*log(n/G) + (p-1)/2 * log(log(n/G)) - log(2/3 * gamma((p-1)/2)) ##CORRECTED
    thresh <- (b+c_alpha)/a #threshold
    thresh <- max(thresh, sqrt(2*log(n)) + c_alpha/sqrt(2*log(n)) )##ASYMPTOTIC
    Reject <- FALSE
  }

  ##Run test-----------------------------
  ind <- index(n,G,0,kap) #fit grid
  if(n - ind[length(ind)] < G) ind <- head(ind, -1)
  stat <- rep(0, n) #initialise statistic vector
  #if(method == "Wald")statW <- vapply(ind, get_Wkn_RCPP, FUN.VALUE = double(1), x=x,p=p,G=G, estim=estim)#get_W_RCPP(x[ind,],p,G,estim)
  C <-t(yX[,2:p]) %*% (yX[,2:p]) / n
  ev <- eigen(C)
  ##t(ev$vectors) %*% diag( (ev$values)^(-.5) ) %*% (ev$vectors)
  ##Score Subsample-------------------------------------
  # if(method == "Score" || method == "Wald" ){
    stat <- rep(0, n)  
    if(method == "Score"){  
      C12 <- (ev$vectors) %*% diag( (ev$values)^(-.5) ) %*% t(ev$vectors)
      for (k in ind) {
        window <- (k-G+1):(k+G)
        mod <- lm(yX[window,1] ~ yX[window,2:p]-1)
        H <- get_H_fixed(mod,yX[window,2:p])
        a <- get_Tk_fixed(mod,k=G,G,H,C12) ##calculate statistic on subsample
        stat[k] <- a ##collect into output vector
        if(a > thresh ){ ##if passes threshold locally
          Reject <- TRUE
          if(k>= 2*G & k <= n-1*G){ #not too close to ends
            #ss <- max(G+1, k-2*G +1); ee <- min(n-G,k+2*G) ##bounds
            newresids <- yX[,1] -  as.vector( yX[,2:p] %*% mod$coefficients )#predict(mod, as.data.frame(X[ss:ee,2:d]), se.fit=F) #rep(0,  ee-ss) #obtain extended residuals
            mod$residuals <- newresids #overwrite residuals
            newH <-  -1 * newresids *yX[,2:p] ##new estimating function series
            max_fill <- min(n-k, 2*G) #for filling in ends
            Tt <- rep(0,  max_fill) #new statistic evaluation
            for (t in  1:max_fill ){
              Tt[t] <- get_Tk_fixed(mod, k =window[t],G,newH, C12)
            }
            
            stat[window[1:max_fill]] <- pmax(stat[window[1:max_fill]], Tt) ##select max at each value
          }
        }  
        
      }
    }
    if(method == "Wald"){  
      statW <- rep(0, length(ind))
      Croot <-  (ev$vectors) %*% diag( (ev$values)^(.5) ) %*% t(ev$vectors) ##matrix root
      for (k in ind) {
        a <- get_Wk_fixed(yX,k,G,p,Croot) ##calculate statistic on subsample
        statW[which(ind==k)] <-a ##collect into output vector
      }
      stat[ind] <- statW
      test_stat <- max(stat)
      cps <- c() #empty changepoint vector
      if(test_stat > thresh){ #compare test stat with threshold
        Reject <- TRUE
        C12 <- (ev$vectors) %*% diag( (ev$values)^(-.5) ) %*% t(ev$vectors) ##need C inverse root for Tk
        times <- which(stat> thresh)#times <- sort(Reduce(union,tlist))
        tlist <- split(times, cumsum(c(1, diff(times) != R))) #split into list of consecutive regions
        
        for (i in 1:length(tlist)) {
          interval <- (max(min(tlist[[i]]) - 2*G, 1)):( min(max(tlist[[i]]) + 2*G, n)) ##extend by +-2G
          #fit var model
          mod <- lm(yX[interval,1] ~ yX[interval,2:p]-1)
          H <- get_H_fixed(mod,yX[interval,2:p])
          
          for (t in (G):(length(interval)-G )) {
            stat[min(interval) -1 + t] <- get_Tk_fixed(mod, t, G, H, C12)
          }
          #stat[interval[(G):(length(interval)-G )]] <- 
          #  vapply( (G):(length(interval)-G ), get_Tk_fixed, FUN.VALUE = 0, mod=mod, G=G,H=H,C12=C12 )  #overwrite statistic
        }
      }
    } 
    
    
    cps <- c()
    if(criterion == "eps"){
      sub_pairs <- get_sub_pairs(stat,thresh,G,kap=kap,nu=nu) #get_sub_pairs
      q <- dim(sub_pairs)[1]
      if(q==0) Reject <- FALSE
      else if (q>0){ ## locate cps
        for (ii in 1:q) {
          interval <- sub_pairs[ii,1]:sub_pairs[ii,2]
          kk <- which.max(stat[interval]) #internal cp location
          cps[ii] <- kk + sub_pairs[ii,1] #- G-p
        }
      }
    }
    if(criterion == "eta"  ){
      cps <- get_local_maxima(stat,thresh,G,nu=2*nu)
      q <- length(cps)
      if(q==0) Reject <- FALSE
    }
  # }
  ##BinSeg----------------------------------------
  # if(method=="BinSeg"){
  #   C12 <- (ev$vectors) %*% diag( (ev$values)^(-.5) ) %*% t(ev$vectors)
  #   cps <- MOSUMBS_fixed(X,1,n,D=thresh,G,p,C12,cps=c(), iter=0)
  #   if(length(cps)==0){ Reject <- F }else{ Reject <- T}
  #   stat <- c()
  # }
  
  ##Plot------------------------------------
  if(do.plot){
    plot.ts(stat, ylab="mosum") # plot test statistic
    abline(h = thresh, col = "blue") #add threshold
    if(Reject==TRUE) abline(v = cps, col = "red")  #if rejecting H0, add estimated cps
    pl <- recordPlot()
  } else pl <- NULL
  # if(method=="BinSeg"){
  #   plot.ts(X[,1], ylab="Y")
  #   if(Reject==TRUE) abline(v = cps, col = "red")  #if rejecting H0, add estimated cps
  # }
  #plot( a*Tn - b); abline(h=c_alpha, col="blue") #rescaled plot
  ##Output------------------------------------
  out <- list(thresh = thresh, mosum = stat, cps = cps, plot = pl)
  attr(out, "class") <- "mosumlm"
  return(out)
}



## binary segmentation

#' Segment data under a linear regression model using MOSUM-Binary Segmentation
#' 
#' @param X matrix of covariates
#' @param y vector of responses
#' @param G integer MOSUM bandwidth
#' @param intercept include intercept in regression
#' @param max.iter maximum number of recursions in Binary Segmentation
#' @param method detector, one of \code{"Wald", "Score"}
#' @param alpha Numeric significance level
#' @param criterion location procedure, one of \code{"eps", "eta"}
#' @param nu Numeric location procedure hyperparameter
#' @param thresh rejection threshold; see reference for default
#' @param do.plot Boolean, return plot
#' @return List containing 
#' \itemize{
#'   \item{\code{thresh}}{ input}
#'   \item{\code{mosum}}{ vector of mosum detector values}
#'   \item{\code{cps}}{ estimated change points}
#'   \item{\code{plot}}{ detector plot}
#' }
#' @examples
#' data(voldata)
#' mosumlm.bs(voldata[,3:5], voldata[,2])
mosumlm.bs <- function(X, y, G = NULL, intercept = TRUE, max.iter = 10, method = c("Wald", "Score"), alpha = 0.05, criterion= c("eps","eta"), nu=0.25, thresh = NULL, do.plot = TRUE){
  X <- as.matrix(X)
  out <- NULL 
  if(intercept) X <- cbind(X, 1) 
  n <- dim(X)[1]
  p <- dim(X)[2]+1 
  yX <- cbind(y,X)
  if(is.null(G)) G <- default.G(n, sqrt(p-1), 1)
  
  method <- match.arg(method, c("Wald","Score"))
  criterion <- match.arg(criterion, c("eps","eta"))
  nu <- min(1, max(0, nu))
  alpha <- min(1, max(0, alpha)) 
  
  ##Test setup----------------------------
  c_alpha <- -log(log( (1-alpha)^(-1/2))) #critical value
  a <- sqrt(2*log(n/G)) #test transform multipliers
  b <- 2*log(n/G) + (p-1)/2 * log(log(n/G)) - log(2/3 * gamma((p-1)/2)) ##CORRECTED
  thresh <- (b+c_alpha)/a #threshold
  thresh <- max(thresh, sqrt(2*log(n)) + c_alpha/sqrt(2*log(n)) )##ASYMPTOTIC
  Reject <- FALSE
  
  ##Run test-----------------------------
  stat <- rep(0, n) #initialise statistic vector
  C <-t(yX[,2:p]) %*% (yX[,2:p]) / n
  ev <- eigen(C)

  ##BinSeg----------------------------------------
  C12 <- (ev$vectors) %*% diag( (ev$values)^(-.5) ) %*% t(ev$vectors)
  cps <- MOSUMBS_fixed(yX,1,n,D=thresh,G,p,C12,cps=c(), iter=0)
  if(length(cps)==0){ Reject <- F }else{ Reject <- T}
  stat <- c()
  
  ##Plot------------------------------------
  if(do.plot){
    plot.ts(yX[,1], ylab="Y")
    if(Reject) abline(v = cps, col = "red")  #if rejecting H0, add estimated cps
    pl <- recordPlot()
  } else pl <- NULL
  
  ##Output------------------------------------
  out <- list(thresh = thresh, mosum = stat, cps = cps, plot = pl)
  attr(out, "class") <- "mosumlm"
  return(out)
}

#' @title mosumlm.bs inner loop
#' @keywords internal
MOSUMBS_fixed <- function(x, s, e, D, G, d, C12, cps, iter=0, max_iter = 10){
  if(e-s > 2*G & iter < max_iter){
    iter <- iter + 1 #control recursion depth
    mod <- lm(x[s:e,1] ~ x[s:e,2:d]-1)
    ##eigen
    H <- get_H_fixed(mod,x[s:e,2:d])
    stat <-  vapply((1*G):(e-s-1*G +1), get_Tk_fixed, FUN.VALUE = 0, mod=mod, G=G,H=H,C12=C12 ) 
    TT <- max(stat)
    if(TT > D){ ##passes test
      k <- which.max(stat)+G +s -1 #shift
      cps <- append(cps, k)
      cps1 <- MOSUMBS_fixed(x, s, k, D, G, d, C12, cps, iter)
      cps2 <- MOSUMBS_fixed(x, k+1, e, D, G, d, C12, cps, iter)
      cps <- union(cps, c(cps1,cps2) )
    }
  }
  return( sort(cps) ) #return cps from current iter
}
