
#' @title get grid indices
#' @keywords internal
index <- function(n, G, order, kap){ #define grid
  a <- 1:n
  R <- floor(kap*G)
  b <- a[seq.int(G+order, n-G, R)]
  return(b)
}
 

#' @title get significant pairs (subsample)
#' @keywords internal
get_sub_pairs <- function(Tn, thresh, G, kap, nu = 0.25){
  n <- length(Tn)
  R <- 1#floor(kap*G)
  rshift <- c(Tn[-(1:R)],rep(0,R)); lshift <- c(rep(0,R),Tn[- ((n-R+1):n)]);
  over <- (Tn >thresh) #indices greater than thresh
  v <- which(over & lshift < thresh )#& lshift >0 ) #lowers
  w <- which(over & rshift < thresh )#& rshift >0) #uppers
  nu_remove <- which(w-v >= nu*R*G) #nu(epsilon) test for distance between
  v_nu <- v[nu_remove]; w_nu <-w[nu_remove] #c(w[nu_remove],n)
  sub_pairs <- cbind(v_nu,w_nu)
  return(sub_pairs)
  #plot(lshift); lines(Tn)
}

#' @title get local maxima
#' @keywords internal
get_local_maxima <- function(Tn, thresh, G, nu = 0.25) {
  n <- length(Tn)
  cps <- c()
  window <- floor(nu*G)
  for(t in (G+1):(n-G)){#all(Tn[(t - window):(t+window)]  > thresh) &
    if( Tn[t] > thresh & Tn[t] == max(Tn[(t - window):(t+window)]) ){cps <- append(cps,t)} ##add to list
  }
  return(cps)
}



#' Segment data under a VAR model on a coarse grid
#'
#' @param x data matrix
#' @param order integer VAR model order
#' @param G integer MOSUM bandwidth; see reference for default
#' @param method detector, one of \code{"Wald", "Score"}
#' @param estim estimator method, one of \code{"C", "H"}
#' @param var.estim variance estimator method, one of \code{"Local", "Global"}
#' @param kap Numeric subsampling resolution constant
#' @param alpha Numeric significance level
#' @param criterion location procedure, one of \code{"eps", "eta"}
#' @param nu Numeric location procedure hyperparameter
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
#' mosumvar.sub(voldata[,2:5], 1, 250)
mosumvar.sub <- function(x, order, G = NULL, method = c("Wald","Score"), estim = c("C","H"), var.estim = c("Local","Global"),
                      kap = 0.1,  alpha = 0.05, criterion= c("eps","eta"), nu=0.25){
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  order <- as.integer(order)
  method <- match.arg(method, c("Wald","Score"))
  estim <- match.arg(estim, c("C","H"))
  var.estim <- match.arg(var.estim, c("Local","Global"))
  alpha <- min(1, max(0, alpha))
  kap <- min(1, max(0, kap))
  criterion <- match.arg(criterion, c("eps","eta"))
  nu <- min(1, max(0, nu))
  
  out <- NULL
  
  if(is.null(G)) G <- default.G(n, p, order)
  dim_warning(n,G,p,order, method)

  R <- floor(kap*G)
  ##Test setup----------------------------
  c_alpha <- -log(log( (1-alpha)^(-1/2))) #critical value
  a <- sqrt(2*log(n/G)) #test transform multipliers
  b <- 2*log(n/G) + p*(p*order+1)/2 * log(log(n/G)) - log(2/3 * gamma(p*(p*order+1)/2)) ##CORRECTED
  thresh <- (b+c_alpha)/a #threshold
  thresh <- max(thresh, sqrt(2*log(n)) + c_alpha/sqrt(2*log(n)) )##ASYMPTOTIC
  Reject <- FALSE
  ##Run test-----------------------------
  ind <- index(n,G,order,kap) #fit grid
  stat <- rep(0, n) #initialise statistic vector
  if(method == "Wald"){
    statW <- vapply(ind, get_Wkn_RCPP, FUN.VALUE = double(1), x=x,p=order,G=G, estim=estim)#get_W_RCPP(x[ind,],p,G,estim)
    stat[ind] <- statW
    #return(stat)
  }  
  if(method == "Score"){
    mod <- ar.ols(x,aic=F,order.max = order,demean = T)
    ##build parameter matrix
    if(order==1) a <- cbind(mod$x.intercept, matrix( mod$ar, nrow=p, ncol=p))
    if(order>1) {
      a <- mod$x.intercept
      for (pp in 1:order) {
        a <- cbind(a,matrix( mod$ar[pp,,], nrow=p, ncol=p) ) ##append lagged parameters
      }
    }
    eps <- mod$resid ; eps[1:order,] <- 1e-4 ##solve NA
    stat <- get_T_RCPP( x, order,G-order-1, Phi = a, eps, PhiList = list(), estim,var_estim = var.estim, R=R) ##calculate statistic on subsample
    # for (k in ind) {
    #   s <- max(k-G, 1); e <- min(n,k+G) #start and end of intervals
    #   subsample <- s:e #set subsample indices 
    #   if(stat[k] > thresh ){ ##if passes threshold locally
    #     Reject <- TRUE
    #     if(k> 2*G & k <= n-1*G){ #not too close to ends   
    #       ss <- max(G+p+1, s-G+1-p); ee <- min(n,e+G) ##bounds
    #       mod_local <- ar.ols(x[ss:ee,],aic=F,order.max = order,demean = T)
    #       newresids <- matrix(0, nrow = ee-ss, ncol = p) #obtain extended residuals
    #       for(t in 1:(ee-ss) ){
    #         newresids[t,] <- predict(mod_local,x[ss - 1 + (t-order):t,], se.fit=F)
    #       }
    #       Tt <-  get_T_RCPP(x[ss:ee,], order,G-order-1,a,newresids, PhiList = list(), estim,var_estim = var.estim) ##evaluate G-window locally
    #       stat[subsample] <- pmax(stat[subsample], Tt[G:(3*G)]) ##select max at each value
    #     }
    #   }
    # } 
  }  
  test_stat <- max(stat)
  cps <- c() #empty changepoint vector
  
  if(test_stat > thresh){ #compare test stat with threshold
      Reject <- TRUE
      times <- which(stat> thresh)#times <- sort(Reduce(union,tlist))
      tlist <- split(times, cumsum(c(1, diff(times) != R))) #split into list of consecutive regions

      for (ii in 1:length(tlist)) {
        interval <- (max(min(tlist[[ii]]) - 2*G, 1)):( min(max(tlist[[ii]]) + 2*G, n)) ##extend by +-2G
        #fit var model
        mod <- ar(x[interval,], order.max = order, demean = T, method = "ols", aic = F)
        mod_a <- mod$x.intercept
        eps <- mod$resid; eps[1:order,] <- 1e-4 ##solve NA
        if(order==1) mod_a <- cbind(mod_a,  matrix( mod$ar, nrow=p, ncol=p))
        if(order>1){
          for (jj in 1:order){ #collect parameters into mat
            mod_a <- cbind(mod_a,  mod$ar[jj,,])
          }
        }
        #interval.sort <- sort( union(interval-p,interval))
        #stat[(min(tlist[[i]])):(max(tlist[[i]]))] <- get_T_RCPP( as.matrix(x[interval,]),p,G,Phi= as.matrix(mod_a), eps=as.matrix(eps),estim = estim) #overwrite statistic
        stat[interval[(1*G):(length(interval)-1*G )]] <- get_T_RCPP(
          as.matrix(x[interval,]),order,G,Phi= as.matrix(mod_a), eps=as.matrix(eps), 
          PhiList = list(), estim = estim)[(1*G):(length(interval)-1*G )] #overwrite statistic
      }
    }
  # }

  

  
  # if(method == "Score"){
  #   #statW <- rep(0, length(ind))
  #   for (k in ind){
  #     s <- max(k-G, 1); e <- min(n,k+G) #start and end of intervals
  #     subsample <- s:e #set subsample indices
  #     mod <- ar.ols(x[subsample,],aic=F,order.max = order,demean = T)
  #     ##build parameter matrix
  #     if(order==1) a <- cbind(mod$x.intercept, matrix( mod$ar, nrow=p, ncol=p))
  #     if(order>1) {
  #       a <- mod$x.intercept
  #       for (pp in 1:order) {
  #         a <- cbind(a,matrix( mod$ar[pp,,], nrow=p, ncol=p) ) ##append lagged parameters
  #       }
  #     }
  #     eps <- mod$resid ; eps[1:order,] <- 1e-4 ##solve NA
  #     Tk <- get_T_RCPP( x[subsample,], order,G-order-1, Phi = a, eps, PhiList = list(), estim,var_estim = var.estim) ##calculate statistic on subsample
  #     stat[k] <- max(Tk) ##collect into output vector
  #     if(stat[k] > thresh ){ ##if passes threshold locally
  #       Reject <- TRUE
  #       if(k> 2*G & k <= n-1*G){ #not too close to ends
  #         ss <- max(G+p+1, s-G+1-p); ee <- min(n,e+G) ##bounds
  #         newresids <- matrix(0, nrow = ee-ss, ncol = p) #obtain extended residuals
  #         for(t in 1:(ee-ss) ){
  #           newresids[t,] <- predict(mod,x[ss - 1 + (t-order):t,], se.fit=F)
  #         }
  #         Tt <-  get_T_RCPP(x[ss:ee,], order,G-order-1,a,newresids, PhiList = list(), estim,var_estim = var.estim) ##evaluate G-window locally
  #         stat[subsample] <- pmax(stat[subsample], Tt[G:(3*G)]) ##select max at each value
  #       }
  #     }
  #   }
  # }
  
  
  
  

  cps1 <- cps2 <- c() ##assign empty cps
  if(criterion == "eps" ){
    sub_pairs <- get_sub_pairs(stat,thresh,G,kap=kap,nu=nu) #get_sub_pairs
    q <- dim(sub_pairs)[1]
    if(q==0) Reject <- FALSE
    else if (q>0){ ## locate cps
      for (ii in 1:q) {
        interval <- sub_pairs[ii,1]:sub_pairs[ii,2]
        kk <- which.max(stat[interval]) #internal cp location
        cps1[ii] <- kk + sub_pairs[ii,1] #- G-p
      }
    }
  }
  if(criterion == "eta" ){
    cps2 <- get_local_maxima(stat,thresh,G,nu=2*nu)
    q <- length(cps2)
    if(q==0) Reject <- FALSE
  }
  cps <- union(cps1,cps2) ##output union

  ##Plot------------------------------------
  plot.ts(stat, ylab="Statistic") # plot test statistic
  abline(h = thresh, col = "blue") #add threshold
  if(Reject==TRUE) abline(v = cps, col = "red")  #if rejecting H0, add estimated cps
  pl <- recordPlot()
  #plot( a*Tn - b); abline(h=c_alpha, col="blue") #rescaled plot
  ##Output------------------------------------
  out <- list(Reject = Reject, Threshold = thresh, mosum = stat, cps = cps, plot = pl, estim=estim)
  attr(out, "class") <- "mosumvar"
  return(out)
} 




