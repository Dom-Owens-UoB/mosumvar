#' @title get change point estimates
#' @keywords internal
get_cps <- function(stat, thresh, G, nu = 1/4, criterion = c("eps","eta") ){
  cps1 <- cps2 <- c() ##assign empty cps
  if(criterion =="eps"){
    sub_pairs <- get_sub_pairs(stat,thresh,G,kap=1,nu=nu) #get_sub_pairs
    q <- dim(sub_pairs)[1]
    if(q==0) Reject <- FALSE
    else if (q>0){ ## locate cps
      for (ii in 1:q) {
        interval <- sub_pairs[ii,1]:sub_pairs[ii,2]
        kk <- which.max(stat[interval]) #internal cp location
        cps1[ii] <- kk + sub_pairs[ii,1] #- G-order
      }
    }
  }
  if(criterion =="eta"){
    cps2 <- get_local_maxima(stat,thresh,G,nu=2*nu)
    q <- length(cps2)
    if(q==0) Reject <- FALSE
  }
  cps <- union(cps1,cps2) ##output union
  return(cps)
}


#' @title Wald test
#' @keywords internal
wald <- function(x, order, G, alpha = 0.05, estim="DiagC", criterion = "eps", nu=1/4, do.bootstrap = FALSE, n.bootstrap = 1000, thresh = NULL, do.plot = TRUE){
  n <- dim(x)[1] #dimensions
  p <- dim(x)[2]
  dim_warning(n,G,p,order,"Wald")
  ##Test setup----------------------------
  if(is.null(thresh)){
  c_alpha <- -log(log( (1-alpha)^(-1/2))) #critical value
  a <- sqrt(2*log(n/G)) #test transform multipliers
  b <- 2*log(n/G) + p*(p*order+1)/2 * log(log(n/G)) - log(2/3 * gamma(p*(p*order+1)/2)) ##CORRECTED
  thresh <- (b+c_alpha)/a #threshold
  thresh <- max(thresh, sqrt(2*log(n)) + c_alpha/sqrt(2*log(n)) )##ASYMPTOTIC
  }

  Reject <- FALSE
  ##Run test-----------------------------
  Wn <- (get_W_RCPP(x,order,G,estim) ) #evaluate statistic at each time k
  test_stat <- max(Wn)
  cps <- c() #empty changepoint vector
  if(do.bootstrap){
    boot <- bootstrap.sim(x, order = order, G = G, uni = FALSE, method = "Wald",
                          estim = estim,
                          n.bootstrap = n.bootstrap)
    thresh <- quantile(boot, 1 - alpha)
  }
  if(test_stat > thresh){ #compare test stat with threshold
    Reject <- TRUE
    cps <- get_cps(Wn,thresh,G, nu=nu, criterion)
    if( is.null(cps) ) Reject <- FALSE #doesn't pass nu-test
  }
  ##Plot------------------------------------
  if(do.plot){
    plot.ts(Wn, ylab= "mosum") # plot test statistic
    abline(h = thresh, col = "blue") #add threshold
    if(Reject==TRUE) abline(v = cps, col = "red")  #if rejecting H0, add estimated cps
    pl <- recordPlot()
    #plot( a*Tn - b); abline(h=c_alpha, col="blue") #rescaled plot
  } else pl <- NULL

  ##Output------------------------------------
  out <- list(thresh = thresh, mosum = Wn, cps = cps, plot = pl, estim=estim)
  return(out)
}

