

#' @title score test
#' @keywords internal
score <- function(x, order, G, Phi=NULL, eps=NULL, alpha = 0.05, estim="C",var.estim = "Local", criterion="eps", nu=0.25, 
                           do.bootstrap = c(F,"multiplier","regression")[1], n.bootstrap=1000, thresh = NULL, do.plot =TRUE){
  n <- dim(x)[1] #dimensions
  p <- dim(x)[2]
  dim_warning(n,G,p,order,"Score")

  # if(is.null(dim(x)) || dim(x)[2] == 1 ) {
    x <- as.matrix(x); 
    if(!is.null(Phi)) Phi <- as.matrix(Phi); 
  # } #handle univariate case
    if(!is.null(eps)) eps <- as.matrix(eps)
  
  ## Estimate model
  if( is.null(Phi) || is.null(eps) ){
    mod <- ar(x, order.max = order, demean = T, method = "ols", aic = F)
    Phi <- mod$x.intercept
    eps <- mod$resid; 
    eps <- as.matrix(eps)
    eps[1:order,] <- rnorm(order)*1e-4 ##solve NA
    if(order==1) Phi <- cbind(Phi,  matrix( mod$ar, nrow=p, ncol=p))
    if(order>1){
      for (jj in 1:order){ #collect parameters into mat
        Phi <- cbind(Phi,  mod$ar[jj,,])
      }
    }
  }
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
  if(p==1){
    Tn <- ts(get_T_univ(z = as.matrix(x), x = as.matrix(x),  p = order, G=G, Phi=as.matrix(Phi), eps=as.matrix(eps),
                        PhiList = list(Phi), estim=estim, var_estim = var.estim))
  } else  Tn <- ts(get_T_RCPP(as.matrix(x),order,G,as.matrix(Phi), as.matrix(eps),PhiList = list(), estim,var.estim)) #evaluate statistic at each time k
  test_stat <- max(Tn)
  cps <- c() #empty changepoint vector
  if(test_stat > thresh){ #compare test stat with threshold
    Reject <- TRUE
    cps <- get_cps(Tn,thresh,G, nu=nu, criterion)
    if( is.null(cps) ) Reject <- FALSE #doesn't pass nu-test
  }
  ##Multiplier Bootstrap--------------------------
  if(do.bootstrap){
    boot <- bootstrap.sim(x, order = order, G = G, uni = FALSE, method = "Score",
                          estim = estim, var.estim = var.estim,
                          n.bootstrap = n.bootstrap)
    thresh <- quantile(boot, 1 - alpha)
  }
  ##Plot------------------------------------
  if(do.plot){
  plot(Tn, ylab = "mosum") # plot test statistic
  abline(h = thresh, col = "blue") #add threshold
  if(Reject==TRUE) abline(v = cps, col = "red")  #if rejecting H0, add estimated cps
  pl <- recordPlot()
  #plot( a*Tn - b); abline(h=c_alpha, col="blue") #rescaled plot
  } else pl <- NULL
  ##Output------------------------------------
  out <- list(thresh = thresh, mosum = Tn, cps = cps, plot = pl, estim = estim)
  return(out)
} 
