

#' Segment data under a VAR model with multiple bandwidths
#'
#' @param x data matrix
#' @param order integer VAR model order
#' @param Gset integer vector of MOSUM bandwidths; see reference for default
#' @param method detector, one of \code{"Wald", "Score","BS"}
#' @param estim estimator method, one of \code{"C", "H"}
#' @param alpha Numeric significance level
#' @param criterion string location procedure
#' @param nu Numeric location procedure hyperparameter
#' @param max.iter integer maximum number of splits (when \code{method = "BS"})
#' @param do.plot Boolean, return plot
#' @return List containing
#' \itemize{  
#'   \item{\code{cps}}{ estimated change points}
#'   \item{\code{plot}}{ detector plot}
#'   \item{\code{estim}}{ input}
#' }
#' @examples
#' data(voldata)
#' mosumvar.ms(voldata[1:3000,2:4], 1)
mosumvar.ms <- function(x, order, Gset = NULL, method = c("Wald","Score","BS"),  estim = c("C","H"),  alpha = 0.05, 
                        criterion= c("eps","eta"), nu = 0.25, max.iter = 3, do.plot = TRUE){
  x <- as.matrix(x)
  p <- ncol(x)
  order <- as.integer(order)
  method <- match.arg(method, c("Wald","Score","BS"))
  estim <- match.arg(estim, c("C","H"))
  alpha <- min(1, max(0, alpha))
  if(is.null(Gset)) {
    G <- default.G(nrow(x), ncol(x), order)
    Gset <- floor(c(G, 4/3*G, 5/3*G))
  }
  out <- NULL
  if(method=="Wald") {
    out <- MFA_Wald(x,order,Gset,estim,alpha,nu)
    cps <- sort(out$ChangePoints)
  }  
  if(method=="Score"){
    mod <- ar(x, order.max = order, demean = T, method = "ols", aic = F)
    Phi <- mod$x.intercept
    eps <- mod$resid; eps[1:order,] <- 1e-4 ##solve NA
    if(p==1) Phi <- cbind(Phi,  matrix( mod$ar, nrow=ncol(x), ncol=ncol(x)))
    if(p>1){
      for (jj in 1:order){ #collect parameters into mat
        Phi <- cbind(Phi,  mod$ar[jj,,])
      }
    }
    out <- MFA_Score(x,order,Gset, Phi, eps, estim,alpha, nu=nu)
    cps <- sort(out$ChangePoints)
  }
  if(method=="BS"){
    for (ii in 1:length(Gset)) {
      cps_temp <- mosumvar.bs(x,order,Gset[[ii]],estim, alpha=alpha, criterion=criterion, nu=nu, max.iter = max.iter)$cps
        if(ii==1) {
          cps <- cps_temp 
        } else if (length(cps_temp) > 0) { 
          for (jj in 1:length(cps_temp)) {
            if (min(abs(cps_temp[jj] - cps)) >= 0.25*Gset[ii]) cps <- c(cps, cps_temp[jj])
          }
        } 
    } 
  }  
  if(do.plot){
    ts.plot(x)
    abline(v = cps, col = "red")  #if rejecting H0, add estimated cps
    pl <- recordPlot()
  } else pl <- NULL

  
  out_processed <- list(cps = cps, plot = pl, estim = estim) 
  attr(out_processed, "class") <- "mosumvar"
  return(out_processed)
}

