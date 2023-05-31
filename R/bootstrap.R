## sim bootstrap

#' @title Simulate bootstrapped threshold
#' @keywords internal
bootstrap.sim <- function(x, order = 1, G = NULL, uni = FALSE, method = c("Wald","Score"), estim = c("C","H"), var.estim = c("Local","Global"),
                          rm.cross.terms = TRUE, global.resids = TRUE, n.bootstrap = 1000){
  n <- nrow(x)
  p <- ncol(x)
  order <- as.integer(order)
  method <- match.arg(method, c("Wald","Score"))
  estim <- match.arg(estim, c("C","H"))
  var.estim <- match.arg(var.estim, c("Local","Global"))
  
  ar <- ar.ols(x, aic = F, order.max = order)
  
  if(order>1){
    coeffs <- list()
    for (ii in 1:order) {
      coeffs[[ii]] <- as.matrix(ar$ar[ii,,])
    }
  } else {
    coeffs <- as.matrix(ar$ar[,,])
  }
  
  out <- rep(0, n.bootstrap)
  Sigma <- ar$var.pred #diag(1, p)
  for (ii in 1:n.bootstrap) {
    if(ii %% floor(n.bootstrap/10) == 0)  print(paste0("Bootstrapping threshold: ", ii, "/", n.bootstrap))
    x_ii <- mosumvar::VAR.sim(n, mu = ar$x.intercept, Sigma = Sigma, coeffs = coeffs)
    if(uni){
      ms_ii <- mosumvar::mosumvar.uni(x_ii, G = G, order = order, method = method, 
                            estim = estim, var.estim = var.estim,
                            rm.cross.terms = rm.cross.terms, global.resids = global.resids, do.plot = FALSE)
    } else {
      ms_ii <- mosumvar::mosumvar(x_ii, G = G, order = order, method = method,  estim = estim, var.estim = var.estim, do.plot = FALSE )
    }
    out[ii] <- max(ms_ii$mosum)
  }
  return(out)
}
 
