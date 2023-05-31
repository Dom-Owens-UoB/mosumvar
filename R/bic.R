# 
#  
# mosumvar.bic <- function(x, order, G, method = c("Wald","Score"), estim = c("C","H"), var.estim = c("Local","Global"),  
#                          uni = FALSE,
#                          alpha = 0.05, nu= 0.25, do.plot = TRUE){
#   criterion <-  "eta"
#   
#   x <- as.matrix(x)
#   order <- as.integer(order)
#   method <- match.arg(method, c("Wald","Score"))
#   estim <- match.arg(estim, c("C","H"))
#   var.estim <- match.arg(var.estim, c("Local","Global"))
#   alpha <- min(1, max(0, alpha))
#   # criterion <- match.arg(criterion, c("eps","eta"))
#   nu <- min(1, max(0, nu))
#   if(is.null(G)) G <- default.G(nrow(x), ncol(x), order, method)
#   if(uni){
#     out <- mosumvar.uni(x, order, G, method, estim, var.estim, alpha = 0, criterion, nu, thresh = 1e-6, do.plot = FALSE )
#   } else {
#     if(method== "Wald") out <- wald(x, order, G, alpha = alpha, estim= estim, 
#                                     criterion = criterion, nu=nu, thresh = 1e-6, do.plot = FALSE)
#     if(method== "Score") out <- score(x, order, G, alpha = alpha, estim= estim, 
#                                       criterion = criterion, nu=nu, thresh = 1e-6, do.plot = FALSE)
#   }
#   ranks <- rank(-out$mosum[out$cps])
#   bic_scores <- fit.model.bic(x, order, out$cps, ranks, floor(nrow(x)/(nu*G)) )
#   
#   min_q <- which.min(bic_scores) - 1
#   out_cps <- out$cps[ranks <= min_q]
#   
#   
#   ##Plot------------------------------------
#   if(do.plot){
#     par(mfrow = c(2, 1))
#     plot(1:length(bic_scores) - 1, bic_scores, xlab = "q", ylab = "BIC")
#     abline(v = min_q, col = "red") #add q
#     
#     
#     plot.ts(out$mosum, ylab="mosum") # plot test statistic
#     abline(h = out$mosum[out$cps[ranks == min_q]], col = "blue") #add threshold
#     abline(v = out_cps, col = "red")  #if rejecting H0, add estimated cps
#     pl <- recordPlot()
#   } else pl <- NULL
#   
#   
#   out_list <- list(call = out, bic_scores = bic_scores, cps = out_cps, pl = pl)
#   return(out_list)
# }
#  
#  
# fit.model.bic <- function(x, order, cps, ranks, max.cps, uni = FALSE){
#   n <- nrow(x)
#   p <- ncol(x)
#   q <- length(cps)
#   max.cps <- min(max.cps, q)
#   
#   model_list <- list()
#   complexity <- log(n)^1.01 * p^uni *(order*p+1) *  0:max.cps #
#   # fidelity <- rep(0, max.cps+1)
#   fid_mat <- array(0, dim = c(p,p,max.cps+1)) #matrix(0, p,p) #
#   cps_ <- c(0, cps, n)
# 
#   model_list[[max.cps+1]] <- list()
#   if(uni){
#     z <- remove_cross_terms(x, order, p)
#     model_list[[max.cps+1]][[1]] <- list()
#     for (ii in 1:(max.cps+1) ) {
#       for (pp in 1:p) {
#         model_list[[max.cps+1]][[ii]][[pp]] <- ar.ols(z[(cps_[ii]+1):cps_[ii+1],pp], aic = F, order.max = order)
#         fid_mat[,,max.cps+1] <- fid_mat[,,max.cps+1] + model_list[[max.cps+1]][[ii]][[pp]]$var.pred * (cps_[ii+1] - cps_[ii]+1)/n
#       }
#     }
#     for (jj in 1:(max.cps+1) ) {
#       cps_jj <- c(0, cps[ranks < jj], n)
#       
#       for (ii in 1:(length(cps_jj)-1) ) {
#         for (pp in 1:p) {
#           model_list[[jj]][[ii]][[pp]] <- ar.ols(x[(cps_[ii]+1):cps_[ii+1],pp], aic = F, order.max = order)
#           fid_mat[,,jj] <- fid_mat[,,jj] + model_list[[jj]][[ii]][[pp]]$var.pred * (cps_[ii+1] - cps_[ii]+1)/n
#         }
#       }
#       
#     }  
#     
#   } else {
#     for (ii in 1:(max.cps+1) ) {
#       model_list[[max.cps+1]][[ii]] <- ar.ols(x[(cps_[ii]+1):cps_[ii+1],], aic = F, order.max = order)
#       fid_mat[,,max.cps+1] <- fid_mat[,,max.cps+1] + model_list[[max.cps+1]][[ii]]$var.pred * (cps_[ii+1] - cps_[ii]+1)/n
#     }
#     for (jj in 1:(max.cps+1) ) {
#       cps_jj <- c(0, cps[ranks < jj], n)
#       
#       for (ii in 1:(length(cps_jj)-1) ) {
#         model_list[[jj]][[ii]] <- ar.ols(x[(cps_jj[ii]+1):cps_jj[ii+1],], aic = F, order.max = order)
#         fid_mat[,,jj] <- fid_mat[,,jj] + model_list[[jj]][[ii]]$var.pred * (cps_jj[ii+1] - cps_jj[ii]+1)/n
#       }
#       
#     }  
#     
#   }
# 
#   
#   
# 
#   
#   # model_list[[jj]] <- model_list[[jj+1]]
#   # added_cp <- ranks[jj]
#   # model_list[[jj]][[added_cp]] <- model_list[[jj]][[added_cp+1]] <- 
#   #   ar.ols(x[(cps_[added_cp]+1):cps_[added_cp+2],], aic = F, order.max = order)
#   #  
#   # fid_mat[,,jj] <- fid_mat[,,jj+1] - 
#   #   model_list[[jj+1]][[added_cp]]$var.pred* (cps_[added_cp+1] - cps_[added_cp]+1) /n - 
#   #   model_list[[jj+1]][[added_cp+1]]$var.pred * (cps_[added_cp+2] - cps_[added_cp+1]+1)/n + 
#   #   model_list[[jj]][[added_cp]]$var.pred * (cps_[added_cp+2] - cps_[added_cp]+1)/n
#   fidelity <- apply(fid_mat, 3, function(x) log(det(x)))
#   par(mfrow = c(2, 1))
#   plot( n/2 * fidelity)
#   lines(complexity)
#   bic_scores <- n/2 * fidelity + complexity
#   return(bic_scores)
# }