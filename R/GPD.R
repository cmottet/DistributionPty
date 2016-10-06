
#' Title
#'
#' @param x
#' @param xi
#' @param beta
#' @param eps
#'
#' @return
#' @export
#' @importFrom QRM dGPD pGPD
#'
#' @examples
gradientpGPD <- function(x,xi,beta,eps = 1e-12)
{
  gradient <- matrix(0, ncol = 1,nrow = 2)
  row.names(gradient) <- c("xi","beta")

  if (is.finite(x) & x>=0 & (1+xi*x/beta)>= 0 )
  {
    # Check with eps rather than 0 for numerical stability
    if (abs(xi)>= eps)
      gradient[1] <- (-1/xi^2*log(1+xi*x/beta) + x/(xi*(beta + xi*x)))*(1-pGPD(x,xi,beta))  # Partial xi

    gradient[2] <- -x/beta*dGPD(x,xi,beta) # Partial Beta
  }

  output <- gradient
  return(output)
}



#' Title
#'
#' @param x
#' @param fitGPD
#' @param h
#' @param hGrad
#' @param alpha
#' @param bonferroni
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
#' library(QRM)
#' sample <- rGPD(500, 1,1)
#'
#' # Find a suitable threshold, and fit a GPD
#' MEplot(sample)
#' u <- 0
#' fitGPD <- fit.GPD(sample, u, type = "ml")
#'
#' # Build a CI on P(4 <= X <= 9), whose true value is 10%
#' h <- function(xi, beta) pGPD(9 -u,xi,beta) - pGPD(4 -u,xi,beta)
#' hGrad <- function(xi,beta) gradientpGPD(9-u,xi,beta) - gradientpGPD(4-u,xi,beta)
#' asymptoticCIforGPDfit(fitGPD,h,hGrad )
asymptoticCIforGPDfit <- function(fitGPD,h,hGrad, alpha = 0.05,bonferroni = 1,verbose = TRUE)
{

  # Getting the parameters with the library ismev
  xi   <- as.numeric(fitGPD$par.ests[1])
  beta <- as.numeric(fitGPD$par.ests[2])
  tmp     <- fitGPD$varcov

  Fnu <- fitGPD$p.less.thresh
  nexc <- fitGPD$n.exceed

  Sigma   <- matrix( (1+xi)*c( (1+xi),-beta,-beta,2*beta^2),ncol =2,nrow = 2)

  # Point estimate of h for the given estimation of xi and beta
  hHat <- h(xi,beta)

  if (xi < -1)   output <- data.frame(lB = NA, hHat = NA, uB = NA)
  if (xi >= -1)
  {

    if (nexc  < 30 && verbose ) print(paste("Unreliable Delta-Method for u =" ,u ,". Nexc < 30",sep=""))
    if (xi< -1/2 && verbose)  print(paste("Unreliable Delta-Method for u =" ,u ,". xi < -1/2",sep=""))

    gradient  <- hGrad(xi,beta)
    sdHat <- sqrt(t(gradient)%*%Sigma%*%gradient/nexc)

    CI <- hHat + qnorm(c(alpha/(2*bonferroni),1-alpha/(2*bonferroni)))*sdHat
  }

  output <- (1 - Fnu)*data.frame(lB = CI[1], hHat = hHat, uB = CI[2])

  return(output)
}



# CI.qGPD = function(prob,fit.gpd,alpha = 0.05,bonferroni = 1)
# {
#   u = fit.gpd$threshold
#
#   # Getting the parameters with the library QRM
#   # xiHat   = fit.gpd$par.ests[1]
#   # betaHat = fit.gpd$par.ests[2]
#   # Sigma   = fit.gpd$varcov
#   #
#   # Getting the parameters with the library ismev
#   xiHat   = fit.gpd$mle[2]
#   betaHat = fit.gpd$mle[1]
#   Sigma   = fit.gpd$cov
#
#   if (0 < prob & prob < 1)
#   {
#     muHat    = qGPD(prob,xiHat,betaHat)
#     gradient = gradientqGPD(prob,xiHat,betaHat)
#
#     sdHat = sqrt(t(gradient)%*%Sigma%*%gradient)
#
#     CI = muHat + qnorm(c(alpha/(2*bonferroni),1-alpha/(2*bonferroni)))*sdHat
#     CI = pmax(0,CI) # Check that CI is smaller than 1
#   } else {muHat = NA ; CI = c(NA,NA)}
#
#   output = as.vector(c(CI[1],muHat,CI[2]))
#   return(output)
# }

# gradientdGPD = function(x,xi,beta)
# {
# # partialXi   = -(-1/xi^2*log(1+xi*x) + (1/xi+1)*x/(1+xi*x))*dGPD(x,xi,beta)
# # partialBeta = -1/beta*dGPD(x,xi,beta)
#
#   partialXi   = -(-1/xi^2*log(1+xi*x/beta) + (1/xi+1)*x/(beta+xi*x))*dGPD(x,xi,beta)
#   partialBeta = (-1/beta + x/beta*(1+xi)/(beta+xi*x))*dGPD(x,xi,beta)
#
#
# gradient = matrix(c(partialXi,partialBeta),nrow = 2)
# return(gradient)
# }
# gradientqGPD = function(prob,xi,beta)
# {
#   partialXi   = -beta/xi^2*( log(1-prob)*xi*(1-prob)^(-xi)  + (1-prob)^(-xi) - 1)
#   partialBeta = 1/xi*((1-prob)^(-xi) - 1)
#
#   gradient = matrix(c(partialXi,partialBeta),nrow = 2)
#   return(gradient)
#
# }
# gradientddGPD = function(x,xi,beta)
# {
#   grad.dGPD   = gradientdGPD(x,xi,beta)
#
# #   partialXi   = -(beta-x)/(beta+beta*x)^2*dGPD(x,xi,beta)- (xi+1)/(beta+xi*x)*grad.dGPD[1,1]
# #   partialBeta =    (xi+1)/(beta+beta*x)^2*dGPD(x,xi,beta)- (xi+1)/(beta+xi*x)*grad.dGPD[2,1]
#
#   partialXi   = -(beta-x)/(beta+xi*x)^2*dGPD(x,xi,beta) - (xi+1)/(beta+xi*x)*grad.dGPD[1,1]
#   partialBeta =    (xi+1)/(beta+xi*x)^2*dGPD(x,xi,beta) - (xi+1)/(beta+xi*x)*grad.dGPD[2,1]
#
#   gradient = matrix(c(partialXi,partialBeta),nrow = 2)
#   return(gradient)
# }
# CI.GPD = function(x,fit.gpd,order = c(-1,0,1),alpha = 0.05,bonferroni = 1)
# {
#   if (order == -1) {gradFunc = gradientpGPD  ; Func = pGPD}
#   if (order == 0)  {gradFunc = gradientdGPD  ; Func = dGPD}
#   if (order == 1)  {gradFunc = gradientddGPD ; Func = ddGPD}
#
#
#   # Getting the parameters with the library ismev
#   xiHat   = fit.gpd$mle[2]
#   betaHat = fit.gpd$mle[1]
#   tmp     = fit.gpd$cov
#
#   Sigma   = matrix(c(tmp[2,2],tmp[1,2],tmp[2,1],tmp[1,1]),ncol =2,nrow = 2)
#   u = as.numeric(fit.gpd$threshold)
#   Fn = ecdf(fit.gpd$xdata)
#
#   #
#   #   # Getting the parameters with the library extRemes
#   #   xiHat   =  as.numeric(fit.gpd$results$par[2])
#   #   betaHat =  as.numeric(fit.gpd$results$par[1])
#   #   tmp     =  solve(fit.gpd$results$hessian)
#   #   Sigma   =  matrix(c(tmp[2,2],tmp[1,2],tmp[2,1],tmp[1,1]),ncol =2,nrow = 2)
#   #
#   #   u = fit.gpd$threshold
#   #   Fn = ecdf(fit.gpd$x)
#   #
#
#
#   if (xiHat < 0)  upperbound = (u - betaHat/xiHat)
#   if (0 <= xiHat) upperbound = Inf
#
#   if (u < x & x<  upperbound)
#   {
#     ## Confidence Interval based on the Delta-Method
#     sdHat = sqrt(t(gradFunc(x-u,xiHat,betaHat))%*%Sigma%*%gradFunc(x-u,xiHat,betaHat))
#     muHat = Func(x-u,xiHat,betaHat)
#
#     Fn = ecdf(fit.gpd$xdata)
#     CI = (1 - Fn(threshold))*(muHat + qnorm(c(alpha/(2*bonferroni),1-alpha/(2*bonferroni)))*sdHat)
#
#     if (order == -1) CI = pmax(0,pmin(1,CI))
#     if (order == 0)  CI = pmax(0,CI)
#
#
#   } else CI = c(0,0)
#
#   return(CI)
# }
