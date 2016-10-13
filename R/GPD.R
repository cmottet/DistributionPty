
#' Gradient of the Generalized Pareto Distribution function w.r.t. to the scale and shape parameter
#'
#' @param x Vector of real values
#' @param xi, beta Shape and scale parameters of the GPD
#' @param eps Numerical tolerance on the shape parameter \code{xi}. If \code{xi < eps},
#' the GPD is considered to be an exponential distribution function with scale \code{beta}.
#' Default value is 1e-12.
#'
#' @export
#' @importFrom QRM dGPD pGPD
#'
#' @examples
#'
#' gradientpGPD(-5:5,1,1)
#'
#' # If the distribution is numerically similar to an exponential distribution
#' gradientpGPD(-5:5,0,1)
#'
gradientpGPD <- function(x,xi,beta,eps = 1e-12)
{
  gradient <- matrix(0, ncol = length(x),nrow = 2, dimnames = list(c("xi","beta")))

  # Check with eps rather than 0 for numerical stability
  if (abs(xi)>= eps){

    supp <- 0 < (1+xi*x/beta)
    xs <- x[supp]

    gradient[1,supp] <- (-1/xi^2*log(1+xi*xs/beta) + xs/(xi*(beta + xi*xs)))*(1-pGPD(xs,xi,beta))
    gradient[2,] <- -x/beta*dGPD(x,xi,beta) # Partial Beta
  }

  # In this case, the GPD fit is numerically an exponential distribution
  if (abs(xi) < eps){
    gradient[1,] <- 0
    gradient[2,] <- x/beta*dexp(x,1/beta)
  }

  output <- gradient
  return(output)
}



#' Build (1-alpha) confidence intervals of functions of GPD parameters
#'
#' This function applies the delta method to derive a (1-alpha) confidence interval
#' of the function h(xi, beta) where xi, and beta are respectively the MLE estimators
#' of the shape and scale parameters of a GPD distribution.
#'
#' @param fitGPD An object obtained using the function \emph{fit.GPD} available in the
#' package \emph{QRM} available on CRAN
#' @param h A function taking for argument the shape \emph{xi}, and scale \emph{beta} parameters
#'  of a GPD distribution
#' @param hGrad The gradient of the function h w.r.t to \emph{xi}, and \emph{beta}, in that order
#' @param alpha The desired level of accuracy for the confidence interval
#' @param verbose A logical value indicating wether messages should be printed on the screen
#'
#' @return a data frame containing three values
#' \item{lB}{the lower bound of the confidence interval}
#' \item{hHat}{the estimated value of h(xi,beta)}
#' \item{uB}{the upper bound of the confidence interval}
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
asymptoticCIforGPDfit <- function(fitGPD,h,hGrad, alpha = 0.05,verbose = TRUE)
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

    CI <- hHat + qnorm(c(alpha/2,1-alpha/2))*sdHat
    output <- (1 - Fnu)*data.frame(lB = CI[1], hHat = hHat, uB = CI[2])
  }


  return(output)
}

