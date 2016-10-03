#' Title
#'
#' @param x
#' @param D
#' @param meanlog
#' @param sdlog
#'
#' @return
#' @export
#'
#' @examples
#'# Simulate a pareto sample
#' n <- 1e3 ; meanlog <- 1 ; sdlog <- 2
#' set.seed(100) ; X <- sort(rlnorm(n, meanlog, sdlog))
#' d <- 1:4
#' # Compare the ECDF with the theoritical CDF
#' plot(X,ecdf(X)(X), type = "l")
#' lines(X,ppareto(X, scale, shape), col = "red")

Dlnorm <- function(x,D,meanlog,sdlog)
{

  output <- rep(0,length(x) )
  supp <- x >0
  xs   <- x[supp]

  Xs <- (log(xs) - meanlog)/sdlog + sdlog

  if (D == 1) output[supp] <- dlnorm(xs,meanlog,sdlog)
  if (D == 2) output[supp] <- -1/(xs*sdlog)*dlnorm(xs,meanlog,sdlog)*Xs
  if (D == 3) output[supp] <- 1/(xs*sdlog)^2*dlnorm(xs,meanlog,sdlog)*( Xs^2 + sdlog*Xs - 1)
  if (D > 3) return("Not available in this package.")
  return(output)
}

#' Title
#'
#' @param k
#' @param meanlog
#' @param sdlog
#' @param lower
#'
#' @return
#' @export
#'
#' @examples
partialExectationlnorm <- function(x, meanlog,sdlog, lower = TRUE)
{
  mean <- exp(meanlog + sdlog^2/2)
  if (lower)  output <- mean*pnorm( (log(x) -meanlog)/sdlog,sdlog,1)
  if (!lower) output <- mean*(1 - pnorm( (log(x) -meanlog)/sdlog,sdlog,1))
  return(output)
}

# # Convex point of the pdf
# convpt.dlnorm = function(meanlog,sdlog) exp(meanlog + sdlog*(sqrt(sdlog^2 + 4) - 3*sdlog)/2 )
#
# UpperTruncMomLogN = function(a,order,meanlog,sdlog) #E[X^iI(x>a)]
# {
#   if (round(order)!= order)
#   {
#     print("The order must be an integer")
#     break
#   }
#
#
#   result = exp(meanlog*order +  1/2*(sdlog*order)^2)*(1-pnorm( (log(a) - meanlog)/sdlog,sdlog*order,1))
#   output = result
#
#   return(output)
# }
