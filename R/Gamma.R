#' Derivative w.r.t. x of the gamma distribution function
#' @param x vector of quantiles
#' @param D an integer between 0 and 3. The derivatives are not available
#'  in this package when D > 3
#' @param shape,rate shape and rate of the gamma distribution
#' with default values of 0 and 1 respectively.
#' @return Dgamma gives the  distribution function when D = 0, the density
#'  when D  = 1, the derivative of the density when D = 2, and so on
#' @export
#'
#' @examples
#' shape <- 2
#' rate <- 1
#' x <- seq(0,10,by = 0.01)
#'
#' # Compare to the thruth when D = 0
#' plot(x,pgamma(x,shape,rate),type="l")
#' lines(x,Dgamma(x,D = 0,shape,rate),col="red")
#'
#' # Compare to the thruth when D = 1
#' plot(x,dgamma(x,shape,rate),type="l")
#' lines(x,Dgamma(x,D = 1,shape,rate),col="red")

#' # Check that the density has the right behavior when D = 2
#' plot(x,dgamma(x,shape,rate),type="l",ylim = c(-1,1))
#' lines(x,Dgamma(x,D = 2,shape,rate),col="red")

Dgamma <- function(x,D,shape,rate)
{
  if (round(D) != D || D < 0) {
    print("D must be a non-negative integer.")
    return(NaN)
  }

  if (D == 0) output <- pgamma(x,shape,rate)
  if (D == 1) output <- dgamma(x,shape,rate)
  if (D == 2) output <- dgamma(x,shape-1,rate)*rate*(1-rate*x/(shape-1))

  if (D > 3)  {
    print("Not available for this package version")
    output[supp]  <- NA
  }

  return(output)
}

