#' Derivative w.r.t. x of the Normal distribution function
#' @param x vector of quantiles
#' @param mean,sd mean and standard deviation of the distribution
#' with default values of 0 and 1 respectively.
#' @return Dnorm gives the  distribution function when D = 0, the density
#'  when D  = 1, the derivative of the density when D = 2, and so on
#' @export
#'
#' @examples
#' mean <- 1
#' sd <- 1/2
#' x <- seq(-4,4,by = 0.01)
#'
#' # Compare to the thruth when D = 0
#' plot(x,pnorm(x,mean,sd),type="l")
#' lines(x,Dnorm(x,D = 0,mean,sd),col="red")
#'
#' # Compare to the thruth when D = 1
#' plot(x,dnorm(x,mean,sd),type="l")
#' lines(x,Dnorm(x,D = 1,mean,sd),col="red")

#' # Check that the density has the right behavior when D = 2
#' plot(x,dnorm(x,mean,sd),type="l",ylim = c(-1,1))
#' lines(x,Dnorm(x,D = 2,mean,sd),col="red")

Dnorm <- function(x,D,mean,sd)
{
  if (round(D) != D || D < 0) {
    print("D must be a non-negative integer.")
    return(NaN)
  }

  if (D == 0) output <- pnorm(x,mean,sd)
  if (D >= 1) output <- 1/sd*(-1)^(D-1)*hermite((x -mean)/sd,D-1)*dnorm( (x -mean)/sd)

  return(output)
}

