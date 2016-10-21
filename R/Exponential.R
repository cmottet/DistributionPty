#' Derivative w.r.t. x of the exponential distribution function
#' @param x vector of quantiles
#' @param rate of the distributio
#' @param D a non-negative integer
#'
#' @return Dexp gives the  distribution function when D = 0, the density
#'  when D  = 1, the derivative of the density when D = 2, and so on
#' @export
#' @details F(x) = 1 - exp(-rate x)
#' @examples
#' Dexp(1,1,0,2)
Dexp <- function(x,D, rate = 1)
{
  if (round(D) != D || D < 0) {
    print("D must be a non-negative integer.")
    return(NaN)
  }
  if (D == 0) output <- pexp(x,rate)
  if (D >= 1 ) output <- ((-rate)^(D - 1))*dexp(x,rate)
  return(output)
}
