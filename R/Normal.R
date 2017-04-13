#' Derivative w.r.t. x of the standard Normal distribution function
#' @param x vector of quantiles
#' @param D an integer between 0 and 3. The derivatives are not available
#'  in this package when D > 3
#'
#' @return Dnorm gives the  distribution function when D = 0, the density
#'  when D  = 1, the derivative of the density when D = 2, and so on
#' @export
#'
#' @examples
#' Dnorm(1,1)
Dlnorm <- function(x,D)
{
  if (round(D) != D || D < 0) {
    print("D must be a non-negative integer.")
    return(NaN)
  }

  if (D == 0) output <- pnorm(x)
  if (D >= 1) output <- (-1)^(D-1)*hermite(x,D-1)*dnorm(x)

  return(output)
}

