#' Derivative w.r.t. x of the log-normal distribution function
#' @param x vector of quantiles
#' @param meanlog,sdlog mean and standard deviation of the distribution
#'  on the log scale with default values of 0 and 1 respectively.
#' @param D an integer between 0 and 3. The derivatives are not available
#'  in this package when D > 3
#'
#' @return Dlnorm gives the  distribution function when D = 0, the density
#'  when D  = 1, the derivative of the density when D = 2, and so on
#' @export
#'
#' @examples
#' Dlnorm(1,1,0,2)
Dlnorm <- function(x,D,meanlog,sdlog)
{
  if (round(D) != D || D < 0) {
    print("D must be a non-negative integer.")
    return(NaN)
    }

  X <- (log(x) - meanlog)/sdlog + sdlog

  if (D == 0) output <- plnorm(x,meanlog,sdlog)
  if (D == 1) output <- dlnorm(x,meanlog,sdlog)
  if (D == 2) output <- -1/(x*sdlog)*dlnorm(x,meanlog,sdlog)*X
  if (D == 3) output <- 1/(x*sdlog)^2*dlnorm(x,meanlog,sdlog)*(X^2 + sdlog*X - 1)
  if (D > 3)  output  <- ("Not available in this package.")
  return(output)
}

#' Compute E[X I(X <= x)] when X follows a log-normal distribution
#'
#' @param x A positive vector
#' @inheritParams Dlnorm
#' @param lower Logical value determining whether the lower or the upper partial
#' expectation should be computed, the function should return E[X I(X <= x)] or
#' E[X I(X > x)].  Default is TRUE for lower expectation.
#'
#' @return
#' @export
#'
#' @examples
#' partialExpectationlnorm(1, 0,1, lower = TRUE)
partialExpectationlnorm <- function(x, meanlog,sdlog, lower = TRUE)
{
  mean <- exp(meanlog + sdlog^2/2)
  if (lower)  output <- mean*pnorm( (log(x) -meanlog)/sdlog,sdlog,1)
  if (!lower) output <- mean*(1 - pnorm( (log(x) -meanlog)/sdlog,sdlog,1))
  return(output)
}
