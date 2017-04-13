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
Dlnorm <- function(x,D,meanlog = 0,sdlog = 1)
{
  if (round(D) != D || D < 0) {
    print("D must be a non-negative integer.")
    return(NaN)
  }


  if (D == 0) output <- plnorm(x,meanlog,sdlog)
  if (D == 1) output <- dlnorm(x,meanlog,sdlog)
  if (D > 1){
    output <- rep(0,length(x))
    supp <- x > 0
    xs <- x[supp]
    Xs <- (log(xs) - meanlog)/sdlog + sdlog


    if (D == 2) output[supp] <- -1/(sdlog*xs)^2*dnorm((log(xs)-meanlog)/sdlog)*Xs
    if (D == 3) output[supp] <- 1/(xs*sdlog)^3*dnorm((log(xs)-meanlog)/sdlog)*(Xs^2 + sdlog*Xs - 1)
    if (D > 3)  {
      print("Not available for this package version")
      output[supp]  <- NA
    }
  }
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
partialExpectationlnorm <- function(x, meanlog = 0,sdlog = 1, lower = TRUE)
{
  mean <- exp(meanlog + sdlog^2/2)
  if (lower)  output <- mean*pnorm( (log(x) -meanlog)/sdlog - sdlog)
  if (!lower) output <- mean*(1 - pnorm( (log(x) -meanlog)/sdlog - sdlog))
  return(output)
}
