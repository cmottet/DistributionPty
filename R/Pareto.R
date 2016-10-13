#' The Pareto Distribution
#'
#' Density, distribution function, distribution function derivatives, quantile function
#' and random generation for the Pareto distribution function with shape \code{shape} and
#' scale \code{scale}.
#' @rdname Pareto
#' @aliases rpareto
#' @aliases ppareto
#' @aliases qpareto
#' @aliases dpareto
#' @param x,q Vector of quantiles
#' @param p  vector of probabilities
#' @param n Number of observations. If length(n) > 1, the length is taken
#' to be the number required. Default value is 1
#' @param scale,shape Positive real values respectively defining the shape and scale
#' parameter of the Pareto distributon. Default value is 1 for both of them
#' @param d An non-negative integer giving the order of the derivative
#' @details The Pareto distribution has the following cumultative distribution function
#' \deqn{F(x) = 1 - (x/scale)^(-\alpha)} for all x > scale
#' @return \code{dpareto} gives the density, \code{ppareto} gives the distribution function,
#' \code{Dpareto} gives the distribution function derivative, \code{qpareto} gives the quantile
#' function, and \code{rpareto} generates random observations.
#'
#' @export
#'
#' @examples
#' # Simulate a pareto sample
#' n <- 1e2 ; scale <- 1 ; shape <- 1
#' set.seed(100) ; X <- sort(rpareto(n, scale , shape))
#'
#' # Compare the ECDF and the distribution function
#' plot(X, ecdf(X)(X))
#' lines(X,ppareto(X,shape,scale) , col = "red")
#'
#' # Compare kernel density and density function
#' with(density(X, from = min(X)), plot(x, y))
#' lines(X,dpareto(X,shape,scale), col = "red")
#'
#' # Visualize the distribution derivatives for 0 <= d <= 3
#' D <- 0:3
#'
#' derivatives <- sapply(D, function(d) Dpareto(X,d,scale, shape))
#' dataPlot <- data.frame(
#'              x = rep(X,length(D)),
#'              y = c(derivatives),
#'              D = rep(paste0("d = ",D ), each = length(X))
#'              )
#'
#' library(ggplot2)
#' ggplot(dataPlot, aes(x = x, y = y)) +
#' geom_line() +
#' facet_grid(D~., scales = "free_y")
rpareto  <- function(n = 1,scale=1,shape=1)
{
  if (scale <= 0 || shape <= 0 ) return("The scale and shape parameters must be positive numbers.")

  U <- runif(n)
  output <- qpareto(U, scale, shape)
  return(output)
}

#' @export
#' @rdname Pareto
ppareto  <- function(x,scale=1,shape=1)
{
  if (scale <= 0 || shape <= 0 ) {
    print("The scale and shape parameters must be positive numbers.")
    return(rep(NA, length(x)))
  }

  p <- rep(0,length(x))
  supp <- scale <= x
  xs <- x[supp]
  p[supp]  <- 1 - (xs/scale)^(-shape)

  output <- p

  return(output)
}

#' @export
#' @rdname Pareto
dpareto  <- function(x,scale=1,shape=1)
{
  if (scale <= 0 || shape <= 0 ) {
    print("The scale and shape parameters must be positive numbers.")
    return(rep(NA, length(x)))
  }

  density <- rep(0,length(x))

  supp <- scale <= x
  xs <- x[supp]
  density[supp] <- shape/scale*(xs/scale)^(-shape-1)

  output <- density

  return(output)
}


#' @export
#' @rdname Pareto
Dpareto <- function(x,d,scale=1,shape=1)
{
  if (round(d) != d || d < 0) return("d must be a non-negative integer.")
  if (scale <= 0 || shape <= 0 ) return("The scale and shape parameters must be positive numbers.")

  derivative <- rep(0,length(x))
  supp <- scale <= x
  xs <- x[supp]
  if (d == 0) derivative[supp] <- ppareto(xs,scale,shape)
  if (d != 0) derivative[supp] <- -(-1/scale)^(d)*(xs/scale)^(-shape-d)*prod(shape + 0:(d-1))

  output <- derivative
  return(output)
}

#' @rdname Pareto
#' @export
qpareto <- function(p,scale=1,shape=1)
{
  q = rep(NA,length(p))

  q[p==1] <- Inf
  q[p==0] <- scale
  q[0<p & p<1] <- scale*(1-p[0<p & p<1])^(-1/shape)

  output <- q

  return(output)
}


#' Calculate E[X^mI(x =>a)] when X follows a pareto distribution
#'
#' @param a A real value
#' @param m A vector of real numbers
#' @inheritParams ppareto
#' @export
#' @examples
#' UpperTruncMomPareto(scale,1,scale=1,shape = 2)
UpperTruncMomPareto <- function(a,m,scale=1,shape = 1)
{
  uppMoments <- rep(NA,length(m))

  # Case when the order is larger or equal than the shape
  uppMoments[(m >= shape)] <- Inf

  # Case when the order is strictly smaller than the shape
  ind <- which(m < shape)
  uppMoments[ind]  <- shape*scale^shape/(shape-m[ind])*max(a,scale)^(m[ind]-shape)

  output <- uppMoments
  return(output)
}

