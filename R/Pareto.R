#' Title
#'
#' @param x
#' @param scale
#'
#' @export
#'
#' @examples
#' # Simulate a pareto sample
#' n <- 1e3 ; scale = 1/2 ; shape = 3
#' set.seed(100) ; X <- sort(rpareto(n, scale, shape))
#'
#' # Compare the ECDF with the theoritical CDF
#' plot(X,ecdf(X)(X), type = "l")
#' lines(X,ppareto(X, scale, shape), col = "red")
rpareto  <- function(n,scale=1,shape=1)
{
  U <- runif(n)
  output <- qpareto(U, scale, shape)
  return(output)
}
#' Title
#'
#' @param shape
#'
#' @return
#' @export
#'
#' @examples
#'
#' x <- seq(-3,15, length = 100)
#' plot(x,ppareto(x), type ="l")
ppareto  <- function(x,scale=1,shape=1)
{
  p <- rep(NA,length(x))
  p[x < scale] <- 0

  xs <- x[x >= scale]
  p[x >= scale]  <- 1 - (xs/scale)^(-shape)

  output <- p

  return(output)
}

#' Title
#'
#' @param x
#' @param scale
#' @param shape
#'
#' @return
#' @export
#'
#' @examples
#' d <- 1:4
#' x <- seq(-3,15, length = 100)
#'
#' # Compute the firts 4 derivatives
#' dataPlot <- NULL
#' for (D in d)
#' dataPlot <- rbind(dataPlot, data.frame(d = as.factor(D), x = x, value = (Dpareto(x,D))))
#'
#' # Plot the derivatives
#' library(ggplot2)
#' ggplot(dataPlot, aes(x = x, y = value, color = d)) + geom_line()

Dpareto  = function(x,d,scale=1,shape=1)
{
  derivative <- rep(NA,length(x))
  derivative[x < scale] <- 0
  xs <- x[x >= scale]
  derivative[x >= scale] <- -(-1/scale)^(d)*(xs/scale)^(-shape-d)*prod(shape + 0:(d-1))

  output <- derivative
  return(output)
}


#' Title
#'
#' @param p
#' @param scale
#' @param shape
#'
#' @return
#' @export
#'
#' @examples
#' x <- seq(1,15, length = 100)
#' p <- ppareto(x)
#' q <- qpareto(p)
#' hist(x-q)

qpareto <- function(p,scale=1,shape=1)
{
  q = rep(NA,length(p))

  q[p==1] <- Inf
  q[p==0] <- scale
  q[0<p & p<1] <- scale*(1-p[0<p & p<1])^(-1/shape)

  output <- q

  return(output)
}


#' Title
#'
#' @param a
#' @param m
#' @param scale
#' @param shape
#'
#' @return
#' @export
#'
#' @examples
UpperTruncMomPareto = function(a,m,scale=1,shape = 1) #E[X^mI(x>a)]
{
  uppMoments = rep(NA,length(m))

  # Case when the order is larger or equal than the shape
  uppMoments[(m >= shape)] = Inf

  # Case when the order is strictly smaller than the shape
  ind = which(m < shape)
  uppMoments[ind]  = shape*scale^shape/(shape-m[ind])*max(a,scale)^(m[ind]-shape)

  output = uppMoments
  return(output)
}

