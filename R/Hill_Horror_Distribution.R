####
#### "Hill Horror Plot" distribution
#### F(x) <- 1 - 1/(x log x)
####
f_inv_horror_dist <- function(d)
{
  # Define the interval in which to search for the root
  lower <- sqrt(d)
  upper <- max(exp(1),d)
  output <- uniroot( function(x,d)x*log(x) -d, d = d, interval = c(lower,upper))$root
  return(output)
}

# Using the inversion CDF method, we simulate
# a sample
rhorror <- function(n)
{
  U <- runif(n)
  d <- 1/(1-U)
  output <- sapply(d,f_inv_horror_dist )
  return(output)
}

# Check that the simulation is accurate
# X = sort(simulate_horror_dist(1e5))
# Fn = ecdf(X)
# plot(Fn,xlim = c(f_inv(1),max(X)),log = "x" )
# lines(X,1- 1/(X*log(X)),col = "red"   )
qhorror <- function(p)
{
  q <- rep(NA,length(p))
  q[p==1] <- Inf
  ind <- 0 <= p & p <1
  q[ind] <-  sapply(1/(1-p[ind]),f_inv_horror_dist)
  output <- q
  return(output)
}
phorror   = function(x)
{
  xmin <- qhorror(0)
  ind_supp <- which(x>= xmin)
  output <- rep(0,length(x))
  z <- x[ind_supp]
  output[ind_supp] <- 1- 1/(z*log(z))
  return(output)
}
dhorror <- function(x){
  xmin <- qhorror(0)
  ind_supp <- which(x>= xmin)
  output <- rep(0,length(x))
  z <- x[ind_supp]
  output[ind_supp] <- (log(z)+1)/(z*log(z))^2
  return(output)
}
ddhorror <- function(x) {
  xmin <- qhorror(0)
  ind_supp <- which(x>= xmin)
  output <- rep(0,length(x))
  z <- x[ind_supp]
  output[ind_supp] <- (log(z) - 2*(log(z)+1)^2)/(z*log(z))^3
  return(output)
}
d2dhorror <- function(x){
  xmin <- qhorror(0)
  ind_supp <- which(x>= xmin)
  output <- rep(0,length(x))
  z <- x[ind_supp]
  output[ind_supp]  <- (6*(log(z)+1)^3 -6*log(z) -7*log(z)^2 )/(z*log(z))^4
  return(output)
}
UpperTruncHorrorMoment <- function(a,order,seed, eps = 1e-4)
{
  xmin <- qhorror(0)
  A <- pmax(a,xmin)
  if (order >= 1) output <- Inf
  if (order < 1)
  {
    if (order == 0) output <-  1 - phorror(A)
    else {
      rate <- (1-order)*log(A)
      if (rate>=1)  output <- as.numeric(A^(order - 1)/log(A) + order*gammainc(rate,0)[2])
      if (rate< 1){
        if (!missing(seed)) set.seed <- seed
        X <- rexp(1e7,rate)
        sig <- sd(1/X*(X>=1))
        nsam <- (2*qnorm(0.975)*sig/(eps*rate))^2
        if (!missing(seed)) set.seed = seed
        X <- rexp(nsam,rate)
        output <-  A^(order - 1)/log(A) + order*mean(1/X*(X>=1))/rate
      }
    }
  }
  return(output)
}
