###
### Horror Hill distribution with log-transform
###
# Check that the simulation is accurate
# X = sort(simulate_horror_dist(1e5))
# Fn = ecdf(X)
# plot(Fn,xlim = c(f_inv(1),max(X)),log = "x" )
# lines(X,1- 1/(X*log(X)),col = "red"   )
qlhorror <- function(p)
{
  output <- log(qhorror(p))
  return(output)
}

plhorror <- function(x)
{
  output <- phorror(exp(x))
  return(output)
}

dlhorror <- function(x){
  output <- exp(x)*dhorror(exp(x))
  return(output)
}

ddlhorror <- function(x) {
  output <-  exp(2*x)*ddhorror(exp(x))+ exp(x)*dhorror(exp(x))
  return(output)
}

d2dlhorror <- function(x){
  output <- exp(3*x)*d2dhorror(exp(x)) + 3*exp(2*x)*ddhorror(exp(x)) + exp(x)*dhorror(exp(x))
  return(output)
}

Dlhorror <- function(x,D){
  if ( !(D %in% 0:3)) return(rep(NA,length(x)))

  if (D == 0) output <- plhorror(x)
  if (D == 1) output <- dlhorror(x)
  if (D == 2) output <- ddlhorror(x)
  if (D == 3) output <- d2dlhorror(x)

  return(output)
}

UpperTrunclHorrorMoment <- function(a,order,seed, eps = 1e-4)
{
  xmin <- qlhorror(0)
  A <- pmax(a,xmin)
  if (order <0){ print("The order must be non-negative") ; return(NA)} else
  {
    if (order == 0) output <- 1 - plhorror(A)
    if (order == 1) {
      if (A > 1)  output <- exp(-A) + as.numeric(gammainc(A,0)[2])
      if (A <= 1) output <- exp(-A) + integrate(function(x)exp(-x)/x,lower = A,upper = Inf)$value
    }
    if ( !(order %in% c(0,1)) ) output <- (1+ 1/(order-1))*as.numeric(gammainc(A,order)[2]) -  A^(order -1)*exp(-A)/(order-1)
  }
  return(output)
}
