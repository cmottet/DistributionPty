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

rlhorror <- function(n) log(rhorror(n))
###
### E[(X-a)_+^j]
###
UpperTrunclHorrorMoment <- function(a,order,seed, n = 1e6)
{
  xmin <- qlhorror(0)
  if (a < xmin){
    print("a must be larger than lower bound of the distribution support, i.e. 0.567")
    return(NaN)
  }
  output <- if (order == 0) 1 - plhorror(a) else{X = rlhorror(n); mean((X-a)^order*I(X >= a)) }
  return(output)
}
