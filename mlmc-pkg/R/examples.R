## Contains a range of different SDE-models and payoff-functions

SDE_geomBM <- function(r = 0.05,s = 0.2,X0=1.0) list(r=r,s=s,X0=X0,mu = function(t,x) r*x,sigma = function(t,x) s*x)

f_europeanOption <- function(r = 0.05, X0 = 1.0,tMax = 1.0){
  function(X_T) return(exp(-r*tMax)*max(0,X_T - X0))
}

f_asianOption <- function(r = 0.05 , X0 = 1.0){
  function(timesteps,X){
    N <- length(X)
    tMax <- timesteps[N]
    Xbar <- 0.5 * sum((X[-1] + X[-N]) * (timesteps[-1] - timesteps[-N]))
    return(exp(-r*tMax)*max(0,X-Xbar))
  }
}

f_loopbackOption <- function(r = 0.05, s = 0.2, X0 = 1.0){
  function(timestep, X){
    N <- length(X)
    tMax <- timesteps[N]
    minX <- min(X)
    h_l <- tMax/(N-1)
    betaStar <- 0.5826
    minX <- minX * (1 - betaStar*s*sqrt(h_l))
    return(exp(-r*tMax) * (X[N] - minX))
  }
}

f_digitalOption <- function(r = 0.05 , X0 = 1.0,tMax = 1.0){
  function(X_T){
    exp(-r*tMax)*as.numeric(X_T > X0)
  }
}
