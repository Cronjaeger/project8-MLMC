# run tests

#' Implements the Euler-Maruyama method to solve the SDE
#' dX_t = drift(t,X_t) dt + dispersion(t,X_t) dW_t
#'
#' @param t_max the time horizon of the simulation.
#' the path of the process will be estimated on the interval
#' [0,floor(t_max/h)*h]
#' @param h the step size taken
#' @param drift the drift term
#'  (should accept arguments of the form (t,X_t))
#' @param dispersion the dispersion term (sigma)
#'  (should accept arguments of the form (t,X_t))
#' @param X0 initial conditions
euler_maruyama <- function(t_max,h,drift,dispersion,X0){
  N_steps <- floor(t_max / h)
  steps <- h * 0:N_steps
  dW <- rnorm(N_steps,sd = sqrt(h)) #precompute all increments
  X <- c(X0,rep(NA,N_steps))
  for(i in 1:N_steps){
    X[i+1] <- X[i] + h*drift(t,X[i]) + dW[i] * dispersion(t,X[i])
  }
  return(list(t=steps,X=X))
}

#' An auxiliary function for computing P_L - P_(L-1)
#'
#' @description simulates N_L independent paths of the SDE
#'
#' dS_t = a(t,S_t) dt + b(t,S_t) dW_t
#'
#' on two differnt grids -- one grid being finer than the
#' other by a factor of M. Each path on the coarse grid is generated
#' by summing noise from the finer grid.
#' The estimates of payoffFunction(S_Tmax) on the two different grids
#' is then compared and returned for each pair of paths
#'
#' @param L The level of the finest grid to be considered.
#' Should be a positive integer
#' @param M The factor of refinement between two levels.
#' Should be a positive integer
#' @param h0 The step size used at level 0. At level L, the step size
#' is given by h0*M^-L
#' @param Tmax the time-horizon of simulations. The underlying random
#' variable is simulated on the time-interval [0,Tmax]
#' @param N_L the number of paths to be sampled.
#' @param a the coefficient of drift in the diffusion
#' @param b the coefficient of dispersion in the diffusion
#' @param S0 the initial condition of the diffusion
#' @param payoffFunction the payoffFunction to be used
#'
#' @return A vector of length \code{N_L}. Each entry containing 1
#' estimate of the difference in discretisation-error between two
#' consecutive grids.
euler_maruyama_multilevel <- function(
    L,
    M,
    h0,
    Tmax,
    N_L = 10^4,
    a,
    b,
    S0,
    payoffFunction = function(S) S)
  {

  #infer implicitly specified parameters
  h_fine <- h0 * M^(-L)
  h_coarse <- h0 * M^(-L + 1)
  n_steps_fine <- Tmax/h_fine
  n_steps_coarse <- round(n_steps_fine/M)

  #precompute increments
  increments <- rnorm(N_L * n_steps_fine , sd = sqrt(h_fine))
  dim(increments) <- c(N_L,M,n_steps_coarse)

  #initialize result-vectors
  S_coarse <- rep(0.0,N_L)
  S_fine <- rep(0.0,N_L)

  #define reduction-steps
  reduction_step_coarse <- function(increments){
    #axis 2 of increments indexes 1:n_steps_coarse
    #axis 1 of increments indexes 1:M
    dW <- apply(increments,MARGIN = 2,FUN = sum)
    t <- 0
    S <- S0
    for(i in 1:n_steps_coarse){
      t <- i*h_coarse
      S <- S + a(t,S)*h_coarse + b(t,S)*dW[i]
    }
    return(S)
  }
  reduction_step_fine <- function(increments){
    dim(increments) <- n_steps_fine
    t <- 0
    S <- S0
    for(i in 1:n_steps_fine){
      t <- i*h_fine
      S <- S + a(t,S)*h_fine + b(t,S)*increments[i]
    }
    return(S)
  }

  #apply reduction-step to "N_L"-axis.
  S_coarse <- apply(X = increments, MARGIN = 1, FUN = reduction_step_coarse)
  S_fine <- apply(X = increments, MARGIN = 1, FUN = reduction_step_fine)

  P_fine <- payoffFunction(S_fine)
  P_coarse <- payoffFunction(S_coarse)

  return(P_fine - P_coarse)
}
  #function(t_max,h_fine,h_coarse,drift,dispersion,X0){
#   #infer implicit parameters
#   M <- round(h_fine/h_coarse)
#   N_steps_fine <- floor(t_max / h_fine)
#   N_steps_coarse <- round(N_steps_fine / M)
#   steps_fine <- h_fine * 0:N_steps_fine
#   steps_coarse <- h_coarse * 0:N_steps_coarse
#
#   dW_fine <- rnorm(N_steps_fine,sd = sqrt(h_fine))
#   #dW_coarse <- rep(0,N_steps_coarse)
#
#   X_fine <- c(X0,rep(NA,N_steps_fine))
#   X_coarse <- c(X0,rep(NA,N_steps))
#
#   #TODO: finiish this
#}

#' A script for testing the EM-mathod
#' @description Simulates a geometric brownian motion, and plots the results.
#' The number of paths sampled is \code{N_paths}; the number of steps taken
#' per path is \code{floor(t_max / h)}.
#'
#' @param N_paths the number of paths to be sampled
#' @param t_max the time-horizon of simulations
#' @param h the step size.
testEM <- function(N_paths = 20,t_max = 1, h = 0.01){
  #test on a geometric brownian motion
  r <- 0.05
  s <- 0.2
  X0 <- 1
  mu <- function(t,x) r*x
  sigma <- function(t,x) s*x

  #set parameters for EM algorithm
  t_max <- 1
  h <- 0.001

  #output is stored as a list of lists
  tests <- list()
  maxVal <- 1 #minimum of plotting
  minVal <- 0.9 #maximum for plotting
  for(i in 1:N_paths){
    result <- euler_maruyama(t_max = t_max,h=h,drift = mu,dispersion = sigma,X0 = X0)
    tests <- c(tests,list(result))
    maxVal <- max(maxVal, max(result$X))
    minVal <- min(minVal, min(result$X))
  }

  ##plot output
  par(bty = "l",family = "HersheySans",font=1)
  colours <- rainbow(N_paths,alpha = 0.7)
#  print(tests)
  for(i in 1:N_paths){
    plot(tests[[i]]$t,tests[[i]]$X,type="l",col = colours[i],xlab = "",ylab = "",ylim=c(minVal,maxVal),xaxs="i")
    par(new=TRUE)
  }
  par(new=FALSE)
#   em_test <- euler_maruyama(t_max = 1,h=0.01,drift = mu,dispersion = sigma,X0 = X0)
#   plot(em_test$t,em_test$X,type = "l")
}
