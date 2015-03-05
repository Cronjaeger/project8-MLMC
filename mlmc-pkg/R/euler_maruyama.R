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

euler_maruyama_multilevel <- function(L,M,h0,Tmax,NL = 10^4,a,b,S0)
  #implement stuff
  return(P_fine -P_coarse)
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
}

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
  minVal <- 0.9 #makimum for plotting
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
