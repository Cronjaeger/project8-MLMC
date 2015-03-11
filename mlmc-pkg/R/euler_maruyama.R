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
#' Should have the form a <- function(t,s) ...
#' @param b the coefficient of dispersion in the diffusion
#' Should have the form b <- function(t,s) ...
#' @param S0 the initial condition of the diffusion
#' @param payoffFunction the payoffFunction to be used
#' #' Should have the form payoffFunction <- function(s_tMax) ...
#' @param payoff_is_a_path_functional is a boolean indicatin weather
#' \code{payoff-function} is a functrion of S[t_max] or S[0:t_max].
#' In the former case, A lot of memmory is saved by not storing whole
#' paths, but only positions mooving as the simulation progresses.
#' @param useParallel a boolean value indicating if paths schoud be
#' generated and computed in parallel. Parallel compuattion depends
#' on the parallel-package.
#' @param nCores the number of cores to use when running in parallel.
#' if set to \code{NA}, the total number of availiable cores is used.
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
    payoffFunction = function(S) S,
    payoff_is_a_path_functional = FALSE,
    useParallel = FALSE,
    nCores = NA)
  {

  if(useParallel) require(parallel)

  #infer implicitly specified parameters
  h_fine <- h0 * M^(-L)
  h_coarse <- h0 * M^(-L + 1)
  n_steps_fine <- round(Tmax/(h_fine))
  n_steps_coarse <- round(Tmax/(h_coarse))

  #print(paste(n_steps_fine,n_steps_coarse,n_steps_fine/n_steps_coarse, M))

  if(!payoff_is_a_path_functional){
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
      return(payoffFunction(S))
    }
    reduction_step_fine <- function(increments){
      dim(increments) <- n_steps_fine
      t <- 0
      S <- S0
      for(i in 1:n_steps_fine){
        t <- i*h_fine
        S <- S + a(t,S)*h_fine + b(t,S)*increments[i]
      }
      return(payoffFunction(S))
    }
  } #end if(!payoff_is_a_path_functional)

  ## we do the above in a more general way
  else{
    #define timesteps
    timesteps_fine = h_fine*0:n_steps_fine
    timesteps_coarse = h_coarse*0:n_steps_coarse

    #define reduction-steps
    reduction_step_coarse <- function(increments){
      #axis 2 of increments indexes 1:n_steps_coarse
      #axis 1 of increments indexes 1:M
      dW <- apply(increments,MARGIN = 2,FUN = sum)
      t <- 0
      S <- c(S0,rep(0.0,n_steps_coarse))
      for(i in 1:n_steps_coarse){
        t <- i*h_coarse
        S[i+1] <- S[i] + a(t,S[i])*h_coarse + b(t,S[i])*dW[i]
      }
      return(payoffFunction(timesteps_coarse,S))
    }
    reduction_step_fine <- function(increments){
      dim(increments) <- n_steps_fine
      t <- 0
      S <- c(S0,rep(0.0,n_steps_fine))
      for(i in 1:n_steps_fine){
        t <- i*h_fine
        S[i+1] <- S[i] + a(t,S[i])*h_fine + b(t,S[i])*increments[i]
      }
      return(payoffFunction(timesteps_fine,S))
    }
  }

  #apply reduction-step to "N_L"-axis.
  if(!useParallel){
    #precompute increments
    if(L > 0){
      increments <- rnorm(N_L * M * max(1,n_steps_coarse) , sd = sqrt(h_fine))
      dim(increments) <- c(N_L,M,max(n_steps_coarse,1))
    } else {
      increments <- rnorm(N_L * n_steps_fine , sd = sqrt(h_fine))
      dim(increments) <- c(N_L,n_steps_fine)
    }
    if(L>0) P_coarse <- apply(X = increments, MARGIN = 1, FUN = reduction_step_coarse)
    else P_coarse = 0
    P_fine <- apply(X = increments, MARGIN = 1, FUN = reduction_step_fine)
  } else {
    if(is.na(nCores)) nCores <- detectCores()

    #precompute increments this time as a list of matrices...
    incrementGenerator <- function(x){
      increments <- rnorm(max(1,n_steps_coarse) * M, sd = sqrt(h_fine))
      dim(increments) <- c(M,max(n_steps_coarse,1))
      return(increments)
    }
    increments <- mclapply(1:N_L,FUN = incrementGenerator,
                           mc.cores = nCores,
                           mc.preschedule = TRUE)
    if(L>0) P_coarse <- mcmapply(FUN = reduction_step_coarse,
                                 increments,
                                 mc.preschedule = TRUE,
                                 mc.cores = nCores)
    else P_coarse = 0

    P_fine <- mcmapply(FUN = reduction_step_fine,
                       increments,
                       mc.preschedule = TRUE,
                       mc.cores = nCores)
  }

  return(P_fine - P_coarse)
}
