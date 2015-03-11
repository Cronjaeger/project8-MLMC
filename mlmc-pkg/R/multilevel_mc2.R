#' Implements the Multilevel Monte Carlo Path Simulation described
#' in the paper by Giles 2008
#'
#' @description obtaines the value of Y_hat
#'
#' Y_hat = N_L^(-1)sum_(i=1)^L (P_l-P_{l-1})
#'
#' obtaining the optimal number of paths N_opt to be sampled of the SDE
#'
#' dS_t = a(t,S_t) dt + b(t,S_t) dW_t
#'
#' using the Euler-Maruyama method in the Rfunction: \code{euler_maruyama}
#'
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
#' @param epsilon r.m.s. accuracy used to determine the optimal
#' number of paths.
#'
#' @return A list containing the optimal Y_hat, the variance of the
#' sequence, the optimal number of paths and the computing times of
#' each path.
#'
#' @example R/financial_options.R
#'
#' @export multilevel_mc

multilevel_mc<-function(
  M,
  h0,
  Tmax,
  N_L = 10^4,
  a,
  b,
  S0,
  #payoffFunction = f_europeanOption(),
  payoffFunction = function(S) S,
  payoff_is_a_path_functional = FALSE,
  useParallel = FALSE,
  nCores = NA,
  epsilon = 0.001){

  L<-0 #It's meant to be ZERO
  P_0<-0
  Y_hat_sequence <- vector()
  N_opt_sequence <- vector()
  var_sequence <- vector()
  times <- vector()
  stopMLMC<-FALSE
  Sum<-0
  while (L<2 || !stopMLMC) {
    t_start <- Sys.time()
#    print(paste("Running loop for L =",L))
    diff<-euler_maruyama_multilevel(L,M,h0,Tmax,N_L,a,b,S0,payoffFunction,payoff_is_a_path_functional,useParallel,nCores)
    #Step 2
    sigma<-var(diff)
    var_sequence <- c(var_sequence,sigma)
    Y_hat<-mean(diff)
    Y_hat_sequence <- c(Y_hat_sequence,Y_hat)
    #Step 3
    Sum <- Sum+sqrt(sigma/(Tmax/(M^L)))

    N_opt_sequence_old <- c(N_opt_sequence,N_L)
    N_opt_sequence <- sapply(X=0:L,FUN = function(l) ceiling(2*epsilon^(-2) * sqrt(var_sequence[l+1]*Tmax/(M^l)) * Sum ))

#     N_opt<-ceiling(2*epsilon^(-2)*sqrt(sigma*Tmax/(M^L))*Sum)
#     N_opt_sequence <- c(N_opt_sequence,N_opt)

#    print(paste("N_opt =",N_opt))
    #Step 4
#     print(N_opt_sequence_old)
#     print(N_opt_sequence)
    for(l in 0:L){
      if(N_opt_sequence_old[l+1] < N_opt_sequence[l+1]){
        alpha <- N_opt_sequence_old[l+1]/N_opt_sequence[l+1]
        diffNew <- euler_maruyama_multilevel(L = l,
                                             M = M,
                                             h0 = Tmax,
                                             Tmax = Tmax,
                                             N_L = N_opt_sequence[l+1] - N_opt_sequence_old[l+1],
                                             a = a,
                                             b = b,
                                             S0 = S0,
                                             payoffFunction = payoffFunction)
        #diff<- c(diff,euler_maruyama_multilevel(L,M,h0,Tmax,N_opt - N_L,a,b,S0,payoffFunction))
        Y_hat_sequence[l+1] <- alpha*Y_hat_sequence[l+1] + (1 - alpha) * mean(diffNew)
      }
    }
    L <- L+1
  #Step 5
    if(L > 1){
    stopMLMC<-if((max((1/M)*abs(Y_hat_old),abs(Y_hat))<1/sqrt(2)*(M-1)*epsilon) || L > 10){TRUE}
            else{FALSE}
#   stopMLMC<-if(L > 4){TRUE}
#   else{FALSE}
    }

    Y_hat_old <- Y_hat
    #Y_hat_sequence <- c(Y_hat_sequence,Y_hat)
#    print(Y_hat_sequence)
    t_end <- Sys.time()
    times <- c(times,difftime(t_end,t_start,units = "secs"))
  }
  return(list( Y = Y_hat_sequence , sigma = var_sequence , N_opt = N_opt_sequence, t_diff = times))
}
