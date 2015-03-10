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
  epsilon = 0.001,
  payoff_is_a_path_functional = FALSE,
  useParallel = FALSE){

  L<-0 #It's meant to be ZERO
  P_0<-0
  Y_hat_sequence <- vector()
  N_opt_sequence <- vector()
  var_sequence <- vector()
  times <- vector()
  stopMLMC<-FALSE
  while (L<2 || !stopMLMC) {
    t_start <- Sys.time()
#    print(paste("Running loop for L =",L))
    diff<-euler_maruyama_multilevel(L,M,h0,Tmax,N_L,a,b,S0,payoffFunction,payoff_is_a_path_functional,useParallel)
    #Step 2
    sigma<-var(diff)
    var_sequence <- c(var_sequence,sigma)
    Y_hat<-mean(diff)
    sum<-0
    for (i in 0:L) {
      sum <- sum+sqrt(sigma/(Tmax/(M^i)))
    }
    #Step 3
    N_opt<-ceiling(2*epsilon^(-2)*sqrt(sigma*Tmax/(M^L))*sum)
    N_opt_sequence <- c(N_opt_sequence,N_opt)
#    print(paste("N_opt =",N_opt))
    #Step 4
    if(N_L<N_opt){
      diff<- c(diff,euler_maruyama_multilevel(L,M,h0,Tmax,N_opt - N_L,a,b,S0,payoffFunction))
      Y_hat<-mean(diff)
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
    Y_hat_sequence <- c(Y_hat_sequence,Y_hat)
#    print(Y_hat_sequence)
    t_end <- Sys.time()
    times <- c(times,difftime(t_end,t_start,units = "secs"))
  }
  return(list( Y = Y_hat_sequence , sigma = var_sequence , N_opt = N_opt_sequence, t_diff = times))
}
