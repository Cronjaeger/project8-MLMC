##perform benchmarking
testParallelization <- function(n_samples){
  require(microbenchmark)
  require(ggplot2)
  #require(mlmcpkg)
  source("R/euler_maruyama.R")
  source("R/examples.R")
  source("R/multilevel_mc2.R")
  set.seed(seed = 1742L)
  Parallel <- function(L,M,N_L) euler_maruyama_multilevel(L = L, M = M, h0 = 1, Tmax = 1, N_L = N_L,a = function(t,s) 0.05*s, b = function(t,s) 0.2*s, S0 = 1,useParallel = TRUE)
  Serial <- function(L,M,N_L) euler_maruyama_multilevel(L = L, M = M, h0 = 1, Tmax = 1, N_L = N_L,a = function(t,s) 0.05*s, b = function(t,s) 0.2*s, S0 = 1,useParallel = FALSE)
  res <- microbenchmark(Serial(L = 1,M = 4,N_L = 10^4),
                        Parallel(L = 1,M = 4,N_L = 10^4),
                        Serial(L = 3,M = 4,N_L = 10^4),
                        Parallel(L = 3,M = 4,N_L = 10^4),
                        Serial(L = 5,M = 4,N_L = 10^3),
                        Parallel(L = 5,M = 4,N_L = 10^3),
                        Serial(L = 6,M = 4,N_L = 94),
                        Parallel(L = 6,M = 4,N_L = 94),
                        times=n_samples)

  print(res)
  pdf(file="benchmarking_serial_vs_parallel_boxplot.pdf")
  #boxplot(res,log=FALSE)
  #vertical boxes
  # par(mar = c(14, 4, 1, 1)+ 0.1)
  # boxplot(res,horizontal = FALSE,las=2,xlab = "")

  #horizontal boxes
  par(mar = c(4, 16, 1, 1)+ 0.1)
  boxplot(res,horizontal = TRUE,log=FALSE,las=1,ylab = "",xlab = "")
  dev.off()
  autoplot(res)
  ggsave(filename = "benchmarking_serial_vs_parallel_autoplot.pdf")
}

testAgainst_euler <- function(epsilon, h_euler, N_euler,trials = 10L){
  P_euler <- vector(mode = "numeric",length = trials)
  P_mlmc <- vector(mode = "numeric",length = trials)
  P_list_euler <- vector(mode = "numeric",length = N_euler)
  runtime_euler <- vector(mode = "numeric",length = trials)
  runtime_mlmc <- vector(mode = "numeric",length = trials)

  modell <- SDE_geomBM(r = 0.05,s = 0.2,X0 = 1)
  f <- f_europeanOption(r = modell$X0, X0 = modell$X0, tMax = 1)

  for(i in 1:trials){
    t1 <- Sys.time()
    for(j in 1:N_euler){
      out_euler <- euler_maruyama(h = h_euler,
                                  t_max = 1.0,
                                  drift = modell$mu,
                                  dispersion =modell$sigma,
                                  X0 = modell$X0)
      P_list_euler[j] <- f(out_euler$X[length(out_euler$X)])
    }
    P_euler[i] <- mean(P_list_euler)
    t2 <- Sys.time()
    out_multilevel_mc <- multilevel_mc(M = 4,
                                       h0 = 1,
                                       Tmax = 1,
                                       N_L = 10^4,
                                       a = modell$mu,
                                       b = modell$sigma,
                                       S0 = modell$X0,
                                       payoffFunction = f,
                                       epsilon = epsilon,
                                       payoff_is_a_path_functional = FALSE,
                                       useParallel = FALSE)
    P_mlmc[i] <- sum(out_multilevel_mc$Y)
    t3 <- Sys.time()
    runtime_euler[i] <- difftime(t2,t1,units = "secs")
    runtime_mlmc[i] <- difftime(t3,t1,units = "secs")
  }
  hist(runtime_euler)
  hist(runtime_mlmc)
  return(list(P_euler = P_euler,P_mlmc = P_mlmc,runtime_euler = runtime_euler,runtime_mlmc = runtime_mlmc))
}
