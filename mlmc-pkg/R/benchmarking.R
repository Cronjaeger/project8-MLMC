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
  P_euler <- vector(0.0,length = trials)
  P_mlmc <- vector(0.0,length = trials)
  runtime_euler <- vector(0.0,length = trials)
  runtime_mlmc <- vector(0.0,length = trials)

  modell <- SDE_geomBM(r = 0.05,s = 0.2,X0 = 1)

  for(i in 1:trials){
    out_euler <- euler_maruyama(h = h_euler,
                                t_max = 1.0,
                                drift = modell$r,
                                dispersion =modell$s,
                                X0 = modell$X0)
    out_multilevel_mc <- multilevel_mc()
  }

}
