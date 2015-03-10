##perform benchmarking

require(microbenchmark)
require(ggplot2)
source("R/euler_maruyama.R")
source("R/examples.R")
source("R/multilevel_mc2.R")
set.seed(seed = 1742L)
Parallel <- function(L,M,N_L) euler_maruyama_multilevel(L = L, M = M, h0 = 1, Tmax = 1, N_L = N_L,a = function(t,s) 0.05*s, b = function(t,s) 0.2*s, S0 = 1,useParallel = TRUE)
Serial <- function(L,M,N_L) euler_maruyama_multilevel(L = L, M = M, h0 = 1, Tmax = 1, N_L = N_L,a = function(t,s) 0.05*s, b = function(t,s) 0.2*s, S0 = 1,useParallel = FALSE)
res <- microbenchmark(Serial(L = 1,M = 2,N_L = 10^4),
                      Parallel(L = 1,M = 2,N_L = 10^4),
                      Serial(L = 3,M = 4,N_L = 10^4),
                      Parallel(L = 3,M = 4,N_L = 10^4),
                      Serial(L = 5,M = 4,N_L = 10^2),
                      Parallel(L = 5,M = 4,N_L = 10^2),
                      times=100L)

print(res)
pdf(file="benchmarking_serial_vs_parallel_boxplot.pdf")
boxplot(res,log=FALSE)
dev.off()
autoplot(res)
ggsave(filename = "benchmarking_serial_vs_parallel_autoplot.pdf")
