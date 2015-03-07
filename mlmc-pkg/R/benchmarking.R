##perform benchmarking

require(microbenchmark)
require(ggplot2)
source("R//euler_maruyama.R")

Parallel <- function(N) euler_maruyama_multilevel(L = 3, M = 4, h0 = 0.5, Tmax = 1, N_L = N,a = function(t,s) 0.05*s, b = function(t,s) 0.2*s, S0 = 1,useParallel = TRUE)
Serial <- function(N) euler_maruyama_multilevel(L = 3, M = 4, h0 = 0.5, Tmax = 1, N_L = N,a = function(t,s) 0.05*s, b = function(t,s) 0.2*s, S0 = 1,useParallel = FALSE)

res <- microbenchmark(Serial(10000),Parallel(10000),times=100L)

print(res)
pdf(file="benchmarking.pdf")
boxplot(res)
dev.off()
autoplot(res)

