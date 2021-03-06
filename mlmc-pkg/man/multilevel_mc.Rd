% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/multilevel_mc.R
\name{multilevel_mc}
\alias{multilevel_mc}
\title{Implements the Multilevel Monte Carlo Path Simulation described
in the paper by Giles 2008}
\usage{
multilevel_mc(M, h0, Tmax, N_L = 10^4, a, b, S0,
  payoffFunction = function(S) S, payoff_is_a_path_functional = FALSE,
  useParallel = FALSE, nCores = NA, epsilon = 0.001)
}
\arguments{
\item{M}{The factor of refinement between two levels.
Should be a positive integer}

\item{h0}{The step size used at level 0. At level L, the step size
is given by h0*M^-L}

\item{Tmax}{the time-horizon of simulations. The underlying random
variable is simulated on the time-interval [0,Tmax]}

\item{N_L}{the number of paths to be sampled.}

\item{a}{the coefficient of drift in the diffusion
Should have the form a <- function(t,s) ...}

\item{b}{the coefficient of dispersion in the diffusion
Should have the form b <- function(t,s) ...}

\item{S0}{the initial condition of the diffusion}

\item{payoffFunction}{the payoffFunction to be used
#' Should have the form payoffFunction <- function(s_tMax) ...}

\item{payoff_is_a_path_functional}{is a boolean indicatin weather
\code{payoff-function} is a functrion of S[t_max] or S[0:t_max].
In the former case, A lot of memmory is saved by not storing whole
paths, but only positions mooving as the simulation progresses.}

\item{useParallel}{a boolean value indicating if paths schoud be
generated and computed in parallel. Parallel compuattion depends
on the parallel-package.}

\item{nCores}{the number of cores to use when running in parallel.
if set to \code{NA}, the total number of availiable cores is used.}

\item{epsilon}{r.m.s. accuracy used to determine the optimal
number of paths.}
}
\value{
A list containing the optimal Y_hat, the variance of the
sequence, the optimal number of paths and the computing times of
each path.
}
\description{
calculates the value of Y_hat

Y_hat = N_L^(-1)sum_(i=1)^L (P_l-P_{l-1})

obtaining the optimal number of paths N_opt to be sampled of the SDE

dS_t = a(t,S_t) dt + b(t,S_t) dW_t

using the Euler-Maruyama method in the Rfunction: \code{euler_maruyama}
}
\examples{
output<-multilevel_mc(4,0.1,1,10^2,a=function(t,s){0.05*s},b=function(t,s){0.2*s},1,p=f_europeanOption(),payoff_is_a_path_functional = FALSE,useParallel = FALSE,nCores = NA,0.1)
}

