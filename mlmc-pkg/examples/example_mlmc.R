output<-multilevel_mc(4,0.1,1,10^2,a=function(t,s){0.05*s},b=function(t,s){0.2*s},1,p=f_europeanOption(),payoff_is_a_path_functional = FALSE,useParallel = FALSE,nCores = NA,0.1)
