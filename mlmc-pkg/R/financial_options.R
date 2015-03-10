#
# #####European Option####
# output<-multilevel_mc(4,0.01,1,10^4,a=function(t,s){0.05*s},b=function(t,s){0.2*s},1,p=f_europeanOption(tMax = 1),payoff_is_a_path_functional = FALSE,useParallel = FALSE,nCores = NA,0.01)
# logs <- log(output$sigma,base=4)
# Sigma<-data.frame(l=seq(1:5),Log_Mvariance=logs)
# ggplot(Sigma,aes(l,Log_Mvariance))+geom_line(colour="turquoise4")
#
# Mean<-data.frame(l=seq(1:5),Log_Mmean=log(abs(output$Y),base=4))
# ggplot(Mean,aes(l,Log_Mmean))+geom_line(colour="turquoise4")
#
# #####Asian Option####
# output.asia<-multilevel_mc(4,0.01,1,10^4,a=function(t,s){0.05*s},b=function(t,s){0.2*s},1,p=f_asianOption(),payoff_is_a_path_functional = TRUE,useParallel = FALSE,nCores = NA,0.01)
# logs.asia <- log(output.asia$sigma,base=4)
# Sigma.asia<-data.frame(l=seq(1:4),Log_Mvariance=logs.asia)
# ggplot(Sigma.asia,aes(l,Log_Mvariance))+geom_line(colour="turquoise4")
#
# Mean.asia<-data.frame(l=seq(1:4),Log_Mmean=log(abs(output.asia$Y),base=4))
# ggplot(Mean.asia,aes(l,Log_Mmean))+geom_line(colour="turquoise4")
#
# #####Loopback Option####
# output.loopback<-multilevel_mc(4,0.01,1,10^4,a=function(t,s){0.05*s},b=function(t,s){0.2*s},1,p=f_loopbackOption(),payoff_is_a_path_functional = TRUE,useParallel = FALSE,nCores = NA,0.01)
# logs.loopback <- log(output.loopback$sigma,base=4)
# Sigma.loopback<-data.frame(l=seq(1:4),Log_Mvariance=logs.loopback)
# ggplot(Sigma.loopback,aes(l,Log_Mvariance))+geom_line(colour="turquoise4")
#
# Mean.loopback<-data.frame(l=seq(1:4),Log_Mmean=log(abs(output.loopback$Y),base=4))
# ggplot(Mean.loopback,aes(l,Log_Mmean))+geom_line(colour="turquoise4")
#
# #####Digital Option####
# output.digital<-multilevel_mc(4,0.01,1,10^4,a=function(t,s){0.05*s},b=function(t,s){0.2*s},1,p=f_digitalOption(),payoff_is_a_path_functional = FALSE,useParallel = FALSE,nCores = NA,0.01)
# logs.digital <- log(output.digital$sigma,base=4)
# Sigma.digital<-data.frame(l=seq(1:4),Log_Mvariance=logs.digital)
# ggplot(Sigma.digital,aes(l,Log_Mvariance))+geom_line(colour="turquoise4")
#
# Mean.digital<-data.frame(l=seq(1:4),Log_Mmean=log(abs(output.digital$Y),base=4))
# ggplot(Mean.digital,aes(l,Log_Mmean))+geom_line(colour="turquoise4")
