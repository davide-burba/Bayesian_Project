#
#   WARNING:  OLD SCRIPT, modello e selezione di covariate trovate senza standardizzare 
#
#
#
#  RNFL_average as response of random intercept and time (random coefficients)
#  


rm(list=ls())
load("../R_object/Glaucoma_better_data.RData")
attach(mydata)


################JAGGAMENTO#######################

library(rjags)


# Compute useful values 
numerosity<-as.integer(table(Patient))
kk=rep(0,length(unique(Patient)))
kk[1]=0
for (i in 2:length(unique(Patient))){
  kk[i]=kk[i-1]+numerosity[i-1]
}

#covariates with random coefficents
Z=cbind(
  rep(1,length(Patient)),  #b1 (intercept)
  visit2                   #b2 
)

# Hyperparameters:



b0=rep(0,dim(Z)[2])
R=diag(rep(1, dim(Z)[2])) 
p=dim(Z)[2]

data=list( y=RNFL_average, Z=Z, npat= length(unique(Patient)), nrow_b=dim(Z)[2], 
             b0=b0, R=R, p=p, numerosity = numerosity, kk=kk)      # dati che passo a jags

#random coefficients initialization (row_i= coefficients b_i):
b= matrix(0,nrow=dim(Z)[2], ncol=length(unique(Patient))) 
invD=p*solve(R)


inits = function() {list(b=b, invD=invD, sigma0=50)} 


modelRegress=jags.model("data6_onlyrandom.bug",data=data,inits=inits,n.adapt=1000,n.chains=1)
update(modelRegress,n.iter=19000)
variable.names=c("b", "sigma0")
n.iter=50000 
thin=10 

library(coda)
library(plotrix)
outputRegress=coda.samples(model=modelRegress,variable.names=variable.names,n.iter=n.iter,thin=thin)


#outputRegress_mcmc <- as.mcmc(outputRegress)


data.out=as.matrix(outputRegress)
data.out=data.frame(data.out)
attach(data.out)
n.chain=dim(data.out)[1] 
n.chain
#summary(data.out)
#head(data.out)

#save.image("../R_object/model_8.RData")

load("../R_object/model_8.RData")

library(coda)
outputRegress_mcmc <- as.mcmc(outputRegress)

quartz()
plot(outputRegress_mcmc)




library(ggplot2)
library(ggmcmc)
library(reshape)



##################### 90% CI RANDOM COEFFICIENTS #####################

# 90% CI random coefficients: random intercept
tmp=data.out[,grep("b.1.", names(data.out), fixed=TRUE)]
names(tmp)
b=tmp
tmp=melt(b)
head(tmp)
colnames(tmp) = c("Parameter", "value")
CI.b = apply(b, 2, quantile, c(0.05, 0.95)) 
CI.b
p=ggs_caterpillar(tmp, thick_ci = c(0.05, 0.95), thin_ci = c(0.025, 0.975))
p  + geom_vline(xintercept=0, col="orange") 
#ggsave('95CI_coefficients_intercept.png', width = 5, height = 6)

# 90% CI random coefficients: macular volume
tmp=data.out[,grep("b.2.", names(data.out), fixed=TRUE)]
names(tmp)
b=tmp
tmp=melt(b)
head(tmp)
colnames(tmp) = c("Parameter", "value")
CI.b = apply(b, 2, quantile, c(0.05, 0.95)) 
CI.b
p=ggs_caterpillar(tmp, thick_ci = c(0.05, 0.95), thin_ci = c(0.025, 0.975))
p  + geom_vline(xintercept=0, col="orange") 

#ggsave('95CI_coefficients_macular.png', width = 5, height = 6)



################################################
####  Credible band for the regression line  ###
################################################

library(ggplot2)

hh=c(kk,1046)# patients i has rows from hh[i]+1 to hh[i+1]

mu=data.out[,grep("mu", names(data.out), fixed=TRUE)]

# plot i-esimi pazienti
i1=50
i2=51
i3=101

mu_pat_i=mu[,(hh[i1]+1):hh[i1+1]]
pred.mean=apply(mu_pat_i,2,mean)

lower=apply(mu_pat_i,2,quantile, 0.05)
upper=apply(mu_pat_i,2,quantile, 0.95)
plotdata <- data.frame(obs=RNFL_average[which(Patient==unique(Patient)[i1])], x=visit2[which(Patient==unique(Patient)[i1])], y=pred.mean, lower = lower, upper = upper)

p1= ggplot(plotdata) + geom_line(aes(y=y, x=x, colour = "sin"))+
  geom_ribbon(aes(ymin=lower, ymax=upper, x=x, fill = "band"), alpha = 0.3)+
  scale_colour_manual("",values="blue")+
  scale_fill_manual("",values="grey12")+
  geom_point(data=plotdata, aes(x=x, y=obs))+
  labs(  y="RNFL_average", x = "Time (Years)")+ theme(legend.position="none")
#p1


mu_pat_i=mu[,(hh[i2]+1):hh[i2+1]]
pred.mean=apply(mu_pat_i,2,mean)

lower=apply(mu_pat_i,2,quantile, 0.05)
upper=apply(mu_pat_i,2,quantile, 0.95)

plotdata <- data.frame(obs=RNFL_average[which(Patient==unique(Patient)[i2])], x=visit2[which(Patient==unique(Patient)[i2])], y=pred.mean, lower = lower, upper = upper)

p2= ggplot(plotdata) + geom_line(aes(y=y, x=x, colour = "sin"))+
  geom_ribbon(aes(ymin=lower, ymax=upper, x=x, fill = "band"), alpha = 0.3)+
  scale_colour_manual("",values="blue")+
  scale_fill_manual("",values="grey12")+
  geom_point(data=plotdata, aes(x=x, y=obs))+
  labs(  y="RNFL_average", x = "Time (Years)")+ theme(legend.position="none")
#p2


mu_pat_i=mu[,(hh[i3]+1):hh[i3+1]]
pred.mean=apply(mu_pat_i,2,mean)

lower=apply(mu_pat_i,2,quantile, 0.05)
upper=apply(mu_pat_i,2,quantile, 0.95)

plotdata <- data.frame(obs=RNFL_average[which(Patient==unique(Patient)[i3])], x=visit2[which(Patient==unique(Patient)[i3])], y=pred.mean, lower = lower, upper = upper)

p3= ggplot(plotdata) + geom_line(aes(y=y, x=x, colour = "sin"))+
  geom_ribbon(aes(ymin=lower, ymax=upper, x=x, fill = "band"), alpha = 0.3)+
  scale_colour_manual("",values="blue")+
  scale_fill_manual("",values="grey12")+
  geom_point(data=plotdata, aes(x=x, y=obs))+
  labs(  y="RNFL_average", x = "Time (Years)")+ theme(legend.position="none")
#p3


source ("multiplot.R")
multiplot(p1,p2,p3) 


######################## RESIDUI BAYESIANI ######################## 

mu=data.out[,grep("mu", names(data.out), fixed=TRUE)]
pred.mean=apply(mu,2,mean)
pred.sd=apply(mu,2,sd)

bres= (pred.mean- RNFL_average)/pred.sd   # residui bayesiani
out2 = (abs(bres) > 2) #as a reference value we take 2
# or 1.8

length(which(out2==TRUE))

# Predictive goodness-of-fit: SUM of the predictive Bayesian residuals
# to compare different models 
# The "best" model is the one with the smallest value for this index

sum(bres^2)/1046





