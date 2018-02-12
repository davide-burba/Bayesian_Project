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


#covariates with fixed coefficents
X=cbind(
  rep(1,length(Patient)),  # b1 (intercept)
  Black,               #beta2
  #Hispanic,            #beta3 
  familiarity_yes,     #beta4
  Sex,                 #beta5
  #POAG,                #beta6
  #Hypertension,        #beta7
  #HyperLipidemia,      #beta8
  #Cardiovascular_Dz,   #beta9
  age65,               #beta10
  #prostaglandin,       #beta11
  brimonidine,         #beta12
  #timolol,             #beta13
  #IOP,                     #beta14 
  MD,                      #beta15
  Vert_integrated_rim_area__vol_,     #beta17
  Horz_integrated_rim_width__area_   #beta18 
  #Rim_area               #beta19
  
)

#covariates with random coefficents
Z=cbind(
  rep(1,length(Patient)),  # b1 (intercept)
  #macular_volume,          # b2 
  visit2                   # time
)

# Hyperparameters:
mu0=rep(0, dim(X)[2])
c=50
S0=c*solve(t(X)%*%X) # Zellner prior 
b0=rep(0,dim(Z)[2])
R=diag(rep(1, dim(Z)[2])) 
#p=dim(Z)[2]
#p=100
p=30

data=list( y=RNFL_average, X=X, Z=Z, npat= length(unique(Patient)), nrow_b=dim(Z)[2], 
           mu0=mu0, S0=S0, b0=b0, R=R, p=p, numerosity = numerosity, kk=kk)      # dati che passo a jags

#fixed coefficients initialization:
beta= rep(0,dim(X)[2])

#random coefficients initialization (row_i= coefficients b_i):
b= matrix(0,nrow=dim(Z)[2], ncol=length(unique(Patient))) 
invD=p*solve(R)


inits = function() {list( beta=beta, b=b, invD=invD, sigma0=50,sigma1=50)} 


modelRegress=jags.model("data4.bug",data=data,inits=inits,n.adapt=1000,n.chains=1)
update(modelRegress,n.iter=19000)
variable.names=c("beta", "b", "sigma0","sigma1","mu")
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

save.image("../R_object/model_13.RData")

load("../R_object/model_13.RData")

#####################
## GODNESS OF MCMC ##
#####################

library(coda)        # pacchetto per analizzare catene
library(plotrix)     # per fare plot CIs



data=as.matrix(outputRegress) # trasformo il dataframe in matrice 
data=data.frame(data)
attach(data)
n.chain=dim(data)[1]   # lunghezza catena (final sample size)

names(data)
# let's check the traceplots of random b1
outputRegress_mcmc = as.mcmc(data[,grep("b.1", names(data), fixed=TRUE)])

quartz()
plot(outputRegress_mcmc)

# some autocorrelations plots
quartz()
acfplot(outputRegress_mcmc) 


# let's check the traceplots of random b2
outputRegress_mcmc = as.mcmc(data[,grep("b.2", names(data), fixed=TRUE)])

quartz()
plot(outputRegress_mcmc[,1:40])

# some autocorrelations plots
quartz()
acfplot(outputRegress_mcmc[,1:10]) 
acfplot(outputRegress_mcmc[,11:20])
acfplot(outputRegress_mcmc[,90:100])

# let's check the traceplots of random b3
outputRegress_mcmc = as.mcmc(data[,grep("b.3", names(data), fixed=TRUE)])

quartz()
plot(outputRegress_mcmc[,1:40])

# some autocorrelations plots
quartz()
acfplot(outputRegress_mcmc[,1:10]) 
acfplot(outputRegress_mcmc[,11:20])
acfplot(outputRegress_mcmc[,90:100])


# let's check the traceplots of fixed coefficients
outputRegress_mcmc = as.mcmc(data[,grep("beta", names(data), fixed=TRUE)])

quartz()
plot(outputRegress_mcmc)

# some autocorrelations plots
quartz()
acfplot(outputRegress_mcmc)




# let's check the traceplots of others
outputRegress_mcmc = as.mcmc(data[,grep("sigma", names(data), fixed=TRUE)])
quartz()
plot(outputRegress_mcmc)

#autocorrelation
quartz()
acfplot(outputRegress_mcmc)



##################### 90% CI FIXED COEFFICIENTS #####################

library(ggplot2)
library(ggmcmc)
library(reshape)

names(data.out[,grep("beta", names(data.out), fixed=TRUE)])
beta=data.out[,grep("beta", names(data.out), fixed=TRUE)]
names(beta)
colnames(beta)=names(data.frame(X))
names(beta)
tmp=melt(beta)
head(tmp)
colnames(tmp) = c("Parameter", "value")
CI.beta = apply(beta, 2, quantile, c(0.05, 0.95)) 
CI.beta
p=ggs_caterpillar(tmp, thick_ci = c(0.05, 0.95), thin_ci = c(0.025, 0.975))
p   + geom_vline(xintercept=0, col="orange") 

#ggsave('95CI_coefficients_random.png', width = 8, height = 4)




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

# 90% CI random coefficients
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



################################
####### GOODNESS OF FIT ########
################################


############# RESIDUI BAYESIANI #############
mu=data[,grep("mu", names(data), fixed=TRUE)]
pred.mean=apply(mu,2,mean)
pred.sd=apply(mu,2,sd)

bres= (pred.mean- RNFL_average)/pred.sd   # residui bayesiani
out2 = (abs(bres) > 2) #as a reference value we take 2, (or 1.8)

length(which(out2==TRUE))

# Predictive goodness-of-fit: SUM (or MEAN) of the predictive Bayesian residuals
# to compare different models 
# The "best" model is the one with the smallest value for this index
MEAN.RES.BAYES=mean(bres^2)
MEAN.RES.BAYES
MEDIAN.RES.BAYES=median(bres^2)
MEDIAN.RES.BAYES

plot(sort(bres^2)) 

## CROSS VALIDATION


#calcolo LPML
n=dim(mydata)[1]
n
cpo=seq(1,n)
names(data)
first.b1.pos=which(names(data)=="b.1.1.")
mcmc.std.dev=data[,grep("sigma0", names(data), fixed=TRUE)]
beta= data[,grep("beta", names(data), fixed=TRUE)]
tmp=X%*%t(beta)
dim(tmp)
npat=length(unique(Patient))
for (i in 1:npat){
  b1=data[,first.b1.pos+2*(i-1)]
  b2=data[,first.b1.pos+1+2*(i-1)]
  tmp_b=cbind(b1,b2)
  for(j in 1:numerosity[i])
  {
    pred=tmp[kk[i]+j,] + Z[kk[i]+j,]%*%t(tmp_b)
    cpo[kk[i]+j]=1/mean(1/dnorm(RNFL_average[kk[i]+j], pred , mcmc.std.dev))
  }
}
LPML=sum(log(cpo))
LPML
sum(cpo)
median(cpo)
mean(cpo)

# BIC & AIC (using the means of coefficients as summary of posterior)

names(data)
first.b1.pos=which(names(data)=="b.1.1.")
beta= data[,grep("beta", names(data), fixed=TRUE)]
tmp=X%*%t(beta)
npat=length(unique(Patient))
for (i in 1:npat){
  b1=data[,first.b1.pos+2*(i-1)]
  b2=data[,first.b1.pos+1+2*(i-1)]
  tmp_b=cbind(b1,b2)
  for(j in 1:numerosity[i])
  {
    pred[kk[i]+j]=tmp[kk[i]+j,] + Z[kk[i]+j,]%*%t(tmp_b)
  }
}

sigma.hat=mean(data[,grep("sigma0", names(data), fixed=TRUE)])
loglik=0
for (k in 1:n){
  loglik=loglik+dnorm(RNFL_average[k],mean=pred[k], sd =sigma.hat , log=TRUE)
}

# loglikelyhood:
loglik

# r is the number of parameters
r=length(names(data[,grep("b.", names(data), fixed=TRUE)])) + length(names(data[,grep("beta", names(data), fixed=TRUE)])) 

#AIC
AIC=2*loglik-2*r
AIC

#BIC
BIC=2*loglik-r*log(n)
BIC


###### OUTLIAR DETECTION

which(cpo<10^-5) 
# those patients are very unlikely to be seen under this model.. measurement error?
