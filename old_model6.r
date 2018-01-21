
#
#   WARNING:  OLD SCRIPT, modello e selezione di covariate trovate senza standardizzare 
#
#
#  RNFL_average as response of several variables, 13 FIXED coefficients
#  attempt of covariate selection; a very simple model
#
# NB:  ZELLNER PRIOR, c=10 (per fixed coefficients)!
#      confrontato con model7: stesse covariate, ma intercetta e macular volume sono RANDOM!


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
        rep(1,length(Patient)),  #beta1 (intercept)
        Black,               #beta2
        White,               #beta3
        familiarity_yes,     #beta4
        Sex,                 #beta5
        POAG,                #beta6
        age65,               #beta7
        PSD,                 #beta8
        macular_volume,      #beta9
        Cup_area,            #beta10
        Rim_area,            #beta11
        MD,                  #beta12
        IOP                  #beta13
        )       



# Hyperparameters:
mu0=rep(0, dim(X)[2])
c=10
S0= c*solve(t(X)%*%X)  #diag(rep(1, dim(X)[2])) 


# NO RANDOM
data=list( y=RNFL_average, X=X,  npat= length(unique(Patient)), 
           mu0=mu0, S0=S0,  numerosity = numerosity, kk=kk)      # dati che passo a jags


#fixed coefficients initialization:
beta= rep(0,dim(X)[2])

inits = function() {list( beta=beta, sigma0=50,sigma1=50)} 


modelRegress=jags.model("data5_norandom.bug",data=data,inits=inits,n.adapt=1000,n.chains=1)
update(modelRegress,n.iter=19000)
variable.names=c("beta", "sigma0","sigma1","mu")
n.iter=50000 
thin=10 

library(coda)
library(plotrix)
outputRegress=coda.samples(model=modelRegress,variable.names=variable.names,n.iter=n.iter,thin=thin)


data.out=as.matrix(outputRegress)
data.out=data.frame(data.out)
attach(data.out)
n.chain=dim(data.out)[1] 
n.chain
#summary(data.out)
#head(data.out)

#save.image("../R_object/model_6.RData")

load("../R_object/model_6.RData")

######### 90% CI  FIXED COEFFICIENTS ########

library(ggplot2)

beta=data.out[,grep("beta", names(data.out), fixed=TRUE)]
names(beta)
colnames(beta)=names(data.frame(X))
colnames(beta)[1]="Intercept"
names(beta)
tmp=melt(beta)
head(tmp)
colnames(tmp) = c("Parameter", "value")
CI.beta = apply(beta, 2, quantile, c(0.05, 0.95)) 
CI.beta
p=ggs_caterpillar(tmp, thick_ci = c(0.05, 0.95), thin_ci = c(0.025, 0.975))
p  + geom_vline(xintercept=0, col="orange") 

# intercept very large, some covariates seems to be uninfluent

#ggsave('something.png', width = 8, height = 4)




############# RESIDUI BAYESIANI #############

mu=data.out[,grep("mu", names(data.out), fixed=TRUE)]
pred.mean=apply(mu,2,mean)
pred.sd=apply(mu,2,sd)

bres= (pred.mean- RNFL_average)/pred.sd   # residui bayesiani
out2 = (abs(bres) > 2) #as a reference value we take 2, (or 1.8)

length(which(out2==TRUE))

# Predictive goodness-of-fit: SUM (or MEAN) of the predictive Bayesian residuals
# to compare different models 
# The "best" model is the one with the smallest value for this index

sum(bres^2)/1046



