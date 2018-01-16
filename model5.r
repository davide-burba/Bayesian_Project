#
#  RNFL_average as response of several variables, 11 fixed and 6 random coefficients
#  1st attempt of a "complete" model; serve a fare SELEZIONE COVARIATE RANDOM (a spanne)
#  
# NOTE: 1) usata ZELLNER PRIOR, c=10  (per fixed coefficients)!
#       2) troppe covariate random!
#       3) intercetta è fixed


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
X=cbind(rep(1,length(Patient)),  #beta1 (intercept)
        Black,               #beta2
        White,               #beta3
        familiarity_yes,     #beta4
        Sex,                 #beta5
        POAG,                #beta6
        Hypertension,        #beta7
        age65,               #beta8
        brimonidine,         #beta9
        timolol,             #beta10
        htnmed              #beta11
        )       

#covariates with random coefficents
Z=cbind(
        IOP,                     #b1 
        MD,                      #b2 
        PSD,                     #b3 
        macular_volume,          #b4 
        Cup_area,                #b5
        Rim_area                 #b6
        )

# Hyperparameters:
mu0=rep(0, dim(X)[2])
c=10
S0= c*solve(t(X)%*%X)   #diag(rep(1, dim(X)[2])) 
b0=rep(0,dim(Z)[2])
R=diag(rep(1, dim(Z)[2])) 
p=dim(Z)[2]

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
variable.names=c("beta", "b", "sigma0","sigma1")
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

#save.image("../R_object/model_5.RData")

load("../R_object/model_5.RData")


################ 90 % CI FIXED COEFFICIENTS ################

library(reshape)
library(ggmcmc)


beta=data.out[,grep("beta", names(data.out), fixed=TRUE)]
names(beta)
colnames(beta)=names(data.frame(X))
colnames(beta)[1]="Intercept"
names(beta)
tmp=melt(beta)
head(tmp)
colnames(tmp) = c("Parameter", "value")
CI.beta = apply(beta, 2, quantile, c(0.05, 0.95))  # 90% 
CI.beta

p=ggs_caterpillar(tmp, thick_ci = c(0.05, 0.95), thin_ci = c(0.025, 0.975))
p  + geom_vline(xintercept=0, col="orange") # una schifezza!


################ 90 % CI RANDOM COEFFICIENTS ##################

##### b1 #######

tmp=data.out[,grep("b.1.", names(data.out), fixed=TRUE)] #; names(tmp)
b=tmp
tmp=melt(b) #;head(tmp)
colnames(tmp) = c("Parameter", "value")
CI.b = apply(b, 2, quantile, c(0.05, 0.95))#; CI.b
p=ggs_caterpillar(tmp, thick_ci = c(0.05, 0.95), thin_ci = c(0.025, 0.975))
p  + geom_vline(xintercept=0, col="orange") 
#ggsave('95CI_coefficients_IOP_model4.png', width = 5, height = 6)
rm(b)

##### b2 #######

tmp=data.out[,grep("b.2.", names(data.out), fixed=TRUE)] #; names(tmp)
b=tmp
tmp=melt(b) #;head(tmp)
colnames(tmp) = c("Parameter", "value")
CI.b = apply(b, 2, quantile, c(0.05, 0.95))#; CI.b
p=ggs_caterpillar(tmp, thick_ci = c(0.05, 0.95), thin_ci = c(0.025, 0.975))
p  + geom_vline(xintercept=0, col="orange") 
rm(b)


##### b3 #######

tmp=data.out[,grep("b.3.", names(data.out), fixed=TRUE)] #; names(tmp)
b=tmp
tmp=melt(b) #;head(tmp)
colnames(tmp) = c("Parameter", "value")
CI.b = apply(b, 2, quantile, c(0.05, 0.95))#; CI.b
p=ggs_caterpillar(tmp, thick_ci = c(0.05, 0.95), thin_ci = c(0.025, 0.975))
p  + geom_vline(xintercept=0, col="orange") 
rm(b)

##### b4 #######

tmp=data.out[,grep("b.4.", names(data.out), fixed=TRUE)] #; names(tmp)
b=tmp
tmp=melt(b) #;head(tmp)
colnames(tmp) = c("Parameter", "value")
CI.b = apply(b, 2, quantile, c(0.05, 0.95))#; CI.b
p=ggs_caterpillar(tmp, thick_ci = c(0.05, 0.95), thin_ci = c(0.025, 0.975))
p  + geom_vline(xintercept=0, col="orange") 
rm(b)
# this (macular volume) is fundamental!


##### b5 #######

tmp=data.out[,grep("b.5.", names(data.out), fixed=TRUE)] #; names(tmp)
b=tmp
tmp=melt(b) #;head(tmp)
colnames(tmp) = c("Parameter", "value")
CI.b = apply(b, 2, quantile, c(0.05, 0.95))#; CI.b
p=ggs_caterpillar(tmp, thick_ci = c(0.05, 0.95), thin_ci = c(0.025, 0.975))
p  + geom_vline(xintercept=0, col="orange") 
rm(b)


##### b6 #######

tmp=data.out[,grep("b.6.", names(data.out), fixed=TRUE)] #; names(tmp)
b=tmp
tmp=melt(b) #;head(tmp)
colnames(tmp) = c("Parameter", "value")
CI.b = apply(b, 2, quantile, c(0.05, 0.95))#; CI.b
p=ggs_caterpillar(tmp, thick_ci = c(0.05, 0.95), thin_ci = c(0.025, 0.975))
p  + geom_vline(xintercept=0, col="orange") 
rm(b)







########### BOXPLOT RANDOM COEFFICIENTS ############

#IOP,   b1 
tmp=stack(as.data.frame(data.out[,grep("b.1.", names(data.out), fixed=TRUE)]))
ggplot(tmp[,]) +  geom_boxplot(aes(x = ind, y = values))

#MD,                      b2 
tmp=stack(as.data.frame(data.out[,grep("b.2.", names(data.out), fixed=TRUE)]))
ggplot(tmp[,]) +  geom_boxplot(aes(x = ind, y = values))

#PSD,                     #b3 
tmp=stack(as.data.frame(data.out[,grep("b.3.", names(data.out), fixed=TRUE)]))
ggplot(tmp[,]) +  geom_boxplot(aes(x = ind, y = values))

#macular_volume,          b4
tmp=stack(as.data.frame(data.out[,grep("b.4.", names(data.out), fixed=TRUE)]))
ggplot(tmp[,]) +  geom_boxplot(aes(x = ind, y = values))

#Cup_area,                b5 
tmp=stack(as.data.frame(data.out[,grep("b.5.", names(data.out), fixed=TRUE)]))
ggplot(tmp[,]) +  geom_boxplot(aes(x = ind, y = values))

#Rim_area,                b6
tmp=stack(as.data.frame(data.out[,grep("b.6.", names(data.out), fixed=TRUE)]))
ggplot(tmp[,]) +  geom_boxplot(aes(x = ind, y = values))





########### means of coefficients ###########
beta.post <- data.out
beta.bayes  <- apply(beta.post,2,"mean")
plot(beta.bayes)
length(beta.bayes)
