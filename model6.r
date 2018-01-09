#
#  RNFL_average as response of several variables, both fixed and random coefficients
#  2nd attempt of covariate selection; a very simple model



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
        IOP,                 #beta8
        PSD,                 #beta9
        macular_volume,      #beta10
        Cup_area,            #beta11
        Rim_area             #beta12
        )       

#covariates with random coefficents
Z=cbind(
        macular_volume,          #b1 
        Cup_area                 #b2
        )

# Hyperparameters:
mu0=rep(0, dim(X)[2])
S0=diag(rep(1, dim(X)[2])) 
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
#summary(data.out)
#head(data.out)


########################### PLOTS COEFFICIENTS (POSTERIOR) ##########################
source ("multiplot.R")
library(ggplot2)

### FIXED COEFFICIENTS ###

# Intercept
ggplot(data.out, aes(data.out$beta.1.))+geom_histogram( binwidth = 5,fill="#366699", col="lightgrey", alpha=I(.8))

# Race
p1=ggplot(data.out, aes(data.out$beta.2.))+geom_histogram( binwidth = 4,fill="#366699", col="lightgrey", alpha=I(.8))+coord_cartesian(xlim = c(-20, 50),ylim=c(0,2000)) 
p2=ggplot(data.out, aes(data.out$beta.3.))+geom_histogram( binwidth = 4,fill="#366699", col="lightgrey", alpha=I(.8))+coord_cartesian(xlim = c(-20, 50),ylim=c(0,2000))
multiplot(p1,p2)

#familiarity
ggplot(data.out, aes(data.out$beta.4.))+geom_histogram( binwidth = 1,fill="#366699", col="lightgrey", alpha=I(.8))

#Sex
ggplot(data.out, aes(data.out$beta.5.))+geom_histogram( binwidth = 1,fill="#366699", col="lightgrey", alpha=I(.8))

#POAG
ggplot(data.out, aes(data.out$beta.6.))+geom_histogram( binwidth = 1,fill="#366699", col="lightgrey", alpha=I(.8))

#age65
ggplot(data.out, aes(data.out$beta.7.))+geom_histogram( binwidth = 1,fill="#366699", col="lightgrey", alpha=I(.8))

#IOP 
ggplot(data.out, aes(data.out$beta.8.))+geom_histogram( binwidth = 0.05,fill="#366699", col="lightgrey", alpha=I(.8))

#PSD
ggplot(data.out, aes(data.out$beta.9.))+geom_histogram( binwidth = 0.1,fill="#366699", col="lightgrey", alpha=I(.8))

#macular_volume
ggplot(data.out, aes(data.out$beta.10.))+geom_histogram( binwidth = 0.5,fill="#366699", col="lightgrey", alpha=I(.8))

#Cup_area
ggplot(data.out, aes(data.out$beta.11.))+geom_histogram( binwidth = 0.5,fill="#366699", col="lightgrey", alpha=I(.8))

#Rim_area 
ggplot(data.out, aes(data.out$beta.12.))+geom_histogram( binwidth = 0.4,fill="#366699", col="lightgrey", alpha=I(.8))







#### RANDOM COEFFICIENTS #####

#macular_volume    b1 
tmp=stack(as.data.frame(data.out[,grep("b.1.", names(data.out), fixed=TRUE)]))
ggplot(tmp[,]) +  geom_boxplot(aes(x = ind, y = values))

#cup_area          b2 
tmp=stack(as.data.frame(data.out[,grep("b.2.", names(data.out), fixed=TRUE)]))
ggplot(tmp[,]) +  geom_boxplot(aes(x = ind, y = values))





# means of coefficients
beta.post <- data.out
beta.bayes  <- apply(beta.post,2,"mean")
plot(beta.bayes)
length(beta.bayes)
