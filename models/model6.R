#
#  Best linear model with only fixed coefficients
#
# NOTE:   Zellner prior, only fixed coefficients, c=50
#
# We found only 



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
        Hispanic,            #beta3 
        #White,               # --> White categoria di riferimento
        familiarity_yes,     #beta4
        Sex,                 #beta5
        POAG,                #beta6
        Hypertension,        #beta7
        HyperLipidemia,      #beta8
        Cardiovascular_Dz,   #beta9
        age65,               #beta10
        prostaglandin,       #beta11
        brimonidine,         #beta12
        timolol,             #beta13
        IOP,                     #beta14 
        MD,                      #beta15
        macular_volume,          #beta16
        Vert_integrated_rim_area__vol_,     #beta17
        Horz_integrated_rim_width__area_,   #beta18 
        Rim_area               #beta19
        
)
# Hyperparameters:
mu0=rep(0, dim(X)[2])
c=50
S0=c*solve(t(X)%*%X) 



# NO RANDOM
data=list( y=RNFL_average, X=X,  npat= length(unique(Patient)), 
           mu0=mu0, S0=S0,  numerosity = numerosity, kk=kk)      # dati che passo a jags


#fixed coefficients initialization:
beta= rep(0,dim(X)[2])


#inits = function() {list( beta=beta, b=b, invD=invD, sigma0=50,sigma1=50)} 
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


save.image("../R_object/model_6.RData")
rm(list=ls())
load("../R_object/model_6.RData")




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
# let's check the traceplots of random slopes
outputRegress_mcmc = as.mcmc(data[,grep("beta.", names(data), fixed=TRUE)])

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

#mcmc is good




########################### PLOTS COEFFICIENTS (POSTERIOR) ##########################


## 90 % CI (fixed) COEFFICIENTS ##

library(reshape)
library(ggmcmc)


names(data.out)
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
p + geom_vline(xintercept=0, col="orange")   #if you don't want colors, use this


# without intercept
p=ggs_caterpillar(tmp[which(tmp$Parameter!="Intercept"),], thick_ci = c(0.05, 0.95), thin_ci = c(0.025, 0.975))
p + geom_vline(xintercept=0, col="orange")   






# If you want colours (EDIT covariates and n)
n=32 # number of betas
pos_neg=rep("grey",n)
tmp2=unique(tmp$Parameter)
for(i in 1:n) # edit looking at CI.beta! blue if both positive, red if both negative
{
  if (tmp2[i] %in% c("Black","familiarity_yes", "Hypertension", "HyperLipidemia", "Cardiovascular_Dz", "prostaglandin","IOP","MD", "cup_disk_horiz_ratio","macular_volume","Horz_integrated_rim_width__area_", "Rim_area") )
  {
    pos_neg[i]="blue"
  }
  if (tmp2[i] %in% c("Hispanic", "Sex", "POAG", "age65", "brimonidine", "timolol", "Vert_integrated_rim_area__vol_") )
  {
    pos_neg[i]="red"
  }
  
}
pos_neg=factor(pos_neg)
pos_neg

#p  +aes(color=pos_neg) + geom_vline(xintercept=0, col="orange") +
# scale_color_manual(values = levels(pos_neg))

#ggsave('95CI_coefficients.png', width = 8, height = 4)


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

plot(sort(bres^2)) # ce n'è qualcuno assurdo, che alza la media


## CROSS VALIDATION

# es CPOi togliendo seconda visita del primo paziente
beta= data[,grep("beta", names(data), fixed=TRUE)]
tmp=X%*%t(beta)
pred= tmp[2,]
1/mean(1/dnorm(RNFL_average[2], pred , data[,grep("sigma0", names(data), fixed=TRUE)])) 
# predizione vs valore reale
mean(pred)
RNFL_average[2]


#calcolo LPML
n=dim(mydata)[1]
n
cpo=seq(1,n)
tmp=X%*%t(beta)
mcmc.std.dev=data[,grep("sigma0", names(data), fixed=TRUE)]

for (i in 1:n){
  pred=tmp[i,]
  cpo[i]=1/mean(1/dnorm(RNFL_average[i], pred , mcmc.std.dev))
}

LPML=sum(log(cpo))
LPML
sum(cpo)
median(cpo)
mean(cpo)
min(cpo)

# BIC & AIC (using the means of coefficients as summary of posterior)
n=dim(mydata)[1]
first.b0.pos=which(names(data)=="b0.1.")
first.slope.pos=which(names(data)=="slope.1.")
mcmc.std.dev=data[,grep("sigma0", names(data), fixed=TRUE)]
tmp=X%*%t(beta)
pred=colMeans(t(tmp))
sigma.hat=mean(data[,grep("sigma0", names(data), fixed=TRUE)])
loglik=0
for (k in 1:n){
  loglik=loglik+dnorm(RNFL_average[k],mean=pred[k], sd =sigma.hat , log=TRUE)
}

# loglikelyhood:
loglik

# r is the number of parameters
r=19

#AIC
AIC=2*loglik-2*r
AIC

#BIC
BIC=2*loglik-r*log(n)
BIC





