#
#  primo tentativo di covariate selection sulla base di 90% CI model5  
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
#summary(data.out)
#head(data.out)

# DIC
dic4<-dic.samples(modelRegress,n.iter=20000,thin=10)

# LPML
#lpml4=...

# RESIDUI BAYESIAN
#b_res4=...

#save.image("../R_object/model_6.RData")
rm(list=ls())
load("../R_object/model_6.RData")



library(coda)
library(plotrix)
outputRegress_mcmc <- as.mcmc(outputRegress)

quartz()
plot(outputRegress_mcmc)

quartz()
acfplot(outputRegress_mcmc)




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



