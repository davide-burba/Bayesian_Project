#
# This model is the same as model 4 but using RNFL scaled come risposta
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

#covariates with fixed coefficents
X=cbind(rep(1,length(Patient)),  #beta1 (intercept)
        Asian,               #beta2
        Black,               #beta3
        Hispanic,            #beta4 
        #White,               # --> White categoria di riferimento
        familiarity_yes,     #beta5
        Sex,                 #beta6
        POAG,                #beta7
        Hypertension,        #beta8
        HyperLipidemia,      #beta9
        Diabetes,            #beta10
        Cardiovascular_Dz,   #beta11
        Smoking,             #beta12
        age65,               #beta13
        Age,                 #beta14
        trusopt,             #beta15
        prostaglandin,       #beta16
        brimonidine,         #beta17
        timolol,             #beta18
        htnmed,              #beta19
        IOP,                     #beta20 
        acuity,                  #beta21
        MD,                      #beta22
        PSD,                     #beta23
        cup_disk_horiz_ratio,    #beta24
        cup_disk_vert_ratio,     #beta25
        macular_volume,          #beta26
        Vert_integrated_rim_area__vol_,     #beta27
        Horz_integrated_rim_width__area_,   #beta28 
        Cup_area,                #beta29
        Rim_area,                #beta30
        cup_disk_area_ratio      #beta31 
        # yearofglaucoma,      # Troppi NA, bisogna pensare a come trattarli (alcuni pazienti non hanno proprio dati)
)       


# Hyperparameters:
mu0=rep(0, dim(X)[2])
c=50
S0=c*diag(rep(1, dim(X)[2])) # usata identità!  solve(t(X)%*%X) non è invertibile!



# NO RANDOM
data=list( y=mydata$RNFL_average_scaled, X=X,  npat= length(unique(Patient)), 
           mu0=mu0, S0=S0,  numerosity = numerosity, kk=kk)      # dati che passo a jags


#fixed coefficients initialization:
beta= rep(0,dim(X)[2])


#inits = function() {list( beta=beta, b=b, invD=invD, sigma0=50,sigma1=50)} 
inits = function() {list( beta=beta, sigma0=50,sigma1=50)} 


modelRegress=jags.model("data5_norandom.bug",data=data,inits=inits,n.adapt=1000,n.chains=1)
update(modelRegress,n.iter=19000)
variable.names=c("beta", "sigma0","sigma1")
n.iter=50000 
thin=10 

library(coda)
library(plotrix)
outputRegress=coda.samples(model=modelRegress,variable.names=variable.names,n.iter=n.iter,thin=thin)

library(coda) 
#



data.out=as.matrix(outputRegress)
data.out=data.frame(data.out)
attach(data.out)
n.chain=dim(data.out)[1] 
n.chain
#summary(data.out)
#head(data.out)

#save.image("../R_object/model_4_c.RData")

rm(list=ls())
load("../R_object/model_4_c.RData")




########################### PLOTS COEFFICIENTS (POSTERIOR) ##########################


## 90 % CI (fixed) COEFFICIENTS ##

library(reshape)
library(ggmcmc)

n=32 # number of betas
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

p=ggs_caterpillar(tmp, thick_ci = c(0.05, 0.95), thin_ci = c(0.025, 0.975))

p + geom_vline(xintercept=0, col="orange")   #if you don't want colors, use this


# without intercept
p=ggs_caterpillar(tmp[which(tmp$Parameter!="Intercept"),], thick_ci = c(0.05, 0.95), thin_ci = c(0.025, 0.975))
p + geom_vline(xintercept=0, col="orange")   

#p  +aes(color=pos_neg) + geom_vline(xintercept=0, col="orange") +
# scale_color_manual(values = levels(pos_neg))

#ggsave('95CI_coefficients.png', width = 8, height = 4)








####### HISTOGRAMS (FIXED) COEFFICIENTS #######
source ("multiplot.R")

# beta0
ggplot(data.out, aes(data.out$beta.1.))+geom_histogram( binwidth = 5,fill="#366699", col="lightgrey", alpha=I(.8))
# beta1,2,3,4,5
p1=ggplot(data.out, aes(data.out$beta.2.))+geom_histogram( binwidth = 4,fill="#366699", col="lightgrey", alpha=I(.8))+coord_cartesian(xlim = c(-20, 50),ylim=c(0,1000)) 
p2=ggplot(data.out, aes(data.out$beta.3.))+geom_histogram( binwidth = 4,fill="#366699", col="lightgrey", alpha=I(.8))+coord_cartesian(xlim = c(-20, 50),ylim=c(0,1000))
p3=ggplot(data.out, aes(data.out$beta.4.))+geom_histogram( binwidth = 4,fill="#366699", col="lightgrey", alpha=I(.8))+coord_cartesian(xlim = c(-20, 50),ylim=c(0,1000))
p4=ggplot(data.out, aes(data.out$beta.5.))+geom_histogram( binwidth = 4,fill="#366699", col="lightgrey", alpha=I(.8))+coord_cartesian(xlim = c(-20, 50),ylim=c(0,1000))
multiplot(p1,p2,p3,p4)
# beta6
ggplot(data.out, aes(data.out$beta.6.))+geom_histogram( binwidth = 1,fill="#366699", col="lightgrey", alpha=I(.8))
# beta7
ggplot(data.out, aes(data.out$beta.7.))+geom_histogram( binwidth = 1,fill="#366699", col="lightgrey", alpha=I(.8))
# beta8
ggplot(data.out, aes(data.out$beta.8.))+geom_histogram( binwidth = 1,fill="#366699", col="lightgrey", alpha=I(.8))
# beta9
ggplot(data.out, aes(data.out$beta.9.))+geom_histogram( binwidth = 1,fill="#366699", col="lightgrey", alpha=I(.8))
# beta10
ggplot(data.out, aes(data.out$beta.10.))+geom_histogram( binwidth = 1,fill="#366699", col="lightgrey", alpha=I(.8))
# beta11
ggplot(data.out, aes(data.out$beta.11.))+geom_histogram( binwidth = 1,fill="#366699", col="lightgrey", alpha=I(.8))
# beta12
ggplot(data.out, aes(data.out$beta.12.))+geom_histogram( binwidth = 1,fill="#366699", col="lightgrey", alpha=I(.8))
# beta13
ggplot(data.out, aes(data.out$beta.13.))+geom_histogram( binwidth = 1,fill="#366699", col="lightgrey", alpha=I(.8))
# beta14
ggplot(data.out, aes(data.out$beta.14.))+geom_histogram( binwidth = 1,fill="#366699", col="lightgrey", alpha=I(.8))
# beta15
ggplot(data.out, aes(data.out$beta.15.))+geom_histogram( binwidth = 0.1,fill="#366699", col="lightgrey", alpha=I(.8))
# beta16
ggplot(data.out, aes(data.out$beta.16.))+geom_histogram( binwidth = 1,fill="#366699", col="lightgrey", alpha=I(.8))
# beta17
ggplot(data.out, aes(data.out$beta.17.))+geom_histogram( binwidth = 1,fill="#366699", col="lightgrey", alpha=I(.8))
# beta18
ggplot(data.out, aes(data.out$beta.18.))+geom_histogram( binwidth = 1,fill="#366699", col="lightgrey", alpha=I(.8))
# beta19
ggplot(data.out, aes(data.out$beta.19.))+geom_histogram( binwidth = 1,fill="#366699", col="lightgrey", alpha=I(.8))
# beta20
ggplot(data.out, aes(data.out$beta.20.))+geom_histogram( binwidth = 1,fill="#366699", col="lightgrey", alpha=I(.8))
# sigma_0
ggplot(data.out, aes(data.out$sigma0))+geom_histogram( binwidth = 1,fill="#366699", col="lightgrey", alpha=I(.8))
# sigma_1
ggplot(data.out, aes(data.out$sigma1))+geom_histogram( binwidth = 50,fill="#366699", col="lightgrey", alpha=I(.8))




##############  BOXPLOT ############# 
require(reshape2)
ggplot(data = melt(beta), aes(x=variable, y=value)) + 
  stat_boxplot(aes(fill=variable)) +
  coord_flip()+
  geom_hline(yintercept=0)

