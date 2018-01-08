#
#  RNFL_average as response of several variables, both fixed and random coefficients
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
        White,               #beta5
        familiarity_yes,     #beta6
        Sex,                 #beta7
        POAG,                #beta8
        Hypertension,        #beta9
        HyperLipidemia,      #beta10
        Diabetes,            #beta11
        Cardiovascular_Dz,   #beta12
        Smoking,             #beta13
        age65,               #beta14
        Age,                 #beta15
        trusopt,             #beta16
        prostaglandin,       #beta17
        brimonidine,         #beta18
        timolol,             #beta19
        htnmed               #beta20
        # yearofglaucoma,      # Troppi NA, bisogna pensare a come trattarli (alcuni pazienti non hanno proprio dati)
        )       

#covariates with random coefficents
Z=cbind(
        IOP,                     #b1 
        acuity,                  #b2 
        MD,                      #b3 
        PSD,                     #b4 
        cup_disk_horiz_ratio,    #b5 
        cup_disk_vert_ratio,     #b6 
        macular_volume,          #b7 
        Vert_integrated_rim_area__vol_,     #b8
        Horz_integrated_rim_width__area_,   #b9 
        Cup_area,                #b10
        Rim_area,                #b11
        cup_disk_area_ratio      #b12 
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

### FIXED COEFFICIENTS ###

# Intercept
ggplot(data.out, aes(data.out$beta.1.))+geom_histogram( binwidth = 5,fill="#366699", col="lightgrey", alpha=I(.8))

# Race
p1=ggplot(data.out, aes(data.out$beta.2.))+geom_histogram( binwidth = 4,fill="#366699", col="lightgrey", alpha=I(.8))+coord_cartesian(xlim = c(-20, 50),ylim=c(0,1000)) 
p2=ggplot(data.out, aes(data.out$beta.3.))+geom_histogram( binwidth = 4,fill="#366699", col="lightgrey", alpha=I(.8))+coord_cartesian(xlim = c(-20, 50),ylim=c(0,1000))
p3=ggplot(data.out, aes(data.out$beta.4.))+geom_histogram( binwidth = 4,fill="#366699", col="lightgrey", alpha=I(.8))+coord_cartesian(xlim = c(-20, 50),ylim=c(0,1000))
p4=ggplot(data.out, aes(data.out$beta.5.))+geom_histogram( binwidth = 4,fill="#366699", col="lightgrey", alpha=I(.8))+coord_cartesian(xlim = c(-20, 50),ylim=c(0,1000))
multiplot(p1,p2,p3,p4)

#familiarity
ggplot(data.out, aes(data.out$beta.6.))+geom_histogram( binwidth = 1,fill="#366699", col="lightgrey", alpha=I(.8))

#Sex
ggplot(data.out, aes(data.out$beta.7.))+geom_histogram( binwidth = 1,fill="#366699", col="lightgrey", alpha=I(.8))

#POAG
ggplot(data.out, aes(data.out$beta.8.))+geom_histogram( binwidth = 1,fill="#366699", col="lightgrey", alpha=I(.8))

#Hypertension
ggplot(data.out, aes(data.out$beta.9.))+geom_histogram( binwidth = 1,fill="#366699", col="lightgrey", alpha=I(.8))

#HyperLipidemia
ggplot(data.out, aes(data.out$beta.10.))+geom_histogram( binwidth = 1,fill="#366699", col="lightgrey", alpha=I(.8))

#Diabetes
ggplot(data.out, aes(data.out$beta.11.))+geom_histogram( binwidth = 1,fill="#366699", col="lightgrey", alpha=I(.8))

#Cardiovascular_Dz
ggplot(data.out, aes(data.out$beta.12.))+geom_histogram( binwidth = 1,fill="#366699", col="lightgrey", alpha=I(.8))

#Smoking
ggplot(data.out, aes(data.out$beta.13.))+geom_histogram( binwidth = 1,fill="#366699", col="lightgrey", alpha=I(.8))

#age65 
ggplot(data.out, aes(data.out$beta.14.))+geom_histogram( binwidth = 1,fill="#366699", col="lightgrey", alpha=I(.8))

#Age
ggplot(data.out, aes(data.out$beta.15.))+geom_histogram( binwidth = 0.1,fill="#366699", col="lightgrey", alpha=I(.8))

#trusopt
ggplot(data.out, aes(data.out$beta.16.))+geom_histogram( binwidth = 1,fill="#366699", col="lightgrey", alpha=I(.8))

#prostaglandin
ggplot(data.out, aes(data.out$beta.17.))+geom_histogram( binwidth = 1,fill="#366699", col="lightgrey", alpha=I(.8))

#brimonidine
ggplot(data.out, aes(data.out$beta.18.))+geom_histogram( binwidth = 1,fill="#366699", col="lightgrey", alpha=I(.8))

#timolol
ggplot(data.out, aes(data.out$beta.19.))+geom_histogram( binwidth = 1,fill="#366699", col="lightgrey", alpha=I(.8))

#htnmed
ggplot(data.out, aes(data.out$beta.20.))+geom_histogram( binwidth = 1,fill="#366699", col="lightgrey", alpha=I(.8))

#sigma0
ggplot(data.out, aes(data.out$sigma0))+geom_histogram( binwidth = 1,fill="#366699", col="lightgrey", alpha=I(.8))

#sigma1
ggplot(data.out, aes(data.out$sigma1))+geom_histogram( binwidth = 50,fill="#366699", col="lightgrey", alpha=I(.8))



#### RANDOM COEFFICIENTS #####

#IOP,   b1 
tmp=stack(as.data.frame(data.out[,grep("b.1.", names(data.out), fixed=TRUE)]))
ggplot(tmp[,]) +  geom_boxplot(aes(x = ind, y = values))

#acuity,                  b2 
tmp=stack(as.data.frame(data.out[,grep("b.2.", names(data.out), fixed=TRUE)]))
ggplot(tmp[,]) +  geom_boxplot(aes(x = ind, y = values))

#MD,                      b3 
tmp=stack(as.data.frame(data.out[,grep("b.3.", names(data.out), fixed=TRUE)]))
ggplot(tmp[,]) +  geom_boxplot(aes(x = ind, y = values))

#PSD,                     #b4 
tmp=stack(as.data.frame(data.out[,grep("b.4.", names(data.out), fixed=TRUE)]))
ggplot(tmp[,]) +  geom_boxplot(aes(x = ind, y = values))

#cup_disk_horiz_ratio,    #b5 
tmp=stack(as.data.frame(data.out[,grep("b.5.", names(data.out), fixed=TRUE)]))
ggplot(tmp[,]) +  geom_boxplot(aes(x = ind, y = values))

#cup_disk_vert_ratio,     b6 
tmp=stack(as.data.frame(data.out[,grep("b.6.", names(data.out), fixed=TRUE)]))
ggplot(tmp[,]) +  geom_boxplot(aes(x = ind, y = values))

#macular_volume,          b7
tmp=stack(as.data.frame(data.out[,grep("b.7.", names(data.out), fixed=TRUE)]))
ggplot(tmp[,]) +  geom_boxplot(aes(x = ind, y = values))


#Vert_integrated_rim_area__vol_,     b8 
tmp=stack(as.data.frame(data.out[,grep("b.8.", names(data.out), fixed=TRUE)]))
ggplot(tmp[,]) +  geom_boxplot(aes(x = ind, y = values))

#Horz_integrated_rim_width__area_,   b9 
tmp=stack(as.data.frame(data.out[,grep("b.9.", names(data.out), fixed=TRUE)]))
ggplot(tmp[,]) +  geom_boxplot(aes(x = ind, y = values))

#Cup_area,                b10   
tmp=stack(as.data.frame(data.out[,grep("b.10.", names(data.out), fixed=TRUE)]))
ggplot(tmp[,]) +  geom_boxplot(aes(x = ind, y = values))

#Rim_area,                b11
tmp=stack(as.data.frame(data.out[,grep("b.11.", names(data.out), fixed=TRUE)]))
ggplot(tmp[,]) +  geom_boxplot(aes(x = ind, y = values))

#cup_disk_area_ratio,     b12
tmp=stack(as.data.frame(data.out[,grep("b.12.", names(data.out), fixed=TRUE)]))
ggplot(tmp[,]) +  geom_boxplot(aes(x = ind, y = values))





# means of coefficients
beta.post <- data.out
beta.bayes  <- apply(beta.post,2,"mean")
plot(beta.bayes)
length(beta.bayes)
