#
#  RNFL as response of RIM, all the obs are considered (random slope)
#





rm(list=ls())
load("../R_object/Glaucoma_data_full.RData")


attach(Gdata_full)
patients<-factor(Patient)
a<-table(patients)
tmp=unique(patients)
tmp_data=tmp[which(a>=6)] # consider only patients with at least 6 visits
length(tmp_data)
Gdata_6vis=Gdata_full[Patient %in% tmp_data,]
rm(tmp)
rm(tmp_data)
rm(patients)
detach(Gdata_full)
attach(Gdata_6vis)
rm(Gdata_full)


patients<-factor(Patient)#recompute values with new data
a<-table(patients)
numerosity<-as.integer(a)

#matrix1<-cbind(Patient,RNFL_average,Rim_area)
#head(matrix1)







################JAGGAMENTO#######################
library(rjags)

kk=rep(0,length(unique(Patient)))
kk[1]=0
for (i in 2:length(unique(Patient))){
  kk[i]=kk[i-1]+numerosity[i-1]
}
data=list( y=RNFL_average, x=Rim_area, npat= length(unique(Patient)), numerosity = numerosity, kk=kk)# dati che passo a jags



inits = function() {list( beta=rep(0,length(unique(Patient))), sigma=5)} # fixed intercept, random slope


modelRegress=jags.model("data2.bug",data=data,inits=inits,n.adapt=1000,n.chains=1)
update(modelRegress,n.iter=19000)
variable.names=c("beta", "sigma")
n.iter=50000 
thin=10 

library(coda)
library(plotrix)
outputRegress=coda.samples(model=modelRegress,variable.names=variable.names,n.iter=n.iter,thin=thin)

 
outputRegress_mcmc <- as.mcmc(outputRegress)
x11()
plot(outputRegress_mcmc)

data.out=as.matrix(outputRegress)
data.out=data.frame(data.out)
attach(data.out)
n.chain=dim(data.out)[1] 
n.chain
summary(data.out)
head(data.out)

beta.post <- data.out
beta.bayes  <- apply(beta.post,2,"mean")
beta2.bayes=beta.bayes[2:103]
hist(beta2.bayes, breaks=5)
qplot(beta2.bayes,
      geom="histogram",
      binwidth = 20,  
      fill=I("blue"), 
      col=I("red"), 
      alpha=I(.2))




