#
#  mean_RNFL_thikness as response of several variables, both fixed and random coefficients
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
Sex
sex=rep(0,dim(Gdata_6vis)[1])
for(i in 1:dim(Gdata_6vis)[1] )
{
  sex[i]=ifelse(Sex[i]=="Male",1,0)
}
sex


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
data=list( y=RNFL_average, rim=Rim_area, sex=sex ,  npat= length(unique(Patient)), numerosity = numerosity, kk=kk)# dati che passo a jags



inits = function() {list( b0=80, b_sex=0, b_rim=rep(0,length(unique(Patient))), sigma=5)} # fixed intercept, random slope


modelRegress=jags.model("data3.bug",data=data,inits=inits,n.adapt=1000,n.chains=1)
update(modelRegress,n.iter=19000)
variable.names=c("b0","b_sex","b_rim", "sigma")
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
qplot(beta2.bayes,
      geom="histogram",
      binwidth = 10,  
      fill=I("blue"), 
      col=I("red"), 
      alpha=I(.2))



qplot(data.out$b_sex,
      geom="histogram",
      binwidth = 1,  
      fill=I("blue"), 
      col=I("red"), 
      alpha=I(.2))

ggplot(data.out, aes(data.out$b_sex))+
  geom_histogram( binwidth = 1,fill="#366699", 
                  col="lightgrey", alpha=I(.8))
