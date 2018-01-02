attach(Gdata)
patients<-factor(Patient)
a<-table(patients)
numerosity<-as.integer(a)
matrix1<-cbind(Patient,RNFL_average,Rim_area)
#flag=0;
#for(ii in 2:3){
#  for(jj in 1:1090){
#    if(matrix[ii,jj]==NA) flag=1;
#  }
#}
#flag
primi<-matrix(0,115,2)
ii=1;
iter=1;
while(iter<1090){
  primi[ii,1]<-matrix1[iter,2]
  primi[ii,2]<-matrix1[iter,3]
  iter<-iter+numerosity[ii];
  ii<-ii+1;
}
head(primi)

################JAGGAMENTO#######################

data=list(y=primi[,1],x=primi[,2],n=115)
inits = function() {list(beta=c(0,0),sigma=5)}

modelRegress=jags.model("data0.bug",data=data,inits=inits,n.adapt=1000,n.chains=1)
update(modelRegress,n.iter=19000)
variable.names=c("beta", "sigma")
n.iter=50000 
thin=10
outputRegress=coda.samples(model=modelRegress,variable.names=variable.names,n.iter=n.iter,thin=thin)

library(coda)
library(plotrix) 

data.out=as.matrix(outputRegress)
data.out=data.frame(data.out)
attach(data.out)
n.chain=dim(data.out)[1] 
n.chain
summary(data.out)
head(data.out)

beta.post <- data.out[,1:2]
beta.bayes  <- apply(beta.post,2,"mean")
beta.bayes

sig.post= data.out[,'sigma']

y1=primi[,1]; x1=primi[,2];
pippo= lm(y1 ~ x1)
summary(pippo)

pippo$coefficients
beta.bayes