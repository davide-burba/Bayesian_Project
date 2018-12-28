#
# This script evaluates the WAIC (and AIC) value; doesn't work well, use
# something else preferably
# 


#EDITABLE VALUES:
n=5000 # number of iterations
names(data.out)[1059]
beta_pos_interval=1:12
sigma_pos=1059

  
  
p=length(unique(Patient))
Mfinal=matrix(0,n,p)
Y=RNFL_average
d=dim(Z)[2]
  
numerosity<-as.integer(table(Patient))
kk=rep(0,length(unique(Patient)+1))
kk[1]=0
for (i in 2:length(unique(Patient))){
  kk[i]=kk[i-1]+numerosity[i-1]
}
kk[length(unique(Patient))+1]  =length(Patient)

library(mnormt)

for(ii in 1:n){
  betas<-data.out[ii,beta_pos_interval]
  for(jj in 1:p){
    Xi<-X[(kk[jj]+1):kk[jj+1],]
    #Zi<-Z[(kk[jj]+1):kk[jj+1],]
    #bi<-data.out[ii,((jj-1)*d+1):(jj*d)]
    Yi<-Y[(kk[jj]+1):kk[jj+1]]
    Sig<-diag(rep(1,numerosity[jj]))*data.out[ii,sigma_pos]
    m<-Xi%*%t(betas)#+Zi%*%t(bi)
    m<-as.vector(m)
    Mfinal[ii,jj]=dmnorm(Yi, m, Sig , log=TRUE )
  }
}

library(loo)
ww<-waic(Mfinal)
ww
loo(Mfinal)

AIC=sum(colMeans(Mfinal))
AIC
