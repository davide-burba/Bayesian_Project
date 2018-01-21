load("../R_object/Glaucoma_better_data.RData")
attach(mydata)

numerosity<-as.integer(table(Patient))
kk=rep(0,length(unique(Patient)))
kk[1]=0
for (i in 2:length(unique(Patient))){
  kk[i]=kk[i-1]+numerosity[i-1]
}
kk<-c(kk,1046)

detection<-c()
flag<-0

for(ii in 1:length(numerosity)){
  for(jj in 1:numerosity[ii]){
    if(Progression[kk[ii]+jj]==1) flag=1
  }
  detection<-c(detection,flag)
  flag=0
}

table(detection)
