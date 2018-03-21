#
# Primo tentativo di mixture di modello di Dirichlet processes: clustering di pazienti
#
# MODEL: RNFL_ij = alpha_i + beta_i*time_ij
#
# alpha_i ~ dnorm
# beta_i |P ~ P
# P ~ DP
#
# alpha=0.5
# n.iter=200000 
# thin=40  
#

rm(list=ls())
load("../R_object/Glaucoma_better_data.RData")
attach(mydata)


npat=length(unique(Patient)) # numero pazienti
npat


################JAGGAMENTO#######################

library(rjags)      # interfaccia R con JAGS
library(plotrix)    # per fare plot CIs
library(coda)
set.seed(1)         # fisso il seed per poter riprodurre

# Compute useful values 
numerosity<-as.integer(table(Patient))
kk=rep(0,length(unique(Patient)))
kk[1]=0
for (i in 2:length(unique(Patient))){
  kk[i]=kk[i-1]+numerosity[i-1]
}

Kmax=40 # troncamento serie Sethuraman

#  genero una lista con i dati da passare a JAGS (tra cui IPERPARAMETRI)
data <- list(y=RNFL_average, visit2=visit2, npat=npat, numerosity = numerosity, kk=kk, 
             b_slope=0,m1=0, Kmax=Kmax, alpha = 0.5,  tau2 = 0.01 )

###########  
# definisco una lista con lo stato iniziale della catena
# fornisco anche il seed e il tipo di generatore di numeri casuali




inits = list(b0=rep(mean(RNFL_average),npat),
             sigma0=50,
             sigma1=50,
             .RNG.seed = 2,
             .RNG.name = 'base::Wichmann-Hill'
)

#########################
### CREATE JAGS MODEL ###
#########################

LMM_DP_1=jags.model("LMM_DP_1.bug",data=data,inits=inits,n.adapt=1000,n.chains=1) 


# Faccio un update del modello per il numero di iterazioni specificato SENZA salvare nulla 
update(LMM_DP_1,19000)

# Monitoro i parametri: 
variable.names=c("b0", "slope", "sigma0", "sigma1","mu")

n.iter=200000 
thin=40  

outputRegress=coda.samples(model=LMM_DP_1,variable.names=variable.names,n.iter=n.iter,thin=thin)
# salvo l'intera catena dei parametri monitorati (si tratta di una lista mcmc)


save.image("../R_object/MDP_2.RData")


load("../R_object/MDP_2.RData")




#####################
## GODNESS OF MCMC ##
#####################

library(coda)        # pacchetto per analizzare catene
library(plotrix)     # per fare plot CIs



data=as.matrix(outputRegress) # trasformo il dataframe in matrice 
data=data.frame(data)
attach(data)
n.chain=dim(data)[1]   # lunghezza catena (final sample size)

# let's check the traceplots of random slopes
outputRegress_mcmc = as.mcmc(data[,grep("slope", names(data), fixed=TRUE)])

quartz()
plot(outputRegress_mcmc[,1:2])

# some autocorrelations plots
quartz()
acfplot(outputRegress_mcmc[,2:5]) # increment lag!
acfplot(outputRegress_mcmc[,11:20])
acfplot(outputRegress_mcmc[,60:63])

# let's check the traceplots of random intercept
outputRegress_mcmc = as.mcmc(data[,grep("b0", names(data), fixed=TRUE)])


# some autocorrelations plots
quartz()
acfplot(outputRegress_mcmc[,1:10]) # increment thinning!
acfplot(outputRegress_mcmc[,40:50])
acfplot(outputRegress_mcmc[,90:100])



# let's check the traceplots of others
outputRegress_mcmc = as.mcmc(data[,grep("sigma", names(data), fixed=TRUE)])
quartz()
plot(outputRegress_mcmc)

#autocorrelation
quartz()
acfplot(outputRegress_mcmc)

##################### 90% CI RANDOM COEFFICIENTS #####################
library(ggplot2)
library(ggmcmc)
library(reshape)

# 90% CI random coefficients: random intercept
tmp=data[,grep("b0", names(data), fixed=TRUE)]
names(tmp)
b=tmp
tmp=melt(b)
head(tmp)
colnames(tmp) = c("Parameter", "value")
CI.b = apply(b, 2, quantile, c(0.05, 0.95)) 
CI.b
p=ggs_caterpillar(tmp, thick_ci = c(0.05, 0.95), thin_ci = c(0.025, 0.975))
p  + geom_vline(xintercept=0, col="orange") 
#ggsave('95CI_coefficients_intercept.png', width = 5, height = 6)

# 90% CI random coefficients
tmp=data[,grep("slope", names(data), fixed=TRUE)]
names(tmp)
b=tmp
tmp=melt(b)
head(tmp)
colnames(tmp) = c("Parameter", "value")
CI.b = apply(b, 2, quantile, c(0.05, 0.95)) 
CI.b
p=ggs_caterpillar(tmp, thick_ci = c(0.05, 0.95), thin_ci = c(0.025, 0.975))
p  + geom_vline(xintercept=0, col="orange") 





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

error=abs(pred.mean- RNFL_average)
mean.error=mean(error)
mean(error)
median.error=median(error)
median.error


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
b0= data[,grep("b0.1.", names(data), fixed=TRUE)]
slope= data[,grep("slope.1.", names(data), fixed=TRUE)]
pred= b0  + slope*visit2[2] 
1/mean(1/dnorm(RNFL_average[2], pred , data[,grep("sigma0", names(data), fixed=TRUE)])) 
# predizione vs valore reale
mean(pred)
RNFL_average[2]


#calcolo LPML
n=dim(mydata)[1]
n
cpo=seq(1,n)
first.b0.pos=which(names(data)=="b0.1.")
first.slope.pos=which(names(data)=="slope.1.")
mcmc.std.dev=data[,grep("sigma0", names(data), fixed=TRUE)]

for (i in 1:npat){
  b0=data[,first.b0.pos-1+i]
  slope=data[,first.slope.pos-1+i]
  for(j in 1:numerosity[i])
  {
    cpo[kk[i]+j]=1/mean(1/dnorm(RNFL_average[kk[i]+j], b0  + slope*visit2[kk[i]+j] , mcmc.std.dev))
  }
}
LPML=sum(log(cpo))
LPML
sum(cpo)
median(cpo)
mean(cpo)

# BIC & AIC (using the means of coefficients as summary of posterior)
n=dim(mydata)[1]
first.b0.pos=which(names(data)=="b0.1.")
first.slope.pos=which(names(data)=="slope.1.")
mcmc.std.dev=data[,grep("sigma0", names(data), fixed=TRUE)]

for (i in 1:npat){
  b0=data[,first.b0.pos-1+i]
  slope=data[,first.slope.pos-1+i]
  for(j in 1:numerosity[i])
  {
    pred[kk[i]+j]= mean(b0  + slope*visit2[kk[i]+j])
  }
}
sigma.hat=mean(data[,grep("sigma0", names(data), fixed=TRUE)])
loglik=0
for (k in 1:n){
  loglik=loglik+dnorm(RNFL_average[k],mean=pred[k], sd =sigma.hat , log=TRUE)
}

# loglikelyhood:
loglik

# r is the number of parameters
r=length(names(data[,grep("b0.", names(data), fixed=TRUE)])) + length(names(data[,grep("b1.", names(data), fixed=TRUE)]))

#AIC
AIC=2*loglik-2*r
AIC

#BIC
BIC=2*loglik-r*log(n)
BIC






##################################
# Cluster estimate di Lau&Green
################################

label.mat = as.matrix(data[,grep("slope", names(data), fixed=TRUE)]) # extract cluster labels
dim(label.mat)
G=n.chain
pihat <- matrix(0,ncol=npat,nrow=npat)
for(i in 1:G){
  ss <- label.mat[i,]
  cij <- outer(ss,ss,'==')*1
  pihat <- pihat+cij
}
pihat <- pihat/G # pihat is the similarity matrix


#####Binder loss function
FF <- vector("numeric")
K <- 0.5 #prova con K=0.5, >=0.7
for(i in 1:G){
  ss <- label.mat[i,] 
  cij <- outer(ss,ss,'==')
  pluto <- (pihat-K)*as.matrix(cij) # nb: il * fa elemento per elemento, non riga per colonna
  pluto <-  pluto[upper.tri(pluto)]
  FF[i] <- sum(pluto)
}

#plot(FF, type="l")

# seleziono come stima dei cluster l'iterazione che minimizza la binder loss
ind.bind <- which.max(FF)[1] 
label.mat[ind.bind,]#leti(ind.bind,L,G,m)
#plot(FF)
ll.bind <- label.mat[ind.bind,] #leti(ind.bind,L,G,m)
unici <- unique(ll.bind)
unici=sort(unici)
unici
l.uni <- length(unici)# numero di gruppi stimato 
l.uni

ncl=l.uni

table(sort(ll.bind))


########################
## "ANALISI" CLUSTER
########################

# Qui si confrontano le caratteristiche dei pazienti nello stesso cluster,
# ad esempio se rispecchiano divisione fatta con classificatore naive 
# (variabile ProgressionStructure all'ultima visita) 

 handy_patient=Patient 
 for (i in 1:npat){
   for(j in 1:numerosity[i])
   {
     handy_patient[kk[i]+j]=i
   }
 }
# handy_patient # pazienti con ID da 1 a 104, più facili da maneggiare



prog=rep(0,npat)
for (i in 1:npat){
  prog[i]=ProgressionStructure[kk[i]+numerosity[i]]
}
prog


mean(prog)*100

# cluster 1
unici[1]
tmp=which(ll.bind==unici[1])
progcl1=prog[tmp]
mean(progcl1)*100 
progcl1
ind.cl1=which(handy_patient %in%tmp)
ind.cl1
unique(Sex[ind.cl1])
Race[ind.cl1]


# cluster 2
unici[2]
tmp=which(ll.bind==unici[2])
progcl2=prog[tmp]
mean(progcl2)*100 
progcl2

# cluster 3
unici[3]
tmp=which(ll.bind==unici[3])
progcl3=prog[tmp]
mean(progcl3)*100 # tutti progrediti
progcl3

# cluster 4
unici[4]
tmp=which(ll.bind==unici[4])
progcl4=prog[tmp]
mean(progcl4)*100 #
progcl4

# cluster 5
unici[5]
tmp=which(ll.bind==unici[5])
progcl5=prog[tmp]
mean(progcl5)*100 #
progcl5


