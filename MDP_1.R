#
# Primo tentativo di mixture di modello di Dirichlet processes: clustering di pazienti
#
# MODEL: RNFL_i = alpha_i + beta_i*time_i
#
# alpha_i ~ dnorm
# beta_i |P ~ P
# P ~ DP
#

rm(list=ls())
load("../R_object/Glaucoma_better_data.RData")
attach(mydata)


npat=length(unique(Patient)) # numero pazienti
npat


################JAGGAMENTO#######################

library(rjags)      # interfaccia R con JAGS
library(plotrix)    # per fare plot CIs
set.seed(1)         # fisso il seed per poter riprodurre

# Compute useful values 
numerosity<-as.integer(table(Patient))
kk=rep(0,length(unique(Patient)))
kk[1]=0
for (i in 2:length(unique(Patient))){
  kk[i]=kk[i-1]+numerosity[i-1]
}

Kmax=40 # troncamento serie Sethuraman

#  genero una lista con i dati da passare a JAGS
data <- list(y=RNFL_average, visit2=visit2, npat=npat, numerosity = numerosity, kk=kk, 
             b_slope=0,m1=0, Kmax=Kmax, alpha = 0.7,  tau2 = 0.01 )

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


variable.names=c("b0", "b1",  "sigma0", "sigma1","mu")
# Monitoro i parametri: 
n.iter=50000 
thin=10  

outGLMM_DP=coda.samples(model=LMM_DP_1,variable.names=variable.names,n.iter=n.iter,thin=thin)
# salvo l'intera catena dei parametri monitorati (si tratta di una lista mcmc)


save(outGLMM_DP,file='GLMMDP_output.Rdata')

