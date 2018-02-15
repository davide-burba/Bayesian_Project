#COVARIATE SELECTION USANDO S&S RELATIVA AL "MODELLO 4"
# NB:  modello 4 è stato modificato!


load("../R_object/Glaucoma_better_data.RData")
attach(mydata)

library("BoomSpikeSlab")
library("coda")


a=cbind(rep(1,length(Patient)),  #beta1 (intercept)
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
)

b=cbind(
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
  cup_disk_area_ratio,      #b12 
  visit2
)

W=cbind(a,b)
W=W[,2:dim(W)[2]]
nomi=colnames(W)
W=data.frame(W)
covariates<-as.matrix(W)
#Y=RNFL_average
L=cbind(W,RNFL_average)
L=data.frame(L)
ypsilon=L[,dim(L)[2]]
#SPIKE N SLAB PRIOR
pr<-SpikeSlabPrior(x=covariates,y=L[,dim(L)[2]],expected.model.size=10,expected.r2 = 0.8)
#Abbiamo 31 covariate, ho messo come expected.model.size 10 (numero di covariate che ci aspettiamo
#di ottenere); questo valore ? una nostra prior belief ed influisce abbastanza sulle inclusion probabilites
#finali.

model1<-lm.spike(ypsilon~covariates,niter=40000)
summary(model1)

#Diagnostics
beta = model1$beta
dim(beta)
beta.mcmc = as.mcmc(beta)
x11()
#plot(beta.mcmc, ask = T)

res = residuals(model1)
dim(res)
hist(res)

#Useful plots
plot(model1, "inclusion", inc = 0.5)
plot(model1, "coefficients", inc = 0.5)
plot(model1, "residuals")
plot(model1, "size")

#CRITERI PER LA SELEZIONE DELLE COVARIATE:
#1)Median Probability (seleziono le covariate con probabilit? di incidenza >0.5)
coeff1 = summary.lm.spike(model1, burn = 5000)$coefficients
colnames(coeff1)
inc_prob1 = coeff1[, 5]
head(inc_prob1)
selected_cov = which(inc_prob1 > 0.5)
names(selected_cov)

#2)Highest Posterior Density
I1 = model1$beta > 0
dim(I1)
labels = numeric(dim(I1)[1])
for(i in 1:dim(I1)[1])
{
  labels[i] = sum( (1:32)*(I1[i,])*2^(1:32) ) # sum( i*gamma_i*2^i)
}
#labels
which.max(table(labels))
idxs = which(labels == 85958115352)
selezionati = which(I1[idxs[1], ]==TRUE)
names(selezionati)
