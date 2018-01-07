#
# This file is used to save a more easy-to-use and "correct" dataset; 
# It also fills some NA values and removes strange observations
#
# Added variables: Asian,Black,Hispanic,White, familiarity_yes, Meds_according_to_Pt_bool,
#                  Ocular_Meds_according_to_Pt_bool, BaselineRNFL1
#
# Riempiti NA di: RNFL_average, RIM_area
#
# Note: - considerati solo pazienti con almeno 6 visite
#       - eliminato paziente 501: valori strani di RNFL (più precisamente di (RNFL-baseline)/baseline)


rm(list=ls())

mydata=read.csv("../data/data.csv", header=T,fill=T,sep=",")


# Consider only patients with at least 6 visits:
tmp=unique(mydata$Patient)
a<-table(mydata$Patient)
tmp_data=tmp[which(a>=6)] 
mydata=mydata[mydata$Patient %in% tmp_data,]
rm(tmp)
rm(tmp_data)
rm(a)


# Modify Sex covariate: Male=1, Female=0
tmp = factor(mydata$Sex, levels=c("Male","Female"), labels=c(1,0))
mydata$Sex=as.numeric(as.character(tmp))
rm(tmp)
mydata$Sex


# Race: build 4 vectors of boolean variables Asian Black Hispanic White 
attach(mydata)
Asian=(Race=="Asian")*1
Black=(Race=="Black")*1
Hispanic=(Race=="Hispanic")*1
White=(Race=="White")*1
detach(mydata)
mydata=cbind(mydata,Asian,Black,Hispanic,White)
mydata$Asian


#  Type_of_glaucoma: identify patients with POAG glaucoma (many NA)
attach(mydata)
POAG=(Type_of_glaucoma=="POAG")*1
detach(mydata)
mydata=cbind(mydata,POAG)
mydata$POAG


#  Allergies: trascuriamolo per ora
unique(mydata$Allergies)


#  Modify Hypertension: Yes=1, No=0
tmp = factor(mydata$Hypertension, levels=c("Yes","No"), labels=c(1,0))
mydata$Hypertension=as.numeric(as.character(tmp))
rm(tmp)
mydata$Hypertension


#  Modify HyperLipidemia: Yes=1, No=0
tmp = factor(mydata$HyperLipidemia, levels=c("Yes","No"), labels=c(1,0))
mydata$HyperLipidemia=as.numeric(as.character(tmp))
rm(tmp)
mydata$HyperLipidemia


#  Modify Diabetes: Yes=1, No=0  (from now on in one row)
mydata$Diabetes=as.numeric(as.character(factor(mydata$Diabetes, levels=c("Yes","No"), labels=c(1,0))))
mydata$Diabetes


#  Modify Cardiovascular_Dz: Yes=1, No=0  (now in one row)
mydata$Cardiovascular_Dz=as.numeric(as.character(factor(mydata$Cardiovascular_Dz, levels=c("Yes","No"), labels=c(1,0))))
mydata$Cardiovascular_Dz


# Other_systemic_dz_1: trascuriamolo per ora
unique(mydata$Other_systemic_dz_1)


# family history: since there are many NA, we made a vector with  1=yes, 0  otherwise
attach(mydata)
familiarity_yes=(family_history=="yes")*1
familiarity_yes
detach(mydata)
mydata=cbind(mydata,familiarity_yes)
mydata$familiarity_yes


# Smoking: Yes=1, No=0  
mydata$Smoking=as.numeric(as.character(factor(mydata$Smoking, levels=c("Yes","No"), labels=c(1,0))))
mydata$Smoking


#Meds_according_to_Pt: many different; for the moment 1="medicine yes", 0="medicine no"
attach(mydata)
Meds_according_to_Pt_bool=(Meds_according_to_Pt==""|Meds_according_to_Pt=="none" )*1
Meds_according_to_Pt_bool
detach(mydata)
mydata=cbind(mydata,Meds_according_to_Pt_bool)
mydata$Meds_according_to_Pt_bool


#Ocular_Meds_according_to_Pt: many different; for the moment 1=" ocular medicine yes", 0="ocul. medicine no"
attach(mydata)
Ocular_Meds_according_to_Pt_bool=(Ocular_Meds_according_to_Pt==""|Ocular_Meds_according_to_Pt=="none" )*1
Ocular_Meds_according_to_Pt_bool
detach(mydata)
mydata=cbind(mydata,Ocular_Meds_according_to_Pt_bool)
mydata$Ocular_Meds_according_to_Pt_bool


# BASELINE RNFL (RNFL AT FIRST VISIT)
attach(mydata)
BaselineRNFL1=rep(0,1050)
for (i in 1:dim(mydata)[1])
{
  if (visit2[i]==0)
  {
    tmp=RNFL_average[i] 
  }
  BaselineRNFL1[i]=tmp
}
detach(mydata)
mydata=cbind(mydata,BaselineRNFL1)
mydata$BaselineRNFL1



##
## RIEMPIMENTO NA: interpolazione lineare più stesse osservazioni alle code
##

library(zoo)

#prova:
tmp=na.approx(c(NA,1,NA,12,3,NA,76,NA,NA),na.rm=FALSE)
tmp
na.fill(tmp,"extend")
rm(tmp)

# riempi NA di RNFL_AVERAGE e RIM_AREA 
numerosity=table(mydata$Patient)
kk=0
for(i in 1:length(unique(mydata$Patient))){
  mydata$RNFL_average[(kk+1):(kk+numerosity[i])]=na.approx(mydata$RNFL_average[(kk+1):(kk+numerosity[i])],na.rm=FALSE)
  mydata$RNFL_average[(kk+1):(kk+numerosity[i])]=na.fill(mydata$RNFL_average[(kk+1):(kk+numerosity[i])],"extend")
  
  mydata$Rim_area[(kk+1):(kk+numerosity[i])]=na.approx(mydata$Rim_area[(kk+1):(kk+numerosity[i])],na.rm=FALSE)
  mydata$Rim_area[(kk+1):(kk+numerosity[i])]=na.fill(mydata$Rim_area[(kk+1):(kk+numerosity[i])],"extend")
  
  kk=kk+numerosity[i]
}
mydata$RNFL_average
mydata$Rim_area


### verifica osservazioni "anomale" RNFL_average

library(ggplot2)

p <- ggplot(data = mydata, aes(x = visit2, y = (RNFL_average-BaselineRNFL1)/BaselineRNFL1,  color=factor(Patient) )  )
p + geom_line() + geom_point() +
  geom_smooth(method = lm, color="black") +
  geom_hline(yintercept = -0.08) 

#evidentemente c'è un paziente "anomalo"; 
mydata$Patient[(mydata$RNFL_average-mydata$BaselineRNFL1)/mydata$BaselineRNFL1==max((mydata$RNFL_average-mydata$BaselineRNFL1)/mydata$BaselineRNFL1)]
#eliminiamo il paziente 501 per ora:
mydata=mydata[mydata$Patient!=501,]






save(mydata, file="../R_object/Glaucoma_better_data.RData")
