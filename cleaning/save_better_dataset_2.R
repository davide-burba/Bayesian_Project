#
# This file is used to save a more easy-to-use and "correct" dataset; 
# uguale a "save_better_dataset.R" ma salva pazienti con almeno 6 visite

rm(list=ls())

mydata_6=read.csv("../data/data.csv", header=T,fill=T,sep=",")


# Consider only patients with at least 6 visits:
tmp=unique(mydata_6$Patient)
a<-table(mydata_6$Patient)
tmp_data=tmp[which(a>=6)] 
mydata_6=mydata_6[mydata_6$Patient %in% tmp_data,]
rm(tmp)
rm(tmp_data)
rm(a)


# Modify Sex covariate: Male=1, Female=0
tmp = factor(mydata_6$Sex, levels=c("Male","Female"), labels=c(1,0))
mydata_6$Sex=as.numeric(as.character(tmp))
rm(tmp)
mydata_6$Sex


# Race: build 4 vectors of boolean variables Asian Black Hispanic White 
attach(mydata_6)
Asian=(Race=="Asian")*1
Black=(Race=="Black")*1
Hispanic=(Race=="Hispanic")*1
White=(Race=="White")*1
detach(mydata_6)
mydata_6=cbind(mydata_6,Asian,Black,Hispanic,White)
mydata_6$Asian


#  Type_of_glaucoma: identify patients with POAG glaucoma (many NA)
attach(mydata_6)
POAG=(Type_of_glaucoma=="POAG")*1
detach(mydata_6)
mydata_6=cbind(mydata_6,POAG)
mydata_6$POAG


#  Allergies: trascuriamolo per ora
unique(mydata_6$Allergies)


#  Modify Hypertension: Yes=1, No=0
tmp = factor(mydata_6$Hypertension, levels=c("Yes","No"), labels=c(1,0))
mydata_6$Hypertension=as.numeric(as.character(tmp))
rm(tmp)
mydata_6$Hypertension


#  Modify HyperLipidemia: Yes=1, No=0
tmp = factor(mydata_6$HyperLipidemia, levels=c("Yes","No"), labels=c(1,0))
mydata_6$HyperLipidemia=as.numeric(as.character(tmp))
rm(tmp)
mydata_6$HyperLipidemia


#  Modify Diabetes: Yes=1, No=0  (from now on in one row)
mydata_6$Diabetes=as.numeric(as.character(factor(mydata_6$Diabetes, levels=c("Yes","No"), labels=c(1,0))))
mydata_6$Diabetes


#  Modify Cardiovascular_Dz: Yes=1, No=0  (now in one row)
mydata_6$Cardiovascular_Dz=as.numeric(as.character(factor(mydata_6$Cardiovascular_Dz, levels=c("Yes","No"), labels=c(1,0))))
mydata_6$Cardiovascular_Dz


# Other_systemic_dz_1: trascuriamolo per ora
unique(mydata_6$Other_systemic_dz_1)


# family history: since there are many NA, we made a vector with  1=yes, 0  otherwise
attach(mydata_6)
familiarity_yes=(family_history=="yes")*1
familiarity_yes
detach(mydata_6)
mydata_6=cbind(mydata_6,familiarity_yes)
mydata_6$familiarity_yes


# Smoking: Yes=1, No=0  
mydata_6$Smoking=as.numeric(as.character(factor(mydata_6$Smoking, levels=c("Yes","No"), labels=c(1,0))))
mydata_6$Smoking


#Meds_according_to_Pt: many different; for the moment 1="medicine yes", 0="medicine no"
attach(mydata_6)
Meds_according_to_Pt_bool=(Meds_according_to_Pt==""|Meds_according_to_Pt=="none" )*1
Meds_according_to_Pt_bool
detach(mydata_6)
mydata_6=cbind(mydata_6,Meds_according_to_Pt_bool)
mydata_6$Meds_according_to_Pt_bool


#Ocular_Meds_according_to_Pt: many different; for the moment 1=" ocular medicine yes", 0="ocul. medicine no"
attach(mydata_6)
Ocular_Meds_according_to_Pt_bool=(Ocular_Meds_according_to_Pt==""|Ocular_Meds_according_to_Pt=="none" )*1
Ocular_Meds_according_to_Pt_bool
detach(mydata_6)
mydata_6=cbind(mydata_6,Ocular_Meds_according_to_Pt_bool)
mydata_6$Ocular_Meds_according_to_Pt_bool


# BASELINE RNFL (RNFL AT FIRST VISIT)
attach(mydata_6)
BaselineRNFL1=rep(0,1050)
for (i in 1:dim(mydata_6)[1])
{
  if (visit2[i]==0)
  {
    tmp=RNFL_average[i] 
  }
  BaselineRNFL1[i]=tmp
}
detach(mydata_6)
mydata_6=cbind(mydata_6,BaselineRNFL1)
mydata_6$BaselineRNFL1



##
## RIEMPIMENTO NA: interpolazione lineare più stesse osservazioni alle code
##

library(zoo)

#prova:
tmp=na.approx(c(NA,1,NA,12,3,NA,76,NA,NA),na.rm=FALSE)
tmp
na.fill(tmp,"extend")
rm(tmp)

# riempi NA di:  RNFL_AVERAGE, RIM_AREA ,IOP, acuity, MD, PSD,  cup_disk_horiz_ratio, cup_disk_vert_ratio,
#                macular_volume, Vert_integrated_rim_area__vol_, Horz_integrated_rim_width__area_,   
#                Cup_area, cup_disk_area_ratio, cup_disk_horiz_ratio, cup_disk_vert_ratio


numerosity=table(mydata_6$Patient)
kk=0
for(i in 1:length(unique(mydata_6$Patient))){
  mydata_6$RNFL_average[(kk+1):(kk+numerosity[i])]=na.approx(mydata_6$RNFL_average[(kk+1):(kk+numerosity[i])],na.rm=FALSE)
  mydata_6$RNFL_average[(kk+1):(kk+numerosity[i])]=na.fill(mydata_6$RNFL_average[(kk+1):(kk+numerosity[i])],"extend")
  
  mydata_6$Rim_area[(kk+1):(kk+numerosity[i])]=na.approx(mydata_6$Rim_area[(kk+1):(kk+numerosity[i])],na.rm=FALSE)
  mydata_6$Rim_area[(kk+1):(kk+numerosity[i])]=na.fill(mydata_6$Rim_area[(kk+1):(kk+numerosity[i])],"extend")
  
  mydata_6$IOP[(kk+1):(kk+numerosity[i])]=na.approx(mydata_6$IOP[(kk+1):(kk+numerosity[i])],na.rm=FALSE)
  mydata_6$IOP[(kk+1):(kk+numerosity[i])]=na.fill(mydata_6$IOP[(kk+1):(kk+numerosity[i])],"extend")
  
  mydata_6$acuity[(kk+1):(kk+numerosity[i])]=na.approx(mydata_6$acuity[(kk+1):(kk+numerosity[i])],na.rm=FALSE)
  mydata_6$acuity[(kk+1):(kk+numerosity[i])]=na.fill(mydata_6$acuity[(kk+1):(kk+numerosity[i])],"extend")
  
  mydata_6$MD[(kk+1):(kk+numerosity[i])]=na.approx(mydata_6$MD[(kk+1):(kk+numerosity[i])],na.rm=FALSE)
  mydata_6$MD[(kk+1):(kk+numerosity[i])]=na.fill(mydata_6$MD[(kk+1):(kk+numerosity[i])],"extend")
  
  mydata_6$PSD[(kk+1):(kk+numerosity[i])]=na.approx(mydata_6$PSD[(kk+1):(kk+numerosity[i])],na.rm=FALSE)
  mydata_6$PSD[(kk+1):(kk+numerosity[i])]=na.fill(mydata_6$PSD[(kk+1):(kk+numerosity[i])],"extend")
  
  mydata_6$cup_disk_horiz_ratio[(kk+1):(kk+numerosity[i])]=na.approx(mydata_6$cup_disk_horiz_ratio[(kk+1):(kk+numerosity[i])],na.rm=FALSE)
  mydata_6$cup_disk_horiz_ratio[(kk+1):(kk+numerosity[i])]=na.fill(mydata_6$cup_disk_horiz_ratio[(kk+1):(kk+numerosity[i])],"extend")
  
  mydata_6$cup_disk_vert_ratio[(kk+1):(kk+numerosity[i])]=na.approx(mydata_6$cup_disk_vert_ratio[(kk+1):(kk+numerosity[i])],na.rm=FALSE)
  mydata_6$cup_disk_vert_ratio[(kk+1):(kk+numerosity[i])]=na.fill(mydata_6$cup_disk_vert_ratio[(kk+1):(kk+numerosity[i])],"extend")
  
  mydata_6$macular_volume[(kk+1):(kk+numerosity[i])]=na.approx(mydata_6$macular_volume[(kk+1):(kk+numerosity[i])],na.rm=FALSE)
  mydata_6$macular_volume[(kk+1):(kk+numerosity[i])]=na.fill(mydata_6$macular_volume[(kk+1):(kk+numerosity[i])],"extend")
  
  mydata_6$Vert_integrated_rim_area__vol_[(kk+1):(kk+numerosity[i])]=na.approx(mydata_6$Vert_integrated_rim_area__vol_[(kk+1):(kk+numerosity[i])],na.rm=FALSE)
  mydata_6$Vert_integrated_rim_area__vol_[(kk+1):(kk+numerosity[i])]=na.fill(mydata_6$Vert_integrated_rim_area__vol_[(kk+1):(kk+numerosity[i])],"extend")
  
  mydata_6$Horz_integrated_rim_width__area_[(kk+1):(kk+numerosity[i])]=na.approx(mydata_6$Horz_integrated_rim_width__area_[(kk+1):(kk+numerosity[i])],na.rm=FALSE)
  mydata_6$Horz_integrated_rim_width__area_[(kk+1):(kk+numerosity[i])]=na.fill(mydata_6$Horz_integrated_rim_width__area_[(kk+1):(kk+numerosity[i])],"extend")
  
  mydata_6$Cup_area[(kk+1):(kk+numerosity[i])]=na.approx(mydata_6$Cup_area[(kk+1):(kk+numerosity[i])],na.rm=FALSE)
  mydata_6$Cup_area[(kk+1):(kk+numerosity[i])]=na.fill(mydata_6$Cup_area[(kk+1):(kk+numerosity[i])],"extend")
  
  mydata_6$cup_disk_area_ratio[(kk+1):(kk+numerosity[i])]=na.approx(mydata_6$cup_disk_area_ratio[(kk+1):(kk+numerosity[i])],na.rm=FALSE)
  mydata_6$cup_disk_area_ratio[(kk+1):(kk+numerosity[i])]=na.fill(mydata_6$cup_disk_area_ratio[(kk+1):(kk+numerosity[i])],"extend")
  
  mydata_6$cup_disk_horiz_ratio[(kk+1):(kk+numerosity[i])]=na.approx(mydata_6$cup_disk_horiz_ratio[(kk+1):(kk+numerosity[i])],na.rm=FALSE)
  mydata_6$cup_disk_horiz_ratio[(kk+1):(kk+numerosity[i])]=na.fill(mydata_6$cup_disk_horiz_ratio[(kk+1):(kk+numerosity[i])],"extend")
  
  mydata_6$cup_disk_vert_ratio[(kk+1):(kk+numerosity[i])]=na.approx(mydata_6$cup_disk_vert_ratio[(kk+1):(kk+numerosity[i])],na.rm=FALSE)
  mydata_6$cup_disk_vert_ratio[(kk+1):(kk+numerosity[i])]=na.fill(mydata_6$cup_disk_vert_ratio[(kk+1):(kk+numerosity[i])],"extend")
  
  mydata_6$AGIS_score[(kk+1):(kk+numerosity[i])]=na.approx(mydata_6$AGIS_score[(kk+1):(kk+numerosity[i])],na.rm=FALSE)
  mydata_6$AGIS_score[(kk+1):(kk+numerosity[i])]=na.fill(mydata_6$AGIS_score[(kk+1):(kk+numerosity[i])],"extend")
  
  
  kk=kk+numerosity[i]
}
mydata_6$RNFL_average
mydata_6$Rim_area
mydata_6$cup_disk_vert_ratio

### verifica osservazioni "anomale" RNFL_average

library(ggplot2)

p <- ggplot(data = mydata_6, aes(x = visit2, y = (RNFL_average-BaselineRNFL1)/BaselineRNFL1,  color=factor(Patient) )  )
p + geom_line() + geom_point() +
  #geom_smooth(method = lm, color="black") +
  #geom_hline(yintercept = -0.08)  +
  labs(  y="(RNFL_average-baseline)/baseline", x = "Time (Years)")+ theme(legend.position="none")


#ggsave('outliar.png', width = 8, height = 4)


#evidentemente c'è un paziente "anomalo"; 
mydata_6$Patient[(mydata_6$RNFL_average-mydata_6$BaselineRNFL1)/mydata_6$BaselineRNFL1==max((mydata_6$RNFL_average-mydata_6$BaselineRNFL1)/mydata_6$BaselineRNFL1)]
#eliminiamo il paziente 501 per ora:
mydata_6=mydata_6[mydata_6$Patient!=501,]



# aggiungo variabile log di RNFL e scale di RNFL

RNFL_average_scaled  =as.vector(scale(mydata_6$RNFL_average))
RNFL_log= log(mydata_6$RNFL_average)

plot(RNFL_log)
plot(RNFL_average_scaled)

mydata_6=cbind(mydata_6,RNFL_average_scaled,RNFL_log)


# standardizzo variabili
mydata_6$RNFL_average  =as.vector(scale(mydata_6$RNFL_average))
mydata_6$Age  =as.vector(scale(mydata_6$Age))
mydata_6$IOP  =as.vector(scale(mydata_6$IOP))
mydata_6$acuity  =as.vector(scale(mydata_6$acuity))
mydata_6$MD  =as.vector(scale(mydata_6$MD))
mydata_6$AGIS_score  =as.vector(scale(mydata_6$AGIS_score))
mydata_6$PSD  =as.vector(scale(mydata_6$PSD))
mydata_6$cup_disk_horiz_ratio  =as.vector(scale(mydata_6$cup_disk_horiz_ratio))
mydata_6$cup_disk_vert_ratio  =as.vector(scale(mydata_6$cup_disk_vert_ratio))
mydata_6$macular_volume  =as.vector(scale(mydata_6$macular_volume))
mydata_6$Vert_integrated_rim_area__vol_  =as.vector(scale(mydata_6$Vert_integrated_rim_area__vol_))
mydata_6$Horz_integrated_rim_width__area_  =as.vector(scale(mydata_6$Horz_integrated_rim_width__area_))
mydata_6$Cup_area     =as.vector(scale(mydata_6$Cup_area))               
mydata_6$Rim_area  =as.vector(scale(mydata_6$Rim_area))
mydata_6$cup_disk_area_ratio        =as.vector(scale(mydata_6$cup_disk_area_ratio))      
mydata_6$yearofglaucoma  =as.vector(scale(mydata_6$yearofglaucoma))


length(unique(mydata_6$Patient))

save(mydata_6, file="../R_object/Glaucoma_over_6_obs.RData")
