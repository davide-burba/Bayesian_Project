


setwd("~/Documents/universita/bayesian_statistics-Guglielmi/GLAUCOMA_PROJECT/src")

##################################################################
####################### LOADING DATA  ############################
##################################################################


Gdata=read.csv("../data/data.csv", header=T,fill=T,sep=",")
dim(Gdata)

attach(Gdata)
length(unique(Patient)) # 115 pazienti



##################################################################
##################### DATA VISUALIZATION  ########################
##################################################################
                           
library(ggplot2)

# patient
ggplot(data= Gdata, aes(factor(Patient), fill=1) ) + 
  geom_bar()+
  theme(legend.position = "none", axis.text.x  = element_text(angle=90, hjust=1, vjust=0.9)) 
#numero visite per paziente:
N=mean(table(Patient))
N
quantile(table(Patient),0.25)
quantile(table(Patient),0.5)
quantile(table(Patient),0.75)

#female vs male
ggplot(data= Gdata, aes(Sex, fill=Sex  ) ) + 
  geom_bar()+
  theme(legend.position = "none", axis.text.x  = element_text(angle=0, hjust=1, vjust=0.9)) 

# type of glaucoma
ggplot(data= Gdata, aes(Type_of_glaucoma, fill=Type_of_glaucoma  ) ) + 
  geom_bar()+
  theme(legend.position = "none", axis.text.x  = element_text(angle=0.45, hjust=1, vjust=0.9)) 



# sup_10thpercent
ggplot(data= Gdata, aes( sup_10thpercent )  ) + 
  geom_histogram(col="white", 
                 fill="blue",
                 binwidth = 10) +
  ylab("relative frequencies") + 
  scale_y_continuous(labels=scales::percent) +
  geom_bar()+
  theme(legend.position = "none", axis.text.x  = element_text(angle=0, hjust=1, vjust=0.9)) 


# sup_25thpercent
ggplot(data= Gdata, aes( sup_25thpercent )  ) + 
  geom_histogram(col="white", 
                 fill="blue",
                 binwidth = 25) +
  ylab("relative frequencies") + 
  scale_y_continuous(labels=scales::percent) +
  geom_bar()+
  theme(legend.position = "none", axis.text.x  = element_text(angle=0, hjust=1, vjust=0.9)) 

# Periph
ggplot(data= Gdata, aes(factor(Periph), fill=factor(Periph) ) ) + 
  geom_bar()+
  theme(legend.position = "none", axis.text.x  = element_text(angle=0, hjust=1, vjust=0.9)) 

###### Proviamo a visualizzare un unico paziente

Paz=590
Gdata_paz=Gdata[Patient==Paz,]
dim(Gdata_paz)

#qualche variabile: 
Gdata_paz$yearofglaucoma
Gdata_paz$Color
Gdata_paz$Periph
Gdata_paz$Visit
Gdata_paz$Progression
Gdata_paz$Progression1
Gdata_paz$ProgressionStructure1
Gdata_paz$ProgressionStructure
Gdata_paz
Gdata_paz$Age # variabile tempo!
Gdata_paz$age65
Gdata_paz$visualizationofON
Gdata_paz$MAP
Gdata_paz$FieldsComment
Gdata_paz$Field2
Gdata_paz$acuity



#MODIFICHE STRUTTURALI (riguardo al nervo ottico)
#RNFL
Gdata_paz$RNFL_average 
Gdata_paz$RNFL_cross_sectinal_area
Gdata_paz$RNFL_Multi_G_center
Gdata_paz$RNFL_Multi_Nasal
Gdata_paz$RNFL_Multi_Temp
Gdata_paz$RNFL_Multi_NI
Gdata_paz$RNFL_Multi_NS
Gdata_paz$RNFL_Multi_TI
Gdata_paz$RNFL_Multi_TS
Gdata_paz$RNFL_quad_Inf
Gdata_paz$RNFL_quad_nasal
Gdata_paz$RNFL_quad_Sup
Gdata_paz$RNFL_quad_temp
Gdata_paz$RNFL_thickness_superior
Gdata_paz$RNFL_thickness_inferior
Gdata_paz$RNFL_thickness_nasal
Gdata_paz$RNFL_thickness_temporal
Gdata_paz$OCT_RNFL_Signal_Strength
Gdata_paz$mean_RNFL_thickness
#RIM
Gdata_paz$Rim_area
Gdata_paz$RimArea
Gdata_paz$rim_volume



#MODIFICHE FUNZIONALI 
Gdata_paz$MD# Deviazione media CAMPO VISIVO!
Gdata_paz$PSD #Pattern Standard Deviation, nel perimetro Humphrey




# Evolution of PSD in time for a single patient
ggplot(data= Gdata_paz, aes(Age, PSD)) +
  geom_line(color = "black") 

# Evolution of RNFL_average in time for a single patient
ggplot(data= Gdata_paz, aes(Age, RNFL_average)) +
  geom_line(color = "black")

#Evaluation of MD
ggplot(data= Gdata_paz, aes(Age, MD)) +
  geom_line(color = "black")

#Let's see if the progression variables ; NB: progression and progressionstructure should be absorbing
ggplot(data= Gdata_paz, aes(Age, Progression)) +
  geom_line(color = "black")

ggplot(data= Gdata_paz, aes(Age, Progression1)) +
  geom_line(color = "black")

ggplot(data= Gdata_paz, aes(Age, ProgressionStructure)) +
  geom_line(color = "black")

ggplot(data= Gdata_paz, aes(Age), ProgressionStructure1)) +
  geom_line(color = "black")




#all the patient-> doesn't work
ggplot(data= Gdata, aes(Age, RNFL_average, color= Patient)) +
  geom_smooth(method = "loess", span = 1/2, se = FALSE)


#Plot of all th patients , RNFL
tmp=Gdata[Patient==502,]
minage=min(tmp$Age)
plot(x=tmp$Age-minage, y=tmp$RNFL_average, col=502, type="l",xlim=c(0,10), ylim=c(0,150))
for (i in factor(Patient))
{
  tmp=Gdata[Patient==i,]
  minage=min(tmp$Age)
  points(x=tmp$Age-minage, y=tmp$RNFL_average, col=i, type="l")
}


#Plot of all th patients , MD
tmp=Gdata[Patient==502,]
minage=min(tmp$Age)
plot(x=tmp$Age-minage, y=tmp$MD, col=502, type="l", ylim=c(-5,5))
for (i in factor(Patient))
{
  tmp=Gdata[Patient==i,]
  minage=min(tmp$Age)
  points(x=tmp$Age-minage, y=tmp$MD, col=i, type="l")
}


# TODO: 
#       Leggere papers e cerca HMM-> guarda che modello (0/1 o RNFL?)
#       modello
#       Riempire dati mancanti

 