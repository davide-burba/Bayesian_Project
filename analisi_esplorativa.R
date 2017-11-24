


setwd("~/Documents/universita/bayesian_statistics-Guglielmi/GLAUCOMA_PROJECT/src")

##################################################################
####################### LOADING DATA  ############################
##################################################################


Gdata=read.table("../data/data.csv", header=T,fill=T,sep=",")
dim(Gdata)

attach(Gdata)
length(unique(Patient)) # 55 pazienti
length(Patient[Patient==0]) # c'è un' osservazione strana, probabile errore
Gdata[ Patient==0,] # quasi tutto NA -> rimuoviamola per ora
Gdata=Gdata[ Patient!=0,]


##################################################################
##################### DATA VISUALIZATION  ########################
##################################################################
                           
library(ggplot2)

# patient
ggplot(data= Gdata, aes(factor(Patient), fill=factor(Patient) ) ) + 
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
Gdata_paz$Age # variabile tempo!
Gdata_paz$age65
Gdata_paz$visualizationofON
Gdata_paz$MAP
Gdata_paz$FieldsComment
Gdata_paz$Field2
Gdata_paz$acuity
Gdata_paz$RNFL_average# questo è importante, riguarda modifiche strutturali del nervo 
Gdata_paz$MD
Gdata_paz$PSD# questo riguarda il campo visivo



# Evolution of PSD in time for a single patient
ggplot(data= Gdata_paz, aes(Age, PSD)) +
  geom_line(color = "black") 

# Evolution of RNFL_average in time for a single patient
ggplot(data= Gdata_paz, aes(Age, RNFL_average)) +
  geom_line(color = "black")

#Let's see if the progression variables mean absorbing states:
ggplot(data= Gdata_paz, aes(Age, Progression)) +
  geom_line(color = "black")

ggplot(data= Gdata_paz, aes(Age, Progression1)) +
  geom_line(color = "black")

ggplot(data= Gdata_paz, aes(Age, ProgressionStructure)) +
  geom_line(color = "black")

ggplot(data= Gdata_paz, aes(Age, ProgressionStructure1)) +
  geom_line(color = "black")

#not always consistent


 