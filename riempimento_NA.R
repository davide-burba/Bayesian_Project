#
# In questo script vengono "riempiti i dati NA di RNFL_average, mean_RNFL_thikness e RIM_AREA"
#

rm(list=ls())
load("../R_object/Glaucoma_data.RData")

Gdata_full=Gdata
rm(Gdata)



# interpolazione con spline
library(zoo)

#prova:
na.spline(c(1,NA,12,3,NA,NA))

# paziente "strano":
Gdata_full$mean_RNFL_thickness[which(Gdata_full$Patient==582)]
#Tutti NA nel paziente 582, togliamolo
Gdata_full=Gdata_full[Patient!=582,]




# riempi NA di RIM_AREA, mean_RNFL_thikness e RNFL_AVERAGE
numerosity=table(Gdata_full$Patient)
kk=0
for(i in 1:length(unique(Gdata_full$Patient))){
  Gdata_full$RNFL_average[(kk+1):(kk+numerosity[i])]=na.spline(Gdata_full$RNFL_average[(kk+1):(kk+numerosity[i])])
  Gdata_full$Rim_area[(kk+1):(kk+numerosity[i])]=na.spline(Gdata_full$Rim_area[(kk+1):(kk+numerosity[i])])
  Gdata_full$Rim_area[(kk+1):(kk+numerosity[i])]=na.spline(Gdata_full$Rim_area[(kk+1):(kk+numerosity[i])])
  Gdata_full$mean_RNFL_thickness[(kk+1):(kk+numerosity[i])]=na.spline(Gdata_full$mean_RNFL_thickness[(kk+1):(kk+numerosity[i])])
  
  kk=kk+numerosity[i]
}

save (Gdata_full, file="../R_object/Glaucoma_data_full.RData")
  


