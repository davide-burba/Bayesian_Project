#
# In questo script vengono "riempiti i dati NA di RNFL e RIM_AREA"
#

rm(list=ls())
load("../R_object/Glaucoma_data.RData")

Gdata_full=Gdata
rm(Gdata)
attach(Gdata_full)


# interpolazione con spline
library(zoo)

#prova:
na.spline(c(1,NA,12,3,NA,NA))


# riempi NA di RIM_AREA e RNFL_AVERAGE
numerosity=table(Patient)
kk=0
for(i in 1:length(unique(Patient))){
  RNFL_average[(kk+1):(kk+numerosity[i])]=na.spline(RNFL_average[(kk+1):(kk+numerosity[i])])
  Rim_area[(kk+1):(kk+numerosity[i])]=na.spline(Rim_area[(kk+1):(kk+numerosity[i])])
  kk=kk+numerosity[i]
}

save (Gdata_full, file="../R_object/Glaucoma_data_full.RData")
  


