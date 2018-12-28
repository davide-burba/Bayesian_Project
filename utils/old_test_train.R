#
# This file saves a training and a testing set; the training set correspond to all 
# the observation for each patient exept the last one which is included in the testing set 
#


rm(list=ls())

load("../R_object/Glaucoma_better_data.RData")

attach(mydata)
length(unique(Patient))

# Compute useful values 
numerosity<-as.integer(table(Patient))



test_set= mydata[cumsum(numerosity),]
test_set$visit2


train_set= mydata[-cumsum(numerosity),]
train_set$visit2

save(train_set, file="../R_object/train_set.RData")
save(test_set, file="../R_object/test_set.RData")

