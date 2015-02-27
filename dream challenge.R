library(caret)
library(randomForest)
setwd("~/Documents/wendy/R/DREAM olfaction")
#read in training files
train.per<-read.table("TrainSet.txt",header=T,sep="\t",comment.char="",strip.white=T)
train.per$Odor<-gsub( " ", "", train.per$Odor)
des<-read.table("molecular_descriptors_data.txt",header=T,sep="\t",comment.char="")

#get molecular descriptors for molecules in training set
CID<-unique(train.per$Compound.Identifier)
train.des<-des[match(CID,des$CID),]

###perception data
#pull out replicates and their original ratings
rep<-train.per[train.per$Replicate=="replicate",]
Odor<-as.character(unique(rep$Odor))
Odor<-strsplit(Odor, "\\(r")
repOdor<-as.character(NA)
for(i in 1:20){
  repOdor<-c(repOdor,Odor[[i]][1])
}
repOdor<-repOdor[-1]
repOdor<-gsub( " ", "", repOdor)
match<-data.frame(match(as.character(train.per$Odor),repOdor),row=c(1:nrow(train.per)))
match<-data.frame(match[complete.cases(match[,1]),])
original<-train.per[match$row,]

#average the scores of originals and replicates (treat NAs as 0)
rep[is.na(rep)]<-0 #replace NaN with 0
original[is.na(original)]<-0 #replace NaN with 0
avg<-cbind(original[,1:6],(original[,7:27]+rep[,7:27])/2)

#replace original and duplicates with the averaged score - total of 338 molecules
train.per2<-rbind(avg,train.per[-as.numeric(c(row.names(original),row.names(rep))),])

#seperating high and low intensity 
high<-subset(train.per2,Intensity=="high")
low<-subset(train.per2,Intensity=="low")

#calculate avg scores for each molecule
train.avgH<-data.frame(matrix(ncol=23))
colnames(train.avgH)<-colnames(train.per)[-c(3:6)]
for(i in 1:length(CID)){
  sub<-subset(high,Compound.Identifier==CID[i])
  train.avgH[i,]<-c(sub[1,1:2],colMeans(sub[,7:27],na.rm=T))
}

train.avgL<-data.frame(matrix(ncol=23))
colnames(train.avgL)<-colnames(train.per)[-c(3:6)]
for(i in 1:length(CID)){
  sub<-subset(low,Compound.Identifier==CID[i])
  train.avgL[i,]<-c(sub[1,1:2],colMeans(sub[,7:27],na.rm=T))
}


###Molecular data
#drop molecules with more than 1000 NAs
nacount<-data.frame(is.na(train.des))
train.des2<-train.des[which(rowSums(nacount)<1000),]
#drop predictors with any NA
nacount<-data.frame(is.na(train.des2))
train.des2<-train.des2[,which(colSums(nacount)==0)]
#drop near zero variables 
nzv.num<-nearZeroVar(train.des2[,-1])
train.des3<-train.des2[,-(nzv.num+1)]


###Build model - randomForest
#High intensity
tb<-cbind(Intensity=train.avgH$INTENSITY.STRENGTH[match(train.des3$CID,train.avgH$Compound.Identifier)], train.des2[,-1]) #nzv removed 
set.seed(500)
modfit<-randomForest(y=as.factor(tb$Intensity), x=tb[,2:length(tb)], ntree=500, importance=TRUE, proximity=TRUE, na.action=na.omit)

###Leaderboard prediction
lbCID<-read.table("CID_leaderboard.txt",header=F,sep="\t")
lb<-des[match(lbCID[,1],des$CID),]
pred<-predict(modfit,lb,type="response")

#Baseline prediction
base<-read.table("CecchiG_DREAM95olf_s2.txt",header=F,sep="\t")

#correlation
cor(pred,as.numeric(base[1:69,3]))

##Need to do more models 