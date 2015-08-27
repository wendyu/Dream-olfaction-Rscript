###load libraries
library(caret)
library(randomForest)
library(plyr)
library(Boruta)
library(ggplot2)
library(reshape2)
#Set working directory
setwd("S:/mainland/Users/Wendy Yu/DREAM olfaction")

###1) Data clean up 
#A) Clean up molecular descriptor 
des<-read.table("molecular_descriptors_data.txt",header=T,sep="\t",comment.char="")
#drop near zero variables
nzv.num<-nearZeroVar(des[,-1])
des2<-des[,-(nzv.num+1)]
#drop predictors with any NA
nacount<-data.frame(is.na(des2))
des3<-des2[,which(colSums(nacount)==0)]
#drop highly correlated variables 
desCor<-cor(des3[,-1])
highlyCorDes<-findCorrelation(desCor, cutoff=0.99)
des4<-des3[,-(highlyCorDes+1)]

#B) Clean up molecular fringerprint 
setwd("S:/mainland/Users/Wendy Yu/DREAM olfaction/Features")#change working directory
finger<-read.table("DREAMset paired fingerprint.txt", header=T,sep="\t", comment.char="")
gsfinger<-read.table("DREAMset vs GoodScent fingerprint.txt",header=T,sep="\t",comment.char="")
colnames(gsfinger)<-paste0("GS",colnames(gsfinger))
finger<-cbind(finger,gsfinger[,-1])
#drop near zero variables
nzv.num<-nearZeroVar(finger[,-1])
finger2<-finger[,-(nzv.num+1)]
#drop highly correlated variables 
fingerCor<-cor(finger2[,-1])
highlyCorDes<-findCorrelation(fingerCor, cutoff=0.99)
finger3<-finger2[,-(highlyCorDes+1)]

#C) Clean up EPI Suite
epi<-read.table("DREAM_EPI.txt",header=T,comment.char="",sep="\t")
#drop near zero variables
nzv.num<-nearZeroVar(epi[,-c(1:4)])
epi2<-epi[,-c(1:3,nzv.num+4)]
#drop predictors with any NA
nacount<-data.frame(is.na(epi2))
epi3<-epi2[,which(colSums(nacount)==0)]
#drop predictors is not numeric 
#(there are some text columns recording the date and time when the file is generated)
epi4<-epi3[,sapply(epi3,is.numeric)]
#drop highly correlated variables
epiCor<-cor(epi4)
highlyCorDes<-findCorrelation(epiCor, cutoff=0.99)
epi5<-epi4[,-(highlyCorDes)]

#D) Clean up perception/phenotype data
train.per<-read.table("TrainSet.txt",header=T,sep="\t",comment.char="",strip.white=T)
train.per$Odor<-gsub( " ", "", train.per$Odor)
levels(train.per$Dilution)
train.per$Dilution<-gsub("1/10,000,000 ",0.0000001,train.per$Dilution) #notice the spac. 
train.per$Dilution<-gsub("1/100,000 ",0.00001,train.per$Dilution)
train.per$Dilution<-gsub("1/1,000 ",0.001,train.per$Dilution)
train.per$Dilution<-gsub("1/10",0.1,train.per$Dilution)
train.per$Dilution<-as.numeric(train.per$Dilution)
colnames(train.per)[1]<-c("CID")
##average preception scores based on dilution.
uCID<-unique(train.per$CID)
train.per2<-data.frame(matrix(nrow=42))
for( i in 1:length(uCID)){
  sub<-subset(train.per,CID==uCID[i])
  uDil<-unique(sub$Dilution)
  for( j in 1:2){
    mean<-data.frame(colMeans(sub[sub$Dilution==uDil[j],c(4:length(sub))],na.rm=T))
    colnames(mean)<-paste0(uCID[i],"_",uDil[j])
    row.names(mean)<-paste0(row.names(mean),".Mean")
    sd<-data.frame(apply(sub[sub$Dilution==uDil[j],c(4:length(sub))],2,sd, na.rm=T))
    colnames(sd)<-paste0(uCID[i],"_",uDil[j])
    row.names(sd)<-paste0(row.names(sd),".sd")
    train.per2<-cbind(train.per2, rbind(mean,sd))}
}
train.per2<-data.frame(t(train.per2[,-1]))
split<-data.frame(t(data.frame(strsplit(row.names(train.per2),"_"))))
train.per2<-cbind(CID=split[,1],Dilution=as.numeric(as.character(split[,2])),train.per2)
train.per2<-train.per2[complete.cases(train.per2),]

###2) prepare data for modeling 
#Models performed better without removing the correlated varaibles 
insert.dragon=des3 
insert.finger=finger2
insert.epi=epi4
trainSet<-cbind(insert.dragon,insert.finger[,-1],insert.epi)
#Standardizing
preObj<-preProcess(trainSet[,-1],method=c("center","scale"))
trainSetCS<-cbind(CID=trainSet[,1],predict(preObj,trainSet[,-1]))
#Combine with perception data.
trainSet2<-merge(train.per2,trainSetCS,by="CID",all.x=T)

###3) feature selection - Recursive feature elimination using RF function
#define rfe control
control <- rfeControl(functions=rfFuncs, method="repeatedcv", number=5,repeats=1)
subsets <- c(1:5, 10, 15, 20, 25,50) #select the number of variables. 
#select variables for prediction on mean intensity. 
tb=trainSet2[,c(2,3,45:length(trainSet2))]
rfRFE<-rfe(tb[,c(1,3:length(tb))],tb[,2],sizes=subsets,rfeControl=control)
rfRFE.var<-varImp(rfRFE)  

###4) Building models and cross validation
#define training control
fitControl <- trainControl(method = "repeatedcv",number = 5,repeats = 0, #5 fold cv, repeat 5 times
                           savePredictions=FALSE) 
#tuning parameter
gbmGrid <-  expand.grid(interaction.depth = c(1, 5, 9),
                        n.trees = (1:4)*50,
                        shrinkage = 0.1)
#set random seed
set.seed(825)
#insert the selected variable table
selectedvar=trainSet2[,c("INTENSITY.STRENGTH.Mean",row.names(rfRFE.var))] #use the variables selected from feature selection
#Gradient boosting selected variables - really fast. good results.  
gbmFit<-train(INTENSITY.STRENGTH.Mean~.,data=selectedvar,
               method = "gbm",
               trControl = fitControl,
               verbose = FALSE,
               tuneGrid = gbmGrid)

###5) External validataion on leaderboard set 
#make prediction on mean intensity in 1/1000 dilution 
lbAns<-read.table("LBs2.txt",header=F)
lbCID<-read.table("CID_leaderboard.txt",header=F,sep="\t")
lbdil<-read.table("dilution_leaderboard.txt",header=T,sep="\t",quote="\"")
lbdil$dilution<-gsub(" '1/100,000'",0.00001,lbdil$dilution)
lbdil$dilution<-gsub(" '1/10'" ,0.1,lbdil$dilution)
lbdil$dilution<-as.numeric(lbdil$dilution)
lbdil$dilution1000<-as.numeric(0.001)
lb.des<-trainSetCS[trainSetCS$CID %in% lbCID[,1],]
lb.des$Dilution<-lbdil$dilution1000 
#make prediction - enter which model
insert.mod= gbmFit
selecteddt<-lb.des[,c(row.names(rfRFE.var))]
pred<-predict(insert.mod,newdata=selecteddt)
cor(as.numeric(as.character(pred)),as.numeric(lbAns[1:69,3]))




