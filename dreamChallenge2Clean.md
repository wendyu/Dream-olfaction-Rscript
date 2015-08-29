Dream Olfaction SubChallenge 2
========================================================
The goal of this script is to build a model that predict the olfactory intensity perception from molecular structures. 
In this script I used the Caret R package to perform machine learning. 


```r
###load libraries
library(caret)
```

```
## Loading required package: lattice
## Loading required package: ggplot2
```

```r
library(randomForest)
```

```
## randomForest 4.6-10
## Type rfNews() to see new features/changes/bug fixes.
```

```r
library(plyr)
library(Boruta)
```

```
## Loading required package: rFerns
```

```r
library(ggplot2)
library(reshape2)
#Set working directory
setwd("S:/mainland/Users/Wendy Yu/DREAM olfaction")
```

Here we have three differnet files that give the physicochemical properties of molecules. Data can contain attributes that are not very informative, for example the zero variance variables that have uniform values. Data can also contain attributes that are highly correlated with each other. Caret package provide methods to remove redundant and non-imformative attributes. 


```r
###1) Data clean up 
#A) Clean up dragon molecular descriptor 
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
```

For the olfactory perception data, every compound is tested in high and low concentraions and each compound is rated twice. Here we averaged the intensity rating according to its diluted concentration. 

```r
#D) Clean up perception/phenotype data
setwd("S:/mainland/Users/Wendy Yu/DREAM olfaction")
train.per<-read.table("TrainSet.txt",header=T,sep="\t",comment.char="",strip.white=T)
train.per$Odor<-gsub( " ", "", train.per$Odor)
#levels(train.per$Dilution)
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
    mean<-data.frame(colMeans(sub[sub$Dilution==uDil[j],c(7:length(sub))],na.rm=T))
    colnames(mean)<-paste0(uCID[i],"_",uDil[j])
    row.names(mean)<-paste0(row.names(mean),".Mean")
    sd<-data.frame(apply(sub[sub$Dilution==uDil[j],c(7:length(sub))],2,sd, na.rm=T))
    colnames(sd)<-paste0(uCID[i],"_",uDil[j])
    row.names(sd)<-paste0(row.names(sd),".sd")
    train.per2<-cbind(train.per2, rbind(mean,sd))}
}
train.per2<-data.frame(t(train.per2[,-1]))
split<-data.frame(t(data.frame(strsplit(row.names(train.per2),"_"))))
train.per2<-cbind(CID=split[,1],Dilution=as.numeric(as.character(split[,2])),train.per2)
train.per2<-train.per2[complete.cases(train.per2),]
```

Next we want to combine the descriptor data and the perception data. After some trial and error we found out that removing correlated attributes decreases the model performance. Here we decide to use the molecular decriptors after removing near-zero-variance attributes but keep the correlated attributes. We then standarized all the descriptors. 

```r
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
```

Selecting the right features is a very importance step in machine learning. Here we use the recursive feature elimination algorithm with random forest model to perform feature selection. 

```r
###3) feature selection - Recursive feature elimination using RF function
#define rfe control
control <- rfeControl(functions=rfFuncs, method="repeatedcv", number=5,repeats=1)
subsets <- c(1:5, 10, 25, 50, 100, 200, 500, 1000, 2000) #select the number of variables. 
#select variables for prediction on mean intensity. 
tb=trainSet2[,c(2,3,45:length(trainSet2))]
rfRFE<-rfe(tb[,c(1,3:length(tb))],tb[,2],sizes=subsets,rfeControl=control)
rfRFE.var<-varImp(rfRFE)  
```


```r
plot(rfRFE, type=c("g","o"))
```

![plot of chunk unnamed-chunk-6](figure/unnamed-chunk-6-1.png) 

After selecting the features, we want to build a machine learning model using those features. We chose gradient boosting model and performed a 5 fold cross validation to exam the model preformance. 

```r
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
               tuneGrid = gbmGrid,
               metric="Rsquared")
```

```
## Loading required package: gbm
## Loading required package: survival
## 
## Attaching package: 'survival'
## 
## The following object is masked from 'package:caret':
## 
##     cluster
## 
## Loading required package: splines
## Loading required package: parallel
## Loaded gbm 2.1.1
```


```r
#plot results
plot(gbmFit, metric="Rsquared")
```

![plot of chunk unnamed-chunk-8](figure/unnamed-chunk-8-1.png) 

To perform external model validation, we applied the gradient boosting model to the leaderboard dataset. 

```r
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
cor(as.numeric(as.character(pred)),as.numeric(lbAns[1:69,3])) #correlation coefficient
```

```
## [1] 0.5627804
```

```r
plot<-data.frame(pred=pred,ans=lbAns[1:69,3])
```


```r
ggplot(plot, aes(x=pred,y=ans))+geom_point()+stat_smooth(method="lm")+
  xlab("Model prediction")+
  ylab("Groud truth")+
  ggtitle("Model prediction on leaderboard data")
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-1.png) 
