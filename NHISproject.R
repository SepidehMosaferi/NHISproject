# ===========================================================================================
# Programming for the NHIS project
# Author: Sepideh Mosaferi
# Date: Sep. 2016
# ============================================================================================

# Required Libraries -------------------------------------------------------------------------

require(PracTools)
require(sampling)
require(reshape)
require(survey)
require(lme4)
require(QRM)
require(reshape2)
require(mi)
require(mice)
require(VIM)

# Preparation of Data -------------------------------------------------------------------------

setwd("C://Users//smosafer//Projects")
getwd()
NHISdata <- read.csv("samadult.csv")
attach(NHISdata)
NHISdata <- NHISdata[order(STRAT_P,PSU_P,HHX,FMX,FPX),] #Sorting the dataset hierarchically

NHISdata$ID <- 1:nrow(NHISdata) #Assigning Unique ID to persons
round(nrow(NHISdata)/15)
NHISdata$STRATUM <- c(rep(1:14,each=round(nrow(NHISdata)/15)),rep(15,nrow(NHISdata)-14*round(nrow(NHISdata)/15)))

set.seed(123)
STR1BG <- rnbinom(7,300,0.5)+c(rep(18,5),rep(19,2))
set.seed(456)
STR2BG <- rnbinom(7,300,0.5)+c(rep(26,5),rep(25,2))
set.seed(789)
STR3BG <- rnbinom(7,300,0.5)+c(rep(42,4),rep(41,3))
set.seed(101112)
STR4BG <- rnbinom(7,300,0.5)+c(rep(27,6),rep(26,1))
set.seed(131415)
STR5BG <- rnbinom(7,300,0.5)+28
set.seed(161718)
STR6BG <- rnbinom(7,300,0.5)+c(rep(16,5),rep(15,2))
set.seed(192021)
STR7BG <- rnbinom(7,300,0.5)+c(rep(22,6),23)
set.seed(222324)
STR8BG <- rnbinom(7,300,0.5)+c(rep(11,5),rep(10,2))
set.seed(252627)
STR9BG <- rnbinom(7,300,0.5)+c(rep(21,3),rep(20,4))
set.seed(282930)
STR10BG <- rnbinom(7,300,0.5)+c(rep(12,3),rep(11,4))
set.seed(313233)
STR11BG <- rnbinom(7,300,0.5)+c(rep(20,3),rep(19,4))
set.seed(343536)
STR12BG <- rnbinom(7,300,0.5)+c(rep(20,4),rep(19,3))
set.seed(373839)
STR13BG <- rnbinom(7,300,0.5)+c(rep(19,5),rep(20,2))
set.seed(404142)
STR14BG <- rnbinom(7,300,0.5)+c(rep(34,5),rep(35,2))
set.seed(434445)
STR15BG <- rnbinom(7,300,0.5)+c(rep(23,4),rep(22,3))

STRBG <- c(STR1BG,STR2BG,STR3BG,STR4BG,STR5BG,
  STR6BG,STR7BG,STR8BG,STR9BG,STR10BG,
  STR11BG,STR12BG,STR13BG,STR14BG,STR15BG)

NHISdata$BG <- unlist(sapply(1:length(STRBG),function(i){rep(i,STRBG[i])}))

# Process of Checking -------------------------------------------------------------------------
sum(aggregate(NHISdata$BG~NHISdata$STRATUM,FUN=function(x){length(unique(x))})[,2])
table(NHISdata$BG)
table(NHISdata$STRATUM)
aggregate(NHISdata$FPX~NHISdata$STRATUM,FUN=function(x){length(x)})[,2]
aggregate(NHISdata$ID~NHISdata$STRATUM,FUN=function(x){length(x)})[,2]
sum(aggregate(NHISdata$ID~NHISdata$HHX+NHISdata$FMX,FUN=function(x){length(x)})[,3])
# Persons in BG and Combination of BG/STRATUM should be the same as BG is nested in STRATUM.
Cper1 <- aggregate(NHISdata$ID~NHISdata$BG,FUN=function(x){length(x)})[,2]
Cper2 <- aggregate(NHISdata$ID~NHISdata$STRATUM+NHISdata$BG,FUN=function(x){length(x)})[,3]
range(Cper1); range(Cper2)

write.csv(NHISdata,"C://Users//smosafer//Projects//FinalNHIS.csv")

# ============================================================================================
# Main Project -------------------------------------------------------------------------
# ============================================================================================

FinalNHIS <- read.csv("FinalNHIS.csv")

# Part A) SAMPLE SELECTION & DESIGN
# sample design: 2-stage stratified desgin, where:
#                                        starta v. = STRATUM (Total:15)
#                                       1-st stage = BG (eah BG has the average of 300 persons)
#                                       2-nd stage = FPX (Person)


MM <- 105 # total number of BGs
M <- 7 #BGs in each STRATUM
m <- 2  #number of selected BGs per STRATUM

Ni <- sapply(1:MM,function(i){length(FinalNHIS$FPX[FinalNHIS$BG==i])}) #BG size 
Ni <- sapply(0:14,function(i){list(Ni[c(((7*i)+1):((7*i)+7))])}) #BG size per stratum

probi <- sapply(1:15,function(i){m*Ni[[i]]/sum(Ni[[i]])}) #1st stage probability (pps): Ni is MOS
probi <- sapply(1:15,function(i){list(probi[,i])})

# selecting BGs from each STRATUM independently:
SAM <- sapply(1:15,function(i){cluster(data=FinalNHIS[FinalNHIS$STRATUM==i,],
                                       clustername="BG",size=m,method="systematic",
                                       pik=probi[[i]],description=TRUE)})

#Extract data from the sampled BGs:
SAMCLUS <- sapply(1:15,function(i){getdata(FinalNHIS[FinalNHIS$STRATUM==i,],as.data.frame(SAM[,i]))})  
SAMCLUS <- sapply(1:15,function(i){rename(as.data.frame(SAMCLUS[,i]),c(Prob="pi1"))})  

#Extract which BGs are selected per STRATUM:
selectedBG <- sapply(1:15,function(i){
  names(which(table(as.data.frame(SAMCLUS[,i])$BG)!=0))})
selectedBG <- sapply(1:15,function(i){list(selectedBG[,i])})

#Process of selecting Persons
ni <- sapply(1:15,function(i){ceiling((10/100)*c(length(FinalNHIS$ID[FinalNHIS$BG==selectedBG[[i]][1]]),
           length(FinalNHIS$ID[FinalNHIS$BG==selectedBG[[i]][2]])))})

Nisam <- sapply(1:15,function(i){c(length(FinalNHIS$ID[FinalNHIS$BG==selectedBG[[i]][1]]),
           length(FinalNHIS$ID[FinalNHIS$BG==selectedBG[[i]][2]]))})

dd <- as.vector(Nisam/sum(Nisam)) #vector of one draw PSU (BG) selection

#treat sample clusters as strata and select SRS sample from each
S <- sapply(1:15,function(i){strata(data=as.data.frame(SAMCLUS[,i]),stratanames="BG",
            size=ni[,i],method="srswor")})

#Extract the Observed data
SAMDATA <- sapply(1:15,function(i){getdata(as.data.frame(SAMCLUS[,i]),as.data.frame(S[,i]))})
SAMDATA <- sapply(1:15,function(i){rename(as.data.frame(SAMDATA[,i]),c(Prob="pi2"))})
# Merging & Final sampled Dataset:
SAMPLEDATA <- as.data.frame(rbind(as.data.frame(SAMDATA[,1]),as.data.frame(SAMDATA[,2]),as.data.frame(SAMDATA[,3]),
                               as.data.frame(SAMDATA[,4]),as.data.frame(SAMDATA[,5]),as.data.frame(SAMDATA[,6]),
                               as.data.frame(SAMDATA[,7]),as.data.frame(SAMDATA[,8]),as.data.frame(SAMDATA[,9]),
                               as.data.frame(SAMDATA[,10]),as.data.frame(SAMDATA[,11]),as.data.frame(SAMDATA[,12]),
                               as.data.frame(SAMDATA[,13]),as.data.frame(SAMDATA[,14]),as.data.frame(SAMDATA[,15])))

#Extract which SSUs (or persons) are selected
selectedFPX <- names(which(table(SAMPLEDATA$FPX)!=0))
SAMPLEDATA$wt <- unlist(1/SAMPLEDATA$pi1/SAMPLEDATA$pi2) #final weight

# Part B) SURVEY DESIGN DEFINE & POINT ESTIMATIONs
# Finding NA's in data set
summary(FinalNHIS)
NotNA <- FinalNHIS[, sapply(FinalNHIS, Negate(anyNA)), drop = FALSE]  #which variables do not contain NA's
sort(names(NotNA))

# Defining Some New Variable for the Estimations
# For the Age-group Domain
min(SAMPLEDATA$AGE_P); max(SAMPLEDATA$AGE_P)
SAMPLEDATA$AGEgrp <- cut(SAMPLEDATA$AGE_P, seq(15,85,10))

# SPECEQ: Have health problem that requires special equipment?
table(SAMPLEDATA$SPECEQ)
table(is.na(SAMPLEDATA$SPECEQ))
SAMPLEDATA$SPECEQ01 <- ifelse(SAMPLEDATA$SPECEQ==1,1,0) #1 stands for Yes

#0-1 Indicator of SEX
SAMPLEDATA$SEX01 <- SAMPLEDATA$SEX-1 #0:Male & 1:Female

#0-1 Indicator for CHDEV
#Have you EVER been told by a doctor to have Coronary heart disease?
SAMPLEDATA$CHDEV01 <- ifelse(SAMPLEDATA$CHLEV==1,1,0) #1 stands for Yes

SAMPLEDATA$CANEV1 <- ifelse(SAMPLEDATA$CANEV==1,1,0)  #indicator for CANEV

# Process of Defining the Designs:
RepSTRATUM <- aggregate(SAMPLEDATA$FPX~SAMPLEDATA$STRATUM,FUN=function(x){length(x)})[,2]
RepBG <- aggregate(SAMPLEDATA$FPX~SAMPLEDATA$STRATUM+SAMPLEDATA$BG,FUN=function(x){length(x)})[,3]

SAMPLEDATA$fpc1 <- c(unlist(sapply(1:14,function(i){rep(2245,RepSTRATUM[i])})),rep(2242,RepSTRATUM[15]))
SAMPLEDATA$fpc2 <- unlist(sapply(1:30,function(i){rep(as.vector(Nisam)[i],RepBG[i])}))

onestage.wfpc <- svydesign(id=~BG, strata=~STRATUM, fpc=~fpc1, weight=~wt, data=SAMPLEDATA, nest=TRUE)
onestage.wofpc <- svydesign(id=~BG, strata=~STRATUM, weight=~wt, data=SAMPLEDATA, nest=TRUE)
twostage.wfpc <- svydesign(id=~BG+FPX,strata=~STRATUM,fpc=~fpc1+fpc2,weight=~wt, data=SAMPLEDATA,nest=TRUE)
twostage.wofpc <- svydesign(id=~BG+FPX,strata=~STRATUM,weight=~wt, data=SAMPLEDATA,nest=TRUE)

# Measure of homogeneity (as an e.g.: for AGE):
BW2stagePPSe(Ni=as.vector(Nisam),ni=as.vector(ni),X=SAMPLEDATA$AGE_P,
             psuID=SAMPLEDATA$BG,w=SAMPLEDATA$wt,m=30,pp=dd)

# Process of Estimation and Design Effect:
# Interested Pointt Estimation for Age:
svymean(~AGE_P,onestage.wfpc,na.rm=TRUE,deff="replace")
svymean(~AGE_P,onestage.wofpc,na.rm=TRUE,deff="replace")
svymean(~AGE_P,twostage.wfpc,na.rm=TRUE,deff="replace")
svymean(~AGE_P,twostage.wofpc,na.rm=TRUE,deff="replace")

svyquantile(~AGE_P,onestage.wfpc,c(0.25,0.5,0.75),ci=TRUE) #Quantiles for AGE

svyby(~AGE_P,by=~SEX,design=onestage.wfpc,na.rm=TRUE,FUN=svymean, deff="replace")
svyby(~AGE_P,by=~SEX,design=onestage.wofpc,na.rm=TRUE,FUN=svymean,deff="replace")
svyby(~AGE_P,by=~SEX,design=twostage.wfpc,na.rm=TRUE,FUN=svymean,deff="replace")
svyby(~AGE_P,by=~SEX,design=twostage.wofpc,na.rm=TRUE,FUN=svymean,deff="replace")

svyby(formula=~factor(SEX),by=~AGEgrp,FUN=svymean,design=onestage.wfpc,na.rm=TRUE,deff="replace")
svyby(formula=~factor(SEX),by=~AGEgrp,FUN=svymean,design=onestage.wofpc,na.rm=TRUE,deff="replace")
svyby(formula=~factor(SEX),by=~AGEgrp,FUN=svymean,design=twostage.wfpc,na.rm=TRUE,deff="replace")
svyby(formula=~factor(SEX),by=~AGEgrp,FUN=svymean,design=twostage.wofpc,na.rm=TRUE,deff="replace")

# Interested Point Estimation for the SPECEQ:
svymean(~SPECEQ01,design=onestage.wfpc,na.rm=TRUE,deff="replace")
svymean(~SPECEQ01,design=onestage.wofpc,na.rm=TRUE,deff="replace")
svymean(~SPECEQ01,design=twostage.wfpc,na.rm=TRUE,deff="replace")
svymean(~SPECEQ01,design=twostage.wofpc,na.rm=TRUE,deff="replace")

svyby(~SPECEQ01,by=~SEX,design=onestage.wfpc,na.rm=TRUE,FUN=svymean,deff="replace")
svyby(~SPECEQ01,by=~SEX,design=onestage.wofpc,na.rm=TRUE,FUN=svymean,deff="replace")
svyby(~SPECEQ01,by=~SEX,design=twostage.wfpc,na.rm=TRUE,FUN=svymean,deff="replace")
svyby(~SPECEQ01,by=~SEX,design=twostage.wofpc,na.rm=TRUE,FUN=svymean,deff="replace")

svyby(formula=~factor(SEX),by=~SPECEQ01,FUN=svymean,design=onestage.wfpc,na.rm=TRUE,deff="replace")
svyby(formula=~factor(SEX),by=~SPECEQ01,FUN=svymean,design=onestage.wofpc,na.rm=TRUE,deff="replace")
svyby(formula=~factor(SEX),by=~SPECEQ01,FUN=svymean,design=twostage.wfpc,na.rm=TRUE,deff="replace")
svyby(formula=~factor(SEX),by=~SPECEQ01,FUN=svymean,design=twostage.wofpc,na.rm=TRUE,deff="replace")

svychisq(~SPECEQ01+SEX,design=onestage.wfpc,statistics="F") #checking the dependency of SEX and SPECEQ


svyby(~SPECEQ01,by=~AGEgrp,design=onestage.wfpc,na.rm=TRUE,FUN=svymean,deff="replace")
svyby(~SPECEQ01,by=~AGEgrp,design=onestage.wofpc,na.rm=TRUE,FUN=svymean,deff="replace")
svyby(~SPECEQ01,by=~AGEgrp,design=twostage.wfpc,na.rm=TRUE,FUN=svymean,deff="replace")
svyby(~SPECEQ01,by=~AGEgrp,design=twostage.wofpc,na.rm=TRUE,FUN=svymean,deff="replace")

svyby(formula=~factor(AGEgrp),by=~SPECEQ01,FUN=svymean,design=onestage.wfpc,na.rm=TRUE,deff="replace")
svyby(formula=~factor(AGEgrp),by=~SPECEQ01,FUN=svymean,design=onestage.wofpc,na.rm=TRUE,deff="replace")
svyby(formula=~factor(AGEgrp),by=~SPECEQ01,FUN=svymean,design=twostage.wfpc,na.rm=TRUE,deff="replace")
svyby(formula=~factor(AGEgrp),by=~SPECEQ01,FUN=svymean,design=twostage.wofpc,na.rm=TRUE,deff="replace")

svychisq(~SPECEQ01+AGEgrp,design=onestage.wfpc,statistics="F") #checking the dependency of SEX and SPECEQ

# Interested Point Estimation for SEX
svymean(~SEX01,onestage.wfpc,na.rm=TRUE,deff="replace")
svymean(~SEX01,onestage.wofpc,na.rm=TRUE,deff="replace")
svymean(~SEX01,twostage.wfpc,na.rm=TRUE,deff="replace")
svymean(~SEX01,twostage.wofpc,na.rm=TRUE,deff="replace")

# Interested Point Estimation for CHDEV
svymean(~CHDEV01,onestage.wfpc,na.rm=TRUE,deff="replace")
svymean(~CHDEV01,onestage.wofpc,na.rm=TRUE,deff="replace")
svymean(~CHDEV01,twostage.wfpc,na.rm=TRUE,deff="replace")
svymean(~CHDEV01,twostage.wofpc,na.rm=TRUE,deff="replace")


# Part C) WEIGHTING & Nonresponse Adjustment

#Base Weight
sum(1/(SAMPLEDATA$pi1*SAMPLEDATA$pi2))
SAMPLEDATA$wt
sum(SAMPLEDATA$wt)

# QC checking
length(unique(SAMPLEDATA$STRATUM))
length(unique(SAMPLEDATA$BG))
aggregate(SAMPLEDATA$BG~SAMPLEDATA$STRATUM,FUN=function(x){length(unique(x))})
aggregate(SAMPLEDATA$BG~SAMPLEDATA$STRATUM,FUN=function(x){unique(x)})
sum(unlist(probi))

ni <- as.vector(ni)
FPX1 <- aggregate(SAMPLEDATA$FPX~SAMPLEDATA$STRATUM+SAMPLEDATA$BG,FUN=function(x){length(x)})[,3]
FPX2 <- aggregate(SAMPLEDATA$FPX~SAMPLEDATA$BG,FUN=function(x){length(x)})[,2]
all.equal(ni,FPX1)
all.equal(ni,FPX2)

aggregate(SAMPLEDATA$pi1~SAMPLEDATA$STRATUM,FUN=function(x){length(unique(x))})
PI1 <- aggregate(SAMPLEDATA$pi1~SAMPLEDATA$STRATUM,FUN=function(x){unique(x)})

aggregate(SAMPLEDATA$pi2~SAMPLEDATA$BG,FUN=function(x){length(unique(x))})
PI2 <- aggregate(SAMPLEDATA$pi2~SAMPLEDATA$BG,FUN=function(x){unique(x)})[,2]
Nisam <- as.vector(Nisam)  
sampsize <- sapply(1:length(Nisam),function(i){sum(rep(PI2[i],Nisam[i]))})
all.equal(sampsize,ni)

# Mising Mechanism
pix <- pGumbel(1,0,1)
SAMPLEDATA$MISS <- rbinom(nrow(SAMPLEDATA),1,
      .5*((SAMPLEDATA$SEX==1 & SAMPLEDATA$AGE_P>=70)+
            (SAMPLEDATA$RACERPI2!=1|2*((SAMPLEDATA$SEX==1)*pix
                                     +(SAMPLEDATA$AGE_P>=70)*(1-pix))/(1+2*pix*(pix -1)))))==1

# Missing Indicator
SAMPLEDATA$M <- ifelse(SAMPLEDATA$MISS==TRUE,1,0)

# Percentage of Nonresponse 
mean(SAMPLEDATA$M)*100
## Response Indicator
SAMPLEDATA$R <- ifelse(SAMPLEDATA$M==0,1,0)

# unweighted logistic model
glm.logit <- glm(R~AGE_P+as.factor(SEX)+as.factor(RACERPI2)+
                   as.factor(CANEV1)+as.factor(APSBSCHK),
                 family=binomial(link=logit),data=SAMPLEDATA)
summary(glm.logit,cor=F)
anova(glm.logit,test = "Chisq")

regTermTest(glm.logit,~as.factor(CANEV1)+as.factor(APSBSCHK), method="Wald",df=NULL)

glm.probit <- glm(R~AGE_P+as.factor(SEX)+as.factor(RACERPI2)+
                    as.factor(CANEV1)+as.factor(APSBSCHK),
                  family=binomial(link=probit),data=SAMPLEDATA)
summary(glm.probit,cor=F)
anova(glm.probit,test = "Chisq")
regTermTest(glm.logit,~as.factor(CANEV1)+as.factor(APSBSCHK), method="Wald",df=NULL)
  
glm.cloglog <- glm(R~AGE_P+as.factor(SEX)+as.factor(RACERPI2)+
                     as.factor(CANEV1)+as.factor(APSBSCHK),
                   family=binomial(link=cloglog),data=SAMPLEDATA) 
summary(glm.cloglog,cor=F)
anova(glm.cloglog,test = "Chisq")
regTermTest(glm.cloglog,~as.factor(CANEV1)+as.factor(APSBSCHK), method="Wald",df=NULL)

par(mfrow=c(2,3))

dat <- data.frame(Logit=glm.logit$fitted.values,Probit=glm.probit$fitted.values,
                  Cloglog=glm.cloglog$fitted.values)
long <- melt(dat)
plot(value ~ variable, data=long)

plot(glm.logit$fitted.values,glm.probit$fitted.values,
     xlab="logit predictions",ylab="probit predictions")
abline(0,1)
plot(glm.logit$fitted.values,glm.cloglog$fitted.values,
     xlab="logit predictions",ylab="cloglog predictions")
abline(0,1)

# weighted logistic model

glm.complex.probit <- svyglm(R~AGE_P+as.factor(SEX)+as.factor(RACERPI2)+
                               as.factor(CANEV1)+as.factor(APSBSCHK),
                             family=binomial(link=probit),data=SAMPLEDATA,
                             design=onestage.wfpc)
summary(glm.complex.probit,cor=F)
anova(glm.complex.probit,test = "Chisq")
regTermTest(glm.complex.probit,~as.factor(CANEV1)+as.factor(APSBSCHK), method="Wald",df=NULL)

glm.complex.logit <- svyglm(R~AGE_P+as.factor(SEX)+as.factor(RACERPI2)+
                              as.factor(CANEV1)+as.factor(APSBSCHK),
                            family=binomial(link=logit),data=SAMPLEDATA,
                            design=onestage.wfpc)
summary(glm.complex.logit,cor=F)
anova(glm.complex.logit,test = "Chisq")
regTermTest(glm.complex.logit,~as.factor(CANEV1)+as.factor(APSBSCHK), method="Wald",df=NULL)

glm.complex.cloglog <- svyglm(R~AGE_P+as.factor(SEX)+as.factor(RACERPI2)+
                                as.factor(CANEV1)+as.factor(APSBSCHK),
                              family=binomial(link=cloglog),data=SAMPLEDATA,
                              design=onestage.wfpc)
summary(glm.complex.cloglog,cor=F)
anova(glm.complex.cloglog,test = "Chisq")
regTermTest(glm.complex.cloglog,~as.factor(CANEV1)+as.factor(APSBSCHK), method="Wald",df=NULL)


dat.wgt <- data.frame(Logit.wgt=glm.complex.logit$fitted.values,
                  Probit.wgt=glm.complex.probit$fitted.values,
                  Cloglog.wgt=glm.complex.cloglog$fitted.values)
long.wgt <- melt(dat.wgt)
plot(value ~ variable, data=long.wgt)

plot(glm.complex.logit$fitted.values,glm.complex.probit$fitted.values,
     xlab="logit weighted predictions",ylab="probit weighted predictions")
abline(0,1)
plot(glm.complex.logit$fitted.values,glm.complex.cloglog$fitted.values,
     xlab="logit weighted predictions",ylab="cloglog weighted predictions")
abline(0,1)

par(mfrow=c(1,3))
plot(glm.probit$fitted.values,glm.complex.probit$fitted.values,
     xlab="Unweighted Predicted Probabilities",
     ylab="Weighted Predicted Probabilities",main="Probit Link Function")
abline(0,1)
plot(glm.logit$fitted.values,glm.complex.logit$fitted.values,
     xlab="Unweighted Predicted Probabilities",
     ylab="Weighted Predicted Probabilities",main="Logit Link Function")
abline(0,1)
plot(glm.cloglog$fitted.values,glm.cloglog$fitted.values,
     xlab="Unweighted Predicted Probabilities",
     ylab="Weighted Predicted Probabilities",main="CLogLog Link Function")
abline(0,1)

# Determine quantiles of response propensities
options(digits=5)
quantiles <- quantile(glm.probit$fitted.values,probs=seq(0,1,0.2))
p.class <- cut(round(glm.probit$fitted.values,20),breaks=quantiles,
               include.lowest=TRUE)
SAMPLEDATA$p.class <- p.class
table(p.class)

#comparing ways of estimating the class response propensity
UA <- by(data=glm.probit$fitted.values,p.class,mean) #unweighted average
UAvec <- as.vector(UA)

fitted <- glm.probit$fitted.values
weight <- SAMPLEDATA$wt
WRP <- by(data = data.frame(fitted,wt = weight),p.class, #weighted response propensity
   function(x) {weighted.mean(x$fitted, x$wt)})
WRPvec <- as.vector(WRP)

URR <- by(as.numeric(SAMPLEDATA[,"R"]),p.class,mean) #unweighted response rate
URRvec <- as.vector(URR)

WRR <- by(data = data.frame(resp = as.numeric(SAMPLEDATA$R), #weighted response rate
                     wt = SAMPLEDATA$wt), p.class,
   function(x) {weighted.mean(x$resp, x$wt)})
WRRvec <- as.vector(WRR)

UMRP <- by(fitted, p.class, median)  #unweighted median response propensity
UMRPvec <- as.vector(UMRP)

SAMPLEDATA$p.class <- p.class
Frameprob <- as.data.frame(cbind(glm.probit$fitted.values,label=p.class))
levels(Frameprob$label) <- c("[0.0055,0.369]", "(0.369,0.599]", "(0.599,0.746]", "(0.746,0.895]", "(0.895,1]")
dev.off()
boxplot(Frameprob$V1~Frameprob$label, main="Boxplot of Predicted Probabilities")


#Illustrative checks
chk1age <- glm(AGE_P~p.class+R+p.class*R,data=SAMPLEDATA)
summary(chk1age)

chk2age <- glm(AGE_P~p.class,data=SAMPLEDATA)
anova(chk1age,chk2age,test="F")


chk1sex <- glm(SEX01~p.class+R+p.class*R,data=SAMPLEDATA)
summary(chk1sex)

chk2sex <- glm(SEX01~p.class,data=SAMPLEDATA)
anova(chk1sex,chk2sex,test="F")

# Poststratification based on Age group

FinalNHIS$AGEgrp <- cut(FinalNHIS$AGE_P,seq(15,85,10))
FREQ <- as.vector(table(FinalNHIS$AGEgrp))
pop.types <- data.frame(AGEgrp=c("(15,25]","(25,35]",
                                   "(35,45]", "(45,55]", "(55,65]",
                                   "(65,75]", "(75,85]"), Freq=FREQ)

onestage.wfpc.p <- postStratify(onestage.wfpc,~AGEgrp,pop.types)

rbind(summary(weights(onestage.wfpc)), summary(weights(onestage.wfpc.p)))

svymean(~WTIA_SA, onestage.wfpc)
svymean(~WTIA_SA, onestage.wfpc.p)
svytotal(~WTIA_SA, onestage.wfpc)
svytotal(~WTIA_SA, onestage.wfpc.p)
c(cv(svymean(~WTIA_SA, onestage.wfpc)), cv(svymean(~WTIA_SA, onestage.wfpc.p)))
c(cv(svytotal(~WTIA_SA, onestage.wfpc)), cv(svytotal(~WTIA_SA, onestage.wfpc.p)))

#Imputation practice
NAvariable <- colnames(FinalNHIS)[colSums(is.na(FinalNHIS)) > 0]  #which variables contain NA's
sort(NAvariable)

SUBSAMPLE <- SAMPLEDATA[1:30,]
SUBSAMPLE$AGE <- replace(SUBSAMPLE$AGE_P,SUBSAMPLE$AGE_P>=70,NA)
SUBSAMPLE <- SUBSAMPLE[,c("CHLYR", "HLTHPROM", "ASISIM", "NIGHTWK","AGE")]

marginplot(SUBSAMPLE[, c("HLTHPROM","AGE")], col = mdc(1:2), cex = 1.2,
           cex.lab = 1.2, cex.numbers = 1.3, pch = 19)  

SUBSAMPLE.imp <- mice(SUBSAMPLE, seed = 23109)
summary(SUBSAMPLE.imp)

fit <- with(SUBSAMPLE.imp, lm(NIGHTWK ~ AGE + HLTHPROM))
summary(fit)





