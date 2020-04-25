R code for the NHIS project ("Data Collection and Analytics Specialization")

Open the R code which contains some commands with the following explanations.


# SOME DESCRIPTIONS


Explanation of Produced Data set:
1. We first order the dataset based on the hierarchical elements of "STRAT_P, PSU_P, HHX, FMX, FPX" in
order not to encounter any problems for the further steps. Then we assign each person a unique id called
ID in the data set.
2. We produce 15 strata with roughly same sizes. This is under a new variable STRATUM.
3. Then in order to produce block groups which are nested in the stratum, we use rnbinom(7,300,0.5).

We exactly have 7 block groups in each stratum.

We need to do manual changes to cover the difference of population size of stratum and
summation of a vector generated by "rnbinom". We assign simple seed that the results of block
group size creation could be replicated.

The name of block group variable is BG. We totally have 105 block groups.
The range of population per block group is from 267 to 358. 


## STEPS

#### Part I) SAMPLE SELECTION & DESIGN

We consider stratified 2-stage sample design, where strata called "STRATUM", 1st stage (PSU) is block-
group called "BG", and 2nd stage (SSU) is person called "FPX" which we have given them unique "ID".

#### Part II) SURVEY DESIGN DEFINITION & POINT ESTIMATIONs

Here we first need to assign a suitable svydesign function to the selected data set in R. In multistage
design, we have a situation where we can only consider the first stage of sample selection and ignore the
rest of stages for the point estimations and variance estimation as well. We call this: ultimate cluster
method. For this project, we assign two survey designs while keeping the strata. These function designs
are as follows:

1. Onestage: when we ignore the second stage of sampling (Ultimate cluster)

With fpc
```
svydesign(id=~BG, strata=~STRATUM, fpc=~fpc1, weight=~wt, data=SAMPLEDATA,
  nest=TRUE)
```
Without fpc
```
svydesign(id=~BG, strata=~STRATUM, weight=~wt, data=SAMPLEDATA, nest=TRUE)
```
2. Twostage: when we consider both stages of sampling

With fpc
```
svydesign(id=~BG+FPX,strata=~STRATUM,fpc=~fpc1+fpc2,weight=~wt,
  data=SAMPLEDATA,nest=TRUE)
```
Without fpc
```
svydesign(id=~BG+FPX,strata=~STRATUM,fpc=~fpc1+fpc2,weight=~wt,
  data=SAMPLEDATA,nest=TRUE)
```  
Then, do the Design Effect Comparison and Point Estimations.  

#### Part III) WEIGHTING & NONRESPONSE ADJUSTMENT

Base Weight:

Base weight is the inverse of selection probability of individual "i". We should keep in mind that the
summation of base weights is equal to the total number of elements in the population.
In our R code the base weight is given by "wt" (See the R code).

Quality Control (QC) Checks of Weighting:

For this part, we need to check the following steps.

1. In the final sample that they extracted from the population, PSUs (BGs) should be nested in the Stratum.
The length of unique BGs per STRATUM should be 2.
2. The length of unique probability selection of BG from the first stage per STRATUM should be 2.
3. The length of unique probability selection of person in the second stage per BG should be 1.
4. The summation of weights for all of selected individuals in the sample should be equal to the size of
population.
5. The summation of first stage selection probabilities for all of BGs in the population should be equal to
the total number of selected BGs, which is 30 based on our work.
6. The summation of second stage selection probabilities for all of individuals within each selected BG
should be equal to the sample size of individuals from each selected BG.

Propensity Score Adjustment:

See the following command from R code
```
pGumbel(1, mu = 0, sigma = 1)= 0.6922006
```
Illustration of Some Quality Checking (QC):

Check1
We illustrate a check on covariate balance by fitting an ANOVA model to AGE, which is continuous. We do
not use the survey weights for this analysis since the interest is in whether balance has been achieved in
the sample that was selected. As an example for our work, we can have the following R code and results:
```
chk1 <- glm(AGE_P~p.class+R+p.class*R,data=SAMPLEDATA)
summary(chk1)
```
Check2
Another check is to fit a second model that includes only "p.class" and to test whether the models are
equivalent? From the following Box, we can realize the models are not equivalent. Therefore, it is
important to keep the interaction term as it was also expected.
```
chk2 <- glm(AGE_P~p.class,data=SAMPLEDATA)
anova(chk1,chk2,test="F")
```

Calibration to Population Control Total:

The idea behind of calibration is using auxiliary information to reduce the variance. Two common ways of
calibration are: Poststratification and Raking. Here based on the Course materials, we recommend the
Poststratification. We consider Age group as the poststrata. The results based on our work are as follows:
```
FinalNHIS$AGEgrp <- cut(FinalNHIS$AGE_P,seq(15,85,10))

FREQ <- as.vector(table(FinalNHIS$AGEgrp))

pop.types <- data.frame(AGEgrp=c("(15,25]","(25,35]",
                                 "(35,45]", "(45,55]", "(55,65]",
                                 "(65,75]", "(75,85]"), Freq=FREQ)
                                   
onestage.wfpc.p <- postStratify(onestage.wfpc,~AGEgrp,pop.types)

rbind(summary(weights(onestage.wfpc)), summary(weights(onestage.wfpc.p)))
```
Then we need to check that weights are calibrated for x’s.

If we consider the weight response variable WTIA_SA , we can have the following results:
```
svymean(~WTIA_SA, onestage.wfpc)

svymean(~WTIA_SA, onestage.wfpc.p)

svytotal(~WTIA_SA, onestage.wfpc)

svytotal(~WTIA_SA, onestage.wfpc.p)

c(cv(svymean(~WTIA_SA, onestage.wfpc)), cv(svymean(~WTIA_SA, onestage.wfpc.p)))

c(cv(svytotal(~WTIA_SA, onestage.wfpc)), cv(svytotal(~WTIA_SA, onestage.wfpc.p)))
```

Impuation Practice:
Most of the time we encounter item nonresponse that a person has not answered to the one or more
than one question in the survey, which is called item non-response. To conduct this part of project,
students need to use "mi" and "mice" packages from R.
Students can consider a small subset of sample that they have drawn from Part I of project and consider
variables with NA values. For example, we have considered a sub-sample of size 30 from the sample drawn
from Part I. The variables that we have considered are:

*CHLYR, HLTHPROM, ASISIM, NIGHTWK, AGE.*

The variable AGE has been created by ourselves in order to have a continuous variable with NA’s values
in the sub-sample. The original AGE variable does not have NA values. The condition that we considered
for the AGE with NA values is: if AGE≥70, then we have non-response for the AGE.


#### Part IV) EXTRA DATA ANALYSIS

After finishing the former parts, we consider the sample extracted from Part I
and a design from Part II. Then we consider an interested
response variable such as SPECEQ as a response variable "y" with following definition:
"SPECEQ: Do you have any health problem that requires you to use special equipment?"
We can recode this variable to 1: Yes and 0: No or other answers.
Then we can find the related covariates such as age, race, sex, etc. using svyglm. Further, we can do
residual analysis, model checking and the odds ratio of the characteristics of having a special health
problem v.s. not having a special health problem requiring special equipment: p/(1 − p).
We also recommend 2 different situations:

Conduct data analysis assuming all of units responded

Conduct data analysis assuming that we have non-respondents following the Part III.

Then, we can investigate the effect of missingness on the estimators when we do not consider the non-
response adjustment and when we consider non-response adjustment.
