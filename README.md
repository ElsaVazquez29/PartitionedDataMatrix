# PartitionedDataMatrix
This repository contains SAS macro %partitionedDataMatrix to get lagged covariates with valid moments for partitioned Bayesian models

# example of use to get partitioned data matrix and then fit a Partitioned method of valid moments model with Bayesian intervals in SAS

Syntax of this macro is the same as that of %partitionedGMM macro, but by running it you only get partitioned data matrix and don't fit any model.
https://www.sas.com/content/dam/SAS/support/en/sas-global-forum-proceedings/2018/2661-2018.pdf



%partitionedDataMatrix(file=add, idvar=ID,outVar=bmi, timeVar=wave, PredVarTD=alcohol feelingscale tvhrs activityscale,
predVarTI=race_, distr=bin, MC=LWY);

The dataset is add and specified in file argument. 
Subjects' ID in dataset "add" is ID and is specified in argument idvar=.
The outVar argument asks the name of the outcome variable in the dataset.
The timeVar= argument asks for variable in the dataset that contains the time variable; In the Add health dataset that variable is wave. 
The PredVarTD argument asks for all time-dependent covariates in the model that we will test for valid moment conditions.
The PredVarTI asks for covariates that are tie-independent (don't change overt iem), in this case just race.
The distr= argument asks if outcome variable is binary (bin) or continous (normal).


The %partitionedDataMatrix outputs a dataset called "mydata3" that contains the outcome and the partitioned data matrix which results for
testing for valid moment conditions

#you can rename "mydata3" outputed dataset
data obesity;
set mydata3;
run;


#use the dataset with partitioned data matrix to fit Partitioned method of valid moments model with Bayes intervals in PROC MCMC

proc mcmc data=obesity NBI=10000 nmc=5000000 THIN=5 seed=7465 stats=all diagnostics=all;
parms b0-b12;
prior b: ~ normal (0, var=10000);
mu=b0 + b1*race_ + b2*alcohol_0 + b3*FEELINGSCALE_0 + b4*TVHRS_0 + b5*ACTIVITYSCALE_0 + b6*ALCOHOL_1
   + b7*FEELINGSCALE_1 + b8*TVHRS_1 + b9*ACTIVITYSCALE_1 + b10*ALCOHOL_2 + b11*ACTIVITYSCALE_2 
   + b12*ACTIVITYSCALE_3;
MU1=EXP(MU)/(1+EXP(MU));
model bmi~binomial(1, p=MU1);
run;



