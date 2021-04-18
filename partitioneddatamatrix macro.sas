/*************************************************************************************************************/
/**	SAS Code: %partitionedDataMatrix				                                    **/
/** Programmer: Elsa Vazquez Arreola					                                    **/
/** Description: Creates partitioned Data matrix for simultaneous modeling                                  **/
/** This code can also be used to run Partitioned MVM model with Bayes intervals                            **/
/** This code is a modification of code for %partitionedGMM macro by Kyle Irimata                           **/
/** The partitioned data matrix is created in %partitioneGMM macro by Kyle Irimata                          **/
/** This macro tests valid moment conditions and creates partitioned data matrix but does not fit any model **/
/*************************************************************************************************************/
%macro partitionedDataMatrix(ds=., file=, timeVar=, outVar=, predVarTD=, idVar=, alpha=0.05, predVarTI=., distr=bin, optim=NLPCG, MC=LWY);

*Check if there is a library;
%if &ds. ^=. %then %do;
	LIBNAME DS &ds.;  
	DATA mydata; 
		SET DS.&file.; 
	RUN;
%end;
%else %do;
	DATA mydata;
		SET &file.;
	RUN;
%end;

TITLE "Partitioned GMM";
PROC SORT DATA=mydata OUT=mydatasorted; 
BY &timeVar.; RUN;

*Check if there is a time independent variable;
%if &predVarTI. ^=.  %then %do;
	%let predVar = &predVarTD &predVarTI;
%end;
%else %do;
	%let predVar = &predVarTD;
%end;

*Use either the LWY or LS moment check;
%if &MC=LWY %then %do; *LWY approach;

*Obtain residuals;
%if &distr.=normal %then %do;
	proc reg data=mydatasorted NOPRINT;
	BY &timeVar.;
	MODEL &outVar. = &predVar.;
	OUTPUT OUT=outpool3 PREDICTED=mu;
	RUN;

	DATA outpool3; SET outpool3;
	wt = 1/mu; 
	rsdraw = &outVar.-mu; 
%end;
%else %if &distr.=bin %then %do;
	PROC logistic DATA=mydatasorted NOPRINT; 
	BY &timeVar.;
	MODEL &outVar. (event='1') = &predVar. / aggregate scale=none;
	OUTPUT OUT=outpool3 P=mu;
	RUN;

	DATA outpool3; SET outpool3;
	wt = mu*(1-mu);
	rsdraw = &outVar.-mu; 
%end;
%else %do;
	%put ERROR must be normal or binomial distributions; 
	%return;
%end;

PROC SORT DATA=outpool3 OUT=outpool3 ;
  BY &idVar. &timeVar.; RUN;
quit;

PROC IML;
use outpool3;                                                           
%if &predVarTI. ^=. %then %do;
	read all VARIABLES {&predVarTD. &predVarTI.} into Zmat; 
	read all var {&predVarTI.} into timeInd;
	read all var {&predVarTD.} into timeDep;
%end;
%else %do;
	read all VARIABLES {&predVarTD.} into Zmat;
%end;
read all var {wt} into wt;
read all var {rsdraw} into rsd;
read all var {&idVar.} INTO ID;
read all var {&timeVar.} INTO time;
close outpool3;


use mydata;
read all var{&predVarTD.} into timeDep2;
read all var {&predVarTI.} into timeInd2;
close mydata;



N=ncol(unique(ID)); *Number of subjects;
T = max(time); *Number of time points;
Np = ncol(Zmat); *Number of predictors, not including intercept;

NpTD2 = ncol(timeDep2);
NpTI2=ncol(timeInd2);

%if &predVarTI. ^=. %then %do;
	NpTI = ncol(timeInd);*The last NPTI columns of Zmat will be time independent;
	NpTD = ncol(timeDep); *The first NPTD columns of Zmat will be time dependent;
%end;

*Find the correlation between X (a) and Y (rsd) across the T time points;
start rho(a,rsd) global(N,T);
abm = j(N,2*T,.);
abm[,1:T] = shape(rsd,N);	* N x T - First T columns are the values of rsd sorted into T columns (by time);
abm[,T+1:2*T] = shape(a,N); * Remaining T columns are the values of a sorted into T columns;
corr = corr(abm);  
rho = corr[1:T,T+1:2*T];    * T x T - Only take the correlations between the X and Y in the T time points;
return(rho);
finish rho;


*Standard deviation for each correlation;
start stddev(a,rsd) global(N,T);
bm = shape(rsd,N);    		 * N x T;
bdev = bm-j(N,1,1)*bm[:,];   * bdev N x T,   bm[:,] Col Mean is a row vector   1 x T - each row of bm minus column means;
bdev2 = bdev#bdev;      
am = shape(a,N);   
adev = am-j(N,1,1)*am[:,];  *N x T;
adev2 = adev#adev;      
stddev = sqrt( (1/N)*t(bdev2)*adev2 );   * T x T;
return(stddev);
finish stddev;

* corrected standardization;
start stdzn(x) global(N,T);
xrows = shape(x,N);   *N x T - by shape default columns are T=nrow(x)/N;
y = xrows - xrows[:,];  *N x T - Each value minus the column mean (mean for that time point);
vcv = (1/(N-1))*t(y)*y; * T x T;
v = diag(vcv); * T x T diagonal elements of vcv;
sinv = sqrt(inv(v));
x2 = y*sinv;   *N x T;
x2  = shape(x2,N*T,1); *N*T x 1 vector of standardized values;
return(x2);
finish stdzn;

pvec = j(Np*T*T,1,.);sevec = j(Np*T*T,1,.);  * pvec   (T*T) x Np;
r4out = j(T,T*Np,.); se4out = j(T,T*Np,.); z4out = j(T,T*Np,.); p4out = j(T,T*Np,.); 

y = rsd;
y_std = stdzn(y);




DO i=1 TO Np;
x = wt#Zmat[,i]; 			 * (N*T) x 1;
x_std = stdzn(x);

*Find p-values for the correlation of X_i and Y;
r = rho(x_std,y_std);		 * T x T;
se = stddev(x_std,y_std);	 * T x T;
z = sqrt(N)*(r/se);
p = 2*(1-cdf("normal",abs(z)));  * T x T;

*Fill the corresponding T columns of each matrix with the calculated values;
r4out[,T*(i-1)+1:T*i] = r;
se4out[,T*(i-1)+1:T*i] = se;
z4out[,T*(i-1)+1:T*i] = z;
p4out[,T*(i-1)+1:T*i] = p;

DO j = 1 TO T;
p[j,j] = 1;  * Not going to test diagonal elements, set to 1;
END;
*Takes the values in p4out and se4out, across row, then down column and creates vectors;
pvec[T*T*(i-1)+1:T*T*i,1] = shape(p,T*T,1);    *(T*T) x 1;
sevec[T*T*(i-1)+1:T*T*i,1] = shape(se,T*T,1);   *(T*T) x 1;
END;


TypeVec2 = (pvec >= &alpha.*j(Np*T*T,1,1) );

Type2 = shape(TypeVec2, Np, T*T);


*Individual test approach;
x = wt; 			 * (N*T) x 1;
x_std = stdzn(x);
r = rho(x_std,y_std);		 * T x T;
se = stddev(x_std,y_std);	 * T x T;
z = sqrt(N)*(r/se);
p = 2*(1-cdf("normal",abs(z)));  * T x T;


T_wt = p>&alpha.;
T_wt = shape(T_wt,1);


TypeMtx3 = j(Np+1,T*T,.);
TypeMtx3[1,] = T_wt;
TypeMtx3[2:Np+1,] = Type2;
Type2[1,] = shape(I(T),1);


*Adjust the last NpTI rows of TypeMtx3 to account for time independent covariates;
%if &predVarTI. ^=. %then %do;
	TypeMtx3= TypeMtx3[1:(NPTD+1),]; *Subset TypeMtx3 to include only time dependent;
	TypeMtxIND = repeat(shape(I(T),1),NPTI); *Create rows of typeMtx for the time independent;
	Np = NpTD;
%end;

print typemtx3; /*typemtx3 contains all valid moments for each of the time-dependent covariates*/
/*TypeMtx3 is not affected by the number of time-independent covariates in the model*/


/* define helper functions ROW and COL */
start row(x);  /* return matrix m such that m[i,j] = i */
   return( repeat( T(1:x), 1, x ));
finish;
start col(x);  /* return matrix m such that m[i,j] = j */
   return( repeat(1:x, x) );
finish;

*Helper matrices;
rt = row(T);
ct = col(T);

PRINT RT;
PRINT CT;


*Indices of upper diagonal for removal;
upperIdx = loc(ct>rt);
lowerIdx= loc(ct<rt);


*Remove backwards conditions;
typeMtxNB = TypeMtx3;
typeMtxNB[,upperIdx] = 0;



/**************************************************************************/
/**************************************************************************/
/**************************************************************************/
/*****************  MATRIX FOR MODIFIED DATASET ***************************/
/**************************************************************************/
/**************************************************************************/
/**************************************************************************/


MatrixCreateNew=j(T,NPTD2,1);
/*creating matrix from first time dependent covariate*/


/*creating matrix for all other time dependent covariate*/

do i=1 to T-1;
	MatrixC1=j(T,NPTD2,0);
	do j=1 to NpTD2;
      MatrixC2=shape(typeMtxNB[j+1,], T, T);
	  MatrixR=j(T-i, T-i,.);
      MatrixR[1:T-i, 1:T-i]=MatrixC2[i+1:T, 1:T-i];
      VecC3=vecdiag(MatrixR);
      MatrixC1[i+1:T, j]=vecC3;
	end;
MatrixCreateNew=MatrixCreateNew||MatrixC1;
end;

cov_id = MatrixCreateNew[+,];

keepcovar_id = loc(cov_id>0);


ncolumns=ncol(keepcovar_id);
matrixCreateNew2=matrixCreateNew[,keepcovar_id];



MatrixCreateNew3 = j(N*T,ncolumns,0);
DO i=1 TO N;
  MatrixCreateNew3[(i-1)*T+1:i*T, 1:ncolumns] = MatrixCreateNew2[1:T, 1:ncolumns];
END;



/*****************************************************************************************/
/*****************************************************************************************/
/*****************************************************************************************/
/*****************************************************************************************/


/*Keep backwards conditions*/
typeMtxfut=typeMtx3;
typeMtxFut[,loweridx]=0;

/*CHECKING HOW TO CHANGE ORDER OF COVARIATES FOR OUTCOMES AS PREDICTORS*/

CREATE TypeMtxFut from TypeMtxFut;
APPEND from TypeMtxFut;
CLOSE TypeMtxFut;



*Create the combination TypeMtx, of maximum size to be pared down later, excluding the intercept;
TypeMtxCombo = j((Np*T),T*T,0); 

TypeMtxCombo2 = j((Np*T),T*T,0);
*Loop across the different types, starting with current;
DO i = 0 TO (T-1);
	*Identify the indices for the appropriate setting;
	Idx = loc(rt-ct = i);
	IdxDif = setdif(1:T*T, Idx);/*THIS RETURNS ALL THOSE INDICES OF THE MATRIX FOR WHICH RT-CT IS DIFFERENT FROM I*/
	temp = typeMtxNB[2:(Np+1),];
	temp[,IdxDif] = 0;
	TypeMtxCombo[(i*Np+1):((i+1)*Np),] = temp;

	*Shift the values to accomodate shifted Y-X relationship;
	if(i>0) then DO; 
		dummyZero = j(Np, i, 0) || temp;
		TypeMtxCombo2[(i*Np+1):((i+1)*Np),] = dummyZero[,1:(T*T)];
	END;
	
	ELSE TypeMtxCombo2[(i*Np+1):((i+1)*Np),] = temp;

END; 



/***************************************************************************************************************/
/********** WILL HAVE TO CREATE SOMETHING SIMILAR TO TYPEMTXCOMBO 2 FOR EACH TIME-DEPENDENT COVARIATE***********/
/**************THIS USING TYPEMTXFUT****************************************************************************/
/***************************************************************************************************************/


*Identify predictors or relationships with no valid moment conditions;

neq = TypeMtxCombo2[,+]; /*ADDING UP NUMBER OF 1'S FOR EACH ROW OF TYPEMTXCOMBO2*/




*Print a note about omitted moment conditions;
if min(neq)=0 then do;
	zeroPred = loc(neq=0);
	Note = "There are no valid moment conditions for " +char(ncol(zeroPred))+" covariate relationship(s).";
	print Note[label="Moment Condition Notes"];
	print "These covariate relationships will be omitted in the analysis.";
end;
else print "All covariate relationships will be evaluated."[label="Moment Condition Notes"];


/*FOR THOSE ROWS RELATED TO LAGS OF THE COVARIATE WITH NO VALID MOMENTS WE DELETE THAT ROW*/
/*THIS BECAUSE WE WILL NOT BE USING THOSE LAGS IN OUR MODELS*/
keepPred = loc(neq>0);

TypeMtxCombo3 = TypeMtxCombo2[keepPred,]; /*WILL ONLY KEEP ROWS OF COVARIATES WITH VALID MOMENTS*/


*Append the intercept conditions to the matrix;
*Append the time independent conditions to the matrix, if available;
%if &predVarTI. ^=. %then %do;
/*	TypeMtxCombo3 = TypeMtxNB[1,] // TypeMtxIND // TypeMtxCombo3; *Tested intercept;*/
	TypeMtxCombo3 = shape(I(T),1) // TypeMtxIND // TypeMtxCombo3; *Type I intercept;
	
%end;
%else %do;
/*	TypeMtxCombo3 = TypeMtxNB[1,] // TypeMtxCombo3; *Tested intercept;*/
	TypeMtxCombo3 = shape(I(T),1) // TypeMtxCombo3; *Type I intercept;
%end;

*If there are time independent variables, they are represented after the intercept, before the time dependent;
CREATE TypeMtxCombo3 from TypeMtxCombo3;
APPEND from TypeMtxCombo3;
CLOSE TypeMtxCombo3;
print "Each row of TypeMtx is the shifted type vector for each of the predictors, by individual test";
print TypeMtxCombo3;

*Create a printable form of the type matrix;
Type2[,upperIdx] = 0;
Type4out = j(T,T*(Np),.);
DO i=1 to Np;
Type4out[,T*(i-1)+1:T*i] = shape(Type2[i,],T,T);
END;

CREATE Type4out from Type4out;
APPEND from Type4out;
CLOSE Type4out;


*Create the modified data set;
USE mydatasorted;
read all var {&predVarTD.} into xnew;
read all var {&outVar.} into Ynew;
read all var {&idVar. &timeVar.} into othersnew;
%if &predVarTI. ^=. %then %do;
	read all VAR {&idVar.  &timeVar. &outVar. &predVarTI.} INTO otherVars;
%end;
%else %do;
	read all VAR {&idVar.  &timeVar. &outVar.} INTO otherVars;
%end;
CLOSE mydatasorted;


*Creating the new  variables;
X2 = j(N*T,Np*(T-1),0);/*X2 has dimensions [total # obs, #of lagged predictors multiplied by number of lags, all values start at 0*/
DO i=1 TO T-1;
	X2[(i*N+1):N*T, ((i-1)*Np+1):i*Np] = Xnew[1:(T-i)*N,];
	/*from row N+1 to row total number of subjects, from column 1 to Np*/ 
    /*creates matrix with first Np columns containing lag 1 measurements, second Np columns containsing lag2 measurements*/

/*X2=[00|00  first two columns represent lag-1 measurments for both time-dependent covariates measured 3 times
      00|00  
	--------
	  11|00  Second two columns represent lag-2 measurements for both time-dependent covariates measured 3 times
	  11|00
	--------
	  22|11
	  22|11]  */
END;


X2 = Xnew || X2;


*Remove the predictors with no valid moments;
X3 = X2[,keepPred];


*Add in the ID, outcome and time (and possible time independent variables;
X3 = otherVars || X3;


*Create the adjusted variable names;
predNames = t(repeat(t({&predVarTD.}),T));
lagName = j(1,Np*T,0);
DO i=1 to (T-1);
	lagName[1,(i*Np+1):(i+1)*Np]=j(1,Np,i);
END;
varNames = catx("_", predNames, char(lagName));
varnamesTD=varnames;

*Retain only the variable names for those with valid moment conditions;
varNames = varNames[,keepPred];
%if &predVarTI. ^=. %then %do;
	varNames = {&predVarTI.} ||  varNames;
%end;

*Add in the names for the ID, outcome and time;
varNames2 = {&idVar.} || {&timeVar.} || {&outVar.} ||  varNames;

*Sort the data;
call sort(X3 , {1 2});

NcolX3=ncol(X3);

X3[1:N*T,4+NpTI2:NcolX3]=X3[1:N*T,4+NpTI2:NcolX3]#MatrixCreateNew3;


CREATE Mydata3 from X3[c=varNames2];
APPEND from X3;
close Mydata3;

%end;

QUIT;
%mend partitionedDataMatrix; *End macro code;



