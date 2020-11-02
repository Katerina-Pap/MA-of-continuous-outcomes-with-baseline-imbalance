*------------------------------------------------------------------------------
* Import data as reported, incl. missing values 
*------------------------------------------------------------------------------;
data TrowmanNA;
	input study Studyname$ Nobs MeanBaseline sdBaseline seBaseline MeanFU sdFU seFU Correlation 
		  MeanCFB sdCFB seCFB group;
	datalines;
1 Chee     91 56.1 8.9 . . 9.4 . . 0.3 . 0.27 1
2 Jensen   24 94.6 .  2.8 89 . 2.6 . . . . 1
3 Lau      95 56.9 7.1 . . . . . 0.52 . 0.27 1
4 Reid    111 66.0 10.0 . . . . . -0.3 1.8 . 1
5 Shapses1 17 84.1 9.4 . . . . . -7 4.6 . 1
6 Shapses2 11 85.9 9.2 . . . . . -6.7 2.6 . 1
7 Shapses3 18 93.7 13.6 . . . . . -6.7 5.5 . 1
8 Winters  13 57.2 4.9 . 56.3 4.3 . . . . . 1
9 Zemel    13 99.8 4.5 . . . . . -8.58 1.6 . 1
1 Chee     82 57.2 9.4 . 57.3 9.4 . . 0.1 . 0.29 0
2 Jensen   24 93.8 .  2.7 89.1 . 3 . . . . 0
3 Lau      90 58.9 7.5 . . . . . -0.26 . 0.28 0
4 Reid 112 68 11.0 . . . . . -0.09 2.4 . 0
5 Shapses1 19 89.4 10.3 . . . . . -7.3 5.3 . 0
6 Shapses2 11 94.2 15.7 . . . . . -7.6 5.7 . 0
7 Shapses3 24 93.5 14.3 . . . . . -4.3 3.5 . 0
8 Winters  10 54.1 7.2 . 54.8 7.2 . . . . . 0
9 Zemel    14 103.1 6.1 . . . . . -6.6 2.58 . 0
;
proc print;
run;

*------------------------------------------------------------------------------
* Start algebraic calculations and imputations of values  ;
*------------------------------------------------------------------------------;
data Trowman; set TrowmanNA;
* calculate post mean from baseline plus difference;
    if (MeanFU=.)then MeanFU=MeanCFB+MeanBaseline;
* calculate mean change from baseline from mean FU and mean baseline;
    if (MeanCFB=.)then MeanCFB=MeanFU-MeanBaseline;
   
* calculate SD from se;
                if (sdBaseline=.) then
                               sdBaseline=(seBaseline)*sqrt(Nobs);
                if (sdFU=.) then
                               sdFU=seFU*sqrt(Nobs);
                if (sdCFB=.) then
                               sdCFB=(seCFB)*sqrt(Nobs);
 * if missing sd post, we use sd post = sd baseline;
    if (sdFU=.) then
                               sdFU=sdBaseline;
 * calculate correlations; 
	if (Correlation=.)then
                               Correlation=(sdBaseline**2+sdFU**2-sdCFB**2)/(2*sdBaseline*sdFU);
run;

* Calculate median correlation per group;
proc means data = Trowman median;
class group;
var correlation; 
run;

data Trowman; set Trowman;
* impute median correlation for those with missing correlation;
    if (correlation = . & group = 1)then correlation = 0.9367901 ;
	if (correlation = . & group = 0)then correlation = 0.9372800 ;
* calculate sd of change from baseline, using imputed correlations;
	if (sdCFB=.) then sdCFB = sqrt(sdBaseline**2+sdFU**2-2*correlation*sdBaseline*sdFU); 
* calculate seCFB if missing;
	if (seCFB=.) then seCFB = sdCFB/sqrt(nobs);
	drop seBaseline seFU; 
run; 
proc print; run; * outputs the final full dataset 

*------------------------------------------------------------------------------
* Perform (standard) AD approaches, Follow-up analysis, Change scores analysis 
* and recovering ANCOVA methods
*------------------------------------------------------------------------------;

* For the first three methods the data needs to be in wide format;
data trowman1; set Trowman;
if group =1;
rename Nobs = Nobs1 MeanBaseline = MeanBaseline1 sdBaseline=sdBaseline1 MeanFU=MeanFU1 sdFU=sdFU1 Correlation=Correlation1 
       MeanCFB=MeanCFB1 sdCFB=sdCFB1 seCFB=seCFB1;
drop group;
run;
data trowman0; set Trowman;
if group =0;
rename Nobs = Nobs0 MeanBaseline = MeanBaseline0 sdBaseline=sdBaseline0 MeanFU=MeanFU0 sdFU=sdFU0 Correlation=Correlation0 
       MeanCFB=MeanCFB0 sdCFB=sdCFB0 seCFB=seCFB0;
drop group Studyname;
run;

data TrowmanWide; 
 	merge trowman0 trowman1;
	by study;
run;
proc print; run; * outputs the data in wide format 

*------------------------------------------------------------------------------
*   METHOD 1: Follow-up (final) scores approach 
*------------------------------------------------------------------------------;
title "Follow-up scores";
data TrowmanFU; set TrowmanWide;

diffFU = MeanFU1 - MeanFU0;
var_diffFU = sdFU1**2/Nobs1 + sdFU0**2/Nobs0;
run;

data effect_est_FU; set TrowmanFU (keep= study diffFU var_diffFU); run;

* variances of the effect estimates, needed to perform random effects meta-analysis;
data effect_est_FU; set effect_est_FU;
row =_n_; col=_n_;value=var_diffFU; run;
* pool the estimated effects in a random effects meta-analysis;
* the gdata command holds variances at known values;
proc mixed  data = effect_est_FU order=data maxiter=5000 maxfunc=1000 nobound;
class study;
model diffFU= /  solution cl ;
ods output covparms=cov;
random study/gdata=effect_est_FU s; 
repeated diag;
run;

*------------------------------------------------------------------------------
*   METHOD 2: Change scores approach 
*------------------------------------------------------------------------------;
title "Change scores";
data TrowmanCS; set TrowmanWide;

diffCS = MeanCFB1 - MeanCFB0;
var_diffCS = (sdFU1**2/Nobs1) + (sdBaseline1**2/Nobs1) + (sdFU0**2/Nobs0)
             + (sdBaseline0**2/Nobs0) - 2*correlation1*sqrt((sdFU1**2/Nobs1) * (sdBaseline1**2/Nobs1))
             -2*correlation0*sqrt((sdFU0**2/Nobs0)*(sdBaseline0**2/Nobs0));
run;

data effect_est_CS; set TrowmanCS (keep= study diffCS var_diffCS); run;

* variances of the effect estimates, needed to perform random effects meta-analysis;
data effect_est_CS; set effect_est_CS;
row =_n_; col=_n_;value=var_diffCS; run;
proc print;run;

* pool the estimated effects in a random effects meta-analysis;
* the gdata command holds variances at known values;
proc mixed data = effect_est_CS order=data;
class study;
model diffCS= /  solution cl ddfm=kenwardroger;
random study/gdata=effect_est_CS s; 
repeated diag;
run;

*------------------------------------------------------------------------------
*   METHOD 3: Recovering ANCOVA estimates approach
*------------------------------------------------------------------------------;
title "Recovering ANCOVA estimates method ";
data TrowmanANCOVA; set TrowmanWide;

sdpooledB = sqrt((((Nobs1 - 1)*(sdBaseline1**2)) + (Nobs0 - 1)*(sdBaseline0**2))/((Nobs1+Nobs0)-2));
sdpooledF = sqrt((((Nobs1 - 1)*(sdFU1**2)) + (Nobs0 - 1)*(sdFU0**2))/((Nobs1+Nobs0)-2));

ripooled = ((Nobs1*correlation1*sdBaseline1*sdFU1 +  Nobs0*correlation0 *sdBaseline0*sdFU0) ) 
                  /((Nobs1+Nobs0)*sdpooledB*sdpooledF);
             
ancova_est = (MeanFU1-MeanFU0)-ripooled*(sdpooledF/sdpooledB)*(MeanBaseline1-MeanBaseline0);

var_ancova_est =  (sdpooledF**2*(1/Nobs1)+sdpooledF**2*(1/Nobs0))*(1-ripooled**2) ;
se_ancova_est = sqrt(var_ancova_est);
run;

data effect_est_ANCOVA; set TrowmanANCOVA (keep= study ancova_est var_ancova_est); run;

* variances of the effect estimates, needed to perform random effects meta-analysis in SAS;
data effect_est_ANCOVA; set effect_est_ANCOVA;
row =_n_; col=_n_;value=var_ancova_est; run;

* pool the estimated effects in a random effects meta-analysis;
* the gdata command holds variances at known values;
proc mixed data = effect_est_ANCOVA order=data;
class study;
model ancova_est= /  solution cl;
random study/gdata=effect_est_ANCOVA s; 
repeated diag;
run;

*------------------------------------------------------------------------------
*   METHOD 5: Modified Trowman approach 
*------------------------------------------------------------------------------;
title "Modified Trowman";
data TrowmanMod; set TrowmanWide;

diffFU = MeanFU1 - MeanFU0;
var_diffFU = sdFU1**2/Nobs1 + sdFU0**2/Nobs0;
diffBaseline = MeanBaseline1 - MeanBaseline0;
run;
proc print;
run; 

data effect_est_mod; set TrowmanMod (keep= study diffFU var_diffFU diffBaseline); run; 

* variances of the effect estimates, needed to perform random effects meta analysis in SAS;
data effect_est_mod; set effect_est_mod;
row =_n_; col=_n_;value=var_diffFU; run;

* pool the estimated effects in a random effects meta analysis;
* the gdata command holds variances at known values;
proc mixed cl data = effect_est_mod order=data maxiter=5000 maxfunc=1000 nobound;
class study;
model diffFU=  diffBaseline/  solution cl; 
random study/gdata=effect_est_mod s; 
repeated diag;
run;

* Modified Trowman with treatment-baseline interaction*;
data effect_est_mod; set TrowmanMod (keep= study diffFU var_diffFU diffBaseline MeanBaseline1); run; 

* variances of the effect estimates, needed to perform random effects-meta analysis;
data effect_est_mod; set effect_est_mod;
row =_n_; col=_n_;value=var_diffFU; run;

* pool the estimated effects in a random effects meta analysis;
* the gdata command holds variances at known values;
proc mixed cl data = effect_est_mod order=data maxiter=5000 maxfunc=1000 nobound;
class study;
model diffFU=  diffBaseline MeanBaseline1/  solution cl; 
random study/gdata=effect_est_mod s; 
repeated diag;
run;

*------------------------------------------------------------------------------
*   			METHOD 6: Pseudo IPD approach 
* First generate pseudo IPD and then fit mixed effects models
*------------------------------------------------------------------------------;
title "Pseudo IPD approach";

* Generate pseudo IPD;

* Simulate two samples Yi1 and Yi2 of size nobs from N(0,1);
data temp; set Trowman;
	do ID=1 to Nobs;
		ytmp1=rannor(123456);
		ytmp2=rannor(7891011);
		output;
	end;
run;

proc sort data=temp;
by study group; run;

* standardize ytmp1 and ytmp2;
proc standard data=temp mean=0 std=1 out=temp2;
  var ytmp1 ytmp2 ;
  by study group;
run;

* regress ytmp2 on ytmp1, save residuals (ytmp22) and the regression coefficient (which is equal to cor(y1tmp, y2tmp));
ods output ParameterEstimates = parms;
proc reg data=temp2 plots=none;
      model ytmp2=ytmp1/noint;
      output out=temp3 r=ytmp22;
	  by study group;
run;
   
* check that correlation of ytmp1 and ytmp2 by group and study is equal to beta from regresion;   
proc corr data=temp2;
var ytmp1 ytmp2;
by study group;
run;

* generate the pseudo IPD;
data ipd;
	merge temp3 parms(keep = study group estimate);* Add the correlation between ytmp1 and ytmp2(estimate) to the ipd dataset;
	by study group;
	* generate ytmp3 with sd(ytmp3) = 1, cor(ytmp3, ytmp1) = observed correlation; 
	ytmp3 = correlation*ytmp1+ sqrt(1-Correlation*Correlation)*ytmp22/sqrt(1-estimate*estimate); 
	y1 = ytmp1*sdBaseline+MeanBaseline ; * y1 now has mean and sd of original data;
    y2 = ytmp3*sdFU + MeanFU; * y2 has mean and sd of original data ;
	drop ytmp1 ytmp2 ytmp22 ytmp3 estimate;
run;

* a check to see if mean pseudo baselines and mean pseudo outcomes are equal to reported mean baseline per group and mean post baseline outcomes;
proc means data=ipd;
	class study group;
	var y1 y2 Meanbaseline MeanFU sdbaseline sdFU Correlation;
run;

proc corr noprint;
var y1 y2;
by study group; run;


* Some data pre-processing before the LMM ;
data ipd; set ipd;
	arm=100*study+group;
	keep study group arm y1 y2;
run;

* data are centered to prevent convergence problems;
proc means data=ipd nway noprint;
	class study;
	var y1;
	output out=means_y1(drop=_TYPE_ _FREQ_) mean=;
run;

data ipd;
	merge ipd means_y1(rename=(y1=mean_y1));
	by study;
run;

data ipd; set ipd;
y1center = y1 - mean_y1;
meanstudy = mean_y1; *needed for interaction models to separate within- and across-studies variability;
groupc=group;
groupcenter = group-0.5;
run;

*---------------------------------------------------------------
*   One-stage meta ANCOVA using the pseudo IPD and LMM              
*---------------------------------------------------------------;

* Stratified study models;

title "Arm/study specific variances, stratified study model";
proc mixed data=ipd method = reml;
class study arm;
model y2= y1center study y1center*study group/s CL ddfm=kenwardroger;
random groupcenter/subject=study type=un s;
repeated/group=arm;
run;

title "Arm/study specific variances, stratified study model with interaction of baseline with treament"; 
proc mixed data=ipd method = reml;
class study arm;
model y2= y1center study group study*y1center group*y1center group*meanstudy/s cl;
random groupcenter/subject=study type=vc s;
repeated/group=arm;
run;

title "study specific variances, stratified study model";
** study specific variances estimated ;
proc mixed data=ipd method = reml;
class study group;
model y2= y1center study y1center*study group/s cl;
random groupcenter/ subject=study type=un s;
repeated/group=study;
run;

title "study specific variances, stratified study model with interaction of baseline with treament";
** study specific variances estimated;* interaction of baseline with treament added; 
proc mixed data=ipd method = reml;
class study groupcenter;
model y2= y1center study group study*y1center group*y1center group*meanstudy/s cl;
random groupcenter/subject=study type=vc s;
repeated/group=study;
run;
* The variance of the random effect assumed for the interaction of treatment group with baseline was estimated at exactly zero thus it was 
* removed from the model.

title "Group specific variances, stratified study model";
* group specific variances estimated;
proc mixed data=ipd method = reml;
class study groupc;
model y2= y1center study study*y1center group/s cl;
random groupcenter/subject=study type=un s;
repeated/group=groupc;
run;

title "group specific variances, stratified study model with interaction of baseline with treament";
* group specific variances estimated; * interaction of baseline with treament added; 
proc mixed data=ipd method = reml;
class study groupc;
model y2= y1center study group study*y1center group*y1center group*meanstudy/s cl;
random groupcenter /subject=study type=vc s;
repeated/group=groupc;
run;

title "One residual variance, stratified study model";
** one residual variance;
proc mixed data=ipd method = reml;
class study group;
model y2= y1center study study*y1center group/s cl;
random groupcenter/subject=study type=un s;
run;

title "One residual variance, stratified study model with interaction of baseline with treament";
proc mixed data=ipd method = reml;
class study groupc;
model y2= y1center study group study*y1center group*y1center group*meanstudy/s cl;
random groupcenter groupcenter*y1center/subject=study type=vc s;
run;

*-----------------------------------------------------
* Two-stage ANCOVA using the pseudo IPD                              
*-----------------------------------------------------;
title "Two stage ANCOVA - main effect";
proc sort data=ipd; by study;
* separate ancovas per study, regress y2 on y1 and group, save regression coefficients;
ods output ParameterEstimates = effect_est;
proc glm data=ipd plots=none;
		 model y2= y1center group/ solution; *ANCOVA;
		 by study;
run;
		 
data effect_est; set effect_est(keep= study Parameter Estimate StdErr); 
if Parameter="group"; 
run;

* variances of the effect estimates, needed to perform random effects-meta analysis;
data effect_est; set effect_est;
row =_n_; col=_n_;value=StdErr*StdErr; run;

* pool the estimated effects in a random effects meta analysis;
* the gdata command holds variances at known values;
proc mixed data = effect_est order=data;
class study;
model estimate= /  solution cl;
random study/gdata=effect_est s; 
repeated diag;
run;

title "Two stage ANCOVA - interaction effect of treatment by baseline";
ods output ParameterEstimates = effect_est_int;
proc glm data=ipd plots=none;
		 model y2= y1center group y1center*group/ solution; *ANCOVA with interaction;
		 by study;
run;		 
data effect_est_int; set effect_est_int(keep= study Parameter Estimate StdErr); 
if Parameter="y1center*group"; 
run;

* variances of the effect estimates, needed to perform random effects meta-analysis;
data effect_est_int; set effect_est_int;
row =_n_; col=_n_;value=StdErr*StdErr; run;
proc print; run;

* pool the estimated effects in a random effects meta analysis;
* the gdata command holds variances at known values;
proc mixed data = effect_est_int order=data;
class study;
model estimate= /  solution cl ddfm=kenwardroger;
random study/gdata=effect_est_int s; 
repeated diag;
run;
* The variance of random effect was estimated at zero and therefore we fitted an FE model using proc glm and weighting by the inverse variance;

data effect_est_int; set effect_est_int;
inverse_variance = 1/StdErr**2;
run;

proc mixed data = effect_est_int order=data;
class study;
weight inverse_variance;
model estimate= /  solution cl;
run;

* or by using proc glm;
proc glm data = effect_est_int order=data;
class study;
weight inverse_variance;
model estimate= / clparm solution;
run;
