*----------------------------------------------------------------------------------------------------------------------------
*     STATA code supplementing "The impact of trial baseline imbalances should be considered in systematic reviews: 
*     a methodological case study revisited
*     Author: Katerina Papadimitropoulou/Saskia le Cessie
*     Date:   November 2020
*----------------------------------------------------------------------------------------------------------------------------

import excel "Trowman_withNAs.xlsx", sheet("Trowman_withNAs") firstrow clear

*----------------------------------------------------------------------------------------------------------------------------
*                       Start algebraic calculations and imputations of values  
*----------------------------------------------------------------------------------------------------------------------------

* Calculate post baseline mean from CFB and baseline 
replace MeanFU = MeanCFB + MeanBaseline if MeanFU==.

* Calculate change score values from baseline and follow-up
replace MeanCFB = MeanFU - MeanBaseline if MeanCFB==.

* Calculate missing standard deviations from standard errors and vice versa
* Calculate SD from SE 
replace sdBaseline =   seBaseline*sqrt(NCFB) if sdBaseline ==.
replace sdFU        =  seFU*sqrt(NCFB) if sdFU ==.
replace sdCFB       =  seCFB*sqrt(NCFB) if sdCFB ==.

* if not possible assume same SD at baseline and follow-up
replace sdFU = sdBaseline if sdFU == .

* Calculate group correlations using Equation (B8) or impute from reported median correlations of the remaining studies
replace Correlation  = (sdBaseline^2+sdFU^2-sdCFB^2)/(2*sdBaseline*sdFU) if Correlation ==. 

* Calculate median Correlation
egen mediancor=median(Correlation), by(group)

replace Correlation = mediancor if Correlation ==.
drop mediancor

* Final calculations of SE from SD
replace seBaseline     =  sdBaseline/sqrt(NCFB) if seBaseline ==.
replace seFU           =  sdFU/sqrt(NCFB) if seFU ==.
replace sdCFB          =  sqrt(sdBaseline^2+sdFU^2-2*Correlation*sdBaseline*sdFU) if sdCFB ==.
replace seCFB          =  sdCFB/sqrt(NCFB) if seCFB ==.

*----------------------------------------------------------------------------------------------------------------------------
*             Perform (standard) AD approaches, Follow-up analysis, Change scores analysis 
*                                  and recovering ANCOVA methods
*----------------------------------------------------------------------------------------------------------------------------

* For the first three methods the data need to be in wide format 
reshape wide  MeanBaseline sdBaseline seBaseline MeanFU sdFU seFU Correlation MeanCFB sdCFB seCFB NCFB, i(ID) j(group) 

*----------------------------------------------------------------------------------------------
*                             Method 1:  Final scores
*----------------------------------------------------------------------------------------------

meta esize NCFB1 MeanFU1 sdFU1 NCFB0 MeanFU0 sdFU0, esize(mdiff, unequal)
meta summarize, se(khartung)

*----------------------------------------------------------------------------------------------
*                             Method 2:  Change scores
*----------------------------------------------------------------------------------------------

meta esize NCFB1 MeanCFB1 sdCFB1 NCFB0 MeanCFB0 sdCFB0, esize(mdiff, unequal)
meta summarize, se(khartung)

*----------------------------------------------------------------------------------------------
*                             Method 3:  Recovering ANCOVA estimates approach
*----------------------------------------------------------------------------------------------

* calculate pooled standard deviations of baseline and follow-up values
generate sdpooledB =  sqrt(((NCFB1 - 1)*sdBaseline1^2 + (NCFB0 - 1)*sdBaseline0^2)/(NCFB1+NCFB0-2))
generate sdpooledF =  sqrt(((NCFB1 - 1)*sdFU1^2 + (NCFB0 - 1)*sdFU0^2)/(NCFB1+NCFB0-2))

* Calculate ancova estimate using formula from Senn et al. 2007
* using the pooled correlation 

generate ripooled = (NCFB1*Correlation1*sdBaseline1*sdFU1 +  NCFB0*Correlation0 *sdBaseline0*sdFU0)  ///
                 /((NCFB1+NCFB0)*sdpooledB*sdpooledF)

generate ancova_est  =  (MeanFU1-MeanFU0)-ripooled*(sdpooledF/sdpooledB)*(MeanBaseline1-MeanBaseline0)

generate se_ancova_est = sqrt( (sdpooledF^2*(1/NCFB1)+sdpooledF^2*(1/NCFB0))*(1-ripooled^2))
* for different sample sizes from McKenzie and from Senn

meta set ancova_est se_ancova_est
meta summarize, se(khartung)

*----------------------------------------------------------------------------------------------
*                             Method 4:  Trowman approach
*----------------------------------------------------------------------------------------------
* Use data in long format:

reshape long  MeanBaseline sdBaseline seBaseline MeanFU sdFU seFU Correlation MeanCFB sdCFB seCFB NCFB, i(ID) j(group) 

* Fitting a simple linear regression weighted by sample size accounting for treatments clustered in trials (ID)
regress MeanFU MeanBaseline group [fweight=NCFB], vce(cluster ID)

* Alternative fitting of the Trowman methods - preferred over fweight
regress MeanFU MeanBaseline group [pweight=NCFB], vce(cluster ID)

* Trowman method with interaction
regress MeanFU MeanBaseline group c.MeanBaseline##group  [fweight=NCFB], vce(cluster ID)

*----------------------------------------------------------------------------------------------
*                              Method 5:  Modified Trowman approach
*----------------------------------------------------------------------------------------------
* Use data in wide format:

reshape wide  MeanBaseline sdBaseline seBaseline MeanFU sdFU seFU Correlation MeanCFB sdCFB seCFB NCFB, i(ID) j(group) 

generate diff = MeanBaseline0-MeanBaseline1
meta esize NCFB1 MeanFU1 sdFU1 NCFB0 MeanFU0 sdFU0, esize(mdiff, unequal)
meta regress diff , se(khartung)

* Modified Trowman with interaction 
meta regress diff MeanBaseline1 , se(khartung)

*----------------------------------------------------------------------------------------------------------------------------
*                               Method 6:  Pseudo IPD approach
*                                First generate the pseudo IPD
*----------------------------------------------------------------------------------------------------------------------------

* Use data in long format:
reshape long  MeanBaseline sdBaseline seBaseline MeanFU sdFU seFU Correlation MeanCFB sdCFB seCFB NCFB, i(ID) j(group) 

* Generate the pseudo baselines and outcomes, make dataset with as many lines as observations
expand NCFB

set seed 123456
generate ytmp1 = rnormal(0,1)
generate ytmp2 = rnormal(0,1)

* Standardize ytmp1 and ytmp2 by study and group
egen mean_ytmp1 = mean(ytmp1), by(ID group)
egen sd_ytmp1 = sd(ytmp1), by(ID group)
replace ytmp1 = (ytmp1 - mean_ytmp1) / sd_ytmp1
egen mean_ytmp2 = mean(ytmp2), by(ID group)
egen sd_ytmp2 = sd(ytmp2), by(ID group)
replace ytmp2 = (ytmp2 - mean_ytmp2) / sd_ytmp2
* check
*by ID group: tabstat ytmp1 ytmp2, stat(mean sd)

* Correlation between y1tmp and y2tmp
bysort ID group : egen corrtmp = corr(ytmp1 ytmp2)
* Check
by ID group: tabstat corrtmp

* Calculate residual of regression of ytmp2 on ytmp1
gen resid = ytmp2-corrtmp*ytmp1
* Generate ytmp3 with sd(ytmp3) = 1, cor(ytmp3, ytmp1) = observed correlation
generate ytmp3 = ytmp1*Correlation + sqrt(1-Correlation^2)*resid/sqrt(1-corrtmp^2)
* Generate pseudo baseline and pseudo follow-up outcomes
generate y1    = ytmp1*sdBaseline + MeanBaseline
* y1 now has mean and sd of original data
generate y2    = ytmp3*sdFU + MeanFU
 *y2 has mean and sd of original data

* Checks
bysort group:tabstat y1 y2,  statistics( mean sd )  by(ID)
bysort group ID: corr(y1 y2)

drop ytmp1 ytmp2 mean_ytmp1 sd_ytmp1 mean_ytmp2 sd_ytmp2 corrtmp resid ytmp3

* Pre-step to calculate centered baseline values by study
egen meany1bystudy =  mean(y1), by(ID)
gen y1center       =  y1 - meany1bystudy
gen groupcenter    =  group - 0.5
gen arm            =  1000*ID + group

*----------------------------------------------------------------------------------------------
*                       One-stage pseudo IPD models 
*----------------------------------------------------------------------------------------------
 
 gen nobs= _n
 gen group0 = 1-group
 
*-----------------------------------------------------------------------------
* Study stratified intercept and random treatment effect ANCOVA - Main effect
*----------------------------------------------------------------------------- 

*arm/study specific residual variances estimated 
mixed y2 group y1center i.ID  i.ID##c.y1center|| ID: groupcenter, covariance(unstructured) resid(ind, by(arm)) noconstant reml dfmethod(satterthwaite)

*study specific residual variances estimated
mixed y2 group y1center i.ID  i.ID##c.y1center|| ID:groupcenter, covariance(unstructured) resid(ind, by(ID)) noconstant reml dfmethod(satterthwaite)

*group specific residual variances estimated
mixed y2 group y1center i.ID  i.ID##c.y1center|| ID:groupcenter, covariance(unstructured) resid(ind, by(group)) noconstant reml dfmethod(satterthwaite)

*one residual variance estimated estimated 
mixed y2 group y1center i.ID  i.ID##c.y1center|| ID: groupcenter, var reml noconstant dfmethod(satterthwaite)

*-----------------------------------------------------------------------------------
* Study stratified intercept and random treatment effect ANCOVA - Interaction effect
*----------------------------------------------------------------------------------- 

*arm/study specific residual variances with treatment-baseline interaction 
mixed y2 group y1center i.ID  i.ID##c.y1center group##c.y1center group##c.meany1bystudy|| ID: groupcenter, covariance(unstructured) resid(ind, by(arm)) noconstant reml dfmethod(satterthwaite)

*study specific residual variances with treatment-baseline interaction
mixed y2 group y1center i.ID  i.ID##c.y1center group##c.y1center group##c.meany1bystudy|| ID:groupcenter, covariance(unstructured) resid(ind, by(ID)) noconstant reml dfmethod(satterthwaite)

*group specific residual variances with treatment-baseline interaction
mixed y2 group y1center i.ID  i.ID##c.y1center group##c.y1center group##c.meany1bystudy|| ID:groupcenter, covariance(unstructured) resid(ind, by(group)) noconstant reml dfmethod(satterthwaite) 

*one residual variance estimated with treatment-baseline interaction
mixed y2 group y1center i.ID  i.ID##c.y1center group##c.y1center group##c.meany1bystudy|| ID: groupcenter, var reml noconstant dfmethod(satterthwaite)

*----------------------------------------------------------------------------------------------
*                       Two-stage pseudo IPD models 
*----------------------------------------------------------------------------------------------

*-----------------------------------------------------------------------------
*  Pooling study-level estimates of treatment effect:  Main effect
*----------------------------------------------------------------------------- 

* Preserve the pseudo IPD 
preserve 
* ANCOVA per study on pseudo IPD for subsequent two-stage MA
* Store coefficients 
sort ID group 
statsby _b _se, by(ID) clear: regress y2 y1center group

meta set  _b_group _se_group
meta summarize,se(khartung)

*-----------------------------------------------------------------------------
*  Pooling study-level estimates of interaction effect:  Interaction effect
*----------------------------------------------------------------------------- 
* Ancovas per study with interaction of baseline and treatment effect on pseudo IPD for subsequent two-stage MA

* restore the pseudo IPD
restore, preserve
 
* ANCOVA per study on pseudo IPD for subsequent two-stage MA
* Store coefficients 
sort ID group 
statsby _b _se, by(ID) clear: regress y2 c.group##c.y1center

meta set  _b_group _se_group
meta summarize,se(khartung)

*interaction terms (have funny names)
meta set  _stat_3 _stat_7
meta summarize,se(khartung)
