/* replicating Tennant Mediator DAG and simulations  */
/* This simulation mirrors Tennant et al Table 1, Scenario 3A, but additionally calculates the target estimands           */
/* for the controlled direct effect of X on Y1 setting Y0 to a specific value and the effect of X on change from Y0 to Y1 */
/* This version of the simulation also calculates potential outcomes varying the value of Y0, to evaluate whether         */
/* Y0 influences change in Y under Tennant's assumptions. WC=X, IC0=Y0, IC1=Y1. "PO" and "po" indicate potential outcomes */
/* call with:
CallResponseTables.do
*/
clear
set obs 1000


* set means and variances: from Tennant as much as possible
local WCmu  = 9.5       
local WCvar = 1.6^2
local IC0mu = 4.0
local IC0var= 0.74^2
local IC1mu = 4.2
local IC1var= 0.74^2       /*constraining variance of IC1 to be same as variance of IC0 entails effects of IC0 on change*/

local pIC0_IC1 =0.65       /*rho IC0 and IC1*/
local betaIC0_IC1=`pIC0_IC1'*(`IC1var'^.5)/(`IC0var'^.5)
local EffSize  = 0.433     /*total effect of X on Y1, or WC on IC1, including via IC0 and via change */
local pWC0_IC0 = 0.5 
local betaWC0_IC0=`pWC0_IC0'*(`IC0var'^.5)/(`WCvar'^.5)
local pWC0_IC1 = `EffSize'-(`pIC0_IC1'*`pWC0_IC0')  /*specified for scenario 3A- effect of X on Y1 subtracting product of path coefficients from X to Y0 and Y0 to Y1*/
*local pWC0_IC1 = 0 /*modified 3a to show the problem w/ ancova */
local betaWC0_IC1=`pWC0_IC1'*(`IC1var'^.5)/(`WCvar'^.5)

* Generate data: in order to calculate the potential outcomes, variables are generated in two pieces, the purely random component and the 
* 'fixed' component that is determined by parent variables 					
gen WC_e=rnormal()*(`WCvar'^.5)           /*random component of WC*/
gen WC=`WCmu'+WC_e                        /*WC is its mean plus each individual's random component */
gen IC0_f=`IC0mu'+(WC-`WCmu')*`betaWC0_IC0'        /*predictable or 'fixed' component of IC0: its mean plus the effect of WC on IC0*/
gen IC0_e=rnormal()*(`IC0var'*(1-`pWC0_IC0'^2))^.5 /*random component of IC0: set so the variance of the sum of fixed and random components matches target */
gen IC0=IC0_f+IC0_e                                /*IC0 is the sum of fixed and random components*/
gen IC0_poWC_0=`IC0mu'+(0-`WCmu')*`betaWC0_IC0' +IC0_e       /*IC0 potential outcome setting WC to 0 */
gen IC0_poWC_1=`IC0mu'+(1-`WCmu')*`betaWC0_IC0' +IC0_e       /*IC0 potential outcome setting WC to 1 */
gen IC1_f=`IC1mu'+(IC0-`IC0mu')*`betaIC0_IC1'+(WC-`WCmu')*`betaWC0_IC1' /*predictable or 'fixed' component of IC1 */
sum IC1_f                                                    /* need this line to calculate the variance of the random component */
gen IC1_e=rnormal()*(`IC1var'-r(Var))^.5                     /* random component of IC1*/
gen IC1_e_WC_0=rnormal()*(`IC1var'-r(Var))^.5  /* random component of IC1 if WC=0*/              /* we make 3 version of this random component but it is not necessary- results the same but w/o uncertainty if you just use 1 random component*/
gen IC1_e_WC_1=rnormal()*(`IC1var'-r(Var))^.5  /* random component of IC1 if WC=1*/

gen IC1=IC1_f+IC1_e /*IC1 is the sum of fixed and random components */
gen IC1_poWC_0=`IC1mu'+(IC0_poWC_0-`IC0mu')*`betaIC0_IC1'+(0-`WCmu')*`betaWC0_IC1'+IC1_e_WC_0    /*IC1 potential outcome setting WC to 0 */
gen IC1_poWC_1=`IC1mu'+(IC0_poWC_1-`IC0mu')*`betaIC0_IC1'+(1-`WCmu')*`betaWC0_IC1'+IC1_e_WC_1    /*IC1 potential outcome setting WC to 1 */

gen IC1_poWC_0_Y0_mean=`IC1mu'+`IC0mu'*`betaIC0_IC1'+(0-`WCmu')*`betaWC0_IC1'+IC1_e_WC_0   /*IC1 potential outcome setting WC to 0 and IC0 to its popln mean*/
gen IC1_poWC_1_Y0_mean=`IC1mu'+`IC0mu'*`betaIC0_IC1'+(1-`WCmu')*`betaWC0_IC1'+IC1_e_WC_1   /*IC1 potential outcome setting WC to 1 and IC0 to its popln mean*/

gen IC1_poWC_0_Y0_p1=`IC1mu'+(`IC0mu'+1)*`betaIC0_IC1'+(0-`WCmu')*`betaWC0_IC1'+IC1_e_WC_0  /*IC1 PO setting WC to 0 and IC0 to its popln mean+1*/
gen IC1_poWC_1_Y0_p1=`IC1mu'+(`IC0mu'+1)*`betaIC0_IC1'+(1-`WCmu')*`betaWC0_IC1'+IC1_e_WC_1  /*IC1 PO setting WC to 1 and IC0 to its popln mean+1*/

* Does this setup entail that IC0 influences change in IC?
* change in IC if IC0=mean and WC=0
sum IC1_poWC_0_Y0_mean
dis r(mean)-`IC0mu'

* change in IC if IC0=mean+1 and WC=0
sum IC1_poWC_0_Y0_p1
dis r(mean)-`IC0mu'+1

* change in IC if IC0=mean and WC=1
sum IC1_poWC_1_Y0_mean
dis r(mean)-`IC0mu'
scalar effectofY0onchange_WC_0=r(mean)-`IC0mu'

* change in IC if IC0=mean+1 and WC=1
sum IC1_poWC_1_Y0_p1
dis r(mean)-`IC0mu'+1
scalar effectofY0onchange_WC_1=r(mean)-`IC0mu'+1

/* Below is just to confirm it delivers consistent with Tennant's estimates */
* Estimate true effects by contrasting potential outcomes
summarize
corr
gen tru_po_diff=(IC1_poWC_1-IC0_poWC_1)-(IC1_poWC_0-IC0_poWC_0)  /* True potential outcome difference in change in IC from 0 to 1 setting WC to 1 vs 0*/
sum tru_po_diff
scalar true_po_diff=r(mean)

gen CDE_est=IC1_poWC_1_Y0_mean-IC1_poWC_0_Y0_mean           /* True CDE on IC1 from setting WC to 1 vs 0 but holding IC0 to its mean*/
sum CDE_est
scalar CDE=r(mean)

reg IC0 WC

* Compare ANCOVA based estimates and CHANGE models
gen dy=IC1-IC0
reg dy WC
matrix betaCHANGE=e(b)             /* estimated effect of WC on change score*/
reg IC1 WC IC0
matrix betaANCOVA=e(b)             /* estimated NDE (should equal CDE)*/
reg IC1 WC                          
matrix betaY1onX=e(b)              /* estimated total effect of WC on IC1*/
scalar betaCHANGE=betaCHANGE[1,1]
scalar betaANCOVA=betaANCOVA[1,1]
scalar betaY1onX=betaY1onX[1,1]
