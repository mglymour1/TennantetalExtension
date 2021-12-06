/* replicating Tennant Mediator DAG and simulations  */
/* This simulation mirrors Tennant et al Table 1, Scenario 3A, but additionally calculates the target estimands           */
/* for the controlled direct effect of X on Y1 setting Y0 to a specific value and the effect of X on change from Y0 to Y1 */
/* Ro1 1 replicates Tennant, row 2 add measurement error to observed variables. WC=X, IC0=Y0, IC1=Y1. "PO" and "po" indicate potential outcomes */
/* call with:
CallResponseTables.do
This code is modified to reproduce Tennant et al case 3a, allowing for a potential outcomes estimation
*/
clear
*set seed 7878
set obs 1000

* set means and variances: from Tennant
local WCmu  = 9.5       
local WCvar = 1.6^2
local IC0mu = 4.0
local IC0var= 0.74^2
local IC1mu = 4.2
local IC1var= 0.74^2

local pIC0_IC1 =0.65   /*correlation of Y0 and Y1 */
local betaIC0_IC1=`pIC0_IC1'*(`IC1var'^.5)/(`IC0var'^.5)
local EffSize  = 0.433
local pWC0_IC0 = 0.5 
local betaWC0_IC0=`pWC0_IC0'*(`IC0var'^.5)/(`WCvar'^.5)
*local pWC0_IC1 =.433  /*as specified in some scenarios */
local pWC0_IC1 = `EffSize'-(`pIC0_IC1'*`pWC0_IC0')  /*specified for scenario 3A*/
*local pWC0_IC1 = 0 /*modified 3a to show the problem w/ ancova */
local betaWC0_IC1=`pWC0_IC1'*(`IC1var'^.5)/(`WCvar'^.5)
/* Settings to add measurement error */
local error_var=.5
						
gen U=rnormal() /*not used in scenario 3A*/
gen WC_e=rnormal()*(`WCvar'^.5)
gen WC=`WCmu'+WC_e
gen IC0_f=`IC0mu'+(WC-`WCmu')*`betaWC0_IC0'  /*predictable or 'fixed' component of IC0 */
gen IC0_e=rnormal()*(`IC0var'*(1-`pWC0_IC0'^2))^.5 /*random component of IC0 */
gen IC0_error=rnormal()*`error_var'  /* error term to add to observed variables */
gen IC0=IC0_f+IC0_e/*IC0 is the sum of fixed and random components*/
gen IC0_obs=IC0+IC0_error /*actual observed IC0 */
gen IC0_poWC_0=`IC0mu'+(0-`WCmu')*`betaWC0_IC0' +IC0_e
gen IC0_poWC_0_obs=IC0_poWC_0+IC0_error /*observed potential outcome setting WC to 0 */
gen IC0_poWC_1=`IC0mu'+(1-`WCmu')*`betaWC0_IC0' +IC0_e
gen IC0_poWC_1_obs=IC0_poWC_1+IC0_error /* observed potential outcome setting WC to 1 */

gen IC1_f=`IC1mu'+(IC0-`IC0mu')*`betaIC0_IC1'+(WC-`WCmu')*`betaWC0_IC1' /*predictable or 'fixed' component of IC1 */
sum IC1_f
gen IC1_e=rnormal()*(`IC1var'-r(Var))^.5       /* random component of IC1*/
gen IC1_e_WC_0=rnormal()*(`IC1var'-r(Var))^.5  /* random component of IC1 if WC=0*/
gen IC1_e_WC_1=rnormal()*(`IC1var'-r(Var))^.5  /* random component of IC1 if WC=1*/
gen IC1_error=rnormal()*`error_var'

gen IC1=IC1_f+IC1_e /*IC1 is the sum of fixed and random components */
gen IC1_obs=IC1+IC1_error
gen IC1_poWC_0=`IC1mu'+(IC0_poWC_0-`IC0mu')*`betaIC0_IC1'+(0-`WCmu')*`betaWC0_IC1'+IC1_e_WC_0
gen IC1_poWC_0_obs=IC1_poWC_0+IC1_error  /* observed potential outcome for iC1 settings WC to 0 */
gen IC1_poWC_1=`IC1mu'+(IC0_poWC_1-`IC0mu')*`betaIC0_IC1'+(1-`WCmu')*`betaWC0_IC1'+IC1_e_WC_1
gen IC1_poWC_1_obs=IC1_poWC_1+IC1_error /* observed potential outcome for iC1 settings WC to 1 */

gen IC1_poWC_0_Y0_mean=`IC1mu'+`IC0mu'*`betaIC0_IC1'+(0-`WCmu')*`betaWC0_IC1'+IC1_e_WC_0
gen IC1_poWC_0_Y0_mean_obs=IC1_poWC_0_Y0_mean+IC1_error
gen IC1_poWC_1_Y0_mean=`IC1mu'+`IC0mu'*`betaIC0_IC1'+(1-`WCmu')*`betaWC0_IC1'+IC1_e_WC_1
gen IC1_poWC_1_Y0_mean_obs=IC1_poWC_1_Y0_mean+IC1_error

summarize
corr IC0 IC1
scalar corr_IC0_IC1=r(rho)
/* row 1 - estimands*/
gen tru_po_diff=(IC1_poWC_1-IC0_poWC_1)-(IC1_poWC_0-IC0_poWC_0)  /*estimand for potential outcome contrast for effect of X on change */
sum tru_po_diff
scalar true_po_diff=r(mean)

gen CDE_est=IC1_poWC_1_Y0_mean-IC1_poWC_0_Y0_mean      /*estimand for potential outcome contrast for effect of X on Y1 setting Y0 to its mean */
sum CDE_est
scalar CDE=r(mean)

/* row 1 */
reg IC1 WC IC0                /*ANCOVA estimate */
matrix betaANCOVA=e(b)

gen dy=IC1-IC0
reg dy WC                     /*Change score estimate */
matrix betaCHANGE=e(b)

reg IC1 WC /* sanity check to confirm consistency with Tennant: total estimated effect of X on Y1 */
matrix betaY1onX=e(b) 

scalar betaCHANGE=betaCHANGE[1,1]
scalar betaANCOVA=betaANCOVA[1,1]
scalar betaY1onX=betaY1onX[1,1]

/* row 2 */
reg IC1_obs WC IC0_obs         /*ANCOVA estimate with measurement error*/
matrix betaANCOVA_obs=e(b)

gen dy_obs=IC1_obs-IC0_obs
reg dy_obs WC                  /*Change score estimate with measurement error*/
matrix betaCHANGE_obs=e(b)

reg IC1_obs WC
matrix betaY1onX_obs=e(b)
scalar betaCHANGE_obs=betaCHANGE_obs[1,1]
scalar betaANCOVA_obs=betaANCOVA_obs[1,1]
scalar betaY1onX_obs=betaY1onX_obs[1,1]

gen po_diff_obsvd=(IC1_poWC_1_obs-IC0_poWC_1_obs)-(IC1_poWC_0_obs-IC0_poWC_0_obs)  /*potential outcomes of observed variables are not actually relevant- estimand is the latent variable */
sum po_diff_obsvd
scalar po_diff_obs=r(mean)

gen CDE_obsvd=IC1_poWC_1_Y0_mean_obs-IC1_poWC_0_Y0_mean_obs  /*potential outcomes of observed variables are not actually relevant- estimand is the latent variable */
sum CDE_obsvd
scalar CDE_obs=r(mean)
sum
