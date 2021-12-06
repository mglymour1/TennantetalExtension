/* replicating Tennant Mediator DAG and simulations  */

/* call with:
CallResponseTables.do
This code is modified to reproduce Tennant et al case 3a, allowing for a potential outcomes estimation
In this code, compared to the paper: WC is X; IC0 is Y0; IC1 is Y1
*/
clear
*set seed 7878
set obs 1000

* set means and variances: from Tennant
local WCmu  = 9.5       
local WCvar = 1.6^2
local IC0mu = 4.0
local IC0var= 0.74^2
*local IC1mu = 4.2     /*IC1 mean not constrained in this setting*/
*local IC1var= 0.74^2  /*IC1 total variance not constrained in this setting*/

* Create a variable "slope" representing biological change in IC 
local slopemu=0
local slopevar=1

*local pIC0_IC1 =0.65   /*correlation of Y0 and Y1 is not an input */
local betaIC0_IC1=1     /*effect of Y0 on Y1 is 1 */
local betaWC_slope  = 0.05 /*key assumption in this simulation */
local pWC0_IC0 = 0.5 
local betaWC0_IC0=`pWC0_IC0'*(`IC0var'^.5)/(`WCvar'^.5)
*local pWC0_IC1 =.433  /*as specified in some scenarios */
*local pWC0_IC1 = `EffSize'-(`pIC0_IC1'*`pWC0_IC0')  /*not constrained in this scenario- derives from effect of WC on slope and IC0*/
*local pWC0_IC1 = 0 /*modified 3a to show the problem w/ ancova */
*local betaWC0_IC1=`pWC0_IC1'*(`IC1var'^.5)/(`WCvar'^.5)  /*not constrained in this scenario- derives from effect of WC on slope and IC0*/
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

gen slope_e=rnormal()*`slopevar'^.5
gen slope=`slopemu'+(WC-`WCmu')*`betaWC_slope'+slope_e
gen IC1=IC0*`betaIC0_IC1'+slope /*IC1 is IC0 plus change*/

gen IC1_error=rnormal()*`error_var'

gen IC1_obs=IC1+IC1_error
gen slope_po_WC_0=`slopemu'+(0-`WCmu')*`betaWC_slope'+slope_e
gen slope_po_WC_1=`slopemu'+(1-`WCmu')*`betaWC_slope'+slope_e
gen IC1_poWC_0=IC0_poWC_0*`betaIC0_IC1'+slope_po_WC_0 /*IC1 is IC0 plus change*/
gen IC1_poWC_0_obs=IC1_poWC_0+IC1_error  /* observed potential outcome for iC1 settings WC to 0 */
gen IC1_poWC_1=IC0_poWC_1*`betaIC0_IC1'+slope_po_WC_1 /*IC1 is IC0 plus change*/
gen IC1_poWC_1_obs=IC1_poWC_1+IC1_error /* observed potential outcome for iC1 settings WC to 1 */

gen IC1_poWC_0_Y0_mean=`IC0mu'*`betaIC0_IC1'+slope_po_WC_0  /*Y0 doesn't affect slope*/
gen IC1_poWC_0_Y0_mean_obs=IC1_poWC_0_Y0_mean+IC1_error
gen IC1_poWC_1_Y0_mean=`IC0mu'*`betaIC0_IC1'+slope_po_WC_1  /*Y0 doesn't affect slope*/
gen IC1_poWC_1_Y0_mean_obs=IC1_poWC_1_Y0_mean+IC1_error

summarize
corr IC0 IC1
scalar corr_IC0_IC1=r(rho)
/*Estimands */
gen tru_po_diff=(IC1_poWC_1-IC0_poWC_1)-(IC1_poWC_0-IC0_poWC_0)   /*estimand for potential outcome contrast for effect of X on change */
sum tru_po_diff
scalar true_po_diff=r(mean)

gen CDE_est=IC1_poWC_1_Y0_mean-IC1_poWC_0_Y0_mean               /*estimand for potential outcome contrast for effect of X on Y1 setting Y0 to its mean */
sum CDE_est
scalar CDE=r(mean)

/* row 3 */
reg IC1 WC IC0            /*ANCOVA estimate */
matrix betaANCOVA=e(b)

gen dy=IC1-IC0           /*Change score estimate */
reg dy WC
matrix betaCHANGE=e(b)

reg IC1 WC /* sanity check, : total estimated effect of X on Y1  */
matrix betaY1onX=e(b) 

scalar betaCHANGE=betaCHANGE[1,1]
scalar betaANCOVA=betaANCOVA[1,1]
scalar betaY1onX=betaY1onX[1,1]

/* row 4 */
reg IC1_obs WC IC0_obs             /*ANCOVA estimate with measurement error*/
matrix betaANCOVA_obs=e(b)

gen dy_obs=IC1_obs-IC0_obs         /*Change score estimate with measurement error*/
reg dy_obs WC
matrix betaCHANGE_obs=e(b)

reg IC1_obs WC
matrix betaY1onX_obs=e(b) 

scalar betaCHANGE_obs=betaCHANGE_obs[1,1]
scalar betaANCOVA_obs=betaANCOVA_obs[1,1]
scalar betaY1onX_obs=betaY1onX_obs[1,1]

/*Below not used- the estimand should be the same regardless of the measurement error */
gen po_diff_obsvd=(IC1_poWC_1_obs-IC0_poWC_1_obs)-(IC1_poWC_0_obs-IC0_poWC_0_obs)
sum po_diff_obsvd
scalar po_diff_obs=r(mean)

gen CDE_obsvd=IC1_poWC_1_Y0_mean_obs-IC1_poWC_0_Y0_mean_obs
sum CDE_obsvd
scalar CDE_obs=r(mean)
sum
