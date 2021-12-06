/* This code calls each of the data simulations, stores the results, and reports the 25th and 97.5th bounds */
/* Demonstrate that Tennant's specification entails that Y0 affects change in Y */
set seed 977
simulate betaCHANGE betaANCOVA betaY1onX true_po_diff CDE effectofY0onchange_WC_0 effectofY0onchange_WC_1, r(1000):do Tennant_IJE_3A_showeffectofY0onchange.do 
rename _sim_1 betaCHANGE
rename _sim_2 betaANCOVA
rename _sim_3 betaY1onX
rename _sim_4 true_po_diff
rename _sim_5 true_CDE
rename _sim_6 effectofY0onchange_WC_0
rename _sim_7 effectofY0onchange_WC_1
sort betaCHANGE
list if _n==25
list if _n==975
sort betaANCOVA
list if _n==25
list if _n==975
sort betaY1onX
list if _n==25
list if _n==975
sort true_po_diff
list if _n==25
list if _n==975
sort true_CDE
list if _n==25
list if _n==975
sort effectofY0onchange_WC_0
list if _n==25
list if _n==975
sort effectofY0onchange_WC_1
list if _n==25
list if _n==975
save  "ShowEffectofY0onchange.dta" , replace

/* Row 1 and 2 of Table 1, replicating Tennant et al and adding (row 2) measurement error to Y0 and Y1 */
set seed 977
clear
simulate betaCHANGE betaANCOVA betaY1onX true_po_diff CDE betaCHANGE_obs /*
*/       betaANCOVA_obs betaY1onX_obs po_diff_obs CDE_obs corr_IC0_IC1 /*
*/       , r(1000):do Tennant_IJE_3A_row1_2.do 
rename _sim_1 betaCHANGE
rename _sim_2 betaANCOVA
rename _sim_3 betaY1onX
rename _sim_4 true_po_diff
rename _sim_5 true_CDE
rename _sim_6 betaCHANGE_obs
rename _sim_7 betaANCOVA_obs
rename _sim_8 betaY1onX_obs
rename _sim_9 po_diff_obs
rename _sim_10 CDE_obs

local varlist betaANCOVA betaCHANGE betaY1onX true_po_diff true_CDE  betaANCOVA_obs betaCHANGE_obs betaY1onX_obs po_diff_obs CDE_obs

foreach i in `varlist' {
sort `i'
dis "*******************************************************************"
mean `i'
dis "2.5TH PERCENTILE LOWER BOUND FOR " 
dis "                                  `i'"
list `i' if _n==25
dis "97.5TH PERCENTILE UPPER BOUND FOR " 
dis "                                   `i'"
list `i' if _n==975
dis "******************************************************************"
}
save rows1_2_k1000v3.csv , replace

/* Row 3 and 4 of Table 1, a setting where CDE and effect on change are identical, without (row 3) and with (row 4) measurement error on Y0 and Y1 */

set seed 977
clear
simulate betaCHANGE betaANCOVA betaY1onX true_po_diff CDE betaCHANGE_obs betaANCOVA_obs betaY1onX_obs po_diff_obs CDE_obs , r(1000):do Tennant_IJE_3A_row3_4.do 
rename _sim_1 betaCHANGE
rename _sim_2 betaANCOVA
rename _sim_3 betaY1onX
rename _sim_4 true_po_diff
rename _sim_5 true_CDE
rename _sim_6 betaCHANGE_obs
rename _sim_7 betaANCOVA_obs
rename _sim_8 betaY1onX_obs
rename _sim_9 po_diff_obs
rename _sim_10 CDE_obs

local varlist  betaANCOVA betaCHANGE betaY1onX true_po_diff true_CDE  betaANCOVA_obs betaCHANGE_obs betaY1onX_obs po_diff_obs CDE_obs

foreach i in `varlist' {
sort `i'
dis "********************************"
mean `i'
dis "2.5TH PERCENTILE LOWER BOUND FOR " 
dis "`i'"
list `i' if _n==25
dis "97.5TH PERCENTILE UPPER BOUND FOR " 
dis "`i'"
list `i' if _n==975
}
save rows3_4_k1000v4.csv , replace

