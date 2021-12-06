# TennantetalExtension
/* This code is to replicate work presented in Tennant PWG, Arnold KF, Ellison GTH, et al.
/* Analyses of ‘change scores’ do not estimate causal effects in observational data. International Journal of Epidemiology. 2021;dyab050. 
/* This simulation mirrors Tennant et al Table 1, Scenario 3A. Tennant's code was in R and used DAGGITY but this is recreated in Stata.
/* I additionally wished to calculate the target estimands by contrasting potential outcomes          */
/* for the controlled direct effect of X on Y1 setting Y0 to a specific value and the effect of X on change from Y0 to Y1 */

/* To create the potential outcomes contrast (i.e., replicate the simulation while setting X or Y0 to specific values),
/* my construction of the data generating model is slightly different. 

/* There are 3 data simulation do files, each called with "CallResponseTables.do", which outputs .csv files and results corresponding
/* to the Table presented in my response to Tennant. 

/* One version of the simulation (Tennant_IJE_3A_showeffectofY0onChange.do) calculates potential outcomes varying the value of Y0, to evaluate whether         */
/* Y0 influences change in Y under Tennant's assumptions. WC=X, IC0=Y0, IC1=Y1. "PO" and "po" indicate potential outcomes */
/* The other two versions of the simulation calculate rows 1 and 2 of the Table in my response and rows 3 and 4 of the Table in my response, respectively.
