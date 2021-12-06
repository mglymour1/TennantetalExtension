# TennantetalExtension
This code is to replicate work presented in Tennant PWG, Arnold KF, Ellison GTH, et al.
"Analyses of ‘change scores’ do not estimate causal effects in observational data". International Journal of Epidemiology. 2021;dyab050. 

This simulation mirrors Tennant et al Table 1, Scenario 3A. Tennant's code was in R and used DAGGITY but this is recreated in Stata.
I additionally wished to calculate the target estimands by contrasting potential outcomes       
for the controlled direct effect of X on Y1 setting Y0 to a specific value and the effect of X on change from Y0 to Y1.

To create the potential outcomes contrast (i.e., replicate the simulation while setting X or Y0 to specific values),
my construction of the data generating model is slightly different than Tennant's, and i do not use DAGGItY, but it preserves
his parameters in the initial scenario. 

There are 3 data simulation do files, each called with "CallResponseTables.do", which outputs .csv files and results corresponding
to the Table presented in my response to Tennant.  
In all, WC=X, IC0=Y0, IC1=Y1. "PO" and "po" indicate potential outcomes

One version of the data generation simulation (Tennant_IJE_3A_showeffectofY0onChange.do) calculates potential outcomes setting the value of Y0 to different
values, to evaluate whether Y0 influences change in Y under Tennant's assumptions. I made this claim in the text and wanted to confirm here.

The other two versions of the simulation calculate rows 1 and 2 of the Table in my response (row 1 precisely mirrors Tennant,
row 2 adds measurement error on Y0 and Y1) and rows 3 and 4 (row 3 examines a situation in which the direct effect and the effect on change
are identical, because Y0 has no effect on change, and row 4 adds measurement error to that situation) of the Table in my response.

Simulations match DAGs shown in my response paper: 

![image](https://user-images.githubusercontent.com/16063827/144901472-2398f16a-d353-47ca-a4fd-293320c6bc29.png)
