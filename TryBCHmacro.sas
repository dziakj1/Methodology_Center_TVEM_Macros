OPTIONS MPRINT MLOGIC;
%INCLUDE "C:\Users\jjd264\Documents\BCHmacro\LCA_Regression_BCH.sas";

PROC IMPORT OUT= WORK.SimData 
            DATAFILE= "C:\Users\jjd264\Documents\BCHmacro\SimpleExample.txt" 
            DBMS=TAB REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN; 

PROC LCA DATA=SimData  
         OUTPOST=Noninclusive_Post  ;
    NCLASS 3;
    ITEMS I1 I2 I3 I4 I5;
    CATEGORIES 2 2 2 2 2;  
    SEED 11111;
    NSTARTS 5;
    ID id X1 X2;
    RHO prior=1; 
RUN;

%LCA_Regression_BCH(postprobs = Noninclusive_Post,
    id = id, 
    Covariates = X1 X2);

/*
%LCA_Regression_BCH(postprobs = Noninclusive_Post,
    id = id, 
    Covariates = X0 X1 X2,
    automatically_add_intercept = 0);

*/
 
