OPTIONS MPRINT MLOGIC;
%INCLUDE "C:\Users\jjd264\Documents\BCHmacro\LCA_Regression_BCH.sas";

PROC IMPORT OUT= OurSimData 
            DATAFILE= "C:\Users\jjd264\Documents\BCHmacro\GroupsExample.txt" 
            DBMS=TAB REPLACE;
     GETNAMES=YES;
     DATAROW=2; 
RUN; 

PROC LCA DATA=OurSimData  
         OUTPOST=Noninclusive_Post  ;
    NCLASS 3;
    ITEMS I1 I2 I3 I4 I5;
    CATEGORIES 2 2 2 2 2;  
    SEED 11111;
    NSTARTS 5;
    ID id X1 X2;
    RHO prior=1; 
	GROUPS G;
RUN;

 
%LCA_Regression_BCH(postprobs = Noninclusive_Post,
    id = id, 
    Covariates = X1 X2,
    Groups = G);
