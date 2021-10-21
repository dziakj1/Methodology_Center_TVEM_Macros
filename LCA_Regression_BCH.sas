   /*-----------------------------------------------------------------*
   | The documentation and code below is supplied by The Methodology              
   | Center, Penn State University.             
   *------------------------------------------------------------------*                          
  /*------------------------------------------------------------------*
   | MACRO NAME  : %LCA_Regression_BCH version 0.1.0 
   | CREATED BY  : John J. Dziak (October 7, 2019)  
   | CONTACT :  Methodology Center Helpdesk (MChelpdesk@psu.edu)
   *------------------------------------------------------------------*
   | Copyright 2020, The Pennsylvania State University
   |
   | This software is distributed "as is" without any warranty.   
   | All warranties, including, but not limited to, implied 
   | warranties of merchantability and fitness for a particular.
   | purpose, are disclaimed.
   *------------------------------------------------------------------*/   
 
%MACRO LCA_Regression_BCH(postprobs,
        id, 
        Covariates,
        Groups=,
        Assignment=Modal,
        Adjustment=BCH,
        ref_class = 1,
        sampling_weight=,
        automatically_add_intercept = 1
    );    
    TITLE "Using LCA_Regression_BCH Macro";
    %IF %TRIM("&groups")="" %THEN %DO;
        %LET there_are_groups = 0;
    %END; %ELSE %DO;
        %LET there_are_groups = 1; 
    %END;      
    %IF %EVAL(&there_are_groups=1) %THEN %DO; 
        DATA BCHMacroAnswers; RUN;
        /* Sort the data according to group, and separate into multiple datasets. */
        PROC IML;
            USE &postprobs;
                READ ALL VAR {&Groups} INTO Group;
            CLOSE &postprobs;
        nGroups = MAX(Group);
        CALL SYMPUT("nGroups",CHAR(nGroups));
                   /* the total number of groups; */  
        %DO this_group=1 %TO &nGroups;            
            DATA _subset_post_group&this_group;
                SET &postprobs;
                IF &Groups=&this_group;
            RUN; 
            %LCA_Regression_BCH_No_Groups(postprobs=_subset_post_group&this_group,
                id=&id, 
                Covariates=&Covariates, 
                Assignment=&Assignment,
                Adjustment=&Adjustment,
                ref_class=&ref_class,
                sampling_weight=&sampling_weight,
                automatically_add_intercept=&automatically_add_intercept
            );   
            DATA BCHMacroAnswersThisGroup; 
                SET BCHMacroAnswersThisGroup;
                Group = &this_group;
            RUN;
            DATA BCHMacroAnswers;
                SET BCHMacroAnswers BCHMacroAnswersThisGroup;
            RUN;                  
        %END;    
    %END; %ELSE %DO;
        %LCA_Regression_BCH_No_Groups(postprobs=&postprobs,
            id=&id, 
            Covariates=&Covariates, 
            Assignment=&Assignment,
            Adjustment=&Adjustment,
            ref_class=&ref_class,
            sampling_weight=&sampling_weight,
            automatically_add_intercept=&automatically_add_intercept
        );   
        DATA BCHMacroAnswers;
            SET BCHMacroAnswersThisGroup;
        RUN;                                  
    %END;
    DATA BCHMacroAnswers;
        SET BCHMacroAnswers;
        WHERE Class IS NOT MISSING;
    RUN;
    TITLE "Estimates from LCA_Regression_BCH Macro";
    PROC PRINT DATA=BCHMacroAnswers NOOBS;
    RUN;
    TITLE;
%MEND LCA_Regression_BCH;

%MACRO LCA_Regression_BCH_No_Groups(postprobs,
    id, 
    Covariates,
    Assignment,
    Adjustment,
    ref_class,
    sampling_weight,
    automatically_add_intercept
    );     
    /******************************************************************************/
    /***** Count the number of classes in the model using the postprobs dataset ***/
    /******************************************************************************/  
    DATA bchtemp1; SET &postprobs; RUN;    
    PROC TRANSPOSE DATA=&postprobs OUT=TemporaryDatasetForBCHmacro;
    RUN;
    PROC IML;
        USE TemporaryDatasetForBCHmacro;
            READ ALL VAR {_NAME_} INTO name;
        CLOSE TemporaryDatasetForBCHmacro;
        nClasses = (SUM(SUBSTR(UPCASE(name),1,6)="POSTLC")); /* number of classes */
        CALL SYMPUT("nClasses",CHAR(nClasses));
    QUIT; 
    PROC DATASETS NOLIST NOWARN NOPRINT;
        DELETE TemporaryDatasetForBCHmacro;
    RUN;
    %put _all_;
    /******************************************************************************/
    /***** Construct BCH-adjusted class weights to represent class membership *****/
    /******************************************************************************/    
    %CalculateClassWeights( PostProbsDataset=&postprobs,
                            SubjectID=&id,
                            nClasses=&nClasses, 
                            Sampling_Weight=&sampling_weight,
                            Assignment=&Assignment,  
                            Adjustment=&Adjustment);
            /* Outputs datasets ClassWeights and ClassWeightsLong. */
    /*******************************************************************************/
    /***** Append class weights to the dataset for analysis. ***********************/
    /*******************************************************************************/
    DATA ForAnalysis;
        MERGE &postprobs ClassWeights;
        BY &id; 
    RUN;
    /*******************************************************************************/
    /***** Add an intercept column to the covariates if appropriate ****************/
    /*******************************************************************************/   
    %IF %EVAL(&automatically_add_intercept>0) %THEN %DO;
        DATA ForAnalysis;
            SET ForAnalysis;
            _Intercept_ = 1;
        RUN;
        %LET Covariates = _Intercept_ &Covariates;
    %END; 
    /*******************************************************************************/
    /***** Create long form of dataset with one row per class per participant. *****/
    /*******************************************************************************/
    DATA ForAnalysisLong;
        SET ForAnalysis; 
        %DO class = 1 %TO &nClasses;
            LC = &class;    ClassWeight = ClassWeightLC&class; OUTPUT; 
        %END;
        DROP %DO class = 1 %TO &nClasses;    
                ClassWeightLC&class 
            %END;
        ;
    RUN;
    /*******************************************************************************/
    /***** Estimate the parameters for the underlying parametric weighted **********/
    /***** multinomial generalized logistic regression.                   **********/
    /*******************************************************************************/     
    %DoWeightedGLogit(InputDataset=ForAnalysisLong,
                        OutcomeClass=LC,
                        ref_class=&ref_class,
                        Weight=ClassWeight,
                        Covariates=&Covariates, 
                        SubjectID=&id,
                        ObservationID=LC );  
            /* Outputs datasets WeightedGLogitCoef and WeightedGLogitSandwichCov. */
    PROC IML;
        USE WeightedGlogitCoef;
            READ ALL;
        CLOSE WeightedGlogitCoef; 
        CovariateIndex = Covariate;
        CovariateNamesString = "&Covariates";
        NumCov = 1 + COUNT(CovariateNamesString," ");
        CovariateName = J(NROW(CovariateIndex),1,"                                  "); 
        DO i = 1 TO NumCov;
            CovariateName[LOC(CovariateIndex=i)] = SCANQ(CovariateNamesString,i);
        END; 
        z_statistic = Estimate / StandardError;
        pvalue = 1 - PROBNORM(ABS(z_statistic));
        CREATE BCHMacroAnswersThisGroup VAR {CovariateName Class Estimate StandardError z_statistic pvalue};
            APPEND;
        CLOSE BCHMacroAnswersThisGroup;
    QUIT; 
%MEND LCA_Regression_BCH_No_Groups;
%MACRO DoWeightedGLogit(InputDataset=,
                    OutcomeClass=,
                         ref_class=1,
                    Weight=,
                    Covariates=,
                    SubjectID=,
                         ObservationID=,
                    RegularizationConstantExp=0,
                    RegularizationConstantForCov=0,
                         RegularizationConstantForSE=0);
    PROC IML;
        USE &InputDataset;
             READ ALL VAR {&OutcomeClass} INTO Y;
            READ ALL VAR {&Weight} INTO W;
            READ ALL VAR {&Covariates} INTO XMatrix;
            READ ALL VAR {&SubjectID} INTO id;
        CLOSE &InputDataset;
        numCovariates  = NCOL(XMatrix);
        numClasses = Y[<>]; 
        CALL SYMPUT("numCovariates",(CHAR(numCovariates)));
        CALL SYMPUT("numClasses",(CHAR(numClasses)));
        DataForNlmixedAsMatrix = id || W || Y || XMatrix;
        VariableNames = "id" || "W" || "Y" || ("X1":"X&numCovariates");  
        CREATE DataForNlmixed FROM DataForNlmixedAsMatrix [COLNAME=VariableNames];
            APPEND FROM DataForNlmixedAsMatrix ;
        CLOSE DataForNlmixed;
    QUIT;    
    ODS EXCLUDE ALL; RUN; 
    ODS NORESULTS; RUN; 
    PROC NLMIXED DATA=DataForNlmixed COV START HESS METHOD=GAUSS TECH=NRRIDG MAXITER=500 
MAXFUNC=1000; 
        PARMS %DO k = 1 %TO %EVAL(&numCovariates);
                %DO c = 1 %TO %EVAL(&numClasses);
                         %IF %EVAL(&c ^= &Ref_Class) %THEN %DO;
                         b&k&c = 0  
                         %END;
                %END;
              %END; ;
        %DO c = 1 %TO %EVAL(&numClasses);
               %IF %EVAL(&c ^= &Ref_Class) %THEN %DO;
            eta&c = b1&c*X1 
                    %IF %EVAL(&numCovariates>1) %THEN %DO;
                        %DO  k = 2 %TO &numCovariates;
                            + b&k&c*X&k 
                        %END;
                    %END;
                  ;
            expEta&c = EXP(eta&c + &RegularizationConstantExp);
               %END; %ELSE %DO;
             eta&c = 0;
               %END;
        %END;
        denominator = 1 
                %DO c = 1 %TO %EVAL(&numClasses);
                         %IF %EVAL(&c ^= &Ref_Class) %THEN %DO;
                          + expEta&c 
                    %END;
            %END; ; 
        %DO c = 1 %TO %EVAL(&numClasses);
               %IF %EVAL(&c ^= &Ref_Class) %THEN %DO;
                 IF Y=&c THEN likelihood = expEta&c/denominator;
               %END; %ELSE %DO;
                  IF Y=&c THEN likelihood = 1/denominator;
               %END;
        %END;
          SumSquaredBetas = 0;
        %DO c = 1 %TO %EVAL(&numClasses);
               %IF %EVAL(&c ^= &Ref_Class) %THEN %DO;
               SumSquaredBetas = SumSquaredBetas + (b1&c)**2 
                    %IF %EVAL(&numCovariates>1) %THEN %DO;
                        %DO  k = 2 %TO &numCovariates;
                            + (b&k&c)**2
                        %END;
                    %END;;
               %END;
        %END;
        loglikelihood = W *LOG(likelihood) ;
        MODEL Y ~ GENERAL(loglikelihood); 
        ODS OUTPUT PARAMETERESTIMATES=WeightedGlogitCoef
                     FITSTATISTICS=WeightedGlogitFitStatistics
                     CONVERGENCESTATUS=WeightedGlogitConvergenceStatus; 
    RUN;  
    PROC PRINT DATA=WeightedGlogitCoef; RUN;
    DATA WeightedGlogitCoef;
        SET WeightedGlogitCoef;
        Row = _N_;
        Covariate = INT((Row-1)/(&numClasses-1))+1;
        Class = Row-(Covariate-1)*(&numClasses-1);
          IF (Class=>&ref_class) THEN Class = Class + 1; /*skip reference class */
        DROP Row;
        StandardError = .;
        KEEP Covariate Class Estimate StandardError;
    RUN;
    PROC PRINT DATA=WeightedGlogitCoef; RUN;
    PROC SORT DATA=WeightedGlogitCoef;
        BY Class Covariate;
    RUN;
    PROC PRINT DATA=WeightedGlogitCoef; RUN;    
    DATA TempDatasetForMacro; SET &InputDataset;
        BY &SubjectID;
        RETAIN BCHMacroIntID 0;
        IF first.&SubjectID THEN BCHMacroIntID = BCHMacroIntID + 1;
            /* This creates a recoded version of the subject ID variable,
               called "BCHMacroIntID," that is composed of consecutive
               integers starting with 1.  This makes it possible for some
               of the code in the macro to be simpler than if subjects could
               have anything for ID's.
             See http://www.ats.ucla.edu/stat/sas/faq/enumerate_id.htm; */
        OUTPUT; 
    RUN;    
    PROC IML;
         /* Might be able to improve memory efficiency using the FREE command 
           or the SYMSIZE and WORKSIZE options? */
        USE  WeightedGlogitCoef;
            READ ALL VAR {Covariate Class Estimate};
        CLOSE WeightedGlogitCoef;
          USE TempDatasetForMacro;
            READ ALL VAR {&OutcomeClass} INTO Y;
            READ ALL VAR {&Weight} INTO W;
            READ ALL VAR {&Covariates} INTO XMatrix;
            READ ALL VAR {BCHMacroIntID} INTO ClusterId;
               READ ALL VAR {&ObservationID} INTO MemberID;
        CLOSE TempDatasetForMacro;
        p  = NCOL(XMatrix); 
        ntotal = NROW(Y); 
        nsub = ClusterID[<>];  /* assumes indexing starts at 1 for this variable */
        PRINT ClusterID;
        PRINT ntotal nsub;
        m = ntotal / nsub;  /* assumes equal cluster sizes */
        n = ntotal / m; /* assumes equal cluster sizes */
        nc = Y[<>];  /* Number of classes.   */
        YAsMatrix = J(n*m,nc,0); 
        DO c = 1 TO nc; 
            IF SUM(Y=c)>0 THEN DO;
                whichrows = LOC(Y=c)`; 
                YAsMatrix[whichrows,c] = 1;
            END; 
        END; 
        CoefAsVector = Estimate;
        ModelInformationMatrix = J((nc-1)*p,(nc-1)*p,0);
        TotalScoreVector = J((nc-1)*p,1,0);
        EmpiricalCovOfScore = J(p*(nc-1),p*(nc-1),0);
        CoefAsMatrix = SHAPE(CoefAsVector,nc-1)`;
        ref_class = &ref_class;          
        IF (ref_class = 1) THEN DO;
             Nonref_classes = 2:nc;
        END; 
        IF (ref_class = nc) THEN DO;
             Nonref_classes = 1:(nc-1);
        END; 
        IF ((ref_class >1) & (ref_class < nc)) THEN DO;
            Nonref_classes = (1:(ref_class-1)) || ((ref_class+1):nc);
        END; 
        DO i = 1 TO n;
            ScoreVector = J((nc-1)*p,1,0);
            DO j = 1 TO m;                 
                wij = w[LOC(ClusterID=i & MemberID=j)]; * Assumes there is only one such entry. ;
                yij = YasMatrix[LOC(ClusterID=i & MemberID=j),]`;
                xij = XMatrix[LOC(ClusterID=i & MemberID=j),]; 
                eta = J(nc,1,0);
                eta[Nonref_classes] =  xij*CoefAsMatrix;
                eta = (eta <> -20) >< +20;
                eta = eta + &RegularizationConstantForSE;
                expEta = EXP(eta);
                piij = expEta / expEta[+]; 
                Mij = DIAG(piij[Nonref_classes]) - piij[Nonref_classes]*piij[Nonref_classes]`;                
                BigXij = (I(nc-1) @ xij)`;
                ScoreVector = ScoreVector + wij*BigXij*(yij[Nonref_classes]-piij[Nonref_classes]);
                ModelInformationMatrix = ModelInformationMatrix + wij*BigXij*Mij*BigXij`;
            END; 
            TotalScoreVector = TotalScoreVector + ScoreVector;
            EmpiricalCovOfScore = EmpiricalCovOfScore + ScoreVector*ScoreVector`;
        END;         
        CoefAsVectorOld = CoefAsVector;  
        lambda = &RegularizationConstantForcov*(VECDIAG(ModelInformationMatrix))[+];
          NaiveCovCoef = INV(ModelInformationMatrix+lambda*I(NROW(ModelInformationMatrix))); 
        CoefAsVector = CoefAsVector + NaiveCovCoef*TotalScoreVector;
        MaxAbsDev = (ABS(CoefAsVector - CoefAsVectorOld))[<>];
        CorrectionFactor = (n/(n-1))*((n*m-1)/(n*m-(nc-1)*p));  
        * In SAS PROC SURVEYLOGISTIC, the option VADJUST=NONE uses n/(n-1).  
        *   The option VADJUST=DF uses (n/(n-1))*((n*m-1)/(n*m-(nc-1)*p)). ; 
        SandwichCov = NaiveCovCoef` * (CorrectionFactor*EmpiricalCovOfScore) * NaiveCovCoef; 
        IF (VECDIAG(SandwichCov)[><]>0) THEN DO;          
            StandardError  = SQRT(VECDIAG(SandwichCov)); 
          END; ELSE DO;
              StandardError = J(NROW(SandwichCov),1,.);
          END;
        EDIT WeightedGlogitCoef;
            REPLACE ALL VAR {StandardError };
        CLOSE WeightedGlogitCoef;
        CREATE WeightedGlogitSandwichCov FROM SandwichCov; 
            APPEND FROM SandwichCov;
        CLOSE WeightedGlogitSandwichCov;
    QUIT;
    ODS EXCLUDE NONE; RUN;
    ODS RESULTS; RUN;
%MEND DoWeightedGLogit;
%MACRO CalculateClassWeights(PostProbsDataset,
                             SubjectID,   
                             nClasses,
                             Sampling_Weight,
                             Assignment,  /* modal or proportional */
                             Adjustment  /* BCH or none */
                             );
    PROC IML;     
         Assignment = LOWCASE(TRIM("&Assignment"));
         Adjustment = LOWCASE(TRIM("&Adjustment"));
          IF (Assignment^="modal" & Assignment^="proportional") THEN DO;
               PRINT ("Assignment option not recognized:");
               PRINT(Assignment);
          END;
          IF (Adjustment^="bch" & Adjustment^="none") THEN DO;
               PRINT ("Adjustment option not recognized:");
               PRINT(Adjustment);
          END;
        /* Calculate BCH weights.*/    
        USE &PostProbsDataset;
            READ ALL VAR {&SubjectID 
                            %DO class = 1 %TO &nClasses;
                                PostLC&class
                            %END;
                          &sampling_weight Best};
        CLOSE &PostProbsDataset;
        PostLC = PostLC1     %DO class = 2 %TO &nClasses;
                        || PostLC&class
                    %END;  ;  /* concatenate the vectors into a matrix */
          IsBestLC = (Best=1)      %DO class = 2 %TO &nClasses;
                              || (Best=&class) 
                         %END; ;
          UnadjustedWeights = J(NROW(IsBestLC),&nClasses,.);
          IF (Assignment="modal") THEN DO;
               UnadjustedWeights = IsBestLC;
          END; ELSE DO;
               UnadjustedWeights = PostLC;
          END;                         
          IF (LENGTH(LOWCASE(TRIM("_&sampling_weight")))>1) THEN DO; 
                DO class = 1 TO &nClasses;
                    UnadjustedWeights[,class] = (&sampling_weight + 0)#UnadjustedWeights[,class];
                            /* adding zero is a weird work-around to get around a macro language limitation*/ 
               END;
          END; ELSE DO; 
          END;
          IF (Adjustment="bch") THEN DO;
                D = J(&nClasses,&nClasses,.);  
               /* Based on Bolck, Croon and Hagenaars (2004, Political
                    Analysis), Vermunt (2010, Political Analysis) and Bakk and Vermunt (2016,
                    Structural Equation Modeling), we will define D so that D(t,s) is the 
                    probability of being assigned to class s, for a subject whose
                    true class is t. */
               DO t = 1 TO &nClasses;
                    DO s = 1 TO &nClasses;
                         D[t,s] = ((PostLC[,t] # UnadjustedWeights[,s])[+]) / (PostLC[+,t]);
                    END;
               END;  
               ClassWeights = UnadjustedWeights * INV(D);  
               /* These are the BCH weights.  */
          END; ELSE DO;
               ClassWeights = UnadjustedWeights;
          END;  
          %DO class = 1 %TO &nClasses;
               ClassWeightLC&class = ClassWeights[,&class];
          %END; 
        CREATE ClassWeights VAR {&SubjectID %DO class = 1 %TO &nClasses;
                ClassWeightLC&class  
                %END;};
            APPEND ;
        CLOSE ClassWeights;
    QUIT; 
     DATA ClassWeightsLong;
          SET ClassWeights;
         %DO class = 1 %TO &nClasses;
               Class = &class;
               ClassWeight = ClassWeightLC&class;
               OUTPUT;
          %END;    
          KEEP &SubjectID Class ClassWeight;
     RUN;
%MEND CalculateClassWeights;


