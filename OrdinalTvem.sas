%MACRO OrdinalTvem(dataset=, /* The name of the SAS dataset 
                                containing the input */
                   id=,   /* The name of the subject ID variable */
                   t=,    /* The name of the observation time variable */
                   cov=,  /* The name or names of the predictor 
				            variables with time-invariant effects,
                            (separated by spaces, not commas, if 
				            there are more than one)   */
                   tcov=, /* The name or names of the predictor 
				             variables with time-varying effects,
                             (separated by spaces, not commas, if 
				             there are more than one)   */
                   y=,    /* The name of the outcome variable */
                   /* Note: The measurement times do not have to be evenly spaced but
                            it is assumed that there are no missing data indicators (.)
                            in any of the variables above */
				   ComputationOption=, /* Which computational method to request 
				                           from SAS PROC GLIMMIX for fitting a 
				                           generalized linear mixed model.
				                           1 for doubly iterative pseudolikelihood (RSPL) 
				                             method with NOBOUND option
				                           2 for doubly iterative pseudolikelihood (RSPL) 
				                             method without NOBOUND option
				                           3 for singly iterative likelihood approximated
				                             with Gauss-Hermite quadrature 
				                       Method 3 is not implemented unless UseRandom > 0 and
				                       UsePropOdds = 1.  In those cases we have a generalized
				                       linear model or linear mixed model, which are simpler
				                       than a generalized linear mixed model.  In the case
				                       of UseRandom > 0 and UsePropOdds = 0, i.e., a normal
				                       linear mixed model, the only difference between options
				                       1 and 2 is the NOBOUND option.  If option 3 is requested
				                       in that case, it is treated as option 1. */
                   MinTToPlot=, /* The first time point to plot on the output graph */
                   MaxTToPlot=, /* The last time point to plot on the output graph */
                   UsePropOdds=1, /* Either 0 or 1.  If 0 is specified, a linear
				                     TVEM will be fit.  If 1 is specified, a 
				                     proportional odds TVEM will be fit. */
                   UseRandom=0,  /* 0 or 1, or 2.  If 1 is specified, a random
				                    intercept will be included.  If 0 is specified,
				                    no random effects are included.  If 2 is specified,
                                    a random intercept and slope are included. */
                   DoPlots=1, /* Either 0 or 1.  If 1 is specified, a plot will
				                    be drawn to show the time-varying coefficient 
				                   over time. */
                   deg=2,  /* 1 for a linear, 2 for a quadratic, 
				              or 3 for a cubic spline */
                   NumInteriorKnots=3 /* The number of interior knots (1 or more */
                   );

    /***************************************************************
    | OrdinalTvem macro Version 2.0
    * By John DZIAK, Runze LI, and Anne BUU
    * Fits a TVEM model to (potentially longitudinal) ordinal data.
    *
    * Uses code which was written by Xianming TAN to construct the spline bases
	* and which is also found in the MixTVEM macro.
	*
    * Copyright:
    * (c) 2014 The Pennsylvania State University
    *
    * License:
    * This program is free software; you can redistribute it and/or
    * modify it under the terms of the GNU General Public License as
    * published by the Free Software Foundation; either version 2 of
    * the License, or (at your option) any later version.
    * This program is distributed in the hope that it will be useful,
    * but WITHOUT ANY WARRANTY; without even the implied warranty of
    * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
    * General Public License for more details.
    *
    * Acknowledgments and references:
    * We fit a varying-coefficient (time-varying effect) model for
    * ordinal data, using either a linear model or a proportional odds
    * logistic model, and an unpenalized B-spline approach.  See
    *    Eilers, P.H.C. and Marx, B.D. (1996). Flexible smoothing using
    *        B-splines and penalized likelihood. Statistical Science
    *        11(2): 89-121.
    *    Hastie, T., & Tibshirani, R. (1993). Varying-coefficient models. Journal
    *        of the Royal Statistical Society, Series B, 55, 757-796.
    *    Shiyko, M. P., Lanza, S. T., Tan, X., Li, R., Shiffman, S. (2012). Using
    *        the Time-Varying Effect Model (TVEM) to Examine Dynamic Associations
    *        between Negative Affect and Self Confidence on Smoking Urges:
    *        Differences between Successful Quitters and Relapsers. Prevention
    *        Science, 13, 288-299.
    *    Ramsay, J., Hooker, G., & Graves, S. (2009). Functional Data Analysis with
    *        R and MATLAB. New York: Springer.
    *    SAS Institute Inc. (2011). SAS/STAT (c) 9.3 user's guide: The GLIMMIX 
    *        procedure (chapter)). Cary, NC: SAS Institute Inc.
    *    Tan, X., Shiyko, M. P., Li, R., Li, Y., & Dierker, L. (2011, November 21).
    *        A time-varying effect model for intensive longitudinal data.
    *        Psychological Methods. Advance online publication. doi: 10.1037/a0025814.
    **************************************************************/
    %LET debug = 0;
    %IF %EVAL(&UsePropOdds = 0) %THEN %DO;
        %LET ModelType = Linear;
    %END; %ELSE %DO;
        %IF %EVAL(&UsePropOdds = 1) %THEN %DO;
            %LET ModelType = ProportionalOdds;
        %END; %ELSE %DO;
            %PUT Error: unrecognized UsePropOdds argument;
        %END;
    %END;
    %IF %EVAL(&UseRandom = 0) %THEN %DO;
        %LET ModelType = &ModelType;
    %END; %ELSE %DO;
        %IF %EVAL(&UseRandom = 1) %THEN %DO;
            %LET ModelType = &ModelType (random intercept);
        %END; %ELSE %DO;
            %IF %EVAL(&UseRandom = 2) %THEN %DO;
                %LET ModelType = &ModelType (random intercept, slope);
            %END; %ELSE %DO;
                %PUT Error: unrecognized UseRandom argument;
            %END;
        %END;
    %END;
    PROC IML;
        START do_Bspline(x, knots, d);
            n_ext = ncol(knots);
            all_knots = (1:(n_ext+2));
            all_knots[1] = knots[1] - (1e-12);
            all_knots[2:(n_ext+1)]= knots;
            all_knots[n_ext+2] = knots[n_ext] + (1e-12);
            n_ext =n_ext+2;
            tmp=0.0*(1:(n_ext-1));
            DO i = 1 TO (n_ext-1);
                IF (x>=all_knots[i] & x<all_knots[i+1]) THEN tmp[i]=1.0;
            END;
            j=1;
            DO WHILE (j<=d);
                  /* the De Boor formula, refer to Eilers and Marx (1996) */
                DO i = 1 TO (n_ext-j-1);
                    w1 = (x-all_knots[i])/(all_knots[i+j] - all_knots[i]);
                    w2 = (all_knots[i+j+1]-x)/(all_knots[i+j+1] - all_knots[i+1]);
                    tmp[i] = w1*tmp[i] +w2*tmp[i+1];
                END;
                j=j+1;
            END;
            RETURN(tmp[1:(n_ext-d-1)]);
        FINISH do_Bspline;
        START vec_Bspline(xx, knots, d);
            n_row = ncol(xx);
            n_col = ncol(knots)+2-d-1;
            out = J(n_row, n_col, .);
            DO i=1 TO n_row;
                out[i, ] = t(do_Bspline(xx[i], knots, d));
            END;
            RETURN(out);
        FINISH vec_Bspline;
        START GenerateTimeBasis(TimeBasis, Time, InteriorKnots, deg);
            * B-Spline basis;
            dx = InteriorKnots[1] - Time[><];
            EarlyKnots = InteriorKnots[><]-dx*((deg):1);
            LateKnots = InteriorKnots[<>]+dx*(1:(deg));
            AllKnots = EarlyKnots || InteriorKnots || LateKnots;
            TimeBasis = vec_Bspline(Time`, AllKnots, deg);
        FINISH GenerateTimeBasis;
        /* Start the main part of the macro */
        USE &dataset;
            READ ALL VAR {&t} INTO t;
            CALL SYMPUT("NumInvariantEffects",CHAR(0));
            %IF %EVAL(%LENGTH(%SUPERQ(cov))>0) %THEN %DO;
                READ ALL VAR {&cov} INTO cov;
                CALL SYMPUT("NumInvariantEffects",CHAR(NCOL(cov)));

            %END;
            %IF %EVAL(%LENGTH(%SUPERQ(tcov))>0) %THEN %DO;
                READ ALL VAR {&tcov} INTO tcov;
            %END;
        CLOSE &dataset;
        MinTime = t[><];
        MaxTime = t[<>]; 
		MeanTime = t[:];
        dt = (MaxTime - MinTime) / (&NumInteriorKnots + 1);
        InteriorKnots = (MinTime + dt*(1:(&NumInteriorKnots))); 
        %IF %EVAL(%LENGTH(%SUPERQ(tcov))>0) %THEN %DO;
            tcov = J(NROW(t),1,1) || tcov;
        %END; %ELSE %DO;
            tcov = J(NROW(t),1,1);
        %END;
        %IF %EVAL(&NumInteriorKnots > 0) %THEN %DO;
            CALL GenerateTimeBasis(TimeBasis,
                           t,
                           InteriorKnots,
                           &deg);
        %END; %ELSE %DO;
            TimeBasis = J(NROW(t),1,1);
            DO d = 1 TO (&deg);
                TimeBasis = TimeBasis || t##d;
            END;
        %END;
        CALL SYMPUT("NumColumns",CHAR(NCOL(TimeBasis)));
        CALL SYMPUT("MeanTime",CHAR(MeanTime));
        CALL SYMPUT("NumTCovsIncludingIntercept",CHAR(NCOL(tcov)));
         
        CALL SYMPUT("NumColumnsTotal",CHAR(NCOL(TimeBasis) * 
                                           &NumTCovsIncludingIntercept  + 
                                           &NumInvariantEffects));
        DO thisTCov = 1 TO (&NumTCovsIncludingIntercept);
            covariateTimesTimeBasis = TCov[,thisTCov]#TimeBasis;
            IF (thisTCov = 1) THEN design = covariateTimesTimeBasis;
            ELSE design = design || covariateTimesTimeBasis; 
        END;
        whichBeta = ((1:&NumTCovsIncludingIntercept)@J(1,&NumColumns,1))`;
        whichThetaWithinBeta = (J(1,&NumTCovsIncludingIntercept,1)@(1:&NumColumns))`;
        %IF %EVAL(&NumInvariantEffects > 0) %THEN %DO;
            design = design || cov;
            whichBeta = whichBeta // J(&NumInvariantEffects,1,0);
            whichThetaWithinBeta = whichThetaWithinBeta // J(&NumInvariantEffects,1,0);
        %END;
        CREATE TCov FROM TCov;
            APPEND FROM TCov; 
        CLOSE TCov;
        CREATE TimeBasis FROM TimeBasis; 
            APPEND FROM TimeBasis; 
        CLOSE TimeBasis;
        CREATE Design FROM Design
                   [COLNAME=(COMPRESS(CONCAT("B",CHAR(1:(&NumColumnsTotal)))))];
            APPEND FROM Design;
        CLOSE Design;
        CREATE WhichBeta FROM whichBeta;
            APPEND FROM whichBeta;
        CLOSE WhichBeta;
        CREATE WhichThetaWithinBeta FROM whichThetaWithinBeta;
             APPEND FROM whichThetaWithinBeta;
        CLOSE WhichThetaWithinBeta;
        CALL SYMPUT("LastB","B%QTRIM(&NumColumnsTotal)");
    QUIT;
    DATA TvemWorkingDataset;
        MERGE &dataset Design;
    RUN;
    DATA TheseFitStats; RUN;
    %IF %EVAL(&UseRandom = 2) %THEN %DO;
		DATA TvemWorkingDataset; 
	        SET TvemWorkingDataset; 
	        tvem_centered_time = &t - &MeanTime; 
	    RUN;
	%END;
    TITLE "Intermediate Calculations for TVEM Model";
    PROC GLIMMIX MAXOPT=5000 NOCLPRINT 
        %IF %EVAL(&ComputationOption = 1) %THEN %DO; NOBOUND %END;
        %IF %EVAL((&ComputationOption = 3)&(&UsePropOdds = 1)) %THEN %DO; METHOD=QUAD %END;
		/* ComputationOption = 2 follows the GLIMMIX defaults */
        DATA=TvemWorkingDataset;
    %IF %EVAL(&UseRandom > 0) %THEN %DO;
        CLASS &id;
    %END;
    %IF %EVAL(&UsePropOdds = 0) %THEN %DO;
        MODEL &y = B1 - &LastB /
                   SOLUTION COVB(DETAILS) NOINT;
    %END; %ELSE %DO;
        MODEL &y(DESCENDING) = B2 - &LastB /
                   SOLUTION COVB(DETAILS) LINK=CUMLOGIT DIST=MULTINOMIAL;
        NLOPTIONS TECH=NRRIDG;
                  /* see support.sas.com/resources/papers/proceedings12/332-2012.pdf */
    %END;
        OUTPUT OUT=DataWithFittedByLevel PRED=predicted;
    %IF %EVAL(&UseRandom = 0) %THEN %DO;
        ODS OUTPUT ParameterEstimates=ParameterEstimates
                   CovB=CovB
                   FitStatistics=TheseFitStats;
    %END;
    %IF %EVAL(&UseRandom = 1) %THEN %DO;
        RANDOM INTERCEPT / SUB = &id;
        ODS OUTPUT ParameterEstimates=ParameterEstimates
                   CovB=CovB
                   CovParms=CovParms
                   FitStatistics=TheseFitStats
                   SolutionR=SolutionR ;
	%END;
    %IF %EVAL(&UseRandom = 2) %THEN %DO;
        RANDOM INTERCEPT tvem_centered_time / SUB = &id;
        ODS OUTPUT ParameterEstimates=ParameterEstimates
                   CovB=CovB
                   CovParms=CovParms
                   FitStatistics=TheseFitStats
                   SolutionR=SolutionR ;
	%END;
    RUN;
    TITLE;
    PROC TRANSPOSE DATA=TheseFitStats OUT=TheseFitStats; ID Descr; RUN; 
    %IF %EVAL(&UsePropOdds = 1) %THEN %DO;
         DATA DataWithFittedByLevel;
             SET DataWithFittedByLevel;
             CALL SYMPUT("SomeLevel",_Level_);
         RUN;
    %END;
    DATA DataWithFitted;
        SET DataWithFittedByLevel;
    %IF %EVAL(&UsePropOdds = 1) %THEN %DO;
        WHERE _Level_ = &SomeLevel;
    %END;
        DROP _Level_;
    RUN;
    DATA BetaFunction0;
        SET ParameterEstimates;
        WHERE ( 
               %IF %EVAL(&UsePropOdds = 0) %THEN %DO;
                  (Effect = "B1") |
               %END; 
               %DO ColumnIndex = 2 %TO %EVAL(&NumInteriorKnots + &deg);
                   (Effect = "B&ColumnIndex") |
               %END;
               %LET ColumnIndex = %EVAL(&NumInteriorKnots + &deg + 1);
               (Effect = "B&ColumnIndex") 
              );
        BetaFunction0 = Estimate;
        KEEP BetaFunction0;
    RUN; 
    %IF %EVAL(&NumTCovsIncludingIntercept>1) /* counting the intercept */
                  %THEN %DO;
        %LET stop = %EVAL(&NumInteriorKnots + &deg + 1);
        %DO ThisTCov = 1 %TO %EVAL(&NumTCovsIncludingIntercept-1);
            %LET increment = %EVAL(&NumInteriorKnots + &deg);
            DATA BetaFunction&ThisTCov;
                SET ParameterEstimates;
                %LET start = %EVAL(&stop + 1);
                %LET stop = %EVAL(&start + &increment);
                %PUT &start;
                %PUT &stop;
                %PUT &increment;
                WHERE (
                %DO ColumnIndex = &start %TO %EVAL(&stop - 1);
                       (Effect = "B&ColumnIndex") |
                %END;
                %LET ColumnIndex = %EVAL(&stop);
                    (Effect = "B&ColumnIndex")
                );
                BetaFunction&ThisTCov = Estimate;
                KEEP BetaFunction&ThisTCov;
            RUN;
        %END;
    %END; 
    PROC IML;
        USE BetaFunction0; READ ALL VAR {BetaFunction0}; CLOSE BetaFunction0;
        %IF %EVAL(&NumTCovsIncludingIntercept>1) %THEN %DO; 
            %DO ThisTCov = 1 %TO %EVAL(&NumTCovsIncludingIntercept-1);
                USE BetaFunction&ThisTCov;
                    READ ALL VAR { BetaFunction&ThisTCov};
                CLOSE BetaFunction&ThisTCov;
            %END;
        %END; 
        USE TimeBasis; READ ALL INTO TimeBasis; CLOSE TimeBasis;
        %IF %EVAL(&UsePropOdds = 0) %THEN %DO;
             BetaFunction0Fit = TimeBasis*BetaFunction0;
        %END; %ELSE %DO;
            BetaFunction0Fit = TimeBasis[,2:(%EVAL(&NumInteriorKnots + &deg + 1))]  *
                               BetaFunction0;
        %END;
        CREATE BetaFunction0Fit FROM BetaFunction0Fit[COLNAME="BetaFunction0Fit"];
            APPEND FROM BetaFunction0Fit;
        CLOSE BetaFunction0Fit;
        %IF %EVAL(&NumTCovsIncludingIntercept>1) %THEN %DO; 
            %DO ThisTCov = 1 %TO %EVAL(&NumTCovsIncludingIntercept-1);
                BetaFunction&ThisTCov.Fit = TimeBasis*BetaFunction&ThisTCov;
                CREATE BetaFunction&ThisTCov.Fit
                      FROM BetaFunction&ThisTCov.Fit[COLNAME="BetaFunction&ThisTCov.Fit"];
                    APPEND FROM BetaFunction&ThisTCov.Fit;
                CLOSE BetaFunction&ThisTCov.Fit;
            %END;
        %END; 
    QUIT;
     PROC IML;
        USE TvemWorkingDataset; READ ALL VAR {&y} INTO y; CLOSE TvemWorkingDataset;
        USE CovB; READ ALL INTO CovBWithExtraColumns; CLOSE CovB;
        %IF %EVAL(&UsePropOdds = 0) %THEN %DO;
            CovB = CovBWithExtraColumns[,2:NCOL(CovBWithExtraColumns)];
            %LET stop = &NumInteriorKnots + &deg + 1;
            %LET increment = %EVAL(&NumInteriorKnots + &deg); 
            B0Indices = 1:(&stop);
            %IF %EVAL(&NumTCovsIncludingIntercept>1) %THEN %DO;
                %DO ThisTCov = 1 %TO %EVAL(&NumTCovsIncludingIntercept-1);  
                     %LET start = %EVAL(&stop + 1);
                     %LET stop = %EVAL(&start + &increment);
                     B&ThisTCov.Indices = &start: &stop;
                %END;
            %END; 
        %END; %ELSE %DO;
            CovB = CovBWithExtraColumns[,3:NCOL(CovBWithExtraColumns)];
            uniquey = UNIQUE(y);
            uniquey = uniquey[LOC(uniquey^=.)];
            numlevels = NROW(uniquey);
            a = (numlevels-1) + 2*&NumInteriorKnots + 2*&deg + 1;
            B0Indices = numlevels:(numlevels+&NumInteriorKnots+&deg-1);
            %IF %EVAL(&debug>0) %THEN %DO;
                PRINT uniquey;
                PRINT numlevels;
                PRINT(CovB);
                PRINT(NROW(CovB));
                PRINT(NCOL(CovB));
                PRINT(B0Indices);
            %END;
            %IF %EVAL(&NumTCovsIncludingIntercept>1) %THEN %DO; 
            %LET stop = &NumInteriorKnots + &deg - 1;
            %LET increment = %EVAL(&NumInteriorKnots + &deg); 
                %DO ThisTCov = 1 %TO %EVAL(&NumTCovsIncludingIntercept-1);
                    %LET start = %EVAL(&stop + 1);
                    %LET stop = %EVAL(&start + &increment);
                    B&ThisTCov.Indices = (numlevels+&start):(numlevels+&stop);
                    %IF %EVAL(&debug>0) %THEN %DO;
                        PRINT(B&ThisTCov.Indices);
                    %END;
                %END;
            %END;
        %END;
        CovBetaFunction0 = CovB[B0Indices,B0Indices];
        CREATE CovBetaFunction0 FROM CovBetaFunction0;
            APPEND FROM CovBetaFunction0;
        CLOSE CovBetaFunction0;
        %IF %EVAL(&NumTCovsIncludingIntercept>1) %THEN %DO; 
            %DO ThisTCov = 1 %TO %EVAL(&NumTCovsIncludingIntercept-1);  
                CovBetaFunction&ThisTCov = CovB[B&ThisTCov.Indices,B&ThisTCov.Indices];
                CREATE CovBetaFunction&ThisTCov FROM CovBetaFunction&ThisTCov; 
                    APPEND FROM CovBetaFunction&ThisTCov;
                CLOSE CovBetaFunction&ThisTCov;
            %END;
        %END;
    QUIT;
    PROC IML;
        USE BetaFunction0Fit; READ ALL VAR {BetaFunction0Fit}; CLOSE BetaFunction0Fit;
        USE CovBetaFunction0; READ ALL INTO CovBetaFunction0; CLOSE CovBetaFunction0; 
        USE TimeBasis; READ ALL INTO TimeBasis; CLOSE TimeBasis; 
        SEBetaFunction0Fit = J(NROW(BetaFunction0Fit),1,0);
        DO i = 1 TO NROW(BetaFunction0Fit);
            %IF %EVAL(&UsePropOdds = 0) %THEN %DO;
                a = TimeBasis[i,];
            %END; %ELSE %DO;
                a = TimeBasis[i,2:(&NumInteriorKnots + &deg + 1)];
            %END;
            m = CovBetaFunction0;
            SEBetaFunction0Fit[i] = SQRT(a*m*a`);
        END;
        CREATE SEBetaFunction0Fit FROM SEBetaFunction0Fit[COLNAME="SEBetaFunction0Fit"]; 
            APPEND FROM SEBetaFunction0Fit;
        CLOSE SEBetaFunction0Fit;
    QUIT;
    %IF %EVAL(&NumTCovsIncludingIntercept>1) %THEN %DO; 
        %DO ThisTCov = 1 %TO %EVAL(&NumTCovsIncludingIntercept-1);  
            PROC IML; 
                USE BetaFunction&ThisTCov.Fit; READ ALL VAR {BetaFunction&ThisTCov.Fit};
                CLOSE BetaFunction&ThisTCov.Fit;
                USE CovBetaFunction&ThisTCov; READ ALL INTO CovBetaFunction&ThisTCov;
                CLOSE CovBetaFunction&ThisTCov;
                USE TimeBasis; READ ALL INTO TimeBasis; CLOSE TimeBasis; 
                SEBetaFunction&ThisTCov.Fit = J(NROW(BetaFunction&ThisTCov.Fit),1,0);
                DO i = 1 TO NROW(BetaFunction&ThisTCov.Fit);
                    a = TimeBasis[i,];
                    m = CovBetaFunction&ThisTCov;
                    SEBetaFunction&ThisTCov.Fit[i] = SQRT(a*m*a`);
                END;
                CREATE SEBetaFunction&ThisTCov.Fit
                       FROM SEBetaFunction&ThisTCov.Fit[COLNAME="SEBetaFunction&ThisTCov.Fit"];
                    APPEND FROM SEBetaFunction&ThisTCov.Fit;
                CLOSE SEBetaFunction&ThisTCov.Fit;
            QUIT; 
        %END;
    %END;
    %IF %EVAL(&NumInvariantEffects>0) %THEN %DO;
        PROC IML;
            USE ParameterEstimates;
                 READ ALL VAR {Estimate} INTO B;
            CLOSE ParameterEstimates;
    		USE CovB;
                 READ ALL INTO CovB;
            CLOSE CovB;
            EstimateInvariantBetas = B[((NROW(B)-&NumInvariantEffects+1):(NROW(B))),1];
            CovMatInvariantBetas = CovB[((NROW(CovB)-&NumInvariantEffects+1):(NROW(CovB))),
                                        ((NCOL(CovB)-&NumInvariantEffects+1):(NCOL(CovB)))];
            Estimate = J(&NumInvariantEffects,1,0);
            StdErr = J(&NumInvariantEffects,1,0);
            z = J(&NumInvariantEffects,1,0);
            p = J(&NumInvariantEffects,1,0);
            DO this = 1 TO &NumInvariantEffects;
    			Estimate[this] = EstimateInvariantBetas[this];
    			StdErr[this] = SQRT(CovMatInvariantBetas[this,this]);
    			z[this] = Estimate[this]/(StdErr[this]+1e-10);
    			p[this] = 2*(1-CDF('Normal',ABS(z[this])));
            END;
            VarName = {&cov}`;
            CREATE InvariantEffects VAR{VarName Estimate StdErr z p};
                APPEND;
            CLOSE InvariantEffects;
    	QUIT;
        PROC PRINT DATA=InvariantEffects;
            TITLE "Invariant Effects Coefficients";
        RUN;TITLE;
    %END;
    DATA DataWithFitted;
        MERGE DataWithFitted
              BetaFunction0Fit 
              SEBetaFunction0Fit 
	    %IF %EVAL(&NumTCovsIncludingIntercept>1) %THEN %DO; 
	        %DO ThisTCov = 1 %TO %EVAL(&NumTCovsIncludingIntercept-1);  
	              BetaFunction&ThisTCov.Fit
	              SEBetaFunction&ThisTCov.Fit
            %END;
            %END; 
              ;
        BetaFunction0Lower = BetaFunction0Fit-1.96*SEBetaFunction0Fit;
        BetaFunction0Upper = BetaFunction0Fit+1.96*SEBetaFunction0Fit;
        %IF %EVAL(&NumTCovsIncludingIntercept>1) %THEN %DO; 
            %DO ThisTCov = 1 %TO %EVAL(&NumTCovsIncludingIntercept-1);  
                BetaFunction&ThisTCov.Lower = BetaFunction&ThisTCov.Fit -
                                               1.96*SEBetaFunction&ThisTCov.Fit;
                BetaFunction&ThisTCov.Upper = BetaFunction&ThisTCov.Fit +
                                               1.96*SEBetaFunction&ThisTCov.Fit;
            %END;
        %END; 
    RUN;
    %IF %EVAL(&DoPlots>0) %THEN %DO;
        %IF %EVAL(&NumTCovsIncludingIntercept>1) %THEN %DO;
             TITLE &ModelType TVEM, y=&y on x=&tcov;
        %END; %ELSE %DO;
             TITLE &ModelType TVEM, y=&y, intercept only;
        %END;
        DATA DataToPlot;
            SET DataWithFitted;
            WHERE ((&t > &MinTToPlot)&(&t < &MaxTToPlot));
        RUN;
        TITLE2 "Beta Function for Intercept";
            PROC GPLOT DATA=DataToPlot;
                SYMBOL1 VALUE=NONE INTERPOL=HILOJ;
                PLOT BetaFunction0Fit*&t
                     BetaFunction0Lower*&t
                     BetaFunction0Upper*&t / VREF=0 OVERLAY;
            RUN;QUIT;
        %IF %EVAL(&NumTCovsIncludingIntercept>1) %THEN %DO; 
            %DO ThisTCov = 1 %TO %EVAL(&NumTCovsIncludingIntercept-1);  
            TITLE2 "Beta Function for Covariate &ThisTCov";
            PROC GPLOT DATA=DataToPlot;
                SYMBOL1 VALUE=NONE INTERPOL=HILOJ;
                PLOT BetaFunction&ThisTCov.Fit*&t
                     BetaFunction&ThisTCov.Lower*&t
                     BetaFunction&ThisTCov.Upper*&t / VREF=0 OVERLAY;
            RUN;QUIT; 
            %END;
        %END;
        TITLE;
        TITLE2;
    %END;
%MEND;
