   /*------------------------------------------------------------------*
   | WeightedTvem.sas version 2.6
   | For TVEM with survey sampling data.
   |
   | Created by John Dziak based on the TVEM 3.1.0 macro 
   | by Li, Dziak, Tan, Yang and Huang (2015), and using
   | the SURVEYREG and SURVEYLOGISTIC procedures in SAS.
   |
   | Copyright 2016 The Pennsylvania State University.
   |
   | Provided by The Methodology Center, The Pennsylvania State University.
   |
   | This program is free software; you can redistribute it and/or
   | modify it under the terms of the GNU General Public License as
   | published by the Free Software Foundation; either version 2 of
   | the License, or (at your option) any later version.
   |
   | This program is distributed in the hope that it will be useful,
   | but WITHOUT ANY WARRANTY; without even the implied warranty of
   | MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
   | General Public License for more details.
   |
   | Example of use : Suppose the dataset data1 has the following variables:  
   |        t:    subjects' measurement time; 
   |              note that this macro version assumes cross-sectional data;
   |        y:    the dependent variable;
   |        x0:   the intercept variable which equals 1 for all observations;  
   |        x1:   a binary covariate;
   |        x2:   a continuous covariate;
   |        w:    a sampling weight;
   |        c:    a clustering variable.
   |   Suppose we want to fit a time-varying coefficient model in which y is the 
   |   dependent variable (binary response), and x0, x1 and x2 are the covariates, 
   |   and the coefficient of x0, x1 and x2 may vary over time.  
   |   We could call the macro as follows:
   |      %WeightedTVEM(      
   |        data = data1,   
   |        time = t,         
   |        dv = y,            
   |        tvary_effect = x0 x1 x2,     
   |        knots = 5 5 5,
   |        dist = binary,
   |        weight = w, 
   |        cluster = c
   |        );
   |    A domain variable can also be used in order to specify that
   |    a TVEM model is to be fit to only part of the dataset.  Suppose 
   |    that there is a variable gender in the dataset and it is
   |    of interest to fit the model on the gender=1 subset.  
   |    In the presence of sampling and survey weights, it is more correct 
   |    to specify 
   |          domain = gender,
   |          which = 1
   |    than to subset the data in a DATA step and then run WeightedTVEM 
   |    on only the gender=1 subset.
   |    If more than one level of clustering is available (e.g., 
   |    observations within students within schools), specify only
   |    the highest-order cluster (e.g., schools).  Robust standard errors
   |    will be applied.
   |    Please leave the normalize_weights option at its default of yes,
   |    unless you are sure that you know what you are doing.  Otherwise
   |    the calculated fit statistics may be very misleading.
   |    For information about penalized pseudolikelihood fit statistics,
   |    see:
   |      Xu, C., Chen, J., & Mantel, H. (2013). Pseudo-likelihood-based 
   |       Bayesian information criterion for variable selection in survey 
   |       data. Survey Methodology, 39, 303-321.
   |      Lumley, T., & Scott, A. (2013). AIC and BIC for survey data. 
   |       Accessed online at 
   |       www.stat.colostate.edu/graybillconference2013/Presentations/Scott.pdf.
   *------------------------------------------------------------------*/
   

/**************************************************************************/
/*                         The Main Macro                                 */
/**************************************************************************/
%MACRO WeightedTVEM(data , 
            time ,
            dv ,
            tvary_effect , 
            knots , 
            dist ,
            weight ,
            cluster ,
            domain ,
            which ,
            degree = 3 ,
            evenly = 0 ,
            normalize_weights = yes,
            invar_effect = ,
			show_all = no,
            output_prefix = tvem_ ,
            outfilename = ,
            plot = full ,
            plot_scale = 100 
); 
    /* Pre-process the input data and look for obvious problems. */
    PROC IML;
        ready_to_go = 1;
        %IF %NRBQUOTE(&data) =  %THEN %DO;
               PRINT("Error in macro:" //
                     "Please specify the dataset name in mydata=.");
            ready_to_go = 0;
        %END;
        %IF %SYSFUNC(EXIST(&data)) %THEN %DO; %END;
        %ELSE %DO; 
            PRINT("Error in macro:" //
                  "The dataset &data was not found.");
            ready_to_go = 0;
        %END; 
        %IF %NRBQUOTE(&time) =  %THEN %DO;
            PRINT("Error in macro:" // 
                     "Please specify the time variable in time=.");
            ready_to_go = 0;
        %END;
        %IF %NRBQUOTE(&dv) =  %THEN %DO;
            PRINT("Error in macro:" // 
                  "Please specify a dependent variable in dep=.");
            ready_to_go = 0;
        %END;
        %IF %NRBQUOTE(&tvary_effect) =  %THEN %DO;
            PRINT("Error in macro:" // 
                  "Please specify at least one time-varying covariate in tvary_effect =.");
            ready_to_go = 0;
        %END;
        %IF %NRBQUOTE(&knots) =  %THEN %DO; 
            PRINT("Error in macro:" //
                   "Please specify an integer for each time-varying covariate in knots=.");
            ready_to_go = 0;
        %END; 
        %IF %NRBQUOTE(&domain) =  %THEN %DO; 
            %IF %NRBQUOTE(&which) =  %THEN %DO; %END; %ELSE %DO;
                PRINT("Error in macro:" //
                   "Please do not specify which level without specifying a domain variable.");
            ready_to_go = 0;
            %END;
        %END; 
        %IF %NRBQUOTE(&domain) =  %THEN %DO; %END; %ELSE %DO;
            %IF %NRBQUOTE(&which) =  %THEN %DO; 
                PRINT("Error in macro:" //
                   "Please do not specify a domain without specifying which level.");
            ready_to_go = 0;
            %END;
        %END; 
        %IF %NRBQUOTE(&which) =  %THEN %DO; %END; %ELSE %DO;
            %LET output_prefix = &output_prefix&which._;
        %END;
        NumberOfTVEArguments = COUNT(COMPBL("&tvary_effect"), " ")+1;
        NumberOfKnotArguments = COUNT(COMPBL("&knots"), " ")+1;
        IF (NumberOfTVEArguments ^= NumberOfKnotArguments) THEN DO;
            PRINT("Error in macro:  The number of arguments for TVE= must equal the" //
                  "number of arguments for knots=.  That is, there must be one value" //
                  "for knots corresponding to each time-varying-effects covariate.");
        END; 
        /* Degree option */
        internal_deg_option = &degree;
        IF ((internal_deg_option^=1) & (internal_deg_option^=2) & (internal_deg_option^=3)) THEN DO;
            PRINT("Error in macro:  The deg option should be 1, 2, or 3, to specify a" //
                  "linear, quadratic or cubic coefficient function between knots.");
            ready_to_go = 0;
        END;
        /* Distribution option */
        internal_dist_option = -1;
        IF (( UPCASE(STRIP("&dist")) = "NORMAL" ) |
            ( UPCASE(STRIP("&dist")) = "GAUSSIAN" )) THEN DO;
            internal_dist_option = 1;
        END; ELSE 
        IF (( UPCASE(STRIP("&dist")) = "LOGISTIC" ) | 
            ( UPCASE(STRIP("&dist")) = "BINARY" ) | 
            ( UPCASE(STRIP("&dist")) = "BINOMIAL" ))
            THEN DO;
            internal_dist_option = 2;
        END; 
        IF (internal_dist_option = -1) THEN DO;
            PRINT "Error in macro:  The distribution option was not recognized.";
            ready_to_go = 0;
        END;
        /* Even distribution of knots option */ 
        internal_evenly_option = -1;
        IF (( UPCASE(STRIP("&evenly")) = "0" ) |
            ( UPCASE(STRIP("&evenly")) = "no" )) THEN DO;
            internal_evenly_option = 0;
        END;
        IF (( UPCASE(STRIP("&evenly")) = "1" ) | 
            ( UPCASE(STRIP("&evenly")) = "yes" ))
            THEN DO;
            internal_evenly_option = 1;
        END;
        IF (internal_evenly_option = -1) THEN DO;
            PRINT "Error in macro:  The option for evenly= was not recognized.";
            ready_to_go = 0;
        END;        
        /* Plot option */
        internal_plot = -1;
        IF (( UPCASE(STRIP("&plot")) = "NO" ) |
            ( UPCASE(STRIP("&plot")) = "NONE" )) THEN DO;
            internal_plot = 1;
        END; 
        IF ( UPCASE(STRIP("&plot")) = "SIMPLE" ) THEN DO;
            internal_plot = 2;
        END; 
        IF (( UPCASE(STRIP("&plot")) = "FULL" ) |
            ( UPCASE(STRIP("&plot")) = "ELEGANT") ) THEN DO;
            internal_plot = 3;
        END; 
        IF (internal_plot = -1) THEN DO;
            PRINT "Error in macro:  The option for plot was not recognized.";
            ready_to_go = 0;
        END;
        /* Weights Normalizing Option */
        internal_normalize_option = -1;
        IF (( UPCASE(STRIP("&normalize_weights")) = "NO" ) |
            ( UPCASE(STRIP("&normalize_weights")) = "FALSE" )|
            ( UPCASE(STRIP("&normalize_weights")) = "0" )) THEN DO;
            internal_normalize_option = 0;
        END; ELSE 
        IF (( UPCASE(STRIP("&normalize_weights")) = "YES" ) | 
            ( UPCASE(STRIP("&normalize_weights")) = "TRUE" ) | 
            ( UPCASE(STRIP("&normalize_weights")) = "1" ))
            THEN DO;
            internal_normalize_option = 1;
        END; 
        IF (internal_normalize_option = -1) THEN DO;
            PRINT "Error in macro:  The option for normalize_weights was not recognized.";
            ready_to_go = 0;
        END;
        /* Show All Output Option -- Default is No*/
        internal_show_all_option = 0; 
        IF (( UPCASE(STRIP("&show_all")) = "YES" ) | 
            ( UPCASE(STRIP("&show_all")) = "TRUE" ) | 
            ( UPCASE(STRIP("&show_all")) = "1" ))
            THEN DO;
            internal_show_all_option = 1;
        END;  
        /* Plot Scale */
        plot_scale = &plot_scale;
        IF ((plot_scale < 100) | (plot_scale > 10000)) THEN DO;
            PRINT("Warning in macro:" // 
                  "We recommend a plot_scale value between 100 and 10000 for the plot."); 
        END;
        /* Check if the output dataset name prefix is too long. */
        IF (LENGTH("&output_prefix")>15) THEN DO;
            PRINT("Error in macro:" // 
                  "Please specify something shorter for the output prefix.");
            ready_to_go = 0;
        END;  
        /* Check if the dependent variable has the correct distribution. */
        USE &data;
            READ ALL VAR {&dv} INTO dep;
                IF SUM(dep^=.)=0 THEN DO;
                    PRINT("Error in macro:" // 
                          "Could not find data for dependent variable.");                
                    ready_to_go = 0;
                END;
                IF internal_dist_option = 2 THEN DO;
                    num_okay = SUM(dep = .)+SUM(dep = 0)+SUM(dep = 1);
                    num_total = NROW(dep);
                    IF (num_okay < num_total) THEN DO;                        
                        PRINT("Error in macro:" // 
                              "The &dv variable should contain only 0's and 1's " //
                              "for a binary logistic regression TVEM.");
                        ready_to_go = 0;
                    END;    
                END;
                nonmissingdep = dep[LOC(dep^=.)];
        CLOSE &data;
        /* Get ready to send information to the back-end macro */
        CALL SYMPUT("internal_deg_option", CHAR(internal_deg_option));
        CALL SYMPUT("internal_dist_option", CHAR(internal_dist_option));
        CALL SYMPUT("internal_evenly_option", CHAR(internal_evenly_option));
        CALL SYMPUT("internal_normalize_option", CHAR(internal_normalize_option));
        CALL SYMPUT("internal_show_all_option", CHAR(internal_show_all_option));
        CALL SYMPUT("internal_plot", CHAR(internal_plot));
        CALL SYMPUT("ready_to_go", CHAR(ready_to_go));
    QUIT;
    /* Unless there is an obvious error, call the appropriate back-end macro. */
    %IF %EVAL(&ready_to_go = 1) %THEN %DO;
            %_TvemW( data  = &data,  
                     time = &time,         
                     dv = &dv,            
                     invar_effect = &invar_effect,            
                     tvary_effect = &tvary_effect,           
                     knots = &knots,   
                     dist_option = &internal_dist_option,
                     weight = &weight,
                     cluster = &cluster,
                     domain = &domain,
                     which = &which,
                     degree = &internal_deg_option,              
                     evenly = &internal_evenly_option,  
                     normalize = &internal_normalize_option,  
                     plot_scale = &plot_scale,
                     plot = &internal_plot, 
                     output_prefix = &output_prefix, 
                     outfilename = &outfilename  
               );
    %END;
    %IF %EVAL(&internal_show_all_option = 1) %THEN %DO; %END; %ELSE %DO;
		PROC DATASETS NOLIST NOWARN;
			DELETE  _InvariantEffectscovariates
					_InvariantEffectscovariatesses
					_t_dense
					_data_grid_distinct
					_data_grid_orig
					_data_grid_orig_spline
					_data_grid_sorted
					_data_grid_spline
					_data_grid_TVE
					_dense_grid_coef
					_dense_grid_covb
					_dense_grid_spline
					_dt_spline 
					_id_dep_cov
					_modelcovariance 
					_plot_data
					_plot_data_or
					_RobustCovariance
					_TVE_spline
					_tmp
					_tmp_mydata
					_tvemp_id
					_tvemp_mu
					_tvemp_temp_id
					_tvemp_y
					_t_dense
					_t_distinct
					_t_orig
					_t_orig_spline
					_t_sorted
					_t_spline
					_t_TVE;
		QUIT;
    %END;
%MEND;

/**************************************************************************/
/*                                 _TvemW                                 */
/**************************************************************************/
%MACRO _TvemW(data ,  
              time ,         
              dv ,            
              invar_effect ,            
              tvary_effect ,           
              knots ,   
              dist_option , 
              weight ,
              cluster ,
              domain ,
              which ,
              degree ,              
              evenly ,  
              normalize , 
              plot_scale ,
              plot , 
   /* 1 = none, 2 = simple, 3 = polished */
              output_prefix ,
              outfilename  
   );
    /*Step 1: Process the input data and create some datasets for 
            intermediate calculations. */
    DATA _data_grid_TVE;   /* covariates with time-varying coefficients */
        SET &data;
        WHERE &time >. ;
        KEEP &tvary_effect; 
    RUN;
    DATA _id_dep_cov; /* Some variables needed for later processing */
        SET &data;
        WHERE &time >. ;
        KEEP &dv &invar_effect &time &weight &cluster &domain; 
    RUN;
    DATA _data_grid_orig;   /* covariates with time-varying coefficients */
        SET &data;
        WHERE &time >. ;
        KEEP &time; 
    RUN;
    PROC SORT DATA=_data_grid_orig OUT=_data_grid_sorted;  
            /*sorted time variable, needed for inner knots calculation */ 
        BY &time;
    RUN;
    PROC SORT DATA=_data_grid_orig OUT=_data_grid_distinct NODUPKEY;  
            /*sorted time variable, needed for inner knots calculation */ 
        BY &time;
    RUN;
    DATA _data_grid_orig;  
        /* Keep _idx to order the observations, 
           be consistent with _id_dep_cov  */
        RETAIN _idx 0;
        SET _data_grid_orig;
        _idx = _idx + 1;
    RUN;
    /* Step 2: Create the design matrix columns to represent the
     B-spline transformation of &tvary_effect */
    PROC IML;
        /* Part 2.1: Define three needed functions: my_Bspline(...),
                     my_vec_Bspline(...), and my_knots(...) */
        START my_Bspline(x, knots, d);
        /*    1.  Evaluate the values of all (= num of knots + d + 1) of the 
              B-spline basis functions at x, where x is a real number. */
        /*    The knots include both exterior and inner knots. */
        /*    2.  Output a vector of (num of knots) + 1 - d real numbers */
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
                /* the De Boor formula, refer to Eilers & Marx (1996) */
                DO i = 1 TO (n_ext-j-1); 
                    w1 = (x-all_knots[i])/(all_knots[i+j] - all_knots[i]); 
                    w2 = (all_knots[i+j+1]-x)/(all_knots[i+j+1] - all_knots[i+1]); 
                    tmp[i] = w1*tmp[i] +w2*tmp[i+1];
                END;
                j=j+1;
            END;    
            RETURN(tmp[1:(n_ext-d-1)]);
        FINISH;
        START my_vec_Bspline(xx, knots, d);
        /*    evaluate the values of the (num of knots) + d +1 B-spline */
        /*        basis functions at xx */ 
        /*    xx is a vector of real numbers */
        /*    output a real matrix with dim(xx) rows,  */
        /* and (num of knots) + d +1 columns */
            n_row = ncol(xx);
            n_col = ncol(knots)+2-d-1;        
            out = J(n_row, n_col, .);  
            DO i=1 TO n_row; 
                out[i, ] = t(my_Bspline(xx[i], knots, d));
            END;            
            RETURN(out);
        FINISH;
        START my_knots(samp, nknots, evenly, d);  
        /* Calculate nknots inner knots and add 2*d exterior knots */
        /* the inner knots are uniform in space (evenly=1) or uniformly */
        /*    on quantiles (evenly != 1) */
        /* samp should be sorted ascendingly, i.e., from small to large */
            len = ncol(samp);
            *len = max(nrow(samp), );
            all_knots = 0*(1:(nknots+2*d));
            all_knots[1:d] = samp[1]-(1e-12)*(d:1);
            all_knots[(nknots+d+1):(nknots+2*d)] = samp[len]+(1e-12)*(1:d);
            IF (evenly=1) THEN DO;
                dist  = (samp[len]-samp[1])/(nknots+1);         
                all_knots[(d+1):(d+nknots)] = samp[1]+dist*(1:nknots);
            END; ELSE DO; 
                 j=1; 
                DO i = 1 TO nknots;
                    DO WHILE ( (j-1)/(len-1) < i/(nknots+1) );
                        j=j+1;
                    END;
                    all_knots[d+i] = samp[j];
                END;
            END;
            RETURN(all_knots);
        FINISH;
        /* Part 2.2: generate Bspline matrix, given t_cov and options 
        /*         for knots, and read data into matrix */
        USE _data_grid_sorted;
        READ ALL INTO mat_data_grid_sorted; /* &time only, sorted */
        USE _data_grid_distinct;
        READ ALL INTO mat_data_grid_distinct; 
                        /* &time only, distinct time points only */
        n_pt = &plot_scale -1;
        mat_t_dense = T( mat_data_grid_sorted[1] +
                    ((0:n_pt)/n_pt)*
                    ( mat_data_grid_sorted[ nrow(mat_data_grid_sorted) ] - 
                        mat_data_grid_sorted[1] ) );
        k_s = {&knots};
            /* a vector of number of inner knots, one for each 
                covariate in &tvary_effect */
        DO i = 1 TO ncol(k_s);
            all_knots = my_knots(t(mat_data_grid_sorted), k_s[i], 
                    &evenly, &degree);
            mat_data_grid_splines = mat_data_grid_splines || 
                    ( my_vec_Bspline( t(mat_data_grid_distinct), 
                                      all_knots, &degree) );
                     /* distinct time points */
            mat_dense_grid_splines = mat_dense_grid_splines || 
                     ( my_vec_Bspline( t(mat_t_dense),
                                       all_knots, &degree) ); 
                     /* uniform distributed time points */
        END;
        mat_data_grid_splines = mat_data_grid_distinct ||
                                mat_data_grid_splines;
        CREATE _data_grid_spline FROM mat_data_grid_splines[colname= "_time"]; 
        APPEND FROM mat_data_grid_splines; 
        /* output mat_dense_grid_splines to dataset dt_spline
            dataset t_dense is useful for plotting  */
        CREATE _dense_grid_spline FROM mat_dense_grid_splines; 
        APPEND FROM mat_dense_grid_splines;
        CREATE _t_dense FROM mat_t_dense [colname= "&time"]; 
                                                /* denser time points */
        APPEND FROM mat_t_dense;
    QUIT;
    PROC SQL;
        CREATE TABLE _data_grid_orig_spline AS
        SELECT A.*, B.*  
        FROM _data_grid_orig A, _data_grid_spline B 
        WHERE A.&time = B._time
        ORDER BY _idx;
    QUIT;
    PROC IML;
        k_s = {&knots}; 
            /* a vector of number of inner knots, one for each */
            /*    covariate in &tvary_effect */
        j = k_s[+] + (&degree+1)*ncol(k_s);
        c_name = ( compress(  CATT( "col", char(j+1, 8, 0)) ) ); 
        USE _data_grid_orig_spline;
        READ ALL VAR("col2":c_name) INTO mat_tspline;  
                /* &time &tvary_effect, unsorted,*/ 
        USE _data_grid_TVE;
        READ ALL VAR{&tvary_effect} INTO mat_TVE;  /* &time &tvary_effect, unsorted,*/ 
              /* We use the syntax VAR{&tvary_effect} to make the order of columns */
              /* in mat_TVE consistent with &tvary_effect (and hence with */
              /* knots) */
        id_u = 0; 
        DO i = 1 TO ncol(k_s);
            id_l = id_u + 1; id_u = id_l + (k_s[i]+&degree+1)-1;
            mat_data_grid_splines = mat_data_grid_splines ||
                            ( (mat_tspline[, id_l: id_u] ) # mat_TVE[,i]);
        END;
        /* Output the matrix mat_TVE_splines to the dataset TVE_spline, and /*
        /* output the matrix mat_dense_grid_splines to the dataset               /*
        /*   dt_spline. */
        cname = CATT("_B",1:NCOL(mat_data_grid_splines));
        CREATE _TVE_spline FROM mat_data_grid_splines[colname=cname];
            APPEND FROM mat_data_grid_splines; 
        CLOSE _TVE_spline;
    QUIT;
    /* Step 3a: Get the internal analysis dataset ready. */ 
    DATA model_data;  /* merge intermediate datasets */
        MERGE _id_dep_cov _TVE_spline;
    RUN; 
    /* In the internal copy of the dataset, remove format labels that 
       may interfere with macro processing */ 
    PROC DATASETS NOLIST NOWARN;
       MODIFY model_data; 
         ATTRIB _all_ label=' '; 
         ATTRIB _all_ format=;
    RUN;
    /* Standardize weights : */
       %IF %NRBQUOTE(&weight) =  %THEN %DO;
            DATA model_data;  SET model_data; 
                SampleWeightForTvem = 1;
            RUN;
        %END; %ELSE %DO;           
            %IF %EVAL(&normalize = 1) %THEN %DO;
                PROC IML;
                    USE model_data; 
                          READ ALL VAR {&weight} INTO weight WHERE (&weight ^= .);
                    CLOSE model_data;
                    mean_weight = weight[:];
                    CALL SYMPUT("MeanWeightForTvem",CHAR(mean_weight));
                QUIT;
                DATA model_data;  SET model_data; 
                       SampleWeightForTvem = &weight / &MeanWeightForTvem;
                RUN;
            %END; %ELSE %DO;
                DATA model_data;  SET model_data; 
                    SampleWeightForTvem = &weight;
                RUN;
            %END;
        %END;
    RUN;
    PROC CONTENTS DATA=_TVE_spline OUT=_tmp NOPRINT;    RUN;
    PROC SQL; /* find out the number of variables reduced from B-splines */
        RESET NOPRINT;
        SELECT (MAX(VARNUM)) INTO :aa FROM _tmp;  
    QUIT;
    %LET size = %TRIM(&aa);	 
    %IF %EVAL(&internal_show_all_option = 1) %THEN %DO;  
		ODS EXCLUDE NONE; RUN;
		ODS RESULTS; RUN;
	%END; %ELSE %DO;
		ODS EXCLUDE ALL; RUN; /* This is like "NOPRINT" in older PROC's */
		ODS NORESULTS; RUN; 
	%END;
    /* Step 3b: Fit the underlying parametric model. */ 
    %IF %EVAL(&dist_option = 1) %THEN %DO; 
        PROC SURVEYREG DATA=model_data ;
            WEIGHT SampleWeightForTvem;
            %IF %NRBQUOTE(&cluster) =  %THEN %DO; %END; %ELSE %DO;           
                CLUSTER &cluster;
            %END;
            %IF %NRBQUOTE(&domain) =  %THEN %DO; %END; %ELSE %DO;           
                DOMAIN &domain;
            %END;
            MODEL &dv = _B1-_B&size &invar_effect / COVB NOINT ; 
            OUTPUT OUT=_survey_out  
                   PRED=predicted 
                   STD=stderr_predicted 
                   RESIDUAL=tvem_residual;                  
            ODS OUTPUT CovB=_survey_CovB
                       FitStatistics=_survey_stats  
                       ParameterEstimates=_survey_est 
                       DataSummary = _survey_nobs
                       %IF %NRBQUOTE(&domain) =  %THEN %DO; %END; %ELSE %DO;  
                       DomainSummary = _survey_domain
                       %END;
                       ;      
        RUN;    
        %IF %NRBQUOTE(&domain) =  %THEN %DO; %END; %ELSE %DO;  
            DATA _survey_stats; SET _survey_stats; 
                WHERE (COMPRESS(LOWCASE(Domain)," ")) = COMPRESS(LOWCASE("&domain=&which")," ");
            RUN;  
            DATA _survey_est; SET _survey_est; 
                WHERE (COMPRESS(LOWCASE(Domain)," ")) = COMPRESS(LOWCASE("&domain=&which")," ");
            RUN;
            DATA _survey_out; SET _survey_out; 
                WHERE (COMPRESS(LOWCASE(Domain)," ")) = COMPRESS(LOWCASE("&domain=&which")," ");
            RUN;
			DATA _survey_all_domains; SET _survey_domain; RUN;
            DATA _survey_domain; SET _survey_domain; 
                WHERE (COMPRESS(LOWCASE(Domain)," ")) = COMPRESS(LOWCASE("&domain=&which")," ");
            RUN;
        %END;
		PROC IML; 
			%IF %NRBQUOTE(&domain) =  %THEN %DO; 
				USE _survey_out;
					READ ALL VAR {tvem_residual SampleWeightForTvem};
				CLOSE _survey_out; 
				GoodRows = (LOC((tvem_residual^=.)&(SampleWeightForTvem^=.)))`; 
			%END; %ELSE %DO;   
				USE _survey_out;
					READ ALL VAR {tvem_residual &domain SampleWeightForTvem};
				CLOSE _survey_out; 
				GoodRows = (LOC((tvem_residual^=.)&(&domain=&which)&(SampleWeightForTvem^=.)))`;
			%END;  
			SampleWeightForTvem = SampleWeightForTvem[GoodRows];
			tvem_residual = tvem_residual[GoodRows];
			NObs = NROW(GoodRows);
			WeightedN = SampleWeightForTvem[+];
			WeightedMSE = ((SampleWeightForTvem#(tvem_residual##2))[+])/WeightedN;
			Neg2Loglik = WeightedN*LOG(2*3.14159)+WeightedN*LOG(WeightedMSE)+WeightedN;
			USE _survey_est;
				READ ALL VAR {Estimate};
			CLOSE _survey_est; 
			P = NROW(Estimate);
			AIC = Neg2Loglik + 2*P;
			BIC = Neg2Loglik + LOG(NObs)*P; 
			CREATE _new_survey_fitstats VAR {NObs WeightedN WeightedMSE P Neg2Loglik AIC BIC};
				APPEND;
			CLOSE _new_survey_fitstats;
		QUIT; 
    %END;
    %IF %EVAL(&dist_option = 2) %THEN %DO; 
        PROC SURVEYLOGISTIC DATA=model_data ;
            MODEL &dv (DESCENDING) = _B1-_B&size &invar_effect / COVB NOINT ; 
            OUTPUT OUT=_survey_out  
                   PRED=predicted  ;             
            WEIGHT SampleWeightForTvem;
            %IF %NRBQUOTE(&cluster) =  %THEN %DO; %END; %ELSE %DO;           
                CLUSTER &cluster;
            %END;
            %IF %NRBQUOTE(&domain) =  %THEN %DO; %END; %ELSE %DO;           
                DOMAIN &domain;
            %END;     
            ODS OUTPUT CovB=_survey_CovB
                       FitStatistics=_survey_stats  
                       ParameterEstimates=_survey_est  
                       NObs = _survey_nobs 
                       %IF %NRBQUOTE(&domain) =  %THEN %DO; %END; %ELSE %DO;  
                       DomainSummary = _survey_domain
                       %END;   
                       ;
        RUN;     
        %IF %NRBQUOTE(&domain) =  %THEN %DO; %END; %ELSE %DO;   
            DATA _survey_stats; SET _survey_stats; 
                WHERE (COMPRESS(LOWCASE(Domain)," ")) = COMPRESS(LOWCASE("&domain=&which")," ");
            RUN;  
            DATA _survey_est; SET _survey_est; 
                WHERE (COMPRESS(LOWCASE(Domain)," ")) = COMPRESS(LOWCASE("&domain=&which")," ");
            RUN;
            DATA _survey_out; SET _survey_out; 
                WHERE (COMPRESS(LOWCASE(Domain)," ")) = COMPRESS(LOWCASE("&domain=&which")," ");
            RUN;
			DATA _survey_all_domains; SET _survey_domain; RUN;
            DATA _survey_domain; SET _survey_domain; 
                WHERE (COMPRESS(LOWCASE(Domain)," ")) = COMPRESS(LOWCASE("&domain=&which")," ");
            RUN;
        %END; 
		 PROC IML;
			%IF %NRBQUOTE(&domain) =  %THEN %DO; 
				USE _survey_out;
					READ ALL VAR {predicted &dv SampleWeightForTvem};
				CLOSE _survey_out;
				GoodRows = LOC((predicted^=.)&
							   (&dv^=.)&
							   (SampleWeightForTvem^=.))`;	
			%END; %ELSE %DO;   
				USE _survey_out;
					READ ALL VAR {predicted &dv &domain SampleWeightForTvem};
				CLOSE _survey_out;
				GoodRows = LOC((predicted^=.)&
							   (&dv^=.)&
							   (&domain=&which)&
							   (SampleWeightForTvem^=.))`;
			%END;
			y = &dv[GoodRows];
			mu = predicted[GoodRows];
			w = SampleWeightForTvem[GoodRows];
			NObs = NROW(GoodRows);
			LogLik = (w#(y#LOG(mu+1e-20)+(1-y)#LOG(1-mu+1e-20)))[+];
			Neg2Loglik = -2*LogLik; 
			USE _survey_est;
				READ ALL VAR {Estimate};
			CLOSE _survey_est; 
			P = NROW(Estimate);
			AIC = Neg2Loglik + 2*P;
			BIC = Neg2Loglik + LOG(NObs)*P; 
			CREATE _new_survey_fitstats VAR {NObs WeightedMSE P Neg2Loglik AIC BIC};
				APPEND;
			CLOSE _new_survey_fitstats;
		QUIT;  
	%END;
    DATA _InvariantEffectsCovariates;
        SET _survey_est;
        %IF %EVAL(&dist_option = 1) %THEN %DO; 
            WHERE ( SUBSTR(Parameter,1,2) ^= "_B");      
        %END;    
        %IF %EVAL(&dist_option = 2) %THEN %DO; 
            WHERE ( SUBSTR(Variable,1,2) ^= "_B");      
        %END;    
    RUN;
    ODS EXCLUDE NONE; RUN;
    ODS RESULTS; RUN;
    /* Step 4: Calculate the coefficient functions and their pointwise 
                confidence bands using Cramer's delta method 
                (Taylor linearization). */
    PROC IML;
        USE _dense_grid_spline;
            READ ALL INTO mat_dense_grid_spline;  
                    /* Spline basis functions at dense grid of time points */
        CLOSE _dense_grid_spline;    
        USE _survey_est;
            READ ALL VAR {estimate} INTO mat_fix_est;
        CLOSE _survey_est;
        k_s = {&knots}; 
            /* a vector of number of inner knots, one for each 
                covariate in &tvary_effect */
        j = k_s[+] + (&degree+1)*NCOL(k_s);
        c_name = ( compress(  CATT( "_B", char(j, 8, 0)) ) ); 
        USE _survey_CovB;
            READ ALL VAR("_B1":c_name) INTO mat_covB;  
            /* only read columns relevant to B-spline basis */
        CLOSE _survey_CovB; 
        idx_u=0;  /* index sets for picking up suitable columns */
        DO i = 1 TO ncol(k_s);
            idx_l=idx_u+1;
            idx_u = idx_l + (k_s[i]+&degree+1) -1 ; 
            /* values at denser points */
            Coef = mat_dense_grid_spline[, idx_l:idx_u]*
                    mat_fix_est[idx_l:idx_u];
            /* pointwise variances */
            cov_Coef = VECDIAG( ( mat_dense_grid_spline[, idx_l:idx_u]*
                                  mat_covB[idx_l:idx_u, (idx_l):(idx_u)])*
                                  T(mat_dense_grid_spline[, idx_l:idx_u]) );
            IF ( MIN(cov_Coef)<0 ) THEN cov_Coef[loc(cov_Coef<0)] = .;  
                        /* if <0, then missing*/
            /* calculate CIs */ 
            Coef_StdErr = sqrt(cov_Coef);
            Coef_L = Coef - 1.96*Coef_StdErr;
            Coef_U = Coef + 1.96*Coef_StdErr;
            mat_coef = mat_coef || 
                        (( (coef_L||coef) || coef_U) || Coef_StdErr);
            var_name = SCAN("&tvary_effect", i);
            var_names = ( CATT(var_name, "_L") || var_name || 
                        CATT(var_name,"_U") || CATT(var_name,"_SE") );
            cname = cname || var_names;
        END;
        CREATE _dense_grid_coef FROM mat_coef[colname=cname]; 
                /* denser time points */
            APPEND FROM mat_coef;
        CLOSE _dense_grid_coef;
        CREATE _dense_grid_covb FROM mat_covb; /* denser time points */
            APPEND FROM mat_covb;
        CLOSE _dense_grid_covb;
    QUIT;
    DATA _plot_data;
        MERGE _t_dense _dense_grid_coef;
    RUN;  
    /* Step Five:  Print the output */
    %_TvemOut(data = &data,  
             time = &time,         
             dv = &dv,            
             invar_effect = &invar_effect,        
             tvary_effect = &tvary_effect,           
             knots = &knots,   
             degree = &degree,              
             evenly = &evenly,       
             plot_scale = &plot_scale,
             dist_option = &dist_option,
             plot = &plot, 
             are_not_using_penalty = 1,
             output_prefix = &output_prefix,
             outfilename = &outfilename );
%MEND;


/**************************************************************************/
/*                                 _TvemOut                               */
/**************************************************************************/
%MACRO _TvemOut(data,    
              time,         
              dv,            
              invar_effect,            
              tvary_effect,           
              knots,   
              degree,              
              evenly,       
              plot_scale,
              dist_option,   
              plot, 
              are_not_using_penalty,
              output_prefix,
              outfilename 
   );
    /* Step 1: Create the dataset for plotting the fitted coefficient  
       functions and their pointwise confidence bands. */   
    PROC CONTENTS DATA=_dense_grid_coef OUT=_tmp NOPRINT;    RUN;
    PROC SQL; 
        /* find out the number of variables, reduced from B-splines */
        RESET NOPRINT;
        SELECT (MAX(VARNUM))/4 INTO :aa FROM _tmp;  
    QUIT; 
    %LET size = %TRIM(&aa);
    /* Step 2: Output this dataset to a file if requested. */ 
    %IF %NRBQUOTE(&outfilename) ^=  %THEN %DO; 
        PROC EXPORT DATA=_plot_data 
                    OUTFILE= "&outfilename" 
                    DBMS=CSV REPLACE;
        RUN;
    %END;
    /* Step 3: Create exponentiated version of this dataset. */                                                                                                       
    %GLOBAL list;                                                                                                                         
    %LET list=;                                                                                                                           
    /** open dataset **/                                                                                                                  
    %LET dsid=%SYSFUNC(OPEN(_plot_data));
    /** cnt will contain the number of variables in the 
        dataset passed in **/
    %LET cnt=%SYSFUNC(ATTRN(&dsid,nvars));
    %DO i = 1 %TO &cnt;
        %LET list=&list %SYSFUNC(varname(&dsid,&i));
    %END;
    /** close dataset **/                                                                                                                 
    %LET rc=%SYSFUNC(CLOSE(&dsid));                                                                                                       
    DATA _plot_data_OR(DROP=&list RENAME=(&time._old=&time));                                                                                                                           
        SET _plot_data;            
        &time._old=&time;  
        &time = 0; 
        %DO i = 1 %TO &cnt;                                                                                                                
        %LET var=%SCAN(&list,&i);                                                                                                         
        exp_&var=EXP(&var);                                                                                                               
        %END;     
        DROP exp_&time;
    RUN;           
    /* Provide output on the screen or listing */
    PROC IML;
        OPTIONS FORMCHAR="|----|+|---+=|-/\<>*";
        FILE PRINT;
        TITLE3 "Weighted Time-Varying Effects Modeling (TVEM) Macro Output"; 
        PUT "================================================================";
        PUT "Weighted Time-Varying Effects Modeling (TVEM) Macro Output";
        PUT "================================================================";
        PUT "Dataset:                         &data";
        PUT "Time variable:                   &time";
        PUT "Response variable:               &dv";
        %IF %TRIM("&weight") ^= "" %THEN %DO;
            PUT "Weighting variable:              &weight";
        %END;
        %IF %TRIM("&cluster") ^= "" %THEN %DO;
            PUT "Clustering variable:             &cluster";
        %END;
        %IF %TRIM("&domain") ^= "" %THEN %DO;
            PUT "Domain:                          &domain = &which.";
        %END;
        %IF %EVAL(&dist_option = 1) %THEN %DO; 
            PUT "Response distribution:           Normal (Gaussian)";
        %END;
        %IF %EVAL(&dist_option = 2) %THEN %DO; 
            PUT "Response distribution:           Binary (Logistic)";
        %END;
        %IF %EVAL(&dist_option = 3) %THEN %DO; 
            PUT "Response distribution:           Poisson";
        %END;
        %IF %TRIM("&invar_effect") ^= "" %THEN %DO;
            PUT "Non-time-varying effects:        &invar_effect";
        %END;
        %IF %TRIM("&tvary_effect") ^= "" %THEN %DO;
            PUT "Time-varying effects:            &tvary_effect";
            PUT "Knots for splines:               &knots";
            PUT "Degree for splines:              &degree";
        %END;
        PUT "================================================================";
        %IF %SYSFUNC(EXIST(_new_survey_fitstats)) %THEN %DO;
            USE _new_survey_fitstats;
                READ ALL VAR {NObs Neg2Loglik AIC BIC};
            CLOSE _new_survey_fitstats;
            PUT "Fit Statistics:";
            PUT "Number of observations used:              " NObs;
            PUT "Negative Two Pseudolikelihood:            " Neg2Loglik;
            PUT "Pseudolikelihood AIC:                     " AIC; 
            PUT "Pseudolikelihood BIC:                     " BIC; 
            PUT "================================================================";
        %END;
        %IF %SYSFUNC(EXIST(_survey_est)) %THEN %DO;
            %IF %SYSFUNC(EXIST(_plot_data)) %THEN %DO;
                PUT "The estimated coefficient functions are stored in the dataset &output_prefix.plot_data.";        
                %IF %NRBQUOTE(&outfilename) ^=  %THEN %DO; 
                    PUT "They are also stored in the file &outfilename..";
                %END;
            %END;
        %END; %ELSE %DO;
            %PUT "Error in macro:  Estimation failed.";
        %END;
        PUT "================================================================";
    QUIT; TITLE1;     
    %IF %SYSFUNC(EXIST(_InvariantEffectsCovariates)) %THEN %DO;
        %IF %NRBQUOTE(&invar_effect) ^=  %THEN %DO;
            PROC PRINT DATA=_InvariantEffectsCovariates; 
                TITLE3 "Fixed Effects Covariates"; 
            RUN; TITLE3;
        %END;
        DATA &output_prefix.invariant_effects; 
            SET _InvariantEffectsCovariates; 
        RUN;
    %END;   
    %IF (&are_not_using_penalty) %THEN %DO;
        %IF %SYSFUNC(EXIST(_survey_covparms)) %THEN %DO;
            PROC PRINT DATA=_survey_covparms;
                TITLE1 "Covariance Parameters";
                WHERE CovParm ^= "Residual (VC)";
            RUN; TITLE1;
            DATA &output_prefix.covparms; SET _survey_covparms; RUN;
        %END;   
    %END;   
    %IF %SYSFUNC(EXIST(_plot_data)) %THEN %DO;
        DATA &output_prefix.plot_data; SET _plot_data; RUN;
    %END;
    %IF %SYSFUNC(EXIST(_plot_data_OR)) %THEN %DO;
        %IF %EVAL(&dist_option = 2) %THEN %DO; 
            DATA &output_prefix.plot_data_OR; SET _plot_data_OR; RUN;
        %END;
    %END;
    %IF %SYSFUNC(EXIST(_survey_converged)) %THEN %DO;
        DATA &output_prefix.converged; SET _survey_converged; RUN;
    %END;
    %IF %SYSFUNC(EXIST(_survey_stats)) %THEN %DO;
        DATA &output_prefix.fitstats; SET _new_survey_fitstats; RUN;
    %END;  
    %IF %SYSFUNC(EXIST(_survey_nobs)) %THEN %DO;
        DATA &output_prefix.nobs; SET _survey_nobs; RUN;
    %END;  
    %IF %SYSFUNC(EXIST(_survey_domain)) %THEN %DO;
        DATA &output_prefix.domain; SET _survey_domain; RUN;
    %END;  
    %IF %SYSFUNC(EXIST(_survey_est)) %THEN %DO;
        DATA &output_prefix.orig_coefs; SET _survey_est; RUN; 
    %END; %ELSE %DO;
        PROC IML;  
            PRINT "Error in macro:  Estimation was not successful."; 
        QUIT;
    %END;  
    %IF %NRBQUOTE(&domain) =  %THEN %DO; %END; %ELSE %DO;
        %IF %SYSFUNC(EXIST(_dense_grid_covb)) %THEN %DO;
            DATA &output_prefix.covar_orig_coefs; SET _dense_grid_covb; RUN;
        %END;
    %END;
    %IF %SYSFUNC(EXIST(_survey_out)) %THEN %DO;
        DATA &output_prefix.predicted; SET _survey_out; RUN;
    %END; 
    /* Plot the dataset. */ 
    %DO i=1 %TO &size;
        %LET var_name = %SCAN(&tvary_effect, &i);
        /* Plot the coefficient function by time.*/     
        %IF %EVAL(&plot = 2) %THEN %DO;  
                    /* Very simple graph style; */ 
            TITLE1 "Coefficient Curve and 95% Confidence Band for Covariate: &var_name";
            AXIS  
                LABEL=('Coefficient' JUSTIFY=RIGHT "of &var_name")
                  WIDTH=3;
            LEGEND LABEL=NONE
                   VALUE=('Lower Bound' 'Estimate' 'Upper Bound')
                   POSITION=(bottom left inside) MODE=SHARE DOWN = 2; 
            SYMBOL1 COLOR=green VALUE=point INTERPOL=join LINE=2 WIDTH=3;
            SYMBOL2 COLOR=red VALUE=point INTERPOL=join LINE=1 WIDTH=3;
            SYMBOL3 COLOR=green VALUE=point INTERPOL=join LINE=2 WIDTH=3;
            PROC GPLOT DATA=_plot_data;
                PLOT (&var_name._L &var_name &var_name._U )* &time / OVERLAY
                     LEGEND = legend VAXIS=axis;
                RUN;
            QUIT;
            TITLE1; SYMBOL1; SYMBOL2; SYMBOL3; AXIS; LEGEND;
        %END; 
        %IF %EVAL(&plot = 3) %THEN %DO;
            PROC TEMPLATE;  /* More elegant graph style, requiring Java; */                                                                                                     
                DEFINE STATGRAPH mygraph_tvem;                                                                                                             
                    BEGINGRAPH;                                                                                                                           
                        ENTRYTITLE 
        "Coefficient Curve and 95% Confidence Band for Covariate: &var_name";                                                                                                             
                        LAYOUT LATTICE / COLUMNDATARANGE=union;
                            /* Overlay a Bandplot and a Seriesplot */                                                                                               
                            LAYOUT OVERLAY /
                                YAXISOPTS =(label="Coefficient of &var_name");                                                                               
                                BANDPLOT X=&time 
                                         LIMITUPPER=&var_name._U 
                                         LIMITLOWER=&var_name._L /
                                             DATATRANSPARENCY=.3 
                                             MODELNAME="fit" 
                                         LEGENDLABEL="95% Confidence Limits";                                                         
                                SERIESPLOT x=&time y=&var_name;  
                                REFERENCELINE y=0/LINEATTRS=(pattern=4);   
                            ENDLAYOUT;                               
                        ENDLAYOUT;                                                                                                                         
                    ENDGRAPH;                                                                                                                             
                END; 
            RUN;   
            PROC SGRENDER DATA=_plot_data TEMPLATE=mygraph_tvem;
            RUN;
        %END; 
        %IF %EVAL(&dist_option = 2) %THEN %DO;
            /* Plot the odds ratio (the exponent of the coefficient
                function) by time.*/ 
            %IF %EVAL(&plot = 2) %THEN %DO;       
                            /* Very simple graph style; */ 
                TITLE1  
        "Odds Ratio and 95% Confidence Band for Covariate: &var_name";
                AXIS  
                    LABEL=('Odds Ratio' JUSTIFY=RIGHT "of &var_name")
                      WIDTH=3;
                LEGEND LABEL=NONE
                       VALUE=('Lower Bound' 'Estimate' 'Upper Bound')
                       POSITION=(bottom left inside) MODE=SHARE DOWN = 2; 
                SYMBOL1 COLOR=green VALUE=point INTERPOL=join LINE=2 WIDTH=3;
                SYMBOL2 COLOR=red VALUE=point INTERPOL=join LINE=1 WIDTH=3;
                SYMBOL3 COLOR=green VALUE=point INTERPOL=join LINE=2 WIDTH=3;
                PROC GPLOT DATA=_plot_data_OR;
                    PLOT (exp_&var_name._L exp_&var_name exp_&var_name._U )*
                            &time / OVERLAY LEGEND = legend VAXIS=axis;
                RUN;
                TITLE1; SYMBOL1; SYMBOL2; SYMBOL3; AXIS; LEGEND;
                QUIT;
            %END;
            %IF %EVAL(&plot = 3) %THEN %DO;
                PROC TEMPLATE; /* More elegant graph style, requiring Java; */                                                                                                     
                    DEFINE STATGRAPH mygraph_tvem;                                                                                                             
                        BEGINGRAPH;                                                                                                                           
                            ENTRYTITLE 
                "Odds Ratio and 95% Confidence Band for Covariate: &var_name";                                                                                                             
                            LAYOUT LATTICE / COLUMNDATARANGE=union ;
                                /* Overlay a Bandplot and a Seriesplot */                                                                                               
                                LAYOUT OVERLAY /
                                    YAXISOPTS=(label="Odds Ratio of &var_name");                                                                               
                                    BANDPLOT X=&time 
                                             LIMITUPPER=exp_&var_name._U 
                                             LIMITLOWER=exp_&var_name._L /
                                                 DATATRANSPARENCY=.3 MODELNAME="fit" 
                                             LEGENDLABEL="95% Confidence Limits";                                                         
                                    SERIESPLOT x=&time y=exp_&var_name;  
                                    REFERENCELINE y=1/LINEATTRS=(pattern=4);   
                                ENDLAYOUT;                               
                            ENDLAYOUT;                                                                                                                         
                        ENDGRAPH;                                                                                                                             
                    END; 
                RUN;   
                PROC SGRENDER DATA=_plot_data_OR TEMPLATE=mygraph_tvem;
                RUN;
            %END;                         
        %END; 
    %END;  /* loop over time-varying effects */
	TITLE3;
%MEND;

