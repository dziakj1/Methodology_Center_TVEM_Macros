   /*------------------------------------------------------------------*
   | Tvem_v311.sas
   | TVEM (Time-Varying Effect Model) SAS Macro Suite for Generalized Linear Models.
   | 
   | Revised by John DZIAK (10/10/2017) 
   |
   | Previous version created by  : Runze LI, John DZIAK, Xianming TAN, 
   |                                Jingyun YANG and Liying HUANG
   |                                (9/24/2015)
   |
   | Copyright 2017 The Methodology Center, The Pennsylvania State University.
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
   | Special thanks are owed to Dr. Deborah Kloska and Dr. Sara Vasilenko 
   | for their assistance in testing and finding limitations in previous versions.
   |
   | Example of use : Suppose the dataset data1 has the following variables: 
   |        subj: subjects' identification variable; 
   |        t:    measurement time for each observation, 
   |              different subjects may have different measurement times;
   |        y:    the dependent variable;
   |        x0:   the intercept variable which equals 1 for all observations;  
   |        x1:   a binary covariate;
   |        x2:   a continuous covariate.  
   |    Suppose some observations in this data set are as follows:
   |                subj   y      t      X0     X1         X2                  
   |                -------------------------------------------------
   |                1      4    0.00591     1     0    -1.71718
   |                1      8    0.01378     1     0    -1.99814
   |                1      0    0.01772     1     1     0.59689
   |                1      3    0.02559     1     1    -0.21331
   |                1      3    0.02953     1     1     1.10735
   |                2      1    0.00197     1     1     1.28186
   |                2      0    0.00591     1     1     1.55535
   |                2      5    0.00984     1     0    -0.84808
   |                2      2    0.01378     1     1     0.80863
   |                2      2    0.02165     1     0     0.60380
   |                .....................................
   |   Suppose we want to fit a time-varying coefficient model in which y is the 
   |   dependent variable (count response), and x0, x1 and x2 are the covariates, 
   |   and the coefficient of x0, x1 and x2 may vary over time.  
   |   We could call the macro as follows:
   |      %TVEM(      
   |        data = data1,  
   |        id = subj,          
   |        time = t,         
   |        dv = y,            
   |        tvary_effect = x0 x1 x2,        
   |        method = b-spline,
   |        knots = 5 5 5,
   |        dist = poisson 
   |        );
   *------------------------------------------------------------------*/


/**************************************************************************/
/*                         The Main Macro                                 */
/**************************************************************************/
%MACRO TVEM(data ,
            id ,
            time ,
            dv ,
            tvary_effect ,
            method ,
            knots , 
            dist = ,
            degree = 3 ,
            evenly = 0 ,
            invar_effect = ,
            output_prefix = tvem_ ,
            outfilename = ,
            plot = full ,
            plot_scale = 100 ,
            random = none ,
            stderr = robust 
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
        %IF %NRBQUOTE(&id) =  %THEN %DO;
            PRINT("Error in macro:" //
                  "Please specify an ID variable in id=.");
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
                  "Please specify at least one time-varying covariate in tvary_effect=.");
            ready_to_go = 0;
        %END;
        %IF %NRBQUOTE(&knots) =  %THEN %DO; 
            PRINT("Error in macro:" //
                   "Please specify an integer for each time-varying covariate in knots=.");
            ready_to_go = 0;
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
        IF (( UPCASE(STRIP("&dist")) = "POISSON" ) )
            THEN DO;
            internal_dist_option = 3;
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
        /* Method option */ 
        internal_method_option = -1;
        IF (( UPCASE(STRIP("&method")) = "B" ) |
            ( UPCASE(STRIP("&method")) = "BSPLINE" ) |
            ( UPCASE(STRIP("&method")) = "B-SPLINE" ) |
            ( UPCASE(STRIP("&method")) = "B_SPLINE" ) |
            ( UPCASE(STRIP("&method")) = "UNPENALIZED" )  ) THEN DO;
            internal_method_option = 1;
        END;
        IF (( UPCASE(STRIP("&method")) = "P" ) | 
            ( UPCASE(STRIP("&method")) = "PSPLINE" ) | 
            ( UPCASE(STRIP("&method")) = "P-SPLINE" ) | 
            ( UPCASE(STRIP("&method")) = "P_SPLINE" ) |
            ( UPCASE(STRIP("&method")) = "PENALIZED" )) THEN DO;
            internal_method_option = 2;
        END;
        IF (internal_method_option =-1)  THEN DO;
            PRINT "Error in macro:  The method option was not recognized.";
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
        /* Random effects option */
        internal_random_option = -1;
        IF (( UPCASE(STRIP("&random")) = "NO" ) |
            ( UPCASE(STRIP("&random")) = "NONE") |
            ( UPCASE(STRIP("&random")) = "INDEPENDENT")) THEN DO;
            internal_random_option = 1;
        END; 
        IF ( UPCASE(STRIP("&random")) = "INTERCEPT" ) THEN DO;
            internal_random_option = 2;
        END; 
        IF ( UPCASE(STRIP("&random")) = "SLOPE" ) THEN DO;
            internal_random_option = 3;
        END; 
        IF (internal_random_option = -1) THEN DO;
            PRINT "Error in macro:  The option for random was not recognized.";
            ready_to_go = 0;
        END;
        /* Plot Scale */
        plot_scale = &plot_scale;
        IF ((plot_scale < 100) | (plot_scale > 10000)) THEN DO;
            PRINT("Warning in macro:" // 
                  "We recommend a plot_scale value between 100 and 10000 for the plot."); 
        END;
        /* Standard error option */ 
        internal_se = -1;
        IF (( UPCASE(STRIP("&stderr")) = "MODEL" ) |
            ( UPCASE(STRIP("&stderr")) = "MODEL-BASED" )) THEN DO;
            internal_se = 1;
        END;
        IF (( UPCASE(STRIP("&stderr")) = "ROBUST" ) | 
            ( UPCASE(STRIP("&stderr")) = "SANDWICH" ))
            THEN DO;
            internal_se = 2;
        END;
        IF (( UPCASE(STRIP("&stderr")) = "MBN" ) | 
            ( UPCASE(STRIP("&stderr")) = "CORRECTED" ))
            THEN DO;
            internal_se = 3;
        END;
        IF (internal_se = -1) THEN DO;
            PRINT("Error in macro:" //
                  "Standard error option was not recognized.");
            ready_to_go = 0;
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
                IF internal_dist_option = 3 THEN DO; 
                    num_negative = SUM(nonmissingdep < 0);
                    num_noninteger = SUM(nonmissingdep - INT(nonmissingdep) ^= 0);
                    IF ((num_negative>0) | (num_noninteger>0)) THEN DO;
                        PRINT("Error in macro:" // 
                              "The &dv variable should contain only nonnegative integers " //
                              "for a Poisson regression TVEM.");
                        ready_to_go = 0;
                    END;
                END;
        CLOSE &data;
        /* Get ready to send information to the back-end macro */
        CALL SYMPUT("internal_deg_option", CHAR(internal_deg_option));
        CALL SYMPUT("internal_dist_option", CHAR(internal_dist_option));
        CALL SYMPUT("internal_evenly_option", CHAR(internal_evenly_option));
        CALL SYMPUT("internal_method_option", CHAR(internal_method_option));
        CALL SYMPUT("internal_plot", CHAR(internal_plot));
        CALL SYMPUT("internal_random_option", CHAR(internal_random_option));
        CALL SYMPUT("internal_se", CHAR(internal_se));
        CALL SYMPUT("internal_plot", CHAR(internal_plot));
        CALL SYMPUT("internal_random_option", CHAR(internal_random_option));
        CALL SYMPUT("internal_se", CHAR(internal_se));
        CALL SYMPUT("ready_to_go", CHAR(ready_to_go));
    QUIT;
    /* Unless there is an obvious error, call the appropriate back-end macro. */
    %IF %EVAL(&ready_to_go = 1) %THEN %DO;
        %IF %EVAL(&internal_method_option = 1) %THEN %DO;
            %_TvemB( data  = &data,  
                     id = &id,              
                     time = &time,         
                     dv = &dv,            
                     invar_effect = &invar_effect,            
                     tvary_effect = &tvary_effect,           
                     knots = &knots,   
                     degree = &internal_deg_option,              
                     evenly = &internal_evenly_option,       
                     plot_scale = &plot_scale,
                     dist_option = &internal_dist_option, 
                     plot = &internal_plot, 
                     random_option = &internal_random_option , 
                     stderr = &internal_se, 
                     output_prefix = &output_prefix, 
                     outfilename = &outfilename  
               );
        %END;
        %IF %EVAL(&internal_method_option = 2) %THEN %DO;
            %_TvemP( data  = &data,  
                     id = &id,              
                     time = &time,         
                     dv = &dv,            
                     invar_effect = &invar_effect,            
                     tvary_effect = &tvary_effect,           
                     knots = &knots,   
                     degree = &internal_deg_option, 
                     evenly = &internal_evenly_option, 
                     plot_scale = &plot_scale,
                     dist_option = &internal_dist_option, 
                     plot = &internal_plot, 
                     output_prefix = &output_prefix, 
                     outfilename = &outfilename );
        %END;
    %END;
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
                _glimmix_converged
                _glimmix_covb
                _glimmix_covparms
                _glimmix_est 
                _glimmix_info 
                _glimmix_nobs
                _glimmix_out
                _glimmix_outdesign
                _glimmix_random_est
                _glimmix_stats
                _id_dep_cov
                _modelcovariance
                _model_data
                _model_data_nomiss
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
%MEND;

/**************************************************************************/
/*                                 _TvemB                                 */
/**************************************************************************/
%MACRO _TvemB(data ,  
              id ,              
              time ,         
              dv ,            
              invar_effect ,            
              tvary_effect ,           
              knots ,   
              degree ,              
              evenly ,       
              plot_scale ,
              dist_option ,
   /* 1 = normal, 2 = logistic, 3 = Poisson */
              plot ,
   /* 1 = none, 2 = simple, 3 = polished */
              random_option ,
   /* 1 = independent, 2 = random intercept, 3 = random slope */
              stderr ,
   /* 1 = model-based, 2 = sandwich, 3 = MBN corrected sandwich */
              output_prefix,
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
        KEEP &id &dv &invar_effect &time; 
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
    /* Step 3: Fit the underlying parametric model. */
    /* The most important results of this step are my_Est (the coefficients) 
        and my_CovB (their covariance matrix) */
    DATA model_data;  /* merge  */
        MERGE _id_dep_cov _TVE_spline;
    RUN;
    PROC CONTENTS DATA=_TVE_spline OUT=_tmp NOPRINT;    RUN;
    PROC SQL; /* find out the number of variables reduced from B-splines */
        RESET NOPRINT;
        SELECT (MAX(VARNUM)) INTO :aa FROM _tmp;  
    QUIT;
    %LET size = %TRIM(&aa);
    ODS EXCLUDE ALL; RUN;
    ODS NORESULTS; RUN;
    PROC GLIMMIX NOITPRINT NOCLPRINT MAXOPT=1000 DATA=model_data METHOD=RSPL IC=PQ
    %IF %EVAL(&stderr = 2) %THEN %DO;  EMPIRICAL  %END;
    %IF %EVAL(&stderr = 3) %THEN %DO;  EMPIRICAL=MBN  %END;
        ;
        CLASS &id;
        MODEL &dv 
        %IF %EVAL(&dist_option = 2) %THEN %DO; (DESCENDING) %END;
               = _B1-_B&size &invar_effect / S COVB NOINT 
        %IF %EVAL(&dist_option = 1) %THEN %DO; 
            LINK=identity DIST=Gaussian 
        %END;
        %IF %EVAL(&dist_option = 2) %THEN %DO; 
            LINK=logit DIST=binary
        %END;
        %IF %EVAL(&dist_option = 3) %THEN %DO; 
            LINK=log DIST=Poisson 
        %END; 
        ; 
        %IF %EVAL(&random_option = 2) %THEN %DO; 
            RANDOM INTERCEPT / S SUB=&id; 
        %END; 
        %IF %EVAL(&random_option = 3) %THEN %DO; 
            RANDOM INTERCEPT &time / S SUB=&id; 
        %END; 
        %IF %EVAL((&dist_option = 2)|(&dist_option = 3)) %THEN %DO;  
             RANDOM _RESIDUAL_;  
        %END;
        OUTPUT OUT=_glimmix_out  
               PRED(NOBLUP ILINK)=predicted 
               PRED(BLUP ILINK)=adjusted_predicted 
               STDERR(NOBLUP ILINK)=stderr_predicted 
               STDERR(BLUP ILINK)=adjusted_stderr_predicted;                  
        ODS OUTPUT ConvergenceStatus=_glimmix_Converged
               CovB=_glimmix_CovB
                %IF %EVAL(  (&random_option = 2)|
                            (&random_option = 3)|
                            (&dist_option = 1)) %THEN %DO; 
                       CovParms=_glimmix_CovParms
                %END;
               FitStatistics=_glimmix_Stats  
               ModelInfo=_glimmix_Info
               NObs=_glimmix_NObs 
               ParameterEstimates=_glimmix_est ;  
    RUN;   
    DATA _InvariantEffectsCovariates;
        SET _glimmix_est;
        WHERE ((SUBSTR(Effect,1,2) ^= "_P") & 
               (SUBSTR(Effect,1,2) ^= "_B") & 
               (StdErr ^= .)); 
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
        USE _glimmix_est;
            READ ALL VAR {estimate} INTO mat_fix_est;
        CLOSE _glimmix_est;
        k_s = {&knots}; 
            /* a vector of number of inner knots, one for each 
                covariate in &tvary_effect */
        j = k_s[+] + (&degree+1)*NCOL(k_s);
        c_name = ( compress(  CATT( "Col", char(j, 8, 0)) ) ); 
        USE _glimmix_CovB;
            READ ALL VAR("Col1":c_name) INTO mat_covB;  
            /* only read columns relevant to B-spline basis */
        CLOSE _glimmix_CovB; 
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
             id = &id,              
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
             random_option = &random_option,
             are_not_using_penalty = 1,
             output_prefix = &output_prefix,
             outfilename = &outfilename );
%MEND;

/**************************************************************************/
/*                                 _TvemP                                 */
/**************************************************************************/

%MACRO _TvemP(  data, 
                id, 
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
                output_prefix, 
                outfilename  );
    /* Step 1: Process the input data and create some datasets for intermediate
        calculations.  Some of this is similar to what is done for the B-spline 
        code. However, an additional step is to temporarily the time variable 
        so that it covers approximately the interval (0,1), before converting 
        it into basis functions. We do this because calculating the truncated 
        power spline involves squaring (for deg >= 2) or cubing (for deg >= 3)
        the values of the time variable, and this might cause numerical problems 
        if time is on an arbitrary scale (i.e., if units are extremely 
        large or small). */
    PROC SQL;
       RESET NOPRINT;
       SELECT MIN(&time) INTO :min_t FROM &data;
          SELECT MAX(&time) INTO :max_t FROM &data;
       RESET PRINT;
    QUIT;
    DATA _tmp_mydata;
        SET &data;
        &time = (&time - &min_t + 0.001)/(&max_t - &min_t + 0.002);
    RUN;
    DATA _t_TVE;   /* covariates with time-varying coefficients */
        SET _tmp_mydata;
        WHERE &time >. ;
        KEEP &tvary_effect; 
    RUN;
    DATA _id_dep_cov; 
            /*id, dependent variable, covariates with constant coefficients, 
            &time may be needed for R specification */
        SET _tmp_mydata;
        WHERE &time >. ;
        KEEP &id &dv &invar_effect &time; 
    RUN;
    DATA _t_orig;   /* covariates with time-varying coefficients */
        SET _tmp_mydata;
        WHERE &time >. ;
        KEEP &time; 
    RUN;
    PROC SORT DATA=_t_orig OUT=_t_sorted;  
        /*sorted time variable, needed for inner knots calculation */ 
        BY &time;
    RUN;
    PROC SORT DATA=_t_orig OUT=_t_distinct NODUPKEY;  
        /*sorted time variable, needed for inner knots calculation */ 
        BY &time;
    RUN;
    DATA _t_orig;   
        /* Keep _idx to order the observations, be consistent with _id_dep_cov  */
        RETAIN _idx 0;
        SET _t_orig;
        _idx = _idx + 1;
    RUN;
    /* Step 2: Create the design matrix columns to represent the
     P-spline transformation of &tvary_effect */
    PROC IML;
        /* Part 2.1: Define three functions, 
                     my_Pspline(...), my_vec_Pspline(...), and my_knots(...) */
        START my_Pspline(x, knots, d);
        /*    Evaluate the values of all of the truncated power spline basis  
           functions (there are number of knots + d + 1 of them) at x, 
           where x is a real number */
        /*    Outputs a vector of (num of knots) + d + 1 real numbers. */
        /*  Unlike B-splines, all of the knots are inner knots. */
            n_ext = ncol(knots);
            tmp=0.0*(1:(n_ext+d+1));
            tmp[1] = 1.0;
            DO i = 2 TO (n_ext+d+1);
                IF i<=(d+1) THEN DO;    /* the power part */
                    tmp[i] = x##(i-1);
                END; ELSE DO;   /* the truncated power part */
                    z = x-knots[i-d-1];
                    IF  z>0 THEN tmp[i] = z##d;
                            ELSE tmp[i] = 0;
                 END;
            END; 
            RETURN(tmp[1:(n_ext+d+1)]);
        FINISH;
        START my_vec_Pspline(xx, knots, d);
        /*    evaluate the values of the (num of knots) + d +1 B-spline basis 
            functions at xx */ 
        /*    xx is a vector of real numbers*/
        /*    output a real matrix with dim(xx) rows,  and (num of knots) + d +1 columns */
            n_row = ncol(xx);
            n_col = ncol(knots)+d+1;        
            out = J(n_row, n_col, .);  
            DO i=1 TO n_row; 
                out[i, ] = t(my_Pspline(xx[i], knots, d));
            END;            
            RETURN(out);
        FINISH;
        START my_knots(samp, nknots, evenly, d);  
        /* Calculate nknots inner knots */
        /* the inner knots are uniform in space (evenly=1) or uniformly on 
            quantiles (evenly != 1) */
        /* samp should be sorted ascendingly, i.e., from small to large */
            len = ncol(samp);
            all_knots = 0*(1:nknots);
            IF (evenly=1) THEN DO;
                dist  = (samp[len]-samp[1])/(nknots+1);         
                all_knots[1:nknots] = samp[1]+dist*(1:nknots);
            END; ELSE DO; 
                 j=1; 
                DO i = 1 TO nknots;
                    DO WHILE ( (j-1)/(len-1) < i/(nknots+1) );
                        j=j+1;
                    END;
                    all_knots[i] = samp[j];
                END;
                /* Following the documentation for B-splines in PROC TRANSREG, 
                    quantiles are considered to be on data points, with no 
                    interpolation; */
            END;
            RETURN(all_knots);
        FINISH; 
        /* Part 2.2: generate Pspline matrix, given t_cov and options */
        /*         for knots, and read data into matrix */
        USE _t_sorted;
        READ ALL INTO mat_t_sorted; /* &time only, sorted */
        USE _t_distinct;
        READ ALL INTO mat_t_distinct; 
                /* &time only, distinct time points only */
        n_pt = &plot_scale -1;
        mat_t_dense = t( mat_t_sorted[1]+ ((0:n_pt)/n_pt)*
                ( mat_t_sorted[ nrow(mat_t_sorted) ] - mat_t_sorted[1] ) );
        k_s = {&knots};
                /* a vector of number of inner knots, one for each 
                    covariate in &tvary_effect */
        DO i = 1 TO ncol(k_s);
            all_knots = my_knots( t(mat_t_sorted), k_s[i], &evenly, &degree );
            mat_t_splines = mat_t_splines || 
                            ( my_vec_Pspline( t(mat_t_distinct), 
                                                all_knots, &degree) );
                            /* distinct time points */
            mat_dt_splines = mat_dt_splines || 
                            ( my_vec_Pspline( t(mat_t_dense), 
                                                all_knots, &degree) );
                            /* uniform separated time points */
        END;
        mat_t_splines = (mat_t_distinct || mat_t_splines) ;
        CREATE _t_spline FROM mat_t_splines[colname= "_time"];
            APPEND FROM mat_t_splines;
        CLOSE _t_spline;
        /* output mat_dt_splines to dataset dt_spline
            dataset t_dense is useful for plotting  */
        CREATE _dt_spline FROM mat_dt_splines; 
            APPEND FROM mat_dt_splines;
        CLOSE _dt_spline;
        CREATE _t_dense FROM mat_t_dense [colname= "&time"]; 
                /* denser time points */
            APPEND FROM mat_t_dense;
        CLOSE _t_dense;
    QUIT;
    PROC SQL;
        CREATE TABLE _t_orig_spline AS
        SELECT A.*, B.*  
        FROM _t_orig A, _t_spline B 
        WHERE A.&time = B._time
        ORDER BY _idx;
    QUIT;
    PROC IML;
        k_s = {&knots}; 
            /* a vector of number of inner knots, one for each 
                covariate in &tvary_effect */
        j = k_s[+] + (&degree + 1)*ncol(k_s);
        cname = ( compress(  CATT( "col", char(j+1, 8, 0)) ) ); 
        USE _t_orig_spline;
        READ ALL VAR("col2":cname) INTO mat_tspline; 
            /* &time &tvary_effect, unsorted,*/ 
        USE _t_TVE;
        READ ALL VAR{&tvary_effect} INTO mat_TVE;  /* &time &tvary_effect, unsorted,*/ 
             /* Using VAR{&tvary_effect} makes sure the order of columns in mat_TVE 
               is consistent with &tvary_effect(and hence knots)*/
        ida = 0;
        id_u = 0; 
        DO i = 1 TO ncol(k_s);
            id_l = id_u + 1; id_u = id_l + (k_s[i]+ &degree +1)-1;
            mat_t_splines = mat_t_splines ||
                            ( (mat_tspline[, id_l: id_u] ) # mat_TVE[,i]);
            DO j = 1 TO (&degree+1);
                ida=ida+1;
                c_name = (c_name || COMPRESS ( CATT("_P", CHAR(ida,8,0) ))  );
            END;
            DO j = 1 TO k_s[i];
                c_name = (c_name || COMPRESS ( CATT( "_TVE", 
                            CHAR(i, 8, 0), "TP", CHAR(j,8, 0) ) )); 
            END;
        END;
        /* output mat_TVE_splines to dataset TVE_spline */
        CREATE _TVE_spline FROM mat_t_splines[colname= c_name];
                        /* col1 - col( sum( k_s[i] + deg + 1 )  )*/
        APPEND FROM mat_t_splines; 
    QUIT;
    /* Step 3: Fit the model */
    DATA _model_data;  /* merge  */
        MERGE _id_dep_cov _TVE_spline;
    RUN; 
    PROC IML;
        NumberOfTVEs = COUNT(COMPBL("&tvary_effect"), " ")+1;
        TotalNumberUnpenCoefsForTVEs = NumberOfTVEs *( &degree +1); 
        CALL SYMPUT("NumberOfTVEs",CHAR(NumberOfTVEs,8));
        CALL SYMPUT("TotalNumberUnpenCoefsForTVEs",
                    CHAR(TotalNumberUnpenCoefsForTVEs,8));
        DO i = 1 TO NumberOfTVEs;
            NumberOfKnotsForThisTVE = SCAN("&knots", i);
            TempString = CATT("NumberOfKnotsForTVE",i);
            CALL SYMPUT(TempString,NumberOfKnotsForThisTVE);
        END;
    QUIT;   
    DATA _model_data_nomiss;
        SET _model_data(WHERE=((&dv IS NOT MISSING) & (&id IS NOT MISSING)
            %DO i = 1 %TO &TotalNumberUnpenCoefsForTVEs;
                       & (_P&i IS NOT MISSING )
            %END;
            %DO i = 1 %TO &NumberOfTVEs;
                %DO j = 1 %TO &&NumberOfKnotsForTVE&i;
                       & (_TVE&i.TP&j IS NOT MISSING)
                %END;
            %END;      ));
    RUN;
    DATA _model_data_nomiss;
        SET _model_data_nomiss;
        %IF %NRBQUOTE(&invar_effect) =  %THEN %DO;
            _Num_Invar_Effects_Missing = 0;
            /* No invariant effect covariates */
        %END; %ELSE %DO;
            %LET list_of_invar_effect =  %SYSFUNC(TRANSLATE(%QCMPRES(&invar_effect),","," "));
            %PUT &invar_effect; 
            %PUT &invar_effect; 
            %PUT &list_of_invar_effect; 
            %PUT &list_of_invar_effect;  
            _Num_Invar_Effects_Missing = NMISS(&list_of_invar_effect);
        %END;
    RUN;
    DATA _model_data_nomiss;
        SET _model_data_nomiss;
        WHERE _Num_Invar_Effects_Missing = 0;
    RUN;
    DATA _model_data_nomiss;
        SET _model_data_nomiss;
        DROP _Num_Invar_Effects_Missing;
    RUN;
    ODS EXCLUDE ALL; RUN;
    ODS NORESULTS; RUN;
    PROC GLIMMIX NOITPRINT NOCLPRINT 
                 OUTDESIGN(NOVAR)=_glimmix_outdesign
                 DATA=_model_data_nomiss METHOD=RSPL  ;
                   /* Note: The METHOD=RSPL is ignored unless there is
                      nonnormal data with random effects. */
        CLASS &id; 
        MODEL &dv 
            %IF %EVAL(&dist_option = 2) %THEN %DO; (DESCENDING) %END; = 
        /* Unpenalized X treated as fixed effects: */
            %DO i = 1 %TO &TotalNumberUnpenCoefsForTVEs;
                        _P&i
            %END;
            &invar_effect / NOINT S COVB 
            %IF %EVAL(&dist_option = 1) %THEN %DO;
                LINK=identity DIST=Gaussian 
            %END;
            %IF %EVAL(&dist_option = 2) %THEN %DO; 
                LINK=logit DIST=binary
            %END;
            %IF %EVAL(&dist_option = 3) %THEN %DO; 
                LINK=log DIST=Poisson 
            %END; 
            DDFM=BW;
        /* Penalized X treated like random effects: */
            %DO i = 1 %TO &NumberOfTVEs;
              RANDOM  _TVE&i.TP1 - _TVE&i.TP&&NumberOfKnotsForTVE&i / G S TYPE=TOEP(1);
            %END;    
        %IF %EVAL((&dist_option = 2)|(&dist_option = 3)) %THEN %DO;  
             RANDOM _RESIDUAL_;  
        %END;
        NLOPTIONS TECHNIQUE=NRRIDG;  
             /* A ridge penalty is needed here to handle possible unstable solutions.
                This was done automatically in PROC MIXED and the GLIMMIX macro. */
        OUTPUT OUT = _glimmix_out   
               PRED(BLUP ILINK) = _mu  ;     
        ODS OUTPUT ConvergenceStatus = _glimmix_Converged 
                   CovB = _glimmix_covb
                   CovParms = _glimmix_CovParms  
                   ModelInfo = _glimmix_Info 
                   NObs = _glimmix_NObs 
                   ParameterEstimates = _glimmix_est 
                   SolutionR = _glimmix_random_est  ;  
    RUN;
    DATA _NULL_;
        SET _glimmix_Converged; 
        CALL SYMPUT("FirstTryStatus",Status);
    RUN;
    %IF %EVAL(&FirstTryStatus ^= 0) %THEN %DO; 
         /* 0 here stands for success, so we are handling the 
            case where REML does not converge.  We will try a    
            cruder method instead. */
        PROC DATASETS NOLIST NOWARN;
            DELETE _glimmix_Converged _glimmix_covb
                   _glimmix_CovParms _glimmix_Info 
                   _glimmix_NObs _glimmix_est 
                   _glimmix_random_est;
        QUIT;           
        DATA _NULL_;
            PUTLOG '-----------------------------------------------' ;
            PUTLOG 'Having difficulty estimating optimal penalty size.;';
            PUTLOG 'Trying again with NOITER.;' ;
            PUTLOG '-----------------------------------------------;' ;
        RUN;
        PROC GLIMMIX NOITPRINT NOCLPRINT 
                     OUTDESIGN(NOVAR)=_glimmix_outdesign
                     DATA=_model_data_nomiss METHOD=RSPL ;
            CLASS &id; 
            MODEL &dv 
                %IF %EVAL(&dist_option = 2) %THEN %DO; (DESCENDING) %END; = 
            /* Unpenalized X treated as fixed effects: */
                %DO i = 1 %TO &TotalNumberUnpenCoefsForTVEs;
                            _P&i
                %END;
                &invar_effect / NOINT S COVB 
                %IF %EVAL(&dist_option = 1) %THEN %DO;
                    LINK=identity DIST=Gaussian 
                %END;
                %IF %EVAL(&dist_option = 2) %THEN %DO; 
                    LINK=logit DIST=binary
                %END;
                %IF %EVAL(&dist_option = 3) %THEN %DO; 
                    LINK=log DIST=Poisson 
                %END; 
                DDFM=BW;
            /* Penalized X treated like random effects: */
                %DO i = 1 %TO &NumberOfTVEs;
                  RANDOM  _TVE&i.TP1 - _TVE&i.TP&&NumberOfKnotsForTVE&i / G S TYPE=TOEP(1);
                %END;
            %IF %EVAL((&dist_option = 2)|(&dist_option = 3)) %THEN %DO;  
                 RANDOM _RESIDUAL_;  
            %END;
            OUTPUT OUT = _glimmix_out   
                   PRED(BLUP ILINK) = _mu  ; 
            PARMS / NOITER ; /* like MIVQUE0 */   
            ODS OUTPUT ConvergenceStatus = _glimmix_Converged 
                       CovB = _glimmix_covb
                       CovParms = _glimmix_CovParms  
                       ModelInfo = _glimmix_Info 
                       NObs = _glimmix_NObs 
                       ParameterEstimates = _glimmix_est 
                       SolutionR = _glimmix_random_est  ;  
        RUN;        
    %END; 
    DATA _glimmix_Converged;
        SET _glimmix_Converged;
        EarlyStatus = &FirstTryStatus;
    RUN;     
    ODS EXCLUDE NONE; RUN;
    ODS RESULTS; RUN;
    /* Step 3.25:  Create some datasets to use in calculating 
        standard errors later */
    DATA _TvemP_Y;
        SET _model_data_nomiss;
        KEEP &dv;
    RUN;
    DATA _TvemP_mu;
        SET _glimmix_out;
        KEEP _mu;
    RUN;
    DATA _TvemP_temp_id;
        SET _model_data_nomiss;
        KEEP &id;
    RUN;
    PROC SORT DATA=_TvemP_temp_id; BY &id; RUN;
    DATA _TvemP_id;
        SET _TvemP_temp_id;
        BY &id;
        RETAIN _integer_id 0;
        IF first.&ID THEN _integer_id = _integer_id + 1;
            /* See http://www.ats.ucla.edu/stat/sas/faq/enumerate_id.htm; */
        OUTPUT; 
        KEEP _integer_id; 
    RUN;
    /*Step 3.5:  Calculate sandwich covariance estimator: */
    PROC IML;
        dist = &dist_option;
        USE _dt_spline;
            READ ALL INTO mat_dt_spline;  
            /* Spline basis functions at denser time points */
        CLOSE _dt_spline;
        USE _glimmix_est; 
            READ ALL VAR {Estimate} INTO b_fixed;
        CLOSE _glimmix_est;
        USE _glimmix_est; 
            READ ALL VAR {StdErr} INTO b_fixed_model_se; 
        CLOSE _glimmix_est;
        USE _glimmix_random_est; 
            READ ALL VAR {Estimate} INTO b_random; 
        CLOSE _glimmix_random_est;
        USE _tvemp_id; 
            READ ALL INTO id; 
        CLOSE _tvemp_id;
        USE _tvemp_mu; 
            READ ALL INTO mu; 
        CLOSE _tvemp_mu;
        USE _dt_spline; 
            READ ALL INTO DenseSplineBases; 
        CLOSE _dt_spline;
        USE _glimmix_OutDesign; 
            READ ALL INTO Design; 
        CLOSE _glimmix_OutDesign;
        USE _tvemp_y; 
            READ ALL INTO y; 
        CLOSE _tvemp_y; 
        USE _glimmix_CovParms; 
            READ ALL VAR {Estimate} INTO PenaltyVariance WHERE 
                (CovParm = "Variance");
        CLOSE _glimmix_CovParms;   
        IF (dist = 1) THEN DO;
            USE _glimmix_CovParms; 
                READ ALL VAR {Estimate} INTO ScaleVariance WHERE 
                (CovParm = "Residual");
            CLOSE _glimmix_CovParms;  
        END; ELSE DO;
            ScaleVariance = 1;
        END;
        p = J(NROW(b_fixed),1,0) %DO i = 1 %TO &NumberOfTVEs; // 
                J(&&NumberOfKnotsForTVE&i,1,ScaleVariance/(PenaltyVariance[&i]+1e-15))
            %END; ;
        nsub = id[<>];
        nobs = NROW(Design);
        npar = NCOL(Design);
        XDX = J(npar,npar,0); 
        b = b_fixed // b_random;
        DO i = 1 TO nsub;
            IF (SUM(id=i)>0) THEN DO;
                these = LOC(id=i); 
                IF (dist = 1) THEN DO;
                    deriv = 1;
                END;
                IF (dist = 2) THEN DO;
                    deriv = (mu[these,])#(1-mu[these]);
                END;
                IF (dist = 3) THEN DO;
                    deriv = mu[these,];
                END;
                XDX = XDX + T(Design[these,])*DIAG(deriv)*Design[these,];  
            END; 
        END;   
        bread = INV(XDX + DIAG(p));    
        meat = J(npar,npar,0);
        DO i = 1 TO nsub;
            IF (SUM(id=i)>0) THEN DO;
                these = LOC(id=i);
                meat = meat + T(Design[these,])*(y[these]-mu[these])*
                              T(y[these]-mu[these])*Design[these,]; 
            END; 
        END; 
        sandwich = bread * meat * bread; 
        CREATE _ModelCovariance FROM bread;     
            APPEND FROM bread;
        CLOSE _ModelCovariance;
        CREATE _RobustCovariance FROM sandwich;     
            APPEND FROM sandwich;
        CLOSE _RobustCovariance; 
    /*step 4: calculate curves and confidence intervals */   
        StartIndexUnpenalizedCoefs = 0;
        StartIndexPenalizedCoefs = 0;
        StartIndexDenseSplineBases = 0;  
        %IF %EVAL(&NumberOfTVEs > 0) %THEN %DO;
            %DO i = 1 %TO &NumberOfTVEs; 
                NumberOfFixedCoefsThisTVE = &degree + 1;
                NumberOfRandomCoefsThisTVE = &&NumberOfKnotsForTVE&i;    
                DesignIndicesThisTVE = ((StartIndexUnpenalizedCoefs +1):
                    (StartIndexUnpenalizedCoefs + 
                        NumberOfFixedCoefsThisTVE)) || 
                     ((NROW(b_fixed)+StartIndexPenalizedCoefs +1):
                        (NROW(b_fixed)+StartIndexPenalizedCoefs + 
                            NumberOfRandomCoefsThisTVE));
                DenseSplineIndicesThisTVE = (StartIndexDenseSplineBases+1):
                    (StartIndexDenseSplineBases + 1 + &degree +
                    &&NumberOfKnotsForTVE&i);
                StartIndexUnpenalizedCoefs = StartIndexUnpenalizedCoefs + 
                        &degree + 1;
                StartIndexPenalizedCoefs = StartIndexPenalizedCoefs +
                        &&NumberOfKnotsForTVE&i;
                StartIndexDenseSplineBases = StartIndexDenseSplineBases + 
                        &degree + 1 + &&NumberOfKnotsForTVE&i;  
                CoefsForThisTVE = b[DesignIndicesThisTVE]; 
                  CovarianceMatrixColumnsThisTVE = 
                        sandwich[DesignIndicesThisTVE,DesignIndicesThisTVE];
                SplineBasesForThisTVE = 
                        DenseSplineBases[,DenseSplineIndicesThisTVE];
                Coef = J(NROW(DenseSplineBases),1,.);
                Coef_StdErr = J(NROW(DenseSplineBases),1,.);
                DO j = 1 TO NROW(DenseSplineBases);
                    Coef[j] = SplineBasesForThisTVE[j,]*CoefsForThisTVE;
                    PointwiseVariance = SplineBasesForThisTVE[j,]*
                                        CovarianceMatrixColumnsThisTVE*
                                        T(SplineBasesForThisTVE[j,]);
                    IF PointwiseVariance > 0 THEN Coef_StdErr[j] = 
                        SQRT(PointwiseVariance);
                    ELSE Coef_StdErr[j] = .;
                END;        
                Coef_L = Coef - 1.96*Coef_StdErr;
                Coef_U = Coef + 1.96*Coef_StdErr; 
                mat_coef = mat_coef || 
                            (( (coef_L || coef) || coef_U) || Coef_StdErr); 
                var_name = SCAN("&tvary_effect", &i);
                var_names = ( CATT(var_name, "_L") || var_name || 
                            CATT(var_name,"_U") || CATT(var_name,"_SE") );
                cname = cname || var_names;
             %END; 
            CREATE _dense_grid_coef FROM mat_coef[colname=cname]; 
                        /* denser time points */
                    APPEND FROM mat_coef;
            CLOSE _dense_grid_coef;
            CREATE _dense_grid_covb FROM sandwich; 
                    /* denser time points */
                APPEND FROM sandwich;
            CLOSE _dense_grid_covb;
        %END;
        /* Step 5:  Show fixed effects covariates */
        %IF %NRBQUOTE(&invar_effect) ^=  %THEN %DO;  /* This line corrects a bug in 3.0.3 */
            StartLeftoversIndex = &NumberOfTVEs * (&degree + 1) + 1;
            IF (StartLeftoversIndex < NROW(b_fixed)) THEN DO;
                temp = VECDIAG(sandwich[StartLeftoversIndex:NROW(b_fixed),
                                        StartLeftoversIndex:NROW(b_fixed)]); 
                IF (SUM(temp<0)>0) THEN DO;
                    temp[LOC(temp<0)] = .;
                END;
            END;        
            IF (StartLeftoversIndex = NROW(b_fixed)) THEN DO;
                temp = sandwich[StartLeftoversIndex:NROW(b_fixed),
                                        StartLeftoversIndex:NROW(b_fixed)]; 
                IF (temp<0) THEN DO;
                    temp = .;
                END;
            END;         
            InvariantEffectsCovariatesSEs = SQRT(temp);
            CREATE _InvariantEffectsCovariatesSEs FROM 
                    InvariantEffectsCovariatesSEs [COLNAME = 'RobustStdErr'];
                APPEND FROM InvariantEffectsCovariatesSEs;
            CLOSE _InvariantEffectsCovariatesSEs;
        %END;
    QUIT;
    %IF %NRBQUOTE(&invar_effect) ^=  %THEN %DO;
        DATA _InvariantEffectsCovariates;
            SET _glimmix_est;
            WHERE ((SUBSTR(Effect,1,2) ^= "_P") AND (StdErr ^= .));
            DROP DF tValue Probt;
        RUN;
        DATA _InvariantEffectsCovariates;
            MERGE _InvariantEffectsCovariates _InvariantEffectsCovariatesSEs;
        RUN;
        DATA _InvariantEffectsCovariates;
            SET _InvariantEffectsCovariates;
            RobustZ = Estimate / RobustStdErr;
            RobustP = 2*( 1 - CDF('Normal',ABS(RobustZ)) );
			DROP StdErr;
        RUN;
    %END;
    /* Step Six:  Print the output */
    DATA _glimmix_est;
        LENGTH Effect $15;  
        SET _glimmix_est _glimmix_random_est;
    RUN;
    DATA _plot_data;
        MERGE _t_dense _dense_grid_coef;
    RUN;  
    DATA _plot_data;
        SET _plot_data;
        &time = &min_t + &time*( &max_t -&min_t + 0.002) - 0.001;
    RUN;
    %_TvemOut(data = &data,  
             id = &id,              
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
             random_option = 0,
             are_not_using_penalty = 0,
             output_prefix = &output_prefix,
             outfilename = &outfilename );
%MEND;

/**************************************************************************/
/*                                 _TvemOut                               */
/**************************************************************************/
%MACRO _TvemOut(data,  
              id,              
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
              random_option,
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
        %DO i = 1 %TO &cnt;                                                                                                                
        %LET var=%SCAN(&list,&i);                                                                                                         
        exp_&var=EXP(&var);                                                                                                               
        %END;     
        &time._old=&time;  
        DROP exp_&time;
    RUN;           
    /* Provide output on the screen or listing */
    PROC IML;
        OPTIONS FORMCHAR="|----|+|---+=|-/\<>*";
        TITLE2 "TVEM Macro Output Summary";
        FILE PRINT;
        PUT "================================================================";
        PUT "Time-Varying Effects Modeling (TVEM) Macro Output";
        PUT "================================================================";
        PUT "Dataset:                  &data";
        PUT "Time variable:            &time";
        PUT "Response variable:        &dv";
        %IF %EVAL(&dist_option = 1) %THEN %DO; 
            PUT "Response distribution:    Normal (Gaussian)";
        %END;
        %IF %EVAL(&dist_option = 2) %THEN %DO; 
            PUT "Response distribution:    Binary (Logistic)";
        %END;
        %IF %EVAL(&dist_option = 3) %THEN %DO; 
            PUT "Response distribution:    Poisson";
        %END;
        %IF %EVAL(&random_option = 1) %THEN %DO; 
            PUT "Random effects:           None";
        %END;
        %IF %EVAL(&random_option = 2) %THEN %DO; 
            PUT "Random effects:           Intercept";
        %END;
        %IF %EVAL(&random_option = 3) %THEN %DO; 
            PUT "Random effects:           Intercept &time";
        %END;
        %IF %TRIM("&invar_effect") ^= "" %THEN %DO;
            PUT "Non-time-varying effects: &invar_effect";
        %END;
        %IF %TRIM("&tvary_effect") ^= "" %THEN %DO;
            PUT "Time-varying effects:     &tvary_effect";
            PUT "Knots for splines:        &knots";
            PUT "Degree for splines:       &degree";
        %END;
        %IF %SYSFUNC(EXIST(_glimmix_nobs)) %THEN %DO;
            USE _glimmix_nobs;
                READ VAR {NObsUsed NObsRead};
            CLOSE _glimmix_nobs;
            CALL SYMPUT("NObsUsed",CHAR(NObsUsed));
            PUT "Number of observations used:  &NObsUsed";
        %END; 
        %IF %SYSFUNC(EXIST(_glimmix_converged)) %THEN %DO;
            USE _glimmix_converged;
                READ VAR {Reason Status };
            CLOSE _glimmix_converged;
            IF (Status = 0) THEN DO;
                PUT "Internal computations in PROC GLIMMIX converged successfully.";
            END; ELSE DO;
                PUT "Internal computations in PROC GLIMMIX did not converge successfully.";
                PUT "Message:" Reason;
            END; 
            %IF (&are_not_using_penalty) %THEN %DO;
            IF (SUM(CONTENTS(_glimmix_converged)="pdG")>0) THEN DO;
                USE _glimmix_converged;
                    READ VAR {Reason pdG};
                CLOSE _glimmix_converged;
                IF (pdG = 0) THEN DO; 
                    PRINT("Error in macro:  The G matrix in PROC GLIMMIX is not" //
                          " positive definite.  Your model might be overspecified" //
                          " (too complicated to interpret correctly given the data).");
                    PRINT("Message from PROC GLIMMIX:" // Reason);      
                END;
            END;
            %END;
        %END; 
        PUT "================================================================";
        %IF %SYSFUNC(EXIST(_glimmix_stats)) %THEN %DO;
            USE _glimmix_stats;
               READ ALL VAR {Descr} INTO Descr;
               READ ALL VAR {Value} INTO Value;
            CLOSE _glimmix_stats; 
            %IF (&are_not_using_penalty) %THEN %DO; 
                IF (SUM((Descr = "-2 Log Likelihood") |
                         (Descr = "-2 Res Log Likelihood"))>0) THEN DO; 
                    NegTwoLogLik = Value[LOC(
                        ((Descr = "-2 Log Likelihood") |
                         (Descr = "-2 Res Log Likelihood")))];
                    LogLik = NegTwoLogLik / (-2);
                    PUT "Log-likelihood (less negative=better fit):        " LogLik 18.3;
                END;
                IF (SUM((Descr = "-2 Log Pseudo-Likelihood") |
                         (Descr = "-2 Res Log Pseudo-Likelihood"))>0) THEN DO; 
                    NegTwoPseudoLogLik = Value[LOC( 
                        ((Descr = "-2 Log Pseudo-Likelihood") |
                         (Descr = "-2 Res Log Pseudo-Likelihood")))];
                    PseudoLogLik = NegTwoPseudoLogLik / (-2);
                    PUT "Pseudo log-likelihood (less negative=better fit): " PseudoLogLik 18.3;
                END;
                IF (SUM((Descr = "AIC  (smaller is better)"))>0) THEN DO;
                    AIC = Value[LOC( 
                        (Descr = "AIC  (smaller is better)"))];
                    PUT "AIC (smaller=better fit):                         " AIC 18.3;
                END;
                IF (SUM((Descr = "BIC  (smaller is better)"))>0) THEN DO;
                    BIC = Value[LOC( 
                        (Descr = "BIC  (smaller is better)"))];
                    PUT "BIC (smaller=better fit):                         " BIC 18.3;
                END; 
                IF (SUM((Descr = "Pseudo-AIC"))>0) THEN DO;
                    AIC = Value[LOC( 
                        (Descr = "Pseudo-AIC"))];
                    PUT "Pseudo-likelihood AIC (smaller=better fit):       " AIC 18.3;
                END;
                IF (SUM((Descr = "Pseudo-BIC"))>0) THEN DO;
                    BIC = Value[LOC( 
                        (Descr = "Pseudo-BIC"))];
                    PUT "Pseudo-likelihood BIC (smaller=better fit):       " BIC 18.3;
                END;
            %END; %ELSE %DO;
                IF (SUM((Descr = "-2 Log Likelihood") |
                         (Descr = "-2 Res Log Likelihood"))>0) THEN DO; 
                    PseudoNegTwoLogLik = Value[LOC( 
                        ((Descr = "-2 Log Likelihood") |
                         (Descr = "-2 Res Log Likelihood")))];
                    PseudoLogLik = PseudoNegTwoLogLik / (-2);
                    PUT "Penalized pseudo log-likelihood (less negative=better fit): " PseudoLogLik 18.3;
                END;
            %END;
            PUT "================================================================";
        %END;
        %IF %SYSFUNC(EXIST(_glimmix_est)) %THEN %DO;
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
    PROC DATASETS NOLIST NOWARN;
        DELETE &output_prefix.invariant_effects 
               &output_prefix.converged
               &output_prefix.covar_orig_coefs
               &output_prefix.covparms
               &output_prefix.fitstats
               &output_prefix.orig_coefs
               &output_prefix.plot_data
               &output_prefix.plot_data_OR
               &output_prefix.predicted;
    QUIT;
    %IF %SYSFUNC(EXIST(_InvariantEffectsCovariates)) %THEN %DO;
        %IF %NRBQUOTE(&invar_effect) ^=  %THEN %DO;
            PROC PRINT DATA=_InvariantEffectsCovariates; 
                TITLE1 "Time-Invariant Effects Covariates"; 
            RUN; TITLE1;
        %END;
        DATA &output_prefix.invariant_effects; 
            SET _InvariantEffectsCovariates; 
        RUN;
    %END;   
    %IF (&are_not_using_penalty) %THEN %DO;
        %IF %SYSFUNC(EXIST(_glimmix_covparms)) %THEN %DO;
            PROC PRINT DATA=_glimmix_covparms;
                TITLE1 "Covariance Parameters";
                WHERE CovParm ^= "Residual (VC)";
            RUN; TITLE1;
            DATA &output_prefix.covparms; SET _glimmix_covparms; RUN;
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
    %IF %SYSFUNC(EXIST(_glimmix_converged)) %THEN %DO;
        DATA &output_prefix.converged; SET _glimmix_converged; RUN;
    %END;
    %IF %SYSFUNC(EXIST(_glimmix_stats)) %THEN %DO;
        DATA &output_prefix.fitstats; SET _glimmix_stats; RUN;
    %END; 
    %IF %SYSFUNC(EXIST(_glimmix_est)) %THEN %DO;
        DATA &output_prefix.orig_coefs; SET _glimmix_est; RUN;
    %END; %ELSE %DO;
        PROC IML;  
            PRINT "Error in macro:  Estimation was not successful."; 
        QUIT;
    %END;
    %IF %SYSFUNC(EXIST(_dense_grid_covb)) %THEN %DO;
        DATA &output_prefix.covar_orig_coefs; SET _dense_grid_covb; RUN;
    %END;
    %IF %SYSFUNC(EXIST(_glimmix_info)) %THEN %DO;
        DATA &output_prefix.glimmix_info; SET _glimmix_info; RUN;
    %END; 
    %IF %SYSFUNC(EXIST(_glimmix_out)) %THEN %DO;
        DATA &output_prefix.predicted; SET _glimmix_out; RUN;
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
%MEND;
