LIBNAME here "C:\Users\jjd264\Documents\Tvem-R\TVEM_v311";
%INCLUDE "C:\Users\jjd264\Documents\Tvem-R\TVEM_v311\Tvem_v311.sas";
 
/* Use simulated data */ 
DATA exampledata;
        SET here.exampledata;
        * Create an intercept column:;
        Intercept = 1;
        * Dummy-code the three-level location variable:;
        Location1 = 0; IF Location=1 THEN Location1 = 1;
        Location2 = 0; IF Location=2 THEN Location2 = 1;
RUN;

 

%TVEM(  dist=normal,
        data = exampledata,
        id = SubjectID,
        time = Time,
        dv = Urge,
        tvary_effect = Intercept,
        method = P-spline,
        knots = 10,    
        invar_effect = Location1 Location2 Male);

%TVEM(  dist=normal,
        data = exampledata,
        id = SubjectID,
        time = Time,
        dv = Urge,
        tvary_effect = Intercept,
        method = B-spline,
        knots = 2,
        random = slope,
        invar_effect = Location1 Location2 Male);

%TVEM(  dist=normal,
        data = exampledata,
        id = SubjectID,
        time = Time,
        dv = Urge,
        tvary_effect = Intercept NegAffect,
        method = P-spline,   
        knots = 10 10,
        invar_effect = Location1 Location2 Male);

%TVEM(  dist = logistic,
        data = exampledata,
        id = SubjectID,
        time = Time,
        dv = ExpectedSuccess,
        tvary_effect = Intercept Urge,
        method = P-spline,
        knots = 10 10 ,
        invar_effect = Location1 Location2 Male);

%TVEM(  dist=logistic,
        data = exampledata,
        id = SubjectID,
        time = Time,
        dv = ExpectedSuccess,
        tvary_effect = Intercept Urge,
        method = B-spline,
        knots = 2 2 ,
        random = slope,
        invar_effect = Location1 Location2 Male);
