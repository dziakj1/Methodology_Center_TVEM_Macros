%MACRO OrdinalTvemLoop( dataset=, /* The name of the SAS dataset containing the input */
	                  id=, /* The name of the subject ID variable */
	                  t=,  /* The name of the observation time variable */
 	                  cov=, /* The name or names of the predictor variables with constant effects*/
 	                  tcov=, /* The name or names of the predictor variables with varying effects */
 	                  y=, /* The name of the outcome variable */
                      ComputationOption=, /* 1, 2, or 3 for different options to be 
					                         sent to OrdinalTVEM. */
                      UsePropOdds=1, /* Either 0 or 1.  If 0 is specified, a linear
				                     TVEM will be fit.  If 1 is specified, a 
				                     proportional odds TVEM will be fit. */
                      UseRandom=0  /* 0, 1 or 2.  If 1 or 2 is specified, a random
				                    intercept will be included.  If 2 is specitife,
                                    a random slope will also be included.*/ 
						);
    /***************************************************************
    * OrdinalTvemLoop macro Version 2.0
    * By John DZIAK, Runze LI, and Anne BUU
    * Repeatedly calls OrdinalTvem for model selection.
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
	**************************************************************/
     %LET MinDeg = 1;
     %LET MaxDeg = 3;
     %LET MinNumInteriorKnots = 0;
     %LET MaxNumInteriorKnots = 10;
     DATA FitStats; RUN;
     DATA FitStats; deg = .; RUN;
     %DO deg = %EVAL(&MinDeg) %TO %EVAL(&MaxDeg);
         %DO NumInteriorKnots = %EVAL(&MinNumInteriorKnots) %TO %EVAL(&MaxNumInteriorKnots);
             %OrdinalTvem(dataset=&dataset,
                          id=&id,
                          t=&t,
                          cov=&cov,
                          tcov=&tcov,
                          y=&y,
                          UsePropOdds=&UsePropOdds,
                          UseRandom=&UseRandom,
                          DoPlots=0,
                          deg=&deg,
                          NumInteriorKnots=&NumInteriorKnots,
                          ComputationOption=&ComputationOption);
             DATA TheseFitStats;
                 SET TheseFitStats;
                 deg=&deg;
                 NumInteriorKnots=&NumInteriorKnots;
             RUN;
             DATA FitStats;
                 SET FitStats TheseFitStats;
                 WHERE deg ^= .;
             RUN;
         %END;
     %END;
%MEND;
