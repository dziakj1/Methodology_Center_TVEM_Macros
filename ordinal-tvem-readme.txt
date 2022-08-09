# OrdinalTvem & OrdinalTvemLoop macros (Version 2.0)

## Authors
John DZIAK, Runze LI, and Anne BUU

## Description 

The OrdinalTvem macro fits a TVEM model to (potentially longitudinal) ordinal data. The OdinalTvemLoop macro repeatedly calls OrdinalTvem for model selection. The underlying model for OrdinalTvem and OrdinalTvemLoop is further explained in the following publication:

> Dziak, J. J., Li, R., Zimmerman, M. A., & Buu, A. (2014). Time-varying effect models for ordinal responses with applications in substance abuse research. *Statistics in Medicine*, 33: 5126-5137. doi: 10.1002/sim.6303

The full text is available here: [https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4227951/](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4227951/)

## Usage 

These macros have been tested on SAS version 9.2. Descriptions of macro parameters are included in the source files. 

## File Manifest / Data Structure
    
- `OrdinalTvem.sas`: The OrdinalTvem macro
- `OrdinalTvemLoop.sas` The OrdinalTvemLoop macro
    
## Acknowledgements & References

`OrdinalTvem` uses code which was written by Xianming TAN to construct the spline bases and which is also found in the MixTVEM macro. We fit a varying-coefficient (time-varying effect) model for ordinal data, using either a linear model or a proportional odds logistic model, and an unpenalized B-spline approach. See:

- Eilers, P.H.C. and Marx, B.D. (1996). Flexible smoothing using B-splines and penalized likelihood. *Statistical Science* 11(2): 89-121.
- Hastie, T., & Tibshirani, R. (1993). Varying-coefficient models. *Journal of the Royal Statistical Society, Series B*, 55, 757-796.
- Shiyko, M. P., Lanza, S. T., Tan, X., Li, R., Shiffman, S. (2012). Using the Time-Varying Effect Model (TVEM) to Examine Dynamic Associations between Negative Affect and Self Confidence on Smoking Urges: Differences between Successful Quitters and Relapsers. *Prevention Science*, 13, 288-299.
- Ramsay, J., Hooker, G., & Graves, S. (2009). Functional Data Analysis with R and MATLAB. New York: Springer.
- SAS Institute Inc. (2011). SAS/STAT (c) 9.3 user's guide: The GLIMMIX procedure (chapter)). Cary, NC: SAS Institute Inc.
- Tan, X., Shiyko, M. P., Li, R., Li, Y., & Dierker, L. (2011, November 21). A time-varying effect model for intensive longitudinal data. *Psychological Methods*. Advance online publication. doi: 10.1037/a0025814.

## Copyright & License 

Copyright (c) 2014 The Pennsylvania State University

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License as
published by the Free Software Foundation; either version 2 of
the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.