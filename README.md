# A-statistical-method-to-model-non-stationarity-in-precipitation-records-changes
This repository contains the necessary code for reproducing the results presented in "A statistical method to model non stationarity in precipitation records changes. We find here 4 files:

1. cmip6_ipsl.zip
It contains the necessary data for reproducing the analysis presentend in the article. This compressed folder has two R object files.
  i) df_pr_cmip6_yearmax_ipsl.rds:  the yearly maxima of daily precipitation from IPSL-CM6A-LR with a spatial resolution of 72 x 36 grid points.
  ii) df_lonlat.rds: latitude and longitude of each gridpoint
  
2. p1rt_functions_for_the study.R
R code with the necessary functions for the estimation of record probability p_{1,r}(t) and its confidence intervals.    

3. p1rt_study_from_the_article.R
R code that reproduces the analysis from the article (Section 3) and produces Figure 3 and Figure 4. The code calls functions from p1rt_functions_for_the_study.R and uses as input data from cmip6_ipsl.zip

4. p1rvarcmip6.zip
 R objects containing the Estimators \hat{p}_{1,r}(t) and variances var(\hat{p}_{1,r}(t)) used in Section 3 of the article. For each gridpoint i, we have two files
  i) p1r_i_cmip6_r1i1p1f1.rds
  ii) var_i_cmip6_r1i1p1f1.rds

One can use the code in p1rt_study_from_the_article.R to estimate \hat{p}_{1,r}(t) and variances var(\hat{p}_{1,r}(t)) or just call objects 4.i) and 4.ii) and direclty plot Figure 3 and Figure 4. The code from p1rt_functions_for_the study.R can be used to for other studies.

