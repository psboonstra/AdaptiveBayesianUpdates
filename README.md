# Adaptive Bayesian Updates

### Current Suggested Citation

Boonstra, Philip S. and Barbaro, Ryan P., "Incorporating Historical Models with Adaptive Bayesian Updates" (March 2018). The University of Michigan Department of Biostatistics Working Paper Series. Working Paper 124.
https://biostats.bepress.com/umichbiostat/paper124

## Executive Summary
The functions <samp>glm_nab</samp> and <samp>glm_sab</samp> contained in the file <samp>Functions.R</samp> represent the primary statistical contribution from this manuscript. With these functions, plus the mean and variance of the coefficients from a historical regression model and the usual ingredients for fitting the current model of interest, a user can fit a Bayesian logistic regression with the adaptive priors that are described in the manuscript.

## Further details

In more detail, there are ten files included in this repository (in addition to this README): one text file (ending in <samp>.txt</samp>), three <samp>R</samp> scripts (ending in  <samp>.R</samp>), and six STAN functions (ending in  <samp>.stan</samp>). The simulation studies reported in Boonstra and Barbaro were run using commit 5; however, all subsequent commits have not changed any methodology and should therefore yield the same results, subject to simulation error. This has not been thoroughly checked. 

### Text file
<samp>runABUSims.txt</samp> is the script for running the simulation study on a cluster using SLURM. The following script does this:

<code> sbatch runABUSims.txt </code>

### <samp>R</samp> files
<samp>runMe.R</samp> provides two means of interfacing with the code. On a local machine, the user may choose a specific scenario (as described in that script) and run the code locally on his/her machine. On a cluster running SLURM, the user can use this script to submit multiple jobs simultaneously. 

<samp>Functions.R</samp> provides all of the necessary functions to fit the methods described in the paper as well as to run the simulation study. 

<samp>GenParams.R</samp> constructs inputs for running the simulation study. As described in the descriptions of this script and <samp>runMe.R</samp>, these inputs can be overwritten by the user.

### STAN files
The STAN files are described below. Note that these currently all implement a logistic link, but changing to a non-logistic link (i.e. log, probit, etc.) will be relatively easy. Upon using these for the first time, <samp>R</samp> will need to compile these programs, creating an R data object (ending in <samp>.rds</samp>) in the current working directory. Recompilation of the STAN files are not necessary as long as they stay the same.

<samp>RegHS_Stable.stan</samp> implements the regularized horseshoe prior, as described in Boonstra and Barbaro, applied to a logistic regression. An <samp>R</samp> user calls this with <samp>glm_standard</samp> in <samp>Functions.R</samp>. 

<samp>NAB_Stable.stan</samp>, <samp>NAB_Dev.stan</samp> both implement the 'naive adaptive Bayesian' prior, as described in Boonstra and Barbaro, applied to a logistic regression. The '<samp>_Dev</samp>' modifier was initially used for testing development versions of the prior against the current stable version. For the results reported in Boonstra and Barbaro, the only difference between the two is in the hyperprior distribution on &eta; (eta): in the former it is distributed as Inv-Gamma(2.5, 2.5), and in the latter it is Inv-Gamma(25, 25). The 'stable' versions are reported in the main manuscript. An <samp>R</samp> user calls this with the function <samp>glm_nab</samp> in <samp>Functions.R</samp>. 

<samp>SAB_Stable.stan</samp>, <samp>SAB_Dev.stan</samp> are analogous versions of the 'sensible adaptive Bayesian' prior. An <samp>R</samp> user calls this with the function <samp>glm_sab</samp> in <samp>Functions.R</samp>. 

<samp>RegStudT.stan</samp> implements a regularized Student-t prior applied to a logistic regression. This is not considered in the simulation study but is used in the exemplar data analysis ('PedRESC2'). A Student-t prior is applied to each regression coefficient using a normal-inverse-gamma distribution, but the latent inverse-gamma scale has a smooth upperbound provided by the user, so as to constrain very large scale values. An <samp>R</samp> user calls this with <samp>glm_studt</samp> in <samp>Functions.R</samp>. 

Built into each <samp>glm_</samp>*** function is a check for divergent iterations. The function will re-run if any divergent transitions are detected, up to a user-specified number of times (<samp>ntries</samp>), and return the results with the fewest. By virtue of the way this check is constructed, the user will see the following warning when divergent transitions are encountered:

<code>
Warning message:
In glm_sab(stan_path = paste0(stan_file_path, sab_stan_filename),  :
  NAs introduced by coercion
</code>


### Note, 10-Jul-2018:

After updating to version 3.5.0, <samp>R</samp> occasionally throws the following 'error':

<code>Error in x$.self$finalize() : attempt to apply non-function</code>

Error is used in quotes because it does not interrupt any processes and does not seem to affect any results. Searching online, this has been asked about by others and seems to be related to garbage collection:

http://discourse.mc-stan.org/t/very-mysterious-debug-error-when-running-rstanarm-rstan-chains-error-in-x-self-finalize-attempt-to-apply-non-function/4746