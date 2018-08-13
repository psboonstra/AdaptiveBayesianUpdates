# Adaptive Bayesian Updates

### Current Suggested Citation

Boonstra, Philip S. and Barbaro, Ryan P., "Incorporating Historical Models with Adaptive Bayesian Updates" (March 2018). The University of Michigan Department of Biostatistics Working Paper Series. Working Paper 124.
https://biostats.bepress.com/umichbiostat/paper124

## Executive Summary
The functions <samp>glm_nab</samp> and <samp>glm_sab</samp> contained in the file <samp>Functions.R</samp> represent the primary statistical contribution from this manuscript. With these functions, plus the mean and variance of the coefficients from a historical regression model and the usual ingredients for fitting the current model of interest, a user can fit a Bayesian logistic regression with the adaptive priors that are described in the manuscript.

## Further details

In more detail, there are twelve files included in this repository (in addition to this README): one text file (ending in <samp>.txt</samp>), five <samp>R</samp> scripts (ending in  <samp>.R</samp>), and six STAN functions (ending in <samp>.stan</samp>). The simulation studies reported in Boonstra and Barbaro were run using commit 22.

### Text file
<samp>runABUSims.txt</samp> is the script for submitting parallel runs of <samp>runABUSims.R</samp> (described below) to a cluster that is running SLURM. The following command does this:

<code> sbatch runABUSims.txt </code>

### <samp>R</samp> files

<samp>Functions.R</samp> provides all of the necessary functions to fit the methods described in the paper. 

<samp>Exemplar.R</samp> creates a single simulated dataset and walks through how to fit the methods described in the manuscript. 

<samp>GenParams.R</samp> constructs inputs for running the simulation study. As described in the script's documentation and the language below, these inputs can be overwritten by the user.

<samp>runABUSims.R</samp> is the script to conduct the large-scale simulation study described in the manuscript. On a local machine, the user may choose a specific <samp>array_id</samp> (as described in this script's documentation) and run the code locally on his/her machine. On a cluster running SLURM, the user can use this script to submit multiple jobs simultaneously (as described above). 

<samp>makeFigures.R</samp> gives the code to create the figures and tables in the manuscript and supplementary material reporting on the simulation study. 


### STAN files
The STAN files are described below. Note that these currently all implement a logistic link, but changing to a non-logistic link (i.e. log, probit, etc.) would be relatively easy. Upon using these for the first time, <samp>R</samp> will need to compile these programs, creating an <samp>R</samp> data object file (ending in <samp>.rds</samp>) in the current working directory. Re-compilation of the STAN files are not necessary as long as they are unchanged.

<samp>RegHS_stable.stan</samp> implements the regularized horseshoe prior, using the settings described described in Boonstra and Barbaro, applied to a logistic regression. An <samp>R</samp> user calls this with <samp>glm_standard</samp> in <samp>Functions.R</samp>. 

<samp>NAB_stable.stan</samp>, <samp>NAB_dev.stan</samp> both implement the 'naive adaptive Bayesian' prior, as described in Boonstra and Barbaro, applied to a logistic regression. The '<samp>_dev</samp>' modifier was initially used for testing development versions of the prior against the current stable version. For the results reported in Boonstra and Barbaro, the only difference between the two is in the hyperprior distribution on &eta; (eta): in the former it is distributed as Inv-Gamma(2.5, 2.5), and in the latter it is Inv-Gamma(25, 25). The 'stable' versions are reported in the main manuscript. An <samp>R</samp> user calls this with the function <samp>glm_nab</samp> in <samp>Functions.R</samp>. 

<samp>SAB_stable.stan</samp>, <samp>SAB_dev.stan</samp> are analogous versions of the 'sensible adaptive Bayesian' prior. An <samp>R</samp> user calls this with the function <samp>glm_sab</samp> in <samp>Functions.R</samp>. 

<samp>RegStudT.stan</samp> implements a regularized Student-t prior applied to a logistic regression. This is not considered in the simulation study but is used in the data analysis ('PedRESC2'). A Student-t prior is applied to each regression coefficient using a normal-inverse-gamma distribution, but the latent inverse-gamma scale has a smooth upperbound provided by the user, so as to constrain very large scale values. An <samp>R</samp> user calls this with <samp>glm_studt</samp> in <samp>Functions.R</samp>. 

### Divergent transitions

Built into each <samp>glm_</samp>*** function is a check for divergent iterations (http://mc-stan.org/users/documentation/case-studies/divergences_and_bias.html), which are simulatenously very helpful, very mysterious, and very frustrating. The function will re-run if any divergent transitions are detected, up to a user-specified number of times (<samp>ntries</samp>), and return the results in which the fewest divergent transitions were encountered. By virtue of the way this check is constructed, the user will see the following warning each time divergent transitions are encountered:

<code>
Warning message:
In glm_sab(stan_path = paste0(stan_file_path, sab_stan_filename),  :
  NAs introduced by coercion
</code>

### Compilation warning

STAN is smart enough to recognize the need for the normalizing constant and so, upon compilation, will give the following warning:

<code>
DIAGNOSTIC(S) FROM PARSER:
Warning (non-fatal):
Left-hand side of sampling statement (~) may contain a non-linear transform of a parameter or local variable.
If it does, you need to include a target += statement with the log absolute determinant of the Jacobian of the transform.
Left-hand-side of sampling statement:
    normalized_beta ~ normal(...)
</code>

This warning can be safely ignored because we do, in fact, calculate the normalizing constant. 

### Note, 10-Jul-2018:

After updating to version 3.5.0, <samp>R</samp> occasionally throws the following 'error':

<code>Error in x$.self$finalize() : attempt to apply non-function</code>

Error is used in quotes because it does not interrupt any processes and does not seem to affect any results. Searching online, this has been asked about by others and seems to be related to garbage collection:

http://discourse.mc-stan.org/t/very-mysterious-debug-error-when-running-rstanarm-rstan-chains-error-in-x-self-finalize-attempt-to-apply-non-function/4746