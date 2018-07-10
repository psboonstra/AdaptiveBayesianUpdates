# Adaptive Bayesian Updates

There are fourteen files included in this repository: one text file (ending in .txt), three R scripts (ending in .R), five STAN functions (ending in .stan), and five R data objects (ending in .rds). 

## Text file
runABUSims.txt is the script for running the simulation study on a cluster using SLURM. The following script does this:

<code> sbatch runABUSims.txt </code>

## R scripts
runMe.R provides two means of interfacing with the code. On a local machine, the user may choose a specific scenario (as described in that script) and run the code locally on his/her machine. On a cluster running SLURM, the user can use this script to submit multiple jobs simultaneously. 

Functions.R provides all of the necessary functions to fit the methods described in the paper as well as to run the simulation study. The functions glm_standard, glm_nab, and glm_sab are the main statistical contributions from this manuscript.

GenParams.R constructs inputs for running the simulation study. As described in the descriptions of this script and runMe.R, these inputs can be overwritten by the user.

## STAN files
The STAN files are described below. Note that these currently all implement a logistic link, but changing to a non-logistic link (i.e. log, probit, etc.) will be relatively easy. 

RegHS_Stable.stan implements the regularized horseshoe prior, as described in Boonstra and Barbaro, applied to a logistic regression. An R user calls this with glm_standard in Functions.R. 

NAB_Stable.stan, NAB_Dev.stan both implement the 'naive adaptive Bayesian' prior, as described in Boonstra and Barbaro, applied to a logistic regression. The '_Dev' suffix was originally used for testing development versions of the prior against the current stable version. For the results reported in Boonstra and Barbaro, the only difference between the two is in the hyperprior distribution on eta: in the former it is IG(2.5, 2.5), and in the latter it is IG(25, 25). The 'stable' versions are reported in the main manuscript. An R user calls this with glm_nab in Functions.R. 

SAB_Stable.stan, SAB_Dev.stan. Analogous versions of the 'sensible adaptive Bayesian' prior. 

## R data objects
These are the compiled versions of the programs described in the STAN file. If R does not find these, it will go ahead and recompile the STAN files and create new versions. 


Note, 10-Jul-2018:

After updating to version 3.5.0, R occasionally throws the following 'error':

<code>Error in x$.self$finalize() : attempt to apply non-function</code>

Error is used in quotes because it does not interrupt any processes and does not seem to affect any results. Searching online, this has been asked about by others and seems to be related to garbage collection:

http://discourse.mc-stan.org/t/very-mysterious-debug-error-when-running-rstanarm-rstan-chains-error-in-x-self-finalize-attempt-to-apply-non-function/4746