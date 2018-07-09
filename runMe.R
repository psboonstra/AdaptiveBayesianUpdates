#Expects that the following files are in your current working directory:
#NullPrior_Stable.stan, 
#SAB_Stable.stan, 
#NAB_Stable.stan,
#GenParams.R
#Functions.R

rm(list=ls(all=TRUE));
library(readstata13);library(mice);library(Hmisc);library(MASS);
library(rstan);library(Matrix);library(mnormt);
library(tidyr);library(plyr);library(dplyr);require(magrittr);
#Flag for whether this is running on a local machine or on a cluster running SLURM
my_computer = F;

if(my_computer) {
  array_id = 1;#Choose from between 1-120 if GenParams.R is used as-is
} else {
  array_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'));
}

#Recommended options from rstan:
rstan_options(auto_write = TRUE);
options(mc.cores = parallel::detectCores());

#If the files ending in .stan are elsewhere, indicate as such with this 
stan_file_path = "";
#Add this number to the saved object. Useful if a new batch of simulations is run and you want to continue the labeling scheme
array_id_offset = 0;
#Should the sim_id labels be randomly permuted across array_ids? Helpful if you want to look at intermediate results 
#and have good representation of all scenarios. Obviously, once all simulations are finished running, the final results will 
#remain the same
permute_sim_ids = F;
#Number of independent replicates to do for each scenario, e.g. 8 replicates times 10 simulated datasets per replicate = 80 
#simulated datasets per scenario
runs_per_scenario = 8;

#Without modification, running 'GenParams.R' produces a list of 120 lists, with each of the 120 lists corresponding to a unique
#scenario, defined as a combination from 6 sample sizes, 10 sets of true regression coefficients, and 2 distributions of predictors. 
source("GenParams.R");
rm(list=setdiff(ls(all=T),c("arglist","array_id","my_computer","permute_sim_ids","runs_per_scenario","array_id_offset")));
source("Functions.R");

if(permute_sim_ids) {
  set.seed(2);
  permute_sim_ids = sample(runs_per_scenario * length(arglist));
} else {
  permute_sim_ids = 1:(runs_per_scenario * length(arglist));
}

curr_args = arglist[[ceiling(permute_sim_ids[array_id]/runs_per_scenario)]];
curr_args[["random_seed"]] = curr_args[["random_seed"]] + array_id;

assign(paste0("sim",array_id_offset + array_id),do.call("run.sim",args = curr_args));
do.call("save",list(paste0("sim",array_id_offset + array_id),file=paste0("Sim",array_id_offset + array_id,".RData")));
