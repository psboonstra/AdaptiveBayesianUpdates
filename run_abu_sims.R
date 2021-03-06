# DESCRIPTION: Running this script without modification produces an object 
# called 'arglist', which is a list of 60 lists, with each of the 60 lists 
# corresponding to a unique scenario, defined as a combination from  sample 
# sizes and 10 sets of true regression coefficients. Then, one of these 
# lists is extracted based upon the value of array_id, and the simulator 
# function called simulator is called on this scenario. Finally, the result 
# is saved as an binary data file (an R workspace) 
# called paste0("Sim",array_id_offset + array_id,".RData")

if(!require(adaptBayes)) {
  library(devtools)
  # may take some time:
  install_github('umich-biostatistics/adaptBayes') 
} 

library(mice);library(Hmisc);library(MASS);
library(rstan);library(Matrix);library(mnormt);
library(tidyverse);

#Flag for whether this is running on a local machine or on a cluster running SLURM
my_computer = T;

# The simulation study in Boonstra and Barbaro was conducted in two batches of 960 
# jobs followed by another 960 (because the cluster only accepts jobs in batches 
# of size up to 1000). 'which_run' indicates whether this is the first batch (= 1 )
# or the second batch ( = 2)
which_run = 1;

if(my_computer) {
  #Choose from between 1-960 if GenParams.R is used as-is
  array_id = 1;
} else {
  array_id = as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'));
}

#Recommended options from rstan:
rstan_options(auto_write = TRUE);
options(mc.cores = parallel::detectCores());


if(which_run == 1) {#first batch
  #'jobs_per_scenario' is the number of parallel independent jobs to send for 
  # each scenario, and 'n_sim' is the number of independent datasets per job
  # to generate. In Boonstra and Barbaro, initially 16 jobs per scenario times 
  # 2 simulated datasets per job = 32 simulated datasets per scenario were run. 
  # Then, an additional 16 jobs times 6 simulated datasets per job = 96 
  # simulated datasets per scenario were run, yielding a total of 100 simulated 
  # datasets per scenario presented in the results. 
  jobs_per_scenario = 16;
  n_sim = 2;#Worst case running time needed is ~23 hours
  #'array_id_offset' is added to the label of the saved object. Useful when a 
  # new batch of jobs is run and you want to continue the labeling scheme. 
  # In Boonstra and Barbaro, 960 jobs were initially submitted (60 scenarios 
  # times 16 jobs per scenario), followed by 960 subsequent jobs (60 scenarios
  # times 16 jobs per scenario). Thus, for the second batch, array_id_offset 
  # was set to be 960, so that the first job from the second batch would be 
  # labeled 961. 
  array_id_offset = 0;
} else {#second batch
  jobs_per_scenario = 16;
  n_sim = 6;#Worst case running time needed is ~75 hours
  array_id_offset = 960;  
}
# Should the sim_id labels be randomly permuted across array_ids? Helpful if 
# you're impatient and intent on looking at intermediate results along the 
# way; you'll have good representation across all scenarios. Obviously, once 
# all jobs are finished running, the final results will 
# remain the same, and this simply results in a permutation of the labels. 
permute_array_ids = T;

# Before calling the next line, the user should specify any parameters that 
# she/he wishes to change from their default values. Any specified values will
# overwrite their default values, which are set by sourcing GenParams.R. 
source("generate_params.R");
source("functions_simulation.R");

# This is the step that permutes the job ids. If FALSE, then all jobs from the 
# same scenario will occur in contiguous blocks. 
if(permute_array_ids) {
  set.seed(2);
  permute_array_ids = sample(jobs_per_scenario * length(arglist));
} else {
  permute_array_ids = 1:(jobs_per_scenario * length(arglist));
}

curr_args = arglist[[ceiling(permute_array_ids[array_id]/jobs_per_scenario)]];
curr_args[["random_seed"]] = curr_args[["random_seed"]] + array_id_offset + array_id;

assign(paste0("sim",array_id_offset + array_id),do.call("simulator",args = curr_args));
do.call("save",list(paste0("sim",array_id_offset + array_id),file=paste0("Sim",array_id_offset + array_id,".RData")));
