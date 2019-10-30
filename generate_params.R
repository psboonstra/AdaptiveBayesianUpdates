#DESCRIPTION: creates a list of lists called 'arglist'. The user should not 
# generally change anything in this script. If it is desired to overwrite
# any of the default values, create that variable in 'run_abu_sims.R' 
# (before sourcing this script on line 74 of that script). The default value 
# below will then be ignored. 

if(!"array_id_offset"%in%ls()){array_id_offset = 0;}
if(!"mc_iter_after_warmup"%in%ls()){mc_iter_after_warmup = 1.5e3;}
if(!"mc_chains"%in%ls()){mc_chains = 2;}
if(!"ntries_per_iter"%in%ls()){ntries_per_iter = 2;}
if(!"mc_warmup"%in%ls()){mc_warmup = 2.5e3;}
if(!"fit_methods"%in%ls()){fit_methods = T;}
if(!"dynamic_run"%in%ls()){dynamic_run = F;}
if(!"local_dof"%in%ls()){local_dof = 1;}
if(!"global_dof"%in%ls()){global_dof = 1;}
if(!"slab_precision"%in%ls()){slab_precision = (1.0/15.0)^2;}
if(!"nab_augmented_scale"%in%ls()){nab_augmented_scale = 0.05;}
if(!"power_prop_nonzero_prior"%in%ls()){power_prop_nonzero_prior = 1/3;}
if(!"phi_params"%in%ls()){
  phi_params = list("Agnostic" = c(mean = 0.50, sd = 2.5),
                    "Optimist" = c(mean = 1, sd = 0.25)
  );
}
if(!"sab_imputes_list"%in%ls()){
  sab_imputes_list = list(c(1,100),
                          c(1,1));
}

if(!"standard_stan_template"%in%ls()){
  standard_stan_template = adaptBayes:::stanmodels$RegHS_Stable;
}
if(!"sab_stan_template"%in%ls()){
  sab_stan_template = adaptBayes:::stanmodels$SAB_Stable;
}
if(!"sab_dev_stan_template"%in%ls()){
  sab_dev_stan_template = adaptBayes:::stanmodels$SAB_Dev;
}
if(!"nab_stan_template"%in%ls()){
  nab_stan_template = adaptBayes:::stanmodels$NAB_Stable;
}
if(!"nab_dev_stan_template"%in%ls()){
  nab_dev_stan_template = adaptBayes:::stanmodels$NAB_Dev;
}

if(!"n_list"%in%ls()){
  n_list = list(n_hist = c(100,100,400,400,1600,1600),#n_hist
                n_curr = c(100,200,100,200,100,200));#n_curr
}

if(!"betas_list"%in%ls()) {
  canonical_beta1 = c(1,1,1,1,1,1)/2;
  canonical_beta2 = c(1,0.5,0,0,0.5,1);
  canonical_beta3 = c(1,-0.5,0,0,-0.5,-1);
  canonical_beta4 = c(0.5,0.5,0,0,1,1);
  canonical_beta5 = c(0.5,0.5,0,0,-1,-1);
  canonical_beta6 = c(rep(0.5,4),rep(0.25,7), c(2, 1, 1, rep(0, 8)) );
  canonical_beta7 = rep(0.20,length=25);
  canonical_beta8 = c(1,1,1,0,0,rep(0.5,10),rep(0,35));
  canonical_beta9 = c(1,1,1,0,0,rep(0.25,20),rep(0,75));
  
  betas_list = list(true_mu_hist = rep(-1.0,10),
                    true_mu_curr = rep(-2.0,10),
                    true_betas_orig = list(canonical_beta1[1:4],
                                           #
                                           canonical_beta2[1:4],
                                           #
                                           canonical_beta3[1:4],
                                           #
                                           canonical_beta4[1:4],
                                           #
                                           canonical_beta5[1:4],
                                           #
                                           canonical_beta6[1:11],
                                           #
                                           canonical_beta7[1:5],
                                           canonical_beta7[1:20],
                                           #
                                           canonical_beta8[1:5],
                                           #
                                           canonical_beta9[1:5]
                    ),
                    true_betas_aug = list(canonical_beta1[5:6],
                                          #
                                          canonical_beta2[5:6],
                                          #
                                          canonical_beta3[5:6],
                                          #
                                          canonical_beta4[5:6],
                                          #
                                          canonical_beta5[5:6],
                                          #
                                          canonical_beta6[12:22],
                                          #
                                          canonical_beta7[6:25],
                                          canonical_beta7[21:25],
                                          #
                                          canonical_beta8[6:50],
                                          #
                                          canonical_beta9[6:100]
                                          #
                    )
  );
}
stopifnot(length(betas_list$true_mu_hist)==length(betas_list$true_mu_curr));
stopifnot(length(betas_list$true_mu_hist)==length(betas_list$true_betas_orig));
stopifnot(length(betas_list$true_mu_hist)==length(betas_list$true_betas_aug));
stopifnot(length(betas_list$true_mu_curr)==length(betas_list$true_betas_orig));
stopifnot(length(betas_list$true_mu_curr)==length(betas_list$true_betas_aug));
stopifnot(length(betas_list$true_betas_orig)==length(betas_list$true_betas_aug));


if(!"covariate_args_list"%in%ls()) {
  covariate_args_list = list(
    list(x_correlation = 0.20, x_orig_binom = "first_half", x_aug_binom = "first_half")
  )
}

if(!"n_sim"%in%ls()) {
  n_sim = 2;
}

all_varying = expand.grid(covariates = 1:length(covariate_args_list),n = 1:length(n_list[[1]]), betas = 1:length(betas_list[[1]]));
all_varying = cbind(all_varying, scenario = array_id_offset + (1:nrow(all_varying)));

random_seeds = sample(.Machine$integer.max - 1e4,nrow(all_varying));

arglist = list();
for(i in 1:nrow(all_varying)) {
  
  true_mu_hist = betas_list$true_mu_hist[all_varying[i,"betas"]];
  true_mu_curr = betas_list$true_mu_curr[all_varying[i,"betas"]];
  true_betas_orig = betas_list$true_betas_orig[[all_varying[i,"betas"]]];
  true_betas_aug = betas_list$true_betas_aug[[all_varying[i,"betas"]]];
  #
  
  covariate_args = covariate_args_list[[all_varying[i,"covariates"]]];
  if(covariate_args$x_orig_binom == "first_half") {
    covariate_args$x_orig_binom = 1:ceiling(length(true_betas_orig)/2);
  } else if(!is.na(covariate_args$x_orig_binom) && !all(covariate_args$x_orig_binom%in% (1:length(true_betas_orig)))) {
    stop("elements of 'x_orig_binom', if not NA, must be a subset of '1:length(true_betas_orig)'")
  }
  if( covariate_args$x_aug_binom == "first_half") {
    covariate_args$x_aug_binom = 1:ceiling(length(true_betas_aug)/2);
  } else if(!is.na(covariate_args$x_orig_aug) && !all(covariate_args$x_aug_binom%in% (1:length(true_betas_aug)))) {
    stop("elements of 'x_aug_binom', if not NA, must be a subset of '1:length(true_betas_aug)'")
  }
  #
  n_hist = n_list$n_hist[all_varying[i,"n"]];
  n_curr = n_list$n_curr[all_varying[i,"n"]];
  
  #
  assign(paste("sim",array_id + array_id_offset,"_params",sep=""),
         list(sim_number = all_varying[i,"scenario"],
              array_id = array_id,
              n_sim = n_sim,
              n_hist = n_hist,
              n_curr = n_curr,
              n_new = 1e3,
              true_mu_hist = true_mu_hist,
              true_mu_curr = true_mu_curr,
              true_betas_orig = true_betas_orig,
              true_betas_aug = true_betas_aug,
              beta_label = all_varying[i,"betas"],
              covariate_args = covariate_args,
              covariate_label = all_varying[i,"covariates"],
              local_dof = local_dof,
              global_dof = global_dof,
              slab_precision = slab_precision,
              nab_augmented_scale = nab_augmented_scale,
              power_prop_nonzero_prior = power_prop_nonzero_prior,
              sab_imputes_list = sab_imputes_list,
              standard_stan_template = standard_stan_template,
              sab_stan_template = sab_stan_template,
              sab_dev_stan_template = sab_dev_stan_template,
              nab_stan_template = nab_stan_template,
              nab_dev_stan_template = nab_dev_stan_template,
              phi_params = phi_params,
              mc_warmup = mc_warmup, 
              mc_iter_after_warmup = mc_iter_after_warmup, 
              mc_chains = mc_chains, 
              mc_thin = 1,
              mc_stepsize = 0.01,
              mc_adapt_delta_relaxed = 0.999,
              mc_adapt_delta_strict = 0.9999999,
              mc_max_treedepth = 20,
              ntries_per_iter = ntries_per_iter,
              random_seed = random_seeds[i],
              fit_marginal = F,
              fit_methods = fit_methods, 
              dynamic_run = dynamic_run))
  arglist = c(arglist,list(get(paste("sim",array_id + array_id_offset,"_params",sep=""))))
}
