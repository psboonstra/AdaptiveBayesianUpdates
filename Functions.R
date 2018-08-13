
#DESCRIPTION: self-explanatory
expit = function(x) { 1/(1+exp(-x));}
logit = function(x) { log(x/(1-x));}

#DESCRIPTION: Error-handling function
#Borrowed from the R package simsalapar 
#https://www.rdocumentation.org/packages/simsalapar/versions/1.0-9/topics/tryCatch.W.E
tryCatch.W.E <- function(expr)
{
  W <- NULL
  w.handler <- function(w){ # warning handler
    W <<- c(W,w)
    invokeRestart("muffleWarning")
  }
  list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                   warning = w.handler),
       warning = W)
}

#DESCRIPTION: Simulator function for drawing binary outcomes (historical, current, and new [for validation]), 
#the probabilities of which are logistic-linear functions of normal and/or bernoulli  distributed 
#predictors. 
#
#
#ARGUMENTS:
#
#n_hist (pos. integer) size of historical data; n_hist in the paper
#
#n_curr (pos. integer) size of current data; n_curr in the paper
#
#n_new (pos. integer) size of testing data (for prediction)
#
#true_mu_hist (real) true intercept for generating historical model. mu_hist in Boonstra and Barbaro 
#
#true_mu_curr (real) true intercept for generating current / new model. mu in Boonstra and Barbaro 
#
#true_betas_orig (vector) true regression coefficients corresponding to original covariates. Beta^o in Boonstra and Barbaro
#
#true_betas_aug (vector) true regression coefficients corresponding to augmented covariates. Beta^a in Boonstra and Barbaro
#
#covariate_args (list) the named arguments are simple ways to govern the distribution of the predictors. 
#x_correlation is the normal correlation between any pair of predictors; x_orig_binom is the 
#integer indices (any subset of 1, ..., length(true_betas_orig)) indicating which of the original 
#covariates should be transformed to binary values based upon being less than or greater than zero; 
#x_aug_binom is the analogous set of indices for the augmented predictors (any subset of 1, ..., length(true_betas_aug))

draw_data = function(n_hist = 150,
                     n_curr = 50,
                     n_new = 100,
                     true_mu_hist = 0,
                     true_mu_curr = 0,
                     true_betas_orig = 0,
                     true_betas_aug = 0,
                     covariate_args = list(x_correlation = 0, 
                                           x_orig_binom = NULL, 
                                           x_aug_binom = NULL)
) {
  
  p = length(true_betas_orig);
  q = length(true_betas_aug);
  stopifnot(length(true_mu_hist) == 1 && length(true_mu_curr) == 1);
  
  x_all = matrix(rnorm((n_hist + n_curr + n_new)*(p+q)),nrow=n_hist + n_curr + n_new)%*%chol(diag(1 - covariate_args$x_correlation,p+q) + covariate_args$x_correlation);#original covariates are N(0,1)
  x_all_orig = x_all[,1:p,drop = F];
  #Binary covariates will be -1 or 1 with equal prevalence (assuming the latent normal is mean zero), 
  #which will result in a random variable with mean zero and variance 1 
  if(length(covariate_args$x_orig_binom)) {
    x_all_orig[,covariate_args$x_orig_binom] = 2 * (x_all_orig[,covariate_args$x_orig_binom,drop = F] > 0) - 1;
  }
  x_all_aug = x_all[,(p + 1):(p + q),drop = F];
  if(length(covariate_args$x_aug_binom)) {
    x_all_aug[,covariate_args$x_aug_binom] = 2 * (x_all_aug[,covariate_args$x_aug_binom,drop = F] > 0) - 1;
  }
  
  #Linear predictors (two for each observation: contribution from original covariates and augmented covariates)
  lin_pred_x_orig = drop(x_all_orig%*%true_betas_orig);
  lin_pred_x_aug = drop(x_all_aug%*%true_betas_aug);
  #Note the difference in intercepts between historical data and current data
  risk_all = 1/(1+exp(-c(rep(true_mu_hist,n_hist),rep(true_mu_curr,n_curr+n_new)) - lin_pred_x_orig - lin_pred_x_aug));
  
  y_all =  rbinom(n_hist + n_curr + n_new, 1, risk_all);
  
  x_hist_orig = x_all_orig[1:n_hist,,drop=F];
  x_curr_orig = x_all_orig[(n_hist + 1):(n_hist + n_curr),,drop=F];
  x_new_orig = x_all_orig[(n_hist + n_curr + 1):(n_hist + n_curr + n_new),,drop=F];
  #
  x_hist_aug = x_all_aug[1:n_hist,,drop=F];
  x_curr_aug = x_all_aug[(n_hist + 1):(n_hist + n_curr),,drop=F];
  x_new_aug = x_all_aug[(n_hist + n_curr + 1):(n_hist + n_curr + n_new),,drop=F];
  #
  y_hist = y_all[1:n_hist];
  y_curr = y_all[(n_hist+1):(n_hist + n_curr)];
  y_new = y_all[(n_hist + n_curr + 1):(n_hist + n_curr + n_new)];
  
  risk_new = risk_all[(n_hist + n_curr + 1):(n_hist + n_curr + n_new)];
  
  list(x_hist_orig = x_hist_orig, 
       x_hist_aug = x_hist_aug, 
       y_hist = y_hist, 
       x_curr_orig = x_curr_orig, 
       x_curr_aug = x_curr_aug, 
       y_curr = y_curr,
       x_new_orig = x_new_orig, 
       x_new_aug = x_new_aug, 
       y_new = y_new,
       lin_pred_x_orig = lin_pred_x_orig,
       lin_pred_x_aug = lin_pred_x_aug,
       risk_new = risk_new)
  
}

#DESCRIPTION: Program for fitting a GLM equipped with the 'standard' prior evaluated 
#in Boonstra and Barbaro, which is the regularized horseshoe. It has two intended uses: 
#compile stan scripts or analyze data. First, if the user provides nothing but a valid 
#'stan_path', then the stan script is compiled. Second, the user provides both a compiled 
#stanfit object as well asvalues for y, x_standardized, #, q, and any other desired 
#arguments to actually fit a regression. 
#
#
#ARGUMENTS:
#
#stan_fit: an R object of class stanfit, which allows the function to run without recompiling the stan code.
#
#stan_path: (character) a path pointing to a .stan file, which indicates the stan code to compile and run. If
#both stan_fit and stan_path are provided, stan_fit takes precedence. 
#
#y (vector) outcomes corresponding to the type of glm desired. This should match whatever datatype is expected 
#by the stan program.
#
#x_standardized (matrix) matrix of numeric values with number of rows equal to the length of y and number of columns
#equal to p+q. It is assumed without verification that each column is standardized to whatever scale the prior 
#expects - in Boonstra and Barbaro, all predictors are marginally generated to have mean zero and unit variance, so no 
#standardization is conducted. In practice, all data should be standardized to have a common scale before model fitting. 
#If regression coefficients on the natural scale are desired, they be easily obtained through unstandardizing. 
#
#p, q (nonneg. integers) numbers, the sum of which add up to the number of columns in x_standardized. For the standard
#prior, this distinction is only needed if a different constant scale parameter (beta_orig_scale, beta_aug_scale), which is 
#the constant 'c' in the notation of Boonstra and Barbaro, is used. 
#
#beta_orig_scale, beta_aug_scale (pos. real) constants indicating the prior scale of the horseshoe. Both values correspond 
#to 'c' in the notation of Boonstra and Barbaro, because that paper never considers beta_orig_scale!=beta_aug_scale
#
#local_dof, global_dof (pos. integer) numbers indicating the degrees of freedom for lambda_j and tau, respectively. Boonstra, 
#et al. never considered local_dof != 1 or global_dof != 1. 
#
#slab_precision (pos. real) the slab-part of the regularized horseshoe, this is equivalent to (1/d)^2 in the notation of
#Boonstra and Barbaro 
#
#intercept_offset (vector) vector of 0's and 1's equal having the same length as y. Those observations with a value of 1 
#have an additional constant offset in their linear predictor, effectively a different intercept. This is useful to jointly 
#regress two datasets in which it is believed that the regression coefficients are the same but not the intercepts and could be 
#useful (but was not used) in the simulation study to compare to a benchmark, namely if both the historical and current datasets 
#were available but there is a desire to adjust for potentially different baseline prevalences. 
#
#only_prior (logical) should all data be ignored, sampling only from the prior?
#
#ntries (pos. integer) the stan function will run up to this many times, stopping either when the number of 
#*divergent transitions* is zero or when ntries has been reached. The reported fit will be that with the fewest number of 
#divergent iterations. 

glm_standard = function(stan_fit = NA, 
                        stan_path,
                        y = c(0,1),
                        x_standardized = matrix(0,length(y),6), 
                        p = 3, 
                        q = 3,
                        beta_orig_scale = 1, 
                        beta_aug_scale = 1, 
                        local_dof = 1, 
                        global_dof = 1, 
                        slab_precision = (1/15)^2, 
                        intercept_offset = NULL,
                        only_prior = F, 
                        mc_warmup = 50, 
                        mc_iter_after_warmup = 50, 
                        mc_chains = 1, 
                        mc_thin = 1, 
                        mc_stepsize = 0.1, 
                        mc_adapt_delta = 0.9,
                        mc_max_treedepth = 15,
                        ntries = 1) {
  
  stopifnot(ncol(x_standardized) == (p+q));
  
  max_divergences = -Inf;
  accepted_divergences = Inf;
  curr_try = 1;
  if(is.null(intercept_offset)) {intercept_offset = numeric(length(y));}
  
  while(curr_try <= ntries) {
    assign("curr_fit",tryCatch.W.E(stan(file = stan_path,
                                        fit = stan_fit,
                                        data = list(n_stan = length(y),
                                                    p_stan = p,
                                                    q_stan = q,
                                                    y_stan = y,
                                                    x_standardized_stan = x_standardized,
                                                    local_dof_stan = local_dof,
                                                    global_dof_stan = global_dof,
                                                    beta_orig_scale_stan = beta_orig_scale,
                                                    beta_aug_scale_stan = beta_aug_scale,
                                                    slab_precision_stan = slab_precision,
                                                    intercept_offset_stan = intercept_offset,
                                                    only_prior = as.integer(only_prior)), 
                                        warmup = mc_warmup, 
                                        iter = mc_iter_after_warmup + mc_warmup, 
                                        chains = mc_chains, 
                                        thin = mc_thin,
                                        control = list(stepsize = mc_stepsize,
                                                       adapt_delta = mc_adapt_delta,
                                                       max_treedepth = mc_max_treedepth)))); 
    if("simpleError"%in%class(curr_fit$value) || "error"%in%class(curr_fit$value)) {
      stop(curr_fit$value);
    }
    if(!"stanfit"%in%class(stan_fit)) {
      break;
    }
    divergent_check = unlist(lapply(curr_fit$warning,grep,pattern="divergent transitions",value=T));
    rhat_check = max(summary(curr_fit$value)$summary[,"Rhat"],na.rm=T);
    #Originally, the break conditions were baesd upon having both no divergent transitions as well as a max Rhat (i.e. gelman-rubin 
    #diagnostic) sufficiently close to 1. I subsequently changed the conditions to be based only upon the first, which is reflected
    #by setting rhat = T immediately below. 
    break_conditions = c(divergence = F, rhat = T);
    if(length(divergent_check) == 0) {#corresponds to zero divergent transitions
      curr_divergences = 0;
      max_divergences = max(max_divergences,curr_divergences,na.rm=T);
      break_conditions["divergence"] = T;
    } else {#corresponds to > zero divergent transitions
      curr_divergences <- max(as.numeric(strsplit(divergent_check," ")$message),na.rm=T);
      max_divergences = max(max_divergences,curr_divergences,na.rm=T);
      curr_try = curr_try + 1;
    }
    #update if fewer divergent transitions were found
    if(curr_divergences < accepted_divergences) {
      accepted_divergences = curr_divergences;
      max_rhat = rhat_check;
      foo = rstan::extract(curr_fit$value);
      hist_beta0 = as.numeric(foo$mu);
      curr_beta0 = as.numeric(foo$mu) + as.numeric(foo$mu_offset);
      curr_beta = foo$beta;
      theta_orig = foo$theta_orig;
      theta_aug = foo$theta_aug;
    }
    if(all(break_conditions)) {
      break;
    }
  }
  if(!"stanfit"%in%class(stan_fit)) {
    curr_fit$value;
  } else {
    list(accepted_divergences = accepted_divergences,
         max_divergences = max_divergences,
         max_rhat = max_rhat,
         hist_beta0 = hist_beta0,
         curr_beta0 = curr_beta0,
         curr_beta = curr_beta,
         theta_orig = theta_orig,
         theta_aug = theta_aug);
  }
}

#DESCRIPTION: Program for fitting a GLM equipped with the 'naive adaptive bayes' prior evaluated in the manuscript
#
#
#ARGUMENTS: (only those distinct from glm_standard are discussed)
#
#alpha_prior_mean (vector) p-length vector giving the mean of alpha from the historical analysis, 
#corresponds to m_alpha in Boonstra and Barbaro 
#
#alpha_prior_cov (matrix) pxp positive definite matrix giving the variance of alpha from the historical
#analysis, corresponds to S_alpha in Boonstra and Barbaro 
#
#phi_mean (real) mean of phi corresponding to a normal distribution, support is truncated to [0,1]
#
#phi_sd (pos. real) sd of phi corresponding to a normal distribution, support is truncated to [0,1]
#
#beta_aug_scale_tilde (pos. real) constant indicating the prior scale of the horseshoe for the augmented
#covariates when phi = 1, i.e. when the historical analysis is fully used. This corresponds to tilde_c in
#Boonstra and Barbaro 

glm_nab = function(stan_fit = NA, 
                   stan_path,
                   y = c(0,1),
                   x_standardized = matrix(0,length(y),6), 
                   alpha_prior_mean = rep(0, 3),
                   alpha_prior_cov = diag(1, 3),
                   phi_mean = 0.5,
                   phi_sd = 2.5,
                   beta_orig_scale = 1, 
                   beta_aug_scale = 1, 
                   beta_aug_scale_tilde = 1,
                   local_dof = 1, 
                   global_dof = 1, 
                   slab_precision = (1/15)^2, 
                   only_prior = F, 
                   mc_warmup = 50, 
                   mc_iter_after_warmup = 50, 
                   mc_chains = 1, 
                   mc_thin = 1, 
                   mc_stepsize = 0.1, 
                   mc_adapt_delta = 0.9,
                   mc_max_treedepth = 15,
                   ntries = 1,
                   eigendecomp_hist_var = NULL,
                   scale_to_variance225 = NULL
) {
  
  if(is.null(eigendecomp_hist_var)) {
    eigendecomp_hist_var = eigen(alpha_prior_cov);
  }
  eigenvec_hist_var = t(eigendecomp_hist_var$vectors);
  sqrt_eigenval_hist_var = sqrt(eigendecomp_hist_var$values);
  
  if(is.null(scale_to_variance225)) {
    scale_to_variance225 = diag(alpha_prior_cov) / 225;
  }
  
  p = length(alpha_prior_mean);
  q = ncol(x_standardized) - p;
  if(p == 1) {
    alpha_prior_mean = array(alpha_prior_mean,dim=1);
    sqrt_eigenval_hist_var = array(sqrt_eigenval_hist_var,dim=1);
    scale_to_variance225 = array(scale_to_variance225,dim=1);
  }
  
  max_divergences = -Inf;
  accepted_divergences = Inf;
  curr_try = 1;
  
  while(curr_try <= ntries) {
    assign("curr_fit",tryCatch.W.E(stan(file = stan_path,
                                        fit = stan_fit,
                                        data = list(n_stan = length(y),
                                                    p_stan = p,
                                                    q_stan = q,
                                                    y_stan = y,
                                                    x_standardized_stan = x_standardized,
                                                    alpha_prior_mean_stan = alpha_prior_mean,
                                                    alpha_prior_cov_stan = alpha_prior_cov,
                                                    sqrt_eigenval_hist_var_stan = sqrt_eigenval_hist_var,
                                                    eigenvec_hist_var_stan = eigenvec_hist_var,
                                                    local_dof_stan = local_dof,
                                                    global_dof_stan = global_dof,
                                                    beta_orig_scale_stan = beta_orig_scale,
                                                    beta_aug_scale_stan = beta_aug_scale,
                                                    beta_aug_scale_tilde_stan = beta_aug_scale_tilde,
                                                    slab_precision_stan = slab_precision,
                                                    scale_to_variance225 = scale_to_variance225,
                                                    phi_mean_stan = phi_mean,
                                                    phi_sd_stan = phi_sd,
                                                    only_prior = as.integer(only_prior)), 
                                        warmup = mc_warmup, 
                                        iter = mc_iter_after_warmup + mc_warmup, 
                                        chains = mc_chains, 
                                        thin = mc_thin,
                                        control = list(stepsize = mc_stepsize,
                                                       adapt_delta = mc_adapt_delta,
                                                       max_treedepth = mc_max_treedepth)))); 
    
    if("simpleError"%in%class(curr_fit$value) || "error"%in%class(curr_fit$value)) {
      stop(curr_fit$value);
    }
    if(!"stanfit"%in%class(stan_fit)) {
      break;
    }
    divergent_check = unlist(lapply(curr_fit$warning,grep,pattern="divergent transitions",value=T));
    rhat_check = max(summary(curr_fit$value)$summary[,"Rhat"],na.rm=T);
    #Originally, the break conditions were baesd upon having both no divergent transitions as well as a max Rhat (i.e. gelman-rubin 
    #diagnostic) sufficiently close to 1. I subsequently changed the conditions to be based only upon the first, which is reflected
    #by setting rhat = T immediately below. 
    break_conditions = c(divergence = F, rhat = T);
    if(length(divergent_check) == 0) {#corresponds to zero divergent transitions
      curr_divergences = 0;
      max_divergences = max(max_divergences,curr_divergences,na.rm=T);
      break_conditions["divergence"] = T;
    } else {#corresponds to > zero divergent transitions
      curr_divergences <- max(as.numeric(strsplit(divergent_check," ")$message),na.rm=T);
      max_divergences = max(max_divergences,curr_divergences,na.rm=T);
      curr_try = curr_try + 1;
    }
    #update if fewer divergent transitions were found
    if(curr_divergences < accepted_divergences) {
      accepted_divergences = curr_divergences;
      max_rhat = rhat_check;
      foo = rstan::extract(curr_fit$value);
      curr_beta0 = as.numeric(foo$mu);
      curr_beta = foo$beta;
      theta_orig = foo$theta_orig;
      theta_aug = foo$theta_aug;
      phi = foo$phi_copy;
      eta = foo$eta;
    }
    if(all(break_conditions)) {
      break;
    }
  }
  if(!"stanfit"%in%class(stan_fit)) {
    curr_fit$value;
  } else {
    list(accepted_divergences = accepted_divergences,
         max_divergences = max_divergences,
         max_rhat = max_rhat,
         curr_beta0 = curr_beta0,
         curr_beta = curr_beta,
         theta_orig = theta_orig,
         theta_aug = theta_aug,
         phi = phi,
         eta = eta);
  }
}

#DESCRIPTION: Program for fitting a GLM equipped with the 'sensible adaptive bayes' prior evaluated in the manuscript
#
#
#ARGUMENTS: (only those distinct from glm_standard and glm_nab are discussed)
#
#aug_projection: (matrix) pxq matrix that approximately projects the regression coefficients of 
#the augmented predictors onto the space of the regression coefficients for the original predictors. 
#This is the matrix P in the notation of Boonstra and Barbaro. It can be calculated using the function 
#'create_projection'. 
#
#eigendecomp_hist_var: R object of class 'eigen' containing a pxp matrix of eigenvectors in each row 
#(equivalent to v_0 in Boonstra and Barbaro) and a p-length vector of eigenvalues. This is by default 
#equal to eigen(alpha_prior_cov)
#
#scale_to_variance225: a vector assumed to be such that, when multiplied by the diagonal elements of 
#alpha_prior_cov, the result is a vector of elements each equal to 225. This is explicitly calculated 
#if it is not provided

glm_sab = function(stan_fit = NA, 
                   stan_path,
                   y = c(0,1),
                   x_standardized = matrix(0,length(y),6), 
                   alpha_prior_mean = rep(0, 3),
                   alpha_prior_cov = diag(1, 3),
                   aug_projection = diag(1, 3),
                   phi_mean = 0.5,
                   phi_sd = 2.5,
                   beta_orig_scale = 1, 
                   beta_aug_scale = 1, 
                   local_dof = 1, 
                   global_dof = 1, 
                   slab_precision = (1/15)^2, 
                   only_prior = F, 
                   mc_warmup = 50, 
                   mc_iter_after_warmup = 50, 
                   mc_chains = 1, 
                   mc_thin = 1, 
                   mc_stepsize = 0.1, 
                   mc_adapt_delta = 0.9,
                   mc_max_treedepth = 15,
                   ntries = 1,
                   eigendecomp_hist_var = NULL,
                   scale_to_variance225 = NULL
) {
  
  if(is.null(eigendecomp_hist_var)) {
    eigendecomp_hist_var = eigen(alpha_prior_cov);
  }
  eigenvec_hist_var = t(eigendecomp_hist_var$vectors);
  sqrt_eigenval_hist_var = sqrt(eigendecomp_hist_var$values);
  
  if(is.null(scale_to_variance225)) {
    scale_to_variance225 = diag(alpha_prior_cov) / 225;
  }
  
  p = length(alpha_prior_mean);
  q = ncol(x_standardized) - p;
  if(p == 1) {
    alpha_prior_mean = array(alpha_prior_mean,dim=1);
    sqrt_eigenval_hist_var = array(sqrt_eigenval_hist_var,dim=1);
    scale_to_variance225 = array(scale_to_variance225,dim=1);
  }
  
  max_divergences = -Inf;
  accepted_divergences = Inf;
  curr_try = 1;
  
  while(curr_try <= ntries) {
    assign("curr_fit",tryCatch.W.E(stan(file = stan_path,
                                        fit = stan_fit,
                                        data = list(n_stan = length(y),
                                                    p_stan = p,
                                                    q_stan = q,
                                                    y_stan = y,
                                                    x_standardized_stan = x_standardized,
                                                    aug_projection_stan = aug_projection,
                                                    alpha_prior_mean_stan = alpha_prior_mean,
                                                    alpha_prior_cov_stan = alpha_prior_cov,
                                                    sqrt_eigenval_hist_var_stan = sqrt_eigenval_hist_var,
                                                    eigenvec_hist_var_stan = eigenvec_hist_var,
                                                    local_dof_stan = local_dof,
                                                    global_dof_stan = global_dof,
                                                    beta_orig_scale_stan = beta_orig_scale,
                                                    beta_aug_scale_stan = beta_aug_scale,
                                                    slab_precision_stan = slab_precision,
                                                    scale_to_variance225 = scale_to_variance225,
                                                    phi_mean_stan = phi_mean,
                                                    phi_sd_stan = phi_sd,
                                                    only_prior = as.integer(only_prior)), 
                                        warmup = mc_warmup, 
                                        iter = mc_iter_after_warmup + mc_warmup, 
                                        chains = mc_chains, 
                                        thin = mc_thin,
                                        control = list(stepsize = mc_stepsize,
                                                       adapt_delta = mc_adapt_delta,
                                                       max_treedepth = mc_max_treedepth)))); 
    
    if("simpleError"%in%class(curr_fit$value) || "error"%in%class(curr_fit$value)) {
      stop(curr_fit$value);
    }
    if(!"stanfit"%in%class(stan_fit)) {
      break;
    }
    divergent_check = unlist(lapply(curr_fit$warning,grep,pattern="divergent transitions",value=T));
    rhat_check = max(summary(curr_fit$value)$summary[,"Rhat"],na.rm=T);
    #Originally, the break conditions were baesd upon having both no divergent transitions as well as a max Rhat (i.e. gelman-rubin 
    #diagnostic) sufficiently close to 1. I subsequently changed the conditions to be based only upon the first, which is reflected
    #by setting rhat = T immediately below.
    break_conditions = c(divergence = F, rhat = T);
    if(length(divergent_check) == 0) {#corresponds to zero divergent transitions
      curr_divergences = 0;
      max_divergences = max(max_divergences,curr_divergences,na.rm=T);
      break_conditions["divergence"] = T;
    } else {#corresponds to > zero divergent transitions
      curr_divergences <- max(as.numeric(strsplit(divergent_check," ")$message),na.rm=T);
      max_divergences = max(max_divergences,curr_divergences,na.rm=T);
      curr_try = curr_try + 1;
    }
    #update if fewer divergent transitions were found
    if(curr_divergences < accepted_divergences) {
      accepted_divergences = curr_divergences;
      max_rhat = rhat_check;
      foo = rstan::extract(curr_fit$value);
      curr_beta0 = as.numeric(foo$mu);
      curr_beta = foo$beta;
      theta_orig = foo$theta_orig;
      theta_aug = foo$theta_aug;
      phi = foo$phi_copy;
      eta = foo$eta;
    }
    if(all(break_conditions)) {
      break;
    }
  }
  if(!"stanfit"%in%class(stan_fit)) {
    curr_fit$value;
  } else {
    list(accepted_divergences = accepted_divergences,
         max_divergences = max_divergences,
         max_rhat = max_rhat,
         curr_beta0 = curr_beta0,
         curr_beta = curr_beta,
         theta_orig = theta_orig,
         theta_aug = theta_aug,
         phi = phi,
         eta = eta);
  }
}

#DESCRIPTION: Helper function to make projection matrix (or list of projection matrices), which is P in the notation 
#of Boonstra and Barbaro
#
#VALUE: A list as long as the the length of imputes_list, with each element containing a different projection matrix 
#using the indices of the imputations specified in the corresponding element of imputes_list.
#
#ARGUMENTS:
#
#x_curr_orig, x_curr_aug (matrices) matrices with equal numbers of rows and p & q columns, 
#respectively. These are used to estimate the joint association between the original 
#and augmented covariates, which is needed for the imputation model
#
#eigenvec_hist_var (matrix) pxp matrix with each row corresponding to an eigenvector. 
#This is v_0 in Boonstra and Barbaro.
#
#imputes_list (list) list of length-2 vectors, with each vector containing the lower and upper indices of the imputations to use 
#for a projection matrix in the SAB method. This is best explained with an example: if imputes_list = list(c(1,15),c(16,100),c(1,100)), 
#then three projection matrices will be returned. One will be based upon the first 15 imputations from a run of MICE, the second based upon
#the last 85 imputations from that same run (i.e. the 16th-100th imputations), and the third will be based upon all 100 imputations from this
#same run. This is coded as such to allow for flexible exploration of the impact of number of imputations or variability due to imputations. 
#
#seed_start (pos. integer) random seed to start each imputation 
#
#predictorMatrix: (matrix) (p+q)x(p+q) matrix equivalent to argument of the same name in
#in the 'mice()' function (type '?mice'). It is specially calculated based upon a monotone
#missingness pattern (x^o is fully observed, x^a is not) Thus it samples from
#[X^a_{1...q}|X^o] = [X^a_1|X^o]*[X^a_2|X^o,X^a_1]*...

create_projection = function(x_curr_orig,
                             x_curr_aug,
                             eigenvec_hist_var,
                             imputes_list = list(c(1,15)),
                             seed_start = sample(.Machine$integer.max,1),
                             predictorMatrix = NULL) {
  require(mice);
  require(magrittr);
  p = ncol(x_curr_orig);
  q = ncol(x_curr_aug);
  orig_covariates = colnames(x_curr_orig);
  aug_covariates = colnames(x_curr_aug);
  x_all = cbind(x_curr_orig,x_curr_aug);
  stopifnot(class(imputes_list) == "list");
  
  n_imputes = max(unlist(lapply(imputes_list,max)));
  
  if(is.null(predictorMatrix)) {  
    predictorMatrix11 = matrix(0,nrow = p,ncol = p);
    predictorMatrix12 = matrix(0,nrow = p,ncol = q);
    predictorMatrix21 = matrix(1,nrow = q,ncol = p);
    predictorMatrix22 = matrix(1,nrow = q,ncol = q);
    predictorMatrix22[upper.tri(predictorMatrix22,diag = T)] = 0;
    predictorMatrix = rbind(cbind(predictorMatrix11,predictorMatrix12),
                            cbind(predictorMatrix21,predictorMatrix22));
    colnames(predictorMatrix) = rownames(predictorMatrix) = c(orig_covariates,aug_covariates);
    rm(predictorMatrix11,predictorMatrix12,predictorMatrix21,predictorMatrix22);
  } else {
    stopifnot(dim(predictorMatrix) == c(p+q,p+q));
  }
  
  #Create data to be marginalized
  #Column 1 identifies each eigenvector, including an additional row for the intercept offset;
  #Columns 2:(p+1) are original covariates;
  #Columns (p+2):(p+q+1) are the augmented covariates to be imputed and averaged;
  dat_for_marg = cbind(c(1:p,nrow(x_all)+1),
                       rbind(eigenvec_hist_var,0),
                       matrix(NA,p+1,q));
  colnames(dat_for_marg) = c("ID",orig_covariates,aug_covariates);
  dat_for_marg = data.matrix(dat_for_marg);
  #Store one copy for each of the different imputations to use
  dat_for_marg = lapply(1:length(imputes_list), function(x) dat_for_marg);
  
  #Now impute the augmented covariates based upon the artificially created original covariates, i.e. the eigenvectors, using
  #the empiric associations in current data. 
  curr_impute = mice(rbind(x_all,
                           #original covariates are constant across the different imputations, so we can just use the first
                           dat_for_marg[[1]][,c(orig_covariates,aug_covariates)]),
                     printFlag = F,
                     predictorMatrix = predictorMatrix,#user-provided, or if missing, created above
                     m = n_imputes,#
                     maxit = 1,#only one iteration is needed, due to the monotone missingness pattern
                     method = "pmm", 
                     seed = seed_start);
  curr_impute = data.matrix(mice::complete(curr_impute,"long"));
  for(k in 1:(p+1)) {
    for(j in 1:length(imputes_list)) {
      #Fill in the data_for marg
      dat_for_marg[[j]][dat_for_marg[[j]][,"ID"] == dat_for_marg[[j]][k,"ID"], aug_covariates] = 
        colMeans(curr_impute[which(curr_impute[,".id"] == k + nrow(x_all))[imputes_list[[j]][1]:imputes_list[[j]][2]],aug_covariates,drop=F]);
    }
  }
  
  #If there is any perfect correlation between predictors, the imputer will return NA. In this case, we fill in the missing data
  #with whatever variable was perfectly correlated with it
  for(j in 1:length(imputes_list)) {
    while(any(is.na(dat_for_marg[[j]][,aug_covariates]))) {
      correlated_column = which(colSums(is.na(dat_for_marg[[j]][,aug_covariates]))>0)[1];
      dat_for_marg[[j]][,names(correlated_column)] =  dat_for_marg[[j]][,setdiff(colnames(x_all)[abs(cor(x_all,x_curr_aug[,correlated_column]) - 1) < sqrt(.Machine$double.eps)],names(correlated_column))[1]];
    }
  }
  #v0_inv is constant across the different imputation sets, so we can just use the first set
  v0_inv = solve(dat_for_marg[[1]][1:p,orig_covariates,drop=F]);
  projections = vector("list",length(imputes_list));
  for(j in 1:length(imputes_list)) {
    projections[[j]] = v0_inv %*% (dat_for_marg[[j]][1:p,aug_covariates,drop=F] - dat_for_marg[[j]][rep(p+1,p),aug_covariates,drop=F]);
  }
  projections;
}

#DESCRIPTION: This is the parent function in the simulation study. For a given data-generating mechanism (characterized by choices 
#of n_hist, n_curr, true_mu_hist, true_mu_curr, true_betas_orig, true_betas_aug, and covariate_args) and modeling choice (characterized 
#by choices of local_dof, global_dof, slab_precision, nab_augmented_scale, power_prop_nonzero_prior, and sab_imputes_list, phi_params),
#all of the methods in Boonstra and Barbaro are run against an arbitrary number of simulated datasets. The user can modify various
#characteristics of the underlying HMC chain. A number of operating characteristics are returned, based both on estimation and prediction. 
#
#
#VALUE: A list of various results. 
#
#
#ARGUMENTS:
#sim_number (arbitrary) This is a label to help the user keep track of multiple different scenario settings
#It is simply returned at the end of the function
#
#array_id (pos. integer) This is intended to be the slurm array id
#
#n_sim (pos. integer) number of simulated datasets to construct
#
#n_hist (pos. integer) size of historical data; n_hist in the paper
#
#n_curr (pos. integer) size of current data; n_curr in the paper
#
#n_new (pos. integer) size of testing data (for prediction)
#
#true_mu_hist (real) true intercept for generating historical model. mu_hist in Boonstra and Barbaro.
#
#true_mu_curr (real) true intercept for generating current / new model. mu in Boonstra and Barbaro.
#
#true_betas_orig (vector) true regression coefficients corresponding to original covariates. Beta^o in Boonstra and Barbaro
#
#true_betas_aug (vector) true regression coefficients corresponding to augmented covariates. Beta^a in Boonstra and Barbaro
#
#covariate_args (list) the named arguments are simple ways to govern the distribution of the predictors. Beta^a in Boonstra and Barbaro.
#
#beta_label (arbitrary) This is a label to help the user keep track of multiple different true regression 
#coefficients. It is simply returned at the end of the function.
#
#covariate_args (list) the named arguments are simple ways to govern the distribution of the predictors. 
#x_correlation is the normal correlation between any pair of predictors; x_orig_binom is the 
#integer indices (any subset of 1, ..., length(true_betas_orig)) indicating which of the original 
#covariates should be transformed to binary values based upon being less than or greater than zero; 
#x_aug_binom is the analogous set of indices for the augmented predictors (any subset of 1, ..., length(true_betas_aug))
#
#covariate_label (arbitrary) This is a label to help the user keep track of multiple different true covariate 
#distributions It is simply returned at the end of the function
#
#local_dof, global_dof (pos. integer) numbers indicating the degrees of freedom for lambda_j and tau, respectively. Boonstra, 
#et al. never considered local_dof != 1 or global_dof != 1. 
#
#slab_precision (pos. real) the slab-part of the regularized horseshoe, this is equivalent to (1/d)^2 in the notation of
#Boonstra and Barbaro 
#
#nab_augmented_scale (pos. real) the scale parameter to accompany the tilde_lambda shrinkage of the augmented covariates
#when phi = 1 under NAB. This is equivalent to tilde_c in the notation of Boonstra and Barbaro.  
#
#power_prop_nonzero_prior (pos. real in [0,1]) exponent for number covariates that are assumed to be non-zero minus 1/2, 
#e.g. p^(1/3) - 0.5
#
#sab_imputes_list (list) list of length-2 vectors, with each vector containing the lower and upper indices of the imputations to use 
#for a projection matrix in the SAB method. This is best explained with an example: if sab_imputes_list = list(c(1,15),c(16,100),c(1,100)), 
#then three projection matrices will be constructed. One will be based upon the first 15 imputations from a run of MICE, the second based upon
#the last 85 imputations from that same run (i.e. the 16th-100th imputations), and the third will be based upon all 100 imputations from this
#same run. This is coded as such to allow for flexible exploration of the impact of number of imputations or variability due to imputations. 
#Note that each of the projection matrices is crossed with each of the hyperpriors on phi (determined by the length of phi_params), i.e. 
#i.e. if phi_params has length 3 and sab_imputes_list has length 2, then there will be six versions of SAB run and reported. Because NAB 
#doesn't depend upon any imputation, there would be just 3 versions of NAB run and reported in this example. 
#
#stan_file_path (character) local path to directory containing stan files
#
#standard_stan_filename (character) file name for standard prior
#
#sab_stan_filename (character) file name for sab prior
#
#sab_dev_stan_filename (character) file name for development version of sab prior (for testing)
#
#nab_stan_filename (character) file name for nab prior
#
#nab_dev_stan_filename  (character) file name for development version of nab prior (for testing)
#
#phi_params (list) list of named lists, each with named components 'mean' and 'sd', corresponding to the mean and standard deviation, respectively, 
#of a normal hyperprior on phi truncated to the [0,1] interval. Its length is the number of different hyperpriors on phi to be explored 
#for each adaptive bayesian update. In Boonstra and Barbaro, 'Agnostic' corresponds to list(mean = 0.5, sd = 2.5), which is approximately uniform, 
#and 'Optimistic' corresponds to list(mean = 1, sd = 0.25). Note that each of the hyperpriors on phi is crossed with each of the constructed 
#projection matrices (determined by the length of sab_imputes_list), i.e. if phi_params has length 3 and sab_imputes_list has length 2, 
#then there will be six versions of SAB run and reported. Because NAB doesn't depend upon any imputation, there would be just 3 versions of
#NAB run and reported in this example. 
#
#mc_warmup (pos. integer) equivalent to warmup in 'stan()' function (type '?stan')
#
#mc_iter_after_warmup (pos. integer) equivalent to iter - warmup in 'stan()' function (type '?stan')
#
#mc_chains (pos. integer) equivalent to chains in 'stan()' function (type '?stan')
#
#mc_thin (pos. integer) equivalent to thin in 'stan()' function (type '?stan')
#
#mc_stepsize (pos. real) equivalent to stepsize in 'stan()' function (type '?stan')
#
#mc_adapt_delta_relaxed, mc_adapt_delta,strict (pos. real in [0,1]) Two alternatives to use as values for adapt_delta 
#in 'stan()' function (type '?stan'). The relaxed version is presumably smaller and will be used for the 
#standard prior, which has fewer numerical issues. The strict version is for the nab and sab priors.
#
#mc_max_treedepth (pos. integer) equivalent to max_treedepth in 'stan()' function (type '?stan')
#
#ntries_per_iter (pos. integer) each method will run up to this many times, stopping either when the number of *divergent transitions*
#is zero or when ntries has been reached. The reported fit will be that with the fewest number of divergent iterations. 
#
#random_seed (pos. integer) where to initialize the simulator
#
#fit_marginal (logical) Set to TRUE to include an asymptotic estimate of the misspecified alpha_orig
#
#fit_methods  (logical) Set to FALSE to do a dry run of the simulator 
#
#skip_methods (vector) any methods falling in the intersection of skip_methods and 
#c("Benchmark","Historical","Standard","NAB","NAB_dev","SAB","SAB_dev") will be skipped
#
#dynamic_run (logical) set to TRUE if stepping through this function and store more results

run.sim <- function(sim_number,
                    array_id,
                    n_sim,
                    n_hist,
                    n_curr,
                    n_new,
                    true_mu_hist,
                    true_mu_curr,
                    true_betas_orig,
                    true_betas_aug,
                    beta_label,
                    covariate_args,
                    covariate_label,
                    local_dof = 1,
                    global_dof = 1,
                    slab_precision = (1.0/15.0)^2,
                    nab_augmented_scale = 0.05,
                    power_prop_nonzero_prior = 1/3,
                    sab_imputes_list = list(c(1,100)),
                    stan_file_path = "",
                    standard_stan_filename = "RegHS_stable.stan",
                    sab_stan_filename = "SAB_stable.stan",
                    sab_dev_stan_filename = "SAB_dev.stan",
                    nab_stan_filename = "NAB_stable.stan",
                    nab_dev_stan_filename = "NAB_dev.stan",
                    phi_params = list("Agnostic" = c(mean = 0.5, sd = 2.5),
                                      "Optimist" = c(mean = 1, sd = 0.25)),
                    mc_warmup = 1e3,
                    mc_iter_after_warmup = 5e3, 
                    mc_chains = 2, 
                    mc_thin = 1,
                    mc_stepsize = 0.01,
                    mc_adapt_delta_relaxed = 0.999,
                    mc_adapt_delta_strict = 0.9999999,
                    mc_max_treedepth = 20,
                    ntries_per_iter = 4,
                    random_seed = sample(.Machine$integer.max,1),
                    fit_marginal = T,#Set to TRUE to use all data to come up with an asymptotic estimate of the misspecified alpha_orig
                    fit_methods = T,#Set to FALSE to do a dry run of the simulator 
                    skip_methods = c("NAB_dev","SAB_dev"),#methods falling in the intersection of skip_methods with c("Benchmark","Historical","Standard","NAB","NAB_dev","SAB","SAB_dev") will be skipped
                    dynamic_run = T) { 
  
  begin_all = Sys.time();
  set.seed(random_seed);
  data_seeds = sample(2^30.999,n_sim);
  informational_messages = list();
  
  stopifnot(class(sab_imputes_list) == "list");
  sab_num_imputes_each = unlist(lapply(sab_imputes_list,diff)) + 1;
  max_sab_index = max(unlist(lapply(sab_imputes_list,max)));
  min_sab_index = min(unlist(lapply(sab_imputes_list,min)));
  if(min_sab_index != 1) {
    stop("The argument 'sab_imputes_list' requires that the smallest index overall be 1");
  }

  base_meth_names = c("Benchmark",
                      "Historical",
                      "Standard",
                      "NAB",
                      "NAB_dev",
                      "SAB",
                      "SAB_dev");
  sab_suffix = expand.grid(paste0(".phi",names(phi_params)),
                           paste0(".imp",1:length(sab_imputes_list)));
  expanded_meth_names = c("Benchmark",
                          "Historical",
                          "Standard",
                          paste0("NAB",names(phi_params)),
                          paste0("NAB_dev",names(phi_params)),
                          paste0("SAB",paste0(sab_suffix[,1],sab_suffix[,2])),
                          paste0("SAB_dev",paste0(sab_suffix[,1],sab_suffix[,2])));
  rm(sab_suffix);
  
  if("Historical" %in% skip_methods) {
    skip_methods = intersect(skip_methods, c("NAB","NAB_dev","SAB","SAB_dev"));
    informational_messages = c(informational_messages, paste0("Skipping all adaptive Bayesian methods because because Historical was skipped"));
  }
  if(("NAB" %in% skip_methods && !"NAB_dev" %in% skip_methods) || ("SAB" %in% skip_methods && !"SAB_dev" %in% skip_methods)) {
    stop("Code assumes that NAB and SAB are to be fit if their corresponding development versions are to be fit");
  }
  
  mse_beta = #mean squared error for all betas
    mse_beta_orig = #mean squared error for original betas
    mse_beta_aug = #mean squared error for augmented betas
    rmspe = #root mean squared prediction error
    mape = #mean absolute predict error
    mins_per_method = #minutes devoted to each method
    max_rhat = #
    matrix(NA,nrow=n_sim,ncol=length(expanded_meth_names),dimnames = list(NULL, expanded_meth_names));
  store_impute_time = rep(NA, n_sim);
  
  posterior_eff_orig = #posterior effective number of original parameters = mean(rowSums(1-kappa[orig]))
    posterior_eff_aug = #posterior effective number of augmented parameters = mean(rowSums(1-kappa[aug]))
    matrix(NA,nrow=n_sim,ncol=length(expanded_meth_names),dimnames = list(NULL,expanded_meth_names));
  #Store characteristics of true risk in the population
  population_params = 
    matrix(NA,
           nrow=n_sim,
           ncol=12,
           dimnames = list(NULL,c("mean_hist",#historical prevalence
                                  "mean_curr",#current prevalence
                                  "signal_comp",#(n_hist+n_curr) * var(lin_pred_all[historical data+current data])
                                  "signal_miss",#(n_hist) * var(lin_pred_aug[historical data])
                                  "historical_signal_obs",#(n_hist) * var(lin_pred_orig[historical data])
                                  "current_signal_obs",#(n_curr) * var(lin_pred_all[current data])
                                  "mean_pop_risk",#mean of true population risk
                                  "sd_pop_risk",#standard deviation of true population risk
                                  "lowerq_pop_risk",#lower and upper quartiles of risk
                                  "upperq_pop_risk",
                                  "log_det_xprod_x_comp",
                                  "log_det_xprod_x_miss")));
  
  num_orig = length(true_betas_orig);
  num_aug = length(true_betas_aug);
  orig_covariates = paste0("orig",1:num_orig);
  aug_covariates = paste0("aug",1:num_aug);
  
  #####For evaluating overall MSE
  matrix_true_beta = matrix(c(true_betas_orig,true_betas_aug),
                            nrow = mc_iter_after_warmup * mc_chains,
                            ncol = num_orig + num_aug,
                            byrow = T);
  #####For evaluating MSE of original betas only
  matrix_true_beta_orig = matrix(true_betas_orig,
                                 nrow = mc_iter_after_warmup * mc_chains,
                                 ncol = num_orig,
                                 byrow = T);
  #####For evaluating MSE of augmented betas only
  matrix_true_beta_aug = matrix(true_betas_aug,
                                nrow = mc_iter_after_warmup * mc_chains,
                                ncol = num_aug,
                                byrow = T);
  
  store_mean_beta = store_sd_beta = vector("list",length(expanded_meth_names));
  names(store_mean_beta) = names(store_sd_beta) = expanded_meth_names;
  for(k in 1:length(store_mean_beta)) {
    store_mean_beta[[k]] = matrix(NA,nrow = n_sim,ncol=num_orig + num_aug,dimnames = list(NULL,c(orig_covariates,aug_covariates)));
    store_sd_beta[[k]] = matrix(NA,nrow = n_sim,ncol=num_orig + num_aug,dimnames = list(NULL,c(orig_covariates,aug_covariates)));
  }
  rm(k);
  #store posterior means and sds of phi (shrinkage weight) and eta (scale factor on S_alpha)
  store_phi_mean =
    store_phi_sd = 
    store_eta_mean =
    store_eta_sd = matrix(NA, nrow = n_sim, 
                          ncol = length(grep("NAB",expanded_meth_names)) + 
                            length(grep("SAB",expanded_meth_names)),
                          dimnames = list(NULL, c(grep("NAB",expanded_meth_names,value=T),
                                                  grep("SAB",expanded_meth_names,value=T))));
  
  #Calculate scale parameters so that E[1-sum kappa_j] = (# of parameters to fit)^(power_prop_nonzero_prior) - 0.5;
  store_hierarchical_scales = 
    prior_eff = #prior effective number of original parameters = mean(rowSums(1-kappa[orig]))
    vector("list",length(base_meth_names));
  names(store_hierarchical_scales) = names(prior_eff) = base_meth_names;
  if(fit_methods) {
    #Benchmark: full access to data (n_hist + n_curr), but agnostic about covariates, so use skeptical hierarchical shrinkage 
    #with most assumed to be zero. 
    foo = solve_for_hiershrink_scale(target_mean1 = -0.5 + (num_orig + num_aug) ^ power_prop_nonzero_prior,
                                     target_mean2 = NA,
                                     npar1 = num_orig + num_aug, 
                                     npar2 = 0,
                                     local_dof = local_dof, 
                                     regional_dof = -Inf, 
                                     global_dof = global_dof,
                                     slab_precision = slab_precision,
                                     n = n_hist + n_curr,
                                     sigma = 2, 
                                     n_sim = round(2e6/(num_orig + num_aug)));
    store_hierarchical_scales$Benchmark = foo$scale1;
    prior_eff$Benchmark = foo$prior_num1;
    rm(foo);
    #Historical: access to historical data alone, with a skeptical prior. This is both a method in and of itself as well as
    #the 'prior' analysis that will be provided to the SAB methods. 
    foo = solve_for_hiershrink_scale(target_mean1 = -0.5 + num_orig ^ power_prop_nonzero_prior,
                                     target_mean2 = NA,
                                     npar1 = num_orig, 
                                     npar2 = 0,
                                     local_dof = local_dof, 
                                     regional_dof = -Inf, 
                                     global_dof = global_dof,
                                     slab_precision = slab_precision,
                                     n = n_hist,
                                     sigma = 2, 
                                     n_sim = round(2e6/(num_orig + num_aug)));
    store_hierarchical_scales$Historical = foo$scale1;
    prior_eff$Historical = foo$prior_num1;
    rm(foo);
    #Standard: access to current data alone, with a skeptical prior. 
    #This also corresponds to the standard scales for SAB and NAB, where standard means that
    #no historical information is used
    foo = solve_for_hiershrink_scale(target_mean1 = -0.5 + (num_orig + num_aug) ^ power_prop_nonzero_prior,
                                     target_mean2 = NA,
                                     npar1 = num_orig + num_aug, 
                                     npar2 = 0,
                                     local_dof = local_dof, 
                                     regional_dof = -Inf, 
                                     global_dof = global_dof,
                                     slab_precision = slab_precision,
                                     n = n_curr,
                                     sigma = 2, 
                                     n_sim = round(2e6/(num_orig + num_aug)));
    store_hierarchical_scales$Standard = 
      store_hierarchical_scales$NAB = 
      store_hierarchical_scales$SAB = 
      foo$scale1;
    prior_eff$Standard = 
      prior_eff$NAB = 
      prior_eff$SAB = 
      foo$prior_num1;
    rm(foo);
    #
    store_hierarchical_scales$NAB_aug_tilde = nab_augmented_scale;
  }
  #This creates a predictor matrix to pass to MICE for creating the projection matrix
  #The assumption is of a monotone missingness pattern (x^o is fully observed, x^a is not)
  #Thus it samples fromm [X^a_{1...q}|X^o] = [X^a_1|X^o]*[X^a_2|X^o,X^a_1]*...
  predictorMatrix11 = matrix(0,nrow = num_orig,ncol = num_orig);
  predictorMatrix12 = matrix(0,nrow = num_orig,ncol = num_aug);
  predictorMatrix21 = matrix(1,nrow = num_aug,ncol = num_orig);
  predictorMatrix22 = matrix(1,nrow = num_aug,ncol = num_aug);
  predictorMatrix22[upper.tri(predictorMatrix22,diag = T)] = 0;
  predictorMatrix = rbind(cbind(predictorMatrix11,predictorMatrix12),
                          cbind(predictorMatrix21,predictorMatrix22));
  colnames(predictorMatrix) = rownames(predictorMatrix) = c(orig_covariates,aug_covariates);
  rm(predictorMatrix11,predictorMatrix12,predictorMatrix21,predictorMatrix22);
  
  ##Store information on divergent transitions
  max_divergences_by_method = matrix(-Inf, n_sim, length(expanded_meth_names),dimnames = list(NULL,expanded_meth_names));
  accepted_divergences_by_method = matrix(Inf, n_sim, length(expanded_meth_names),dimnames = list(NULL,expanded_meth_names));
  
  #Store asymptotically misspecified estimates of alpha_orig
  if(fit_marginal) {
    complete_dat = draw_data(n_hist = 5, 
                             n_curr = 5,
                             n_new = round(max(1e5, 2e7 / (num_orig + num_aug))),
                             true_mu_hist = true_mu_hist,
                             true_mu_curr = true_mu_curr,
                             true_betas_orig = true_betas_orig,
                             true_betas_aug = true_betas_aug,
                             covariate_args = covariate_args);
    true_alphas_orig =  round(as.numeric(coef(glm(complete_dat$y_new ~ complete_dat$x_new_orig, family="binomial"))[-1]),2);
    rm(complete_dat);
  } else {
    true_alphas_orig = rep(NA,num_orig);
  }
  
  stan_compiled = F;
  begin_sim = Sys.time();
  i=1;
  
  for(i in 1:n_sim) {
    set.seed(data_seeds[i]);
    
    stopifnot("x_correlation"%in%names(covariate_args) && covariate_args$x_correlation >= 0 && covariate_args$x_correlation < 1);
    
    complete_dat = draw_data(n_hist = n_hist, 
                             n_curr = n_curr,
                             n_new = n_new,
                             true_mu_hist = true_mu_hist,
                             true_mu_curr = true_mu_curr,
                             true_betas_orig = true_betas_orig,
                             true_betas_aug = true_betas_aug,
                             covariate_args = covariate_args);
    
    #####RMSPE
    matrix_risk_new = matrix(complete_dat$risk_new, nrow=mc_iter_after_warmup*mc_chains, ncol=n_new, byrow=T);
    
    y_hist = complete_dat$y_hist;
    y_curr = complete_dat$y_curr;
    x_hist_orig = as.matrix(complete_dat$x_hist_orig);
    x_curr_orig = as.matrix(complete_dat$x_curr_orig);
    colnames(x_hist_orig) = colnames(x_curr_orig) = orig_covariates;
    x_hist_aug = as.matrix(complete_dat$x_hist_aug);#will be unobserved, used for benchmarking performance
    x_curr_aug = as.matrix(complete_dat$x_curr_aug);
    colnames(x_hist_aug) = colnames(x_curr_aug) = aug_covariates;
    population_params[i,"mean_hist"] = mean(y_hist);
    population_params[i,"mean_curr"] = mean(y_curr);
    population_params[i,"signal_comp"] = (n_hist+n_curr) * var((complete_dat$lin_pred_x_orig + complete_dat$lin_pred_x_aug)[1:(n_hist+n_curr)]);
    population_params[i,"signal_miss"] = n_hist * var((complete_dat$lin_pred_x_aug)[1:n_hist]);
    population_params[i,"historical_signal_obs"] = n_hist * var((complete_dat$lin_pred_x_orig)[1:n_hist]);
    population_params[i,"current_signal_obs"] = n_curr * var((complete_dat$lin_pred_x_orig + complete_dat$lin_pred_x_aug)[(n_hist+1):(n_hist+n_curr)]);
    population_params[i,"mean_pop_risk"] = mean(complete_dat$risk_new);
    population_params[i,"sd_pop_risk"] = sd(complete_dat$risk_new);
    population_params[i,c("lowerq_pop_risk","upperq_pop_risk")] = quantile(complete_dat$risk_new,probs=c(0.25,0.75));
    population_params[i,"log_det_xprod_x_comp"] = determinant(crossprod(rbind(cbind(x_hist_orig, x_hist_aug),
                                                                              cbind(x_curr_orig, x_curr_aug))))$modulus;
    population_params[i,"log_det_xprod_x_miss"] = determinant(crossprod(x_hist_aug))$modulus;
    
    if(fit_methods) {
      
      ## Compile templates ====##########################################################################
      if(!stan_compiled) {
        begin_compile = Sys.time();
        
        assign("Standard_template",glm_standard(stan_path = paste0(stan_file_path,standard_stan_filename)));
        assign("NAB_template",glm_nab(stan_path = paste0(stan_file_path,nab_stan_filename)));
        if(!is.null(nab_dev_stan_filename) && (!"NAB_dev" %in% skip_methods)) {
          assign("NAB_dev_template",glm_nab(stan_path = paste0(stan_file_path,nab_dev_stan_filename)));
        }
        assign("SAB_template",glm_sab(stan_path = paste0(stan_file_path,sab_stan_filename)));
        if(!is.null(sab_dev_stan_filename) && (!"SAB_dev" %in% skip_methods)) {
          assign("SAB_dev_template",glm_sab(stan_path = paste0(stan_file_path,sab_dev_stan_filename)));
        }
        end_compile = Sys.time();  
        stan_compiled = T;
        
      } 
      
      only_prior = F;
      p = num_orig;
      q = num_aug;
      
      ## Benchmark ====##########################################################################
      #Full knowledge of all data, accounting for different intercepts between studies
      curr_method = curr_base_method = "Benchmark";
      if(!curr_base_method %in% skip_methods) {
        begin_curr_method = Sys.time();
        y = c(y_hist, y_curr);
        #Allow for different intercepts between cohorts
        intercept_offset = c(rep(0,n_hist),
                             rep(1,n_curr));
        x_standardized = rbind(cbind(x_hist_orig, x_hist_aug),
                               cbind(x_curr_orig, x_curr_aug));
        beta_orig_scale  = 
          beta_aug_scale = 
          store_hierarchical_scales[[curr_method]];
        
        foo = glm_standard(stan_path = paste0(stan_file_path,standard_stan_filename), 
                           stan_fit = Standard_template,
                           y = y, 
                           x_standardized = x_standardized, 
                           p = p,
                           q = q,
                           beta_orig_scale = beta_orig_scale, 
                           beta_aug_scale = beta_aug_scale, 
                           local_dof = local_dof, 
                           global_dof = global_dof, 
                           slab_precision = slab_precision,
                           intercept_offset = intercept_offset,
                           only_prior = only_prior, 
                           mc_warmup = mc_warmup, 
                           mc_iter_after_warmup = mc_iter_after_warmup, 
                           mc_chains = mc_chains, 
                           mc_thin = mc_thin, 
                           mc_stepsize = mc_stepsize, 
                           mc_adapt_delta = mc_adapt_delta_relaxed,
                           mc_max_treedepth = mc_max_treedepth,
                           ntries = ntries_per_iter);
        
        max_divergences_by_method[i,curr_method] = foo$max_divergences;
        accepted_divergences_by_method[i,curr_method] = foo$accepted_divergences;
        max_rhat[i,curr_method] = foo$max_rhat;
        
        curr_beta0 = foo$curr_beta0;
        curr_beta = foo$curr_beta;
        
        ##
        store_mean_beta[[curr_method]][i,] = colMeans(curr_beta);
        store_sd_beta[[curr_method]][i,] = apply(curr_beta,2,sd);
        curr_kappa = 1 / (1 + 0.25 * length(y) * cbind(foo$theta_orig,foo$theta_aug)^2);
        posterior_eff_orig[i,curr_method] = mean(rowSums(1 - curr_kappa[,1:p,drop=F]));
        posterior_eff_aug[i,curr_method] = mean(rowSums(1 - curr_kappa[,(p+1):(p+q),drop=F]));
        
        mse_beta[i,curr_method] = mean(rowSums((curr_beta - matrix_true_beta)^2));
        mse_beta_orig[i,curr_method] = mean(rowSums((curr_beta[,1:p,drop=F] - matrix_true_beta_orig)^2));
        mse_beta_aug[i,curr_method] = mean(rowSums((curr_beta[,(p+1):(q+p),drop=F] - matrix_true_beta_aug)^2));
        
        fitted_new = expit(cbind(curr_beta0,curr_beta)%*%t(cbind(1,complete_dat$x_new_orig,complete_dat$x_new_aug)));
        rmspe[i,curr_method] = mean(sqrt(colMeans((fitted_new - matrix_risk_new)^2)));
        mape[i,curr_method] = mean(abs(fitted_new - matrix_risk_new));
        
        if(dynamic_run) {
          #Keep copy of values for further study 
          assign(paste0("beta0_",curr_method),curr_beta0);
          assign(paste0("beta_",curr_method),curr_beta);
        }
        
        mins_per_method[i,curr_method] = difftime(Sys.time(), begin_curr_method, units = "mins");
        rm(begin_curr_method,beta_orig_scale,beta_aug_scale,y,x_standardized,intercept_offset,foo,curr_beta,curr_beta0,fitted_new,curr_kappa);
      } else {
        if(i == 1) {
          informational_messages = c(informational_messages,paste0("Skipping ", curr_base_method, " because it is in 'skip_methods' arg"));
        }
      }
      rm(curr_method,curr_base_method);
      
      ## Historical ====##########################################################################
      #Full knowledge of historical data
      curr_base_method = curr_method = "Historical";
      if(!curr_base_method %in% skip_methods) {
        begin_curr_method = Sys.time();
        y = y_hist;
        x_standardized = x_hist_orig;
        intercept_offset = rep(0,n_hist);
        beta_orig_scale = 
          beta_aug_scale = 
          store_hierarchical_scales[[curr_method]];
        
        foo = glm_standard(stan_path = paste0(stan_file_path,standard_stan_filename), 
                           stan_fit = Standard_template,
                           y = y, 
                           x_standardized = x_standardized, 
                           p = p,
                           q = 0,#Note this difference from the other models
                           beta_orig_scale = beta_orig_scale, 
                           beta_aug_scale = beta_aug_scale, 
                           local_dof = local_dof, 
                           global_dof = global_dof, 
                           slab_precision = slab_precision,
                           intercept_offset = intercept_offset,
                           only_prior = only_prior, 
                           mc_warmup = mc_warmup, 
                           mc_iter_after_warmup = mc_iter_after_warmup, 
                           mc_chains = mc_chains, 
                           mc_thin = mc_thin, 
                           mc_stepsize = mc_stepsize, 
                           mc_adapt_delta = mc_adapt_delta_relaxed,
                           mc_max_treedepth = mc_max_treedepth,
                           ntries = ntries_per_iter);
        
        max_divergences_by_method[i,curr_method] = foo$max_divergences;
        accepted_divergences_by_method[i,curr_method] = foo$accepted_divergences;
        max_rhat[i,curr_method] = foo$max_rhat;
        
        curr_beta0 = foo$hist_beta0;
        #Historical model does not estimate all betas; for purposes of performance, we 
        #interpret this as infinite shrinkage of beta_aug
        curr_beta = cbind(foo$curr_beta,matrix(0, nrow = nrow(foo$curr_beta), ncol = q));
        
        ##
        store_mean_beta[[curr_method]][i,] = colMeans(curr_beta);
        store_sd_beta[[curr_method]][i,] = apply(curr_beta,2,sd);
        curr_kappa = 1 / (1 + 0.25 * length(y) * foo$theta_orig^2);
        posterior_eff_orig[i,curr_method] = mean(rowSums(1 - curr_kappa[,1:p,drop=F]));
        posterior_eff_aug[i,curr_method] = 0;
        
        mse_beta[i,curr_method] = mean(rowSums((curr_beta - matrix_true_beta)^2));
        mse_beta_orig[i,curr_method] = mean(rowSums((curr_beta[,1:p,drop=F] - matrix_true_beta_orig)^2));
        mse_beta_aug[i,curr_method] = mean(rowSums((curr_beta[,(p+1):(q+p),drop=F] - matrix_true_beta_aug)^2));
        
        fitted_new = expit(cbind(curr_beta0,curr_beta)%*%t(cbind(1,complete_dat$x_new_orig,complete_dat$x_new_aug)));
        rmspe[i,curr_method] = mean(sqrt(colMeans((fitted_new - matrix_risk_new)^2)));
        mape[i,curr_method] = mean(abs(fitted_new - matrix_risk_new));
        
        ##Keep copy of values for use in Adaptive Bayesian Updates
        assign(paste0("beta0_",curr_method),curr_beta0);
        assign(paste0("beta_",curr_method),curr_beta[,1:p,drop = F]);
        
        mins_per_method[i,curr_method] = difftime(Sys.time(), begin_curr_method, units = "mins");
        rm(begin_curr_method,beta_orig_scale,beta_aug_scale,y,x_standardized,intercept_offset,foo,curr_beta,curr_beta0,fitted_new,curr_kappa);
      } else {
        if(i == 1) {
          informational_messages = c(informational_messages,paste0("Skipping ", curr_base_method, " because it is in 'skip_methods' arg"));
        }
      }
      rm(curr_method,curr_base_method);
      
      ## Standard ====##########################################################################
      #Standard analysis of current data, ignoring historical model
      curr_base_method = curr_method = "Standard";
      if(!curr_base_method %in% skip_methods) {
        begin_curr_method = Sys.time();
        y = y_curr;
        x_standardized = cbind(x_curr_orig,x_curr_aug);
        intercept_offset = rep(0,n_curr);
        beta_orig_scale = 
          beta_aug_scale = 
          store_hierarchical_scales[[curr_method]];
        
        foo = glm_standard(stan_path = paste0(stan_file_path,standard_stan_filename), 
                           stan_fit = Standard_template,
                           y = y, 
                           x_standardized = x_standardized, 
                           p = p,
                           q = q,
                           beta_orig_scale = beta_orig_scale, 
                           beta_aug_scale = beta_aug_scale, 
                           local_dof = local_dof, 
                           global_dof = global_dof, 
                           slab_precision = slab_precision,
                           intercept_offset = intercept_offset,
                           only_prior = only_prior, 
                           mc_warmup = mc_warmup, 
                           mc_iter_after_warmup = mc_iter_after_warmup, 
                           mc_chains = mc_chains, 
                           mc_thin = mc_thin, 
                           mc_stepsize = mc_stepsize, 
                           mc_adapt_delta = mc_adapt_delta_relaxed,
                           mc_max_treedepth = mc_max_treedepth,
                           ntries = ntries_per_iter);
        
        max_divergences_by_method[i,curr_method] = foo$max_divergences;
        accepted_divergences_by_method[i,curr_method] = foo$accepted_divergences;
        max_rhat[i,curr_method] = foo$max_rhat;
        
        #The curr_beta0 assignment below is not a typo: 'hist_beta0' corresponds to those for whom 'intercept_offset = 0', 
        #and 'curr_beta0' corresponds to those for whom 'intercept_offset = 1'. Because nobody falls in the latter category, 
        #the parameter that controls this offset is not identified and should not be included. 
        curr_beta0 = foo$hist_beta0;
        curr_beta = foo$curr_beta;
        
        ##
        store_mean_beta[[curr_method]][i,] = colMeans(curr_beta);
        store_sd_beta[[curr_method]][i,] = apply(curr_beta,2,sd);
        curr_kappa = 1 / (1 + 0.25 * length(y) * cbind(foo$theta_orig,foo$theta_aug)^2);
        posterior_eff_orig[i,curr_method] = mean(rowSums(1 - curr_kappa[,1:p,drop=F]));
        posterior_eff_aug[i,curr_method] = mean(rowSums(1 - curr_kappa[,(p+1):(p+q),drop=F]));
        
        mse_beta[i,curr_method] = mean(rowSums((curr_beta - matrix_true_beta)^2));
        mse_beta_orig[i,curr_method] = mean(rowSums((curr_beta[,1:p,drop=F] - matrix_true_beta_orig)^2));
        mse_beta_aug[i,curr_method] = mean(rowSums((curr_beta[,(p+1):(q+p),drop=F] - matrix_true_beta_aug)^2));
        
        fitted_new = expit(cbind(curr_beta0,curr_beta)%*%t(cbind(1,complete_dat$x_new_orig,complete_dat$x_new_aug)));
        rmspe[i,curr_method] = mean(sqrt(colMeans((fitted_new - matrix_risk_new)^2)));
        mape[i,curr_method] = mean(abs(fitted_new - matrix_risk_new));
        
        if(dynamic_run) {
          #Keep copy of values for further study 
          assign(paste0("beta0_",curr_method),curr_beta0);
          assign(paste0("beta_",curr_method),curr_beta);
        }
        
        mins_per_method[i,curr_method] = difftime(Sys.time(), begin_curr_method, units = "mins");
        rm(begin_curr_method,beta_orig_scale,beta_aug_scale,y,x_standardized,intercept_offset,foo,curr_beta,curr_beta0,fitted_new,curr_kappa);
      } else {
        if(i == 1) {
          informational_messages = c(informational_messages,paste0("Skipping ", curr_base_method, " because it is in 'skip_methods' arg"));
        }
      }
      rm(curr_method,curr_base_method);
      
      ## NAB ====##########################################################################
      #Naive Adaptive Bayes: apply Historical analysis directly as a prior on beta_orig. 
      if(exists("beta_Historical")) {#Can only do these methods if beta_Historical was estimated
        #
        curr_base_method = "NAB";
        if(!curr_base_method %in% skip_methods) {
          y = y_curr;
          x_standardized = cbind(x_curr_orig,x_curr_aug);
          alpha_prior_mean = colMeans(beta_Historical);
          alpha_prior_cov = var(beta_Historical);
          
          beta_orig_scale = 
            beta_aug_scale = store_hierarchical_scales[[curr_base_method]];
          beta_aug_scale_tilde = store_hierarchical_scales[[paste0(curr_base_method,"_aug_tilde")]];

          scale_to_variance225 = diag(alpha_prior_cov) / 225;
          eigendecomp_hist_var = eigen(alpha_prior_cov);
          prior_type = names(phi_params)[1];
          
          for(prior_type in names(phi_params)) {
            begin_curr_method = Sys.time();
            curr_method = paste0(curr_base_method,prior_type);
            phi_mean = eval(phi_params[[prior_type]][["mean"]]);
            phi_sd = eval(phi_params[[prior_type]][["sd"]]);
            
            foo = glm_nab(stan_path = paste0(stan_file_path,nab_stan_filename),
                          stan_fit = NAB_template,
                          y = y, 
                          x_standardized = x_standardized, 
                          alpha_prior_mean = alpha_prior_mean,
                          alpha_prior_cov = alpha_prior_cov,
                          phi_mean = phi_mean,
                          phi_sd = phi_sd,
                          beta_orig_scale = beta_orig_scale, 
                          beta_aug_scale = beta_aug_scale, 
                          beta_aug_scale_tilde = beta_aug_scale_tilde,
                          local_dof = local_dof, 
                          global_dof = global_dof, 
                          slab_precision = slab_precision,
                          only_prior = only_prior, 
                          mc_warmup = mc_warmup, 
                          mc_iter_after_warmup = mc_iter_after_warmup, 
                          mc_chains = mc_chains, 
                          mc_thin = mc_thin, 
                          mc_stepsize = mc_stepsize, 
                          mc_adapt_delta = mc_adapt_delta_strict,
                          mc_max_treedepth = mc_max_treedepth,
                          ntries = ntries_per_iter,
                          eigendecomp_hist_var = eigendecomp_hist_var,
                          scale_to_variance225 = scale_to_variance225);
            
            max_divergences_by_method[i,curr_method] = foo$max_divergences;
            accepted_divergences_by_method[i,curr_method] = foo$accepted_divergences;
            max_rhat[i,curr_method] = foo$max_rhat;
            
            curr_beta0 = foo$curr_beta0;
            curr_beta = foo$curr_beta;
            curr_phi = foo$phi;
            curr_eta = foo$eta;
            
            ##
            store_mean_beta[[curr_method]][i,] = colMeans(curr_beta);
            store_sd_beta[[curr_method]][i,] = apply(curr_beta,2,sd);
            curr_kappa = 1 / (1 + 0.25 * length(y) * cbind(foo$theta_orig,foo$theta_aug)^2);
            posterior_eff_orig[i,curr_method] = mean(rowSums(1 - curr_kappa[,1:p,drop=F]));
            posterior_eff_aug[i,curr_method] = mean(rowSums(1 - curr_kappa[,(p+1):(p+q),drop=F]));
            
            mse_beta[i,curr_method] = mean(rowSums((curr_beta - matrix_true_beta)^2));
            mse_beta_orig[i,curr_method] = mean(rowSums((curr_beta[,1:p,drop=F] - matrix_true_beta_orig)^2));
            mse_beta_aug[i,curr_method] = mean(rowSums((curr_beta[,(p+1):(q+p),drop=F] - matrix_true_beta_aug)^2));
            
            fitted_new = expit(cbind(curr_beta0,curr_beta)%*%t(cbind(1,complete_dat$x_new_orig,complete_dat$x_new_aug)));
            rmspe[i,curr_method] = mean(sqrt(colMeans((fitted_new - matrix_risk_new)^2)));
            mape[i,curr_method] = mean(abs(fitted_new - matrix_risk_new));
            
            store_phi_mean[i,curr_method] = mean(curr_phi);
            store_phi_sd[i,curr_method] = sd(curr_phi);
            
            store_eta_mean[i,curr_method] = mean(curr_eta);
            store_eta_sd[i,curr_method] = sd(curr_eta);
            
            if(dynamic_run) {
              #Keep copy of values for further study 
              assign(paste0("beta0_",curr_method),curr_beta0);
              assign(paste0("beta_",curr_method),curr_beta);
            }
            
            mins_per_method[i,curr_method] = difftime(Sys.time(), begin_curr_method, units = "mins");
            rm(begin_curr_method,curr_method,phi_mean,phi_sd,foo,curr_beta,curr_beta0,fitted_new,curr_kappa,curr_phi,curr_eta);
          }
          rm(prior_type,y,x_standardized,alpha_prior_mean,alpha_prior_cov,beta_orig_scale,beta_aug_scale,beta_aug_scale_tilde,eigendecomp_hist_var,scale_to_variance225);
        } else {
          if(i == 1) {
            informational_messages = c(informational_messages,paste0("Skipping ", curr_base_method, " because it is in 'skip_methods' arg"));
          }
        }
        rm(curr_base_method);
      } 
      
      ## NAB_dev ====##########################################################################
      #Naive Adaptive Bayes (development version): for testing development version against stable version
      if(exists("beta_Historical")) {#Can only do these methods if beta_Historical was estimated
        #
        curr_base_method = "NAB";
        if(!paste0(curr_base_method,"_dev") %in% skip_methods) {
          y = y_curr;
          x_standardized = cbind(x_curr_orig,x_curr_aug);
          alpha_prior_mean = colMeans(beta_Historical);
          alpha_prior_cov = var(beta_Historical);
          
          beta_orig_scale = 
            beta_aug_scale = store_hierarchical_scales[[curr_base_method]];
          beta_aug_scale_tilde = store_hierarchical_scales[[paste0(curr_base_method,"_aug_tilde")]];
          prior_type = names(phi_params)[1];
          
          scale_to_variance225 = diag(alpha_prior_cov) / 225;
          eigendecomp_hist_var = eigen(alpha_prior_cov);
          prior_type = names(phi_params)[1];
          
          for(prior_type in names(phi_params)) {
            begin_curr_method = Sys.time();
            curr_method = paste0(curr_base_method,"_dev",prior_type);
            phi_mean = eval(phi_params[[prior_type]][["mean"]]);
            phi_sd = eval(phi_params[[prior_type]][["sd"]]);
            
            foo = glm_nab(stan_path = paste0(stan_file_path,nab_dev_stan_filename),
                          stan_fit = NAB_dev_template,
                          y = y, 
                          x_standardized = x_standardized, 
                          alpha_prior_mean = alpha_prior_mean,
                          alpha_prior_cov = alpha_prior_cov,
                          phi_mean = phi_mean,
                          phi_sd = phi_sd,
                          beta_orig_scale = beta_orig_scale, 
                          beta_aug_scale = beta_aug_scale, 
                          beta_aug_scale_tilde = beta_aug_scale_tilde,
                          local_dof = local_dof, 
                          global_dof = global_dof, 
                          slab_precision = slab_precision,
                          only_prior = only_prior, 
                          mc_warmup = mc_warmup, 
                          mc_iter_after_warmup = mc_iter_after_warmup, 
                          mc_chains = mc_chains, 
                          mc_thin = mc_thin, 
                          mc_stepsize = mc_stepsize, 
                          mc_adapt_delta = mc_adapt_delta_strict,
                          mc_max_treedepth = mc_max_treedepth,
                          ntries = ntries_per_iter,
                          eigendecomp_hist_var = eigendecomp_hist_var,
                          scale_to_variance225 = scale_to_variance225);
            
            max_divergences_by_method[i,curr_method] = foo$max_divergences;
            accepted_divergences_by_method[i,curr_method] = foo$accepted_divergences;
            max_rhat[i,curr_method] = foo$max_rhat;
            
            curr_beta0 = foo$curr_beta0;
            curr_beta = foo$curr_beta;
            curr_phi = foo$phi;
            curr_eta = foo$eta;
            
            ##
            store_mean_beta[[curr_method]][i,] = colMeans(curr_beta);
            store_sd_beta[[curr_method]][i,] = apply(curr_beta,2,sd);
            curr_kappa = 1 / (1 + 0.25 * length(y) * cbind(foo$theta_orig,foo$theta_aug)^2);
            posterior_eff_orig[i,curr_method] = mean(rowSums(1 - curr_kappa[,1:p,drop=F]));
            posterior_eff_aug[i,curr_method] = mean(rowSums(1 - curr_kappa[,(p+1):(p+q),drop=F]));
            
            mse_beta[i,curr_method] = mean(rowSums((curr_beta - matrix_true_beta)^2));
            mse_beta_orig[i,curr_method] = mean(rowSums((curr_beta[,1:p,drop=F] - matrix_true_beta_orig)^2));
            mse_beta_aug[i,curr_method] = mean(rowSums((curr_beta[,(p+1):(q+p),drop=F] - matrix_true_beta_aug)^2));
            
            fitted_new = expit(cbind(curr_beta0,curr_beta)%*%t(cbind(1,complete_dat$x_new_orig,complete_dat$x_new_aug)));
            rmspe[i,curr_method] = mean(sqrt(colMeans((fitted_new - matrix_risk_new)^2)));
            mape[i,curr_method] = mean(abs(fitted_new - matrix_risk_new));
            
            store_phi_mean[i,curr_method] = mean(curr_phi);
            store_phi_sd[i,curr_method] = sd(curr_phi);
            
            store_eta_mean[i,curr_method] = mean(curr_eta);
            store_eta_sd[i,curr_method] = sd(curr_eta);
            
            if(dynamic_run) {
              #Keep copy of values for further study 
              assign(paste0("beta0_",curr_method),curr_beta0);
              assign(paste0("beta_",curr_method),curr_beta);
            }
            
            mins_per_method[i,curr_method] = difftime(Sys.time(), begin_curr_method, units = "mins");
            rm(begin_curr_method,curr_method,phi_mean,phi_sd,foo,curr_beta,curr_beta0,fitted_new,curr_kappa,curr_phi,curr_eta);
          }
          rm(prior_type,y,x_standardized,alpha_prior_mean,alpha_prior_cov,beta_orig_scale,beta_aug_scale,beta_aug_scale_tilde,eigendecomp_hist_var,scale_to_variance225);
        } else {
          if(i == 1) {
            informational_messages = c(informational_messages,paste0("Skipping ", paste0(curr_base_method,"_dev"), " because it is in 'skip_methods' arg"));
          }
        }
        rm(curr_base_method);
      } 
      #all_current_objects = ls(all=T);
      
      ## SAB ====##########################################################################
      #Sensible Adaptive Bayes: apply Historical analysis as a prior on beta_orig + projection%*%beta_aug 
      if(exists("beta_Historical")) {#Can only do these methods if beta_Historical was estimated
        curr_base_method = "SAB";
        if(!curr_base_method %in% skip_methods) {
          ##
          impute_time = Sys.time();
          cat("Creating projection matrix for SAB...\n\n");
          eigendecomp_hist_var = eigen(var(beta_Historical));
          eigenvec_hist_var = t(eigendecomp_hist_var$vectors);
          aug_projection = create_projection(x_curr_orig = x_curr_orig,
                                             x_curr_aug = x_curr_aug,
                                             eigenvec_hist_var = eigenvec_hist_var,
                                             imputes_list = sab_imputes_list,
                                             seed_start = data_seeds[i],
                                             predictorMatrix = predictorMatrix);
          impute_time = difftime(Sys.time(), impute_time, units = "mins");
          store_impute_time[i] = impute_time;
          ##
          y = y_curr;
          x_standardized = cbind(x_curr_orig,x_curr_aug);
          alpha_prior_mean = colMeans(beta_Historical);
          alpha_prior_cov = var(beta_Historical);
          
          scale_to_variance225 = diag(alpha_prior_cov) / 225;
          beta_orig_scale = 
            beta_aug_scale = store_hierarchical_scales[[curr_base_method]];
          m = 1;
          
          for(m in 1:length(sab_imputes_list)) {
            
            prior_type = names(phi_params)[1];
            
            for(prior_type in names(phi_params)) {
              curr_method = paste0(curr_base_method,".phi",prior_type,".imp",m);
              begin_curr_method = Sys.time();
              phi_mean = eval(phi_params[[prior_type]][["mean"]]);
              phi_sd = eval(phi_params[[prior_type]][["sd"]]);
              
              foo = glm_sab(stan_path = paste0(stan_file_path,sab_stan_filename),
                            stan_fit = SAB_template,
                            y = y, 
                            x_standardized = x_standardized, 
                            alpha_prior_mean = alpha_prior_mean,
                            alpha_prior_cov = alpha_prior_cov,
                            aug_projection = aug_projection[[m]],#This is the only thing that varies with m
                            phi_mean = phi_mean,
                            phi_sd = phi_sd,
                            beta_orig_scale = beta_orig_scale, 
                            beta_aug_scale = beta_aug_scale, 
                            local_dof = local_dof, 
                            global_dof = global_dof, 
                            slab_precision = slab_precision, 
                            only_prior = only_prior, 
                            mc_warmup = mc_warmup, 
                            mc_iter_after_warmup = mc_iter_after_warmup, 
                            mc_chains = mc_chains, 
                            mc_thin = mc_thin, 
                            mc_stepsize = mc_stepsize, 
                            mc_adapt_delta = mc_adapt_delta_strict,
                            mc_max_treedepth = mc_max_treedepth,
                            ntries = ntries_per_iter,
                            eigendecomp_hist_var = eigendecomp_hist_var,
                            scale_to_variance225 = scale_to_variance225);
              max_divergences_by_method[i,curr_method] = foo$max_divergences;
              accepted_divergences_by_method[i,curr_method] = foo$accepted_divergences;
              max_rhat[i,curr_method] = foo$max_rhat;
              
              curr_beta0 = foo$curr_beta0;
              curr_beta = foo$curr_beta;
              curr_phi = foo$phi;
              curr_eta = foo$eta;
              
              ##
              store_mean_beta[[curr_method]][i,] = colMeans(curr_beta);
              store_sd_beta[[curr_method]][i,] = apply(curr_beta,2,sd);
              curr_kappa = 1 / (1 + 0.25 * length(y) * cbind(foo$theta_orig,foo$theta_aug)^2);
              posterior_eff_orig[i,curr_method] = mean(rowSums(1 - curr_kappa[,1:p,drop=F]));
              posterior_eff_aug[i,curr_method] = mean(rowSums(1 - curr_kappa[,(p+1):(p+q),drop=F]));
              
              mse_beta[i,curr_method] = mean(rowSums((curr_beta - matrix_true_beta)^2));
              mse_beta_orig[i,curr_method] = mean(rowSums((curr_beta[,1:p,drop=F] - matrix_true_beta_orig)^2));
              mse_beta_aug[i,curr_method] = mean(rowSums((curr_beta[,(p+1):(q+p),drop=F] - matrix_true_beta_aug)^2));
              
              fitted_new = expit(cbind(curr_beta0,curr_beta)%*%t(cbind(1,complete_dat$x_new_orig,complete_dat$x_new_aug)));
              rmspe[i,curr_method] = mean(sqrt(colMeans((fitted_new - matrix_risk_new)^2)));
              mape[i,curr_method] = mean(abs(fitted_new - matrix_risk_new));
              
              store_phi_mean[i,curr_method] = mean(curr_phi);
              store_phi_sd[i,curr_method] = sd(curr_phi);
              
              store_eta_mean[i,curr_method] = mean(curr_eta);
              store_eta_sd[i,curr_method] = sd(curr_eta);
              
              if(dynamic_run) {
                #Keep copy of values for further study 
                assign(paste0("beta0_",curr_method),curr_beta0);
                assign(paste0("beta_",curr_method),curr_beta);
              }
              mins_per_method[i,curr_method] = 
                difftime(Sys.time(), begin_curr_method, units = "mins") + #Stan running time
                (sab_num_imputes_each[[m]]/max_sab_index) * impute_time; #Prorated imputation running time
              rm(begin_curr_method,curr_method,foo,curr_beta,curr_beta0,fitted_new,curr_kappa,curr_phi,phi_mean, phi_sd,curr_eta);
            }
          }
          rm(m,prior_type,y,x_standardized,alpha_prior_mean,alpha_prior_cov,beta_orig_scale,beta_aug_scale,scale_to_variance225);
          #Keep this if the development version will also be run
          if(paste0(curr_base_method,"_dev") %in% skip_methods) {
            rm(impute_time,aug_projection,eigenvec_hist_var,eigendecomp_hist_var);
          }
        } else {
          if(i == 1) {
            informational_messages = c(informational_messages,paste0("Skipping ", curr_base_method, " because it is in 'skip_methods' arg"));
          }
        }
        rm(curr_base_method);
      } 
      #setdiff(ls(), all_current_objects);
      
      ## SAB_dev ====##########################################################################
      #Sensible Adaptive Bayes (development version): for testing development version against stable version
      if(exists("beta_Historical")) {#Can only do these methods if beta_Historical was estimated
        curr_base_method = "SAB";
        if(!paste0(curr_base_method,"_dev") %in% skip_methods) {
          y = y_curr;
          x_standardized = cbind(x_curr_orig,x_curr_aug);
          alpha_prior_mean = colMeans(beta_Historical);
          alpha_prior_cov = var(beta_Historical);
          
          scale_to_variance225 = diag(alpha_prior_cov) / 225;
          beta_orig_scale = 
            beta_aug_scale = store_hierarchical_scales[[curr_base_method]];
          m = 1;
          
          for(m in 1:length(sab_imputes_list)) {
            
            prior_type = names(phi_params)[1];
            
            for(prior_type in names(phi_params)) {
              #SAB
              curr_method = paste0(curr_base_method,"_dev.phi",prior_type,".imp",m);
              begin_curr_method = Sys.time();
              phi_mean = eval(phi_params[[prior_type]][["mean"]]);
              phi_sd = eval(phi_params[[prior_type]][["sd"]]);
              
              foo = glm_sab(stan_path = paste0(stan_file_path,sab_dev_stan_filename),
                            stan_fit = SAB_dev_template,
                            y = y, 
                            x_standardized = x_standardized, 
                            alpha_prior_mean = alpha_prior_mean,
                            alpha_prior_cov = alpha_prior_cov,
                            aug_projection = aug_projection[[m]],#This is the only thing that varies with m
                            phi_mean = phi_mean,
                            phi_sd = phi_sd,
                            beta_orig_scale = beta_orig_scale, 
                            beta_aug_scale = beta_aug_scale, 
                            local_dof = local_dof, 
                            global_dof = global_dof, 
                            slab_precision = slab_precision,
                            only_prior = only_prior, 
                            mc_warmup = mc_warmup, 
                            mc_iter_after_warmup = mc_iter_after_warmup, 
                            mc_chains = mc_chains, 
                            mc_thin = mc_thin, 
                            mc_stepsize = mc_stepsize, 
                            mc_adapt_delta = mc_adapt_delta_strict,
                            mc_max_treedepth = mc_max_treedepth,
                            ntries = ntries_per_iter,
                            eigendecomp_hist_var = eigendecomp_hist_var,
                            scale_to_variance225 = scale_to_variance225);
              max_divergences_by_method[i,curr_method] = foo$max_divergences;
              accepted_divergences_by_method[i,curr_method] = foo$accepted_divergences;
              max_rhat[i,curr_method] = foo$max_rhat;
              
              curr_beta0 = foo$curr_beta0;
              curr_beta = foo$curr_beta;
              curr_phi = foo$phi;
              curr_eta = foo$eta;
              
              ##
              store_mean_beta[[curr_method]][i,] = colMeans(curr_beta);
              store_sd_beta[[curr_method]][i,] = apply(curr_beta,2,sd);
              curr_kappa = 1 / (1 + 0.25 * length(y) * cbind(foo$theta_orig,foo$theta_aug)^2);
              posterior_eff_orig[i,curr_method] = mean(rowSums(1 - curr_kappa[,1:p,drop=F]));
              posterior_eff_aug[i,curr_method] = mean(rowSums(1 - curr_kappa[,(p+1):(p+q),drop=F]));
              
              mse_beta[i,curr_method] = mean(rowSums((curr_beta - matrix_true_beta)^2));
              mse_beta_orig[i,curr_method] = mean(rowSums((curr_beta[,1:p,drop=F] - matrix_true_beta_orig)^2));
              mse_beta_aug[i,curr_method] = mean(rowSums((curr_beta[,(p+1):(q+p),drop=F] - matrix_true_beta_aug)^2));
              
              fitted_new = expit(cbind(curr_beta0,curr_beta)%*%t(cbind(1,complete_dat$x_new_orig,complete_dat$x_new_aug)));
              rmspe[i,curr_method] = mean(sqrt(colMeans((fitted_new - matrix_risk_new)^2)));
              mape[i,curr_method] = mean(abs(fitted_new - matrix_risk_new));
              
              store_phi_mean[i,curr_method] = mean(curr_phi);
              store_phi_sd[i,curr_method] = sd(curr_phi);
              
              store_eta_mean[i,curr_method] = mean(curr_eta);
              store_eta_sd[i,curr_method] = sd(curr_eta);
              
              if(dynamic_run) {
                #Keep copy of values for further study 
                assign(paste0("beta0_",curr_method),curr_beta0);
                assign(paste0("beta_",curr_method),curr_beta);
              }
              mins_per_method[i,curr_method] = difftime(Sys.time(), begin_curr_method, units = "mins");
              rm(begin_curr_method,curr_method,foo,curr_beta,curr_beta0,fitted_new,curr_kappa,curr_phi,phi_mean,phi_sd,curr_eta);
            }
          }
          rm(m,prior_type,y,x_standardized,alpha_prior_mean,alpha_prior_cov,beta_orig_scale,beta_aug_scale,scale_to_variance225);
          rm(impute_time,aug_projection,eigenvec_hist_var,eigendecomp_hist_var);
        } else {
          if(i == 1) {
            informational_messages = c(informational_messages,paste0("Skipping ", paste0(curr_base_method,"_dev"), " because it is in 'skip_methods' arg"));
          }
        }
        rm(curr_base_method);
      } 
      #
      rm(matrix_risk_new);     
      cat(array_id,"-",i," proj. run time:",round(proj_run_time <- difftime(end_compile,begin_compile,units="hours") + (n_sim * (difftime(Sys.time(),begin_sim,units="hours") - difftime(end_compile,begin_compile,units="hours"))/i),2),"hours. proj. end time: ",format.Date(begin_sim  + proj_run_time),"\n\n",file="tracktime.txt",append=T);
    }
  }
  ## Summarize results ====##########################################################################
  
  mins_per_sim = difftime(Sys.time(),begin_sim,units="mins")/n_sim;
  
  average_accepted_divergences_by_method = colMeans(accepted_divergences_by_method);
  informational_messages = c(informational_messages, paste0("The range of the average divergences across methods was {", paste0(range(average_accepted_divergences_by_method,finite = T),collapse=","),"}"));
  
  mc_params = list(mc_warmup = mc_warmup, 
                   mc_iter_after_warmup = mc_iter_after_warmup, 
                   mc_chains = mc_chains, 
                   mc_thin = mc_thin,
                   mc_stepsize = mc_stepsize,
                   mc_adapt_delta_relaxed = mc_adapt_delta_relaxed,
                   mc_adapt_delta_strict = mc_adapt_delta_strict,
                   mc_max_treedepth = mc_max_treedepth,
                   ntries_per_iter = ntries_per_iter);
  
  model_params = list(standard_stan_filename = standard_stan_filename,
                      sab_stan_filename = sab_stan_filename,
                      sab_dev_stan_filename = sab_dev_stan_filename,
                      nab_stan_filename = nab_stan_filename,
                      nab_dev_stan_filename = nab_dev_stan_filename,
                      hierarchical_scale = store_hierarchical_scales,
                      local_dof = local_dof,
                      global_dof = global_dof,
                      slab_precision = slab_precision,
                      phi_params = phi_params,
                      sab_imputes_list = sab_imputes_list,
                      power_prop_nonzero_prior = power_prop_nonzero_prior);
  
  sim_params = list(sim_number = sim_number,
                    array_id = array_id,
                    n_sim = n_sim,
                    random_seed = random_seed,
                    data_seeds = data_seeds,
                    mins_toImpute_per_sim = store_impute_time,
                    mins_per_sim = mins_per_sim,
                    mins_per_method_per_sim = mins_per_method,
                    mins_total_runtime = difftime(Sys.time(),begin_all,units="mins"), 
                    time_start = begin_all,
                    time_end = Sys.time());
  list(sim_params = sim_params,
       max_divergences_by_method = max_divergences_by_method,
       accepted_divergences_by_method = accepted_divergences_by_method,
       max_rhat = max_rhat,
       generating_params = list(n_hist = n_hist,
                                n_curr = n_curr,
                                true_mu_hist = true_mu_hist,
                                true_mu_curr = true_mu_curr,
                                true_betas_orig = true_betas_orig,
                                true_betas_aug = true_betas_aug,
                                true_alphas_orig = true_alphas_orig,
                                beta_label = beta_label),
       covariate_args = c(covariate_args,covariate_label = covariate_label),
       model_params = model_params,
       mc_params = mc_params,
       population_params = population_params,
       mean_beta = store_mean_beta,
       sd_beta = store_sd_beta,
       prior_eff = prior_eff,
       posterior_eff_orig = posterior_eff_orig,
       posterior_eff_aug = posterior_eff_aug,
       phi_mean = store_phi_mean,
       phi_sd = store_phi_sd,
       eta_mean = store_eta_mean,
       eta_sd = store_eta_sd,
       mse_beta = mse_beta,
       mse_beta_orig = mse_beta_orig,
       mse_beta_aug = mse_beta_aug, 
       rmspe = rmspe,
       mape = mape,
       informational_messages = informational_messages);
}

#DESCRIPTION: This function calculates a numerical-based solution to the scale parameter c in the the equation three lines 
#from the top of page 7 in Section 2 of Boonstra and Barbaro. If desired, the user may request regional scale values
#for a partition of the covariates into two regions, defined by the first npar1 covariates and the second npar2 covariates.
#This functionality was not used in Boonstra and Barbaro. 
#
#
#ARGUMENTS:
#target_mean1, target_mean2 (pos. reals): the desired prior number of effective parameters (tilde xi_eff in Boonstra and Barbaro). 
#If one scale parameter is desired, leave target_mean2 = NA. An error will be thrown if target_mean1 > npar1 or if 
#target_mean2 > npar2. 
#
#npar1, npar2 (pos. integers): the number of covariates. If one scale parameter is required, then leave npar2 = 0. 
#
#local_dof, global_dof (pos. integer) numbers indicating the degrees of freedom for lambda_j and tau, respectively. Boonstra
#and Barbaro never considered local_dof != 1 or global_dof != 1. 
#
#regional_dof (pos. integer) Not used in Boonstra and Barbaro. If 
#
#n (pos. integer) sample size 
#
#sigma (pos. real) square root of the assumed dispersion. In Boonstra and Barbaro, this was always 2, corresponding to the
#maximum possible value: sqrt(1/[0.5 * (1 - 0.5)]). 
#
#tol (pos. real) numerical tolerance for convergence of solution
#
#max_iter (pos. integer) maximum number of iterations to run without convergence before giving up
#
#n_sim (pos. integer) number of simulated draws from the underlying student-t hyperpriors to calculate the Monte Carlo-based
#approximation of the expectation. 
#
#slab_precision (pos. real) the slab-part of the regularized horseshoe, this is equivalent to (1/d)^2 in the notation of
#Boonstra and Barbaro 


solve_for_hiershrink_scale = function(target_mean1,
                                      target_mean2 = NA,
                                      npar1,
                                      npar2 = 0,
                                      local_dof = 1,
                                      regional_dof = -Inf,
                                      global_dof = 1,
                                      slab_precision = (1/15)^2,
                                      n,
                                      sigma = 2,
                                      tol = .Machine$double.eps^0.5,
                                      max_iter = 100, 
                                      n_sim = 2e5
) {
  
  npar = npar1 + npar2;
  stopifnot(isTRUE(all.equal(npar%%1,0)));#Ensure integers
  do_local = (local_dof > 0);
  do_regional = (regional_dof > 0);
  do_global = (global_dof > 0);
  if(do_local) {
    lambda = matrix(rt(n_sim * npar,df = local_dof), nrow = n_sim);
  } else {
    lambda = matrix(1, nrow = n_sim, ncol = npar);
  }
  if(do_regional) {
    #Now do local-regional
    stopifnot(npar2 > 0 && isTRUE(all.equal(npar1%%1,0)) && isTRUE(all.equal(npar2%%1,0)));
    gamma_p = 1 + abs(rt(n_sim, df = regional_dof));
    gamma_q = 1 + abs(rt(n_sim, df = regional_dof));
    lambda[,1:npar1] = lambda[,1:npar1,drop=F] * gamma_p;
    lambda[,(npar1+1):(npar1+npar2)] = lambda[,(npar1+1):(npar1+npar2),drop=F] * gamma_q;
    rm(gamma_p,gamma_q);
  }
  if(do_global) {
    tau = rt(n_sim,df=global_dof);
    lambda[,1:npar1] = lambda[,1:npar1,drop=F] * tau;
    if(npar2 > 0) {
      lambda[,(npar1+1):(npar1+npar2)] = lambda[,(npar1+1):(npar1+npar2),drop=F] * tau;
    }
    rm(tau);
  } 
  if(is.na(target_mean2)) {
    npar1 = npar;
  }
  stopifnot(target_mean1 > 0 && target_mean1 < npar1);#Ensure proper bounds
  log_scale1 = diff_target = numeric(max_iter);
  log_scale1[1] = log(target_mean1/(npar1 - target_mean1)*sigma/sqrt(n));
  random_scales = 1 / (slab_precision + 1/(exp(2*log_scale1[1]) * lambda^2));
  kappa = 1/(1+n*random_scales[,1:npar1,drop=F]/sigma^2);
  diff_target[1] = mean(rowSums(1-kappa)) - target_mean1;
  log_scale1[2] = 0.02 + log_scale1[1];
  random_scales = 1 / (slab_precision + 1/(exp(2*log_scale1[2]) * lambda^2));
  kappa = 1/(1+n*random_scales[,1:npar1,drop=F]/sigma^2);
  diff_target[2] = mean(rowSums(1-kappa)) - target_mean1;
  i=2;
  while(T) {
    i = i+1;
    if(i > max_iter) {i = i-1; break;}
    log_scale1[i] = log_scale1[i-1] - diff_target[i-1]*(log_scale1[i-1]-log_scale1[i-2])/(diff_target[i-1]-diff_target[i-2]);
    random_scales = 1 / (slab_precision + 1/(exp(2*log_scale1[i]) * lambda^2));
    kappa = 1/(1+n*random_scales[,1:npar1,drop=F]/sigma^2);
    diff_target[i] = mean(rowSums(1-kappa)) - target_mean1;
    if(abs(diff_target[i]-diff_target[i-1]) < tol) {break;}
  }
  scale1 = exp(log_scale1[i]);
  diff1 = abs(diff_target[i]);
  iter1 = i;
  prior_num1 = diff_target[i] + target_mean1;
  
  if(!is.na(target_mean2)) {
    stopifnot(target_mean2 > 0 && target_mean2 < npar2);
    log_scale2 = diff_target = numeric(max_iter);
    log_scale2[1] = log(target_mean2/(npar2 - target_mean2)*sigma/sqrt(n));
    random_scales = 1 / (slab_precision + 1/(exp(2*log_scale2[1]) * lambda^2));
    kappa = 1/(1+n*random_scales[,(npar1+1):(npar1+npar2),drop=F]/sigma^2);
    diff_target[1] = mean(rowSums(1-kappa)) - target_mean2;
    log_scale2[2] = 0.02 + log_scale2[1];
    random_scales = 1 / (slab_precision + 1/(exp(2*log_scale2[2]) * lambda^2));
    kappa = 1/(1+n*random_scales[,(npar1+1):(npar1+npar2),drop=F]/sigma^2);
    diff_target[2] = mean(rowSums(1-kappa)) - target_mean2;
    i=2;
    while(T) {
      i = i+1;
      if(i > max_iter) {i = i-1; break;}
      log_scale2[i] = log_scale2[i-1] - diff_target[i-1]*(log_scale2[i-1]-log_scale2[i-2])/(diff_target[i-1]-diff_target[i-2]);
      random_scales = 1 / (slab_precision + 1/(exp(2*log_scale2[i]) * lambda^2));
      kappa = 1/(1+n*random_scales[,(npar1+1):(npar1+npar2),drop=F]/sigma^2);
      diff_target[i] = mean(rowSums(1-kappa)) - target_mean2;
      if(abs(diff_target[i]-diff_target[i-1]) < tol) {break;}
    }
    scale2 = exp(log_scale2[i]);
    diff2 = abs(diff_target[i]);
    iter2 = i;
    prior_num2 = diff_target[i] + target_mean2;
  } else {
    scale2 = NA;
    diff2 = NA;
    iter2 = NA;
    prior_num2 = NA;
  }
  
  return(list(scale1 = scale1,
              diff_from_target1 = diff1,
              iter1 = iter1,
              prior_num1 = prior_num1,
              scale2 = scale2,
              diff_from_target2 = diff2,
              iter2 = iter2,
              prior_num2 = prior_num2));
}

#DESCRIPTION: This function is the inverse of 'solve_for_hiershrink_scale'. Instead of providing a desired effective number of 
#parameters, the user provides the scale value(s), which is c in the notation of Boonstra and Barbaro, and the the function gives
#the implied prior number of effective parameters based upon this. As with 'solve_for_hiershrink_scale', the user can provide
#one global scale parameter (scale1, leaving scale2 = NA) that applies to all parameters, or two regional scale parameters (scale1, 
#scale2), that applies to a partition of the parameters as defined by the first npar1 parameters and the second npar2 parameters.
#
#
#ARGUMENTS:
#target_mean1, target_mean2 (pos. reals): the desired prior number of effective parameters (tilde xi_eff in Boonstra and Barbaro). 
#If one scale parameter is desired, leave target_mean2 = NA. An error will be thrown if target_mean1 > npar1 or if 
#target_mean2 > npar2. 

calculate_m_eff = function(scale1,
                           scale2 = NA,
                           npar1,
                           npar2 = 0,
                           local_dof = 1,
                           regional_dof = -Inf,
                           global_dof = 1,
                           slab_precision = (1/15)^2,
                           n,
                           sigma = 2,
                           tol = .Machine$double.eps^0.5,
                           max_iter = 100, 
                           n_sim = 2e5
) {
  do_local = (local_dof > 0);
  do_regional =(regional_dof > 0);
  do_global = (global_dof > 0);
  npar = npar1 + npar2;
  stopifnot(isTRUE(all.equal(npar%%1,0)));#Ensure integers
  if(do_local) {
    lambda = matrix(rt(n_sim * npar,df = local_dof), nrow = n_sim);
  } else {
    lambda = matrix(1, nrow = n_sim, ncol = npar);
  }
  if(do_regional) {
    #Now do local-regional
    stopifnot(npar2 > 0 && isTRUE(all.equal(npar1%%1,0)) && isTRUE(all.equal(npar2%%1,0)));
    gamma_p = 1 + abs(rt(n_sim, df = regional_dof));
    gamma_q = 1 + abs(rt(n_sim, df = regional_dof));
    lambda[,1:npar1] = lambda[,1:npar1,drop=F] * gamma_p;
    lambda[,(npar1+1):(npar1+npar2)] = lambda[,(npar1+1):(npar1+npar2),drop=F] * gamma_q;
    rm(gamma_p,gamma_q);
  }
  if(do_global) {
    tau = rt(n_sim,df=global_dof);
    lambda[,1:npar1] = lambda[,1:npar1,drop=F] * tau;
    if(npar2 > 0) {
      lambda[,(npar1+1):(npar1+npar2)] = lambda[,(npar1+1):(npar1+npar2),drop=F] * tau;
    }
    rm(tau);
  }
  random_scales = 1 / (slab_precision + 1/(scale1^2*lambda^2));
  if(is.na(scale2)) {
    npar1 = npar;
  }
  kappa1 = 1/(1+n*random_scales[,1:npar1,drop=F]/sigma^2);
  prior_num1 = mean(rowSums(1-kappa1))
  
  if(!is.na(scale2)) {
    random_scales = 1 / (slab_precision + 1/(scale2^2*lambda^2));
    kappa2 = 1/(1+n*random_scales[,(npar1+1):(npar1+npar2),drop=F]/sigma^2);
    prior_num2  = mean(rowSums(1-kappa2));
  } else {
    prior_num2 = NA;
  }
  
  return(list(prior_num1 = prior_num1,
              prior_num2 = prior_num2));
}


as.dummy = function(x,full_rank=T) {
  single.as.dummy <- function(x,full_rank) {
    levels_x = levels(x);
    1*matrix(rep(x,nlevels(x)) == rep(levels_x,each=length(x)),nrow=length(x),ncol=nlevels(x),dimnames=list(NULL,levels_x))[,(1+full_rank):length(levels_x),drop=F];
  }
  if("factor"%in%class(x)) {
    if(length(full_rank)>1) {warning("ignoring all but first element of 'full_rank'");}
    result = single.as.dummy(x,full_rank[1]);
  } else if("logical"%in%class(x)) {
    if(length(full_rank)>1) {warning("ignoring all but first element of 'full_rank'");}
    result = single.as.dummy(factor(x,levels=c(F,T)),full_rank[1]);
  } else if("integer"%in%class(x)) {
    if(length(full_rank)>1) {warning("ignoring all but first element of 'full_rank'");}
    result = single.as.dummy(factor(x,levels=sort(unique(x),decreasing = T)),full_rank[1]);
  } else if(class(x)=="data.frame") {
    result = NULL;
    full_rank = rep(full_rank,length=ncol(x));
    for(i in 1:ncol(x)) {
      if("factor"%in%class(x[,i])) {
        foo = single.as.dummy(x[,i],full_rank[i]);
        colnames(foo) = paste0(colnames(x)[i],colnames(foo));
        result = cbind(result,foo);
      } else if("logical"%in%class(x[,i])) {
        foo = single.as.dummy(factor(x[,i],levels=c(F,T)),full_rank[i]);
        colnames(foo) = paste0(colnames(x)[i],colnames(foo));
        result = cbind(result,foo);
      } else if("integer"%in%class(x[,i])) {
        foo = single.as.dummy(factor(x[,i],levels=sort(unique(x[,i]),decreasing = T)),full_rank[i]);
        colnames(foo) = paste0(colnames(x)[i],colnames(foo));
        result = cbind(result,foo);
      } else {
        stop("x must be either a factor, logical, or a dataframe comprised of factors/logicals/integers");
      }
    } 
  } else {
    stop("x must be either a factor, logical, or a dataframe comprised of factors/logicals/integers");
  }
  data.frame(result);
}

#DESCRIPTION: Program for fitting a GLM equipped with a regularized student-t prior on the regression coefficients, 
#parametrized using the normal-inverse-gamma distribution. The 'regularization' refers to the fact that the inverse-gamma
#scale is has a finite upper bound that it smoothly approaches. This method was not used in the simulation study but was used 
#in the data analysis. Specifically, it corresponds to 'PedRESC2'. 
#
#
#ARGUMENTS: (only those distinct from glm_standard are discussed)
#
#beta_scale (pos. real) constants indicating the prior scale of the student-t prior. 
#
#dof (pos. integer) degrees of freedom for the student-t prior

glm_studt = function(stan_fit = NA, 
                     stan_path,
                     y = c(0,1),
                     x_standardized = matrix(0,length(y),3), 
                     beta_scale = 1, 
                     dof = 1, 
                     slab_precision = (1/15)^2, 
                     only_prior = F, 
                     mc_warmup = 50, 
                     mc_iter_after_warmup = 50, 
                     mc_chains = 1, 
                     mc_thin = 1, 
                     mc_stepsize = 0.1, 
                     mc_adapt_delta = 0.9,
                     mc_max_treedepth = 15,
                     ntries = 1) {
  
  max_divergences = -Inf;
  accepted_divergences = Inf;
  curr_try = 1;

  while(curr_try <= ntries) {
    assign("curr_fit",tryCatch.W.E(stan(file = stan_path,
                                        fit = stan_fit,
                                        data = list(n_stan = length(y),
                                                    p_stan = ncol(x_standardized),
                                                    y_stan = y,
                                                    x_standardized_stan = x_standardized,
                                                    dof_stan = dof,
                                                    beta_scale_stan = beta_scale,
                                                    slab_precision_stan = slab_precision,
                                                    only_prior = as.integer(only_prior)), 
                                        warmup = mc_warmup, 
                                        iter = mc_iter_after_warmup + mc_warmup, 
                                        chains = mc_chains, 
                                        thin = mc_thin,
                                        control = list(stepsize = mc_stepsize,
                                                       adapt_delta = mc_adapt_delta,
                                                       max_treedepth = mc_max_treedepth)))); 
    if("simpleError"%in%class(curr_fit$value) || "error"%in%class(curr_fit$value)) {
      stop(curr_fit$value);
    }
    if(!"stanfit"%in%class(stan_fit)) {
      break;
    }
    divergent_check = unlist(lapply(curr_fit$warning,grep,pattern="divergent transitions",value=T));
    rhat_check = max(summary(curr_fit$value)$summary[,"Rhat"],na.rm=T);
    #Originally, the break conditions were baesd upon having both no divergent transitions as well as a max Rhat (i.e. gelman-rubin 
    #diagnostic) sufficiently close to 1. I subsequently changed the conditions to be based only upon the first, which is reflected
    #by setting rhat = T immediately below. 
    break_conditions = c(divergence = F, rhat = T);
    if(length(divergent_check) == 0) {#corresponds to zero divergent transitions
      curr_divergences = 0;
      max_divergences = max(max_divergences,curr_divergences,na.rm=T);
      break_conditions["divergence"] = T;
    } else {#corresponds to > zero divergent transitions
      curr_divergences <- max(as.numeric(strsplit(divergent_check," ")$message),na.rm=T);
      max_divergences = max(max_divergences,curr_divergences,na.rm=T);
      curr_try = curr_try + 1;
    }
    #update if fewer divergent transitions were found
    if(curr_divergences < accepted_divergences) {
      accepted_divergences = curr_divergences;
      max_rhat = rhat_check;
      foo = rstan::extract(curr_fit$value);
      curr_beta0 = as.numeric(foo$mu);
      curr_beta = foo$beta;
      theta = foo$theta;
    }
    if(all(break_conditions)) {
      break;
    }
  }
  if(!"stanfit"%in%class(stan_fit)) {
    curr_fit$value;
  } else {
    list(accepted_divergences = accepted_divergences,
         max_divergences = max_divergences,
         max_rhat = max_rhat,
         curr_beta0 = curr_beta0,
         curr_beta = curr_beta,
         theta = theta);
  }
}
