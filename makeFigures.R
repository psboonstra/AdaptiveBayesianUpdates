#DESCRIPTION: run this script after completing both stages of simulations in 'runABUSims.R'. This will create
#the figures reported in the simulation study. Note that I don't know where you've saved your simulation results, 
#i.e. the 1920 .RData files that were created, so you will need to provide that path if different from your current
#working directory. 

rm(list=ls(all=TRUE));
require(rstan);
require(lattice);require(RColorBrewer);require(data.table);
require(plyr);require(dplyr);require(tidyr);
phils_computer = T;
if(phils_computer) {
  read_path = "/Users/philb/Desktop/Work/Barbaro/PedRescuersV2/MethodsPaper/Biostatistics/revisionSims/";
  write_path = "/Users/philb/Desktop/Work/Barbaro/PedRescuersV2/MethodsPaper/Biostatistics/";
} else {
  read_path = 
    write_path = "";
}

sim_set = 1:1920;
file_name = "Sim";
max_run_time = -Inf;

#NOTE: It will appear that R is nonresponsive upon running this for-loop. It hasn't, but apparently 
#loadingthis many binary data files is fairly compute-intensive. R just needs a minute to catch up. 
for(curr_sim in sim_set) {
  foo = try(load(file=paste(read_path,file_name,curr_sim,".RData",sep="")),silent=T);
  if(class(foo)=="try-error") {
    sim_set = setdiff(sim_set,curr_sim);
  } else if ((get(foo)$sim_params$mins_total_runtime) > max_run_time) {
    print(curr_sim);
    max_run_time = max(max_run_time,(get(foo)$sim_params$mins_total_runtime));
  } 
}

sims_to_look_at = sim_set;

sim_params_names = "sim_number";
generating_params_names = c("n_hist","n_curr","beta_label");
covariate_args_names = c("x_correlation","x_orig_binom","x_aug_binom");
population_params_names = c("signal_comp","signal_miss","historical_signal_obs","current_signal_obs","mean_pop_risk");
meth_names = colnames(get(paste0("sim",sim_set[1]))$mse_beta);
meth_names = setdiff(meth_names,"Benchmark");

meth_names_toprint =  c(
  "Historical" = "Historical", 
  "Standard" = "Standard",
  "NABAgnostic" = "NAB_agn",
  "NABOptimist" = "NAB_opt",
  "NAB_devAgnostic" = "NAB_agn_eta2",
  "NAB_devOptimist" = "NAB_opt_eta2",
  "SAB.phiAgnostic.imp1" = "SAB_agn",
  "SAB.phiOptimist.imp1" = "SAB_opt",
  "SAB.phiAgnostic.imp2" = "SAB_agn_1imp",
  "SAB.phiOptimist.imp2" = "SAB_opt_1imp",
  "SAB_dev.phiAgnostic.imp1" = "SAB_agn_eta2",
  "SAB_dev.phiOptimist.imp1" = "SAB_opt_eta2",
  "SAB_dev.phiAgnostic.imp2" = "SAB_agn_eta2_1imp",
  "SAB_dev.phiOptimist.imp2" = "SAB_opt_eta2_1imp"
);

#I had to hard code the pretty names; this is a check to make sure that 
#I didn't change the method names and forget to change the corresponding pretty names
stopifnot(all(names(meth_names_toprint)==meth_names));
#Check to make sure that I've ordered the method names in accrodance with the way they appear in a result
stopifnot(all(names(meth_names_toprint)==setdiff(colnames(get(paste0("sim",sim_set[1]))$mse_beta),"Benchmark")));

individual_results = NULL;

curr_scenario = NA;
array_id = sim_set[1];
for(array_id in sim_set) {
  curr_n_sim = get(paste0("sim",array_id))$sim_params$n_sim;
  
  foo = get(paste0("sim",array_id));
  foo$covariate_args["x_orig_binom"] = ifelse(is.na(foo$covariate_args["x_orig_binom"]), F, T);
  foo$covariate_args["x_aug_binom"] = ifelse(is.na(foo$covariate_args["x_aug_binom"]), F, T);
  curr_params = c(foo$generating_params[c(generating_params_names)],foo$covariate_args[covariate_args_names]);
  curr_scenario = NA;
  curr_dataset_ids = curr_n_sim*(array_id-1) + 1:curr_n_sim;
  curr_p = length(foo$generating_params$true_betas_orig);
  curr_q = length(foo$generating_params$true_betas_aug);
  
  orig_covariates = paste0("orig",1:curr_p);
  aug_covariates = paste0("aug",1:curr_q);
  attach(foo);
  curr_mape = matrix(NA, nrow = curr_n_sim, ncol = length(meth_names),dimnames = list(NULL,paste0("mape_._",meth_names)));
  curr_rmspe = matrix(NA, nrow = curr_n_sim, ncol = length(meth_names),dimnames = list(NULL,paste0("rmspe_._",meth_names)));
  curr_mse_beta = matrix(NA, nrow = curr_n_sim, ncol = length(meth_names),dimnames = list(NULL,paste0("mse_beta_._",meth_names)));
  curr_mse_beta_orig = matrix(NA, nrow = curr_n_sim, ncol = length(meth_names),dimnames = list(NULL,paste0("mse_beta_orig_._",meth_names)));
  curr_mse_beta_aug = matrix(NA, nrow = curr_n_sim, ncol = length(meth_names),dimnames = list(NULL,paste0("mse_beta_aug_._",meth_names)));
  curr_divergences = matrix(NA, nrow = curr_n_sim, ncol = length(meth_names),dimnames = list(NULL,paste0("accepted_divergences_by_method_._",meth_names)));
  curr_phi = matrix(NA, nrow = curr_n_sim, ncol = length(setdiff(meth_names,c("Historical","Standard"))),dimnames = list(NULL,paste0("phi_._",setdiff(meth_names,c("Historical","Standard")))));
  curr_runtime = matrix(NA, nrow = curr_n_sim, ncol = length(meth_names),dimnames = list(NULL,paste0("runtime_._",meth_names)));
  #
  curr_mape[,paste0("mape_._",meth_names)] = mape[,meth_names];
  curr_rmspe[,paste0("rmspe_._",meth_names)] = rmspe[,meth_names];
  curr_mse_beta[,paste0("mse_beta_._",meth_names)] = mse_beta[,meth_names];
  curr_mse_beta_orig[,paste0("mse_beta_orig_._",meth_names)] = mse_beta_orig[,meth_names];
  curr_mse_beta_aug[,paste0("mse_beta_aug_._",meth_names)] = mse_beta_aug[,meth_names];
  curr_divergences[,paste0("accepted_divergences_by_method_._",meth_names)] = accepted_divergences_by_method[,meth_names];
  curr_phi[,paste0("phi_._",setdiff(meth_names,c("Historical","Standard")))] = phi_mean[,setdiff(meth_names,c("Historical","Standard"))];
  curr_runtime[,paste0("runtime_._",meth_names)] = sim_params$mins_per_method_per_sim[,meth_names];
  
  #
  individual_results = rbind(individual_results,
                             cbind(array_id,
                                   curr_dataset_ids,
                                   curr_scenario,
                                   matrix(rep(c(unlist(sim_params[sim_params_names]),
                                                curr_p,
                                                curr_q,
                                                unlist(generating_params[generating_params_names]),
                                                unlist(covariate_args[covariate_args_names])),
                                              each=curr_n_sim),
                                          nrow=curr_n_sim),
                                   population_params[,population_params_names,drop=F],
                                   curr_mape,
                                   curr_rmspe,
                                   curr_mse_beta,
                                   curr_mse_beta_orig,
                                   curr_mse_beta_aug,
                                   curr_divergences,
                                   curr_phi,
                                   curr_runtime));
  
  detach(foo);
  old_params = curr_params;
}
colnames(individual_results) =
  c("array_id",
    "dataset_id",
    "scenario",
    sim_params_names,
    "length_orig",
    "length_aug",
    generating_params_names,
    covariate_args_names,
    population_params_names,
    paste0("mape_._",meth_names_toprint),
    paste0("rmspe_._",meth_names_toprint),
    paste0("mse_beta_._",meth_names_toprint),
    paste0("mse_beta_orig_._",meth_names_toprint),
    paste0("mse_beta_aug_._",meth_names_toprint),
    paste0("accepted_divergences_by_method_._",meth_names_toprint),
    paste0("phi_._",meth_names_toprint[setdiff(meth_names,c("Historical","Standard"))]),
    paste0("runtime_._",meth_names_toprint));
individual_results = data.frame(individual_results);
individual_results$percent_signal_missing = individual_results$signal_miss/individual_results$signal_comp;
individual_results$beta_group = cut(individual_results$beta_label, breaks = c(0,1,2,3,4,8,9,10),right = T,labels = F)
individual_results$label_aug = factor(individual_results$length_aug);
individual_results$label_orig = factor(individual_results$length_orig);
individual_results$pretty_beta_label = factor(individual_results$beta_label,
                                              levels = sort(unique(individual_results$beta_label)),
                                              labels = unique(eval(bquote(paste("b[",individual_results$beta_label,"]:group('{',list(p,q),'}') == group('{',list(",individual_results$length_orig,",",individual_results$length_aug,"),'}')",sep=""))))[order(unique(individual_results$beta_label))],
                                              ordered = T);
individual_results$pretty_n_curr_label = factor(eval(bquote(paste("n[curr] ==", individual_results$n_curr))),ordered = T);

#
individual_results = 
  individual_results %>% 
  gather(v, value, contains("_._")) %>% 
  separate(v, c("Metric", "Method"),sep="_._") %>% 
  arrange(beta_label,x_correlation,x_orig_binom,n_hist,n_curr,array_id,dataset_id) %>% 
  spread(Method,value);

foo = which(!duplicated(individual_results[,c("beta_label","x_correlation","x_orig_binom","n_hist","n_curr")]));
foo2 = as.numeric(rep(unlist(lapply(strsplit(as.character(individual_results[foo,"sim_number"]), "[.]"),"[[",1)), diff(c(foo,nrow(individual_results)+1))))
individual_results$scenario = foo2;
rm(foo,foo2);
#
#Make ratios before collapsing over scenarios
all_comb = expand.grid(as.character(meth_names_toprint),
                       as.character(meth_names_toprint),stringsAsFactors = F)
foo = individual_results[,all_comb[,1]] / individual_results[,all_comb[,2]];
colnames(foo) = paste0(all_comb[,1],"_over_",all_comb[,2]);
individual_results = cbind(individual_results,foo);
rm(foo,all_comb);
#Collapse by scenario
columns_to_summarize = setdiff(colnames(individual_results), c("array_id","dataset_id","scenario","sim_number","Metric","label_orig","label_aug","pretty_beta_label","pretty_n_curr_label"))
scenario_results =
  individual_results %>%
  group_by(scenario,Metric,label_orig,label_aug,pretty_beta_label) %>%
  summarize_at(columns_to_summarize,mean,na.rm=T)
######################################################

######################################################
figure1_results = 
  individual_results %>%
  filter(Metric == "mse_beta") %>%
  dplyr::select("Metric","scenario","n_hist","n_curr","beta_label","beta_group","x_correlation","x_orig_binom","label_orig","label_aug","pretty_beta_label","pretty_n_curr_label","percent_signal_missing",
                "SAB_agn_over_Standard",
                "SAB_opt_over_Standard",
                "NAB_agn_over_Standard",
                "NAB_opt_over_Standard") %>% 
  reshape2::melt(id.vars=c("Metric","scenario","n_hist","n_curr","beta_label","beta_group","x_correlation","x_orig_binom","label_orig","label_aug","pretty_beta_label","pretty_n_curr_label","percent_signal_missing"))
figure1_results$variable = factor(figure1_results$variable,
                                  labels = c(
                                    "SAB(agnostic)/Standard",
                                    "SAB(optimist)/Standard",
                                    "NAB(agnostic)/Standard",
                                    "NAB(optimist)/Standard"));
#

ggplot(filter(figure1_results,beta_label < 6, x_orig_binom == T), 
       aes(x = factor(n_hist), 
           y = sqrt(value))) + 
  geom_hline(yintercept = 1, linetype = 2, color = "#303030") + 
  geom_boxplot(aes(group = interaction(n_hist,variable),
                   color = factor(variable),
                   fill = factor(variable)),
               position = position_dodge(width =  0.85), 
               alpha = 1,
               outlier.size = 0.5,
               size = 0.2) + 
  facet_grid(pretty_n_curr_label~pretty_beta_label, labeller = label_parsed) +
  scale_fill_grey() + 
  scale_color_grey(start = 0, end = 0) +  
  scale_y_continuous(trans = "log2",
                     breaks = 2^seq(-3, 1, by = 0.5),
                     minor_breaks = 2^seq(-3, 1, by = 0.25),
                     labels = function(x) {formatC(x,digits=2,format="f");}) + 
  labs(x = expression(n[hist]), 
       y = "RMSE Ratio",
       fill = "",
       color = "") + 
  theme(axis.text = element_text(colour = 'black', angle = 0, size = 12, hjust = 0.5, vjust = 0.5),
        axis.title= element_text(size=12),
        strip.text = element_text(size = 12),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.position = "top",
        legend.text.align = 0);
ggsave(paste0(write_path,"/fig1a.eps"),height = 5., width = 9);
#
ggplot(filter(figure1_results,beta_label %in% c(6:10), x_orig_binom == T), 
       aes(x = factor(n_hist), 
           y = sqrt(value))) + 
  geom_hline(yintercept = 1, linetype = 2, color = "#303030") + 
  geom_boxplot(aes(group = interaction(n_hist,variable),
                   color = factor(variable),
                   fill = factor(variable)),
               position = position_dodge(width =  0.85), 
               alpha = 1,
               outlier.size = 0.5,
               size = 0.2) + 
  facet_grid(pretty_n_curr_label~pretty_beta_label, labeller = label_parsed) +
  scale_fill_grey() + 
  scale_color_grey(start = 0, end = 0) +  
  scale_y_continuous(trans = "log2",
                     breaks = 2^seq(-4,0,by=1), 
                     minor_breaks = 2^seq(-4, 0.5, by = 0.5),
                     labels = function(x) {formatC(x,digits=2,format="f");}) + 
  labs(x = expression(n[hist]), 
       y = "RMSE Ratio",
       fill = "",
       color = "") +
  guides(color = F, fill = F) + 
  theme(axis.text = element_text(colour = 'black', angle = 0, size = 12, hjust = 0.5, vjust = 0.5),
        axis.title= element_text(size=12),
        strip.text = element_text(size = 10),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.position = "top",
        legend.text.align = 0);
ggsave(paste0(write_path,"/fig1b.eps"),height = 4.75, width = 9);


summary(filter(figure1_results,beta_label <= 5, n_hist == 100, n_curr == 100, x_orig_binom == T, variable %in% c("SAB(agnostic)/Standard"))[,"value"]);
summary(filter(figure1_results,beta_label <= 5, n_hist == 1600, n_curr == 100, x_orig_binom == T, variable %in% c("SAB(agnostic)/Standard"))[,"value"]);
#
summary(filter(figure1_results,beta_label == 5, n_hist == 100, n_curr == 100, x_orig_binom == T, variable %in% c("NAB(agnostic)/Standard"))[,"value"]);
summary(filter(figure1_results,beta_label == 5, n_hist == 1600, n_curr == 100, x_orig_binom == T, variable %in% c("NAB(agnostic)/Standard"))[,"value"]);
#
summary(filter(figure1_results,beta_label < 6, x_orig_binom == T, variable %in% c("NAB(agnostic)/Standard"))[,"value"]);
summary(filter(figure1_results,beta_label < 6, x_orig_binom == T, variable %in% c("SAB(agnostic)/Standard"))[,"value"]);
summary(filter(figure1_results,beta_label >= 6, x_orig_binom == T, variable %in% c("NAB(agnostic)/Standard"))[,"value"]);
summary(filter(figure1_results,beta_label >= 6, x_orig_binom == T, variable %in% c("SAB(agnostic)/Standard"))[,"value"]);
######################################################


######################################################
figureS1_results = 
  individual_results %>%
  filter(Metric == "mse_beta_orig") %>%
  dplyr::select("Metric","scenario","n_hist","n_curr","beta_label","beta_group","x_correlation","x_orig_binom","label_orig","label_aug","pretty_beta_label","pretty_n_curr_label","percent_signal_missing",
                "SAB_agn_over_Standard",
                "SAB_opt_over_Standard",
                "NAB_agn_over_Standard",
                "NAB_opt_over_Standard") %>% 
  reshape2::melt(id.vars=c("Metric","scenario","n_hist","n_curr","beta_label","beta_group","x_correlation","x_orig_binom","label_orig","label_aug","pretty_beta_label","pretty_n_curr_label","percent_signal_missing"))
figureS1_results$variable = factor(figureS1_results$variable,
                                   labels = c(
                                     "SAB(agnostic)/Standard",
                                     "SAB(optimist)/Standard",
                                     "NAB(agnostic)/Standard",
                                     "NAB(optimist)/Standard"));
#

ggplot(filter(figureS1_results,beta_label < 6, x_orig_binom == T),
       aes(x = factor(n_hist), 
           y = sqrt(value))) + 
  geom_hline(yintercept = 1, linetype = 2, color = "#303030") + 
  geom_boxplot(aes(group = interaction(n_hist,variable),
                   color = factor(variable),
                   fill = factor(variable)),
               position = position_dodge(width =  0.85), 
               alpha = 1,
               outlier.size = 0.5,
               size = 0.2) + 
  facet_grid(pretty_n_curr_label~pretty_beta_label, labeller = label_parsed) +
  #scale_fill_grey() + 
  scale_color_grey(start = 0, end = 0) +  
  scale_y_continuous(trans = "log2",
                     breaks = 2^seq(-4, 1, by = 0.5),
                     minor_breaks = 2^seq(-4, 1, by = 0.25),
                     labels = function(x) {formatC(x,digits=2,format="f");}) + 
  labs(x = expression(n[hist]), 
       y = "RMSE Ratio (Original)",
       fill = "",
       color = "") + 
  theme(axis.text = element_text(colour = 'black', angle = 0, size = 12, hjust = 0.5, vjust = 0.5),
        axis.title= element_text(size=12),
        strip.text = element_text(size = 12),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.position = "top",
        legend.text.align = 0);
ggsave(paste0(write_path,"/figS1a.eps"),height = 5., width = 9);
#
ggplot(filter(figureS1_results,beta_label %in% c(6:10), x_orig_binom == T),
       aes(x = factor(n_hist), 
           y = sqrt(value))) + 
  geom_hline(yintercept = 1, linetype = 2, color = "#303030") + 
  geom_boxplot(aes(group = interaction(n_hist,variable),
                   color = factor(variable),
                   fill = factor(variable)),
               position = position_dodge(width =  0.85), 
               alpha = 1,
               outlier.size = 0.5,
               size = 0.2) + 
  facet_grid(pretty_n_curr_label~pretty_beta_label, labeller = label_parsed) +
  #scale_fill_grey() + 
  scale_color_grey(start = 0, end = 0) +  
  scale_y_continuous(trans = "log2",
                     breaks = 2^seq(-5,0,by=1), 
                     minor_breaks = 2^seq(-5, 0.5, by = 0.5),
                     labels = function(x) {formatC(x,digits=2,format="f");}) + 
  labs(x = expression(n[hist]), 
       y = "RMSE Ratio (Original)",
       fill = "",
       color = "") +
  guides(color = F, fill = F) + 
  theme(axis.text = element_text(colour = 'black', angle = 0, size = 12, hjust = 0.5, vjust = 0.5),
        axis.title= element_text(size=12),
        strip.text = element_text(size = 10),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.position = "top",
        legend.text.align = 0);
ggsave(paste0(write_path,"/figS1b.eps"),height = 4.75, width = 9);
######################################################

######################################################
figureS2_results = 
  individual_results %>%
  filter(Metric == "mse_beta_aug") %>%
  dplyr::select("Metric","scenario","n_hist","n_curr","beta_label","beta_group","x_correlation","x_orig_binom","label_orig","label_aug","pretty_beta_label","pretty_n_curr_label","percent_signal_missing",
                "SAB_agn_over_Standard",
                "SAB_opt_over_Standard",
                "NAB_agn_over_Standard",
                "NAB_opt_over_Standard") %>% 
  reshape2::melt(id.vars=c("Metric","scenario","n_hist","n_curr","beta_label","beta_group","x_correlation","x_orig_binom","label_orig","label_aug","pretty_beta_label","pretty_n_curr_label","percent_signal_missing"))
figureS2_results$variable = factor(figureS2_results$variable,
                                   labels = c(
                                     "SAB(agnostic)/Standard",
                                     "SAB(optimist)/Standard",
                                     "NAB(agnostic)/Standard",
                                     "NAB(optimist)/Standard"));
#

ggplot(filter(figureS2_results,beta_label < 6, x_orig_binom == T),
       aes(x = factor(n_hist), 
           y = sqrt(value))) + 
  geom_hline(yintercept = 1, linetype = 2, color = "#303030") + 
  geom_boxplot(aes(group = interaction(n_hist,variable),
                   color = factor(variable),
                   fill = factor(variable)),
               position = position_dodge(width =  0.85), 
               alpha = 1,
               outlier.size = 0.5,
               size = 0.2) + 
  facet_grid(pretty_n_curr_label~pretty_beta_label, labeller = label_parsed) +
  #scale_fill_grey() + 
  scale_color_grey(start = 0, end = 0) +  
  scale_y_continuous(trans = "log2",
                     breaks = 2^seq(-3, 1, by = 0.5),
                     minor_breaks = 2^seq(-3, 1, by = 0.25),
                     labels = function(x) {formatC(x,digits=2,format="f");}) + 
  labs(x = expression(n[hist]), 
       y = "RMSE Ratio (Added)",
       fill = "",
       color = "") + 
  theme(axis.text = element_text(colour = 'black', angle = 0, size = 12, hjust = 0.5, vjust = 0.5),
        axis.title= element_text(size=12),
        strip.text = element_text(size = 12),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.position = "top",
        legend.text.align = 0);
ggsave(paste0(write_path,"/figS2a.eps"),height = 5., width = 9);
#
ggplot(filter(figureS2_results,beta_label %in% c(6:10), x_orig_binom == T),
       aes(x = factor(n_hist), 
           y = sqrt(value))) + 
  geom_hline(yintercept = 1, linetype = 2, color = "#303030") + 
  geom_boxplot(aes(group = interaction(n_hist,variable),
                   color = factor(variable),
                   fill = factor(variable)),
               position = position_dodge(width =  0.85), 
               alpha = 1,
               outlier.size = 0.5,
               size = 0.2) + 
  facet_grid(pretty_n_curr_label~pretty_beta_label, labeller = label_parsed) +
  #scale_fill_grey() + 
  scale_color_grey(start = 0, end = 0) +  
  scale_y_continuous(trans = "log2",
                     breaks = 2^seq(-4,0,by=1), 
                     minor_breaks = 2^seq(-4, 0.5, by = 0.5),
                     labels = function(x) {formatC(x,digits=2,format="f");}) + 
  labs(x = expression(n[hist]), 
       y = "RMSE Ratio (Added)",
       fill = "",
       color = "") +
  guides(color = F, fill = F) + 
  theme(axis.text = element_text(colour = 'black', angle = 0, size = 12, hjust = 0.5, vjust = 0.5),
        axis.title= element_text(size=12),
        strip.text = element_text(size = 10),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.position = "top",
        legend.text.align = 0);
ggsave(paste0(write_path,"/figS2b.eps"),height = 4.75, width = 9);
######################################################

######################################################
figureS3_results = 
  individual_results %>%
  filter(Metric == "mse_beta") %>%
  dplyr::select("Metric","scenario","n_hist","n_curr","beta_label","beta_group","x_correlation","x_orig_binom","label_orig","label_aug","pretty_beta_label","pretty_n_curr_label","percent_signal_missing",
                "SAB_agn_over_SAB_agn_eta2",
                "SAB_opt_over_SAB_opt_eta2",
                "NAB_agn_over_NAB_agn_eta2",
                "NAB_opt_over_NAB_opt_eta2") %>% 
  reshape2::melt(id.vars=c("Metric","scenario","n_hist","n_curr","beta_label","beta_group","x_correlation","x_orig_binom","label_orig","label_aug","pretty_beta_label","pretty_n_curr_label","percent_signal_missing"))
figureS3_results$variable = factor(figureS3_results$variable,
                                   labels = c(
                                     "SAB(agnostic)",
                                     "SAB(optimist)",
                                     "NAB(agnostic)",
                                     "NAB(optimist)"));
#

ggplot(filter(figureS3_results,beta_label < 6, x_orig_binom == T), 
       aes(x = factor(n_hist), 
           y = sqrt(value))) + 
  geom_hline(yintercept = 1, linetype = 2, color = "#303030") + 
  geom_boxplot(aes(group = interaction(n_hist,variable),
                   color = factor(variable),
                   fill = factor(variable)),
               position = position_dodge(width =  0.85), 
               alpha = 1,
               outlier.size = 0.5,
               size = 0.2) + 
  facet_grid(pretty_n_curr_label~pretty_beta_label, labeller = label_parsed) +
  #scale_fill_grey() + 
  scale_color_grey(start = 0, end = 0) +  
  scale_y_continuous(trans = "log2",
                     breaks = 2^seq(-0.5, 0.5, by = 0.5/4),
                     minor_breaks = 2^seq(-0.5, 0.5, by = 0.25/4),
                     labels = function(x) {formatC(x,digits=2,format="f");}) + 
  labs(x = expression(n[hist]), 
       y = "RMSE Ratio",
       fill = "",
       color = "") + 
  theme(axis.text = element_text(colour = 'black', angle = 0, size = 12, hjust = 0.5, vjust = 0.5),
        axis.title= element_text(size=12),
        strip.text = element_text(size = 12),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.position = "top",
        legend.text.align = 0);
ggsave(paste0(write_path,"/figS3a.eps"),height = 5., width = 9);
#
ggplot(filter(figureS3_results,beta_label %in% c(6:10), x_orig_binom == T), 
       aes(x = factor(n_hist), 
           y = sqrt(value))) + 
  geom_hline(yintercept = 1, linetype = 2, color = "#303030") + 
  geom_boxplot(aes(group = interaction(n_hist,variable),
                   color = factor(variable),
                   fill = factor(variable)),
               position = position_dodge(width =  0.85), 
               alpha = 1,
               outlier.size = 0.5,
               size = 0.2) + 
  facet_grid(pretty_n_curr_label~pretty_beta_label, labeller = label_parsed) +
  #scale_fill_grey() + 
  scale_color_grey(start = 0, end = 0) +  
  scale_y_continuous(trans = "log2",
                     breaks = 2^seq(-0.75, 0.75, by = 0.5/3),
                     minor_breaks = 2^seq(-0.75, 0.75, by = 0.25/3),
                     labels = function(x) {formatC(x,digits=2,format="f");}) + 
  labs(x = expression(n[hist]), 
       y = "RMSE Ratio",
       fill = "",
       color = "") +
  guides(color = F, fill = F) + 
  theme(axis.text = element_text(colour = 'black', angle = 0, size = 12, hjust = 0.5, vjust = 0.5),
        axis.title= element_text(size=12),
        strip.text = element_text(size = 10),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.position = "top",
        legend.text.align = 0);
ggsave(paste0(write_path,"/figS3b.eps"),height = 4.75, width = 9);


summary(filter(figureS3_results,beta_label < 6, x_orig_binom == T, variable %in% c("NAB(agnostic)"))[,"value"])
summary(filter(figureS3_results,beta_label < 6, x_orig_binom == T, variable %in% c("SAB(agnostic)"))[,"value"])
summary(filter(figureS3_results,beta_label >= 6, x_orig_binom == T, variable %in% c("NAB(agnostic)"))[,"value"])
summary(filter(figureS3_results,beta_label >= 6, x_orig_binom == T, variable %in% c("SAB(agnostic)"))[,"value"])
######################################################

######################################################
runtime_summaries = 
  individual_results %>%
  filter(Metric == "runtime", x_orig_binom == 1) %>%
  dplyr::select("Metric","scenario","n_hist","n_curr","beta_label","beta_group","label_orig","label_aug","pretty_beta_label","pretty_n_curr_label","percent_signal_missing",
                "Standard",
                "NAB_agn",
                "NAB_opt",
                "SAB_agn",
                "SAB_opt",
                "SAB_agn_over_Standard",
                "SAB_opt_over_Standard",
                "NAB_agn_over_Standard",
                "NAB_opt_over_Standard",
                "NAB_agn_over_NAB_opt",
                "SAB_agn_over_SAB_opt") %>% 
  reshape2::melt(id.vars=c("Metric","scenario","n_hist","n_curr","beta_label","beta_group","label_orig","label_aug","pretty_beta_label","pretty_n_curr_label","percent_signal_missing"))
runtime_summaries = 
  runtime_summaries %>% 
  group_by(scenario,variable,n_hist,n_curr,beta_label) %>% 
  dplyr::summarize(median_runtime = median(value));
runtime_summaries$variable = factor(runtime_summaries$variable,
                                    labels = c(
                                      "Standard",
                                      "NAB(agn)",
                                      "NAB(opt)",
                                      "SAB(agn)",
                                      "SAB(opt)",
                                      "NAB(agn)/Standard",
                                      "NAB(opt)/Standard",
                                      "SAB(agn)/Standard",
                                      "SAB(opt)/Standard",
                                      "NAB(agn)/NAB(opt)",
                                      "SAB(agn)/SAB(opt)"));


filter(runtime_summaries, variable == "Standard") %>%
  arrange(median_runtime);
filter(runtime_summaries, variable == "Standard") %>%
  arrange(desc(median_runtime));
filter(runtime_summaries, variable == "NAB(agn)") %>%
  arrange(median_runtime)
filter(runtime_summaries, variable == "NAB(agn)") %>%
  arrange(desc(median_runtime))
filter(runtime_summaries, variable == "SAB(agn)") %>%
  arrange(median_runtime)
filter(runtime_summaries, variable == "SAB(agn)") %>%
  arrange(desc(median_runtime))


file_name = paste0(write_path,"median_runtimes.txt");
methods_to_show = c("Standard",
                    "NAB(agn)",
                    "NAB(opt)",
                    "SAB(agn)",
                    "SAB(opt)");
write.table(matrix(c("Coef. Set",methods_to_show),nrow=1),file=file_name,row.names=F,col.names=F,sep="&",na="",quote=F,eol="\\\\\\hline\n",append=F)

table_vals = NULL;
for(curr_beta in 1:10) {
  foo = filter(runtime_summaries, beta_label == curr_beta, variable %in% methods_to_show) %>%
    group_by(variable) %>%
    dplyr::summarize(shortest = min(median_runtime),
                     longest = max(median_runtime));
  curr_vals = c(paste0("$b_{",curr_beta,"}$"),
                paste0("(",formatC(pull(foo,shortest),format="f",digits = 2),",",
                       formatC(pull(foo,longest),format="f",digits = 2),")"));
   
  table_vals = rbind(table_vals,curr_vals);
  if(curr_beta == "10") {
    table_vals[nrow(table_vals),ncol(table_vals)] = paste(table_vals[nrow(table_vals),ncol(table_vals)]," \\\\ \\hline ",sep="")
  } else {
    table_vals[nrow(table_vals),ncol(table_vals)] = paste(table_vals[nrow(table_vals),ncol(table_vals)]," \\\\ ",sep="")
  }
  rm(foo);
}
write.table(table_vals,file=file_name,row.names=F,col.names=F,sep="&",na="",quote=F,eol="\n",append=T);
######################################################


######################################################
figureA1_results = 
  individual_results %>%
  filter(Metric == "mse_beta") %>%
  dplyr::select("Metric","scenario","n_hist","n_curr","beta_label","beta_group","x_correlation","x_orig_binom","label_orig","label_aug","pretty_beta_label","pretty_n_curr_label","percent_signal_missing",
                "SAB_agn_over_Standard",
                "SAB_agn_1imp_over_Standard") %>% 
  reshape2::melt(id.vars=c("Metric","scenario","n_hist","n_curr","beta_label","beta_group","x_correlation","x_orig_binom","label_orig","label_aug","pretty_beta_label","pretty_n_curr_label","percent_signal_missing"))
figureA1_results$variable = factor(figureA1_results$variable,
                                   labels = c(
                                     "SAB(agn)/Standard",
                                     "SAB1Imp(agn)/Standard"));
#
ggplot(filter(figureA1_results,beta_label < 6, x_orig_binom == T), 
       aes(x = factor(n_hist), 
           y = sqrt(value))) + 
  geom_hline(yintercept = 1, linetype = 2, color = "#303030") + 
  geom_boxplot(aes(group = interaction(n_hist,variable),
                   color = factor(variable),
                   fill = factor(variable)),
               position = position_dodge(width =  0.85), 
               alpha = 1,
               outlier.size = 0.5,
               size = 0.2) + 
  facet_grid(pretty_n_curr_label~pretty_beta_label, labeller = label_parsed) +
  #scale_fill_grey() + 
  scale_color_grey(start = 0, end = 0) +  
  scale_y_continuous(trans = "log2",
                     breaks = 2^seq(-3, 1, by = 0.5),
                     minor_breaks = 2^seq(-3, 1, by = 0.25),
                     labels = function(x) {formatC(x,digits=2,format="f");}) + 
  labs(x = expression(n[hist]), 
       y = "RMSE Ratio",
       fill = "",
       color = "") + 
  theme(axis.text = element_text(colour = 'black', angle = 0, size = 12, hjust = 0.5, vjust = 0.5),
        axis.title= element_text(size=12),
        strip.text = element_text(size = 12),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.position = "top",
        legend.text.align = 0);
ggsave(paste0(write_path,"/figA1a.eps"),height = 5., width = 9);
#
ggplot(filter(figureA1_results,beta_label %in% c(6:10), x_orig_binom == T), 
       aes(x = factor(n_hist), 
           y = sqrt(value))) + 
  geom_hline(yintercept = 1, linetype = 2, color = "#303030") + 
  geom_boxplot(aes(group = interaction(n_hist,variable),
                   color = factor(variable),
                   fill = factor(variable)),
               position = position_dodge(width =  0.85), 
               alpha = 1,
               outlier.size = 0.5,
               size = 0.2) + 
  facet_grid(pretty_n_curr_label~pretty_beta_label, labeller = label_parsed) +
  #scale_fill_grey() + 
  scale_color_grey(start = 0, end = 0) +  
  scale_y_continuous(trans = "log2",
                     breaks = 2^seq(-4,0,by=1), 
                     minor_breaks = 2^seq(-4, 0.5, by = 0.5),
                     labels = function(x) {formatC(x,digits=2,format="f");}) + 
  labs(x = expression(n[hist]), 
       y = "RMSE Ratio",
       fill = "",
       color = "") +
  guides(color = F, fill = F) + 
  theme(axis.text = element_text(colour = 'black', angle = 0, size = 12, hjust = 0.5, vjust = 0.5),
        axis.title= element_text(size=12),
        strip.text = element_text(size = 10),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.position = "top",
        legend.text.align = 0);
ggsave(paste0(write_path,"/figA1b.eps"),height = 4.75, width = 9);

######################################################
figureA2_results = 
  individual_results %>%
  filter(Metric == "phi") %>%
  dplyr::select("Metric","scenario","n_hist","n_curr","beta_label","beta_group","x_correlation","x_orig_binom","label_orig","label_aug","pretty_beta_label","pretty_n_curr_label","percent_signal_missing",
                "SAB_agn",
                "SAB_agn_1imp") %>% 
  reshape2::melt(id.vars=c("Metric","scenario","n_hist","n_curr","beta_label","beta_group","x_correlation","x_orig_binom","label_orig","label_aug","pretty_beta_label","pretty_n_curr_label","percent_signal_missing"))
figureA2_results$variable = factor(figureA2_results$variable,
                                   labels = c(
                                     "SAB(agn)",
                                     "SAB1Imp(agn)"));

ggplot(filter(figureA2_results,beta_label < 6, x_orig_binom == T), 
       aes(x = factor(n_hist), 
           y = sqrt(value))) + 
  geom_hline(yintercept = 1, linetype = 2, color = "#303030") + 
  geom_boxplot(aes(group = interaction(n_hist,variable),
                   color = factor(variable),
                   fill = factor(variable)),
               position = position_dodge(width =  0.85), 
               alpha = 1,
               outlier.size = 0.5,
               size = 0.2) + 
  facet_grid(pretty_n_curr_label~pretty_beta_label, labeller = label_parsed) +
  scale_fill_grey() + 
  scale_color_grey(start = 0, end = 0) +  
  scale_y_continuous(trans = "log2",
                     breaks = 2^seq(-8, 8, by = 2),
                     minor_breaks = 2^seq(-8, 8, by = 1),
                     labels = function(x) {formatC(x,digits=2,format="f");}) + 
  labs(x = expression(n[hist]), 
       y = "RMSE Ratio",
       fill = "",
       color = "") + 
  theme(axis.text = element_text(colour = 'black', angle = 0, size = 12, hjust = 0.5, vjust = 0.5),
        axis.title= element_text(size=12),
        strip.text = element_text(size = 12),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.position = "top",
        legend.text.align = 0);
ggsave(paste0(write_path,"/figA2a.eps"),height = 5., width = 9);
#
ggplot(filter(figureA2_results,beta_label %in% c(6:10), x_orig_binom == T), 
       aes(x = factor(n_hist), 
           y = sqrt(value))) + 
  geom_hline(yintercept = 1, linetype = 2, color = "#303030") + 
  geom_boxplot(aes(group = interaction(n_hist,variable),
                   color = factor(variable),
                   fill = factor(variable)),
               position = position_dodge(width =  0.85), 
               alpha = 1,
               outlier.size = 0.5,
               size = 0.2) + 
  facet_grid(pretty_n_curr_label~pretty_beta_label, labeller = label_parsed) +
  scale_fill_grey() + 
  scale_color_grey(start = 0, end = 0) +  
  scale_y_continuous(trans = "log2",
                     breaks = 2^seq(-8, 8, by = 2),
                     minor_breaks = 2^seq(-8, 8, by = 1),
                     labels = function(x) {formatC(x,digits=2,format="f");}) + 
  labs(x = expression(n[hist]), 
       y = "RMSE Ratio",
       fill = "",
       color = "") +
  guides(color = F, fill = F) + 
  theme(axis.text = element_text(colour = 'black', angle = 0, size = 12, hjust = 0.5, vjust = 0.5),
        axis.title= element_text(size=12),
        strip.text = element_text(size = 10),
        legend.title = element_text(size=12),
        legend.text = element_text(size=12),
        legend.position = "top",
        legend.text.align = 0);
ggsave(paste0(write_path,"/figA2b.eps"),height = 4.75, width = 9);


######################################################
#WNAR plots
wnar_results = filter(figure1_results,
                      beta_label %in% c(1,5,6,7,8,10), 
                      x_orig_binom == T,
                      n_curr == 200, 
                      variable %in% c("SAB(agnostic)/Standard", "NAB(agnostic)/Standard"));
wnar_results$wnar_variable = factor(as.character(wnar_results$variable),labels = c("Naive/Standard","Sensible/Standard"));
wnar_results$wnar_pretty_beta_label = factor(as.character(wnar_results$pretty_beta_label),
                                             levels = intersect(as.character(wnar_results$pretty_beta_label),levels(wnar_results$pretty_beta_label)),
                                             labels = c("b[1]","b[2]","b[3]","b[4]","b[5]","b[6]"));

ggplot(filter(wnar_results, beta_label %in% c(1,5,6)),
       aes(x = factor(n_hist), 
           y = sqrt(value))) + 
  geom_hline(yintercept = 1, linetype = 2, color = "#303030") + 
  geom_boxplot(aes(group = interaction(n_hist,wnar_variable),
                   color = factor(wnar_variable),
                   fill = factor(wnar_variable)),
               position = position_dodge(width =  0.85), 
               alpha = 1,
               outlier.size = 0.5,
               size = 0.2) + 
  facet_grid(~wnar_pretty_beta_label, labeller = label_parsed) +
  scale_color_grey(start = 0, end = 0) +  
  scale_y_continuous(trans = "log2",
                     breaks = 2^seq(-3, 1, by = 0.5),
                     minor_breaks = 2^seq(-3, 1, by = 0.25),
                     labels = function(x) {formatC(x,digits=2,format="f");}) + 
  labs(x = expression(n[hist]), 
       y = "",
       fill = "",
       color = "") + 
  theme(text = element_text(size = 20),
        axis.text = element_text(colour = 'black', angle = 0, size = 20, hjust = 0.5, vjust = 0.5),
        axis.title= element_text(size=20),
        strip.text = element_text(size = 20),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20),
        legend.position = "bottom",
        legend.text.align = 0);
ggsave(paste0(write_path,"/fig1a_WNAR.eps"),height = 5., width = 9);
#
ggplot(filter(wnar_results, beta_label %in% c(7,8,10)),
       aes(x = factor(n_hist), 
           y = sqrt(value))) + 
  geom_hline(yintercept = 1, linetype = 2, color = "#303030") + 
  geom_boxplot(aes(group = interaction(n_hist,wnar_variable),
                   color = factor(wnar_variable),
                   fill = factor(wnar_variable)),
               position = position_dodge(width =  0.85), 
               alpha = 1,
               outlier.size = 0.5,
               size = 0.2) + 
  facet_grid(~wnar_pretty_beta_label, labeller = label_parsed) +
  #scale_fill_grey() + 
  scale_color_grey(start = 0, end = 0) +  
  scale_y_continuous(trans = "log2",
                     breaks = 2^seq(-3, 1, by = 0.5),
                     minor_breaks = 2^seq(-3, 1, by = 0.25),
                     #labels = c(expression(2^-1.5),expression(2^-1),expression(2^-0.5),expression(2^0),expression(2^0.5),expression(2^1)))+
                     labels = function(x) {formatC(x,digits=2,format="f");}) + 
  labs(x = expression(n[hist]), 
       y = "",
       fill = "",
       color = "") + 
  theme(text = element_text(size = 20),
        axis.text = element_text(colour = 'black', angle = 0, size = 20, hjust = 0.5, vjust = 0.5),
        axis.title= element_text(size=20),
        strip.text = element_text(size = 20),
        legend.title = element_text(size=20),
        legend.text = element_text(size=20),
        legend.position = "bottom",
        legend.text.align = 0);
ggsave(paste0(write_path,"/fig1b_WNAR.eps"),height = 5., width = 9);
#
