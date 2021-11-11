######### Create plots for the priors #########

#### 1. Load Bayes factors ##### 

priors <- c("beta_prime", 
            "fixed_g", 
            "hyper_gn", 
            "intrinsic", 
            "robust", 
            "unit_information")

# make sure that this is in the same order as the prior list above!
prior_functions <- list(beta.prime(),
                     g.prior(100),
                     hyper.g.n(),
                     intrinsic(),
                     robust(),
                     "Unit Information")


column_names_large <-  colMeans(structure(vapply(sim_data_u0_large, function(x) sum(x$arr),
                                                 numeric(1)), dim=dim(sim_data_u0_large)))

column_names_small <-  colMeans(structure(vapply(sim_data_u0_small, function(x) sum(x$arr),
                                                 numeric(1)), dim=dim(sim_data_u0_small)))


# read in bayes factors for large n
bayes_factors_nlarge <- lapply(priors, function(x) {
  # read logit-based Bayes factors
  file_name_log <- paste0("results/Bayesian/priors/", 
         x, "/",  
         "logistic_n1500.RData")
  print(file_name_log)
  
  # read poisson-based Bayes factors
  file_name_poi <- paste0("results/Bayesian/priors/", 
                          x, "/", 
                          "poisson_n1500.RData")
  
  list(logistic = readRDS(file_name_log),
       poisson = readRDS(file_name_poi))
})

names(bayes_factors_nlarge) <- priors

# create plots save them in the plots folder
lapply(1:length(bayes_factors_nlarge),function(i) {
  print(priors[i])
  create_bf_plots(input_prior = prior_functions[[i]],
                  input_prior_name = priors[i],
                  input_logistic_data = bayes_factors_nlarge[[i]]$logistic,
                  input_poisson_data = bayes_factors_nlarge[[i]]$poisson,
                  input_column_names = column_names_large)
})


# read in bayes factors for small n
bayes_factors_nsmall <- lapply(priors, function(x) {
  # read logit-based Bayes factors
  file_name_log <- paste0("results/Bayesian/priors/", 
                          x, "/", 
                          "logistic_n100.RData")
  
  # read poisson-based Bayes factors
  file_name_poi <- paste0("results/Bayesian/priors/", 
                          x, "/", 
                          "poisson_n100.RData")
  
  list(logistic = readRDS(file_name_log),
       poisson = readRDS(file_name_poi))
})

names(bayes_factors_nsmall) <- priors


# create plots for small n
lapply(1:length(bayes_factors_nsmall),function(i) {
  print(priors[i])
  create_bf_plots(input_prior = prior_functions[[i]],
                  input_prior_name = priors[i],
                  input_logistic_data = bayes_factors_nsmall[[i]]$logistic,
                  input_poisson_data = bayes_factors_nsmall[[i]]$poisson,
                  input_column_names = column_names_small)
})


