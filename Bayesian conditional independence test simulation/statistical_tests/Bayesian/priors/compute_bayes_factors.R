#### 1. Create Bayes factors ##### 

# priors for the simulation 
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

# in case, the data sets are not yet loaded 
sim_data_u0_large <- readRDS("data/sim_data_u0_large.RData")
sim_data_u1_large <- readRDS("data/sim_data_u1_large.RData")
sim_data_u2_large <- readRDS("data/sim_data_u2_large.RData")
sim_data_u3_large <- readRDS("data/sim_data_u3_large.RData")

# number of columns (sample size specifications), for which  the Bayes factor should be computed 
#sim_data_large <- list(sim_data_u0_large[,1:min(15,ncol(sim_data_u0_large))], 
#                       sim_data_u1_large[,1:min(15,ncol(sim_data_u0_large))])

sim_data_large <- list(sim_data_u0_large, 
                       sim_data_u1_large)

sim_data_u0_small <- readRDS("data/sim_data_u0_small.RData")
sim_data_u1_small <- readRDS("data/sim_data_u1_small.RData")
sim_data_u2_small <- readRDS("data/sim_data_u2_small.RData")
sim_data_u3_small <- readRDS("data/sim_data_u3_small.RData")


sim_data_small<- list(sim_data_u0_small,
                      sim_data_u1_small,
                      sim_data_u2_small,
                      sim_data_u3_small)
                
# set to TRUE if results should be overwritten 
overwrite_bf <- TRUE

lapply(1:length(prior_functions), function(i) {
  
  print(paste("Computing Bayes factors for", priors[i]))
  
  print("computing for poisson model (small n)")
  # compute Bayes factors for small n and poisson model
  poisson_n100 <- lapply(sim_data_small, 
                           function(x) {
                             get_bayes_factors(sim_data = x,
                                               prior = prior_functions[[i]],
                                               M = "poisson") %>% 
                               as.data.frame() %>%
                               mutate(prior = priors[i])
                           })
  
  
  save_file(poisson_n100,path = paste0("results/Bayesian/priors/", 
                                          priors[i], "/"), 
            overwrite = overwrite_bf)
  
  print("computing for logistic model (small n)")
  
  # compute Bayes factors for small n and logistic model
  logistic_n100 <- lapply(sim_data_small, 
                                   function(x) {
                                     get_bayes_factors(sim_data = x,
                                                       prior = prior_functions[[i]],
                                                       M = "binomial_expanded")%>% 
                                       as.data.frame() %>%
                                       mutate(prior = priors[i])
                                   })
  
  save_file(logistic_n100, path = paste0("results/Bayesian/priors/", priors[i], "/"),
            overwrite = overwrite_bf)
  
  print("computing for poisson model (large n)")
  # compute Bayes factors for large n and poisson model
  poisson_n1500 <- lapply(sim_data_large, 
                           function(x) {
                             get_bayes_factors(sim_data = x,
                                               prior = prior_functions[[i]],
                                               M = "poisson")%>% 
                               as.data.frame() %>%
                               mutate(prior = priors[i])
                           })
  
  save_file(poisson_n1500, path = paste0("results/Bayesian/priors/", priors[i], "/"),
            overwrite = overwrite_bf)
  
  print("computing for logistic model (large n)")
  # compute Bayes factors for small n and logistic model
  logistic_n1500 <- lapply(sim_data_large, 
                            function(x) {
                              get_bayes_factors(sim_data = x,
                                                prior = prior_functions[[i]],
                                                M = "binomial_expanded") %>%
                                as.data.frame() %>%
                                mutate(prior = priors[i])
                            })
  
  save_file(logistic_n1500, path = paste0("results/Bayesian/priors/", priors[i], "/"),
            overwrite = overwrite_bf)
  
})
  

