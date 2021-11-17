u1_fixed_data_medium <- readRDS("data/u1_fixed_data_medium.RData")


u1_fixed_04_log_g100 <- lapply(u1_fixed_data_medium, function(x) {
  get_bayes_factors(sim_data = x, 
                    prior = g.prior(100),
                    M = "binomial_expanded")
})

u1_fixed_04_poi_g100 <- lapply(u1_fixed_data_medium, function(x) {
  get_bayes_factors(sim_data = x, 
                    prior = g.prior(100))
})



u1_fixed_04_log_hgn <- lapply(u1_fixed_data_medium, function(x) {
  get_bayes_factors(sim_data = x, 
                    prior = hyper.g.n(), 
                    M = "binomial_expanded")
                    
})

u1_fixed_04_poi_hgn <- lapply(u1_fixed_data_medium, function(x) {
  get_bayes_factors(sim_data = x, 
                    prior = hyper.g.n())
})

u1_fixed_04_log_beta_prime <- lapply(u1_fixed_data_medium, function(x) {
  get_bayes_factors(sim_data = x, 
                    prior = beta.prime(), 
                    M = "binomial_expanded")
  
})

u1_fixed_04_poi_beta_prime <- lapply(u1_fixed_data_medium, function(x) {
  get_bayes_factors(sim_data = x, 
                    prior = beta.prime())
})

