########### Compute p-values for the simulated data sets ############
#####################################################################


##### 1. load data #####

# sample size range 20 to 100
if (!"sim_data_small" %in% ls()) {
  sim_data_u0_small <- readRDS("data/sim_data_u0_small.RData")
  sim_data_u1_small <- readRDS("data/sim_data_u1_small.RData")
  sim_data_u2_small <- readRDS("data/sim_data_u2_small.RData")
  sim_data_u3_small <- readRDS("data/sim_data_u3_small.RData")
  
  
  sim_data_small <- list(sim_data_u0_small, 
                         sim_data_u1_small, 
                         sim_data_u2_small, 
                         sim_data_u3_small)
  
  
}

# verify sample sizes for the simulated data set
column_names_small <-  colMeans(structure(vapply(sim_data_u0_small, function(x) sum(x$arr),
                                                 numeric(1)), dim=dim(sim_data_u0_small)))


# (sample size range 100 to 1500)

sim_data_u0_large <- readRDS("data/sim_data_u0_large.RData")
sim_data_u1_large <- readRDS("data/sim_data_u1_large.RData")
sim_data_u2_large <- readRDS("data/sim_data_u0_large.RData")
sim_data_u3_large <- readRDS("data/sim_data_u1_large.RData")
  


####. 2. calculate p-values ####

print("Compute p-values for small n")

p_values_u0_small <- get_pvals(sim_data_u0_small)
p_values_u1_small <- get_pvals(sim_data_u1_small)
p_values_u2_small <- get_pvals(sim_data_u2_small)
p_values_u3_small <- get_pvals(sim_data_u3_small)


print("Compute p-values for large n")

p_values_u0_large <- get_pvals(sim_data_u0_large)
p_values_u1_large <- get_pvals(sim_data_u1_large)
p_values_u2_large <- get_pvals(sim_data_u2_large)
p_values_u3_large <- get_pvals(sim_data_u3_large)

#### 3. save calculated p-values ####

overwrite_p <- TRUE

save_file(p_values_u0_small, path = "results/frequentist/", overwrite = overwrite_p)
save_file(p_values_u1_small, path = "results/frequentist/", overwrite = overwrite_p)
save_file(p_values_u2_small, path = "results/frequentist/", overwrite = overwrite_p)
save_file(p_values_u3_small, path = "results/frequentist/", overwrite = overwrite_p)

save_file(p_values_u0_large, path = "results/frequentist/", overwrite = overwrite_p)
save_file(p_values_u1_large, path = "results/frequentist/", overwrite = overwrite_p)
save_file(p_values_u2_large, path = "results/frequentist/", overwrite = overwrite_p)
save_file(p_values_u3_large, path = "results/frequentist/", overwrite = overwrite_p)







