########## Create data using the log-linear representation of the CMH test ##############

##### 1. Setup ######


# Note that in this script
# u0 refers to no effect, u1 refers to small effect, 
# u2 refers to medium effect, u3 refers to large effect size

# Cohen's D
small_d <- 0.2
medium_d <- 0.5
large_d <- 0.8

# conditiona log odds ratio effect size
small_u12 <- small_d*pi/sqrt(3)
medium_u12 <- medium_d*pi/sqrt(3)
large_u12 <- large_d*pi/sqrt(3)


save_file <- function(data, overwrite=FALSE, path) {
  str_name <- as.character(paste0(path, deparse(substitute(data)), ".RData"))
  if (!file.exists(str_name)|overwrite) {
    saveRDS(data, str_name)
    if (file.exists(str_name)) {
      print("file will be overwritten")
    }
    print(paste(deparse(substitute(data)), "was saved in", paste0(getwd(), path)))
  } else {
    print(paste(str_name, "already exists."))
  }
}


# function to save files locally 
if (!"overwrite_data" %in% ls()) {
  overwrite_data <- FALSE
}

if (!"u_per_n_input" %in% ls()) {
  u_per_n_input <- 25
}

if (!"samples_per_u_input" %in% ls()) {
  samples_per_u_input <- 25
}


#### 2. Simulate data for small sample size ranges #####

if (overwrite_data | !("sim_data_u0_small.RData" %in% list.files("data"))) {

  set.seed(20210919)
  # simulate data for u0 and small sample size ranges 
  sim_data_u0_small <- simulate_data(u12 = 0,  
                                     s_range = c(20, 100),
                                     step_size = 10,
                                     u_per_n  = u_per_n_input,
                                     samples_per_u = samples_per_u_input, 
                                     k = 2)
  
  save_file(sim_data_u0_small, path = "data/", 
            overwrite = overwrite_data)
  
}

if (overwrite_data | !("sim_data_u1_small.RData" %in% list.files("data"))) {
  
  set.seed(20200919)
  # simulate data for u1 and small sample size ranges 
  sim_data_u1_small <- simulate_data(u12 = small_u12, 
                                     s_range = c(20, 100),
                                     step_size = 10,
                                     u_per_n  = u_per_n_input,
                                     samples_per_u = samples_per_u_input, 
                                     k = 2)
  
  save_file(sim_data_u1_small, path = "data/", 
            overwrite = overwrite_data)
}



if (overwrite_data | !("sim_data_u2_small.RData" %in% list.files("data"))) {
 
   set.seed(20200920)
  # simulate data for u2 and small sample size ranges 
  sim_data_u2_small <- simulate_data(u12 = medium_u12, 
                                     s_range = c(20, 100),
                                     step_size = 10,
                                     u_per_n  = u_per_n_input,
                                     samples_per_u = samples_per_u_input, 
                                     k = 2)
 
   save_file(sim_data_u2_small, path = "data/", 
            overwrite = overwrite_data)
}


if (overwrite_data | !("sim_data_u2_small.RData" %in% list.files("data"))) {
  set.seed(20200921)
  # simulate data for u3 and small sample size ranges 
  sim_data_u3_small <- simulate_data(u12 = large_u12, 
                                     s_range = c(20, 100),
                                     step_size = 10,
                                     u_per_n  = u_per_n_input,
                                     samples_per_u = samples_per_u_input, 
                                     k = 2)
  
    save_file(sim_data_u3_small, path = "data/", 
            overwrite = overwrite_data)

}

####### 3. Simulate data for large sample size ranges ########

if (overwrite_data | !("sim_data_u0_large.RData" %in% list.files("data"))) {

  set.seed(20210920)
  # simulate data for u0 and large sample size ranges 
  sim_data_u0_large <- simulate_data(u12 = 0,  
                                     s_range = c(100, 3000),
                                     step_size = 100,
                                     u_per_n  = u_per_n_input,
                                     samples_per_u = samples_per_u_input, 
                                     k = 2)
  
  # save files for large sample size ranges
  save_file(sim_data_u0_large, path = "data/", 
            overwrite = overwrite_data)

}


 if (overwrite_data | !("sim_data_u1_large.RData" %in% list.files("data"))) {
  set.seed(20210921)
  # simulate data for u1 and large sample size ranges 
  sim_data_u1_large <- simulate_data(u12 = small_u12,  
                                     s_range = c(100, 3000),
                                     step_size = 100,
                                     u_per_n  = u_per_n_input,
                                     samples_per_u = samples_per_u_input, 
                                     k = 2)
  # save files for large sample size ranges
  save_file(sim_data_u1_large, path = "data/", 
            overwrite = overwrite_data)
  
}

if (overwrite_data | !("sim_data_u2_large.RData" %in% list.files("data"))) {
  
  set.seed(20210922)
  # simulate data for u2 and large sample size ranges 
  sim_data_u2_large <- simulate_data(u12 = medium_u12, 
                                     s_range = c(100, 3000),
                                     step_size = 100,
                                     u_per_n  = u_per_n_input,
                                     samples_per_u = samples_per_u_input, 
                                     k = 2)
  
  # save files for large sample size ranges
  save_file(sim_data_u2_large, path = "data/", 
            overwrite = overwrite_data)
  

}

if (overwrite_data | !("sim_data_u3_large.RData" %in% list.files("data"))) {
  set.seed(20230922)
  # simulate data for u3 and large sample size ranges 
  sim_data_u3_large <- simulate_data(u12 = large_u12, 
                                     s_range = c(100, 3000),
                                     step_size = 100,
                                     u_per_n  = u_per_n_input,
                                     samples_per_u = samples_per_u_input, 
                                     k = 2)
  
  # save files for large sample size ranges
  save_file(sim_data_u3_large, path = "data/", 
            overwrite = overwrite_data)
}





