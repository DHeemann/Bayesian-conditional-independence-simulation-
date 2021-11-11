########### Simulation data creation functions ###########
##########################################################


# function to save files locally 
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

# get u-parameters for simulating data for the poisson GLM representation of the CMH test
get_u_params <- function(u12, k = 2, sampling = "uniform",
                         starting_avg = 300) {
  
  # Args:
  #      u12: test-relevant parameter for the CMH test
  #      k: number of layers in the 2x2xk design
  #      starting_avg: the log() of this value will determine the average value per cell 
  #    
  
  # Apart from u0 and u12, we sample the other u-parameter values 
  # from an uniform distribution between 0 and 1
  # 
  u_parameters <- runif(5,-1, 1)
  
  if (sampling == "fixed_u") {
    u_parameters[u_index] <- u_value # fix u value
  }
  
  u_params <- c(log(starting_avg/(2*2*k)),   
                u_parameters,
                u12)
  
  
  # assigning names to the values
  names(u_params) <-  c("b0", "b1", "b2",
                        paste("b3", 1:(k-1), sep = "_"), 
                        paste("b13", 1:(k-1), sep = "_"), 
                        paste("b23", 1:(k-1), sep = "_"), "b12")
  
  return(u_params)
}



create_data <- function(k = 2,expected_n,u_params) {
  # function that creates a data.frame of frequency values based on the Possion regression model
  # Args:
  #   k: the number of layers (default is 2)
  #   expected_n: The expected sum of frequency values. 
  #               This will affect the size of the of lambda parameters while not affecting the ratio of the 
  #               expected cell counts. 
  #   u_params: the u parameters which will define the ratios of the expected cell counts
  # 
  # Returns: list of u-parameter values, 
  #          the simulateddata.frame of frequency and an array of the same data
  
  # at first, we define a model matrix representing the parameters needed for the CMH test
  M <- model.matrix(Freq ~ Var1+Var2+Var3+Var1*Var3+Var2*Var3+Var1*Var2, 
                    data = data.frame(as.table(array(0, dim = c(2,2,k)))))
  
  lambda_vals <- exp(M%*%u_params) # the expected value of each cell count
  lambda_vals <- lambda_vals * (expected_n/sum(lambda_vals)) # resizing lambdas based on expected value
  
  # we sample observed cell values based on the above defined lambda values
  # from a poisson distribution 
  
  y <- rpois(n = nrow(lambda_vals), lambda = lambda_vals[,1])  
  #y <- round(y * (expected_n/sum(y))) # resizing y based on expected value
  
  # formatting values according to the 2x2xk design for the CMH test
  tbl <- array(y, dim = c(2,2,k))
  ex.2 <- data.frame(as.table(tbl))
 
  return(list(u_params = u_params, df = ex.2, arr = tbl))
}


fixed_n <- function(n, u_per_n  = 10, samples_per_u = 5, u12, 
                    u_sampling = "uniform", n_fixed = TRUE) {
  # creates data sets for a fixed sample size (sum of frequncy values)
  # Args:
  #     n: fixed sum of frequency values (datasets with other values will be rejected)
  #     u_per_n: how many times we want to sample different u-values per sample size value
  #     samples_per_u: how many samples we want to draw per vector of sampled u-values
  #     u12: the specified effect size parameter
  # 
  # for each n value, the function will create u_per_n*samples_per_u data sets
  # we iniate the sampling u-par values multiple times to cancel out potential  
  # differences due to the nuasance parameters 
  print(paste("Creating samples for n=", n))
  
  # set to TRUE if n has to be fixed
  n_fixed <- TRUE
  
  c(replicate(u_per_n, {
    # create u parameter values
    u_params <- get_u_params(u12 = u12, sampling = u_sampling, 
                             starting_avg = n)
    
    replicate(samples_per_u, {
      
      if (n_fixed == FALSE) {
        create_data(k = 2, expected_n = n, u_params = u_params)
      } else {
        while (TRUE) {
        temp_data <- create_data(k = 2, expected_n = n, u_params = u_params)
        freq_temp <- temp_data$df[,"Freq"]
        sample_s <- sum(freq_temp)
        zeros_allowed <- FALSE
        # resample until sample size is reached and no value is 0
        if (sample_s == n & (zeros_allowed|all(freq_temp > 0))) {
          return(temp_data)
          break 
        }
        }
      }
    }, simplify = FALSE)
  }))
}

simulate_data <- function(u12, s_range = c(100, 3000), nuisance_sampling = "uniform",
                          fixed_n = TRUE, step_size = 50,u_per_n  = 10, samples_per_u = 20,k = 2) {
  # function that loops through multiple fixed sample size specifications to create the data sets
  
  # Args: 
  #     u12: test-relevant parameter for the CMH test
  #     s_range: min and max sample size to be simulated
  #     step_size: how many steps there should be between the sample sizes
  #     u_per_n: how many times we want to sample different u-values per sample size value
  #     samples_per_u: how many samples we want to draw per vector of sampled u-values
  
  # returns: a matrix of lists containing the simulated data sets 
  # each column represents one fixed sample size
  sim_df <- pbsapply(seq(s_range[1], s_range[2], by = step_size), 
                     function(x) fixed_n(n = x, u12 = u12, 
                                         u_sampling = nuisance_sampling,
                                         u_per_n  = u_per_n, 
                                         samples_per_u = samples_per_u,
                                         n_fixed = fixed_n))
  
  return(sim_df)
}

# example
# simulate data
# sim_data_u0 <- simulate_data(u12 = 0,
#                              s_range = c(100, 3000),
#                              step_size = 100,
#                              u_per_n  = 100,
#                              samples_per_u = 100, k = 2)
# save(sim_data_u0, file = "sim_data_u0.RData")

