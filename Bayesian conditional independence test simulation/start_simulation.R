##### Script to initiate the simulation #####

# The simulation studies the effect of different effect sizes and g-priors
# on the Bayes factors for conditional independence tests in 
# 2x2x2 designs


#### 1.Set up #### 

# set path for simulation project 
user_path  <- "C:/Users/Daniel/Desktop/Bayesian conditional independence test simulation"

setwd(user_path)


# load and install relevant packages 
# if all packages are already installed, they will not be installed again
source("install_packages.R")

#### 2. Run all the set-up scripts needed to start the simulation ####

sapply(list.files("setup functions"), function(x) {
  source(paste0("setup functions/",x))
  })

#### 3. Create data ####

# If over_write_data = FALSE, new data will be only simulated if not present in 
# in the data folder
overwrite_data <- TRUE

source("create_data_loglinear.R")


### 4. Compute p-values #### 

source("statistical_tests/Frequentist/p_value_computations.R")

#### 4. Compute Bayes factors ####

source("statistical_tests/Bayesian/priors/compute_bayes_factors.R")

#### 5. Create figures ##### 

# for p_values 
source("statistical_tests/Frequentist/p_value_comparisons.R")

# for Bayes factors
source("statistical_tests/Bayesian/priors/create_plots.R")



