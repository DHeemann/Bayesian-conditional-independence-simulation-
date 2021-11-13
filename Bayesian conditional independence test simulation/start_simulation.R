##### Script to initiate the simulation #####

# The simulation studies the effect of different effect sizes and g-priors
# on the Bayes factors for conditional independence tests in 
# 2x2x2 designs


#### 1.Set up #### 

# set path for simulation project 
# user_path would e.g.  "C:/Users/your_name/Desktop"
user_path <- "path_where_folder_is_stored"

setwd(paste0(user_path, "/Bayesian conditional independence test simulation"))

# create subfolders where Bayes factor results per prior will be stored.
priors <- c("beta_prime", 
            "fixed_g", 
            "hyper_gn", 
            "intrinsic", 
            "robust", 
            "unit_information")

sapply(priors, function(x) dir.create(file.path("results/Bayesian/priors", x)))


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
       
# determines how many samples will be simulated
if (!"u_per_n_input" %in% ls()) {
            u_per_n_input <- 25
            }
if (!"samples_per_u_input" %in% ls()) {
            samples_per_u_input <- 25
            }
    
samples_per_u_input <- 25

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



