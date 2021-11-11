# Bayesian conditional independence test - simulation study 
Simulation study in R on Bayesian conditional independence tests in 2x2x2 contingency tables using different (g-)prior variants. 

The Bayes factors are computed using the R Package BAS (Bayesian Variable Selection and Model Averaging using Bayesian Adaptive Sampling).

The figure below gives an overview of how the simulation process works. 
<img width="617" alt="simulation_process" src="https://user-images.githubusercontent.com/36103689/141364834-365e2e5a-e869-4e4c-ac9a-ad293a6c7cf7.PNG">


The script "start_simulation.R" will initiate the simulation by: 
(1) Loading the needed packages and installing them if not already installed 
(2) Running all the set-up scripts (for the functions needed later on) 
(3) Creating the simulated data sets (sourcing the script "create_data_loglinear.R" -> numbers of samples can be adjusted here to reduce computing time)
(4) Computing p-values (CMH, logistic and log-linear model comparison test)
(5) Computing Bayes factors 
(6) Creating plots for the p-values
(7) Creating plots for the Bayes factors 
