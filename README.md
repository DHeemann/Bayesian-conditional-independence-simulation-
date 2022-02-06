# Bayesian conditional independence test - simulation study 
Simulation study in R on Bayesian conditional independence tests in 2x2x2 contingency tables using different (g-)prior variants. 

The whole project with subfolders can be downloaded from github with the following link:

The Bayes factors are computed using the R Package BAS (Bayesian Variable Selection and Model Averaging using Bayesian Adaptive Sampling).

The conditional independence test (and the data creation) is based on a log-linear representation of the CMH test as discussed in the paper "Log-linear representations of the mantel-haenszel and the breslow-day tests" (von Eye and Indurkhya, 2000). The relevant effect size (conditional log odds ratio) is referred to as u_12 in the figure below, while the other u-paramters are nuisance parameters. 

The figure below gives an overview of how the simulation process works. 
<img width="617" alt="simulation_process" src="https://user-images.githubusercontent.com/36103689/141364834-365e2e5a-e869-4e4c-ac9a-ad293a6c7cf7.PNG">


The script "start_simulation.R" will initiate the simulation by: 
- Loading the needed packages and installing them if not already installed 
- Running all the set-up scripts (for the functions needed later on) 
- Creating the simulated data sets (sourcing the script "create_data_loglinear.R" -> numbers of samples can be adjusted here to reduce computing time)
- Computing p-values (CMH, logistic and log-linear model comparison test)
- Computing Bayes factors 
- Creating plots for the p-values
- Creating plots for the Bayes factors 
