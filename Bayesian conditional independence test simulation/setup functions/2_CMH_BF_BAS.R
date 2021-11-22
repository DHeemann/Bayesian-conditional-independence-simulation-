#######  CMH Bayes factor functions #############
#################################################



get_bf <- function(data,prior, model) {
  # compute CMH Bayes factor based on the BAS package
  # Args:
  #     data: data.frame of frequency values 
  #     prior: (G-)prior to be used in the bas.glm function
  #     model: either Poisson or Binomial model
  # 
  # Returns: the CMH Bayes factor as numerical value
  #  
  # if (any(data$Freq == 0)){
  #   return(NA)
  # }
  if (model == "poisson") {
    b <-  bas.glm(Freq ~ Var1 + Var2 + Var3 + Var1*Var3 +Var2*Var3 + Var1*Var2, 
                  data=data,
                  family=poisson(),
                  betaprior= prior, modelprior=uniform(),
                  include.always = "~ 1 + Var1 + Var2 + Var3 + Var1*Var3 + Var2*Var3",
                  n.models=2^10, MCMC.iterations=1000,
                  prob.rw=.95)
  }
  if (model == "binomial_expanded") {
    
    data_expanded <- data[rep(seq_len(nrow(data)), data$Freq),]
    data_expanded$Var1 <- as.numeric(data_expanded$Var1)-1
    
    b <-  bas.glm(Var1 ~ Var3 + Var2, 
                  data= data_expanded,
                  family=binomial(),
                  betaprior= prior, modelprior=uniform(),
                  include.always = "~ 1 + Var3", n.models=2^10, 
                  MCMC.iterations=1000,
                  prob.rw=.95)
  }
  z <- summary(b)
  bf <- tryCatch({
    # in case both BF's in the summary are 1, return 1 
    if (sum(summary(b)["BF",], na.rm = T) == 2) {
      bf <- 1
    } 
    # in case the model with the u12 term has the highest Bayes factor, return 1/(the other bayes factor)
    else if (which(summary(b)["BF",] != 1) == 2) {
      bf <- 1/(z["BF", 2])
    } 
    
    else{
      bf <- z["BF", 3]
    }
  }, warning = function(cond) {
    print("something is wrong with this case, check summary of b: ")
    print(summary(b))
    #print(x)
    break 
  })
  
  return(bf)
}


get_bayes_factors <- function(sim_data, prior = EB.local(), M = "poisson") {
  # apply BF function to multiple simulated datasets
  # 
  # Args:
  #     sim_datadata: multiple data.frames containing frequency values 
  #     prior: (G-)prior to be used in the bas.glm function
  #     model: either Poisson or Binomial model
  # 
  # Returns: Bayes factors in the same dimensions structure as the sim_data object
  # 
  bf <- structure(pbsapply(sim_data, function(x) tryCatch({ 
  
  if (class(prior) != "prior") {
    
    # unit information prior with n = sum(m)
    if (prior == "Unit Information") {
      prior <- g.prior(g = sum(x$arr))
    } else if  (prior == "Beta_prime_n_flex") {
      prior <- beta.prime(sum(x$arr))
    }
    
  }
    
    get_bf(x$df, prior = prior, model = M)
  }, error = function(cond) { return(NA)
  })), dim = dim(sim_data))
  
  return(list(bf = bf))
}


#get_bayes_factors(sim_data = sim_data_u0_large[1], prior = "Unit Information")
#get_bayes_factors(sim_data = sim_data_u0_large[1], prior = g.prior(100))

                           
                           
# function to compute modified hyper g bayes factor across many data sets
# based on Poisson regression models
# Requires the function hyper_gn_bf.R 
get_mod_hgn_bf <- function(sim_data, n_flex = FALSE,k_param) {
  # apply modified Hyper G/N Bayes factor to multiple datasets
  # 
  # Args:
  #     sim_data: multiple data.frames containing frequency values 
  # 
  # Returns: Bayes factors in the same dimensions structure as the sim_data object
  # 
  bf <- structure(pbsapply(sim_data, function(x) #tryCatch({
    if (n_flex) {
      manual_hyper_gn(sim_data = x$df, 
                      k1 = 1/sum(x$arr))
    
    } else {
      manual_hyper_gn(sim_data = x$df, 
                      k1 = k_param)
    }
 # }, 
#error = function(cond) { return(NA)
#  }
#)
), dim = dim(sim_data))
  
  return(list(bf = bf))
}


