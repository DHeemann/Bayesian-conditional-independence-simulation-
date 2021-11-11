################### Hyper g/n prior Bayes factor #######################

# In this script, the Bayes factor using the hyper-g/n prior for the 
# conditional independence test in 2x2x2 contigency tables 
# will be implemented from scratch with most parts 
# following the specifictions in the paper 
# "Mixtures of g-Priors in Generalized Linear Models"

# The difference to the original implementation in the BAS R-Package
# is that this script accomodates for large sample sizes which lead to NA 
# values in the BAS implementation.

# While the script uses the hyper-g/n prior, changing its parameters, other prior
# variants can be recreated as well. 

# M7 is the model including the conditional interaction term between rows and colums
# M4 is the restricted model 

#### 1. Setup #####

library(BAS)
library(survival)
library(dplyr)
library(tidyr)
library(aod)
library(numDeriv)


#### 2. Set-up functions and objects #####

# these likelihood functions will be used in the main function below
L7_poisson <- function(b0) {
  # loglikelihood for poisson model (m7)
  return(-sum(exp(M7_poisson%*%c(b0,beta_parameters_m7))) +sum(y* (M7_poisson%*%c(b0, beta_parameters_m7))) - sum(lfactorial(y)))
}

L0_poisson <- function(b0) {
  # loglikelihood for poisson model (m0)
  return(-sum(exp(M0_poisson%*%b0)) +sum(y* (M0_poisson%*%b0)) - sum(lfactorial(y)))
}

L4_poisson <- function(b0) {
  # loglikelihood for poisson model (m4)
  return(-sum(exp(M4_poisson%*%c(b0,beta_parameters_m4))) +sum(y*(M4_poisson%*%c(b0,beta_parameters_m4))) - sum(lfactorial(y)))
}

# creating relevat model matrices 
M7_poisson <- model.matrix(Freq ~ Var1+Var2+Var3+Var1*Var3+Var2*Var3+Var1*Var2, 
                   data = data.frame(as.table(array(0, dim = c(2,2,2)))))

M0_poisson <- model.matrix(Freq ~ 1, 
                   data = data.frame(as.table(array(0, dim = c(2,2,2)))))

M4_poisson <- model.matrix(Freq ~ Var1+Var2+Var3+Var1*Var3+Var2*Var3, 
                   data = data.frame(as.table(array(0, dim = c(2,2,2)))))

# parameters for the hyper g/n function
a1 = 1
b1 = 2
r1 = 1.5
s1 = 0
v1 = 1


#### 3. Hyper g/n prior Bayes factor function ####


manual_hyper_gn <- function(sim_data, no_zero = TRUE, k1) {
  # Computes the conditional independence Bayes factor using the hyper g/n prior
  #
  # Args:
  #  sim_data: a data.frame object with the frequency table
  #  no_zero: Boolean indicating if 0 values are permitted as cell values
  #  k1: k parameter, which is set to  (1/n) in the hyper g/n
  #
  # Returns: 
  #     Numeric Bayes factor of M7 over M4. 
  #
  #
  freq_df <- as.data.frame(sim_data)
  
  # stop here if 0s are not permitted by the user
  y <<- freq_df$Freq 
  if (no_zero & any(y == 0)) {
    return(NA)
    break 
  }
  
  
  ##  fit models to get MLE estimates 
  # model 7 (with u12 term)
  m7 <- glm(Freq ~ Var1 + Var2 + Var3 + Var1*Var3 +Var2*Var3 + Var1*Var2 , 
            data = freq_df, family = poisson)
  
  # model 4 (without u12 term)
  m4 <- glm(Freq ~ Var1 + Var2 + Var3 + Var1*Var3 +Var2*Var3, 
            data = freq_df, family = poisson)
  
  # empty model as reference 
  m0 <- glm(Freq ~ 1, data = freq_df, family = poisson)
  
  ## get observed fisher information of fitted models
  # observed fisher information for MLE Intercept of Null model
  obfish_m0 <- -hessian(func=L0_poisson, x=coef(m0)[1])[1] 
  
  # observed fisher information for MLE Intercept of m7
  beta_parameters_m7 <<- coef(m7)[-1] #!
  obfish_m7 <- -hessian(func=L7_poisson, x=coef(m7)[1])[1] # same as result from manual calculations
  
  # observed fisher information for MLE Intercept of m4
  beta_parameters_m4 <<- coef(m4)[-1]
  obfish_m4 <- -hessian(func=L4_poisson, x=coef(m4)[1])[1] 
  
  # observed fisher information for MLE betas (except intercept) for m7 and m4
  In7 <- solve(vcov(m7)[-1,-1])
  In4 <- solve(vcov(m4)[-1,-1])
  
  # 2*log likelihood ratio 
  zm7 <- 2*(logLik(m7) - logLik(m0))[1]
  # same as anova(m7, m0, test = "Chisq")
  
  zm4 <- 2*(logLik(m4) - logLik(m0))[1]
  # same as anova(m4, m0, test = "Chisq")
  
  # Wald statistic for m7 and m4
  Qm7 <- (t(coef(m7)[-1])%*%In7%*%coef(m7)[-1])[1,1]
  Qm4 <- (t(coef(m4)[-1])%*%In4%*%coef(m4)[-1])[1,1]
  
  Pm_7 <- qr(M7_poisson)$rank - 1 # number of predictors
  Pm_4 <- qr(M4_poisson)$rank - 1 # number of predictors
  
  # set-up for humbert function
  a <- b1/2
  b <- r1
  x <- 1-k1
  # c and y for model 7
  c7 <- (a1+b1+Pm_7)/2
  y7 <- (s1+Qm7)/(2*v1)
  # c and y for model 4 
  c4 <- (a1+b1+Pm_4)/2
  y4 <- (s1+Qm4)/(2*v1)
  
  # Here we are choosing a small enough constant that allows us to compute the 
  # Bayes factor even if the Wald statistic is large
  # 
  max_value <- 300 # max_value should be defined such that exp(max_value) can still be computed in R
  max_y <- max(y7, y4) # we need to accomodate for the model yielding the larger Wald statistic
  rep_y <- ceiling(max_y/max_value) # defines how many times we split up the exp(y*t) term
  
  # small_constant defines by which value we multiply exp(y*t)
  small_constant <- 1/exp(max(0, (max_y - max_value)/rep_y))
  
  # Humbert function for model 7 to be integrated 
  integral_humbert7 <- Vectorize(function(t) {
    t^(a-1)*(1-t)^(c7-a-1)*(1-x*t)^(-b)*
      prod(small_constant*exp(rep((y7)*t/rep_y, rep_y)))
  })
  
  # Humbert function for model 4 to be integrated
  integral_humbert4 <- Vectorize(function(t) {
    t^(a-1)*(1-t)^(c4-a-1)*(1-x*t)^(-b)*
      prod(small_constant*exp(rep((y4)*t/rep_y, rep_y)))
  })
  
  # used in the Bayes factor formula 
  C_2 <- BAS::phi1(b1/2, 
                   r1, 
                   (a1+b1)/2,
                   s1/(2*v1),
                   1-k1)
  
  # Output of the Phi function for model 7
  C7_1 <- gamma(c7)/prod(gamma(a), gamma(c7-a)) * 
    integrate(f = integral_humbert7, 
              lower = 0, 
              upper = 1)[1]$value
  
  # Output of the Phi function for model 4 
  C4_1 <- gamma(c4)/prod(gamma(a), gamma(c4-a)) *
    integrate(f = integral_humbert4, 
              lower = 0, 
              upper = 1)[1]$value
  
  # this part of the Bayes factor formula is rewritten from exp(a)/exp(b)  
  # to exp(a-b) to make large values in the Wald statistic computationally feasible  
  md <- exp((zm7/2 - Qm7/(2*v1))- (zm4/2 - Qm4/(2*v1)))
  
  # numerator of the Bayes factor 
  chh_bf_7_0 <- prod(sqrt(obfish_m0/obfish_m7)*v1^(-Pm_7/2) * md,
                     (beta((a1+Pm_7)/2, b1/2)*C7_1)/(beta(a1/2,b1/2)*C_2))
  
  # denominator of the Bayes factor 
  chh_bf_4_0 <- prod(sqrt(obfish_m0/obfish_m4)*v1^(-Pm_4/2),
                     (beta((a1+Pm_4)/2, b1/2)*C4_1)/(beta(a1/2,b1/2)*C_2))
  

  
  # obtaining the Bayes factor
  (chh_bf_7_0/chh_bf_4_0)
  
}

#### 5. Example and validity check ####

example_df_small <- structure(list(Var1 = structure(c(1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L),
                                .Label = c("A", "B"), class = "factor"), 
               Var2 = structure(c(1L,  1L, 2L, 2L, 1L, 1L, 2L, 2L), 
                                .Label = c("A", "B"), class = "factor"), 
               Var3 = structure(c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L), 
                                .Label = c("A", "B"), class = "factor"), 
               Freq = c(12L, 12L, 9L, 15L, 7L, 9L, 13L, 23L)), 
          class = "data.frame", row.names = c(NA, -8L))


example_df_large <- structure(list(Var1 = structure(c(1L, 2L, 1L, 2L, 1L, 2L, 1L, 2L), 
                                .Label = c("A", "B"), class = "factor"), 
               Var2 = structure(c(1L, 1L, 2L, 2L, 1L, 1L, 2L, 2L),
                                .Label = c("A", "B"), class = "factor"), 
               Var3 = structure(c(1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L),
                                .Label = c("A", "B"), class = "factor"), 
               Freq = c(120L, 85L, 266L, 301L,  101L, 146L, 523L, 958L)), 
          class = "data.frame", row.names = c(NA,-8L))


# function to return only Bayes factor 
get_bf <- function(data,prior, model) {
  if (model == "poisson") {
    b <-  bas.glm(Freq ~ Var1 + Var2 + Var3 + Var1*Var3 +Var2*Var3 + Var1*Var2, 
                  data=data,
                  family=poisson(),
                  betaprior= prior, modelprior=uniform(),
                  include.always = "~ 1 + Var1 + Var2 + Var3 + Var1*Var3 + Var2*Var3",
                  n.models=2^10, MCMC.iterations=10,
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
                  MCMC.iterations=10,
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
  },warning = function(cond) {
    print("something is wrong with this case, check summary of b: ")
    print(summary(b))
    #print(x)
    break 
  })
  
  return(bf)
}

# BAS results
get_bf(example_df_small, prior = hyper.g.n(),
       model = "poisson")

# comparing own implementation with BAS package
manual_hyper_gn(sim_data = example_df_small, k1 = 1/8)


# setting k1 to 1, we create the hyper g prior
manual_hyper_gn(sim_data = example_df_small, k1 = 1)
get_bf(example_df_small, prior = hyper.g(),
       model = "poisson")

# below is an example where the BF is NA using BAS

# BAS results: summary of bas.glm() will return NA in this case
# get_bf(example_df_large, prior = hyper.g.n(),
#        model = "poisson")


# comparing own implementation with BAS package (works)
manual_hyper_gn(sim_data = example_df_large,k1 = 1/8)





