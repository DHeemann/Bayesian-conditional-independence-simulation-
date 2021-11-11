#######  CMH P-Value  functions #################
#################################################

# Note that all functions are set to return NA in case of 0 count cells.

get_p_ll <- function(data) {
  
  # retrieve p-value based on Poisson presentation of the CMH test
  # Args:
  #    data: data.frame containing the frequencies
  
  if (any(data$Freq == 0)) {
    return(NA)
    break
  }
  
  M7 <- glm(Freq ~ Var1 + Var2 + Var3 + Var1*Var3 +Var2*Var3 + Var1*Var2 , 
            data = data, family = poisson)
  M4 <- glm(Freq ~ Var1 + Var2 + Var3 + Var1*Var3 +Var2*Var3 , 
            data = data, family = poisson)
  
  
  return(anova(M4,M7, test = "Chisq")[2,5])
}


get_p_log <- function(data) {
  # retrieve p-value based on Logistic presentation of the CMH test
  # Args:
  #    data: data.frame containing the frequencies
  if (any(data$Freq == 0)) {
    return(NA)
    break
  }
  rows <- cbind(A = data[data$Var1 == "A", "Freq"],B= data[data$Var1 == "B", "Freq"])
  columns <- c("A", "B", "A", "B")
  layers <- c("A", "A", "B", "B")
  #    data: data.frame containing the frequencies
  M7 <- glm(rows ~ columns + layers, family = binomial(link = "logit"))
  M4 <- glm(rows ~ layers, family = binomial(link = "logit"))
  
  return(anova(M4,M7, test = "Chisq")[2,5])
}


get_p_CMH <- function(data, correct_cmh = TRUE) {
  # retrieve p-value based on Logistic presentation of the CMH test
  # Args:
  #    data: data.frame containing the frequencies
  if (any(data == 0)) {
    return(NA)
    break
  }
  
  return(mantelhaen.test(data, correct = correct_cmh)[[3]])
}


get_pvals <- function(sim_df) {
  # obtain CMH p-values based on simulated data sets
  # Computes (1) the Poisson CMH test, 
  #          (2) the CMH test (from the stats package) with correction
  #          (3) the CMH test without correction 
  # Args:
  #    sim_df: Simulated data
  print("computing p-values from log-linear model")
  p_vals_ll <- structure(vapply(sim_df, function(x) get_p_ll(x$df), 
                                numeric(1)), dim=dim(sim_df))
  #colnames(p_vals_ll) <-  paste("n_",seq(s_range[1],, s_range[2], by = step_size),sep = "")
  
  print("computing p-values from logistic model")
  p_vals_logistic <- structure(vapply(sim_df, function(x) get_p_log(x$df), 
                                      numeric(1)), dim=dim(sim_df))
  
  print("computing p-values based on original CMH test (stats:mantelhaen.test)")
  p_vals_cmh <- structure(vapply(sim_df, function(x) get_p_CMH(x$arr), 
                                 numeric(1)), dim=dim(sim_df))
  #colnames(p_vals_cmh) <-  paste("n_",seq(s_range[1], s_range[2], by = step_size),sep = "")
  
  print("computing p-values based on original CMH test (stats:mantelhaen.test) without any correction")
  p_vals_cmh_no_correction <-  structure(vapply(sim_df, function(x) get_p_CMH(x$arr, correct_cmh = FALSE), 
                                                numeric(1)), dim=dim(sim_df))
  #colnames(p_vals_cmh_no_correction) <-  paste("n_",seq(s_range[1], s_range[2], by = step_size),sep = "")
  
  return(list(log_linear = p_vals_ll, 
              logistic = p_vals_logistic,
              CMH = p_vals_cmh, 
              CMH_no_correct = p_vals_cmh_no_correction))
}
