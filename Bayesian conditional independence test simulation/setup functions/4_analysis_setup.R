# u12 values (no effect, small effect , medium effect, large effect)

# Cohen's D
small_d <- 0.2
medium_d <- 0.5
large_d <- 0.8

# effect size for the conditional log odds ratio 
small_u12 <- small_d*pi/sqrt(3)
medium_u12 <- medium_d*pi/sqrt(3)
large_u12 <- large_d*pi/sqrt(3)

u12 <- round(c(0, small_u12, medium_u12, large_u12),3)

##### 1. load data #####

# read in data generated from the null, small, medium and large effect size 
# (sample size range 10 to 100)


#### 2. Create minimal data set #####

if (!"minimal_dataset" %in% ls()) {
#### 3.  Minimal sample size #####
  minimal_dataset <- simulate_data(u12 = 0, 
                                   s_range = c(100,100),
                                   u_per_n = 1, 
                                   step_size = 1,
                                   samples_per_u = 1,k = 2)
  
  minimal_dataset[[1]]$arr <- array(1, dim= c(2,2,2))
  minimal_dataset[[1]]$df$Freq <- rep(1, 8)
}


#### 3. function to create cleaned bayes factors for logistic and log-linear results #### 

create_cleaned_bf <- function(poisson_data, logistic_data,
                              prior_name, column_names) {
  
  # Args: 
  # poisson_data: Matrix containing the poisson bayes factors.
  #               Each column represents one sample size.
  #               Number of rows represent the number of samples per n
  # logistic_data: same as poisson_data but for the logistic test results
  # prior_name: the name of the prior (e.g. "Beta Prime")
  # column_names: the a vector containing the ordered sample size values for the input data 
  
  poisson_nlarge_df <- lapply(seq_along(poisson_data), function(i) {
    
    df_bf_large <- as.data.frame(poisson_data[[i]]) %>%
      select(-prior)
    colnames(df_bf_large) <- column_names
    
    df_bf_large %>% 
      pivot_longer(cols = everything()) %>%
      rename(n = name, bf = value) %>%
      mutate(u12 = u12[i], 
             method = "Poisson", 
             prior = prior_name)
    
  }) %>% 
    do.call(what = rbind)
  
  logistic_nlarge_df <- lapply(seq_along(logistic_data), function(i) {
    
    df_bf_large <- as.data.frame(logistic_data[[i]]) %>%
      select(-prior)
    colnames(df_bf_large) <- column_names
    df_bf_large %>% 
      pivot_longer(cols = everything()) %>%
      rename(n = name, bf = value) %>%
      mutate(u12 = u12[i], 
             method = "Logistic", 
             prior = prior_name)
  }) %>% 
    do.call(what = rbind)
  
  # combine poisson and log-linear 
  bf_df <- rbind(poisson_nlarge_df,
                 logistic_nlarge_df)
  
  bf_df
}

summarise_bf <- function(cleaned_bf, prior_name) {
  # compute median per combination of u12 and sample size
  cleaned_bf %>% 
    mutate(prior = prior_name) %>%
    mutate(n = extract_numeric(n)) %>%
    group_by(n,u12, method, prior) %>%
    summarise(median_bf = median(bf, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(u12 = as.character(u12),
           n = as.numeric(n)) %>%
    mutate(method = case_when(method == "poisson" ~ "Poisson",
                              method == "Logistic expanded" ~ "Logistic",
                              TRUE ~ method)) %>%
    filter(method %in% c("Poisson", "Logistic"))
  
}


create_bf_plots <- function(input_prior, input_prior_name,
                            input_logistic_data, input_poisson_data, 
                            input_column_names) 
  {
  
  
  #### 3. Create bayes factor for minimal n #####
  minimal_logistic <- get_bayes_factors(sim_data = minimal_dataset, 
                                        prior = input_prior,
                                        M = "binomial_expanded")
  
  minimal_poisson <- get_bayes_factors(sim_data = minimal_dataset, 
                                       prior = input_prior)
  
  minimal_df <- data.frame(n = 8, 
                           method = c("Poisson","Logistic"),
                           prior = input_prior_name,
                           median_bf = c(minimal_poisson$bf, 
                                         minimal_logistic$bf)) %>%
    full_join(data.frame(u12 = as.character(u12), 
                         prior= input_prior_name), 
              by = "prior")
  
  
  
  # calculate cleaned BF for large data sets 
  bf_nlarge_df <- create_cleaned_bf(poisson_data = input_poisson_data,
                                    logistic_data = input_logistic_data, 
                                    column_names = input_column_names, 
                                    prior_name = "hyper_gn")
  

  
  # calculate medians 
  bf_median_nlarge_df <- summarise_bf(cleaned_bf = bf_nlarge_df, 
                                      prior_name = input_prior_name)
  
  
  if (max(input_column_names) > 100) {
    plot_breaks <-  seq(min(input_column_names), 
                        max(input_column_names), 
                        200)
    plot_limits <-  range(input_column_names)
  } else {
    plot_breaks <-  plot_breaks <-  seq(min(input_column_names), 
                                        max(input_column_names), 
                                        10)
    plot_limits <-  range(input_column_names)
  }
  
  # plot results
  p_log_poi <- theme_bayesfactor +
    geom_line(data= bf_median_nlarge_df,
              aes(group = u12, 
                  x = n, y = median_bf,
                  col = u12), size = 0.7) +
    ylab("\n\nMedian Bayes factor") +
    xlab("Sample size") +
    facet_wrap(. ~ method, nrow = 2) +
    scale_x_continuous(breaks = plot_breaks, 
                       limits = plot_limits)+
    theme(axis.text = element_text(size = 7), 
          axis.title = element_text(size = 9))
  
  
  # save plot 
  file_name <- paste0("results/Bayesian/plots/",
                     input_prior_name,
                     "_n",
                     max(input_column_names),
                     ".pdf")
 
  
  ggsave(filename = file_name, plot =  p_log_poi)
  print(paste(file_name, "was saved"))
  
}







