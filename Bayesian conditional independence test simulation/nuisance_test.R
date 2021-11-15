##### test of nuisance parameters ####

# we need to keep n constant and then loop through 
# u1 from -1 to 1 in steps of 0.1


sapply(list.files("setup functions"), function(x) source(paste0("setup functions/",x)))

check_pattern <- function(bayes_factors, data) {
  
  poisson_df_fixed <- lapply(seq_along(u), function(i) {
    df <- as.data.frame(bayes_factors[[i]])
    n_u0_v2 <- structure(vapply(data[[i]], function(x) sum(x$arr),
                                numeric(1)), dim=dim(data[[i]]))
    df_n <- as.data.frame(n_u0_v2)
    data.frame(df %>% pivot_longer(cols = everything()),
               df_n %>% pivot_longer(cols = everything())) %>%
      dplyr::select(bf = value, n = value.1) %>%
      mutate(u = u[i])
  }) %>% do.call(what = rbind) %>%
    filter(!is.na(bf))
  
  
  plot_df <- poisson_df_fixed %>%
    mutate(n = as.numeric(n)) %>% 
    group_by(u, n_group = floor(n/100)*100) %>%
    summarise(median_bf = median(bf, na.rm = TRUE)) %>%
    filter(n_group <= end) 
  
  
  plot_df %>% 
    ungroup() %>%
    mutate(u = as.factor(u)) %>%
    ggplot(aes(x = n_group, y = median_bf)) +
    #scale_y_continuous(trans = "log10")+
    geom_line(aes(col = u, group = u)) +
    geom_point(aes(col = u)) +
    #geom_hline(yintercept = 1)+
    ggpubr::theme_pubr()+ 
    theme(legend.position = "right")+
    xlab("Sample size") + ylab("Median BF") +
    scale_color_discrete(name = unname(TeX("u_{1}")))
  
}


#u <- seq(-1, 1, by = 0.5)

# here we fix u1
u_index <- 1

u <- c(-1, 0, 1)

start <- 50
end <- 1000
step <- 50

# create data

set.seed(201918)
u1_fixed_data_medium <- lapply(u, function(u1) {
  u_value <<- u1
  simulate_data(u12 = 0.4, s_range = c(start, end), step_size = step, u_per_n = 50,
                fixed_n = FALSE,samples_per_u = 20, nuisance_sampling = "fixed_u")
})

set.seed(201917)
u1_0_fixed_data_medium <- lapply(u, function(u1) {
  u_value <<- u1
  simulate_data(u12 = 0, s_range = c(start, end), step_size = step, u_per_n = 50,
                fixed_n = FALSE,samples_per_u = 20, nuisance_sampling = "fixed_u")
})

set.seed(201916)
u1_data_medium <- lapply(u, function(u1) {
  simulate_data(u12 = 0.4, s_range = c(start, end), step_size = step, u_per_n = 100,
                fixed_n = FALSE,samples_per_u = 50, nuisance_sampling = "uniform")
})


# create Bayes factors

u1_fixed_log_exp_hyper <- lapply(u1_fixed_data_medium, function(x) {
  get_bayes_factors(sim_data = x, 
                    prior = g(100),
                    M = "binomial_expanded")
})

u1_fixed_poisson_g100 <- lapply(u1_fixed_data_medium, function(x) {
  get_bayes_factors(sim_data = x, 
                    prior = g(100))
})

u1_0_fixed_poisson_g100 <- lapply(u1_0_fixed_data_medium, function(x) {
  get_bayes_factors(sim_data = x, 
                    prior = g(100))
})

u1_fixed_poisson_hgn <- lapply(u1_fixed_data_medium, function(x) {
  get_mod_hgn_bf(sim_data = x, n_flex = TRUE)
})

u1_0_fixed_poisson_hgn <- lapply(u1_0_fixed_data_medium, function(x) {
  get_mod_hgn_bf(sim_data = x, n_flex = TRUE)
})


u1_fixed_hgn_exp_log_medium <- lapply(u1_fixed_data_medium, function(x) {
  get_mod_hgn_bf_expanded_logistic(sim_data = x)
  
})

u1_0_fixed_hgn_exp_log_medium <- lapply(u1_0_fixed_data_medium, function(x) {
  get_mod_hgn_bf_expanded_logistic(sim_data = x)

})


u1_exp_bi_robust_medium <- lapply(u1_data_medium, function(x) {
  get_bayes_factors(sim_data = x, 
                    prior = robust(),
                    M = "binomial_expanded")
})

u1_exp_bi_robust_medium<- lapply(u1_data_medium, function(x) {
  get_bayes_factors(sim_data = x, 
                    prior = robust(),
                    M = "binomial_expanded")
})


u1_fixed_log_exp_intrinsic_medium <- lapply(u1_fixed_data_medium, function(x) {
  get_bayes_factors(sim_data = x, 
                    prior = intrinsic(),
                    M = "binomial_expanded")
})


u1_fixed_binomial_beta_prime_medium <- lapply(u1_fixed_data_medium, function(x) {
  get_bayes_factors(sim_data = x, 
                    prior = beta.prime(),
                    M = "binomial_expanded")
})



check_pattern(u1_fixed_hgn_exp_log_medium, data = u1_fixed_data_medium) +
  ggtitle(unname(TeX("u_{12} = 0.4"))) +
  scale_y_continuous(trans = "log10") +
  ylab("Median BF")

check_pattern(u1_0_fixed_hgn_exp_log_medium, data = u1_0_fixed_data_medium) +
  ggtitle(unname(TeX("u_{12} = 0"))) +
  scale_y_continuous(trans = "log10") 


