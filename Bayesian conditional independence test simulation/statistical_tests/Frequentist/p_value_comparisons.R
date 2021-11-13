######### P-value comparisons #############
###########################################


# Set-up 

# Cohen's D
small_d <- 0.2
medium_d <- 0.5
large_d <- 0.8

# conditiona log odds ratio effect size
small_u12 <- small_d*pi/sqrt(3)
medium_u12 <- medium_d*pi/sqrt(3)
large_u12 <- large_d*pi/sqrt(3)

#### 1. Read in p-value data ##### 

# for small sample size range (10-100)
p_values_u0_small <- readRDS("results/frequentist/p_values_u0_small.RData")
p_values_u1_small <- readRDS("results/frequentist/p_values_u1_small.RData")
p_values_u2_small <- readRDS("results/frequentist/p_values_u2_small.RData")
p_values_u3_small <- readRDS("results/frequentist/p_values_u3_small.RData")

# read in data generated from the null, small, medium and large effect size 
# (sample size range 110 to 3000)
sim_data_u0_large <- readRDS("data/sim_data_u0_large.RData")
sim_data_u1_large <- readRDS("data/sim_data_u1_large.RData")
sim_data_u2_large <- readRDS("data/sim_data_u2_large.RData")
sim_data_u3_large <- readRDS("data/sim_data_u3_large.RData")


p_values_small <- list(u0 = p_values_u0_small, 
                       u1 = p_values_u1_small,
                       u2 = p_values_u2_small,
                       u3 = p_values_u3_small)


p_values_u0_large <- readRDS("results/frequentist/p_values_u0_large.RData")
p_values_u1_large <- readRDS("results/frequentist/p_values_u1_large.RData")
p_values_u2_large <- readRDS("results/frequentist/p_values_u2_large.RData")
p_values_u3_large <- readRDS("results/frequentist/p_values_u3_large.RData")

p_values_large <- list(u0 = p_values_u0_large, 
                       u1 = p_values_u1_large,
                       u2 = p_values_u2_large,
                       u3 = p_values_u3_large)
  

#### 2. Comparison of log-linear and logistic approach ####

glm_diff_u0 <- abs(p_values_u0_small$log_linear - p_values_u0_small$logistic)
glm_diff_u1 <- abs(p_values_u1_small$log_linear - p_values_u1_small$logistic)
glm_diff_u2 <- abs(p_values_u2_small$log_linear - p_values_u2_small$logistic)
glm_diff_u3 <- abs(p_values_u3_small$log_linear - p_values_u3_small$logistic)

# overall mean difference in p-values
# mean(glm_diff_u0, na.rm = TRUE)
# mean(glm_diff_u1, na.rm = TRUE)
# mean(glm_diff_u2, na.rm = TRUE)
# mean(glm_diff_u3, na.rm = TRUE)

# -> practically identical, might even be rounding differences

glm_diff_u0 <- abs(p_values_u0_large$log_linear - p_values_u0_large$CMH_no_correct)
glm_diff_u1 <- abs(p_values_u1_large$log_linear - p_values_u1_large$CMH_no_correct)
glm_diff_u2 <- abs(p_values_u2_large$log_linear - p_values_u2_large$CMH_no_correct)
glm_diff_u3 <- abs(p_values_u3_large$log_linear - p_values_u3_large$CMH_no_correct)

# for the comparison using the correction
# glm_diff_u0 <- abs(p_values_u0_large$log_linear - p_values_u0_large$CMH)
# glm_diff_u1 <- abs(p_values_u1_large$log_linear - p_values_u1_large$CMH)
# glm_diff_u2 <- abs(p_values_u2_large$log_linear - p_values_u2_large$CMH)
# glm_diff_u3 <- abs(p_values_u3_large$log_linear - p_values_u3_large$CMH)


column_names_large <-  colMeans(structure(vapply(sim_data_u0_large, function(x) sum(x$arr),
                                                 numeric(1)), dim=dim(sim_data_u0_large)))

column_names_small <-  colMeans(structure(vapply(sim_data_u0_small, function(x) sum(x$arr),
                                                 numeric(1)), dim=dim(sim_data_u0_small)))


glm_diff_u0 <- as.data.frame(glm_diff_u0)
colnames(glm_diff_u0) <- column_names_large

glm_diff_df <- glm_diff_u0 %>% 
  pivot_longer(cols = everything()) %>%
  group_by(n = name) %>%
  summarise(mean_diff = mean(value, na.rm = T)) 

glm_diff_df %>%
  mutate(n = as.numeric(n)) %>%
  ggplot(aes(x = n, y = mean_diff)) + 
  geom_line() + 
  geom_point()+
  ggpubr::theme_pubr()+ 
  ylab("(Abs) difference of p-values") + 
  #scale_y_continuous(trans = "log", 
  #                   breaks = seq(0, 0.1, 0.01)) +
  xlab("Sample size") + 
  geom_hline(yintercept = 0, linetype = "dashed")

# -> difference decreases with larger sample size 
# Having established that the two approaches are almost identical, we can 
# continue with the comparison of either with the original CMH test 
# We will focus on comparing the log-linear model with the CMH test

#### 3. Comparison of log-linear and CMH test (corrected) ####

# check for uniform distribution for p-values for  data generated from the Null
cmh_ll_u0 <- as.data.frame(p_values_u0_small$log_linear) %>% 
  pivot_longer(starts_with("V")) %>%
  group_by(n = extract_numeric(name)*10) %>% 
  rename(p_value = value) %>%
  mutate(test = "Log-linear test") %>%
  bind_rows(as.data.frame(p_values_u0_small$CMH_no_correct) %>% 
              pivot_longer(starts_with("V")) %>%
              group_by(n = extract_numeric(name)*10) %>% 
              rename(p_value = value) %>%
              mutate(test = "CMH test (without correction)")) %>%
  select( -name)

p_cmh_hist <- cmh_ll_u0 %>% 
  ggplot(aes(x =p_value))+
  geom_histogram(col = "black", alpha = 0.8) + 
  facet_wrap(test~ ., nrow = 2) + 
  coord_cartesian(xlim = c(0,1)) +
  xlab("p-value")+ ylab("Count")+
  ggpubr::theme_pubr() +
  theme(strip.text = element_text(size = 10), 
        axis.title = element_text(size = 12), 
        axis.text = element_text(size = 10))

ggsave("results/frequentist/plots/histogram_u0_CMH_loglinear.pdf", 
       plot = p_cmh_hist)


#### 4. Median p-value pattern
u12 <- as.character(round(c(0, small_u12, medium_u12, large_u12), 3))

p_values_nlarge <- lapply(1:length(p_values_large), function(u) {
  u_par <- u12[u]
   p_u_df <- lapply(1:length(p_values_large[[u]]), function(t) {
     method <- names(p_values_large[[u]])[t]
     # print(method)
     p_df <- as.data.frame(p_values_large[[u]][[t]])
     colnames(p_df) <- column_names_large
     p_df_long <- p_df %>% pivot_longer(cols = everything()) %>%
       rename(n = name, p_val = value) %>% 
       mutate(method = method)
     return(p_df_long)
      }) %>% do.call(what = rbind) 
   

   p_u_df %>% mutate(u12 = u_par)
}) %>% do.call(what = rbind)


p_values_nlarge_median <- p_values_nlarge %>% 
  filter(method %in% c("CMH", "log_linear")) %>%
  mutate(method = case_when(method == "CMH_no_correct" ~
                              "CMH (without correction)", 
                            method == "CMH" ~
                              "CMH (with correction)", 
                            TRUE ~ method)) %>%
  mutate(u12 = as.character(u12)) %>%
  mutate(method = case_when(method =="log_linear" ~
                              "Log-linear",
                            TRUE ~ method)) %>%
  mutate(n = as.numeric(n)) %>%
  group_by(n, method, u12) %>%
  summarise(median_pval = median(p_val, na.rm = TRUE))

p_values_nlarge_median_plot <- p_values_nlarge_median %>%
  ggplot(aes(x = n, y = median_pval)) + 
  geom_point(aes(col = u12)) +
  geom_line(aes(group = u12, col = u12)) + 
  facet_wrap(. ~ method, nrow = 4) + 
  ggpubr::theme_pubr() + 
  theme(legend.position = "right") + 
  scale_x_continuous(limits = c(100, 1500), 
                     breaks = seq(100, 1500, 200)) +
  ylab("Median p-value") + xlab("Sample size") +
  scale_color_discrete(name = unname(TeX("u_{12}")))


ggsave("results/frequentist/plots/p_values_nlarge_median_plot.pdf", 
       plot = p_values_nlarge_median_plot)


#### analysis for small n ####  

p_values_nsmall <- lapply(1:length(p_values_small), function(u) {
  u_par <- u12[u]
  p_u_df <- lapply(1:length(p_values_small[[u]]), function(t) {
    method <- names(p_values_small[[u]])[t]
    p_df <- as.data.frame(p_values_small[[u]][[t]])
    colnames(p_df) <- column_names_small
    p_df_long <- p_df %>% pivot_longer(cols = everything()) %>%
      rename(n = name, p_val = value) %>% 
      mutate(method = method)
    return(p_df_long)
  }) %>% do.call(what = rbind) 
  
  
  p_u_df %>% mutate(u12 = u_par)
}) %>% do.call(what = rbind)


p_values_nsmall_median_plot <- p_values_nsmall %>% 
  filter(method %in% c("CMH_no_correct", "log_linear")) %>%
  mutate(method = case_when(method == "CMH_no_correct" ~
                              "CMH (without correction)", 
                            method == "CMH" ~
                              "CMH (with correction)", 
                            TRUE ~ method)) %>%
  mutate(u12 = as.character(u12)) %>%
  mutate(method = case_when(method =="log_linear" ~
                              "Log-linear",
                            TRUE ~ method)) %>%
  #filter(u12 %in% c(0, "0.363")) %>%
  mutate(n = as.numeric(n)) %>%
  group_by(n, method, u12) %>%
  summarise(median_pval = median(p_val, na.rm = TRUE)) %>%
  ggplot(aes(x = n, y = median_pval)) + 
  geom_point(aes(col = u12)) +
  geom_line(aes(group = u12, col = u12)) + 
  facet_wrap(. ~ method, nrow = 4) + 
  ggpubr::theme_pubr() + 
  theme(legend.position = "right") + 
  scale_x_continuous(limits = range(column_names_small), 
                     breaks = seq(min(column_names_small),
                                  max(column_names_small),
                                      10))+
  ylab("Median p-value") + xlab("Sample size") +
  scale_color_discrete(name = unname(TeX("u_{12}")))

ggsave("results/frequentist/plots/p_values_nsmall_median_plot.pdf", 
       plot = p_values_nsmall_median_plot)



######## Agreement matrix ############
p_val_table <- table(p_values_u0_large$log_linear < 0.05, p_values_u0_large$CMH < 0.05)
round(p_val_table/sum(p_val_table), 3)





