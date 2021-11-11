#### install and load list of packages #### 

required_packages <- c("dplyr",
                       "tidyr",
                       "stringr", 
                       "ggplot2",
                       "ggpubr",
                       "latex2exp",
                       "BAS",
                       "readr", 
                       "stargazer",
                       "survival",
                       "aod",
                       "pbapply",
                       "numDeriv")

installed_packages <- installed.packages()[,"Package"]

# check if required packages are installed, 
# if not, install them 
sapply(required_packages, function(p) {
  if(!p %in% installed_packages) {
    install.packages(p)
  }
})

# load packages
lapply(required_packages, require, character.only = TRUE)

