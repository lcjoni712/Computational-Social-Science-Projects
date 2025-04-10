install.packages(c("tidyverse", "MatchIt", "cobalt", "optmatch", "gridExtra"))

library(tidyverse)
library(MatchIt)
library(cobalt)
library(optmatch)
library(gridExtra)
data <- read_csv("C:/Users/lande/Downloads/ypsps.csv")
baseline_vars <- names(data)[!str_detect(names(data), "1973|1976|1979|1980")]
names(data)
data_base <- data %>%
  select(all_of(baseline_vars)) %>%
  drop_na(college, student_ppnscal) 
set.seed(123)
covariate_name <- "student_vote"
head(data_base)  
colnames(data_base)  
table(data_base[[covariate_name]])
set.seed(123)
get_simple_balance <- function(data, covariate_name) {
  tryCatch({
    
    treat <- sample(c(0, 1), nrow(data), replace = TRUE)
    
    
    if(length(unique(data[[covariate_name]])) <= 2) {
      prop_treat <- mean(data[[covariate_name]][treat == 1], na.rm = TRUE)
      prop_control <- mean(data[[covariate_name]][treat == 0], na.rm = TRUE)
      return(abs(prop_treat - prop_control))
    } 
    
    else {
      mean_treat <- mean(data[[covariate_name]][treat == 1], na.rm = TRUE)
      mean_control <- mean(data[[covariate_name]][treat == 0], na.rm = TRUE)
      sd_all <- sd(data[[covariate_name]], na.rm = TRUE)
      return(abs(mean_treat - mean_control) / sd_all)
    }
  }, error = function(e) {
    message("Error in simulation: ", e$message)
    return(NA)
  })}
sim_results <- replicate(1000, get_simple_balance(data_base, covariate_name))
sim_results <- sim_results[!is.na(sim_results)]
summary(sim_results)
sim_df <- data.frame(balance_diff = sim_results)
hist(sim_results, 
     main = "Covariate Balance Simulation", 
     xlab = paste("Balance Difference in", covariate_name),
     col = "green",
     breaks = 30)
library(ggplot2)
ggplot(sim_df, aes(x = balance_diff)) +
  geom_histogram(bins = 30, fill = "green", color = "white") +
  labs(
    title = "Covariate Balance Across Random Assignments",
    x = paste("Balance Difference in", covariate_name),
    y = "Count"
  )
## Questions
#What do you see across your simulations? Why does independence of treatment assignment and baseline covariates not guarantee balance of treatment assignment and baseline covariates?
#Your Answer:





