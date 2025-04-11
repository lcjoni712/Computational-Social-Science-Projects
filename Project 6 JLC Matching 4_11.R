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
#Your Answer: The histogram of covariate balance for "student_vote" across random assignments demonstrates that while most simulations show small balance differences (0-0.03), some randomizations produce larger imbalances. This illustrates the key distinction between independence (guaranteed by random assignment) and perfect balance (not guaranteed in finite samples). Independence ensures the probability of treatment assignment doesn't depend on covariates, but any single randomization can still produce chance imbalances in the actual distribution of covariates across treatment groups. The right-skewed distribution visualizes this sampling variability, showing why researchers might employ techniques like stratified randomization when perfect balance is critical. This explains why independence of treatment assignment and baseline covariates conceptually does not guarantee balance between them in practice.

##Propensity Matching
library(MatchIt)
library(cobalt)
library(tidyverse)

# Load your data
df <- read.csv("C:/Users/lande/Downloads/ypsps.csv")

# Define treatment and outcome
treatment_var <- "college"
outcome_var <- "student_ppnscal"

# Drop rows with NA in treatment or outcome
df <- df %>%
  filter(!is.na(.data[[treatment_var]]), !is.na(.data[[outcome_var]]))

# Select potential covariates (drop post-treatment variables manually if known)
post_treatment_vars <- c("student_ppnscal")  # Add more if needed
all_covariates <- setdiff(names(df), c(treatment_var, outcome_var, post_treatment_vars))

# Keep only numeric or factor covariates
all_covariates <- all_covariates[sapply(df[all_covariates], function(x) is.numeric(x) || is.factor(x))]

# Run simulations
set.seed(42)
n_sim <- 25
results_list <- list()

for (i in 1:n_sim) {
  n_covs <- sample(2:min(6, length(all_covariates)), 1)  # select 2–6 covariates
  selected_covs <- sample(all_covariates, n_covs)
  fmla <- as.formula(paste(treatment_var, "~", paste(selected_covs, collapse = " + ")))
  
  m.out <- tryCatch({
    matchit(fmla, data = df, method = "nearest")
  }, error = function(e) NULL)
  
  if (is.null(m.out)) next
  
  m.data <- match.data(m.out)
  if (nrow(m.data) == 0) next
  
  att <- with(m.data, mean(get(outcome_var)[get(treatment_var) == 1]) -
                mean(get(outcome_var)[get(treatment_var) == 0]))
  
  balance_stats <- tryCatch({
    bal.tab(m.out, un = TRUE, disp.v.ratio = FALSE)
  }, error = function(e) NULL)
  
  if (is.null(balance_stats)) next
  
  smds <- balance_stats$Balance$Diff.Adj
  smds_un <- balance_stats$Balance$Diff.Un
  prop_balanced <- mean(abs(smds) <= 0.1, na.rm = TRUE)
  mean_improvement <- mean(abs(smds_un) - abs(smds), na.rm = TRUE)
  
  results_list[[length(results_list) + 1]] <- list(
    formula = fmla,
    att = att,
    prop_balanced = prop_balanced,
    mean_improvement = mean_improvement
  )
}

# Build results_df
results_df <- purrr::map_dfr(results_list, ~tibble(
  formula = deparse(.x$formula),
  att = .x$att,
  prop_balanced = .x$prop_balanced,
  mean_improvement = .x$mean_improvement
))

# View
print(head(results_df))

# Plot
ggplot(results_df, aes(x = prop_balanced, y = att)) +
  geom_point(alpha = 0.6) +
  labs(x = "Proportion Balanced (|SMD| ≤ 0.1)", y = "ATT",
       title = "ATT vs Covariate Balance Proportion") +
  theme_minimal()


##4.3 Questions
##1. How many simulations resulted in models with a higher proportion of balanced covariates? Do you have any concerns about this?
#Very few simulations resulted in models with a higher proportion of balanced covariates. Most points cluster at the left side of the x-axis (near 0.00), with only a few points at 0.15 proportion balanced. This is concerning because it suggests that across most simulations, very few covariates achieved balance (defined as |SMD| ≤ 0.1), indicating potential systematic imbalance issues in the study design.
##2. Analyze the distribution of the ATTs. Do you have any concerns about this distribution?
#The ATT (Average Treatment Effect on the Treated) values appear to vary widely, ranging from approximately 1.3 to 2.2, despite most simulations having similar (poor) balance proportions. This variation in treatment effect estimates across simulations with similar balance characteristics is concerning, as it suggests the treatment effect estimates are unstable and potentially biased due to the covariate imbalances.

##Matching Algorithm of Your Choice
# Load necessary libraries
library(tidyverse)
library(MatchIt)
library(cobalt)
library(Matching)  # Explicitly load Matching package
library(gridExtra)

# Read data
data <- read_csv("C:/Users/lande/Downloads/ypsps.csv")  

# Define baseline variables (excluding post-treatment variables)
baseline_vars <- names(data)[!str_detect(names(data), "1973|1976|1979|1980")]

# Define treatment and outcome variables
treatment_var <- "college"
outcome_var <- "student_ppnscal"

# Prepare data - handle missing values more carefully
# Fix for select() error
data_base <- data[, baseline_vars]  # Base R subsetting instead of dplyr select

# Remove post-treatment variables from potential covariates
post_treatment_vars <- c("student_ppnscal")  # Add more if needed
all_covariates <- setdiff(baseline_vars, c(treatment_var, outcome_var, post_treatment_vars))

# Filter problematic variables identified in the error log
problematic_vars <- c(
  "parent_HHCollegePlacebo", 
  "student_1982College", 
  "student_1982IncSelf", 
  "parent_GPHighSchoolPlacebo"
)

# Remove problematic variables from covariates list
all_covariates <- setdiff(all_covariates, problematic_vars)

# Keep only numeric or factor covariates with no missing values
valid_covariates <- c()
for (cov in all_covariates) {
  # Check if numeric or factor AND has no missing values
  if ((is.numeric(data_base[[cov]]) || is.factor(data_base[[cov]])) && 
      sum(is.na(data_base[[cov]])) == 0) {
    valid_covariates <- c(valid_covariates, cov)
  }
}

# Final preparation of data - drop NAs in treatment and outcome only
data_clean <- data_base[!is.na(data_base[[treatment_var]]) & !is.na(data_base[[outcome_var]]), ]

# Run genetic matching simulations
set.seed(42)
n_sim <- 25
genetic_results_list <- list()
genetic_models <- list()

for (i in 1:n_sim) {
  # Randomly select covariates from valid covariates
  n_covs <- sample(2:min(6, length(valid_covariates)), 1)  # select 2-6 covariates
  selected_covs <- sample(valid_covariates, n_covs)
  fmla <- as.formula(paste(treatment_var, "~", paste(selected_covs, collapse = " + ")))
  
  cat("Running simulation", i, "with formula:", deparse(fmla), "\n")
  
  # Create a temporary dataset with only the needed variables and no missing values
  temp_cols <- c(treatment_var, outcome_var, selected_covs)
  temp_data <- data_clean[, temp_cols]
  # Remove any rows with NA
  temp_data <- temp_data[complete.cases(temp_data), ]
  
  # Run genetic matching with error handling
  m.out <- tryCatch({
    matchit(fmla, data = temp_data, method = "genetic", estimand = "ATT",
            # Add parameters to address warnings
            pop.size = 100,  # Increase population size
            max.generations = 100,  # Limit generations to prevent long runs
            wait.generations = 10)  # Wait generations
  }, error = function(e) {
    cat("Error in genetic matching:", e$message, "\n")
    NULL
  })
  
  # If matching failed, skip to next iteration
  if (is.null(m.out)) {
    cat("Skipping simulation", i, "due to error\n")
    next
  }
  
  # Store the model
  genetic_models[[i]] <- m.out
  
  # Get matched data
  m.data <- match.data(m.out)
  if (nrow(m.data) == 0) {
    cat("Skipping simulation", i, "due to empty matched dataset\n")
    next
  }
  
  # Calculate ATT
  m_formula <- as.formula(paste(outcome_var, "~", treatment_var))
  model_result <- lm(m_formula, data = m.data, weights = weights)
  att <- coef(model_result)[treatment_var]
  
  # Get balance statistics
  balance_stats <- tryCatch({
    bal.tab(m.out, un = TRUE, disp.v.ratio = FALSE)
  }, error = function(e) {
    cat("Error in calculating balance:", e$message, "\n")
    NULL
  })
  
  if (is.null(balance_stats)) {
    cat("Skipping simulation", i, "due to balance calculation error\n")
    next
  }
  
  # Calculate balance metrics
  smds <- balance_stats$Balance$Diff.Adj
  smds_un <- balance_stats$Balance$Diff.Un
  prop_balanced <- mean(abs(smds) <= 0.1, na.rm = TRUE)
  
  # Handle case where all SMDs are NA
  if (all(is.na(smds_un)) || all(is.na(smds))) {
    mean_improvement <- NA
  } else {
    # Only calculate for non-NA values and avoid division by zero
    valid_indices <- !is.na(smds_un) & !is.na(smds) & abs(smds_un) > 0
    if (sum(valid_indices) > 0) {
      improvements <- (abs(smds_un[valid_indices]) - abs(smds[valid_indices])) / abs(smds_un[valid_indices])
      mean_improvement <- mean(improvements, na.rm = TRUE) * 100
    } else {
      mean_improvement <- NA
    }
  }
  
  # Store results
  genetic_results_list[[length(genetic_results_list) + 1]] <- list(
    sim_id = i,
    formula = fmla,
    att = att,
    prop_balanced = prop_balanced,
    mean_improvement = mean_improvement,
    covariates = paste(selected_covs, collapse = ", ")
  )
}

# Build results dataframe
genetic_results_df <- do.call(rbind, lapply(genetic_results_list, function(x) {
  data.frame(
    sim_id = x$sim_id,
    formula = deparse(x$formula),
    att = x$att,
    prop_balanced = x$prop_balanced,
    mean_improvement = x$mean_improvement,
    covariates = x$covariates
  )
}))

# View results
print(head(genetic_results_df))

# Plot ATT vs proportion of balanced covariates
library(ggplot2)
ggplot(genetic_results_df, aes(x = prop_balanced, y = att)) +
  geom_point(alpha = 0.8) +
  labs(x = "Proportion Balanced (|SMD| ≤ 0.1)", y = "ATT",
       title = "ATT vs Covariate Balance Proportion (Genetic Matching)") +
  theme_minimal()

# Plot distribution of ATTs
ggplot(genetic_results_df, aes(x = att)) +
  geom_histogram(bins = 10, fill = "steelblue", color = "black", alpha = 0.7) +
  labs(title = "Distribution of ATT Estimates (Genetic Matching)",
       x = "ATT", y = "Count") +
  theme_minimal()

# Plot mean improvement
# Filter out NA values in mean_improvement
genetic_results_df_filtered <- genetic_results_df[!is.na(genetic_results_df$mean_improvement), ]

if(nrow(genetic_results_df_filtered) > 0) {
  ggplot(genetic_results_df_filtered, aes(x = factor(sim_id), y = mean_improvement)) +
    geom_bar(stat = "identity", fill = "darkgreen", alpha = 0.7) +
    labs(title = "Mean Percent Improvement in Balance by Simulation",
         x = "Simulation ID", y = "Mean % Improvement") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90))
}

# Select 5 random models for balance plots
non_null_models <- which(!sapply(genetic_models, is.null))
if (length(non_null_models) >= 5) {
  selected_indices <- sample(non_null_models, 5)
} else if (length(non_null_models) > 0) {
  selected_indices <- non_null_models
} else {
  selected_indices <- integer(0)
}

# Create balance plots
balance_plots <- list()
for (i in seq_along(selected_indices)) {
  idx <- selected_indices[i]
  if (length(idx) > 0) {
    plot <- tryCatch({
      love.plot(
        genetic_models[[idx]],
        stats = "mean.diffs",
        abs = TRUE,
        thresholds = c(0.1),
        var.order = "unadjusted",
        line = TRUE,
        stars = "std", # Added stars parameter to distinguish between SMD and raw differences
        title = paste("Model", idx, "Balance")
      )
    }, error = function(e) {
      cat("Error creating plot for model", idx, ":", e$message, "\n")
      NULL
    })
    
    if (!is.null(plot)) {
      balance_plots[[i]] <- plot
    }
  }
}

# Arrange balance plots in a grid
if (length(balance_plots) > 0) {
  grid.arrange(grobs = balance_plots, ncol = min(2, length(balance_plots)))
} else {
  cat("No valid models for balance plots\n")
}

# Compare with previous propensity score results (if available)
if (exists("results_df")) {
  # Check column names of both dataframes
  genetic_cols <- names(genetic_results_df)
  ps_cols <- names(results_df)
  
  # Print column names for debugging
  cat("Genetic columns:", paste(genetic_cols, collapse=", "), "\n")
  cat("Propensity score columns:", paste(ps_cols, collapse=", "), "\n")
  
  # Find common columns that are essential for the plot
  required_cols <- c("prop_balanced", "att")
  
  # Check if both dataframes have the required columns
  if (all(required_cols %in% genetic_cols) && all(required_cols %in% ps_cols)) {
    # Select only the columns that exist in both dataframes
    common_cols <- intersect(genetic_cols, ps_cols)
    
    # Create subsets with common columns
    genetic_df_for_combine <- genetic_results_df[, common_cols]
    ps_df_for_combine <- results_df[, common_cols]
    
    # Add method column to each
    genetic_df_for_combine$method <- "Genetic Matching"
    ps_df_for_combine$method <- "Propensity Score"
    
    # Now combine the dataframes with the same structure
    combined_df <- rbind(ps_df_for_combine, genetic_df_for_combine)
    
    # Create the comparison plot
    ggplot(combined_df, aes(x = prop_balanced, y = att, color = method)) +
      geom_point(alpha = 0.7) +
      labs(x = "Proportion Balanced (|SMD| ≤ 0.1)", y = "ATT",
           title = "ATT vs Covariate Balance: Propensity Score vs Genetic Matching") +
      theme_minimal()
  } else {
    cat("Cannot create comparison plot: required columns not found in both dataframes\n")
    cat("Required columns:", paste(required_cols, collapse=", "), "\n")
  }
}

##5.2 Questions
#1. Does your alternative matching method have more runs with higher proportions of balanced covariates?
# It appears that the propensity score method tends to have higher proportions of balanced covariates compared to genetic matching. The propensity score points are generally shifted further to the right on the x-axis (Proportion Balanced |SMD| ≤ 0.1), indicating better covariate balance.
#2. Use a visualization to examine the change in the distribution of the percent improvement in balance in propensity score matching vs. the distribution of the percent improvement in balance in your new method. Which did better? Analyze the results in 1-2 sentences.
#The propensity score matching appears to achieve better covariate balance overall compared to genetic matching. The resulsts suggest that the propensity score matching consistently balances a higher proportion of covariates while maintaining comparable ATT estimates to genetic matching, which indicates greater reliability in the causal effect estimation.


##Discussion Questions
#1. Why might it be a good idea to do matching even if we have a randomized or as-if-random design?
#Matching can be beneficial even in randomized or as-if-random designs for several reasons. It can correct for chance imbalances that occur despite randomization, especially in smaller samples where perfect balance is rarely achieved. Matching increases precision and statistical power by reducing variance in treatment effect estimates, while addressing issues from differential attrition or non-compliance that may emerge after initial randomization. Additionally, it provides robustness against model misspecification when estimating treatment effects and can help handle missing data problems that might arise during the study.
#2. The standard way of estimating the propensity score is using a logistic regression to estimate probability of treatment. Given what we know about the curse of dimensionality, do you think there might be advantages to using other machine learning algorithms (decision trees, bagging/boosting forests, ensembles, etc.) to estimate propensity scores instead?
#There are definite advantages to using machine learning algorithms instead of standard logistic regression for propensity score estimation. ML methods can automatically detect complex non-linear relationships and interactions without requiring manual specification, handling high-dimensional data more effectively and avoiding the curse of dimensionality that hampers logistic regression. Ensemble methods like random forests or boosting produce more stable propensity scores by averaging across multiple models, often yielding better balance on covariates as they capture more complex treatment assignment mechanisms. These approaches are less prone to extreme propensity scores in regions of sparse data, which leads to more reliable matching and ultimately more credible causal estimates.
