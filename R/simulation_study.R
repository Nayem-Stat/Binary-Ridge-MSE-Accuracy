# Clear environment
rm(list = ls())

# Load necessary libraries
library(MASS)
library(glmnet)
library(pls)
library(caret)
library(openxlsx) # For writing results to an Excel file
library(psych)

# Set parameters
sample_sizes <- c(100, 200, 300) # Different sample sizes
p_values <- c(3, 5, 10) # Number of predictors
rho_values <- c(0.80, 0.90, 0.95, 0.99) # Different collinearity levels
sigma_squared_values <- c(1,25) # Error variance
num_simulations <- 5000 # Number of simulations
lambda_names <- c("k1", "k2", "k3", "k4", "k5", "k6", "k7", "k8", "k9", "k10", "k11", "k12", "k13",
                  "k14", "k15", "k16", "k17", "k18", "k19", "k20", "k21", "k22", "k23") # New lambda (k) names
model_names <- c(lambda_names, "OLS", "Lasso", "ElasticNet", "OptimalRidge")

set.seed(12345678) # For reproducibility

# Generate X
generate_X <- function(n, p, rho) {
  Z <- matrix(rnorm(n * p), nrow = n, ncol = p)
  X <- sqrt(1 - rho^2) * Z + rho * matrix(rnorm(n), nrow = n, ncol = p)
  return(X)
}

# Generate binary y (using logistic regression link function)
generate_y_binary <- function(X, beta, sigma_squared) {
  n <- nrow(X)
  linear_predictor <- X %*% beta
  probabilities <- 1 / (1 + exp(-linear_predictor)) # Logistic link
  y <- rbinom(n, 1, probabilities) # Generate binary outcomes
  return(y)
}

# Function to compute different lambda values
compute_lambda <- function(k, sigma_squared, eigenvalues, beta_hat, beta_hat_squared, n, p) {
  q_i <- max(eigenvalues) / (((n - p) * sigma_squared) + (max(eigenvalues) * beta_hat_squared))

  lambda <- switch(k,
                   "k1" = max(1 / q_i),
                   "k2" = max((sigma_squared / beta_hat)),
                   "k3" = min(sigma_squared / ((sigma_squared / beta_hat_squared) + (1 / eigenvalues))),
                   "k4" = sqrt(p * sum(max(eigenvalues) * sigma_squared / ((n - p) * sigma_squared + max(eigenvalues) * beta_hat_squared))),
                   "k5" = median((sigma_squared / beta_hat_squared) + (1 / eigenvalues)),
                   "k6" = min((sigma_squared / beta_hat_squared) + (1 / eigenvalues)),
                   "k7" = max((sigma_squared / beta_hat_squared) + (1 / eigenvalues)),
                   "k8" = exp(mean(log(((sigma_squared / beta_hat_squared) + (1 / eigenvalues))))),
                   "k9" = harmonic.mean(((sigma_squared / beta_hat_squared) + (1 / eigenvalues))),
                   "k10" = median(max(eigenvalues) * sigma_squared / ((n - p) * sigma_squared + eigenvalues * beta_hat_squared)),
                   "k11" = max(max(eigenvalues) * sigma_squared / ((n - p) * sigma_squared + max(eigenvalues) * beta_hat_squared)),
                   "k12" = sqrt(p * sum(eigenvalues * sigma_squared / ((n - p) * sigma_squared + eigenvalues * beta_hat_squared))),
                   "k13" = sum(max((eigenvalues) * sigma_squared) / ((n - p) * sigma_squared + eigenvalues * beta_hat_squared)) / p,
                   "k14" = max(1 / ((eigenvalues * sigma_squared) / ((n - p) * sigma_squared + eigenvalues * beta_hat_squared))),
                   "k15" = prod((max(eigenvalues) * sigma_squared) / ((n - p) * sigma_squared + eigenvalues * beta_hat_squared))^(1/p),
                   "k16" = max(1 / sqrt((eigenvalues * sigma_squared) / ((n - p) * sigma_squared + max(eigenvalues) * beta_hat_squared))),
                   "k17" = max(1 / ((eigenvalues * sigma_squared) / ((n - p) * sigma_squared + eigenvalues * max(beta_hat_squared)))),
                   "k18" = prod(1 / q_i)^(1/p),
                   "k19" = max(1 / ((eigenvalues * sigma_squared) / ((n - p) * sigma_squared + eigenvalues * max(beta_hat_squared)))),
                   "k20" = 1 / (sqrt((2 * sigma_squared) / (max(eigenvalues) * max(beta_hat_squared)))),
                   "k21" = 1 / (sqrt((2 * p * sigma_squared) / (max(eigenvalues) * sum(beta_hat_squared)))),
                   "k22" = (sqrt(sigma_squared)) * max((p^(1+(p/n))), (p^(1+(1/p)))),
                   "k23"= (sqrt(sigma_squared)) * (p^(1+(p/n))),
                   stop("Unsupported k value")
  )

  # Ensure lambda is non-negative
  lambda <- ifelse(lambda < 0, abs(lambda), lambda)

  return(lambda)
}

# Function to calculate AMSE (Corrected)
calc_AMSE <- function(beta_hat_values, beta) {
  mse_values <- apply(beta_hat_values, 1, function(beta_hat) {
    mean((beta_hat - beta)^2)  # Simplified and more accurate
  })
  mean_mse <- mean(mse_values)
  return(mean_mse)
}


# Function to calculate accuracy
calc_accuracy <- function(y_true, y_pred) {
  mean(y_true == y_pred)
}

# Create a list to store all results
results_all <- list()

# Loop over each combination of p, rho, and sigma_squared
for (p_idx in seq_along(p_values)) {
  p <- p_values[p_idx]

  for (rho in rho_values) {
    for (sigma_squared in sigma_squared_values) {

      # Run simulations for each sample size
      for (n in sample_sizes) {
        AMSE_values <- setNames(vector("list", length(model_names)), model_names)
        accuracy_values <- setNames(vector("list", length(model_names)), model_names)
        lambda_values <- setNames(vector("list", length(lambda_names)), lambda_names)

        # Initialize the AMSE_values and accuracy_values lists with matrices for each model
        for (model_name in model_names) {
          AMSE_values[[model_name]] <- matrix(0, nrow = num_simulations, ncol = p) # Use p directly
          accuracy_values[[model_name]] <- vector("numeric", length = num_simulations)
        }

        for (sim in 1:num_simulations) {
          X <- generate_X(n, p, rho)
          eigen_decomp <- eigen(t(X) %*% X)
          beta <- eigen_decomp$vectors[, 1] # True beta (from eigenvector)
          beta <- beta / sqrt(sum(beta^2))
          y <- generate_y_binary(X, beta, sigma_squared)

          # Split data into training and testing sets
          train_index <- sample(1:n, 0.7 * n)
          X_train <- X[train_index, ]
          y_train <- y[train_index]
          X_test <- X[-train_index, ]
          y_test <- y[-train_index]


          # Convert to data frames:  <--- This is the crucial change
          X_train_df <- as.data.frame(X_train)
          X_test_df <- as.data.frame(X_test)

          # Fit OLS model
          model_ols <- glm(y_train ~ ., data = X_train_df, family = binomial) # Use data frame
          beta_hat_ols <- coef(model_ols)[-1] # Remove intercept
          AMSE_values[["OLS"]][sim, ] <- beta_hat_ols
          y_pred_ols <- predict(model_ols, newdata = X_test_df, type = "response") # Use data frame
          y_pred_ols_class <- ifelse(y_pred_ols > 0.5, 1, 0)
          accuracy_values[["OLS"]][sim] <- calc_accuracy(y_test, y_pred_ols_class)

          # Fit initial ridge regression model with a small lambda
          initial_lambda <- 0.01
          model_ridge <- glmnet(X_train, y_train, alpha = 0, lambda = initial_lambda, family = "binomial", intercept = FALSE)
          beta_hat <- as.vector(coef(model_ridge, s = initial_lambda))[-1] # Use a small lambda to get initial betas
          beta_hat_squared <- beta_hat^2

          eigenvalues <- eigen(cor(X))$values

          for (lambda_name in lambda_names) {
            lambda <- compute_lambda(lambda_name, sigma_squared, eigenvalues, beta_hat, beta_hat_squared, n, p)
            # Ensure lambda is non-negative and valid
            if (lambda < 0 || is.na(lambda) || is.infinite(lambda)) {
              lambda <- abs(lambda)
              cat("Warning: Negative or invalid lambda value for", lambda_name, "in simulation", sim, "\n")
            }
            lambda_values[[lambda_name]] <- lambda # Store the computed lambda value
            model_ridge <- glmnet(X_train, y_train, alpha = 0, lambda = lambda, family = "binomial", intercept = FALSE)
            beta_hat <- as.vector(coef(model_ridge, s = lambda))[-1]
            AMSE_values[[lambda_name]][sim, ] <- beta_hat
            y_pred_ridge <- predict(model_ridge, newx = X_test, s = lambda, type = "response") # glmnet can use matrix
            y_pred_ridge_class <- ifelse(y_pred_ridge > 0.5, 1, 0)
            accuracy_values[[lambda_name]][sim] <- calc_accuracy(y_test, y_pred_ridge_class)

          }

          # Lasso
          model_lasso <- cv.glmnet(X_train, y_train, alpha = 1, family = "binomial", intercept = FALSE)
          beta_hat_lasso <- as.vector(coef(model_lasso, s = model_lasso$lambda.min))[-1]
          AMSE_values[["Lasso"]][sim, ] <- beta_hat_lasso
          y_pred_lasso <- predict(model_lasso, newx = X_test, s = model_lasso$lambda.min, type = "response") # glmnet can use matrix
          y_pred_lasso_class <- ifelse(y_pred_lasso > 0.5, 1, 0)
          accuracy_values[["Lasso"]][sim] <- calc_accuracy(y_test, y_pred_lasso_class)


          # ElasticNet (0.5)
          model_enet <- cv.glmnet(X_train, y_train, alpha = 0.5, family = "binomial", intercept = FALSE)
          beta_hat_enet <- as.vector(coef(model_enet, s = model_enet$lambda.min))[-1]
          AMSE_values[["ElasticNet"]][sim, ] <- beta_hat_enet
          y_pred_enet <- predict(model_enet, newx = X_test, s = model_enet$lambda.min, type = "response")
          y_pred_enet_class <- ifelse(y_pred_enet > 0.5, 1, 0)
          accuracy_values[["ElasticNet"]][sim] <- calc_accuracy(y_test, y_pred_enet_class)

          # Optimal Ridge (using CV)
          cv_ridge <- cv.glmnet(X_train, y_train, alpha = 0, family = "binomial", intercept = FALSE)
          beta_hat_opt_ridge <- as.vector(coef(cv_ridge, s = cv_ridge$lambda.min))[-1]
          AMSE_values[["OptimalRidge"]][sim, ] <- beta_hat_opt_ridge
          y_pred_opt_ridge <- predict(cv_ridge, newx = X_test, s = cv_ridge$lambda.min, type = "response")
          y_pred_opt_ridge_class <- ifelse(y_pred_opt_ridge > 0.5, 1, 0)
          accuracy_values[["OptimalRidge"]][sim] <- calc_accuracy(y_test, y_pred_opt_ridge_class)

        }

        # Store results for the current sample size and sigma_squared
        results_all[[paste("n =", n, "_p =", p, "_rho =", rho, "_sigma_squared =", sigma_squared)]] <- list(AMSE_values, accuracy_values, lambda_values)
      }
    }
  }
}


# Aggregate and store results in data frames
all_results_df <- data.frame()


for (p_idx in seq_along(p_values)) {
  p <- p_values[p_idx]
  for (rho in rho_values) {
    for (sigma_squared in sigma_squared_values) {
      for (n in sample_sizes) {
        results <- results_all[[paste("n =", n, "_p =", p, "_rho =", rho, "_sigma_squared =", sigma_squared)]]
        AMSE_values <- results[[1]]
        accuracy_values <- results[[2]]

        for (model_name in model_names) {
          mean_amse <- calc_AMSE(AMSE_values[[model_name]], beta) # Use the correct beta
          # ... (rest of your aggregation code)

          mean_accuracy <- mean(accuracy_values[[model_name]])

          all_results_df <- rbind(all_results_df, data.frame(
            Model = model_name,
            P = p,
            Rho = rho,
            n = n,
            Sigma_squared = sigma_squared,
            Mean_AMSE = mean_amse,
            Mean_Accuracy = mean_accuracy
          ))
        }
      }
    }
  }
}

write.csv(all_results_df, "C:\\Users\\hnayem\\OneDrive - Florida International University\\Documents\\Private\\Phd work\\Folder for Simulation work\\Binary Ridge - 4rth Paper\\ridge_mse_accuarcy21925.csv", row.names = FALSE)
