# Clear environment
rm(list = ls())

# Load necessary libraries
library(MASS)
library(glmnet)
library(pls)
library(caret)
library(openxlsx)
library(psych)

# Load real data from CSV
real_data <- read.csv("C:\\Users\\hnayem\\OneDrive - Florida International University\\Documents\\Private\\Phd work\\Folder for Simulation work\\Binary Ridge - 4rth Paper\\cancer_data.csv")

str(real_data)

# Extract outcome (DEATH_EVENT) and predictors
y <- real_data$remiss
X <- real_data[, names(real_data) != "DEATH_EVENT"] # Exclude the outcome

# Convert X to a matrix
X <- as.matrix(X)


cor(X)

# Set parameters
n <- nrow(X)
p <- ncol(X)
sigma_squared <- var(y) # Estimate sigma_squared from the data variance.

lambda_names <- c("k1", "k2", "k3", "k4", "k5", "k6", "k7", "k8", "k9", "k10", "k11", "k12", "k13",
                  "k14", "k15", "k16", "k17", "k18", "k19", "k20", "k21", "k22", "k23")
model_names <- c(lambda_names, "OLS", "Lasso", "ElasticNet", "OptimalRidge")

set.seed(12345678) # For reproducibility

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
                   "k17" = max(1 / ((eigenvalues) * sigma_squared / ((n - p) * sigma_squared + eigenvalues * max(beta_hat_squared)))),
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

# Function to calculate accuracy
calc_accuracy <- function(y_true, y_pred) {
  mean(y_true == y_pred)
}

# Function to calculate MSE (beta vs beta_hat)
calc_beta_mse <- function(beta_true, beta_hat) {
  mean((beta_true - beta_hat)^2)
}

# Estimate true beta
eigen_decomp <- eigen(t(X) %*% X)
beta_true <- eigen_decomp$vectors[, 1] # True beta (from eigenvector)
beta_true <- beta_true / sqrt(sum(beta_true^2))

# Split data into training and testing sets (70% train, 30% test)
train_index <- sample(1:n, 0.8 * n)
X_train <- X[train_index, ]
y_train <- y[train_index]
X_test <- X[-train_index, ]
y_test <- y[-train_index]

# Convert to data frames
X_train_df <- as.data.frame(X_train)
X_test_df <- as.data.frame(X_test)

# Fit OLS model
model_ols <- glm(y_train ~ ., data = X_train_df, family = binomial)
beta_hat_ols <- coef(model_ols)[-1] # Remove intercept

# Fit initial ridge regression model with a small lambda
initial_lambda <- 0.01
model_ridge <- glmnet(X_train, y_train, alpha = 0, lambda = initial_lambda, family = "binomial", intercept = FALSE)
beta_hat <- as.vector(coef(model_ridge, s = initial_lambda))[-1]
beta_hat_squared <- beta_hat^2

eigenvalues <- eigen(cor(X))$values

mse_ridge <- setNames(vector("numeric", length(lambda_names)), lambda_names)
beta_hats_ridge <- setNames(vector("list", length(lambda_names)), lambda_names) #store beta_hats

for (lambda_name in lambda_names) {
  lambda <- compute_lambda(lambda_name, sigma_squared, eigenvalues, beta_hat, beta_hat_squared, n, p)
  if (lambda < 0 || is.na(lambda) || is.infinite(lambda)) {
    lambda <- abs(lambda)
    cat("Warning: Negative or invalid lambda value for", lambda_name, "\n")
  }
  model_ridge <- glmnet(X_train, y_train, alpha = 0, lambda = lambda, family = "binomial", intercept = FALSE)
  beta_hat_ridge <- as.vector(coef(model_ridge, s = lambda))[-1]
  beta_hats_ridge[[lambda_name]] <- beta_hat_ridge #store beta_hat
}

# Lasso
model_lasso <- cv.glmnet(X_train, y_train, alpha = 1, family = "binomial", intercept = FALSE)
beta_hat_lasso <- as.vector(coef(model_lasso, s = model_lasso$lambda.min))[-1]

# ElasticNet (0.5)
model_enet <- cv.glmnet(X_train, y_train, alpha = 0.5, family = "binomial", intercept = FALSE)
beta_hat_enet <- as.vector(coef(model_enet, s = model_enet$lambda.min))[-1]

# Optimal Ridge (using CV)
cv_ridge <- cv.glmnet(X_train, y_train, alpha = 0, family = "binomial", intercept = FALSE)
beta_hat_opt_ridge <- as.vector(coef(cv_ridge, s = cv_ridge$lambda.min))[-1]

# Calculate MSEs (beta vs beta_hat)
mse_ols <- calc_beta_mse(beta_true, beta_hat_ols)
for (lambda_name in lambda_names) {
  mse_ridge[lambda_name] <- calc_beta_mse(beta_true, beta_hats_ridge[[lambda_name]])
}
mse_lasso <- calc_beta_mse(beta_true, beta_hat_lasso)
mse_enet <- calc_beta_mse(beta_true, beta_hat_enet)
mse_opt_ridge <- calc_beta_mse(beta_true, beta_hat_opt_ridge)

# Calculate accuracy
y_pred_ols <- predict(model_ols, newdata = X_test_df, type = "response")
y_pred_ols_class <- ifelse(y_pred_ols > 0.5, 1, 0)
accuracy_ols <- calc_accuracy(y_test, y_pred_ols_class)

accuracy_ridge <- setNames(vector("numeric", length(lambda_names)), lambda_names)
for (lambda_name in lambda_names) {
  model_ridge <- glmnet(X_train, y_train, alpha = 0, lambda = compute_lambda(lambda_name, sigma_squared, eigenvalues, beta_hat, beta_hat_squared, n, p), family = "binomial", intercept = FALSE)
  y_pred_ridge <- predict(model_ridge, newx = X_test, s = compute_lambda(lambda_name, sigma_squared, eigenvalues, beta_hat, beta_hat_squared, n, p), type = "response")
  y_pred_ridge_class <- ifelse(y_pred_ridge > 0.5, 1, 0)
  accuracy_ridge[lambda_name] <- calc_accuracy(y_test, y_pred_ridge_class)
}

y_pred_lasso <- predict(model_lasso, newx = X_test, s = model_lasso$lambda.min, type = "response")
y_pred_lasso_class <- ifelse(y_pred_lasso > 0.5, 1, 0)
accuracy_lasso <- calc_accuracy(y_test, y_pred_lasso_class)

y_pred_enet <- predict(model_enet, newx = X_test, s = model_enet$lambda.min, type = "response")
y_pred_enet_class <- ifelse(y_pred_enet > 0.5, 1, 0)
accuracy_enet <- calc_accuracy(y_test, y_pred_enet_class)

y_pred_opt_ridge <- predict(cv_ridge, newx = X_test, s = cv_ridge$lambda.min, type = "response")
y_pred_opt_ridge_class <- ifelse(y_pred_opt_ridge > 0.5, 1, 0)
accuracy_opt_ridge <- calc_accuracy(y_test, y_pred_opt_ridge_class)

# Store results in a data frame
results_df <- data.frame(
  Model = model_names,
  MSE = c(mse_ridge, mse_ols, mse_lasso, mse_enet, mse_opt_ridge),
  Accuracy = c(accuracy_ridge, accuracy_ols, accuracy_lasso, accuracy_enet, accuracy_opt_ridge)
)

# Print or save the results
print(results_df)

#If you want to save the result as csv
write.csv(results_df, "heart_failure_results_mse.csv", row.names = FALSE)

#If you want to save the lambdas
lambdas_df <- data.frame(
  Lambda_name = lambda_names,
  Lambda_value = unlist(lapply(lambda_names, function(x) compute_lambda(x, sigma_squared, eigenvalues, beta_hat, beta_hat_squared, n, p)))
)

write.csv(lambdas_df, "heart_failure_lambdas_mse.csv", row.names = FALSE)
