############################################################
# LOAD PACKAGES
############################################################

library(readr)
library(ggplot2)
library(dplyr)

############################################################
# LOAD DATA
############################################################

data <- read_csv("C:\\Users\\lebom\\OneDrive\\Documents\\Masters\\MVA 880\\Exam 2026\\gmrdat2026.csv")

X <- data$x
Y <- data$y
n <- length(Y)

############################################################
# 1. POLYNOMIAL DESIGN MATRIX (FIRST PRINCIPLES)
############################################################

poly_design <- function(x, degree) {
  Xmat <- matrix(1, nrow = length(x), ncol = degree + 1)
  for (d in 1:degree) {
    Xmat[, d + 1] <- x^d
  }
  return(Xmat)
}

gmr_em <- function(X, Y, K, degree,
                   max_iter = 100,
                   tol = 1e-6) {
  
  n <- length(Y)
  Xmat <- poly_design(X, degree)
  p <- ncol(Xmat)
  
  ########################################################
  # INIT (k-means)
  ########################################################
  
  set.seed(123)
  clusters <- kmeans(Y, centers = K)$cluster
  
  pi_k <- as.numeric(table(clusters) / n)
  
  beta_k <- list()
  sigma_k <- numeric(K)
  
  for (k in 1:K) {
    idx <- which(clusters == k)
    
    # FIX: avoid empty clusters
    if (length(idx) < p) {
      idx <- sample(1:n, p + 1)
    }
    
    fit <- lm(Y[idx] ~ Xmat[idx, ] - 1)
    
    beta_k[[k]] <- as.numeric(coef(fit))
    
    sigma_k[k] <- sd(residuals(fit))
    if (is.na(sigma_k[k]) || sigma_k[k] <= 1e-6) sigma_k[k] <- 1
  }
  
  loglik_old <- -Inf
  
  ########################################################
  # EM LOOP
  ########################################################
  
  for (iter in 1:max_iter) {
    
    gamma <- matrix(0, n, K)
    
    for (k in 1:K) {
      
      mu_k <- as.vector(Xmat %*% beta_k[[k]])
      
      sigma_k[k] <- max(sigma_k[k], 1e-6)
      
      gamma[, k] <- pi_k[k] * dnorm(Y, mean = mu_k, sd = sigma_k[k])
    }
    
    row_sums <- rowSums(gamma)
    
    # FIX: avoid division by zero
    row_sums[row_sums == 0] <- 1e-12
    
    gamma <- gamma / row_sums
    
    Nk <- colSums(gamma)
    
    pi_k <- Nk / n
    
    for (k in 1:K) {
      
      W <- diag(gamma[, k] + 1e-8)
      
      XtWX <- t(Xmat) %*% W %*% Xmat
      XtWY <- t(Xmat) %*% W %*% Y
      
      # FIX: ridge regularisation for stability
      beta_k[[k]] <- solve(XtWX + diag(1e-6, ncol(XtWX)), XtWY)
      
      residuals <- Y - Xmat %*% beta_k[[k]]
      
      sigma_k[k] <- sqrt(sum(gamma[, k] * residuals^2) / Nk[k])
      sigma_k[k] <- max(sigma_k[k], 1e-6)
    }
    
    ########################################################
    # LOG-LIKELIHOOD (SAFE)
    ########################################################
    
    ll_vec <- numeric(n)
    
    for (i in 1:n) {
      
      dens_i <- numeric(K)
      
      for (k in 1:K) {
        mu_ik <- sum(Xmat[i, ] * beta_k[[k]])
        
        dens_i[k] <- pi_k[k] * dnorm(Y[i],
                                     mean = mu_ik,
                                     sd = sigma_k[k])
      }
      
      ll_vec[i] <- sum(dens_i)
    }
    
    ll_vec[ll_vec <= 1e-12] <- 1e-12
    
    loglik <- sum(log(ll_vec))
    
    if (abs(loglik - loglik_old) < tol) break
    
    loglik_old <- loglik
  }
  
  return(list(pi = pi_k,
              beta = beta_k,
              sigma = sigma_k,
              gamma = gamma,
              loglik = loglik))
}


poly_design <- function(x, degree) {
  Xmat <- matrix(1, nrow = length(x), ncol = degree + 1)
  for (d in 1:degree) {
    Xmat[, d + 1] <- x^d
  }
  return(Xmat)
}
############################################################
# 3. MODEL SELECTION (K + DEGREE)
############################################################

Ks <- 1:5
degrees <- 1:3

results <- data.frame()

for (d in degrees) {
  for (K in Ks) {
    
    cat("Fitting K =", K, "Degree =", d, "\n")
    
    fit <- tryCatch(
      gmr_em(X, Y, K, d),
      error = function(e) NULL
    )
    
    if (is.null(fit)) next
    
    p <- K * (d + 2)
    
    bic <- -2 * fit$loglik + log(length(Y)) * p
    
    results <- rbind(results,
                     data.frame(K = K,
                                degree = d,
                                loglik = fit$loglik,
                                BIC = bic))
  }
}
############################################################
# 4. SELECT BEST MODEL
############################################################

best <- results[which.min(results$BIC), ]
best
cat("Selected polynomial degree:", best$degree, "\n")
############################################################
# 5. FIT FINAL MODEL
############################################################

final_fit <- gmr_em(X, Y, best$K, best$degree)

############################################################
# 6. PREDICTIONS
############################################################

x_grid <- seq(min(X), max(X), length.out = 200)

Xgrid <- poly_design(x_grid, best$degree)

gmr_comp <- matrix(0, nrow = length(x_grid), ncol = best$K)

for (k in 1:best$K) {
  gmr_comp[, k] <- Xgrid %*% final_fit$beta[[k]]
}

gmr_mean <- rowSums(gmr_comp * final_fit$pi)

############################################################
# 7. VISUALISATION
############################################################
library(ggplot2)
library(tidyr)
library(dplyr)

df_data <- data.frame(X = X, Y = Y)

df_pred <- data.frame(
  X = x_grid,
  MixtureMean = gmr_mean
)

# add components
for (k in 1:best$K) {
  df_pred[[paste0("Component_", k)]] <- gmr_comp[, k]
}
df_long <- df_pred %>%
  pivot_longer(
    cols = -X,
    names_to = "Series",
    values_to = "Yhat"
  )
ggplot() +
  
  # data
  geom_point(data = df_data,
             aes(x = X, y = Y),
             color = "grey70",
             alpha = 0.7,
             size = 1.5) +
  
  # fitted curves
  geom_line(data = df_long,
            aes(x = X, y = Yhat, color = Series, linetype = Series),
            linewidth = 1.2) +
  
  scale_color_manual(values = c(
    "MixtureMean" = "black",
    setNames(rainbow(best$K), paste0("Component_", 1:best$K))
  )) +
  
  scale_linetype_manual(values = c(
    "MixtureMean" = "solid",
    setNames(rep("dashed", best$K), paste0("Component_", 1:best$K))
  )) +
  
  labs(
    title = "Polynomial Gaussian Mixture Regression (First Principles)",
    x = "X",
    y = "Y",
    color = "Model Component",
    linetype = "Model Component"
  ) +
  
  theme_minimal() +
  
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )



### Scree plot
ggplot(results, aes(x = K, y = BIC)) +
  geom_line() +
  geom_point() +
  facet_wrap(~ degree) +
  labs(
    title = "BIC Scree Curves by Polynomial Degree",
    x = "Number of Components (K)",
    y = "BIC"
  ) +
  theme_minimal()



### Degree
library(ggplot2)

ggplot(results, aes(x = degree, y = BIC, group = factor(K), color = factor(K))) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  labs(
    title = "BIC vs Polynomial Degree for Different K",
    x = "Polynomial Degree (d)",
    y = "BIC",
    color = "Number of Components (K)"
  ) +
  theme_minimal()



cat("====================================\n")
cat("FINAL SELECTED GMR MODEL\n")
cat("====================================\n\n")

cat("Selected polynomial degree:", best$degree, "\n")
cat("Selected number of components K:", best$K, "\n\n")

for (k in 1:best$K) {
  
  cat("Component", k, ":\n")
  cat("pi =", final_fit$pi[k], "\n")
  cat("sigma =", final_fit$sigma[k], "\n")
  cat("beta coefficients:\n")
  print(final_fit$beta[[k]])
  cat("\n------------------------\n")
}









############################################################
# 3.2 SHANNON ENTROPY UNCERTAINTY ANALYSIS (FINAL CLEAN)
############################################################
library(ggplot2)
library(tidyr)
library(dplyr)

############################################################
# 1. DATA PREPARATION
############################################################

# Main data frame (observations)
df <- data.frame(
  X = X,
  Y = Y,
  uncertain = uncertain
)

df$uncertainty_label <- ifelse(df$uncertain,
                               "High uncertainty",
                               "Low uncertainty")

# Mixture component curves
comp_df <- as.data.frame(gmr_comp)
comp_df$x_grid <- x_grid

comp_long <- pivot_longer(
  comp_df,
  cols = -x_grid,
  names_to = "component",
  values_to = "y"
)

# Mixture mean
mean_df <- data.frame(
  x_grid = x_grid,
  gmr_mean = gmr_mean
)

############################################################
# 2. PLOT
############################################################

ggplot() +
  
  ########################################################
# Observations (uncertainty only)
########################################################
geom_point(
  data = df,
  aes(x = X, y = Y, color = uncertainty_label),
  size = 2
) +
  
  scale_color_manual(
    values = c(
      "Low uncertainty" = "grey70",
      "High uncertainty" = "red"
    )
  ) +
  
  ########################################################
# Mixture mean (black solid line)
########################################################
geom_line(
  data = mean_df,
  aes(x = x_grid, y = gmr_mean),
  color = "black",
  linewidth = 1.3
) +
  
  ########################################################
# Component curves (fixed colors: blue, green, red)
########################################################
geom_line(
  data = comp_long,
  aes(x = x_grid, y = y, group = component, color = component),
  linetype = "dashed",
  linewidth = 1
) +
  
  scale_color_manual(
    values = c(
      "Low uncertainty" = "grey70",
      "High uncertainty" = "red",
      "Component 1" = "blue",
      "Component 2" = "green",
      "Component 3" = "red"
    )
  ) +
  
  ########################################################
# Labels
########################################################
labs(
  title = "GMR Model with Mixture Components and Uncertainty",
  x = "X",
  y = "Y",
  color = "Legend"
) +
  
  theme_minimal() +
  
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold")
  )




library(ggplot2)

gamma <- final_fit$gamma
entropy <- -rowSums(gamma * log(gamma + 1e-12))

cutoff <- quantile(entropy, 0.75)

entropy_df <- data.frame(entropy = entropy)

ggplot(entropy_df, aes(x = entropy)) +
  geom_histogram(bins = 30,
                 fill = "grey80",
                 color = "white") +
  geom_vline(xintercept = cutoff,
             color = "red",
             linetype = "dashed",
             linewidth = 1.2) +
  labs(
    title = "Shannon Entropy of GMR Classification Uncertainty",
    x = "Entropy",
    y = "Frequency"
  ) +
  theme_minimal()

############################################################
# SHANNON ENTROPY PER OBSERVATION (Hi)
############################################################

# posterior probabilities from GMR
gamma <- final_fit$gamma

# ensure numerical stability (important)
gamma <- gamma / rowSums(gamma)

# compute entropy for each observation
entropy <- -rowSums(gamma * log(gamma + 1e-12))

# attach to dataset for inspection
df_entropy <- data.frame(
  X = X,
  Y = Y,
  entropy = entropy
)

# show first few values
head(df_entropy)

# summary of entropy
summary(entropy)

# optional: print max/min
cat("Min entropy:", min(entropy), "\n")
cat("Max entropy:", max(entropy), "\n")


library(ggplot2)

df_entropy$uncertainty_level <- cut(
  df_entropy$entropy,
  breaks = quantile(df_entropy$entropy, probs = c(0, 0.5, 0.75, 1)),
  labels = c("Low", "Medium", "High"),
  include.lowest = TRUE
)

ggplot() +
  
  # data points coloured by entropy level
  geom_point(
    data = df_entropy,
    aes(x = X, y = Y, color = uncertainty_level),
    size = 2,
    alpha = 0.8
  ) +
  
  # mixture mean
  geom_line(
    data = mean_df,
    aes(x = x_grid, y = gmr_mean),
    color = "black",
    linewidth = 1.3
  ) +
  
  # component curves
  geom_line(
    data = comp_long,
    aes(x = x_grid, y = y, group = component),
    linetype = "dashed",
    linewidth = 1,
    color = "grey40"
  ) +
  
  scale_color_manual(
    values = c(
      "Low" = "darkgreen",
      "Medium" = "orange",
      "High" = "red"
    )
  ) +
  
  labs(
    title = "GMR Model with Shannon Entropy-Based Uncertainty",
    x = "X",
    y = "Y",
    color = "Entropy Level"
  ) +
  
  theme_minimal()















############################################################
# 3.3 R-SQUARED FOR GMR (CORRECT FIRST PRINCIPLES)
############################################################

# Design matrix (must match model)
Xmat <- poly_design(X, best$degree)

K <- best$K
n <- length(Y)

# extract parameters
gamma <- final_fit$gamma
beta <- final_fit$beta

# store component means
mu_mat <- matrix(0, nrow = n, ncol = K)

for (k in 1:K) {
  mu_mat[, k] <- as.vector(Xmat %*% beta[[k]])
}

# proper GMR prediction: weighted expectation
y_hat <- rowSums(gamma * mu_mat)

# compute R2
y_bar <- mean(Y)

SST <- sum((Y - y_bar)^2)
SSE <- sum((Y - y_hat)^2)

R2 <- 1 - SSE / SST

cat("Corrected GMR R-squared:", R2, "\n")








############################################################
# 3.4 BOOTSTRAP 95% CI FOR R^2 (FIXED FIRST PRINCIPLES)
############################################################

set.seed(123)

B <- 200
n <- length(Y)

bootstrap_R2 <- numeric(0)

for (b in 1:B) {
  
  ########################################################
  # Step 1: bootstrap sample
  ########################################################
  
  idx <- sample(1:n, size = n, replace = TRUE)
  
  X_b <- X[idx]
  Y_b <- Y[idx]
  
  ########################################################
  # Step 2: fit GMR model with RANDOM initialization
  ########################################################
  
  fit_b <- tryCatch({
    
    # IMPORTANT: refit model
    gmr_em(X_b, Y_b, best$K, best$degree)
    
  }, error = function(e) NULL)
  
  if (is.null(fit_b)) next
  
  ########################################################
  # Step 3: design matrix
  ########################################################
  
  Xmat_b <- poly_design(X_b, best$degree)
  
  gamma_b <- fit_b$gamma
  beta_b <- fit_b$beta
  
  K <- best$K
  n_b <- length(Y_b)
  
  ########################################################
  # Step 4: compute component means
  ########################################################
  
  mu_mat_b <- matrix(0, nrow = n_b, ncol = K)
  
  for (k in 1:K) {
    mu_mat_b[, k] <- as.vector(Xmat_b %*% beta_b[[k]])
  }
  
  ########################################################
  # Step 5: GMR prediction (first principles)
  ########################################################
  
  y_hat_b <- rowSums(gamma_b * mu_mat_b)
  
  ########################################################
  # Step 6: compute R^2
  ########################################################
  
  y_bar_b <- mean(Y_b)
  
  SST_b <- sum((Y_b - y_bar_b)^2)
  SSE_b <- sum((Y_b - y_hat_b)^2)
  
  R2_b <- 1 - SSE_b / SST_b
  
  ########################################################
  # Step 7: store ONLY valid values
  ########################################################
  
  if (is.finite(R2_b)) {
    bootstrap_R2 <- c(bootstrap_R2, R2_b)
  }
}

############################################################
# SAFETY CHECK
############################################################

cat("Valid bootstrap samples:", length(bootstrap_R2), "\n")

summary(bootstrap_R2)
sd(bootstrap_R2)

############################################################
# 95% CONFIDENCE INTERVAL
############################################################

CI_lower <- quantile(bootstrap_R2, 0.025)
CI_upper <- quantile(bootstrap_R2, 0.975)

cat("\nBootstrap 95% CI for R2:\n")
cat("Lower bound:", CI_lower, "\n")
cat("Upper bound:", CI_upper, "\n")


############################################################
# BOOTSTRAP DISTRIBUTION PLOT FOR R^2
############################################################

library(ggplot2)

df_boot <- data.frame(R2 = bootstrap_R2)

ggplot(df_boot, aes(x = R2)) +
  
  geom_histogram(
    bins = 30,
    fill = "grey80",
    color = "white"
  ) +
  
  geom_vline(
    xintercept = quantile(bootstrap_R2, 0.025),
    color = "red",
    linetype = "dashed",
    linewidth = 1.1
  ) +
  
  geom_vline(
    xintercept = quantile(bootstrap_R2, 0.975),
    color = "red",
    linetype = "dashed",
    linewidth = 1.1
  ) +
  
  labs(
    title = "Bootstrap Distribution of R² (GMR Model)",
    x = expression(R^2),
    y = "Frequency"
  ) +
  
  theme_minimal() +
  
  theme(
    plot.title = element_text(face = "bold")
  )