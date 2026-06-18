# Load necessary libraries
library(readr)
library(ggplot2)
library(dplyr)
library(readxl)


# Read the data
data <- read_excel("C:\\Users\\lebom\\OneDrive\\Documents\\Masters\\MVA 880\\Exam 2026\\randmixa.xlsx")


###################################################
################ Question 1.1
###################################################
# Histogram with ggplot2
ggplot(data, aes(x = X)) +
  geom_histogram(binwidth = 5, fill = "skyblue", color = "black") +
  labs(title = "Histogram of Variable X", x = "X", y = "Frequency") +
  theme_minimal()
ggplot(data, aes(x = y)) +
  geom_histogram(binwidth = 5, fill = "lightgreen", color = "black") +
  labs(title = "Histogram of Variable y", x = "y", y = "Frequency") +
  theme_minimal()



###################################################
################ Question 1.2
##################################################
# 0. Libraries
###################################################
library(mvtnorm)

###################################################
# 1. Log-likelihood function
###################################################

gmm_loglik_1d <- function(data, pi, mu, sigma, K) {
  n <- length(data)
  dens <- matrix(0, nrow = n, ncol = K)
  
  for (k in 1:K) {
    dens[, k] <- pi[k] * dnorm(data, mean = mu[k], sd = sigma[k])
  }
  
  sum(log(rowSums(dens)))
}

###################################################
# 2. EM algorithm with K-means initialisation
###################################################

gmm_em_1d <- function(data, K, max_iter = 200, tol = 1e-6) {
  n <- length(data)
  
  # ---- K-means init ----
  set.seed(123)
  km <- kmeans(data, centers = K, nstart = 20)
  
  mu <- km$centers
  
  sigma <- numeric(K)
  pi <- numeric(K)
  
  for (k in 1:K) {
    cluster_k <- data[km$cluster == k]
    
    if (length(cluster_k) < 2) {
      sigma[k] <- sd(data)
    } else {
      sigma[k] <- sd(cluster_k)
    }
    
    if (is.na(sigma[k]) || sigma[k] == 0) {
      sigma[k] <- sd(data)
    }
    
    pi[k] <- length(cluster_k) / n
  }
  
  loglik_old <- -Inf
  
  # ---- EM loop ----
  for (iter in 1:max_iter) {
    
    # E-step
    gamma <- matrix(0, n, K)
    
    for (k in 1:K) {
      gamma[, k] <- pi[k] * dnorm(data, mu[k], sigma[k])
    }
    
    gamma <- gamma / rowSums(gamma)
    
    # M-step
    Nk <- colSums(gamma)
    
    for (k in 1:K) {
      mu[k] <- sum(gamma[, k] * data) / Nk[k]
      sigma[k] <- sqrt(sum(gamma[, k] * (data - mu[k])^2) / Nk[k])
      pi[k] <- Nk[k] / n
    }
    
    # log-likelihood
    loglik <- sum(log(rowSums(sapply(1:K, function(k)
      pi[k] * dnorm(data, mu[k], sigma[k])))))
    
    if (abs(loglik - loglik_old) < tol) break
    loglik_old <- loglik
  }
  
  list(pi = pi, mu = mu, sigma = sigma, loglik = loglik)
}

###################################################
# 3. BIC function
###################################################

compute_bic <- function(loglik, K, n) {
  p <- (K - 1) + K + K
  -2 * loglik + p * log(n)
}

###################################################
# 4. Data
###################################################

X <- data$X
Y <- data$y

K_values <- 1:5

###################################################
# 5. MODEL SELECTION (BIC)
###################################################

bic_X <- numeric(length(K_values))
bic_Y <- numeric(length(K_values))

models_X <- list()
models_Y <- list()

for (i in seq_along(K_values)) {
  
  K <- K_values[i]
  
  # ---- X ----
  models_X[[i]] <- gmm_em_1d(X, K)
  bic_X[i] <- compute_bic(models_X[[i]]$loglik, K, length(X))
  
  # ---- Y ----
  models_Y[[i]] <- gmm_em_1d(Y, K)
  bic_Y[i] <- compute_bic(models_Y[[i]]$loglik, K, length(Y))
}

best_K_X <- K_values[which.min(bic_X)]
best_K_Y <- K_values[which.min(bic_Y)]

###################################################
# 6. CONFIRM BEST K
###################################################

cat("Best K for X:", best_K_X, "\n")
cat("Best K for Y:", best_K_Y, "\n")

###################################################
# 7. FINAL FIT USING K = 3 (CONFIRMED)
###################################################

K_final <- 3

fit_X <- gmm_em_1d(X, K_final)
fit_Y <- gmm_em_1d(Y, K_final)

###################################################
# 8. PRINT FINAL ESTIMATES
###################################################

cat("\n====================================\n")
cat("FINAL GMM (X) - K = 3\n")
cat("====================================\n")

cat("\nMixing proportions (pi):\n")
print(fit_X$pi)

cat("\nMeans (mu):\n")
print(fit_X$mu)

cat("\nStandard deviations (sigma):\n")
print(fit_X$sigma)


cat("\n====================================\n")
cat("FINAL GMM (Y) - K = 3\n")
cat("====================================\n")

cat("\nMixing proportions (pi):\n")
print(fit_Y$pi)

cat("\nMeans (mu):\n")
print(fit_Y$mu)

cat("\nStandard deviations (sigma):\n")
print(fit_Y$sigma)

###################################################
# 9. OPTIONAL: BIC plots
###################################################

plot(K_values, bic_X, type = "b", pch = 19,
     xlab = "K", ylab = "BIC", main = "BIC for X")

plot(K_values, bic_Y, type = "b", pch = 19,
     xlab = "K", ylab = "BIC", main = "BIC for Y")







###################################################
################ Question 1.3
#####################################################
# 1. HARD CLUSTERING FUNCTION (from fitted GMM)
###################################################

hard_cluster_gmm <- function(data, fit, K) {
  
  n <- length(data)
  
  ###################################################
  # Step 1: Compute responsibilities
  ###################################################
  
  responsibilities <- matrix(0, nrow = n, ncol = K)
  
  for (k in 1:K) {
    responsibilities[, k] <- fit$pi[k] * dnorm(
      data,
      mean = fit$mu[k],
      sd = fit$sigma[k]
    )
  }
  
  responsibilities <- responsibilities / rowSums(responsibilities)
  
  ###################################################
  # Step 2: Hard cluster assignment
  ###################################################
  
  hard_cluster <- apply(responsibilities, 1, which.max)
  
  ###################################################
  # Step 3: Re-estimate parameters from hard clusters
  ###################################################
  
  hard_pi <- as.numeric(table(hard_cluster) / n)
  
  hard_mu <- tapply(data, hard_cluster, mean)
  hard_sigma <- tapply(data, hard_cluster, sd)
  
  ###################################################
  # Step 4: Organise results
  ###################################################
  
  hard_estimates <- data.frame(
    Component = 1:K,
    PI_hard = round(hard_pi, 4),
    mu_hard = round(hard_mu, 4),
    sigma_hard = round(hard_sigma, 4)
  )
  
  em_estimates <- data.frame(
    Component = 1:K,
    PI_EM = round(fit$pi, 4),
    mu_EM = round(as.numeric(fit$mu), 4),
    sigma_EM = round(fit$sigma, 4)
  )
  
  comparison <- merge(hard_estimates, em_estimates, by = "Component")
  
  return(list(
    hard_cluster = hard_cluster,
    comparison = comparison
  ))
}

###################################################
# 2. APPLY TO X AND Y (USING K = 3 FITS)
###################################################

K_final <- 3

hard_X <- hard_cluster_gmm(X, fit_X, K_final)
hard_Y <- hard_cluster_gmm(Y, fit_Y, K_final)

###################################################
# 3. OUTPUT RESULTS
###################################################

cat("\n==============================\n")
cat("HARD CLUSTER COMPARISON - X\n")
cat("==============================\n")
print(hard_X$comparison)

cat("\n==============================\n")
cat("HARD CLUSTER COMPARISON - Y\n")
cat("==============================\n")
print(hard_Y$comparison)

###################################################
# 4. Cluster sizes
###################################################

cat("\nCluster counts (X):\n")
print(table(hard_X$hard_cluster))

cat("\nCluster counts (Y):\n")
print(table(hard_Y$hard_cluster))

###################################################
# 5. VISUALISATION
###################################################

library(ggplot2)

###################################################
# (a) X with clusters
###################################################

df_X <- data.frame(
  X = X,
  cluster = as.factor(hard_X$hard_cluster)
)

ggplot(df_X, aes(x = X, fill = cluster)) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  labs(title = "Hard GMM Clustering on X", x = "X", y = "Density")

###################################################
# (b) Y with clusters
###################################################

df_Y <- data.frame(
  Y = Y,
  cluster = as.factor(hard_Y$hard_cluster)
)

ggplot(df_Y, aes(x = Y, fill = cluster)) +
  geom_histogram(bins = 30, alpha = 0.6, position = "identity") +
  labs(title = "Hard GMM Clustering on Y", x = "Y", y = "Density")

###################################################
# (c) X vs Y coloured by X clusters
###################################################

df_XY <- data.frame(
  X = X,
  Y = Y,
  cluster = as.factor(hard_X$hard_cluster)
)

ggplot(df_XY, aes(x = X, y = Y, color = cluster)) +
  geom_point(alpha = 0.7) +
  labs(title = "X vs Y (Clustered using X-GMM)", x = "X", y = "Y")












###################################################
################ Question 1.4
##################################################
###################################################
# 1. SIMULATE FROM FITTED GMM (X and Y)
###################################################

gmm_rsample <- function(n, fit) {
  K <- length(fit$pi)
  
  # sample component labels
  z <- sample(1:K, size = n, replace = TRUE, prob = fit$pi)
  
  # generate observations
  x <- numeric(n)
  
  for (k in 1:K) {
    nk <- sum(z == k)
    x[z == k] <- rnorm(nk, mean = fit$mu[k], sd = fit$sigma[k])
  }
  
  return(x)
}

###################################################
# 2. GENERATE LARGE SAMPLES FROM BOTH MODELS
###################################################

set.seed(123)

n_sim <- 10000

X_sim <- gmm_rsample(n_sim, fit_X)
Y_sim <- gmm_rsample(n_sim, fit_Y)

###################################################
# 3. KS TEST: COMPARE DISTRIBUTIONS
###################################################
###################################################
# KS STATISTIC (FIRST PRINCIPLES)
###################################################

# Ensure sorted samples
X_sorted <- sort(X_sim)
Y_sorted <- sort(Y_sim)

n1 <- length(X_sorted)
n2 <- length(Y_sorted)

# All unique points across both samples
all_points <- sort(unique(c(X_sorted, Y_sorted)))

# Empirical CDFs evaluated at each point
F_X <- ecdf(X_sorted)
F_Y <- ecdf(Y_sorted)

# Compute absolute differences in CDFs
ks_values <- abs(F_X(all_points) - F_Y(all_points))

# KS statistic
ks_statistic <- max(ks_values)

###################################################
# OUTPUT
###################################################

cat("\n====================================\n")
cat("KS TEST (FIRST PRINCIPLES)\n")
cat("====================================\n")

cat("\nKS Statistic (D):", ks_statistic, "\n")

###################################################
# 4. OUTPUT RESULTS
###################################################

cat("\n====================================\n")
cat("KS TEST: GMM(X) vs GMM(Y)\n")
cat("====================================\n")

print(ks_result)

cat("\nKS Statistic:\n")
print(ks_result$statistic)

###################################################
# 5. OPTIONAL VISUAL CHECK
###################################################

plot(density(X_sim),
     col = "blue",
     lwd = 2,
     main = "GMM-Based Simulated Distributions",
     xlab = "Value")

lines(density(Y_sim),
      col = "red",
      lwd = 2)

legend("topright",
       legend = c("X (GMM)", "Y (GMM)"),
       col = c("blue", "red"),
       lwd = 2)










###################################################
################ Question 1.5
###################################################
# 1. SETUP
###################################################

X <- X
Y <- y

n1 <- length(X)
n2 <- length(Y)

pooled <- c(X, Y)

set.seed(123)
B <- 2000

###################################################
# 2. STORAGE
###################################################

boot_ks   <- numeric(B)
boot_mean <- numeric(B)
boot_var  <- numeric(B)

p_track_ks <- numeric(B)

###################################################
# 3. OBSERVED STATISTICS
###################################################

obs_ks   <- ks.test(X, Y)$statistic
obs_mean <- mean(X) - mean(Y)
obs_var  <- var(X) - var(Y)

###################################################
# 4. BOOTSTRAP LOOP (NULL: same distribution)
###################################################

for (b in 1:B) {
  
  sample_b <- sample(pooled, replace = TRUE)
  
  Xb <- sample_b[1:n1]
  Yb <- sample_b[(n1 + 1):(n1 + n2)]
  
  # ---------------- KS statistic ----------------
  boot_ks[b] <- ks.test(Xb, Yb)$statistic
  
  # running KS p-value (convergence tracking)
  p_track_ks[b] <- mean(boot_ks[1:b] >= obs_ks)
  
  # ---------------- mean difference ----------------
  boot_mean[b] <- mean(Xb) - mean(Yb)
  
  # ---------------- variance difference ----------------
  boot_var[b] <- var(Xb) - var(Yb)
}

###################################################
# 5. FINAL BOOTSTRAP RESULTS
###################################################

# KS p-value
p_ks <- mean(boot_ks >= obs_ks)

# Mean test
p_mean <- mean(abs(boot_mean) >= abs(obs_mean))
ci_mean <- quantile(boot_mean, c(0.025, 0.975))

# Variance test
p_var <- mean(abs(boot_var) >= abs(obs_var))
ci_var <- quantile(boot_var, c(0.025, 0.975))

# KS CI
ci_ks <- quantile(boot_ks, c(0.025, 0.975))

###################################################
# 6. OUTPUT RESULTS
###################################################

cat("\n====================================\n")
cat("BOOTSTRAP TEST: X vs Y\n")
cat("====================================\n")

cat("\nKS Statistic:\n")
cat("Observed:", obs_ks, "\n")
cat("p-value:", round(p_ks, 4), "\n")
cat("95% CI:", ci_ks, "\n")

cat("\nMean Difference:\n")
cat("Observed:", obs_mean, "\n")
cat("p-value:", round(p_mean, 4), "\n")
cat("95% CI:", ci_mean, "\n")

cat("\nVariance Difference:\n")
cat("Observed:", obs_var, "\n")
cat("p-value:", round(p_var, 4), "\n")
cat("95% CI:", ci_var, "\n")

###################################################
# 7. KS CONVERGENCE PLOT
###################################################

plot(p_track_ks, type = "l",
     col = "blue",
     lwd = 1,
     xlab = "Bootstrap Iterations",
     ylab = "Estimated p-value",
     main = "Bootstrap Convergence of KS Test")

abline(h = p_ks, col = "red", lwd = 2)

legend("bottomright",
       legend = c("Running p-value", "Final p-value"),
       col = c("blue", "red"),
       lty = 1,
       lwd = c(1, 2))

###################################################
# 8. DISTRIBUTION VISUALISATION
###################################################

plot(density(X),
     col = "blue",
     lwd = 2,
     main = "Density Comparison: X vs Y",
     xlab = "Value")

lines(density(Y),
      col = "red",
      lwd = 2)

legend("topright",
       legend = c("X", "Y"),
       col = c("blue", "red"),
       lwd = 2)

###################################################
# 9. KS BOOTSTRAP DISTRIBUTION
###################################################

plot(density(boot_ks),
     main = "Bootstrap Distribution of KS Statistic",
     xlab = "KS statistic")

abline(v = obs_ks, col = "red", lwd = 2)

###################################################
# 7. VISUALISATION
###################################################

par(mfrow = c(1,3))

hist(boot_ks,
     main = "Bootstrap KS Statistic",
     col = "lightblue",
     xlab = "KS")

abline(v = obs_ks, col = "red", lwd = 2)

hist(boot_mean_diff,
     main = "Bootstrap Mean Difference",
     col = "lightgreen",
     xlab = "Mean Diff")

abline(v = obs_mean_diff, col = "red", lwd = 2)

hist(boot_var_diff,
     main = "Bootstrap Variance Difference",
     col = "lightpink",
     xlab = "Variance Diff")

abline(v = obs_var_diff, col = "red", lwd = 2)








###################################################
################ Question 1.6
###################################################
###################################################
# 1. GMM density function (for any fitted model)
###################################################

gmm_rsample <- function(n, model) {
  K <- length(model$pi)
  
  # Step 1: sample component labels
  z <- sample(1:K, size = n, replace = TRUE, prob = model$pi)
  
  # Step 2: sample from corresponding Gaussian
  x <- numeric(n)
  
  for (k in 1:K) {
    nk <- sum(z == k)
    x[z == k] <- rnorm(nk, mean = model$mu[k], sd = model$sigma[k])
  }
  
  return(x)
}

###################################################
# 2. Monte Carlo simulation
###################################################

set.seed(123)

N <- 100000

X_sim <- gmm_rsample(N, fit_X)
Y_sim <- gmm_rsample(N, fit_Y)

###################################################
# 3. Probability estimates
###################################################

# (a) P(X > 100)
p_X_gt_100 <- mean(X_sim > 100)

# (b) P(Y < 120)
p_Y_lt_120 <- mean(Y_sim < 120)

###################################################
# 4. OUTPUT RESULTS
###################################################

cat("\n====================================\n")
cat("MONTE CARLO PROBABILITIES (GMM)\n")
cat("====================================\n")

cat("\nP(X > 100):", round(p_X_gt_100, 4), "\n")
cat("P(Y < 120):", round(p_Y_lt_120, 4), "\n")

###################################################
# 5. OPTIONAL VISUALISATION
###################################################

par(mfrow = c(1, 2))

hist(X_sim, breaks = 50, col = "lightblue",
     main = "Simulated X (GMM)", xlab = "X")

abline(v = 100, col = "red", lwd = 2)

hist(Y_sim, breaks = 50, col = "lightgreen",
     main = "Simulated Y (GMM)", xlab = "Y")

abline(v = 120, col = "red", lwd = 2)
