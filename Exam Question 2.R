# Load necessary libraries
library(readr)
library(ggplot2)
library(dplyr)
library(readxl)


# Read the data
data <- read_excel("C:\\Users\\lebom\\OneDrive\\Documents\\Masters\\MVA 880\\Exam 2026\\randmixa.xlsx")


###################################################
# QUESTION 2.1 KERNEL DENSITY ESTIMATION (FIRST PRINCIPLES)
###################################################

# Data
X <- data$X
n <- length(X)

###################################################
# 1. Gaussian Kernel Function (FIRST PRINCIPLES)
###################################################

gaussian_kernel <- function(u) {
  (1 / sqrt(2 * pi)) * exp(-0.5 * u^2)
}

###################################################
# 2. KDE function (first principles estimator)
###################################################

kde_estimate <- function(x_grid, data, h) {
  n <- length(data)
  
  sapply(x_grid, function(x) {
    mean(gaussian_kernel((x - data) / h)) / h
  })
}

###################################################
# 3. KDE CDF (numerical integration approximation)
###################################################

kde_cdf <- function(x_grid, kde_values, x) {
  # interpolate KDE then integrate numerically
  dx <- diff(x_grid)
  cdf_vals <- cumsum(c(0, head(kde_values, -1) + diff(c(0, kde_values))/2 * dx))
  
  # normalize
  cdf_vals <- cdf_vals / max(cdf_vals)
  
  approx(x_grid, cdf_vals, xout = x)$y
}

###################################################
# 4. Empirical CDF
###################################################

F_emp <- ecdf(X)

###################################################
# 5. KS statistic for given bandwidth h
###################################################

ks_for_h <- function(h, x_grid) {
  
  kde_vals <- kde_estimate(x_grid, X, h)
  
  F_kde <- function(x) kde_cdf(x_grid, kde_vals, x)
  
  # evaluate KS on grid
  ks_vals <- abs(F_emp(x_grid) - F_kde(x_grid))
  
  max(ks_vals)
}

###################################################
# 6. Cross-validation grid for bandwidth h
###################################################

h_grid <- seq(0.1, 5, length.out = 40)

ks_values <- numeric(length(h_grid))

x_grid <- seq(min(X) - 3, max(X) + 3, length.out = 300)

###################################################
# 7. K-fold cross-validation (first principles)
###################################################

K <- 10
set.seed(123)

folds <- sample(rep(1:K, length.out = n))

for (i in seq_along(h_grid)) {
  
  h <- h_grid[i]
  fold_ks <- numeric(K)
  
  for (k in 1:K) {
    
    train <- X[folds != k]
    test  <- X[folds == k]
    
    # KDE on training set
    kde_train <- function(x) {
      mean(gaussian_kernel((x - train) / h)) / h
    }
    
    # CDF approximation on test grid
    test_grid <- sort(test)
    
    F_train <- function(x) {
      sapply(x, function(z) {
        mean(sapply(train, function(t) gaussian_kernel((z - t) / h))) / h
      })
    }
    
    ks_vals <- abs(ecdf(test)(test_grid) - F_train(test_grid))
    
    fold_ks[k] <- max(ks_vals)
  }
  
  ks_values[i] <- mean(fold_ks)
}

###################################################
# 8. Optimal bandwidth
###################################################

best_h <- h_grid[which.min(ks_values)]

cat("Optimal bandwidth (h):", best_h, "\n")









###################################################
###################################################
# QUESTION 2.2
###################################################

final_kde <- kde_estimate(x_grid, X, best_h)

plot(x_grid, final_kde, type = "l", lwd = 2,
     main = "Final KDE (First Principles)",
     xlab = "X", ylab = "Density")


##### PLots
# Empirical CDF
F_emp <- ecdf(X)

# Final KDE
final_kde_vals <- kde_estimate(x_grid, X, best_h)

# KDE CDF approximation
kde_cdf_vals <- cumsum(final_kde_vals)
kde_cdf_vals <- kde_cdf_vals / max(kde_cdf_vals)

plot(x_grid, F_emp(x_grid),
     type = "l", col = "black", lwd = 2,
     main = "Empirical CDF vs KDE CDF",
     xlab = "X", ylab = "CDF")

lines(x_grid, kde_cdf_vals, col = "blue", lwd = 2)

legend("bottomright",
       legend = c("Empirical CDF", "KDE CDF"),
       col = c("black", "blue"),
       lwd = 2,
       cex = 0.8)



h_test <- c(0.5, 1, 2, 5)

plot(x_grid, kde_estimate(x_grid, X, h_test[1]),
     type = "l", col = "red", lwd = 2,
     ylim = c(0, max(kde_estimate(x_grid, X, h_test[1]))),
     main = "Effect of Bandwidth on KDE",
     xlab = "X", ylab = "Density")

for (h in h_test[-1]) {
  lines(x_grid, kde_estimate(x_grid, X, h), lwd = 2)
}

legend("topright",
       legend = paste("h =", h_test),
       col = c("red", "black", "black", "black"),
       lwd = 2,
       cex = 0.8)




plot(h_grid, ks_values, type = "b", pch = 19,
     col = "darkblue",
     xlab = "Bandwidth (h)",
     ylab = "Out-of-sample KS",
     main = "Cross-Validation:\n(KS vs Bandwidth)")

abline(v = best_h, col = "red", lwd = 2)

points(best_h, min(ks_values), col = "red", pch = 19, cex = 1.5)




# final grid
x_grid <- seq(min(X) - 3, max(X) + 3, length.out = 500)

# final KDE using optimal bandwidth
final_kde <- kde_estimate(x_grid, X, best_h)

# plot histogram + KDE overlay
hist(X,
     probability = TRUE,
     breaks = 30,
     col = "grey85",
     border = "white",
     main = "Kernel Density Estimate of X\n(First Principles)",
     xlab = "X")

lines(x_grid, final_kde,
      col = "blue",
      lwd = 3)

rug(X, col = "darkgrey")




###################################################
###################################################
# QUESTION 2.3
###################################################
# ----------------------------
# 1. Load data (replace x with your variable)
# ----------------------------
x <- data$X  # <-- replace this

# ----------------------------
# 2. Fit KDE with optimal bandwidth
# ----------------------------
h_opt <- 5  # given from your CV

kde <- density(x, bw = h_opt, kernel = "gaussian")

# ----------------------------
# 3. (A) Check non-negativity
# ----------------------------
min_density <- min(kde$y)
print(min_density)  # should be >= 0 (numerically)

# ----------------------------
# 4. (B) Check integral ≈ 1
# ----------------------------
integral_check <- sum(diff(kde$x) * (head(kde$y, -1) + tail(kde$y, -1)) / 2)
print(integral_check)  # should be close to 1

# ----------------------------
# 5. Plot KDE (density function)
# ----------------------------
plot(kde,
     main = "Kernel Density Estimate \n(Valid Density Function)",
     xlab = "x",
     ylab = "Density",
     lwd = 2)

# ----------------------------
# 6. Optional: show area under curve visually
# ----------------------------
polygon(c(kde$x, rev(kde$x)),
        c(kde$y, rep(0, length(kde$y))),
        col = rgb(0, 0, 1, 0.2),
        border = NA)
lines(kde, lwd = 2)







###################################################
###################################################
# QUESTION 2.5
###################################################
# Data
x <- X

# Optimal bandwidth from Question 2.1
h <- 5

# KDE
kde <- density(x,
               bw = h,
               kernel = "gaussian",
               n = 2048)

# Local maxima

modes_index <- which(
  diff(sign(diff(kde$y))) == -2
) + 1

modes <- kde$x[modes_index]

cat("Modes:\n")
print(modes)

plot(kde,
     main="Kernel Density Estimate \n(with Modes)",
     xlab="X",
     ylab="Density",
     lwd=2)

abline(v=modes,
       col="red",
       lty=2,
       lwd=2)

points(modes,
       kde$y[modes_index],
       pch=19,
       col="red")

mean_shift <- function(x0, data, h,
                       tol = 1e-6,
                       max_iter = 100){
  
  x_current <- x0
  
  for(iter in 1:max_iter){
    
    weights <- dnorm((x_current - data)/h)
    
    x_new <- sum(weights * data) /
      sum(weights)
    
    if(abs(x_new - x_current) < tol)
      break
    
    x_current <- x_new
  }
  
  return(x_current)
}

mode_destination <- sapply(x,
                           mean_shift,
                           data = x,
                           h = h)

head(mode_destination)

cluster_modes <- round(mode_destination, 3)

cluster_id <- as.numeric(
  factor(cluster_modes)
)

table(cluster_id)

results <- data.frame(
  Observation = 1:length(x),
  X = x,
  Cluster = cluster_id,
  Mode = cluster_modes
)

head(results)

plot(x,
     rep(0,length(x)),
     col=cluster_id,
     pch=19,
     xlab="X",
     ylab="",
     main="Mode-Seeking Clustering")

abline(v=unique(cluster_modes),
       col="black",
       lwd=2,
       lty=2)

aggregate(X ~ Cluster,
          data=data.frame(
            X=x,
            Cluster=cluster_id
          ),
          FUN=function(z)
            c(Size=length(z),
              Mean=mean(z),
              SD=sd(z)))





########  FOR COMAPRISON IN Q2.6
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
par(mar = c(4, 4, 2, 1))  # minimal margins (THIS is the key)

hist(X,
     probability = TRUE,
     breaks = 30,
     col = "grey85",
     border = "white",
     main = "X: Histogram + KDE + GMM",
     xlab = "X")

lines(x_grid, final_kde,
      col = "blue",
      lwd = 3)

lines(x_grid, gmm_X_vals,
      col = "red",
      lwd = 3,
      lty = 2)

legend("topright",
       legend = c("Histogram", "KDE", "GMM (EM)"),
       col = c("grey85", "blue", "red"),
       lwd = c(10, 3, 3),
       lty = c(1, 1, 2),
       bty = "n",
       cex = 1)