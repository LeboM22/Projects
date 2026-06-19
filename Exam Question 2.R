###################################################
# QUESTION 2.1
# KDE FROM FIRST PRINCIPLES
# KS-CV USING KDE CDF
###################################################

library(readxl)
library(ggplot2)

###################################################
# DATA
###################################################

data <- read_excel(
  "C:/Users/lebom/OneDrive/Documents/Masters/MVA 880/Exam 2026/randmixa.xlsx"
)

X <- data$X
n <- length(X)

###################################################
# GAUSSIAN KERNEL
###################################################

gaussian_kernel <- function(u){
  dnorm(u)
}

###################################################
# KDE DENSITY
###################################################

kde_density <- function(x_eval, data, h){
  
  sapply(x_eval, function(x){
    
    mean(
      gaussian_kernel((x-data)/h)
    )/h
    
  })
  
}

###################################################
# KDE CDF
###################################################

kde_cdf <- function(grid, density_vals){
  
  dx <- diff(grid)
  
  cdf <- c(
    0,
    cumsum(
      (head(density_vals,-1) +
         tail(density_vals,-1))/2 * dx
    )
  )
  
  cdf / max(cdf)
  
}

###################################################
# BANDWIDTH GRID
###################################################

h_grid <- seq(0.1,5,length.out=40)

###################################################
# 10-FOLD CV
###################################################

Kfold <- 10

set.seed(123)

folds <- sample(
  rep(1:Kfold,length.out=n)
)

ks_values <- numeric(length(h_grid))

###################################################
# CV LOOP
###################################################

for(i in seq_along(h_grid)){
  
  h <- h_grid[i]
  
  fold_ks <- numeric(Kfold)
  
  for(k in 1:Kfold){
    
    train <- X[folds != k]
    test  <- X[folds == k]
    
    #################################################
    # KDE ON TRAINING DATA
    #################################################
    
    grid <- seq(
      min(train)-3*h,
      max(train)+3*h,
      length.out = 500
    )
    
    dens_train <- kde_density(
      grid,
      train,
      h
    )
    
    cdf_train <- kde_cdf(
      grid,
      dens_train
    )
    
    #################################################
    # KDE-CDF AT TEST POINTS
    #################################################
    
    kde_cdf_test <- approx(
      x = grid,
      y = cdf_train,
      xout = sort(test),
      rule = 2
    )$y
    
    #################################################
    # TEST EMPIRICAL CDF
    #################################################
    
    F_test <- ecdf(test)
    
    emp_test <- F_test(sort(test))
    
    #################################################
    # KS
    #################################################
    
    fold_ks[k] <- max(
      abs(emp_test - kde_cdf_test)
    )
  }
  
  ks_values[i] <- mean(fold_ks)
}

###################################################
# BEST BANDWIDTH
###################################################

best_h <- h_grid[
  which.min(ks_values)
]

cat(
  "Optimal bandwidth =",
  best_h,
  "\n"
)

###################################################
# PLOT
###################################################

df_h <- data.frame(
  h = h_grid,
  KS = ks_values
)

ggplot(df_h,
       aes(h,KS))+
  
  geom_line(
    linewidth=1
  )+
  
  geom_point(
    size=2
  )+
  
  geom_vline(
    xintercept = best_h,
    colour = "red",
    linetype = "dashed"
  )+
  
  annotate(
    "text",
    x = best_h,
    y = max(ks_values),
    label = paste0(
      "h = ",
      round(best_h,2)
    ),
    colour = "red",
    hjust = -0.1
  )+
  
  labs(
    title = "10-Fold CV: KS Statistic vs Bandwidth",
    x = "Bandwidth (h)",
    y = "Average Out-of-Sample KS"
  )+
  
  theme_minimal()







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
h_opt <- best_h  # given from your CV

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
# QUESTION 2.5 (MODE-SEEKING CLUSTERING USING KDE)
###################################################

x <- X
h <- 1  # bandwidth from Q2.1

# -----------------------------
# 1. KDE ESTIMATION (from Q2.1)
# -----------------------------
kde <- density(x,
               bw = h,
               kernel = "gaussian",
               n = 2048)

# -----------------------------
# 2. IDENTIFY KDE MODES
# -----------------------------
modes_index <- which(diff(sign(diff(kde$y))) == -2) + 1
modes <- kde$x[modes_index]

cat("Estimated KDE modes:\n")
print(modes)

# -----------------------------
# 3. VISUALISE KDE + MODES
# -----------------------------
plot(kde,
     main = "Kernel Density Estimate with Modes",
     xlab = "X",
     ylab = "Density",
     lwd = 2)

abline(v = modes,
       col = "red",
       lty = 2,
       lwd = 2)

points(modes,
       kde$y[modes_index],
       pch = 19,
       col = "red")

# -----------------------------
# 4. MEAN-SHIFT (MODE-SEEKING STEP)
# -----------------------------
mean_shift <- function(x0, data, h,
                       tol = 1e-6,
                       max_iter = 100){
  
  x_current <- x0
  
  for(iter in 1:max_iter){
    
    weights <- dnorm((x_current - data)/h)
    x_new <- sum(weights * data) / sum(weights)
    
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

# -----------------------------
# 5. ASSIGN EACH POINT TO NEAREST KDE MODE
# -----------------------------
assign_mode <- function(val, modes){
  which.min(abs(val - modes))
}

cluster_id <- sapply(mode_destination,
                     assign_mode,
                     modes = modes)

cluster_id <- as.factor(cluster_id)

# -----------------------------
# 6. CLUSTER SUMMARY (REQUIRED FOR REPRESENTATION)
# -----------------------------
cluster_summary <- aggregate(
  X ~ cluster_id,
  data = data.frame(X = x, cluster_id = cluster_id),
  FUN = function(z) c(
    Size = length(z),
    Mean = mean(z),
    SD = sd(z)
  )
)

cluster_summary <- do.call(data.frame, cluster_summary)
colnames(cluster_summary) <- c("Cluster", "Size", "Mean", "SD")

cat("\nCluster Summary:\n")
print(cluster_summary)

# -----------------------------
# 7. CLUSTERING VISUALISATION (FINAL REPRESENTATION)
# -----------------------------
plot(x,
     rep(0, length(x)),
     col = cluster_id,
     pch = 19,
     xlab = "X",
     ylab = "",
     main = "Mode-Seeking Clustering based on KDE")

# KDE mode locations (final cluster representatives)
abline(v = modes,
       col = "red",
       lwd = 2,
       lty = 2)

points(modes,
       rep(0, length(modes)),
       pch = 8,
       col = "red",
       cex = 1.5)

# Legend
legend("bottomright",
       legend = c("KDE Modes (Clusters)"),
       col = c("black", "red"),
       pch = c(19, 8),
       lty = c(NA, 2),
       bty = "n")

# -----------------------------
# 8. CLUSTER FREQUENCY TABLE
# -----------------------------
table(cluster_id)








########  FOR COMAPRISON IN Q2.6
# 1. HARD CLUSTERING FUNCTION (from fitted GMM)
###################################################
############################################################
# GMM DENSITY CURVE (FIRST PRINCIPLES)
############################################################

# function to compute Gaussian mixture density on grid
gmm_density <- function(x_grid, fit) {
  
  K <- length(fit$pi)
  n_grid <- length(x_grid)
  
  dens <- numeric(n_grid)
  
  for (k in 1:K) {
    
    dens <- dens + fit$pi[k] * dnorm(
      x_grid,
      mean = fit$mu[k],
      sd = fit$sigma[k]
    )
  }
  
  return(dens)
}

############################################################
# COMPUTE GMM CURVE FOR X
############################################################

gmm_X_vals <- gmm_density(x_grid, fit_X)

############################################################
# PLOT
############################################################
# Compute a safe upper y-limit
y_max <- max(
  hist(X, plot = FALSE, probability = TRUE, breaks = 30)$density,
  final_kde,
  gmm_X_vals
) * 1.2   # add 20% headroom


par(mar = c(5, 5, 3, 2),
    pin = c(3.1, 2.0),
    cex = 0.85)

y_max <- max(
  hist(X, plot = FALSE, probability = TRUE, breaks = 30)$density,
  final_kde,
  gmm_X_vals
) * 1.2

hist(X,
     probability = TRUE,
     breaks = 30,
     col = "grey85",
     border = "white",
     main = "X: Histogram + KDE + GMM",
     xlab = "X",
     ylab = "Density",
     ylim = c(0, y_max))

# KDE
lines(x_grid, final_kde,
      col = "blue",
      lwd = 3)

# GMM (EM)
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
       cex = 0.85)
