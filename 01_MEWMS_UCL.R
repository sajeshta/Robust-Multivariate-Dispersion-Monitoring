# ============================================================
# 01_MEWMS_UCL.R
# Estimation of UCLs for MEWMS Statistics
# ============================================================

library(MASS)
library(robustbase)
library(rrcov)

# ------------------------------------------------------------
# Function: RobBivDispindUCL
#
# m : number of Monte Carlo trials
# n : number of observations in each trial
# p : dimension
# w : smoothing parameter
#
# Returns empirical UCLs for:
# D1 = sum |diag(S_t) - 1|
# D2 = sum (diag(S_t) - 1)^2
# D3 = trace(S_t)
# ------------------------------------------------------------

RobBivDispindUCL <- function(m, n, p, w) {
  
  single_trial <- function() {
    
    # Generate an in-control sample
    g <- mvrnorm(
      n = n,
      mu = rep(0, p),
      Sigma = diag(p)
    )
    
    # --------------------------------------------------------
    # Classical estimator
    # --------------------------------------------------------
    eig <- eigen(cov(g))
    
    Si <- eig$vectors %*%
      diag(1 / sqrt(eig$values)) %*%
      t(eig$vectors)
    
    x <- sweep(
      g,
      2,
      colMeans(g)
    ) %*% Si
    
    # --------------------------------------------------------
    # GK(S) estimator
    # --------------------------------------------------------
    K1 <- gkcov(g, "S")
    
    eig <- eigen(K1)
    
    Si <- eig$vectors %*%
      diag(1 / sqrt(eig$values)) %*%
      t(eig$vectors)
    
    x <- sweep(
      g,
      2,
      colMedians(g)
    ) %*% Si
    
    # --------------------------------------------------------
    # MEWMS statistics
    # --------------------------------------------------------
    
    S <- array(0, dim = c(p, p, n))
    
    D1 <- numeric(n)
    D2 <- numeric(n)
    D3 <- numeric(n)
    
    for (i in 1:n) {
      
      xi_xt <- tcrossprod(x[i, ])
      
      if (i == 1) {
        S[, , i] <- w * xi_xt +
          (1 - w) * xi_xt
      } else {
        S[, , i] <- w * xi_xt +
          (1 - w) * S[, , i - 1]
      }
      
      d <- diag(S[, , i]) - 1
      
      D1[i] <- sum(abs(d))
      D2[i] <- sum(d^2)
      D3[i] <- sum(diag(S[, , i]))
    }
    
    list(
      D1 = D1,
      D2 = D2,
      D3 = D3
    )
  }
  
  # ----------------------------------------------------------
  # Monte Carlo simulation
  # ----------------------------------------------------------
  
  results <- lapply(
    1:m,
    function(i) single_trial()
  )
  
  all_D1 <- unlist(
    lapply(results, `[[`, "D1")
  )
  
  all_D2 <- unlist(
    lapply(results, `[[`, "D2")
  )
  
  all_D3 <- unlist(
    lapply(results, `[[`, "D3")
  )
  
  # ----------------------------------------------------------
  # UCLs
  # ----------------------------------------------------------
  
  QD1 <- quantile(
    all_D1,
    probs = c(0.95, 0.99, 0.9973)
  )
  
  QD2 <- quantile(
    all_D2,
    probs = c(0.95, 0.99, 0.9973)
  )
  
  QD3 <- quantile(
    all_D3,
    probs = c(0.95, 0.99, 0.9973)
  )
  
  QD <- list(
    QD1 = QD1,
    QD2 = QD2,
    QD3 = QD3
  )
  
  return(QD)
}