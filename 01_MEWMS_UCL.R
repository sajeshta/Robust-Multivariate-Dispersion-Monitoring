# ============================================================
# 01_MEWMS_UCL.R
# Estimation of UCLs for MEWMS Statistics
# ============================================================

library(MASS)
library(robustbase)
library(rrcov)
library(pcaPP)

# ------------------------------------------------------------
# Function: RobBivDispindUCL
#
# m : number of Monte Carlo trials
# n : number of observations in each trial
# p : dimension
# w : smoothing parameter
#
# Methods:
# m1 = Comedian
# m2 = Sn
# m3 = Classical
# m4 = MCD
# m5 = MVE
# m6 = GK(Q)
# m7 = GK(S)
#
# Returns empirical UCLs for:
# D1 = sum |diag(S_t) - 1|
# D2 = sum (diag(S_t) - 1)^2
# D3 = trace(S_t)
# ------------------------------------------------------------

RobBivDispindUCL <- function(m, n, p, w, method) {
  
  single_trial <- function() {
    
    # --------------------------------------------------------
    # Generate an in-control sample
    # --------------------------------------------------------
    
    g <- mvrnorm(
      n = n,
      mu = rep(0, p),
      Sigma = diag(p)
    )
    
    # --------------------------------------------------------
    # Robust standardization
    # --------------------------------------------------------
    
    switch(
      method,
      
      # ------------------------------------------------------
      # m1 = Comedian
      # ------------------------------------------------------
      
      m1 = {
        
        K1 <- Comed(g)
        
        eig <- eigen(K1)
        
        Si <- eig$vectors %*%
          diag(1 / sqrt(eig$values)) %*%
          t(eig$vectors)
        
        x <- sweep(
          g,
          2,
          colMedians(g)
        ) %*% Si
      },
      
      # ------------------------------------------------------
      # m2 = Sn
      # ------------------------------------------------------
      
      m2 = {
        
        K1 <- covSn(g)
        
        eig <- eigen(K1)
        
        Si <- eig$vectors %*%
          diag(1 / sqrt(eig$values)) %*%
          t(eig$vectors)
        
        x <- sweep(
          g,
          2,
          colMedians(g)
        ) %*% Si
      },
      
      # ------------------------------------------------------
      # m3 = Classical
      # ------------------------------------------------------
      
      m3 = {
        
        K1 <- cov(g)
        
        eig <- eigen(K1)
        
        Si <- eig$vectors %*%
          diag(1 / sqrt(eig$values)) %*%
          t(eig$vectors)
        
        x <- sweep(
          g,
          2,
          colMeans(g)
        ) %*% Si
      },
      
      # ------------------------------------------------------
      # m4 = MCD
      # ------------------------------------------------------
      
      m4 = {
        
        K1 <- covMcd(g)
        
        eig <- eigen(K1$cov)
        
        Si <- eig$vectors %*%
          diag(1 / sqrt(eig$values)) %*%
          t(eig$vectors)
        
        x <- sweep(
          g,
          2,
          K1$center
        ) %*% Si
      },
      
      # ------------------------------------------------------
      # m5 = MVE
      # ------------------------------------------------------
      
      m5 = {
        
        K1 <- CovMve(g)
        
        eig <- eigen(K1$cov)
        
        Si <- eig$vectors %*%
          diag(1 / sqrt(eig$values)) %*%
          t(eig$vectors)
        
        x <- sweep(
          g,
          2,
          K1$center
        ) %*% Si
      },
      
      # ------------------------------------------------------
      # m6 = GK(Q)
      # ------------------------------------------------------
      
      m6 = {
        
        K1 <- gkcov(g, "Q")
        
        eig <- eigen(K1)
        
        Si <- eig$vectors %*%
          diag(1 / sqrt(eig$values)) %*%
          t(eig$vectors)
        
        x <- sweep(
          g,
          2,
          colMedians(g)
        ) %*% Si
      },
      
      # ------------------------------------------------------
      # m7 = GK(S)
      # ------------------------------------------------------
      
      m7 = {
        
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
      },
      
      stop(
        "Invalid method. Use m1, m2, m3, m4, m5, m6, or m7."
      )
    )
    
    # --------------------------------------------------------
    # MEWMS statistics
    # --------------------------------------------------------
    
    S <- array(
      0,
      dim = c(p, p, n)
    )
    
    D1 <- numeric(n)
    D2 <- numeric(n)
    D3 <- numeric(n)
    
    for (i in 1:n) {
      
      xi_xt <- tcrossprod(x[i, ])
      
      if (i == 1) {
        
        S[, , i] <-
          w * xi_xt +
          (1 - w) * xi_xt
        
      } else {
        
        S[, , i] <-
          w * xi_xt +
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
  # Empirical UCLs
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


# ============================================================
# Supporting robust covariance functions
# ============================================================


# ------------------------------------------------------------
# Sn covariance estimator
# ------------------------------------------------------------

covSn <- function(w) {
  
  x <- w[, 1]
  y <- w[, 2]
  
  n <- length(x)
  
  v <- rep(0, n)
  
  for (i in 1:n) {
    
    v[i] <-
      ((x[i] - x) * (y[i] - y))[-c(i)][
        order(
          ((x[i] - x) * (y[i] - y))[-c(i)]
        )
      ][n %/% 2]
  }
  
  sn_cov <-
    1.422295 *
    v[order(v)][(n + 1) %/% 2]
  
  Snc <- matrix(
    c(
      Sn(x)^2,
      sn_cov,
      sn_cov,
      Sn(y)^2
    ),
    2, 2
  )
  
  return(Snc)
}


# ------------------------------------------------------------
# Comedian covariance estimator
# ------------------------------------------------------------

Comed <- function(w) {
  
  x <- w[, 1]
  y <- w[, 2]
  
  comedian <-
    (1.4826^2) *
    median(
      (x - median(x)) *
        (y - median(y))
    )
  
  Com <- matrix(
    c(
      mad(x)^2,
      comedian,
      comedian,
      mad(y)^2
    ),
    2, 2
  )
  
  return(Com)
}


# ------------------------------------------------------------
# Gnanadesikan-Kettenring covariance estimator
#
# method = "Q" : GK(Q)
# method = "S" : GK(S)
# method = "M" : GK(M)
# method = "T" : GK(T)
# ------------------------------------------------------------

gkcov <- function(w, method) {
  
  x <- w[, 1] + w[, 2]
  y <- w[, 1] - w[, 2]
  
  switch(
    method,
    
    # --------------------------------------------------------
    # GK(Q)
    # --------------------------------------------------------
    
    Q = {
      
      gk <-
        0.25 *
        (qn(x)^2 - qn(y)^2)
      
      gkc <- matrix(
        c(
          qn(w[, 1])^2,
          gk,
          gk,
          qn(w[, 2])^2
        ),
        2, 2
      )
    },
    
    # --------------------------------------------------------
    # GK(S)
    # --------------------------------------------------------
    
    S = {
      
      gk <-
        0.25 *
        (Sn(x)^2 - Sn(y)^2)
      
      gkc <- matrix(
        c(
          Sn(w[, 1])^2,
          gk,
          gk,
          Sn(w[, 2])^2
        ),
        2, 2
      )
    },
    
    # --------------------------------------------------------
    # GK(M)
    # --------------------------------------------------------
    
    M = {
      
      gk <-
        0.25 *
        (mad(x)^2 - mad(y)^2)
      
      gkc <- matrix(
        c(
          mad(w[, 1])^2,
          gk,
          gk,
          mad(w[, 2])^2
        ),
        2, 2
      )
    },
    
    # --------------------------------------------------------
    # GK(T)
    # --------------------------------------------------------
    
    T = {
      
      gk <-
        0.25 *
        (
          scaleTau2(x)^2 -
            scaleTau2(y)^2
        )
      
      gkc <- matrix(
        c(
          scaleTau2(w[, 1])^2,
          gk,
          gk,
          scaleTau2(w[, 2])^2
        ),
        2, 2
      )
    },
    
    stop(
      "Invalid method. Use Q, S, M, or T."
    )
  )
  
  return(gkc)
}