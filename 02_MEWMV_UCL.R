# ============================================================
# 02_MEWMV_UCL.R
# Estimation of UCLs for MEWMV Statistics
# ============================================================

library(MASS)
library(robustbase)
library(rrcov)
library(pcaPP)

# ------------------------------------------------------------
# Function: RobBivDispind_Meanshift_UCL
#
# m      : number of Monte Carlo trials
# n      : number of observations in each trial
# p      : dimension
# w      : smoothing parameter for dispersion
# lambda : smoothing parameter for mean shift
# method : covariance estimator
#
# Methods:
# m1 = Comedian
# m2 = Sn
# m3 = Classical
# m4 = MCD
# m5 = MVE
# m6 = GK(Q)
# m7 = GK(S)
# ------------------------------------------------------------

RobBivDispind_Meanshift_UCL <- function(
    m, n, p, w, lambda, method) {
  
  ct <- (2 - lambda) /
    (2 * (1 - lambda)^2)
  
  # ----------------------------------------------------------
  # One Monte Carlo trial
  # ----------------------------------------------------------
  
  single_trial <- function() {
    
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
      # Comedian
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
      # Sn
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
      # Classical
      # ------------------------------------------------------
      m3 = {
        eig <- eigen(cov(g))
        
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
      # MCD
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
      # MVE
      # ------------------------------------------------------
      m5 = {
        K2 <- CovMve(g)
        
        eig <- eigen(K2$cov)
        
        Si <- eig$vectors %*%
          diag(1 / sqrt(eig$values)) %*%
          t(eig$vectors)
        
        x <- sweep(
          g,
          2,
          K2$center
        ) %*% Si
      },
      
      # ------------------------------------------------------
      # GK(Q)
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
      # GK(S)
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
      }
    )
    
    # --------------------------------------------------------
    # MEWMV statistics
    # --------------------------------------------------------
    
    S <- array(
      0,
      dim = c(p, p, n)
    )
    
    S1 <- array(
      0,
      dim = c(p, p, n)
    )
    
    D1 <- numeric(n)
    D2 <- numeric(n)
    D3 <- numeric(n)
    
    y <- matrix(
      0,
      n,
      p
    )
    
    # --------------------------------------------------------
    # First observation
    # --------------------------------------------------------
    
    y[1, ] <- lambda * x[1, ]
    
    xi_yt <- tcrossprod(
      x[1, ] - y[1, ]
    )
    
    S[, , 1] <-
      ct * (
        w * xi_yt +
          (1 - w) * xi_yt
      )
    
    S1[, , 1] <-
      w * xi_yt +
      (1 - w) * xi_yt
    
    d <- diag(xi_yt) - 1
    
    D1[1] <- sum(abs(d))
    D2[1] <- sum(d^2)
    D3[1] <- sum(diag(S1[, , 1]))
    
    # --------------------------------------------------------
    # Remaining observations
    # --------------------------------------------------------
    
    for (i in 2:n) {
      
      y[i, ] <-
        lambda * x[i, ] +
        (1 - lambda) * y[i - 1, ]
      
      xi_yt <- tcrossprod(
        x[i, ] - y[i, ]
      )
      
      S[, , i] <-
        ct * (
          w * xi_yt +
            (1 - w) * S[, , i - 1]
        )
      
      S1[, , i] <-
        w * xi_yt +
        (1 - w) * S1[, , i - 1]
      
      d <- diag(S[, , i]) - 1
      
      D1[i] <- sum(abs(d))
      D2[i] <- sum(d^2)
      D3[i] <- sum(diag(S1[, , i]))
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