# ============================================================
# 05_ARL_Calculation.R
# ARL Simulation for MEWMS-Based Control Charts
# ============================================================
library(MASS)
library(pcaPP)
library(Matrix)
library(robustbase)
library(rrcov)
library(LaplacesDemon)
library(mvtnorm)


# m      - Maximum run length
# n      - Number of Phase I samples
# iter   - Number of simulation replications
# a1,a2  - Phase I standard deviations
# b1,b2  - Phase II standard deviations
# rho1   - Phase I correlation
# rho2   - Phase II correlation
# per    - Percentage of contaminated observations
# w      - Smoothing constant
# cv     - Control limits in the order:
#          Comedian, Sn, Classical, MCD, MVE


BivIndCCsim_combinedARL1 <- function(
    m, n, iter,
    a1, a2, b1, b2,
    rho1, rho2, per, w, cv
) {
  
  om1 <- 0
  om2 <- 0
  
  p <- 2
  
  RLS  <- numeric(iter)
  RLCl <- numeric(iter)
  RLC  <- numeric(iter)
  RLMc <- numeric(iter)
  RLMv <- numeric(iter)
  
  c1 <- diag(p)
  
  c2 <- matrix(
    c(
      b1^2,
      b1 * b2 * rho2,
      b1 * b2 * rho2,
      b2^2
    ),
    2, 2
  )
  
  C <- matrix(
    c(
      a1^2,
      a1 * a2 * rho1,
      a1 * a2 * rho1,
      a2^2
    ),
    2, 2
  )
  
  for (j in 1:iter) {
    
    n1 <- round(n * per / 100)
    
    wt <- numeric(n)
    
    y1 <- mvrnorm(
      (n - n1),
      mu = rep(0, p),
      Sigma = c1
    )
    
    if (n1 == 0) {
      
      y <- y1
      
    } else {
      
      # y2 <- mvrnorm(
      #   n1,
      #   rep(om1, p),
      #   Sigma = C
      # )
      
      # y2 <- rmvt(
      #   n1,
      #   sigma = C,
      #   df = 3,
      #   delta = rep(om1, p),
      #   type = "shifted"
      # )
      
      # y2 <- rmvpe(
      #   n1,
      #   mu = rep(0, p),
      #   Sigma = C,
      #   kappa = 0.5
      # )
      
      # Multivariate skewed normal
      y2 <- rmsn(
        n = n1,
        xi = c(0, 0),
        Omega = C,
        alpha = c(5, 5)
      )
      
      y <- rbind(y1, y2)
    }
    
    MR <- colMedians(y)
    MC <- colMeans(y)
    
    CCC <- Comed(y)
    CSn <- covSn(y)
    
    kk1 <- CovMcd(y)
    CMc <- kk1$cov
    
    kk2 <- CovMve(y)
    CMv <- kk2$cov
    
    CCl <- cov(y)
    
    g <- mvrnorm(
      m,
      mu = rep(om2, p),
      Sigma = c2
    )
    
    
    # --------------------------------------------------------
    # Comedian
    # --------------------------------------------------------
    
    eig <- eigen(CCC)
    
    Si <- eig$vectors %*%
      diag(1 / sqrt(eig$values)) %*%
      t(eig$vectors)
    
    x <- sweep(
      g,
      2,
      MR
    ) %*% Si
    
    cv1 <- cv[1]
    
    S <- array(
      0,
      dim = c(p, p, m)
    )
    
    D1 <- numeric(m)
    D2 <- numeric(m)
    D3 <- numeric(m)
    
    r1 <- 1
    
    xi_xt <- tcrossprod(x[1, ])
    
    S[, , 1] <- xi_xt
    
    d <- diag(xi_xt) - 1
    
    D3[1] <- sum(
      diag(S[, , 1])
    )
    
    RLC[j] <- r1
    
    if (is.finite(D3[1]) &&
        D3[1] <= cv1) {
      
      for (i in 2:m) {
        
        xi_xt <- tcrossprod(x[i, ])
        
        S[, , i] <-
          w * xi_xt +
          (1 - w) * S[, , (i - 1)]
        
        d <- diag(S[, , i]) - 1
        
        D3[i] <- sum(
          diag(S[, , i])
        )
        
        r1 <- r1 + 1
        
        RLC[j] <- r1
        
        if (is.finite(D3[i]) &&
            D3[i] > cv1) {
          break
        }
      }
    }
    
    
    # --------------------------------------------------------
    # Sn
    # --------------------------------------------------------
    
    eig <- eigen(CSn)
    
    Si <- eig$vectors %*%
      diag(1 / sqrt(eig$values)) %*%
      t(eig$vectors)
    
    x <- sweep(
      g,
      2,
      MR
    ) %*% Si
    
    cv2 <- cv[2]
    
    S <- array(
      0,
      dim = c(p, p, m)
    )
    
    D1 <- numeric(m)
    D2 <- numeric(m)
    D3 <- numeric(m)
    
    r2 <- 1
    
    xi_xt <- tcrossprod(x[1, ])
    
    S[, , 1] <- xi_xt
    
    d <- diag(xi_xt) - 1
    
    D3[1] <- sum(
      diag(S[, , 1])
    )
    
    RLS[j] <- r2
    
    if (is.finite(D3[1]) &&
        D3[1] <= cv2) {
      
      for (i in 2:m) {
        
        xi_xt <- tcrossprod(x[i, ])
        
        S[, , i] <-
          w * xi_xt +
          (1 - w) * S[, , (i - 1)]
        
        d <- diag(S[, , i]) - 1
        
        D3[i] <- sum(
          diag(S[, , i])
        )
        
        r2 <- r2 + 1
        
        RLS[j] <- r2
        
        if (is.finite(D3[i]) &&
            D3[i] > cv2) {
          break
        }
      }
    }
    
    
    # --------------------------------------------------------
    # Classical
    # --------------------------------------------------------
    
    eig <- eigen(CCl)
    
    Si <- eig$vectors %*%
      diag(1 / sqrt(eig$values)) %*%
      t(eig$vectors)
    
    x <- sweep(
      g,
      2,
      MC
    ) %*% Si
    
    cv3 <- cv[3]
    
    S <- array(
      0,
      dim = c(p, p, m)
    )
    
    D1 <- numeric(m)
    D2 <- numeric(m)
    D3 <- numeric(m)
    
    r3 <- 1
    
    xi_xt <- tcrossprod(x[1, ])
    
    S[, , 1] <- xi_xt
    
    d <- diag(xi_xt) - 1
    
    D3[1] <- sum(
      diag(S[, , 1])
    )
    
    RLCl[j] <- r3
    
    if (is.finite(D3[1]) &&
        D3[1] <= cv3) {
      
      for (i in 2:m) {
        
        xi_xt <- tcrossprod(x[i, ])
        
        S[, , i] <-
          w * xi_xt +
          (1 - w) * S[, , (i - 1)]
        
        d <- diag(S[, , i]) - 1
        
        D3[i] <- sum(
          diag(S[, , i])
        )
        
        r3 <- r3 + 1
        
        RLCl[j] <- r3
        
        if (is.finite(D3[i]) &&
            D3[i] > cv3) {
          break
        }
      }
    }
    
    
    # --------------------------------------------------------
    # MCD
    # --------------------------------------------------------
    
    eig <- eigen(CMc)
    
    Si <- eig$vectors %*%
      diag(1 / sqrt(eig$values)) %*%
      t(eig$vectors)
    
    x <- sweep(
      g,
      2,
      MR
    ) %*% Si
    
    cv4 <- cv[4]
    
    S <- array(
      0,
      dim = c(p, p, m)
    )
    
    D1 <- numeric(m)
    D2 <- numeric(m)
    D3 <- numeric(m)
    
    r4 <- 1
    
    xi_xt <- tcrossprod(x[1, ])
    
    S[, , 1] <- xi_xt
    
    d <- diag(xi_xt) - 1
    
    D3[1] <- sum(
      diag(S[, , 1])
    )
    
    RLMc[j] <- r4
    
    if (is.finite(D3[1]) &&
        D3[1] <= cv4) {
      
      for (i in 2:m) {
        
        xi_xt <- tcrossprod(x[i, ])
        
        S[, , i] <-
          w * xi_xt +
          (1 - w) * S[, , (i - 1)]
        
        d <- diag(S[, , i]) - 1
        
        D3[i] <- sum(
          diag(S[, , i])
        )
        
        r4 <- r4 + 1
        
        RLMc[j] <- r4
        
        if (is.finite(D3[i]) &&
            D3[i] > cv4) {
          break
        }
      }
    }
    
    
    # --------------------------------------------------------
    # MVE
    # --------------------------------------------------------
    
    eig <- eigen(CMv)
    
    Si <- eig$vectors %*%
      diag(1 / sqrt(eig$values)) %*%
      t(eig$vectors)
    
    x <- sweep(
      g,
      2,
      MR
    ) %*% Si
    
    cv5 <- cv[5]
    
    S <- array(
      0,
      dim = c(p, p, m)
    )
    
    D1 <- numeric(m)
    D2 <- numeric(m)
    D3 <- numeric(m)
    
    r5 <- 1
    
    xi_xt <- tcrossprod(x[1, ])
    
    S[, , 1] <- xi_xt
    
    d <- diag(xi_xt) - 1
    
    D3[1] <- sum(
      diag(S[, , 1])
    )
    
    RLMv[j] <- r5
    
    if (is.finite(D3[1]) &&
        D3[1] <= cv5) {
      
      for (i in 2:m) {
        
        xi_xt <- tcrossprod(x[i, ])
        
        S[, , i] <-
          w * xi_xt +
          (1 - w) * S[, , (i - 1)]
        
        d <- diag(S[, , i]) - 1
        
        D3[i] <- sum(
          diag(S[, , i])
        )
        
        r5 <- r5 + 1
        
        RLMv[j] <- r5
        
        if (is.finite(D3[i]) &&
            D3[i] > cv5) {
          break
        }
      }
    }
  }
  
  ARL <- c(
    ARLC  = mean(RLC),
    ARLS  = mean(RLS),
    ARLCl = mean(RLCl),
    ARLMc = mean(RLMc),
    ARLMv = mean(RLMv)
  )
  
  return(ARL)
}


# ============================================================
# Supporting robust covariance functions
# ============================================================

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


Comed <- function(w) {
  
  x <- w[, 1]
  y <- w[, 2]
  
  n <- length(x)
  
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


gkcov <- function(w, method) {
  
  x <- w[, 1] + w[, 2]
  y <- w[, 1] - w[, 2]
  
  switch(
    method,
    
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
    }
  )
  
  return(gkc)
}