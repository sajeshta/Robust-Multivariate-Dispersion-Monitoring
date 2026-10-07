# ============================================================
# 05_ARL_Calculation.R
# ARL Simulation for MEWMV-Based Control Charts
# ============================================================

# ------------------------------------------------------------
# Required packages
# ------------------------------------------------------------

library(MASS)
library(pcaPP)
library(Matrix)
library(robustbase)
library(beepr)
library(rrcov)


# ============================================================
# 1. MEWMV-based ARL calculation
# ============================================================

# m      - Maximum run length
# n      - Number of Phase I samples
# iter   - Number of simulation replications
# a1,a2  - Phase I standard deviations
# b1,b2  - Phase II standard deviations
# rho1   - Phase I correlation
# rho2   - Phase II correlation
# per    - Percentage of contaminated observations
# w      - Smoothing constant
# lambda - Mean-shift smoothing parameter
# om1    - Phase I location parameter
# om2    - Phase II location parameter
# cv     - Control limits in the order:
#          Comedian, Sn, Classical, MCD, MVE, GK_Qn, GK_Sn


BivIndCC_meanshift_combARL1 <- function(
    m, n, iter, a1, a2, b1, b2,
    rho1, rho2, per, w, lambda,
    om1, om2, cv
) {
  
  if (length(cv) < 7 ||
      any(is.na(cv)) ||
      any(!is.finite(cv))) {
    
    stop(
      "cv must contain 7 finite control-limit values in the order: ",
      "Comedian, Sn, Classical, MCD, MVE, GK_Qn, GK_Sn"
    )
  }
  
  p <- 2
  
  ct <- (2 - lambda) /
    (2 * ((1 - lambda)^2))
  
  RLS  <- numeric(iter)
  RLCl <- numeric(iter)
  RLC  <- numeric(iter)
  RLMc <- numeric(iter)
  RLMv <- numeric(iter)
  RLGq <- numeric(iter)
  RLGs <- numeric(iter)
  
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
      
      # Multivariate Exponential
      y2 <- rmvpe(
        n1,
        mu = rep(0, p),
        Sigma = C,
        kappa = 0.5
      )
      
      y <- rbind(y1, y2)
    }
    
    MR  <- colMedians(y)
    MC  <- colMeans(y)
    
    CCC <- Comed(y)       # Comedian
    CSn <- covSn(y)       # Sn
    
    kk1 <- CovMcd(y)      # MCD
    CMc <- kk1$cov
    
    kk2 <- CovMve(y)      # MVE
    CMv <- kk2$cov
    
    CGq <- gkcov(y, "Q")  # GK(Q)
    CGs <- gkcov(y, "S")  # GK(S)
    
    CCl <- cov(y)         # Classical
    
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
    
    S <- array(
      0,
      dim = c(p, p, m)
    )
    
    D1 <- numeric(m)
    D2 <- numeric(m)
    D3 <- numeric(m)
    
    r1 <- 1
    
    z <- matrix(
      0,
      m,
      p
    )
    
    z[1, ] <- lambda * x[1, ]
    
    xi_zt <- tcrossprod(
      x[1, ] - z[1, ]
    )
    
    S[, , 1] <-
      w * xi_zt +
      (1 - w) * xi_zt
    
    d <- diag(xi_zt) - 1
    
    D3[1] <- sum(
      diag(S[, , 1])
    )
    
    RLC[j] <- r1
    
    if (is.finite(D3[i]) &&
        D3[i] <= cv[1]) {
      
      for (i in 2:m) {
        
        z[i, ] <-
          lambda * x[i, ] +
          (1 - lambda) * z[(i - 1), ]
        
        xi_zt <- tcrossprod(
          x[i, ] - z[i, ]
        )
        
        S[, , i] <-
          w * xi_zt +
          (1 - w) * S[, , (i - 1)]
        
        d <- diag(S[, , i]) - 1
        
        D3[i] <- sum(
          diag(S[, , i])
        )
        
        r1 <- r1 + 1
        
        RLC[j] <- r1
        
        if (is.finite(D3[i]) &&
            D3[i] > cv[1]) {
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
    
    S <- array(
      0,
      dim = c(p, p, m)
    )
    
    D1 <- numeric(m)
    D2 <- numeric(m)
    D3 <- numeric(m)
    
    r2 <- 1
    
    z <- matrix(
      0,
      m,
      p
    )
    
    z[1, ] <- lambda * x[1, ]
    
    xi_zt <- tcrossprod(
      x[1, ] - z[1, ]
    )
    
    S[, , 1] <-
      w * xi_zt +
      (1 - w) * xi_zt
    
    d <- diag(xi_zt) - 1
    
    D3[1] <- sum(
      diag(S[, , 1])
    )
    
    RLS[j] <- r2
    
    if (is.finite(D3[i]) &&
        D3[i] <= cv[2]) {
      
      for (i in 2:m) {
        
        z[i, ] <-
          lambda * x[i, ] +
          (1 - lambda) * z[(i - 1), ]
        
        xi_zt <- tcrossprod(
          x[i, ] - z[i, ]
        )
        
        S[, , i] <-
          w * xi_zt +
          (1 - w) * S[, , (i - 1)]
        
        d <- diag(S[, , i]) - 1
        
        D3[i] <- sum(
          diag(S[, , i])
        )
        
        r2 <- r2 + 1
        
        RLS[j] <- r2
        
        if (is.finite(D3[i]) &&
            D3[i] > cv[2]) {
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
    
    S <- array(
      0,
      dim = c(p, p, m)
    )
    
    D1 <- numeric(m)
    D2 <- numeric(m)
    D3 <- numeric(m)
    
    r3 <- 1
    
    z <- matrix(
      0,
      m,
      p
    )
    
    z[1, ] <- lambda * x[1, ]
    
    xi_zt <- tcrossprod(
      x[1, ] - z[1, ]
    )
    
    S[, , 1] <-
      w * xi_zt +
      (1 - w) * xi_zt
    
    d <- diag(xi_zt) - 1
    
    D3[1] <- sum(
      diag(S[, , 1])
    )
    
    RLCl[j] <- r3
    
    if (is.finite(D3[i]) &&
        D3[i] <= cv[3]) {
      
      for (i in 2:m) {
        
        z[i, ] <-
          lambda * x[i, ] +
          (1 - lambda) * z[(i - 1), ]
        
        xi_zt <- tcrossprod(
          x[i, ] - z[i, ]
        )
        
        S[, , i] <-
          w * xi_zt +
          (1 - w) * S[, , (i - 1)]
        
        d <- diag(S[, , i]) - 1
        
        D3[i] <- sum(
          diag(S[, , i])
        )
        
        r3 <- r3 + 1
        
        RLCl[j] <- r3
        
        if (is.finite(D3[i]) &&
            D3[i] > cv[3]) {
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
    
    S <- array(
      0,
      dim = c(p, p, m)
    )
    
    D1 <- numeric(m)
    D2 <- numeric(m)
    D3 <- numeric(m)
    
    r4 <- 1
    
    z <- matrix(
      0,
      m,
      p
    )
    
    z[1, ] <- lambda * x[1, ]
    
    xi_zt <- tcrossprod(
      x[1, ] - z[1, ]
    )
    
    S[, , 1] <-
      w * xi_zt +
      (1 - w) * xi_zt
    
    d <- diag(xi_zt) - 1
    
    D3[1] <- sum(
      diag(S[, , 1])
    )
    
    RLMc[j] <- r4
    
    if (is.finite(D3[i]) &&
        D3[i] <= cv[4]) {
      
      for (i in 2:m) {
        
        z[i, ] <-
          lambda * x[i, ] +
          (1 - lambda) * z[(i - 1), ]
        
        xi_zt <- tcrossprod(
          x[i, ] - z[i, ]
        )
        
        S[, , i] <-
          w * xi_zt +
          (1 - w) * S[, , (i - 1)]
        
        d <- diag(S[, , i]) - 1
        
        D3[i] <- sum(
          diag(S[, , i])
        )
        
        r4 <- r4 + 1
        
        RLMc[j] <- r4
        
        if (is.finite(D3[i]) &&
            D3[i] > cv[4]) {
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
    
    S <- array(
      0,
      dim = c(p, p, m)
    )
    
    D1 <- numeric(m)
    D2 <- numeric(m)
    D3 <- numeric(m)
    
    r5 <- 1
    
    z <- matrix(
      0,
      m,
      p
    )
    
    z[1, ] <- lambda * x[1, ]
    
    xi_zt <- tcrossprod(
      x[1, ] - z[1, ]
    )
    
    S[, , 1] <-
      w * xi_zt +
      (1 - w) * xi_zt
    
    d <- diag(xi_zt) - 1
    
    D3[1] <- sum(
      diag(S[, , 1])
    )
    
    RLMv[j] <- r5
    
    if (is.finite(D3[i]) &&
        D3[i] <= cv[5]) {
      
      for (i in 2:m) {
        
        z[i, ] <-
          lambda * x[i, ] +
          (1 - lambda) * z[(i - 1), ]
        
        xi_zt <- tcrossprod(
          x[i, ] - z[i, ]
        )
        
        S[, , i] <-
          w * xi_zt +
          (1 - w) * S[, , (i - 1)]
        
        d <- diag(S[, , i]) - 1
        
        D3[i] <- sum(
          diag(S[, , i])
        )
        
        r5 <- r5 + 1
        
        RLMv[j] <- r5
        
        if (is.finite(D3[i]) &&
            D3[i] > cv[5]) {
          break
        }
      }
    }
    
    
    # --------------------------------------------------------
    # GK(Q)
    # --------------------------------------------------------
    
    eig <- eigen(CGq)
    
    Si <- eig$vectors %*%
      diag(1 / sqrt(eig$values)) %*%
      t(eig$vectors)
    
    x <- sweep(
      g,
      2,
      MR
    ) %*% Si
    
    S <- array(
      0,
      dim = c(p, p, m)
    )
    
    D1 <- numeric(m)
    D2 <- numeric(m)
    D3 <- numeric(m)
    
    r6 <- 1
    
    z <- matrix(
      0,
      m,
      p
    )
    
    z[1, ] <- lambda * x[1, ]
    
    xi_zt <- tcrossprod(
      x[1, ] - z[1, ]
    )
    
    S[, , 1] <-
      w * xi_zt +
      (1 - w) * xi_zt
    
    d <- diag(xi_zt) - 1
    
    D3[1] <- sum(
      diag(S[, , 1])
    )
    
    RLGq[j] <- r6
    
    if (is.finite(D3[i]) &&
        D3[i] <= cv[6]) {
      
      for (i in 2:m) {
        
        z[i, ] <-
          lambda * x[i, ] +
          (1 - lambda) * z[(i - 1), ]
        
        xi_zt <- tcrossprod(
          x[i, ] - z[i, ]
        )
        
        S[, , i] <-
          w * xi_zt +
          (1 - w) * S[, , (i - 1)]
        
        d <- diag(S[, , i]) - 1
        
        D3[i] <- sum(
          diag(S[, , i])
        )
        
        r6 <- r6 + 1
        
        RLGq[j] <- r6
        
        if (is.finite(D3[i]) &&
            D3[i] > cv[6]) {
          break
        }
      }
    }
    
    
    # --------------------------------------------------------
    # GK(S)
    # --------------------------------------------------------
    
    eig <- eigen(CGs)
    
    Si <- eig$vectors %*%
      diag(1 / sqrt(eig$values)) %*%
      t(eig$vectors)
    
    x <- sweep(
      g,
      2,
      MR
    ) %*% Si
    
    S <- array(
      0,
      dim = c(p, p, m)
    )
    
    D1 <- numeric(m)
    D2 <- numeric(m)
    D3 <- numeric(m)
    
    r7 <- 1
    
    z <- matrix(
      0,
      m,
      p
    )
    
    z[1, ] <- lambda * x[1, ]
    
    xi_zt <- tcrossprod(
      x[1, ] - z[1, ]
    )
    
    S[, , 1] <-
      w * xi_zt +
      (1 - w) * xi_zt
    
    d <- diag(xi_zt) - 1
    
    D3[1] <- sum(
      diag(S[, , 1])
    )
    
    RLGs[j] <- r7
    
    if (is.finite(D3[i]) &&
        D3[i] <= cv[7]) {
      
      for (i in 2:m) {
        
        z[i, ] <-
          lambda * x[i, ] +
          (1 - lambda) * z[(i - 1), ]
        
        xi_zt <- tcrossprod(
          x[i, ] - z[i, ]
        )
        
        S[, , i] <-
          w * xi_zt +
          (1 - w) * S[, , (i - 1)]
        
        d <- diag(S[, , i]) - 1
        
        D3[i] <- sum(
          diag(S[, , i])
        )
        
        r7 <- r7 + 1
        
        RLGs[j] <- r7
        
        if (is.finite(D3[i]) &&
            D3[i] > cv[7]) {
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
    ARLMv = mean(RLMv),
    ARLGq = mean(RLGq),
    ARLGs = mean(RLGs)
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


