# ============================================================
# 04_MEWMV_Plot.R
# Plotting MEWMV Statistics
# ============================================================

library(MASS)
library(robustbase)
library(rrcov)
library(pcaPP)

# ------------------------------------------------------------
# Function: RobDisplot_MShift
#
# x1     : Phase I data
# x2     : Phase II data
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
#
# st     : monitoring statistic
#          1 = D1
#          2 = D2
#          3 = D3
#
# cv     : estimated UCL
# ------------------------------------------------------------

RobDisplot_MShift <- function(
    x1, x2, w, lambda, method, st, cv) {
  
  g <- x1
  
  p <- ncol(x1)
  N <- nrow(x2)
  
  # ----------------------------------------------------------
  # Standardization using Phase I estimates
  # ----------------------------------------------------------
  
  switch(
    method,
    
    # --------------------------------------------------------
    # Comedian
    # --------------------------------------------------------
    m1 = {
      title <- "Comedian control chart"
      
      K1 <- Comed(g)
      
      eig <- eigen(K1)
      
      Si <- eig$vectors %*%
        diag(1 / sqrt(eig$values)) %*%
        t(eig$vectors)
      
      x <- sweep(
        x2,
        2,
        colMedians(g)
      ) %*% Si
    },
    
    # --------------------------------------------------------
    # Sn
    # --------------------------------------------------------
    m2 = {
      title <- "Sn control chart"
      
      K1 <- covSn(g)
      
      eig <- eigen(K1)
      
      Si <- eig$vectors %*%
        diag(1 / sqrt(eig$values)) %*%
        t(eig$vectors)
      
      x <- sweep(
        x2,
        2,
        colMedians(g)
      ) %*% Si
    },
    
    # --------------------------------------------------------
    # Classical
    # --------------------------------------------------------
    m3 = {
      title <- "Classical control chart"
      
      eig <- eigen(cov(g))
      
      Si <- eig$vectors %*%
        diag(1 / sqrt(eig$values)) %*%
        t(eig$vectors)
      
      x <- sweep(
        x2,
        2,
        colMeans(g)
      ) %*% Si
    },
    
    # --------------------------------------------------------
    # MCD
    # --------------------------------------------------------
    m4 = {
      title <- "MCD control chart"
      
      K1 <- covMcd(g)
      
      eig <- eigen(K1$cov)
      
      Si <- eig$vectors %*%
        diag(1 / sqrt(eig$values)) %*%
        t(eig$vectors)
      
      x <- sweep(
        x2,
        2,
        K1$center
      ) %*% Si
    },
    
    # --------------------------------------------------------
    # MVE
    # --------------------------------------------------------
    m5 = {
      title <- "MVE control chart"
      
      K2 <- CovMve(g)
      
      eig <- eigen(K2$cov)
      
      Si <- eig$vectors %*%
        diag(1 / sqrt(eig$values)) %*%
        t(eig$vectors)
      
      x <- sweep(
        x2,
        2,
        K2$center
      ) %*% Si
    },
    
    # --------------------------------------------------------
    # GK(Q)
    # --------------------------------------------------------
    m6 = {
      title <- "GK_Qn control chart"
      
      K1 <- gkcov(g, "Q")
      
      eig <- eigen(K1)
      
      Si <- eig$vectors %*%
        diag(1 / sqrt(eig$values)) %*%
        t(eig$vectors)
      
      x <- sweep(
        x2,
        2,
        colMedians(g)
      ) %*% Si
    },
    
    # --------------------------------------------------------
    # GK(S)
    # --------------------------------------------------------
    m7 = {
      title <- "GK_Sn control chart"
      
      K1 <- gkcov(g, "S")
      
      eig <- eigen(K1)
      
      Si <- eig$vectors %*%
        diag(1 / sqrt(eig$values)) %*%
        t(eig$vectors)
      
      x <- sweep(
        x2,
        2,
        colMedians(g)
      ) %*% Si
    },
    
    stop("Invalid method. Use m1, m2, m3, m4, m5, m6, or m7.")
  )
  
  # ----------------------------------------------------------
  # MEWMV statistics
  # ----------------------------------------------------------
  
  S <- array(
    0,
    dim = c(p, p, N)
  )
  
  S1 <- array(
    0,
    dim = c(p, p, N)
  )
  
  D1 <- numeric(N)
  D2 <- numeric(N)
  D3 <- numeric(N)
  
  y <- matrix(
    0,
    N,
    p
  )
  
  # ----------------------------------------------------------
  # First observation
  # ----------------------------------------------------------
  
  y[1, ] <- lambda * x[1, ]
  
  xi_yt <- tcrossprod(
    x[1, ] - y[1, ]
  )
  
  S[, , 1] <-
    w * xi_yt +
    (1 - w) * xi_yt
  
  S1[, , 1] <-
    w * xi_yt +
    (1 - w) * xi_yt
  
  d <- diag(xi_yt) - 1
  
  D1[1] <- sum(abs(d))
  D2[1] <- sum(d^2)
  
  # IMPORTANT:
  # D3 is based on S1, as in the UCL-estimation code.
  D3[1] <- sum(diag(S1[, , 1]))
  
  # ----------------------------------------------------------
  # Remaining observations
  # ----------------------------------------------------------
  
  for (i in 2:N) {
    
    y[i, ] <-
      lambda * x[i, ] +
      (1 - lambda) * y[i - 1, ]
    
    xi_yt <- tcrossprod(
      x[i, ] - y[i, ]
    )
    
    S[, , i] <-
      w * xi_yt +
      (1 - w) * S[, , i - 1]
    
    S1[, , i] <-
      w * xi_yt +
      (1 - w) * S1[, , i - 1]
    
    d <- diag(S[, , i]) - 1
    
    D1[i] <- sum(abs(d))
    D2[i] <- sum(d^2)
    
    # IMPORTANT:
    # D3 is based on S1, not S.
    D3[i] <- sum(diag(S1[, , i]))
  }
  
  # ----------------------------------------------------------
  # Plot selected MEWMV statistic
  # ----------------------------------------------------------
  
  if (st == 1) {
    
    plot(
      seq_len(N),
      D1,
      type = "l",
      main = title,
      xlab = "Sample number",
      ylab = expression("MEWMV" * L[1]),
      cex.main = 0.9
    )
    
    abline(
      h = cv,
      col = "red",
      lty = 1
    )
    
    return(D1)
    
  } else if (st == 2) {
    
    plot(
      seq_len(N),
      D2,
      type = "l",
      main = title,
      xlab = "Sample number",
      ylab = expression("MEWMV" * L[2]),
      cex.main = 0.9
    )
    
    abline(
      h = cv,
      col = "red",
      lty = 1
    )
    
    return(D2)
    
  } else if (st == 3) {
    
    plot(
      seq_len(N),
      D3,
      type = "l",
      main = title,
      xlab = "Sample number",
      ylab = "MEWMV",
      cex.main = 0.9
    )
    
    abline(
      h = cv,
      col = "red",
      lty = 1
    )
    
    return(D3)
    
  } else {
    
    stop("st must be 1, 2, or 3.")
  }
}