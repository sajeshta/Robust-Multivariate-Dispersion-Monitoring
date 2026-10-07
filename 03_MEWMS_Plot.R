# ============================================================
# 03_MEWMS_Plot.R
# Plotting MEWMS Statistics
# ============================================================

library(MASS)
library(robustbase)
library(rrcov)

# ------------------------------------------------------------
# Function: RobDisplot_MEWMS
#
# x1 : Phase I data
# x2 : Phase II data
# w  : smoothing parameter
# method : covariance estimator
# st : monitoring statistic
#       1 = D1
#       2 = D2
#       3 = D3
# cv : estimated UCL
# ------------------------------------------------------------

RobDisplot_MEWMS <- function(
    x1,
    x2,
    w,
    method,
    st,
    cv) {
  
  g <- x1
  
  n <- nrow(x1)
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
    }
  )
  
  # ----------------------------------------------------------
  # MEWMS statistics
  # ----------------------------------------------------------
  
  S <- array(
    0,
    dim = c(p, p, N)
  )
  
  D1 <- numeric(N)
  D2 <- numeric(N)
  D3 <- numeric(N)
  
  for (i in 1:N) {
    
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
  
  # ----------------------------------------------------------
  # Plot selected MEWMS statistic
  # ----------------------------------------------------------
  
  if (st == 1) {
    
    plot(
      seq_len(N),
      D1,
      type = "l",
      main = title,
      xlab = "Sample number",
      ylab = expression("MEWMS" * L[1]),
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
      ylab = expression("MEWMS" * L[2]),
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
      ylab = "MEWMS",
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