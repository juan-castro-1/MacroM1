VAR1.coef <- function(VAR, m, p) {
  if (p == 1) {
    Pi <- Bcoef(VAR)[, (m * p + 1)]
    Phi <- Bcoef(VAR)[, 1:(m * p)]
  } else if (p > 1) {
    Pi <- c(Bcoef(VAR)[, (m * p + 1)], rep(0, m * (p - 1)))
    Phi <- rbind(Bcoef(VAR)[, 1:(m * p)], cbind(diag(m * (p - 1)), matrix(0, nrow = m * (p - 1), ncol = m)))
  }
  coefficients <- list(Pi = Pi, Phi = Phi)
  return(coefficients)
}

sirf <- function(Phi, P, m, H, cumulative = FALSE) {
  Theta <- array(NA, c(m, m, H + 1))
  Theta[, , 1] <- P
  tmp_Phi <- Phi
  for (h in 1:H) {
    Theta[, , h + 1] <- tmp_Phi[1:m, 1:m] %*% P
    tmp_Phi <- tmp_Phi %*% Phi
  }
  if (cumulative == FALSE) {
    return(Theta)
  } else {
    return(aperm(apply(Theta, c(1, 2), cumsum), c(2, 3, 1)))
  }
}

fevd <- function(Phi, P, m, H) {
  Theta <- sirf(Phi, P, m, H)
  Omega <- array(NA, c(m, m, H + 1))
  for (h in 1:(H + 1)) {
    for (j in 1:m) {
      for (i in 1:m) {
        Omega[i, j, h] <- sum(Theta[i, j, 1:h] ^ 2)
      }
    }
  }
  MSE <- apply(Omega, c(1, 3), sum)
  for (h in 1:(H + 1)) {
    for (j in 1:m) {
      for (i in 1:m) {
        Omega[i, j, h] <- Omega[i, j, h] / MSE[i, h]
      }
    }
  }
  return(100 * Omega)
}

hd <- function(Phi, P, Y, E, m, N) {
  Theta <- sirf(Phi, P, m, N)
  U <- t(solve(P, t(E)))
  HD <- array(NA, c(m, m + 1, N))
  for (t in 1:N) {
    for (j in 1:m) {
      for (i in 1:m) {
        HD[i, j, t] <- c(crossprod(Theta[i, j, 1:t], U[t:1, j]))
      }
    }
  }
  NS <- Y - t(apply(HD[, 1:m, ], c(1, 3), sum))
  for (i in 1:m) {
    HD[i, m + 1, ] <- NS[, i]
  }
  return(HD)
}

erpt <- function(Phi, P, m, H, vx, vy, cumulative = FALSE) {
  if (cumulative == FALSE) {
    I <- sirf(Phi, P, m, H)
  } else {
    I <- sirf(Phi, P, m, H, cumulative = TRUE)
  }
  return(100 * I[vx, vy, ] / I[vy, vy, ])
}

SVAR.sirf <- function(SVAR, H, cumulative = FALSE) {
  RF <- VAR1.coef(SVAR$var, SVAR$var$K, SVAR$var$p)
  if (cumulative == FALSE) {
    SIRF <- sirf(RF$Phi, solve(SVAR$A, SVAR$B), SVAR$var$K, H)
  } else {
    SIRF <- sirf(RF$Phi, solve(SVAR$A, SVAR$B), SVAR$var$K, H, cumulative = TRUE)
  }
  dimnames(SIRF)[[1]] <- toupper(colnames(SVAR$var$y))
  dimnames(SIRF)[[2]] <- paste("S.", 1:m, sep = "")
  dimnames(SIRF)[[3]] <- 0:H
  return(SIRF)
}

SVAR.fevd <- function(SVAR, H) {
  RF <- VAR1.coef(SVAR$var, SVAR$var$K, SVAR$var$p)
  FEVD <- fevd(RF$Phi, solve(SVAR$A, SVAR$B), SVAR$var$K, H)
  dimnames(FEVD)[[1]] <- toupper(colnames(SVAR$var$y))
  dimnames(FEVD)[[2]] <- paste("S.", 1:m, sep = "")
  dimnames(FEVD)[[3]] <- 0:H
  return(FEVD)
}

SVAR.hd <- function(SVAR) {
  RF <- VAR1.coef(SVAR$var, SVAR$var$K, SVAR$var$p)
  HD <- hd(RF$Phi, solve(SVAR$A, SVAR$B), SVAR$var$datamat[, 1:SVAR$var$K], resid(SVAR$var), SVAR$var$K, SVAR$var$obs)
  dimnames(HD)[[1]] <- toupper(colnames(SVAR$var$y))
  dimnames(HD)[[2]] <- c(paste("S.", 1:m, sep = ""), "NS")
  dimnames(HD)[[3]] <- 1:N
  return(HD)
}

SVAR.erpt <- function(SVAR, H, vx, vy, cumulative = FALSE) {
  RF <- VAR1.coef(SVAR$var, SVAR$var$K, SVAR$var$p)
  if (cumulative == FALSE) {
    ERPT <- erpt(RF$Phi, solve(SVAR$A, SVAR$B), SVAR$var$K, H, vx, vy)
  } else {
    ERPT <- erpt(RF$Phi, solve(SVAR$A, SVAR$B), SVAR$var$K, H, vx, vy, cumulative = TRUE)
  }
  return(ERPT)
}