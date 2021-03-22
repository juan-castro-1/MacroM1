boot.ci.efron <- function(x, a) {
  lb <- (1 - a) / 2
  ub <- 1 - (1 - a) / 2
  ci <- quantile(x,  probs = c(lb, ub))
  return(ci)
}

boot.rb.resample <- function(Y0, E, Pi, Phi, m, p, pmax, N) {
  Yr <- array(NA, c(pmax + N, m))
  Yr[1:pmax, ] <- Y0
  for (t in (pmax + 1):(pmax + N)) {
    i <- sample(1:N, 1, replace = TRUE)
    e <- E[i, ]
    if (p == 1) {
      y1 <- Yr[(t - 1), ]
    } else {
      y1 <- t(Yr[(t - p):(t - 1), ])
      y1 <- c(y1[, rev(1:p)])
    }
    Yr[t, ] <- Pi[1:m] + Phi[1:m, 1:(m * p)] %*% y1 + e
  }
  return(Yr)
}

boot.rb.replicate <- function(VAR, Y0, pmax, R) {
  Yb <- array(NA, c(pmax + VAR$obs, VAR$K, R))
  RF <- VAR1.coef(VAR, VAR$K, VAR$p)
  RF.E <- resid(VAR)
  for (r in 1:R) {
    Yb[, , r] <- boot.rb.resample(Y0, RF.E, RF$Pi, RF$Phi, VAR$K, VAR$p, pmax, VAR$obs)
  }
  return(Yb)
}

boot.rb.VAR <- function(nvars, Yr, m, p, pmax) {
  # pr <- VARselect(Yr, lag.max = pmax, type = "const")
  # pr <- pr$selection[1] # AIC
  pr <- p # Fixed
  Yr <- Yr[(pmax - pr + 1):dim(Yr)[1], ]
  colnames(Yr) <- nvars
  VAR.rep <- VAR(Yr, p = pr, type = "const")
  VAR.rep <- restrict(VAR.rep, method = "man", resmat = matC(m, pr, 1))
  return(VAR.rep)
}

SVAR.sirf.boot <- function(SVAR, Amat, Bmat, Yb, pmax, H, a, R, cumulative = FALSE) {
  Ib <- array(NA, c(SVAR$var$K, SVAR$var$K, H + 1, R))
  nvars <- colnames(SVAR$var$y)
  for (r in 1:R) {
    VAR.rep.R <- boot.rb.VAR(nvars, Yb[, , r], SVAR$var$K, SVAR$var$p, pmax)
    VAR.rep.S <- SVAR(VAR.rep.R, Amat = Amat, Bmat = Bmat, lrtest = FALSE)
    if (cumulative == FALSE) {
      Ib[, , , r] <- SVAR.sirf(VAR.rep.S, H)
    } else {
      Ib[, , , r] <- SVAR.sirf(VAR.rep.S, H, cumulative = TRUE)
    }
  }
  Il <- array(NA, c(SVAR$var$K, SVAR$var$K, H + 1))
  Iu <- array(NA, c(SVAR$var$K, SVAR$var$K, H + 1))
  if (cumulative == FALSE) {
    Ic <- SVAR.sirf(SVAR, H)
  } else {
    Ic <- SVAR.sirf(SVAR, H, cumulative = TRUE)
  }
  for (h in 1:(H + 1)) {
    for (j in 1:m) {
      for (i in 1:m) {
        ci <- boot.ci.efron(Ib[i, j, h, ], a)
        Il[i, j, h] <- ci[1]
        Iu[i, j, h] <- ci[2]
      }
    }
  }
  dimnames(Il) <- dimnames(Ic)
  dimnames(Iu) <- dimnames(Ic)
  SIRF <- list(lb = Il, pe = Ic, ub = Iu)
  return(SIRF)
}

SVAR.fevd.boot <- function(SVAR, Amat, Bmat, Yb, pmax, H, a, R) {
  Fb <- array(NA, c(SVAR$var$K, SVAR$var$K, H + 1, R))
  nvars <- colnames(SVAR$var$y)
  for (r in 1:R) {
    VAR.rep.R <- boot.rb.VAR(nvars, Yb[, , r], SVAR$var$K, SVAR$var$p, pmax)
    VAR.rep.S <- SVAR(VAR.rep.R, Amat = Amat, Bmat = Bmat, lrtest = FALSE)
    Fb[, , , r] <- SVAR.fevd(VAR.rep.S, H)
  }
  Fl <- array(NA, c(SVAR$var$K, SVAR$var$K, H + 1))
  Fu <- array(NA, c(SVAR$var$K, SVAR$var$K, H + 1))
  Fc <- SVAR.fevd(SVAR, H)
  for (h in 1:(H + 1)) {
    for (j in 1:m) {
      for (i in 1:m) {
        ci <- boot.ci.efron(Fb[i, j, h, ], a)
        Fl[i, j, h] <- ci[1]
        Fu[i, j, h] <- ci[2]
      }
    }
  }
  dimnames(Fl) <- dimnames(Fc)
  dimnames(Fu) <- dimnames(Fc)
  FEVD <- list(lb = Fl, pe = Fc, ub = Fu)
  return(FEVD)
}

SVAR.erpt.boot <- function(SVAR, Amat, Bmat, Yb, pmax, H, vx, vy, a, R, cumulative = FALSE) {
  Eb <- array(NA, c(H + 1, R))
  nvars <- colnames(SVAR$var$y)
  for (r in 1:R) {
    VAR.rep.R <- boot.rb.VAR(nvars, Yb[, , r], SVAR$var$K, SVAR$var$p, pmax)
    VAR.rep.S <- SVAR(VAR.rep.R, Amat = Amat, Bmat = Bmat, lrtest = FALSE)
    if (cumulative == FALSE) {
      Eb[, r] <- SVAR.erpt(VAR.rep.S, H, vx, vy)
    } else {
      Eb[, r] <- SVAR.erpt(VAR.rep.S, H, vx, vy, cumulative = TRUE)
    }
  }
  El <- array(NA, c(H + 1))
  Eu <- array(NA, c(H + 1))
  if (cumulative == FALSE) {
    Ec <- SVAR.erpt(SVAR, H, vx, vy)
  } else {
    Ec <- SVAR.erpt(SVAR, H, vx, vy, cumulative = TRUE)
  }
  for (h in 1:(H + 1)) {
    ci <- boot.ci.efron(Eb[h, ], a)
    El[h] <- ci[1]
    Eu[h] <- ci[2]
  }
  ERPT <- list(lb = El, pe = Ec, ub = Eu)
  return(ERPT)
}