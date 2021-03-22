remove(list = ls(all.names = TRUE))
gc()

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Data ####

source("PS1_Data.R")

Yl.f <- cbind(pcom, er, pc) # Raw data in log
Yl.f <- log(Yl.f) # log transformation
Yd.f <- 100 * diff(Yl.f) # Raw data in log-differences

Yl <- window(Yl.f, start = c(2003, 01), end = c(2019, 12))
Yd <- window(Yd.f, start = c(2003, 01), end = c(2019, 12))

# VAR estimation ####

library(vars)

# Lag order selection
pmax <- 12

popt <- VARselect(Yd, lag.max = pmax, type = "const")
popt
p <- popt$selection[1] # AIC

Yd0 <- Yd[1:pmax, ] # Initial values
Ydt <- Yd[(pmax - p + 1):nrow(Yd), ] # Starting in Jan-04

# Estimation
VAR <- VAR(Ydt, p = p, type = "const")

m <- VAR$K # No. of variables in the VAR
N <- VAR$obs # No. of effective sample observations, excluding "p" starting values

# Ad hoc function
matC <- function(m, p, vx) {
  vy <- setdiff(1:m, vx)
  Cm <- matrix(1, m, m * p + 1)
  for (i in vx) {
    for (l in 1:p) {
      for (j in vy) {
        Cm[i, m * (l - 1) + j] <- 0
      }
    }
  }
  return(Cm)
}

# Re-estimate VAR (no feedback from local vars. to pcom)
VAR <- restrict(VAR, method = "man", resmat = matC(m, p, 1))
VAR

# Model checking
roots(VAR, modulus = TRUE)
serial.test(VAR, lags.bg = 12, type = "ES")

# SVAR estimation ####

# A Matrix
Amat <- diag(m)
for (i in 2:m) {
  for (j in 1:(i - 1)) {
    Amat[i, j] <- NA
  }
}

# B Matrix
Bmat <- matrix(0, m, m)
for (i in 1:m) {
  Bmat[i, i] <- NA
}

# SVAR estimation (AB model configuration)
SVAR <- SVAR(VAR, Amat = Amat, Bmat = Bmat, lrtest = FALSE)
SVAR

# SVAR t0 impact matrix (Cholesky decomposition)
S <- t(resid(VAR)) %*% resid(VAR) / N
P.chol <- t(chol(S)) # Cholesky decomposition
S

# SVAR t0 impact matrix (implied by AB model)
P <- solve(SVAR$A, SVAR$B) # inv(A) %% B
S.SVAR <- P %*% t(P)
S.SVAR

# Other SVAR parameters
pars.R <- Bcoef(VAR) # VAR
pars.S <- solve(P, pars.R) # SVAR
pars.R
pars.S

# SVAR analysis ####

source("PS2_SVAR_Analysis.R")
source("PS2_SVAR_Bootstrap.R")
source("PS2_SVAR_Plots.R")

H <- 18 # Horizon
H_ERPT <- 120 # Horizon for ERPT

# IRF
SVAR.SIRF <- SVAR.sirf(SVAR, H)
plot.sirf(SVAR.SIRF, m, H)

# Cumulative IRF
SVAR.SIRF.c <- SVAR.sirf(SVAR, H, cumulative = TRUE)
plot.sirf(SVAR.SIRF.c, m, H)

# FEVD
SVAR.FEVD <- SVAR.fevd(SVAR, H)
plot.fevd(SVAR.FEVD, m, H)

# HD
SVAR.HD <- SVAR.hd(SVAR)
plot.hd(Yd, SVAR.HD, m, pmax)

# ERPT
SVAR.ERPT <- SVAR.erpt(SVAR, H_ERPT, 3, 2, cumulative = TRUE)
plot.erpt(SVAR.ERPT, H_ERPT)

# Bootstrap inference ####

a <- 0.95 # Confidence level
R <- 500 # No. of bootstrap replications

# Bootstrap replications
Yb <- boot.rb.replicate(VAR, Yd0, pmax, R)

# IRF (bootstrap)
SVAR.SIRF.boot <- SVAR.sirf.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R)
plot.sirf.boot(SVAR.SIRF.boot, m, H)

# Cumulative IRF (bootstrap)
SVAR.SIRF.c.boot <- SVAR.sirf.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R, cumulative = TRUE)
plot.sirf.boot(SVAR.SIRF.c.boot, m, H)

# FEVD (bootstrap)
SVAR.FEVD.boot <- SVAR.fevd.boot(SVAR, Amat, Bmat, Yb, pmax, H, a, R)
plot.fevd.boot(SVAR.FEVD.boot, m, H)

# ERPT (bootstrap)
SVAR.ERPT.boot <- SVAR.erpt.boot(SVAR, Amat, Bmat, Yb, pmax, H_ERPT, 3, 2, a, R, cumulative = TRUE)
plot.erpt.boot(SVAR.ERPT.boot, H_ERPT)

# ERPT in log-levels ####

# Ad-hoc function
assess.erpt <- function(Amat, Bmat, Y, pmax, r, H, R, a, cumulative = FALSE) {
  popt <- VARselect(Y, lag.max = pmax, type = "const")
  if (cumulative == FALSE) {
    p <- popt$selection[3] + 1 # SC + 1, Killian & Lütkepohl (pp. 374)
  } else {
    p <- popt$selection[3] # SC
  }
  Y0 <- Y[1:pmax, ]
  Yt <- Y[(pmax - p + 1):nrow(Y), ]
  VAR <- VAR(Yt, p = p, type = "const")
  VAR <- restrict(VAR, method = "man", resmat = matC(m, p, 1))
  print(roots(VAR, modulus = TRUE))
  print(serial.test(VAR, lags.bg = r, type = "ES"))
  SVAR <- SVAR(VAR, Amat = Amat, Bmat = Bmat, lrtest = FALSE)
  Yb <- boot.rb.replicate(VAR, Y0, pmax, R)
  if (cumulative == FALSE) {
    SVAR.ERPT.boot <- SVAR.erpt.boot(SVAR, Amat, Bmat, Yb, pmax, H, 3, 2, a, R, cumulative = FALSE)
  } else {
    SVAR.ERPT.boot <- SVAR.erpt.boot(SVAR, Amat, Bmat, Yb, pmax, H, 3, 2, a, R, cumulative = TRUE)
  }
  return(SVAR.ERPT.boot)
}

# ERPT in log-levels
SVAR.ERPT.boot.lvl <- assess.erpt(Amat, Bmat, Yl, pmax, 12, H_ERPT, R, a)
plot.erpt.boot(SVAR.ERPT.boot.lvl, H_ERPT)