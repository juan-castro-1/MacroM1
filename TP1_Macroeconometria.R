############################################################ MACROECONOMETRIA - TP 1 ############################################################ 

rm(list = ls())

#Para cargar directo el directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#setwd("D:\\Nicolas Der\\Desktop\\UdeSA\\Cursada\\Tercer Trimestre\\Macroeconometría\\TPs\\TP1")
setwd("C:\Users\juan_\Desktop\MacroMetrics\TP_1")
source("PS1_Data.R")
source("PS2_SVAR_Analysis.R")
source("PS2_SVAR_Bootstrap.R")
source("PS2_SVAR_Plots.R")

# 1----
# Data 

Yl.f <- cbind(pcom, er, pc) # Raw data in log
Yl.f <- log(Yl.f) # log transformation
Yd.f <- 100 * diff(Yl.f) # Raw data in log-differences

Yl <- window(Yl.f, start = c(2004, 01), end = c(2019, 12))
Yd <- window(Yd.f, start = c(2004, 01), end = c(2019, 12))

# VAR estimation 

library(vars)

# Lag order selection
pmax <- 12

popt <- VARselect(Yd, lag.max = pmax, type = "const")
popt
p <- popt$selection[1] # AIC

Yd0 <- Yd[1:pmax, ] # Initial values
Ydt <- Yd[(pmax - p + 1):nrow(Yd), ] # Starting in Jan-04

# Estimation (VAR irrestricto... Restringue mas abajo)
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
summary(VAR)

# Model checking (Autovalores en modulo menor a 1 => Variable Estacionaria)
roots(VAR, modulus = TRUE)
serial.test(VAR, lags.bg = 12, type = "PT.asymptotic")

# SVAR estimation 

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

# SVAR analysis

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

# Bootstrap inference

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

# ERPT in log-levels
# Ad-hoc function
assess.erpt <- function(Amat, Bmat, Y, pmax, r, H, R, a, cumulative = FALSE) {
  popt <- VARselect(Y, lag.max = pmax, type = "const")
  if (cumulative == FALSE) {
    p <- popt$selection[3] + 1 # SC + 1, Killian & L?tkepohl (pp. 374)
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

############################################################################################################################################
########################################################## PUNTO 2 #########################################################################
############################################################################################################################################


rm(list=ls())
#$setwd("D:\\Nicolas Der\\Desktop\\UdeSA\\Cursada\\Tercer Trimestre\\Macroeconometría\\TPs\\TP1")
source("PS1_Data.R")
source("PS2_SVAR_Analysis.R")
source("PS2_SVAR_Bootstrap2.R")
source("PS2_SVAR_Plots.R")
# Data
Yl.f <- cbind(er, pc) # Raw data in log
Yl.f <- log(Yl.f) # log transformation
Yd.f <- 100 * diff(Yl.f) # Raw data in log-differences

Yl <- window(Yl.f, start = c(2004, 01), end = c(2019, 12))
Yd <- window(Yd.f, start = c(2004, 01), end = c(2019, 12))

# VAR estimation

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
#VAR <- restrict(VAR, method = "man", resmat = matC(m, p, 1))
summary(VAR)

# Model checking
roots(VAR, modulus = TRUE)
serial.test(VAR, lags.bg = 12, type = "ES")

# SVAR estimation

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

# SVAR analysis
source("PS2_SVAR_Analysis.R")
source("PS2_SVAR_Bootstrap2.R")
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
SVAR.ERPT <- SVAR.erpt(SVAR, H_ERPT, 2, 1, cumulative = TRUE)
plot.erpt(SVAR.ERPT, H_ERPT)

# Bootstrap inference

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
SVAR.ERPT.boot <- SVAR.erpt.boot(SVAR, Amat, Bmat, Yb, pmax, H_ERPT, 2, 1, a, R, cumulative = TRUE)
plot.erpt.boot(SVAR.ERPT.boot, H_ERPT)

# ERPT in log-levels 
# Ad-hoc function
assess.erpt <- function(Amat, Bmat, Y, pmax, r, H, R, a, cumulative = FALSE) {
  popt <- VARselect(Y, lag.max = pmax, type = "const")
  if (cumulative == FALSE) {
    p <- popt$selection[3] + 1 # SC + 1, Killian & L?tkepohl (pp. 374)
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
    SVAR.ERPT.boot <- SVAR.erpt.boot(SVAR, Amat, Bmat, Yb, pmax, H, 2, 1, a, R, cumulative = FALSE)
  } else {
    SVAR.ERPT.boot <- SVAR.erpt.boot(SVAR, Amat, Bmat, Yb, pmax, H, 2, 1, a, R, cumulative = TRUE)
  }
  return(SVAR.ERPT.boot)
}

# ERPT in log-levels
SVAR.ERPT.boot.lvl <- assess.erpt(Amat, Bmat, Yl, pmax, 12, H_ERPT, R, a)
plot.erpt.boot(SVAR.ERPT.boot.lvl, H_ERPT) 


############################################################################################################################################
########################################################## PUNTO 3 #########################################################################
############################################################################################################################################

rm(list = ls())
#setwd("D:\\Nicolas Der\\Desktop\\UdeSA\\Cursada\\Tercer Trimestre\\Macroeconometría\\TPs\\TP1")
source("PS1_Data.R")
source("PS2_SVAR_Analysis.R")
source("PS2_SVAR_Bootstrap.R")
source("PS2_SVAR_Plots.R")

# (a) ene-2005 a jul-2011, ----

# Data

Yl.f <- cbind(pcom, er, pc) # Raw data in log
Yl.f <- log(Yl.f) # log transformation
Yd.f <- 100 * diff(Yl.f) # Raw data in log-differences

Yl <- window(Yl.f, start = c(2004, 01), end = c(2011, 07))
Yd <- window(Yd.f, start = c(2004, 01), end = c(2011, 07))

# VAR estimation 

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

# SVAR estimation 

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

# SVAR analysis

source("practica/PS2//PS2_SVAR_Analysis.R")
source("practica/PS2/PS2_SVAR_Bootstrap.R")
source("practica/PS2/PS2_SVAR_Plots.R")

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

# Bootstrap inference

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

# ERPT in log-levels
# Ad-hoc function
assess.erpt <- function(Amat, Bmat, Y, pmax, r, H, R, a, cumulative = FALSE) {
  popt <- VARselect(Y, lag.max = pmax, type = "const")
  if (cumulative == FALSE) {
    p <- popt$selection[3] + 1 # SC + 1, Killian & L?tkepohl (pp. 374)
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

# (b) ago-2011 a sep-2015 -----

# Data 
source("practica/PS1/PS1_Data.R")
source("practica/PS2//PS2_SVAR_Analysis.R")
source("practica/PS2/PS2_SVAR_Bootstrap.R")
source("practica/PS2/PS2_SVAR_Plots.R")

Yl.f <- cbind(pcom, er, pc) # Raw data in log
Yl.f <- log(Yl.f) # log transformation
Yd.f <- 100 * diff(Yl.f) # Raw data in log-differences

Yl <- window(Yl.f, start = c(2010, 08), end = c(2015, 09))
Yd <- window(Yd.f, start = c(2010, 08), end = c(2015, 09))

# VAR estimation 

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

# SVAR estimation 

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

# SVAR analysis

source("practica/PS2//PS2_SVAR_Analysis.R")
source("practica/PS2/PS2_SVAR_Bootstrap.R")
source("practica/PS2/PS2_SVAR_Plots.R")

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

# Bootstrap inference

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

# ERPT in log-levels
# Ad-hoc function
assess.erpt <- function(Amat, Bmat, Y, pmax, r, H, R, a, cumulative = FALSE) {
  popt <- VARselect(Y, lag.max = pmax, type = "const")
  if (cumulative == FALSE) {
    p <- popt$selection[3] + 1 # SC + 1, Killian & L?tkepohl (pp. 374)
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


#(c) oct-2015 a ago-2019; ----

# Data 
source("practica/PS1/PS1_Data.R")
source("practica/PS2//PS2_SVAR_Analysis.R")
source("practica/PS2/PS2_SVAR_Bootstrap.R")
source("practica/PS2/PS2_SVAR_Plots.R")

Yl.f <- cbind(pcom, er, pc) # Raw data in log
Yl.f <- log(Yl.f) # log transformation
Yd.f <- 100 * diff(Yl.f) # Raw data in log-differences

Yl <- window(Yl.f, start = c(2014, 10), end = c(2019, 08))
Yd <- window(Yd.f, start = c(2014, 10), end = c(2019, 08))

# VAR estimation 

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

# SVAR estimation 

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

# SVAR analysis

source("practica/PS2//PS2_SVAR_Analysis.R")
source("practica/PS2/PS2_SVAR_Bootstrap.R")
source("practica/PS2/PS2_SVAR_Plots.R")

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

# Bootstrap inference

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

# ERPT in log-levels
# Ad-hoc function
assess.erpt <- function(Amat, Bmat, Y, pmax, r, H, R, a, cumulative = FALSE) {
  popt <- VARselect(Y, lag.max = pmax, type = "const")
  if (cumulative == FALSE) {
    p <- popt$selection[3] + 1 # SC + 1, Killian & L?tkepohl (pp. 374)
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


############################################################################################################################################
########################################################## PUNTO 4 #########################################################################
############################################################################################################################################



############################################################################################################################################
########################################################## PUNTO 5 #########################################################################
############################################################################################################################################

setwd("D:\\Nicolas Der\\Desktop\\UdeSA\\Cursada\\Tercer Trimestre\\Macroeconometría\\TPs\\TP1")
source("PS1_Data.R")
source("PS2_SVAR_Analysis.R")
source("PS2_SVAR_Bootstrap.R")
source("PS2_SVAR_Plots.R")

# Data 
Y_aux<-cbind(sal ,pc)
Yl.faux <- log(Y_aux) # log transformation
Yd.faux <- 100 * diff(Yl.faux) # Raw data in log-differences
Ylaux <- window(Yl.faux, start = c(2004, 01), end = c(2019, 12))
Ydaux <- window(Yd.faux, start = c(2004, 01), end = c(2019, 12))
pmax <- 12
p_VAR_aux<- VARselect(Ydaux, lag.max = pmax, type = "const")
p_VAR_aux
p <- p_VAR_aux$selection[1]
Yd0_aux <- Ydaux[1:pmax, ] # Initial values
Ydt_aux <- Ydaux[(pmax - p + 1):nrow(Ydaux), ] # Starting in Jan-04

VAR_aux<- VAR(Ydt_aux,lag.max = pmax, type = "const")
causality(VAR_aux, cause = "sal")
causality(VAR_aux, cause = "pc") # PARECE SER QUE EN UN VAR1 IPC causa en Sentido de Granger a SALARIOS--> PONER PRIMERO IPC EN EL SVAR


Yl.f <- cbind(pcom, er,pc,sal) # Raw data in log
Yl.f <- log(Yl.f) # log transformation
Yd.f <- 100 * diff(Yl.f) # Raw data in log-differences

Yl <- window(Yl.f, start = c(2004, 01), end = c(2019, 12))
Yd <- window(Yd.f, start = c(2004, 01), end = c(2019, 12))

# VAR estimation 

library(vars)

# Lag order selection
pmax <- 12

popt <- VARselect(Yd, lag.max = pmax, type = "const")
popt
p <- popt$selection[1] # AIC

Yd0 <- Yd[1:pmax, ] # Initial values
Ydt <- Yd[(pmax - p + 1):nrow(Yd), ] # Starting in Jan-04

# Estimation (VAR irrestricto... Restringue mas abajo)
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
summary(VAR)

# Model checking (Autovalores en modulo menor a 1 => Variable Estacionaria)
roots(VAR, modulus = TRUE)
serial.test(VAR, lags.bg = 12, type = "PT.asymptotic")

# SVAR estimation 

# A Matrix
Amat <- diag(m)
for (i in 2:m) {
  for (j in 1:(i - 1)) {
    Amat[i, j] <- NA
  }
}

Amat[4,3]<- 0

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

# SVAR analysis

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

# Bootstrap inference

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

# ERPT in log-levels
# Ad-hoc function
assess.erpt <- function(Amat, Bmat, Y, pmax, r, H, R, a, cumulative = FALSE) {
  popt <- VARselect(Y, lag.max = pmax, type = "const")
  if (cumulative == FALSE) {
    p <- popt$selection[3] + 1 # SC + 1, Killian & L?tkepohl (pp. 374)
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

############################################################################################################################################
########################################################## PUNTO 6 #########################################################################
############################################################################################################################################

setwd("D:\\Nicolas Der\\Desktop\\UdeSA\\Cursada\\Tercer Trimestre\\Macroeconometría\\TPs\\TP1")
source("PS1_Data.R")
source("PS2_SVAR_Analysis.R")
source("PS2_SVAR_Bootstrap.R")
source("PS2_SVAR_Plots.R")

# Data 
dlpc<- diff(log(pc))
Y_aux<-cbind(dlog.emae.adj,dlpc)
Ydaux <- window(Y_aux, start = c(2004, 02), end = c(2019, 12))

pmax <- 12
p_VAR_aux<- VARselect(Ydaux, lag.max = pmax, type = "const")
p_VAR_aux
p <- p_VAR_aux$selection[1]
Yd0_aux <- Ydaux[1:pmax, ] # Initial values
Ydt_aux <- Ydaux[(pmax - p + 1):nrow(Ydaux), ] # Starting in Jan-04

VAR_aux<- VAR(Ydt_aux,lag.max = pmax, type = "const")
causality(VAR_aux, cause = "dlog.emae.adj")
causality(VAR_aux, cause = "dlpc") # PARECE SER QUE EN UN VAR3 hay causalidad instantanea entre si


Yl.f <- cbind(pcom, er,pc,sal) # Raw data in log
Yl.f <- log(Yl.f) # log transformation
Yd.f <- 100 * diff(Yl.f) # Raw data in log-differences

Yl <- window(Yl.f, start = c(2004, 01), end = c(2019, 12))
Yd <- window(Yd.f, start = c(2004, 01), end = c(2019, 12))

# VAR estimation 

library(vars)

# Lag order selection
pmax <- 12

popt <- VARselect(Yd, lag.max = pmax, type = "const")
popt
p <- popt$selection[1] # AIC

Yd0 <- Yd[1:pmax, ] # Initial values
Ydt <- Yd[(pmax - p + 1):nrow(Yd), ] # Starting in Jan-04

# Estimation (VAR irrestricto... Restringue mas abajo)
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
summary(VAR)

# Model checking (Autovalores en modulo menor a 1 => Variable Estacionaria)
roots(VAR, modulus = TRUE)
serial.test(VAR, lags.bg = 12, type = "PT.asymptotic")

# SVAR estimation 

# A Matrix
Amat <- diag(m)
for (i in 2:m) {
  for (j in 1:(i - 1)) {
    Amat[i, j] <- NA
  }
}

Amat[4,3]<- 0

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

# SVAR analysis

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

# Bootstrap inference

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

# ERPT in log-levels
# Ad-hoc function
assess.erpt <- function(Amat, Bmat, Y, pmax, r, H, R, a, cumulative = FALSE) {
  popt <- VARselect(Y, lag.max = pmax, type = "const")
  if (cumulative == FALSE) {
    p <- popt$selection[3] + 1 # SC + 1, Killian & L?tkepohl (pp. 374)
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


############################################################################################################################################
########################################################## PUNTO 7 #########################################################################
############################################################################################################################################

setwd("D:\\Nicolas Der\\Desktop\\UdeSA\\Cursada\\Tercer Trimestre\\Macroeconometría\\TPs\\TP1")
source("PS1_Data.R")
source("PS2_SVAR_Analysis.R")
source("PS2_SVAR_Bootstrap.R")
source("PS2_SVAR_Plots.R")

# Data 
Y_aux<-cbind(m2 ,pc)
Yl.faux <- log(Y_aux) # log transformation
Yd.faux <- 100 * diff(Yl.faux) # Raw data in log-differences
Ylaux <- window(Yl.faux, start = c(2004, 01), end = c(2019, 12))
Ydaux <- window(Yd.faux, start = c(2004, 01), end = c(2019, 12))
pmax <- 12
p_VAR_aux<- VARselect(Ydaux, lag.max = pmax, type = "const")
p_VAR_aux
p <- p_VAR_aux$selection[1]
Yd0_aux <- Ydaux[1:pmax, ] # Initial values
Ydt_aux <- Ydaux[(pmax - p + 1):nrow(Ydaux), ] # Starting in Jan-04

VAR_aux<- VAR(Ydt_aux,lag.max = p, type = "const")
causality(VAR_aux, cause = "m2")
causality(VAR_aux, cause = "pc") # PARECE SER QUE EN UN VAR3 IPC causa en Sentido de Granger a M2--> PONER PRIMERO IPC EN EL SVAR


Yl.f <- cbind(pcom, er,pc,m2) # Raw data in log
Yl.f <- log(Yl.f) # log transformation
Yd.f <- 100 * diff(Yl.f) # Raw data in log-differences

Yl <- window(Yl.f, start = c(2004, 01), end = c(2019, 12))
Yd <- window(Yd.f, start = c(2004, 01), end = c(2019, 12))

# VAR estimation 

library(vars)

# Lag order selection
pmax <- 12

popt <- VARselect(Yd, lag.max = pmax, type = "const")
popt
p <- popt$selection[1] # AIC

Yd0 <- Yd[1:pmax, ] # Initial values
Ydt <- Yd[(pmax - p + 1):nrow(Yd), ] # Starting in Jan-04

# Estimation (VAR irrestricto... Restringue mas abajo)
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
summary(VAR)

# Model checking (Autovalores en modulo menor a 1 => Variable Estacionaria)
roots(VAR, modulus = TRUE)
serial.test(VAR, lags.bg = 12, type = "ES")

# SVAR estimation 

# A Matrix
Amat <- diag(m)
for (i in 2:m) {
  for (j in 1:(i - 1)) {
    Amat[i, j] <- NA
  }
}

Amat[4,3]<-0

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

# SVAR analysis

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

# Bootstrap inference

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

# ERPT in log-levels
# Ad-hoc function
assess.erpt <- function(Amat, Bmat, Y, pmax, r, H, R, a, cumulative = FALSE) {
  popt <- VARselect(Y, lag.max = pmax, type = "const")
  if (cumulative == FALSE) {
    p <- popt$selection[3] + 1 # SC + 1, Killian & L?tkepohl (pp. 374)
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


############################################################################################################################################
########################################################## PUNTO 8 #########################################################################
############################################################################################################################################

#setwd("D:\\Nicolas Der\\Desktop\\UdeSA\\Cursada\\Tercer Trimestre\\Macroeconometría\\TPs\\TP1")
source("PS1_Data.R")
source("PS2_SVAR_Analysis.R")
source("PS2_SVAR_Bootstrap.R")
source("PS2_SVAR_Plots.R")

# Data 
Y_aux<-cbind(td ,pc)
Yl.faux <- log(Y_aux) # log transformation
Yd.faux <- 100 * diff(Yl.faux) # Raw data in log-differences
Ylaux <- window(Yl.faux, start = c(2004, 01), end = c(2019, 12))
Ydaux <- window(Yd.faux, start = c(2004, 01), end = c(2019, 12))
pmax <- 12
p_VAR_aux<- VARselect(Ydaux, lag.max = pmax, type = "const")
p_VAR_aux
p <- p_VAR_aux$selection[2]
Yd0_aux <- Ydaux[1:pmax, ] # Initial values
Ydt_aux <- Ydaux[(pmax - p + 1):nrow(Ydaux), ] # Starting in Jan-04

VAR_aux<- VAR(Ydt_aux,lag.max = 3, type = "const")
summary(VAR_aux)
serial.test(VAR_aux,type="ES")
causality(VAR_aux, cause = "td")
causality(VAR_aux, cause = "pc") # PARECE SER QUE EN UN VAR2 no hayI causalidAD en Sentido de Granger ENTRE IPC Y TASA


Yl.f <- cbind(pcom, er,td,pc) # Raw data in log
Yl.f <- log(Yl.f) # log transformation
Yd.f <- 100 * diff(Yl.f) # Raw data in log-differences

Yl <- window(Yl.f, start = c(2004, 01), end = c(2019, 12))
Yd <- window(Yd.f, start = c(2004, 01), end = c(2019, 12))

# VAR estimation 

library(vars)

# Lag order selection
pmax <- 12

popt <- VARselect(Yd, lag.max = pmax, type = "const")
popt
p <- popt$selection[2] # AIC

Yd0 <- Yd[1:pmax, ] # Initial values
Ydt <- Yd[(pmax - p + 1):nrow(Yd), ] # Starting in Jan-04

# Estimation (VAR irrestricto... Restringue mas abajo)
VAR <- VAR(Ydt, p = 4, type = "const")

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
VAR <- restrict(VAR, method = "man", resmat = matC(m, p, 1)) #ACA ME ESTA TRAYENDO ERROR EN LA MATc!!!!!!!!!!!!!!!
summary(VAR)# CON EL ERROR DE LA MATC, se CHOTEA LAS IRF Y PCOM TIENE EFECTO DE LAS VARIABLES CUANDO EN REALIDAD DEBE SER UNA RESTRICCION

# Model checking (Autovalores en modulo menor a 1 => Variable Estacionaria)
roots(VAR, modulus = TRUE)
serial.test(VAR, lags.bg = 12, type = "ES")

# SVAR estimation 

# A Matrix
Amat <- diag(m)
for (i in 2:m) {
  for (j in 1:(i - 1)) {
    Amat[i, j] <- NA
  }
}

Amat[4,3]<-0

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

# SVAR analysis

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

# Bootstrap inference

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

# ERPT in log-levels
# Ad-hoc function
assess.erpt <- function(Amat, Bmat, Y, pmax, r, H, R, a, cumulative = FALSE) {
  popt <- VARselect(Y, lag.max = pmax, type = "const")
  if (cumulative == FALSE) {
    p <- popt$selection[3] + 1 # SC + 1, Killian & L?tkepohl (pp. 374)
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

############################################################################################################################################
########################################################## PUNTO 9 #########################################################################
############################################################################################################################################

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("PS1_Data.R")
source("PS2_SVAR_Analysis.R")
source("PS2_SVAR_Bootstrap.R")
source("PS2_SVAR_Plots.R")


# Data pcom, er, pc
Y_aux<-cbind(sal,m2,td,pcom,er,pc)
Y_aux<-cbind(pcom,er,td,pc)

Yl.faux <- log(Y_aux) # log transformation
Yd.faux <- 100 * diff(Yl.faux) # Raw data in log-differences
Ylaux <- window(Yl.faux, start = c(2004, 01), end = c(2019, 12))
Ydaux <- window(Yd.faux, start = c(2004, 01), end = c(2019, 12))
pmax <- 12
p_VAR_aux<- VARselect(Ydaux, lag.max = pmax, type = "const")
p_VAR_aux
p <- p_VAR_aux$selection[2]
Yd0_aux <- Ydaux[1:pmax, ] # Initial values
Ydt_aux <- Ydaux[(pmax - p + 1):nrow(Ydaux), ] # Starting in Jan-04

VAR_aux<- VAR(Ydt_aux,lag.max = 4, type = "const")
summary(VAR_aux)
serial.test(VAR_aux,type="ES")
causality(VAR_aux, cause = "td")
causality(VAR_aux, cause = "pc") # PARECE SER QUE EN UN VAR2 no hayI causalidAD en Sentido de Granger ENTRE IPC Y TASA


Yl.f <- cbind(pcom,sal,m2,td,er,pc) # Raw data in log
Yl.f <- na.omit(Yl.f)
Yl.f <- log(Yl.f) # log transformation
Yd.f <- 100 * diff(Yl.f) # Raw data in log-differences

Yl.f<-cbind(Yl.f,dlog.emae)
Yd.f<-cbind(Yd.f,dlog.emae.adj)

Yl.f <- Yl.f[,c(1,5,7,2,3,4,6)]
Yd.f <- Yd.f[,c(1,5,7,2,3,4,6)]

Yl <- window(Yl.f, start = c(2004, 01), end = c(2019, 12))
Yd <- window(Yd.f, start = c(2004, 02), end = c(2019, 12))

# VAR estimation 

library(vars)

# Lag order selection
pmax <- 12

popt <- VARselect(Yd, lag.max = pmax, type = "const")
popt
p <- popt$selection[1] # AIC

Yd0 <- Yd[1:pmax, ] # Initial values
Ydt <- Yd[(pmax - p + 1):nrow(Yd), ] # Starting in Jan-04

# Estimation (VAR irrestricto... Restringue mas abajo)
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
VAR <- restrict(VAR, method = "man", resmat = matC(m, p, 1)) #ACA ME ESTA TRAYENDO ERROR EN LA MATc!!!!!!!!!!!!!!!
summary(VAR)# CON EL ERROR DE LA MATC, se CHOTEA LAS IRF Y PCOM TIENE EFECTO DE LAS VARIABLES CUANDO EN REALIDAD DEBE SER UNA RESTRICCION

# Model checking (Autovalores en modulo menor a 1 => Variable Estacionaria)
roots(VAR, modulus = TRUE)
serial.test(VAR, lags.bg = 12, type = "ES")

# SVAR estimation 

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

# SVAR analysis

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

# Bootstrap inference

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

# ERPT in log-levels
# Ad-hoc function
assess.erpt <- function(Amat, Bmat, Y, pmax, r, H, R, a, cumulative = FALSE) {
  popt <- VARselect(Y, lag.max = pmax, type = "const")
  if (cumulative == FALSE) {
    p <- popt$selection[3] + 1 # SC + 1, Killian & L?tkepohl (pp. 374)
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




############################################################################################################################################
########################################################## PUNTO 12 #########################################################################
############################################################################################################################################


#############################
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(readxl)
ccl<-read_excel("ccl.xlsx",sheet = "ccl")
ccl<-ccl[, -c(1:2)]
ccl<-ts(ccl, frequency=12, start = c(2003,1))
ccl[1:96]<-0
ccl[157:204]<-0
brecha_ccl <- ccl        #creo la brecha
#############################

source("PS1_Data.R")
source("PS2_SVAR_Analysis.R")
source("PS2_SVAR_Bootstrap.R")
source("PS2_SVAR_Plots.R")



# 1----
# Data 

Yl.f <- cbind(pcom, er, brecha_ccl, pc) # Raw data in log
Yl.f[is.na(brecha_ccl)]<-0
Yl.f <- log(Yl.f) # log transformation
Yd.f <- 100 * diff(Yl.f) # Raw data in log-differences
Yd.f[is.na(brecha_ccl)]<-0
is.na(Yd.f)<-0

Yl <- window(Yl.f, start = c(2004, 01), end = c(2019, 12))
Yd <- window(Yd.f, start = c(2004, 01), end = c(2019, 12))
Yd.f <- window(Yd.f, start = c(2004, 01), end = c(2019, 12))


# VAR estimation 

library(vars)

# Lag order selection
pmax <- 12

popt <- VARselect(Yd.f, lag.max = pmax, type = "const")
popt
p <- popt$selection[1] # AIC

Yd0 <- Yd[1:pmax, ] # Initial values
Ydt <- Yd[(pmax - p + 1):nrow(Yd), ] # Starting in Jan-04

# Estimation (VAR irrestricto... Restringue mas abajo)
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
summary(VAR)

# Model checking (Autovalores en modulo menor a 1 => Variable Estacionaria)
roots(VAR, modulus = TRUE)
serial.test(VAR, lags.bg = 12, type = "PT.asymptotic")

# SVAR estimation 

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

# SVAR analysis

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

# Bootstrap inference

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

# ERPT in log-levels
# Ad-hoc function
assess.erpt <- function(Amat, Bmat, Y, pmax, r, H, R, a, cumulative = FALSE) {
  popt <- VARselect(Y, lag.max = pmax, type = "const")
  if (cumulative == FALSE) {
    p <- popt$selection[3] + 1 # SC + 1, Killian & L?tkepohl (pp. 374)
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






