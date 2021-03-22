remove(list = ls(all.names = TRUE))
gc()

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Data ####

source("PS1_Data.R") # corre el archivo R sin abrirlo

Yl <- cbind(pcom, er, pc) # Raw data in log
Yl <- log(Yl) # log transformation
Yl <- window(Yl, start = c(2004, 01), end = c(2019, 12)) # Keep specific observations
Yd <- 100 * diff(Yl) # Raw data in log-differences

plot(Yl) # Plot the variables (log-levels)
plot(Yd) # Plot the variables (log-differences)

# VAR estimation ####

library(vars)

# Lag order selection
pmax <- 6 # Maximum lag length

popt <- VARselect(Yd, lag.max = pmax, type = "const")
popt
p <- popt$selection[2] # HQIC

# Estimation
VAR <- VAR(Yd, p = p, type = "const") # Inclusion of exogenous variables is also possible
summary(VAR)

plot(VAR)

# Manual plotting of residuals
plot(resid(VAR)[, 1], type = "l", ylab = "Resid. pcom")
plot(resid(VAR)[, 2], type = "l", ylab = "Resid. er")
plot(resid(VAR)[, 3], type = "l", ylab = "Resid. pc")

m <- VAR$K # Number of variables in the VAR
N <- VAR$obs # Number of effective sample observations, excluding "p" starting values

# Granger causality ####

# Manual computation of Wald Test via SUR model

library(systemfit)

# Data and regressors
y_pcom <-  as.numeric(Yd[, "pcom"])
y_er <- as.numeric(Yd[, "er"])
y_pc <- as.numeric(Yd[, "pc"])

y_pcom_l1 <- y_pcom[1:N]
y_er_l1 <- y_er[1:N]
y_pc_l1 <- y_pc[1:N]

y_pcom <-  y_pcom[(p + 1):(N + p)]
y_er <- y_er[(p + 1):(N + p)]
y_pc <- y_pc[(p + 1):(N + p)]

# SUR equations
eqpcom <- y_pcom ~ y_pcom_l1 + y_er_l1 + y_pc_l1
eqer <- y_er ~ y_pcom_l1 + y_er_l1 + y_pc_l1
eqpc <- y_pc ~ y_pcom_l1 + y_er_l1 + y_pc_l1

# SUR Estimation
system <- list(pcomreg = eqpcom, erreg = eqer, pcreg = eqpc)
fitsur <- systemfit(system, method = "OLS")

# Constraints for Granger Causality Test
Rmat <- matrix(0, nrow = 1, ncol = 12)
Rmat[1, 11] <- 1
qvec <- c(0)

# Granger Causality Test, ER -> PC (Manual)
linearHypothesis(fitsur, Rmat, qvec, test = "F" )

# Granger causality test, Local Vars. -> PCOM (asymptotic) 
VAR.test.gc.as <- causality(VAR, cause = c("er", "pc"))
VAR.test.gc.as

# Granger causality test, Local Vars. -> PCOM (bootstrap)
VAR.test.gc.boot <- causality(VAR, cause = c("er", "pc"), boot = TRUE, boot.runs = 2000)
VAR.test.gc.boot

# Model checking ####

q <- min(10, trunc(N / 5)) # Rule of thumb for Portmanteau tests (Rob Hyndman) # https://robjhyndman.com/hyndsight/ljung-box-test/
r <- q # Lag order for BG test auxiliary regression
s <- q # Lag order for ARCH test auxiliary regression

# VAR stability

# Eigenvalues
VAR.roots <- roots(VAR, modulus = TRUE)
VAR.roots

# Residual serial correlation

# Portmanteau test
VAR.test.serial.pt <- serial.test(VAR, lags.pt = q, type = "PT.asymptotic")
VAR.test.serial.pt

# Adjusted Portmanteau test
VAR.test.serial.pt.ss <- serial.test(VAR, lags.pt = q, type = "PT.adjusted") # Small sample correc.
VAR.test.serial.pt.ss

# Breusch-Godfrey test
VAR.test.serial.bg <- serial.test(VAR, lags.bg = r, type = "BG")
VAR.test.serial.bg

# Adjusted Breusch-Godfrey test (Edgerton & Shukur)
VAR.test.serial.bg.ss <- serial.test(VAR, lags.bg = r, type = "ES") # Small sample correc.
VAR.test.serial.bg.ss

# Residual normality

# Multivariate Jarque-Bera test
VAR.test.norm <- normality.test(VAR, multivariate.only = FALSE)
VAR.test.norm

# Residual heteroskedasticity

# Multivariate ARCH test
VAR.test.arch <- arch.test(VAR, lags.multi = s, multivariate.only = FALSE)
VAR.test.arch

# Ad hoc function (checks everything at once for different lag orders)
fast.check <- function(Y, p, pmax, N, q, r, s, X = NULL) {
  CHK <- data.frame(matrix(NA, pmax, 8))
  CHK[, 1] <- 1:pmax
  if (is.null(X)) {
    for (l in 1:pmax) {
      # Yt <- Y
      Yt <- Y[(pmax - l + 1):N + p, ]
      VARt <- VAR(Yt, p = l, type = "const")
      PT.a <- serial.test(VARt, lags.pt = q, type = "PT.adjusted")
      BG.a <- serial.test(VARt, lags.bg = r, type = "ES")
      JB <- normality.test(VARt, multivariate.only = TRUE)
      ARCH <- arch.test(VARt, lags.multi = s, multivariate.only = TRUE)
      CHK[l, 2] <- max(roots(VARt, modulus = TRUE))
      CHK[l, 3] <- PT.a$serial[3]
      CHK[l, 4] <- BG.a$serial[3]
      CHK[l, 5] <- JB$jb.mul[1]$JB[3]
      CHK[l, 6] <- JB$jb.mul[2]$Skewness[3]
      CHK[l, 7] <- JB$jb.mul[3]$Kurtosis[3]
      CHK[l, 8] <- ARCH$arch.mul[3]
    }
  } else {
    for (l in 1:pmax) {
      Yt <- Y
      # Yt <- Y[(pmax - l + 1):N + p, ]
      VARt <- VAR(Yt, p = l, type = "const", exogen = X)
      PT.a <- serial.test(VARt, lags.pt = q, type = "PT.adjusted")
      BG.a <- serial.test(VARt, lags.bg = r, type = "ES")
      JB <- normality.test(VARt, multivariate.only = TRUE)
      ARCH <- arch.test(VARt, lags.multi = s, multivariate.only = TRUE)
      CHK[l, 2] <- max(roots(VARt, modulus = TRUE))
      CHK[l, 3] <- PT.a$serial[3]
      CHK[l, 4] <- BG.a$serial[3]
      CHK[l, 5] <- JB$jb.mul[1]$JB[3]
      CHK[l, 6] <- JB$jb.mul[2]$Skewness[3]
      CHK[l, 7] <- JB$jb.mul[3]$Kurtosis[3]
      CHK[l, 8] <- ARCH$arch.mul[3]
    }
  }
  CHK <- round(CHK, 3)
  CHK <- data.frame(CHK)
  colnames(CHK) <- c("Lag", "MEV", "PT.adj", "BG.adj", "JB", "JB-S", "JB-K", "ARCH")
  return(CHK)
}

fast.check(Yd, p, pmax, N, q, r, s)

# Other ####

# Stability analysis (possible structural change in parameters)
stab.test <- stability(VAR, type = "fluctuation")
plot(stab.test)

# Formal test
sctest(stab.test$stability$pcom) # https://cran.r-project.org/web/packages/strucchange/

# Model simplification ####

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

# Estimate VAR with zero constraints
VAR <- restrict(VAR, method = "man", resmat = matC(m, p, 1)) # GLS estimation
VAR

# Forecasting ####

a <- 0.95 # confidence level
H <- 12 # forecast horizon

# VAR forecasts
VAR.fcst <- predict(VAR, n.head = H, ci = a) # ***
VAR.fcst
# *** Normality assumption is relevant for the confidence bands
# se deberia hacer con bootst

# Plot all forecasts
plot(VAR.fcst)
fanchart(VAR.fcst)

# Plot only inflation
plot(VAR.fcst, names = "pc")
fanchart(VAR.fcst, names = "pc")

# Unit root testing ####

# library(urca) # "vars" is dependent on this

plot(Yl)

# Augmented-Dickey-Fuller Unit Root Test
adf <- ur.df(y = Yl[, "pcom"], type = "drift", lags = 2, selectlags = "Fixed") # ***
summary(adf)
# https://stats.stackexchange.com/questions/24072/interpreting-rs-ur-df-dickey-fuller-unit-root-test-results

# Elliott, Rothenberg & Stock Unit Root Test
# ers <- ur.ers(y = Yl[, "pcom"], type = "DF-GLS", model = "const", lag.max = 4)
# sumary(ers)

# Kwiatkowski-Phillips-Schmidt-Shin Unit Root Test
# kpss <- ur.kpss(y = Yl[, "pcom"], type = "tau", lags = "short")
# summary(kpss)

# Phillips & Perron Unit Root Test
# pp <- ur.pp(y = Yl[, "pcom"], type = "Z-tau", model = "trend", lags = "short")
# summary(pp)

# Schmidt \& Phillips Unit Root Test
# sp <- ur.sp(y = Yl[, "pcom"], type = "tau", pol.deg = 1, signif = 0.01)
# summary(sp)

# Zivot & Andrews Unit Root Test
# za <- ur.za(y = Yl[, "pcom"], model = "both", lag = 2)
# summary(za)