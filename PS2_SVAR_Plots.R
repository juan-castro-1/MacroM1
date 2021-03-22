plot.sirf <- function(X, m, H) {
  for (i in 1:m) {
    par(mfrow = c(m, 1))
    for (j in 1:m) {
      plot(0:H, X[i, j, ],
           main = paste("Response of", dimnames(X)[[1]][i], "to", dimnames(X)[[2]][j], "shock", sep = " "),
           xlab = "Horizon", ylab = "",
           type = "o", lwd = 2)
      grid(NULL, NULL, lty = 1)
    }
  }
}

plot.fevd <- function(X, m, H) {
  par(mfrow = c(1, 1))
  for (i in 1:m) {
    barplot(X[i, , ], names.arg = as.character(0:H),
            main = paste("FEVD of ", dimnames(X)[[1]][i], " (log. differences)", sep = ""),
            xlab = "Horizon", ylab = "%",
            legend.text = paste("Shock:", dimnames(X)[[2]], sep = " "))
  }
}

plot.hd <- function(Y, X, m, pmax) {
  for (i in 1:m) {
    par(mfrow = c(m + 1, 1))
    for (j in 1:(m + 1)) {
      plot.ts(ts(X[i, j, ], end = end(Y), frequency = 12),
           main = paste("Influence of", dimnames(X)[[2]][j], "in", dimnames(X)[[1]][i], sep = " "),
           xlab = "Time", ylab = "",
           ylim = c(min(Y[(p + 1):nrow(Y), i] - 1), max(Y[(p + 1):nrow(Y), i] + 1)),
           type = "l", lwd = 2)
      lines(ts(Y[(pmax - p + 1):nrow(Y), i], end = end(Y), frequency = 12), col = "red", lwd = 2)
      abline(h = 0, col = "grey", lwd = 2)
      grid(NULL, NULL, lty = 1)
    }
  }
}

plot.erpt <- function(X, H) {
  par(mfrow = c(1, 1))
  plot(0:H, X,
       main = "Exchange-rate pass-through to consumer prices",
       xlab = "Horizon", ylab = "%",
       ylim = c(0, 150),
       type = "o", lwd = 2)
  grid(NULL, NULL, lty = 1)
}

plot.sirf.boot <- function(X, m, H) {
  for (i in 1:m) {
    par(mfrow = c(m, 1))
    for (j in 1:m) {
      plot(0:H, X$pe[i, j, ],
           main = paste("Response of", dimnames(X$pe)[[1]][i], "to", dimnames(X$pe)[[2]][j], "shock", sep = " "),
           xlab = "Horizon", ylab = "",
           ylim = c(min(X$lb[i, j, ]), max(X$ub[i, j, ])),
           type = "o", lwd = 2)
      grid(NULL, NULL, lty = 1)
      xx <- c(0:H, H:0)
      yy <- c(c(X$lb[i, j, ]), rev(c(X$ub[i, j, ])))
      polygon(xx, yy, col = adjustcolor("grey", alpha.f = 0.5), border = NA)
    }
  }
}

plot.fevd.boot <- function(X, m, H) {
  for (i in 1:m) {
    par(mfrow = c(m, 1))
    for (j in 1:m) {
      plot(0:H, X$pe[i, j, ],
           main = paste("Contribution of", dimnames(X$pe)[[2]][j], "shock to variance of", dimnames(X$pe)[[1]][i], sep = " "),
           xlab = "Horizon", ylab = "%",
           ylim = c(0, 100),
           type = "o", lwd = 2)
      grid(NULL, NULL, lty = 1)
      xx <- c(0:H, H:0)
      yy <- c(c(X$lb[i, j, ]), rev(c(X$ub[i, j, ])))
      polygon(xx, yy, col = adjustcolor("grey", alpha.f = 0.5), border = NA)
    }
  }
}

plot.erpt.boot <- function(X, H) {
  par(mfrow = c(1, 1))
  plot(0:H, X$pe,
       main = "Exchange-rate pass-through to consumer prices",
       xlab = "Horizon", ylab = "%",
       ylim = c(0, 150),
       type = "o", lwd = 2)
  grid(NULL, NULL, lty = 1)
  xx <- c(0:H, H:0)
  yy <- c(c(X$lb), rev(c(X$ub)))
  polygon(xx, yy, col = adjustcolor("grey", alpha.f = 0.5), border = NA)
}