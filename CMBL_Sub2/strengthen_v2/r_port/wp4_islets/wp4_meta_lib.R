## wp4_meta_lib.R - the pooling functions of wp4_islets/meta_analysis.py (protocol 4.4), as a library for the scripts that pool (wp4_meta_analysis.R, positive_control.R, wider_search_analyze.R). Source pycompat.R first (np_sum).
## Estimators (fixed in advance): tau-squared by REML; Hartung-Knapp standard error with t(k-1) quantile and no floor = PRIMARY interval; sensitivity: Knapp-Hartung factor floored at 1, DerSimonian-Laird tau-squared with the
## Hartung-Knapp interval, DerSimonian-Laird with the normal interval, fixed effect. Sums follow numpy (np_sum: sequential below 8 values), so that the iteration of the REML estimate gives the same digits as the Python run.
MARGIN <- 0.10
ssum <- np_sum
reml_tau2 <- function(y, v, tol = 1e-12, maxit = 500) {
  n <- length(y); tau2 <- if (n > 1) max(0.0, (ssum((y - ssum(y) / n)^2) / (n - 1)) - ssum(v) / n) else 0.0
  for (i in seq_len(maxit)) {
    w <- 1 / (v + tau2); mu <- ssum(w * y) / ssum(w); new <- max(0.0, ssum(w^2 * ((y - mu)^2 - v)) / ssum(w^2) + 1 / ssum(w)); if (abs(new - tau2) < tol) return(new); tau2 <- new
  }
  tau2
}
dl_tau2 <- function(y, v) { w <- 1 / v; mu <- ssum(w * y) / ssum(w); Q <- ssum(w * (y - mu)^2); k <- length(y); max(0.0, (Q - (k - 1)) / (ssum(w) - ssum(w^2) / ssum(w))) }
pool <- function(y, se, tau = "reml", ci = "hk", floor = FALSE) {      # list(mu, se, lo, hi, tau2, I2, Q, k) in the key order of the Python dict
  y <- as.numeric(y); v <- as.numeric(se)^2; k <- length(y); w0 <- 1 / v; mu0 <- ssum(w0 * y) / ssum(w0); Q <- ssum(w0 * (y - mu0)^2); I2 <- if (Q > 0) max(0.0, (Q - (k - 1)) / Q) else 0.0
  if (ci == "fixed") return(list(mu = mu0, se = sqrt(1 / ssum(w0)), lo = mu0 - 1.959964 / sqrt(ssum(w0)), hi = mu0 + 1.959964 / sqrt(ssum(w0)), tau2 = 0.0, I2 = I2, Q = Q, k = as.integer(k)))
  tau2 <- if (tau == "reml") reml_tau2(y, v) else dl_tau2(y, v); w <- 1 / (v + tau2); mu <- ssum(w * y) / ssum(w)
  if (ci == "normal") { se_ <- sqrt(1 / ssum(w)); q <- 1.959964 } else { fac <- ssum(w * (y - mu)^2) / (k - 1); if (floor) fac <- max(1.0, fac); se_ <- sqrt(fac / ssum(w)); q <- qt(0.975, k - 1) }
  list(mu = mu, se = se_, lo = mu - q * se_, hi = mu + q * se_, tau2 = tau2, I2 = I2, Q = Q, k = as.integer(k))
}
verdict <- function(lo, hi, m = MARGIN) if (lo > m) "consensus pairs are co-expressed more than random pairs" else if (lo >= -m && hi <= m) "no meaningful excess co-expression" else "inconclusive"
