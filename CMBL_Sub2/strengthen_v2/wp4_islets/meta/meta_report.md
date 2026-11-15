# WP4 meta-analysis (protocol 4.4)

Primary = depth-adjusted E, all frozen beta-cells, REML tau-squared + Hartung-Knapp 95% CI (t, k - 1 df), margin +/-0.10 on the E (Fisher z) scale.

## primary, adjusted mode: 5 qualifying dataset(s): gse81608, gse84133_human, gse86469, gse81547, emtab5061
* REML+HK (primary)              pooled E +0.013  95% CI [-0.052, +0.078]  tau2 0.0011  I2 71%  -> no meaningful excess co-expression
* REML+HK, factor floored at 1   pooled E +0.013  95% CI [-0.052, +0.078]  tau2 0.0011  I2 71%  -> no meaningful excess co-expression
* DL+HK                          pooled E +0.014  95% CI [-0.053, +0.080]  tau2 0.0012  I2 71%  -> no meaningful excess co-expression
* DL, normal CI                  pooled E +0.014  95% CI [-0.025, +0.053]  tau2 0.0012  I2 71%  -> no meaningful excess co-expression
* fixed effect                   pooled E +0.003  95% CI [-0.016, +0.022]  tau2 0.0000  I2 71%  -> no meaningful excess co-expression
  leave out gse81608        pooled E +0.026  [-0.062, +0.115] -> inconclusive
  leave out gse84133_human  pooled E +0.005  [-0.049, +0.059] -> no meaningful excess co-expression
  leave out gse86469        pooled E +0.026  [-0.071, +0.123] -> inconclusive
  leave out gse81547        pooled E +0.015  [-0.095, +0.125] -> inconclusive
  leave out emtab5061       pooled E +0.000  [-0.082, +0.083] -> no meaningful excess co-expression

## primary, unadjusted mode: 5 qualifying dataset(s): gse81608, gse84133_human, gse86469, gse81547, emtab5061
* REML+HK (primary)              pooled E -0.001  95% CI [-0.065, +0.063]  tau2 0.0010  I2 68%  -> no meaningful excess co-expression
* REML+HK, factor floored at 1   pooled E -0.001  95% CI [-0.065, +0.063]  tau2 0.0010  I2 68%  -> no meaningful excess co-expression
* DL+HK                          pooled E +0.000  95% CI [-0.068, +0.068]  tau2 0.0013  I2 68%  -> no meaningful excess co-expression
* DL, normal CI                  pooled E +0.000  95% CI [-0.041, +0.041]  tau2 0.0013  I2 68%  -> no meaningful excess co-expression
* fixed effect                   pooled E -0.009  95% CI [-0.030, +0.011]  tau2 0.0000  I2 68%  -> no meaningful excess co-expression
  leave out gse81608        pooled E +0.013  [-0.065, +0.091] -> no meaningful excess co-expression
  leave out gse84133_human  pooled E -0.006  [-0.063, +0.051] -> no meaningful excess co-expression
  leave out gse86469        pooled E +0.009  [-0.090, +0.109] -> inconclusive
  leave out gse81547        pooled E -0.004  [-0.107, +0.099] -> inconclusive
  leave out emtab5061       pooled E -0.015  [-0.092, +0.061] -> no meaningful excess co-expression

## nondiabetic, adjusted mode: 4 qualifying dataset(s): gse81608, gse86469, gse81547, emtab5061
* REML+HK (primary)              pooled E +0.003  95% CI [-0.085, +0.090]  tau2 0.0024  I2 81%  -> no meaningful excess co-expression
* REML+HK, factor floored at 1   pooled E +0.003  95% CI [-0.085, +0.090]  tau2 0.0024  I2 81%  -> no meaningful excess co-expression
* DL+HK                          pooled E +0.003  95% CI [-0.084, +0.090]  tau2 0.0022  I2 81%  -> no meaningful excess co-expression
* DL, normal CI                  pooled E +0.003  95% CI [-0.049, +0.054]  tau2 0.0022  I2 81%  -> no meaningful excess co-expression
* fixed effect                   pooled E +0.003  95% CI [-0.020, +0.025]  tau2 0.0000  I2 81%  -> no meaningful excess co-expression
  leave out gse81608        pooled E +0.004  [-0.163, +0.171] -> inconclusive
  leave out gse86469        pooled E +0.025  [-0.049, +0.100] -> no meaningful excess co-expression
  leave out gse81547        pooled E -0.005  [-0.166, +0.156] -> inconclusive
  leave out emtab5061       pooled E -0.015  [-0.138, +0.107] -> inconclusive

## nondiabetic, unadjusted mode: 4 qualifying dataset(s): gse81608, gse86469, gse81547, emtab5061
* REML+HK (primary)              pooled E -0.004  95% CI [-0.075, +0.067]  tau2 0.0012  I2 64%  -> no meaningful excess co-expression
* REML+HK, factor floored at 1   pooled E -0.004  95% CI [-0.075, +0.067]  tau2 0.0012  I2 64%  -> no meaningful excess co-expression
* DL+HK                          pooled E -0.004  95% CI [-0.075, +0.067]  tau2 0.0011  I2 64%  -> no meaningful excess co-expression
* DL, normal CI                  pooled E -0.004  95% CI [-0.045, +0.037]  tau2 0.0011  I2 64%  -> no meaningful excess co-expression
* fixed effect                   pooled E -0.002  95% CI [-0.027, +0.022]  tau2 0.0000  I2 64%  -> no meaningful excess co-expression
  leave out gse81608        pooled E -0.007  [-0.144, +0.131] -> inconclusive
  leave out gse86469        pooled E +0.015  [-0.025, +0.054] -> no meaningful excess co-expression
  leave out gse81547        pooled E -0.012  [-0.142, +0.119] -> inconclusive
  leave out emtab5061       pooled E -0.016  [-0.130, +0.097] -> inconclusive
