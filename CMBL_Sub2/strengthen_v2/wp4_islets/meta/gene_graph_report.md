# WP4 add-on (post hoc): meta-analysis with a null that keeps the gene graph of the consensus pairs

Same estimators and margin as meta_analysis.py; E_graph and SE_graph replace E and SE of the pre-specified analysis. The rule of 4.4 stays defined on the pre-specified SE; this is a sensitivity analysis.

* REML+HK                        pooled E +0.012  95% CI [-0.050, +0.074]  tau2 0.0010  I2 67%  -> no meaningful excess co-expression
* REML+HK, factor floored at 1   pooled E +0.012  95% CI [-0.050, +0.074]  tau2 0.0010  I2 67%  -> no meaningful excess co-expression
* DL+HK                          pooled E +0.013  95% CI [-0.051, +0.077]  tau2 0.0012  I2 67%  -> no meaningful excess co-expression
* DL, normal CI                  pooled E +0.013  95% CI [-0.026, +0.052]  tau2 0.0012  I2 67%  -> no meaningful excess co-expression
* fixed effect                   pooled E +0.003  95% CI [-0.017, +0.023]  tau2 0.0000  I2 67%  -> no meaningful excess co-expression

  leave out gse81608        pooled E +0.025  [-0.058, +0.108] -> inconclusive
  leave out gse84133_human  pooled E +0.005  [-0.049, +0.059] -> no meaningful excess co-expression
  leave out gse86469        pooled E +0.024  [-0.068, +0.116] -> inconclusive
  leave out gse81547        pooled E +0.012  [-0.089, +0.112] -> inconclusive
  leave out emtab5061       pooled E -0.004  [-0.076, +0.069] -> no meaningful excess co-expression