# WP4 extension A1-A3: positive control of the islet cross-check (reading rules fixed in PRESPEC_WP4b_positive_control_and_wider_search.md before this ran)

## A1. Pooled verdict over the five human datasets, per nominal rho (200 replicates; no injection = the primary result, pooled E +0.013 [-0.052, +0.078])

| mode | nominal rho | no meaningful excess | inconclusive | co-expressed more | pooled E mean (SD) | pooled CI mean |
|---|---|---|---|---|---|---|
| detected | 0.05 | 0.99 | 0.01 | 0.00 | +0.020 (0.003) | [-0.042, +0.083] |
| detected | 0.10 | 0.87 | 0.13 | 0.00 | +0.028 (0.004) | [-0.034, +0.090] |
| detected | 0.20 | 0.30 | 0.70 | 0.00 | +0.044 (0.006) | [-0.018, +0.106] |
| detected | 0.30 | 0.03 | 0.97 | 0.00 | +0.063 (0.007) | [-0.003, +0.129] |
| detected | 0.50 | 0.00 | 1.00 | 0.00 | +0.108 (0.008) | [+0.035, +0.182] |
| all | 0.05 | 0.00 | 1.00 | 0.00 | +0.119 (0.007) | [+0.071, +0.168] |
| all | 0.10 | 0.00 | 0.08 | 0.92 | +0.186 (0.009) | [+0.126, +0.246] |
| all | 0.20 | 0.00 | 0.00 | 1.00 | +0.290 (0.010) | [+0.225, +0.355] |
| all | 0.30 | 0.00 | 0.00 | 1.00 | +0.381 (0.009) | [+0.315, +0.447] |
| all | 0.50 | 0.00 | 0.00 | 1.00 | +0.582 (0.010) | [+0.523, +0.641] |

**Reading rule R1 (mode "detected", nominal rho 0.3, false-reassurance rate <= 20%): MET: the cross-check is SENSITIVE** (false-reassurance rate at 0.3: 0.03).
R2 not needed.

## A1. Per dataset (recovered excess and attenuation; mode "detected" unless stated)

| dataset | mode | nominal rho | E primary | E recovered mean (SD) | attenuation | share with empirical P < 0.05 |
|---|---|---|---|---|---|---|
| gse81608 | detected | 0.05 | -0.023 | -0.013 (0.005) | 0.20 | 0.00 |
| gse81608 | detected | 0.10 | -0.023 | -0.005 (0.007) | 0.18 | 0.00 |
| gse81608 | detected | 0.20 | -0.023 | +0.012 (0.009) | 0.18 | 0.09 |
| gse81608 | detected | 0.30 | -0.023 | +0.030 (0.011) | 0.17 | 0.66 |
| gse81608 | detected | 0.50 | -0.023 | +0.071 (0.013) | 0.17 | 1.00 |
| gse81608 | all | 0.05 | -0.023 | +0.135 (0.013) | 3.16 | 1.00 |
| gse81608 | all | 0.10 | -0.023 | +0.216 (0.018) | 2.39 | 1.00 |
| gse81608 | all | 0.20 | -0.023 | +0.337 (0.015) | 1.78 | 1.00 |
| gse81608 | all | 0.30 | -0.023 | +0.439 (0.018) | 1.49 | 1.00 |
| gse81608 | all | 0.50 | -0.023 | +0.642 (0.018) | 1.21 | 1.00 |
| gse84133_human | detected | 0.05 | +0.192 | +0.191 (0.003) | -0.03 | 1.00 |
| gse84133_human | detected | 0.10 | +0.192 | +0.190 (0.004) | -0.02 | 1.00 |
| gse84133_human | detected | 0.20 | +0.192 | +0.189 (0.006) | -0.01 | 1.00 |
| gse84133_human | detected | 0.30 | +0.192 | +0.192 (0.008) | 0.00 | 1.00 |
| gse84133_human | detected | 0.50 | +0.192 | +0.205 (0.012) | 0.02 | 1.00 |
| gse84133_human | all | 0.05 | +0.192 | +0.319 (0.012) | 2.55 | 1.00 |
| gse84133_human | all | 0.10 | +0.192 | +0.365 (0.012) | 1.72 | 1.00 |
| gse84133_human | all | 0.20 | +0.192 | +0.421 (0.014) | 1.13 | 1.00 |
| gse84133_human | all | 0.30 | +0.192 | +0.478 (0.013) | 0.92 | 1.00 |
| gse84133_human | all | 0.50 | +0.192 | +0.620 (0.015) | 0.78 | 1.00 |
| gse86469 | detected | 0.05 | -0.020 | -0.011 (0.007) | 0.18 | 0.00 |
| gse86469 | detected | 0.10 | -0.020 | -0.002 (0.009) | 0.18 | 0.00 |
| gse86469 | detected | 0.20 | -0.020 | +0.017 (0.013) | 0.18 | 0.07 |
| gse86469 | detected | 0.30 | -0.020 | +0.035 (0.015) | 0.18 | 0.47 |
| gse86469 | detected | 0.50 | -0.020 | +0.089 (0.018) | 0.20 | 1.00 |
| gse86469 | all | 0.05 | -0.020 | +0.094 (0.018) | 2.28 | 1.00 |
| gse86469 | all | 0.10 | -0.020 | +0.158 (0.020) | 1.78 | 1.00 |
| gse86469 | all | 0.20 | -0.020 | +0.263 (0.023) | 1.40 | 1.00 |
| gse86469 | all | 0.30 | -0.020 | +0.355 (0.023) | 1.21 | 1.00 |
| gse86469 | all | 0.50 | -0.020 | +0.559 (0.021) | 1.05 | 1.00 |
| gse81547 | detected | 0.05 | +0.025 | +0.026 (0.005) | 0.02 | 0.01 |
| gse81547 | detected | 0.10 | +0.025 | +0.028 (0.007) | 0.03 | 0.12 |
| gse81547 | detected | 0.20 | +0.025 | +0.033 (0.010) | 0.04 | 0.38 |
| gse81547 | detected | 0.30 | +0.025 | +0.042 (0.011) | 0.06 | 0.75 |
| gse81547 | detected | 0.50 | +0.025 | +0.065 (0.014) | 0.07 | 0.98 |
| gse81547 | all | 0.05 | +0.025 | +0.109 (0.017) | 1.67 | 1.00 |
| gse81547 | all | 0.10 | +0.025 | +0.165 (0.020) | 1.40 | 1.00 |
| gse81547 | all | 0.20 | +0.025 | +0.257 (0.021) | 1.14 | 1.00 |
| gse81547 | all | 0.30 | +0.025 | +0.341 (0.020) | 1.02 | 1.00 |
| gse81547 | all | 0.50 | +0.025 | +0.539 (0.023) | 0.94 | 1.00 |
| emtab5061 | detected | 0.05 | +0.045 | +0.057 (0.006) | 0.24 | 1.00 |
| emtab5061 | detected | 0.10 | +0.045 | +0.069 (0.009) | 0.24 | 1.00 |
| emtab5061 | detected | 0.20 | +0.045 | +0.093 (0.013) | 0.23 | 1.00 |
| emtab5061 | detected | 0.30 | +0.045 | +0.121 (0.014) | 0.24 | 1.00 |
| emtab5061 | detected | 0.50 | +0.045 | +0.182 (0.017) | 0.25 | 1.00 |
| emtab5061 | all | 0.05 | +0.045 | +0.112 (0.014) | 1.34 | 1.00 |
| emtab5061 | all | 0.10 | +0.045 | +0.168 (0.016) | 1.23 | 1.00 |
| emtab5061 | all | 0.20 | +0.045 | +0.266 (0.021) | 1.09 | 1.00 |
| emtab5061 | all | 0.30 | +0.045 | +0.359 (0.020) | 1.01 | 1.00 |
| emtab5061 | all | 0.50 | +0.045 | +0.571 (0.023) | 0.96 | 1.00 |
| gse84133_mouse | detected | 0.05 | +0.049 | +0.051 (0.003) | 0.05 | 0.00 |
| gse84133_mouse | detected | 0.10 | +0.049 | +0.053 (0.004) | 0.04 | 0.00 |
| gse84133_mouse | detected | 0.20 | +0.049 | +0.055 (0.006) | 0.03 | 0.00 |
| gse84133_mouse | detected | 0.30 | +0.049 | +0.057 (0.009) | 0.03 | 0.01 |
| gse84133_mouse | detected | 0.50 | +0.049 | +0.068 (0.011) | 0.04 | 0.12 |
| gse84133_mouse | all | 0.05 | +0.049 | +0.242 (0.015) | 3.86 | 1.00 |
| gse84133_mouse | all | 0.10 | +0.049 | +0.293 (0.017) | 2.44 | 1.00 |
| gse84133_mouse | all | 0.20 | +0.049 | +0.372 (0.017) | 1.59 | 1.00 |
| gse84133_mouse | all | 0.30 | +0.049 | +0.443 (0.019) | 1.27 | 1.00 |
| gse84133_mouse | all | 0.50 | +0.049 | +0.621 (0.018) | 1.04 | 1.00 |

## A2. Ribosomal protein pairs (100 random sets of 11 pairs per dataset; excess E_RP over matched random pairs)

| dataset | RP genes eligible | mean E_RP (SD) | E_RP in SE units of the consensus E | share of sets with E_RP >= 0.10 | mean rho of the RP pairs |
|---|---|---|---|---|---|
| gse81608 | 77 | +0.162 (0.026) | 10.5 | 1.00 | +0.169 |
| gse84133_human | 77 | +0.215 (0.020) | 2.5 | 1.00 | +0.222 |
| gse86469 | 77 | +0.101 (0.028) | 4.6 | 0.52 | +0.194 |
| gse81547 | 76 | +0.444 (0.037) | 20.5 | 1.00 | +0.512 |
| emtab5061 | 75 | +0.153 (0.028) | 7.0 | 0.97 | +0.196 |
| gse84133_mouse | 72 | +0.177 (0.027) | 3.4 | 1.00 | +0.193 |

**Reading rule (mean E_RP > +0.10 in at least 4 of the 5 human datasets): MET** (5 of 5).

## A3. Margin sensitivity of the primary pooled interval [-0.052, +0.078] (the pre-set margin stays +/-0.10)

| margin | verdict of the primary interval |
|---|---|
| 0.050 | inconclusive |
| 0.075 | inconclusive |
| 0.100 | no meaningful excess co-expression |
| 0.150 | no meaningful excess co-expression |
| 0.200 | no meaningful excess co-expression |

The smallest margin at which the verdict is "no meaningful excess" is 0.078: the data exclude an excess in Fisher z above that value.
