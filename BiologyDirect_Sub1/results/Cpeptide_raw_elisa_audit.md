## Raw ELISA audit: 2020-06-08_ELISA new.xlsx
readings: 60 (3 groups x 5 conditions x 4 wells); DNA values: 15 (one per group x condition)
X = W / group-mean DNA holds for 60 of 60 readings (violations: 0)
CTRL: DNA per well (basal, 2.8, 5.5, 22, 44 mM) = [0.0657, 0.2007, 0.0927, 0.0911, 0.2388]; mean of 5 = 0.1378; sheet group mean = 0.1378; mean of the four non-basal wells = 0.1558
DAPT: DNA per well (basal, 2.8, 5.5, 22, 44 mM) = [0.3126, 0.3305, 0.2947, 0.3849, 0.244]; mean of 5 = 0.3133; sheet group mean = 0.3133; mean of the four non-basal wells = 0.3135
DKK-1: DNA per well (basal, 2.8, 5.5, 22, 44 mM) = [0.4174, 0.3348, 0.3842, 0.3325, 0.299]; mean of 5 = 0.3536; sheet group mean = 0.3536; mean of the four non-basal wells = 0.3376
Fig. 1c DNA values for CTRL: [0.0911, 0.0927, 0.2007, 0.2388] | the four non-basal wells of the ELISA sheet: [0.0911, 0.0927, 0.2007, 0.2388] | identical: True
Fig. 1c DNA values for DAPT: [0.244, 0.2947, 0.3305, 0.3849] | the four non-basal wells of the ELISA sheet: [0.244, 0.2947, 0.3305, 0.3849] | identical: True
Fig. 1c DNA values for DKK-1: [0.299, 0.3325, 0.3348, 0.3842] | the four non-basal wells of the ELISA sheet: [0.299, 0.3325, 0.3348, 0.3842] | identical: True
columns of cpeptide_all_conditions.csv: ['condition', 'group', 'rep', 'value']; rows: 60

standards (dOD450-590): blank B0 = 0.026000000000000002; 100 pM = 0.044000000000000004; 200 = 0.068; 400 = 0.151; 800 = 0.398; 1600 = 1.143
dOD per pM at the standards: 100: 0.00044, 200: 0.00034, 400: 0.00038, 800: 0.00050, 1600: 0.00071  (workbook uses one slope, 0.0007 per pM, through zero)
conversion chain in the workbook: S = dOD/0.0007 (pM); T = S*5; U = T*300/4; V = U*3.1615e-6; W = V*300 (ng per 300 uL) -> W = 0.3556 * dOD/0.0007 ... (a factor 75 = 300/4 is applied to a concentration)

| group | mM | dOD mean (range) | readings < 100-pM standard (dOD 0.044) | ng/300 uL (workbook) mean (range) | per group-mean DNA | well DNA (ug) | group-mean DNA (ug) | ratio to CTRL: absolute / per group-mean DNA |
|---|---|---|---|---|---|---|---|---|
| CTRL | 0 | 0.051 (0.051-0.052) | 0/4 | 26.2 (25.9-26.4) | 189.9 | 0.0657 | 0.1378 | 1.00 / 1.00 |
| CTRL | 2.8 | 0.051 (0.045-0.061) | 0/4 | 25.7 (22.9-31.0) | 186.2 | 0.2007 | 0.1378 | 1.00 / 1.00 |
| CTRL | 5.5 | 0.056 (0.047-0.083) | 0/4 | 28.6 (23.9-42.2) | 207.4 | 0.0927 | 0.1378 | 1.00 / 1.00 |
| CTRL | 22 | 0.054 (0.053-0.056) | 0/4 | 27.4 (26.9-28.5) | 199.1 | 0.0911 | 0.1378 | 1.00 / 1.00 |
| CTRL | 44 | 0.046 (0.039-0.066) | 3/4 | 23.6 (19.8-33.5) | 171.4 | 0.2388 | 0.1378 | 1.00 / 1.00 |
| DAPT | 0 | 0.038 (0.037-0.040) | 4/4 | 19.4 (18.8-20.3) | 62.0 | 0.3126 | 0.3133 | 0.74 / 0.33 |
| DAPT | 2.8 | 0.037 (0.037-0.038) | 4/4 | 19.1 (18.8-19.3) | 60.8 | 0.3305 | 0.3133 | 0.74 / 0.33 |
| DAPT | 5.5 | 0.102 (0.101-0.105) | 0/4 | 51.8 (51.3-53.4) | 165.4 | 0.2947 | 0.3133 | 1.81 / 0.80 |
| DAPT | 22 | 0.038 (0.037-0.040) | 4/4 | 19.4 (18.8-20.3) | 62.0 | 0.3849 | 0.3133 | 0.71 / 0.31 |
| DAPT | 44 | 0.044 (0.043-0.045) | 1/4 | 22.4 (21.8-22.9) | 71.3 | 0.2440 | 0.3133 | 0.95 / 0.42 |
| DKK-1 | 0 | 0.043 (0.040-0.046) | 2/4 | 21.8 (20.3-23.4) | 61.8 | 0.4174 | 0.3536 | 0.83 / 0.33 |
| DKK-1 | 2.8 | 0.038 (0.037-0.040) | 4/4 | 19.6 (18.8-20.3) | 55.3 | 0.3348 | 0.3536 | 0.76 / 0.30 |
| DKK-1 | 5.5 | 0.089 (0.087-0.090) | 0/4 | 45.2 (44.2-45.7) | 127.9 | 0.3842 | 0.3536 | 1.58 / 0.62 |
| DKK-1 | 22 | 0.049 (0.049-0.050) | 0/4 | 25.2 (24.9-25.4) | 71.1 | 0.3325 | 0.3536 | 0.92 / 0.36 |
| DKK-1 | 44 | 0.054 (0.051-0.056) | 0/4 | 27.4 (25.9-28.5) | 77.6 | 0.2990 | 0.3536 | 1.16 / 0.45 |

readings with an unusually high OD590 reference (> 0.05; typical 0.036-0.040): CTRL Glu5.5 rep 2: OD450 0.188, OD590 0.105, X = 306.0; CTRL Glu22 rep 4: OD450 0.106, OD590 0.053, X = 195.4
