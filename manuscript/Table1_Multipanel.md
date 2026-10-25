# Table 1. Cohort characteristics, MTI phenotype, validation, comparator models, and added clinical value

## Panel A. Pretreatment TDI analysis cohort

| Characteristic | Overall (n=270) | CTRCD (n=27) | No CTRCD (n=243) |
|---|---:|---:|---:|
| Age, y | 58 [49-66] | 61 [54-64] | 57 [49-66] |
| Heart rate, beats/min | 72.5 [64.3-80.0] | 78 [71-87] | 72 [64-79] |
| LVEF, % | 65.7 [61.2-70.9] | 61.4 [58.4-66.9] | 65.8 [61.6-72.0] |
| LV end-diastolic diameter, cm | 4.29 [3.97-4.63] | 4.22 [3.81-4.64] | 4.30 [3.97-4.63] |
| LV end-systolic diameter, cm | 2.70 [2.44-2.94] | 2.85 [2.37-3.12] | 2.69 [2.45-2.93] |
| Follow-up, d | 437 [220-874] | 316 [209-720] | 453 [224-879] |
| Anthracycline therapy | 163/248 (65.7%) | 20/24 (83.3%) | 143/224 (63.8%) |
| Anti-HER2 therapy | 158/248 (63.7%) | 17/24 (70.8%) | 141/224 (62.9%) |
| Hypertension | 62/247 (25.1%) | 6/24 (25.0%) | 56/223 (25.1%) |
| Dyslipidemia | 53/248 (21.4%) | 6/24 (25.0%) | 47/224 (21.0%) |
| Diabetes mellitus | 14/248 (5.6%) | 0/24 (0%) | 14/224 (6.2%) |
| Prior anthracycline exposure | 29/248 (11.7%) | 7/24 (29.2%) | 22/224 (9.8%) |
| Prior anti-HER2 exposure | 7/248 (2.8%) | 3/24 (12.5%) | 4/224 (1.8%) |

Values are median [IQR] or n/N (%). Denominators are shown for variables with missing values.

## Panel B. Four-component Mechanical Transition Instability phenotype

| MTI component | Clinical interpretation | CTRCD median [IQR] | No CTRCD median [IQR] | AUROC (95% CI) | P |
|---|---|---:|---:|---:|---:|
| Transition slope deficit | Slower resolution of systolic velocity | 120.5 [103.3-164.1] | 151.0 [122.4-181.9] | 0.650 (0.528-0.773) | 0.011 |
| Transition curvature deficit | Blunted systolic-to-relaxation bend | 742.0 [540.3-930.9] | 909.9 [717.7-1146.5] | 0.652 (0.531-0.766) | 0.010 |
| Early-relaxation latency | Later E′ after systolic termination | 0.175 [0.150-0.200] | 0.152 [0.132-0.175] | 0.632 (0.510-0.742) | 0.024 |
| E′ attenuation | Smaller early-diastolic velocity | 9.15 [7.89-10.48] | 10.51 [8.11-12.59] | 0.621 (0.519-0.714) | 0.039 |

P values are exploratory two-sided Mann-Whitney tests. Components were oriented so that higher standardized values indicated greater abnormality in the MTI score.

## Panel C. Primary cross-fitted performance and time-to-event association

| Analysis | Estimate |
|---|---|
| MTI AUROC | 0.695 (95% CI, 0.590-0.794) |
| MTI AUPRC | 0.192 (95% CI, 0.141-0.321) |
| MTI Harrell C-index | 0.677 (95% CI, 0.570-0.784) |
| Baseline LVEF AUROC | 0.657 (95% CI, 0.549-0.759) |
| Published E′ timing/velocity comparator AUROC | 0.663 (95% CI, 0.574-0.747) |
| MTI + LVEF AUROC | 0.712 (95% CI, 0.608-0.807) |
| MTI + LVEF C-index | 0.715 (95% CI, 0.611-0.810) |
| Paired C-index difference, MTI + LVEF vs LVEF | +0.037 (95% CI, -0.059 to 0.137) |
| MTI HR per 1 SD, unadjusted | 2.37 (95% CI, 1.49-3.78); P<0.001 |
| MTI HR per 1 SD, adjusted for LVEF | 2.20 (95% CI, 1.36-3.55); P=0.001 |
| MTI HR per 1 SD, adjusted for age + LVEF + treatment | 1.95 (95% CI, 1.19-3.18); P=0.008 |
| Preserved LVEF sensitivity cohort, MTI HR per 1 SD | 2.22 (95% CI, 1.41-3.51); P<0.001 |

Repeated prediction metrics used 100 stratified 5-fold cross-validation repeats. All preprocessing and model fitting were restricted to training folds.

## Panel D. Baseline model comparisons

| Model | AUROC (95% CI) | AUPRC (95% CI) |
|---|---:|---:|
| MTI + LVEF | 0.705 (0.598-0.804) | 0.213 (0.150-0.345) |
| Clinical + MTI | 0.704 (0.590-0.806) | 0.235 (0.159-0.392) |
| MTI | 0.684 (0.566-0.790) | 0.186 (0.134-0.306) |
| Clinical model | 0.659 (0.554-0.764) | 0.181 (0.123-0.321) |
| Fractional shortening | 0.658 (0.557-0.747) | 0.156 (0.120-0.237) |
| Published E′ timing/velocity metric | 0.648 (0.560-0.732) | 0.141 (0.114-0.211) |
| LVEF | 0.645 (0.536-0.745) | 0.156 (0.118-0.254) |
| Age + LVEF | 0.633 (0.518-0.739) | 0.152 (0.114-0.244) |
| E′ amplitude | 0.613 (0.505-0.704) | 0.138 (0.104-0.228) |
| Expanded TDI landmarks | 0.601 (0.511-0.693) | 0.130 (0.104-0.203) |
| E′ timing | 0.592 (0.465-0.707) | 0.147 (0.107-0.233) |
| S′ + E′ + A′ landmarks | 0.592 (0.483-0.696) | 0.135 (0.103-0.231) |
| FPCA-5* | 0.619 | 0.138 |
| Waveform ridge-25* | 0.533 | 0.106 |

Clinical model: age, LVEF, anthracycline therapy, and anti-HER2 therapy. *FPCA-5 and waveform ridge-25 used 20 repeated 5-fold cross-validation because of computational cost; all other benchmark models used 100 repeats.

## Panel E. Exploratory 1-year added-value and preserved-LVEF analyses

| Analysis | Result |
|---|---|
| Known 1-year status | 172 women / 16 events |
| LVEF AUROC at 1 y | 0.652 |
| MTI + LVEF AUROC at 1 y | 0.698 |
| Continuous NRI, LVEF → MTI + LVEF | +0.378 (95% CI, -0.122 to 0.869) |
| IDI, LVEF → MTI + LVEF | +0.015 (95% CI, -0.003 to 0.034) |
| Clinical model AUROC at 1 y | 0.618 |
| Clinical + MTI AUROC at 1 y | 0.653 |
| Clinical model AUPRC → Clinical + MTI AUPRC | 0.144 → 0.175 |
| Net benefit at 10% 1-y threshold, LVEF | 0.014 |
| Net benefit at 10% 1-y threshold, MTI + LVEF | 0.032 |
| Baseline LVEF ≥50%, highest MTI tertile | 15/83 events (18.1%) |
| Baseline LVEF ≥50%, lower two MTI tertiles | 11/178 events (6.2%) |
| Highest vs lower MTI tertiles, LVEF-adjusted HR | 3.38 (95% CI, 1.55-7.38); P=0.002 |

NRI, IDI, decision-curve, and tertile analyses are exploratory. There are no validated absolute CTRCD risk categories or MTI cutoffs in this cohort.

Abbreviations: AUPRC, area under the precision-recall curve; AUROC, area under the receiver-operating-characteristic curve; CTRCD, cancer therapy-related cardiac dysfunction; FPCA, functional principal component analysis; HR, hazard ratio; IDI, integrated discrimination improvement; LVEF, left ventricular ejection fraction; MTI, Mechanical Transition Instability; NRI, net reclassification improvement; TDI, Tissue Doppler imaging.
