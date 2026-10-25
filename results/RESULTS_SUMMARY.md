# Analysis results summary

## Cohort

- Full BC_cardiotox cohort: **531 women, 54 CTRCD events**
- TDI analysis cohort: **270 women, 27 CTRCD events (10.0%)**
- Baseline LVEF >=50% TDI sensitivity cohort: **261 women, 26 events**
- Median follow-up in the TDI cohort: **437 days**
- Median time to CTRCD among events: **316 days**

## Primary MTI performance

| Model | AUROC (95% CI) | AUPRC | C-index (95% CI) | 1-year AUROC |
|---|---:|---:|---:|---:|
| MTI | 0.695 (0.590–0.794) | 0.192 | 0.677 (0.570–0.784) | 0.657 |
| LVEF | 0.657 (0.549–0.759) | 0.166 | 0.678 (0.578–0.767) | 0.672 |
| Published vE′ comparator | 0.663 (0.574–0.747) | 0.150 | 0.660 (0.570–0.750) | 0.637 |
| MTI + LVEF | 0.712 (0.608–0.807) | 0.211 | 0.715 (0.611–0.810) | 0.707 |

The paired C-index difference for MTI + LVEF versus LVEF was **+0.037** with 95% CI **−0.059 to 0.137**. The corresponding paired AUROC difference was **+0.055** with 95% CI **−0.040 to 0.153**. Therefore, discrimination superiority is not established.

## Cox associations

- MTI unadjusted: HR per 1 SD **2.37** (95% CI 1.49–3.78; P<0.001)
- MTI adjusted for LVEF: HR **2.20** (95% CI 1.36–3.55; P=0.0013)
- MTI adjusted for age and LVEF: HR **2.19** (95% CI 1.34–3.58; P=0.0017)
- MTI adjusted for age, LVEF, anthracycline therapy, and anti-HER2 therapy: HR **1.95** (95% CI 1.19–3.18; P=0.008), n=248, 24 events
- Preserved baseline LVEF >=50% sensitivity cohort: MTI HR **2.22** (95% CI 1.41–3.51; P<0.001)

## MTI component results

| Component | Event-oriented AUROC | P value |
|---|---:|---:|
| Transition slope deficit | 0.650 | 0.0105 |
| Transition curvature deficit | 0.652 | 0.0096 |
| Early-relaxation latency | 0.632 | 0.0243 |
| E′ attenuation | 0.621 | 0.0386 |

## Tertiles

- Low MTI: 5/90 events
- Middle MTI: 6/90 events
- High MTI: 16/90 events

## Robustness

Across Savitzky-Golay smoothing windows from 21 to 61 samples, MTI AUROC ranged from **0.686 to 0.710**, C-index from **0.673 to 0.694**, and the unadjusted HR per 1 SD from **2.10 to 2.66**. Leave-one-component-out analysis did not identify a single component responsible for the full score.

## Important limitations

The study is internally validated but not externally validated. Only 27 CTRCD events are available. The same source cohort informed MTI development and evaluation. TDI was processed in only 270/531 women, and those with TDI differed from those without processed TDI in age, anti-HER2 exposure, and follow-up duration. No GLS endpoint is available.
