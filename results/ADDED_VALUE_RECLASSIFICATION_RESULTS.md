# Added-value and reclassification analyses

## Fixed 1-year analysis cohort
- 172 women had known 1-year status; 16 developed CTRCD by 1 year.
- Predictions were generated with 100 repeated stratified 5-fold cross-validation. MTI standardization, imputation, scaling, and model fitting occurred in training folds only.

## 1-year discrimination and calibration
| model        |   n |   events |   AUROC_1y |   AUPRC_1y |   Brier_1y |
|:-------------|----:|---------:|-----------:|-----------:|-----------:|
| LVEF         | 172 |       16 |   0.652244 |   0.143325 |  0.0843062 |
| MTI+LVEF     | 172 |       16 |   0.697516 |   0.170657 |  0.0837779 |
| Clinical     | 172 |       16 |   0.618189 |   0.143908 |  0.0850633 |
| Clinical+MTI | 172 |       16 |   0.653045 |   0.174846 |  0.0841744 |

## Continuous reclassification
- LVEF -> MTI+LVEF: continuous NRI +0.378 (95% bootstrap CI -0.122 to +0.869); IDI +0.015 (-0.003 to +0.034).
- Clinical -> Clinical+MTI: continuous NRI +0.378 (95% bootstrap CI -0.135 to +0.853); IDI +0.016 (-0.007 to +0.039).
- These estimates are directionally favorable but imprecise and are exploratory.

## Decision-curve interpretation
- At a 10% 1-year risk threshold, MTI+LVEF net benefit was 0.0323 versus 0.0142 for LVEF alone. Clinical+MTI net benefit was 0.0239 versus 0.0084 for the clinical model.
- These threshold analyses are exploratory and should not be treated as validated surveillance cutoffs.

## Preserved-LVEF phenotype
- Among 261 women with baseline LVEF >=50%, the highest MTI tertile had 15/83 events (18.1%) versus 11/178 (6.2%) in the lower two tertiles.
- After adjustment for baseline LVEF, the highest MTI tertile had HR 3.38 (95% CI 1.55-7.38; P=0.0022) versus the lower two tertiles. This tertile analysis is exploratory and the threshold is cohort-relative.
