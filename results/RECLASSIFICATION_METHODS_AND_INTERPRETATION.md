# Reclassification and incremental utility analyses

These analyses were motivated by the risk-reclassification framework used in the European Heart Journal CV-Age study, but were adapted to the BC_cardiotox data and MTI question.

## Why a 1-year horizon?
The source dataset contains variable follow-up. A fixed 1-year horizon provides an interpretable clinical time point. Women were included in the fixed-horizon analysis only if they developed CTRCD within 365 days or remained under observation beyond 365 days. This yielded 172 women with known 1-year status, including 16 events. Women censored before 1 year were not treated as event-free.

## Cross-fitted risk prediction
Risk probabilities were estimated using 100 repeated stratified 5-fold cross-validation partitions. Within each fold, MTI standardization, missing-value imputation, scaling, and logistic model fitting were estimated using training data only and applied to held-out women. The four compared models were LVEF, MTI + LVEF, a clinical model (age + LVEF + anthracycline + anti-HER2 therapy), and the clinical model + MTI.

## Reclassification
Continuous NRI was used as the primary reclassification measure because no validated absolute 1-year CTRCD probability categories exist for this cohort. Event NRI quantifies the balance of upward vs downward risk movement among women who developed CTRCD; non-event NRI quantifies the balance of downward vs upward movement among event-free women. IDI quantifies change in the separation of mean predicted risk between events and non-events. Patient-level stratified bootstrap resampling (5,000 replicates) was used for confidence intervals.

Risk-quartile matrices are descriptive only. Quartiles were defined from the cross-fitted baseline clinical-model risk distribution and are not proposed as clinical cutoffs.

## Decision curve analysis
Decision curves compare net benefit across 1-year risk thresholds from 2% to 30%. Thresholds are shown as a continuum and no single threshold is proposed for clinical use. At a 10% threshold, MTI + LVEF had net benefit 0.0323 vs 0.0142 for LVEF; clinical + MTI had 0.0239 vs 0.0084 for the clinical model.

## Preserved LVEF analysis
Among 261 women with baseline LVEF >=50%, the highest cohort-relative MTI tertile contained 15/83 CTRCD events (18.1%), compared with 11/178 (6.2%) in the lower two tertiles. In a Cox model additionally adjusting for baseline LVEF, the highest MTI tertile had HR 3.38 (95% CI 1.55-7.38; P=0.0022). This is an exploratory risk-enrichment analysis; the tertile threshold is not a validated clinical cutoff.

## Interpretation guardrails
- Continuous NRI estimates were positive (+0.378) but confidence intervals crossed zero, so improved reclassification is not established.
- IDI estimates were small and imprecise.
- Decision-curve results are exploratory and are best used to support potential surveillance utility, not to define treatment thresholds.
- The preserved-LVEF result is clinically intuitive but uses a cohort-relative MTI tertile and requires external validation.
