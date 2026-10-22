# Mechanical Transition Instability on Pretreatment Tissue Doppler Imaging Identifies Risk of Cancer Therapy-Related Cardiac Dysfunction in Women With Breast Cancer

## 1. Introduction

Cardio-oncology surveillance is designed to identify cardiovascular vulnerability early enough to preserve both cardiac health and effective cancer treatment. Yet pretreatment risk assessment remains difficult. Left ventricular ejection fraction (LVEF) is central to routine evaluation, but it summarizes chamber-level systolic performance and may not capture subtle abnormalities in the timing of myocardial motion. Global longitudinal strain improves sensitivity to evolving dysfunction, but its principal clinical role is surveillance after therapy begins.¹ These limitations motivate complementary baseline markers that describe how the myocardium moves, rather than only how much the ventricle ejects.

Tissue Doppler imaging (TDI) provides a beat-level velocity trajectory spanning systole and diastole. In the public BC_cardiotox cohort, later and lower-amplitude pretreatment E′ velocity was associated with subsequent cancer therapy-related cardiac dysfunction (CTRCD).²,³ We hypothesized that susceptibility may be encoded more broadly in the transition from systolic motion into early relaxation. We therefore developed Mechanical Transition Instability (MTI), an interpretable 4-component phenotype derived from the pretreatment TDI trajectory, and evaluated its association with subsequent CTRCD, its relationship to baseline LVEF, and its potential for surveillance risk enrichment.

## 2. Methods

We performed a secondary analysis of BC_cardiotox, a publicly available cohort of 531 women with HER2-positive breast cancer treated with potentially cardiotoxic therapy at the University Hospital of A Coruña.² The source study was approved by the Drug Research Ethical Committee of Galicia (CEIm-G; 2022/492), conducted in accordance with the Declaration of Helsinki, and received a waiver of consent for the retrospective analysis. The dataset is distributed under a Creative Commons Attribution 4.0 license. A processed pretreatment TDI waveform was available for 270 women and defined the imaging cohort. Each waveform represented one cardiac cycle sampled at 1,001 normalized phases. CTRCD was defined in the source dataset as LVEF <50% with an absolute decrease of at least 10 percentage points from baseline. Time was measured from treatment initiation to CTRCD or censoring.

Waveforms were smoothed with a third-order Savitzky-Golay filter. The systolic S peak, first positive-to-negative zero crossing after S, and subsequent E′ nadir were identified automatically. MTI combined 4 study-defined measures: transition slope deficit, reduced steepness from peak systolic velocity toward zero; transition curvature deficit, blunting of the systolic-to-relaxation bend; early-relaxation latency, delayed E′ after systolic termination; and E′ attenuation, reduced early-diastolic velocity magnitude (Figure 1). Components were oriented so that higher values indicated greater abnormality, standardized within training folds, and averaged with equal weight.

The primary association analysis used Cox proportional hazards regression and reported hazard ratios (HRs) per 1-SD increase in MTI. Models evaluated MTI alone and after adjustment for baseline LVEF, age, and cancer-therapy exposures. Discrimination used 100 repeated stratified 5-fold cross-validation partitions. MTI standardization and model coefficients were estimated in training folds only and applied to held-out women. Repeated held-out predictions were averaged to calculate AUROC, AUPRC, and Harrell C-index; 95% CIs used patient bootstrap resampling. Comparators were baseline LVEF and the previously published E′ timing/velocity ratio.³ Prespecified sensitivity analyses restricted the cohort to baseline LVEF ≥50%, varied signal smoothing, and omitted MTI components individually.

Exploratory analyses evaluated whole-waveform morphology with functional permutation testing, the shape of the continuous MTI-risk relation with restricted cubic splines, and surveillance efficiency using cross-fitted risk rankings. The latter asked what proportion of subsequent events would fall within a prespecified fraction of patients selected for intensified surveillance. These analyses were considered hypothesis-generating and did not define a clinical threshold.

## 3. Results

Among 270 women with pretreatment TDI, 27 (10.0%) developed CTRCD. Median follow-up was 437 days, and median time to CTRCD among event cases was 316 days. Women who developed CTRCD had lower baseline LVEF and higher heart rate (Table 1A). All 4 MTI components were directionally associated with outcome, with event-oriented AUROCs from 0.621 to 0.652 (Table 1B). Pretreatment waveforms showed a distributed difference across the systolic-to-relaxation transition; the global functional permutation test was P=0.047 and was attenuated after accounting for baseline LVEF (P=0.056), without a narrow phase interval surviving cluster-wise family-error correction (Figure 1).

MTI yielded a cross-fitted AUROC of 0.695 (95% CI, 0.590 to 0.794), AUPRC of 0.192 (95% CI, 0.141 to 0.321), and C-index of 0.677 (95% CI, 0.570 to 0.784). Baseline LVEF yielded an AUROC of 0.657 and C-index of 0.678; the published E′ timing/velocity comparator yielded an AUROC of 0.663 and C-index of 0.660. MTI plus LVEF yielded an AUROC of 0.712 (95% CI, 0.608 to 0.807) and C-index of 0.715 (95% CI, 0.611 to 0.810). The paired C-index difference versus LVEF was 0.037 (95% CI, -0.059 to 0.137), so discrimination superiority was not established (Table 1C; Figure 2).

Each 1-SD higher MTI was associated with greater hazard of CTRCD (HR, 2.37; 95% CI, 1.49 to 3.78; P<0.001). The association persisted after adjustment for LVEF (HR, 2.20; 95% CI, 1.36 to 3.55; P=0.001), age and LVEF (HR, 2.19; 95% CI, 1.34 to 3.58; P=0.002), and, among 248 complete cases, anthracycline and anti-HER2 exposure in addition to age and LVEF (HR, 1.95; 95% CI, 1.19 to 3.18; P=0.008). Among 261 women with baseline LVEF ≥50%, MTI remained associated with CTRCD (adjusted HR, 2.22; 95% CI, 1.41 to 3.51; P<0.001). CTRCD occurred in 5/90, 6/90, and 16/90 women across increasing MTI tertiles. A spline analysis showed no evidence of nonlinearity (P=0.999). In an exploratory surveillance-efficiency analysis, the highest-risk 20% by cross-fitted MTI plus LVEF prediction contained 13 of 27 events (48.1%), with a positive predictive value of 24.1% versus 10.0% cohort prevalence (Figure 2).

## 4. Discussion and Conclusion

Pretreatment myocardial motion contained prognostic information before overt systolic dysfunction. A slower, smoother, later, and attenuated transition from systolic motion into early relaxation was associated with subsequent CTRCD, independent of baseline LVEF and among women with preserved pretreatment LVEF. These findings extend prior observations centered on the E′ landmark³ by treating the systolic-to-relaxation transition as an integrated mechanical phenotype.

MTI is not intended to replace LVEF, GLS, biomarkers, or established clinical risk assessment. Its potential value is different: a single pretreatment TDI cycle may provide an interpretable measure of temporal mechanical vulnerability before therapy-related decline is measurable. The continuous risk gradient and event-concentration analysis support a risk-enrichment rather than diagnostic role. The numerically higher discrimination of MTI plus LVEF also suggests complementarity, but paired uncertainty crossed zero and incremental superiority was not established.

Limitations include 27 events, single-center derivation, internal rather than external validation, absence of GLS, and processed TDI availability in only 270 of 531 women. Women with available TDI differed from those without TDI in age, follow-up duration, anti-HER2 exposure, and hypertension, although CTRCD prevalence was similar. Cardiac-cycle normalization removes absolute timing information, and MTI depends on automated landmark extraction and smoothing. External validation in prospectively acquired TDI and contemporary cardio-oncology cohorts is required.

Mechanical Transition Instability from a single pretreatment TDI cycle was independently associated with subsequent CTRCD. Temporal myocardial mechanics may provide a practical layer of baseline risk information for identifying patients who warrant more intensive cardiovascular surveillance.

## Ethics Statement

The source study was approved by the Drug Research Ethical Committee of Galicia (CEIm-G; 2022/492), was conducted in accordance with the Declaration of Helsinki, and received a waiver of consent for the retrospective analysis.²

## Data and Code Availability

BC_cardiotox is publicly available through Figshare (doi:10.6084/m9.figshare.22650748.v4) under CC BY 4.0. Analysis code, aggregate results, figure-generation scripts, and the interactive research companion are provided in the accompanying repository. Patient-level data are not redistributed in the repository.

## Declaration of generative AI and AI-assisted technologies in the manuscript preparation process

During preparation of this work, the authors used ChatGPT (OpenAI) to assist with language editing, organization, and clarity. The authors reviewed and edited the content and take full responsibility for the article.

## References

1. Lyon AR, López-Fernández T, Couch LS, et al. 2022 ESC Guidelines on cardio-oncology developed in collaboration with the EHA, ESTRO and IC-OS. *Eur Heart J*. 2022;43:4229-4361. doi:10.1093/eurheartj/ehac244.
2. Piñeiro-Lamas B, López-Cheda A, Cao R, et al. A cardiotoxicity dataset for breast cancer patients. *Sci Data*. 2023;10:527. doi:10.1038/s41597-023-02419-1.
3. Piñeiro-Lamas B, López-Cheda A, Cao R, Lesta-Mellid R, Bouzas-Mosquera A, Barbeito-Caamaño C. Tissue Doppler Imaging as a predictor of therapy-related cardiac dysfunction in breast cancer patients. *Int J Cardiovasc Imaging*. 2025;41:1505-1512. doi:10.1007/s10554-025-03439-1.
4. Sawaya H, Sebag IA, Plana JC, et al. Early detection and prediction of cardiotoxicity in chemotherapy-treated patients. *Am J Cardiol*. 2011;107:1375-1380.
5. Negishi K, Negishi T, Hare JL, Haluska BA, Plana JC, Marwick TH. Independent and incremental value of deformation indices for prediction of trastuzumab-induced cardiotoxicity. *J Am Soc Echocardiogr*. 2013;26:493-498.
