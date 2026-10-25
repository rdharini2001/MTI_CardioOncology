Pretreatment Tissue Doppler Mechanical Transition Instability and Subsequent Cancer Therapy-Related Cardiac Dysfunction in HER2-Positive Breast Cancer
1. Introduction
Cardiovascular surveillance during potentially cardiotoxic cancer therapy remains challenging because early myocardial dysfunction may develop while left ventricular ejection fraction (LVEF) is still preserved. LVEF and global longitudinal strain (GLS) are recommended for baseline echocardiographic assessment, but neither fully captures pretreatment susceptibility to cancer therapy-related cardiac dysfunction (CTRCD).¹,² This raises an important clinical question: can subtle abnormalities in myocardial mechanics identify vulnerability before conventional systolic function declines?
Tissue Doppler imaging (TDI) provides a continuous record of myocardial velocity across systole and diastole. Prior studies have shown that a later and lower-amplitude pretreatment E′ wave is associated with subsequent CTRCD.³,⁴ Deformation imaging studies similarly suggest that abnormalities in myocardial mechanics can precede a fall in LVEF.⁵,⁶ These observations suggest that risk may be encoded not only in individual velocity peaks, but also in how myocardial motion transitions from systolic contraction into early relaxation. We therefore developed Mechanical Transition Instability (MTI), an interpretable 4-component phenotype of the pretreatment TDI trajectory and evaluated whether it was associated with subsequent CTRCD and provided information complementary to conventional baseline assessment.
2. Methods
We analyzed data from the Cardiotox cohort, which included 531 women with HER2-positive breast cancer treated with potentially cardiotoxic therapy at the University Hospital of A Coruña between 2007 and 2021.³ Pretreatment clinical characteristics, echocardiographic measurements, cancer treatment exposures, and longitudinal follow-up were available. The present study included 270 women with an available pretreatment TDI waveform and follow-up for cardiac outcomes. The primary outcome was subsequent CTRCD, defined as an absolute decrease in LVEF of at least 10 percentage points to <50%.
Each waveform represented one cardiac cycle sampled at 1,001 normalized phases. After smoothing, we identified the systolic S peak, the first positive-to-negative zero crossing after S, and the subsequent E′ nadir. MTI combined 4 transient measurements: transition slope deficit, a slower fall in systolic velocity toward zero; transition curvature deficit, a blunted systolic-to-relaxation bend; early-relaxation latency, a later E′ after systolic termination; and E′ attenuation, a smaller early-diastolic velocity (Figure 1). Components were oriented so that higher values indicated greater abnormality and averaged after training-fold standardization.
We assessed the association between MTI and time to CTRCD using Cox proportional hazards regression and reported hazard ratios (HRs) per 1-SD increase in MTI. We evaluated discrimination using 100 repeated stratified 5-fold cross-validation, with all standardization and model fitting performed exclusively within the training folds. We averaged held-out predictions across repetitions and calculated the area under the receiver-operating-characteristic curve (AUROC), area under the precision-recall curve (AUPRC), and Harrell C-index, with 95% CIs estimated by bootstrap resampling. We compared MTI with baseline LVEF, fractional shortening, conventional TDI landmarks, the previously published E′ timing/velocity metric, clinical models, functional principal components, and a regularized whole-waveform model.
We also evaluated the incremental value of MTI at a fixed 1-year horizon among women with known 1-year outcome status. Using cross-fitted predictions, we assessed changes in discrimination, Brier score, continuous net reclassification improvement (NRI), integrated discrimination improvement (IDI), and decision-curve net benefit after adding MTI to baseline models. We performed sensitivity analyses in women with preserved baseline LVEF and examined the robustness of our findings to signal smoothing, omission of individual MTI components, and outcome-agnostic whole-waveform analyses.
3. Results
Among the 270 women with pretreatment TDI, 27 (10.0%) developed CTRCD. Median follow-up was 437 days, and the median time to CTRCD among women who developed an event was 316 days. Women who subsequently developed CTRCD had a lower baseline LVEF than those who did not (61.4% vs 65.8%) (Table 1A). Each of the 4 MTI components showed a directionally consistent association with outcome, with individual AUROCs ranging from 0.621 to 0.652 (Table 1B). Analysis of the complete TDI waveform supported a broader difference across the systolic-to-relaxation transition rather than a single isolated phase (global permutation P=0.047; after adjustment for LVEF, P=0.056) (Figure 1).
In repeated held-out cross-validation, MTI achieved an AUROC of 0.695 (95% CI, 0.590 to 0.794), an AUPRC of 0.192, and a C-index of 0.677 (95% CI, 0.570 to 0.784). By comparison, the AUROC was 0.657 for baseline LVEF and 0.663 for the previously published E′ metric. The combination of MTI and LVEF produced the highest point estimates, with an AUROC of 0.712 (95% CI, 0.608 to 0.807) and a C-index of 0.715 (95% CI, 0.611 to 0.810). However, the paired difference in C-index compared with LVEF alone was not statistically definitive (difference, 0.037; 95% CI, −0.059 to 0.137). In the harmonized held-out benchmark, MTI also showed higher discrimination than conventional TDI measures, clinical models, functional principal components, and the regularized whole-waveform model (Table 1C; Figure 2).
Higher MTI was consistently associated with subsequent CTRCD. For each 1-SD increase in MTI, the unadjusted HR was 2.37 (95% CI, 1.49 to 3.78; P<0.001). The association remained after adjustment for baseline LVEF (HR, 2.20; 95% CI, 1.36 to 3.55; P=0.001) and after further adjustment for age, LVEF, anthracycline exposure, and anti-HER2 therapy (HR, 1.95; 95% CI, 1.19 to 3.18; P=0.008). Among the 261 women with a baseline LVEF ≥50%, MTI remained associated with CTRCD after adjustment for LVEF (HR per 1-SD increase, 2.22; 95% CI, 1.41 to 3.51; P<0.001). In this preserved-LVEF subgroup, 15 of 83 women (18.1%) in the highest MTI tertile developed CTRCD compared with 11 of 178 (6.2%) in the lower two tertiles, corresponding to an LVEF-adjusted HR of 3.38 (95% CI, 1.55 to 7.38; P=0.002).
We next examined whether MTI added information at a fixed 1-year horizon. Among 172 women with known 1-year outcome status, including 16 events, adding MTI to LVEF increased the AUROC from 0.652 to 0.698. Adding MTI to the clinical model increased the AUROC from 0.618 to 0.653 and the AUPRC from 0.144 to 0.175. For LVEF plus MTI, the continuous NRI was +0.378 (95% CI, −0.122 to 0.869) and the IDI was +0.015 (95% CI, −0.003 to 0.034). At a 10% predicted-risk threshold, net benefit was 0.032 for MTI plus LVEF compared with 0.014 for LVEF alone. Although these findings consistently favored the addition of MTI, the reclassification confidence intervals crossed zero and the added-value analyses should therefore be considered exploratory.
4. Discussion and Conclusion
Several findings are important. Pretreatment TDI contained prognostic information before the development of a treatment-related decline in LVEF. This information was not confined to the E′ wave, but extended across the broader transition from systolic contraction to early relaxation and remained evident when MTI was compared with conventional echocardiographic measures and data-driven waveform models. Importantly, MTI was also associated with subsequent CTRCD among women with preserved baseline LVEF, suggesting that abnormalities in temporal myocardial mechanics may be present even when conventional systolic function appears normal.
These findings raise the possibility that MTI could be used for baseline risk enrichment in cardio-oncology if validated prospectively. TDI is acquired during standard echocardiography and does not require an additional imaging modality. A pretreatment TDI waveform could therefore provide complementary information about myocardial mechanical vulnerability before therapy begins. Among women with preserved baseline LVEF, those in the highest MTI tertile had nearly 3 times the event rate of women in the lower two tertiles. Analyses of reclassification and clinical net benefit also showed favorable point estimates when MTI was added to LVEF, although these findings remain exploratory.
MTI is not intended to replace LVEF, GLS, or established clinical risk assessment. Rather, it captures a different aspect of cardiac function. LVEF reflects overall ventricular systolic performance, whereas MTI characterizes the timing and coordination of myocardial motion as contraction transitions into relaxation. Models combining MTI with LVEF showed higher discrimination than LVEF alone, but the incremental improvement was not statistically definitive. These results therefore support MTI as a potentially complementary marker rather than a superior alternative to established measures.
Several limitations should be considered. The number of CTRCD events was modest, and the data were derived from a single-center cohort. GLS was not available for direct comparison. Complete processed TDI waveforms were available in 270 of 531 women, and differences between women with and without available TDI may have introduced selection bias. The TDI signals were also platform dependent. Normalization of each cardiac cycle removes information about absolute timing, and MTI depends on automated landmark detection and signal smoothing. In addition, the fixed 1-year reclassification analysis included only 16 events, limiting the precision of estimates of incremental clinical value. External prospective validation will therefore be necessary before MTI can be considered for clinical implementation.
In women with HER2-positive breast cancer, higher pretreatment Mechanical Transition Instability was associated with subsequent CTRCD, including among those with preserved baseline LVEF. These findings suggest that temporal features of myocardial motion may provide an interpretable layer of baseline risk information that complements conventional assessment of systolic function.
Ethics Statement
The source study was approved by the Drug Research Ethical Committee of Galicia (CEIm-G; 2022/492), was conducted in accordance with the Declaration of Helsinki, and received a waiver of consent for the retrospective analysis.3
Data and Code Availability
The Cardiotox data is publicly available through Figshare (doi:10.6084/m9.figshare.22650748.v4) under CC BY 4.0. Analysis code, aggregate results, figure-generation scripts, and the interactive research website are provided at https://github.com/rdharini2001/MTI_CardioOncology. Patient-level data are not redistributed in the repository.
Declaration of generative AI and AI-assisted technologies in the manuscript preparation process
During preparation of this work, the authors used ChatGPT (OpenAI) to assist with language editing, organization, and clarity. The authors reviewed and edited the content and take full responsibility for the article.
References
1. Lyon AR, López-Fernández T, Couch LS, et al. 2022 ESC Guidelines on cardio-oncology developed in collaboration with the EHA, ESTRO and IC-OS. Eur Heart J. 2022;43:4229-4361. doi:10.1093/eurheartj/ehac244.
2. Raisi-Estabragh Z, Murphy AC, Ramalingam S, et al. Cardiovascular considerations before cancer therapy: gaps in evidence and JACC: CardioOncology Expert Panel recommendations. JACC CardioOncol. 2024;6:631-654. doi:10.1016/j.jaccao.2024.07.017.
3. Piñeiro-Lamas B, López-Cheda A, Cao R, et al. A cardiotoxicity dataset for breast cancer patients. Sci Data. 2023;10:527. doi:10.1038/s41597-023-02419-1.
4. Piñeiro-Lamas B, López-Cheda A, Cao R, Lesta-Mellid R, Bouzas-Mosquera A, Barbeito-Caamaño C. Tissue Doppler Imaging as a predictor of therapy-related cardiac dysfunction in breast cancer patients. Int J Cardiovasc Imaging. 2025;41:1505-1512. doi:10.1007/s10554-025-03439-1.
5. Negishi K, Negishi T, Hare JL, Haluska BA, Plana JC, Marwick TH. Independent and incremental value of deformation indices for prediction of trastuzumab-induced cardiotoxicity. J Am Soc Echocardiogr. 2013;26:493-498.
6. Sawaya H, Sebag IA, Plana JC, et al. Early detection and prediction of cardiotoxicity in chemotherapy-treated patients. Am J Cardiol. 2011;107:1375-1380.


 


Table 1. Cohort characteristics, MTI phenotype, validation, comparator models, and added clinical value
Panel A. Pretreatment TDI analysis cohort
Characteristic	Overall (n=270)	CTRCD (n=27)	No CTRCD (n=243)
Age, y	58 [49-66]	61 [54-64]	57 [49-66]
Heart rate, beats/min	72.5 [64.3-80.0]	78 [71-87]	72 [64-79]
LVEF, %	65.7 [61.2-70.9]	61.4 [58.4-66.9]	65.8 [61.6-72.0]
LV end-diastolic diameter, cm	4.29 [3.97-4.63]	4.22 [3.81-4.64]	4.30 [3.97-4.63]
LV end-systolic diameter, cm	2.70 [2.44-2.94]	2.85 [2.37-3.12]	2.69 [2.45-2.93]
Follow-up, d	437 [220-874]	316 [209-720]	453 [224-879]
Anthracycline therapy	163/248 (65.7%)	20/24 (83.3%)	143/224 (63.8%)
Anti-HER2 therapy	158/248 (63.7%)	17/24 (70.8%)	141/224 (62.9%)
Hypertension	62/247 (25.1%)	6/24 (25.0%)	56/223 (25.1%)
Dyslipidemia	53/248 (21.4%)	6/24 (25.0%)	47/224 (21.0%)
Diabetes mellitus	14/248 (5.6%)	0/24 (0%)	14/224 (6.2%)
Prior anthracycline exposure	29/248 (11.7%)	7/24 (29.2%)	22/224 (9.8%)
Prior anti-HER2 exposure	7/248 (2.8%)	3/24 (12.5%)	4/224 (1.8%)
Values are median [IQR] or n/N (%). Denominators are shown for variables with missing values.
Panel B. Four-component Mechanical Transition Instability phenotype
MTI component	Clinical interpretation	CTRCD median [IQR]	No CTRCD median [IQR]	AUROC (95% CI)	P
Transition slope deficit	Slower resolution of systolic velocity	120.5 [103.3-164.1]	151.0 [122.4-181.9]	0.650 (0.528-0.773)	0.011
Transition curvature deficit	Blunted systolic-to-relaxation bend	742.0 [540.3-930.9]	909.9 [717.7-1146.5]	0.652 (0.531-0.766)	0.010
Early-relaxation latency	Later E′ after systolic termination	0.175 [0.150-0.200]	0.152 [0.132-0.175]	0.632 (0.510-0.742)	0.024
E′ attenuation	Smaller early-diastolic velocity	9.15 [7.89-10.48]	10.51 [8.11-12.59]	0.621 (0.519-0.714)	0.039
P values are exploratory two-sided Mann-Whitney tests. Components were oriented so that higher standardized values indicated greater abnormality in the MTI score.
Panel C. Primary cross-fitted performance and time-to-event association
Analysis	Estimate
MTI AUROC	0.695 (95% CI, 0.590-0.794)
MTI AUPRC	0.192 (95% CI, 0.141-0.321)
MTI Harrell C-index	0.677 (95% CI, 0.570-0.784)
Baseline LVEF AUROC	0.657 (95% CI, 0.549-0.759)
Published E′ timing/velocity comparator AUROC	0.663 (95% CI, 0.574-0.747)
MTI + LVEF AUROC	0.712 (95% CI, 0.608-0.807)
MTI + LVEF C-index	0.715 (95% CI, 0.611-0.810)
Paired C-index difference, MTI + LVEF vs LVEF	+0.037 (95% CI, -0.059 to 0.137)
MTI HR per 1 SD, unadjusted	2.37 (95% CI, 1.49-3.78); P<0.001
MTI HR per 1 SD, adjusted for LVEF	2.20 (95% CI, 1.36-3.55); P=0.001
MTI HR per 1 SD, adjusted for age + LVEF + treatment	1.95 (95% CI, 1.19-3.18); P=0.008
Preserved LVEF sensitivity cohort, MTI HR per 1 SD	2.22 (95% CI, 1.41-3.51); P<0.001
Repeated prediction metrics used 100 stratified 5-fold cross-validation repeats. All preprocessing and model fitting were restricted to training folds.
Panel D. Baseline model comparisons 
Model	AUROC (95% CI)	AUPRC (95% CI)
MTI + LVEF	0.705 (0.598-0.804)	0.213 (0.150-0.345)
Clinical + MTI	0.704 (0.590-0.806)	0.235 (0.159-0.392)
MTI	0.684 (0.566-0.790)	0.186 (0.134-0.306)
Clinical model	0.659 (0.554-0.764)	0.181 (0.123-0.321)
Fractional shortening	0.658 (0.557-0.747)	0.156 (0.120-0.237)
Published E′ timing/velocity metric	0.648 (0.560-0.732)	0.141 (0.114-0.211)
LVEF	0.645 (0.536-0.745)	0.156 (0.118-0.254)
Age + LVEF	0.633 (0.518-0.739)	0.152 (0.114-0.244)
E′ amplitude	0.613 (0.505-0.704)	0.138 (0.104-0.228)
Expanded TDI landmarks	0.601 (0.511-0.693)	0.130 (0.104-0.203)
E′ timing	0.592 (0.465-0.707)	0.147 (0.107-0.233)
S′ + E′ + A′ landmarks	0.592 (0.483-0.696)	0.135 (0.103-0.231)
FPCA-5*	0.619	0.138
Waveform ridge-25*	0.533	0.106
Clinical model: age, LVEF, anthracycline therapy, and anti-HER2 therapy. *FPCA-5 and waveform ridge-25 used 20 repeated 5-fold cross-validation because of computational cost; all other benchmark models used 100 repeats.
Panel E. Exploratory 1-year added-value and preserved-LVEF analyses
Analysis	Result
Known 1-year status	172 women / 16 events
LVEF AUROC at 1 y	0.652
MTI + LVEF AUROC at 1 y	0.698
Continuous NRI, LVEF -> MTI + LVEF	+0.378 (95% CI, -0.122 to 0.869)
IDI, LVEF -> MTI + LVEF	+0.015 (95% CI, -0.003 to 0.034)
Clinical model AUROC at 1 y	0.618
Clinical + MTI AUROC at 1 y	0.653
Clinical model AUPRC -> Clinical + MTI AUPRC	0.144 -> 0.175
Net benefit at 10% 1-y threshold, LVEF	0.014
Net benefit at 10% 1-y threshold, MTI + LVEF	0.032
Baseline LVEF >=50%, highest MTI tertile	15/83 events (18.1%)
Baseline LVEF >=50%, lower two MTI tertiles	11/178 events (6.2%)
Highest vs lower MTI tertiles, LVEF-adjusted HR	3.38 (95% CI, 1.55-7.38); P=0.002
NRI, IDI, decision-curve, and tertile analyses are exploratory. There are no validated absolute CTRCD risk categories or MTI cutoffs in this cohort.
Abbreviations: AUPRC, area under the precision-recall curve; AUROC, area under the receiver-operating-characteristic curve; CTRCD, cancer therapy-related cardiac dysfunction; FPCA, functional principal component analysis; HR, hazard ratio; IDI, integrated discrimination improvement; LVEF, left ventricular ejection fraction; MTI, Mechanical Transition Instability; NRI, net reclassification improvement; TDI, Tissue Doppler imaging.

