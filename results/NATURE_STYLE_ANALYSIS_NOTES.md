# Nature-style figure and extended analysis notes

## Recommended JCF figure set

### Figure 1. Pretreatment mechanical-transition phenotype
Use `results/nature_figures/Figure1_JCF_MTI_mechanism.png` (or TIFF).

Panels are generated directly from the public BC_cardiotox data and show:
1. Outcome-associated pretreatment TDI waveforms.
2. MTI construction from four interpretable waveform features.
3. A patient-level waveform atlas ordered by MTI.
4. A phase-resolved event-minus-nonevent waveform difference with bootstrap uncertainty and functional permutation testing.

**Interpretation guardrail:** the global unadjusted functional waveform test was P=0.047. After residualizing the waveform for baseline LVEF, P=0.056. The strongest contiguous phase cluster was approximately phase 0.499–0.600, but did not survive cluster-wise family-error correction (FWER P=0.094). Therefore, describe the waveform morphology as a distributed transition phenotype rather than claiming a statistically isolated E-prime time window.

### Figure 2. Clinical association, discrimination, and surveillance enrichment
Use `results/nature_figures/Figure2_JCF_MTI_clinical.png` (or TIFF).

Panels show:
1. Cross-fitted ROC curves.
2. CTRCD-free survival by MTI tertile.
3. Cox-model hazard ratios per 1-SD higher MTI.
4. Event-capture / risk-concentration curves for surveillance triage.

At the top 20% of the cross-fitted MTI+LVEF risk ranking, 13 of 27 future CTRCD events were captured (48.1%; bootstrap 95% CI approximately 25.9%–63.0%). The PPV in that selected group was 24.1% (95% CI approximately 13.0%–31.5%), compared with 10.0% event prevalence in the full TDI cohort. This is an exploratory resource-allocation analysis, not a validated clinical threshold.

## New statistically rigorous analyses

### 1. Continuous dose-response rather than threshold hunting
A natural-cubic-spline Cox model adjusted for baseline LVEF showed essentially no evidence of nonlinearity (`P_nonlin=0.999`). The observed relation is consistent with a progressive risk gradient rather than a discrete MTI threshold. This is preferable to choosing a post-hoc cutoff.

### 2. MTI-LVEF complementarity
A continuous MTI-by-LVEF interaction term was not significant (`P=0.891`). This argues against claiming synergy. MTI and LVEF are better described as providing complementary additive information. The point estimate for cross-fitted MTI+LVEF discrimination is higher than LVEF alone, but paired bootstrap uncertainty crosses zero.

### 3. Event concentration for targeted surveillance
The cross-fitted risk ranking was converted to an event-capture curve. This asks a clinically tangible question: if only X% of patients could receive intensified monitoring, what fraction of future CTRCD events would fall within that group? The ranking is based on held-out predictions, which avoids using in-sample fitted risk for the triage analysis.

### 4. Outcome-agnostic functional PCA
Functional PCA was performed without outcome labels. PC1 explained 33.2% of waveform variance and was associated with CTRCD after LVEF adjustment (HR 1.71 per 1 SD, nominal P=0.016), but did not meet BH-FDR significance across the first five components (q=0.079). Treat this as exploratory orthogonal support, not confirmatory evidence.

### 5. TDI availability selection audit
Compared with women without complete TDI waveforms, the TDI subset was older (SMD 0.45), had shorter follow-up (SMD -0.41), had less anti-HER2 exposure (SMD -0.34), and more hypertension (SMD 0.22). The CTRCD event rate was essentially identical (SMD -0.01). These differences should remain in the limitations because waveform availability is not missing completely at random.

## Manuscript-safe wording

Suggested Results sentence:

> In an exploratory surveillance-efficiency analysis using cross-fitted risk rankings, the highest-risk 20% by MTI plus LVEF contained 13 of 27 subsequent CTRCD events (48.1%), with a positive predictive value of 24.1% compared with 10.0% event prevalence in the full TDI cohort. A spline analysis showed no evidence of departure from a linear MTI-risk association (P for nonlinearity=0.999).

Suggested Discussion sentence:

> The continuous risk gradient and event-concentration analysis support MTI as a candidate risk-enrichment measure rather than a diagnostic cutoff; these analyses are exploratory and require external validation before use for surveillance allocation.

## Code

- `src/nature_theme.py`: publication styling and restrained Nature-family palette.
- `src/nature_figures_extended.py`: extended statistics, numeric output tables, all panel plots, and composite figures.
- All plotted values are saved under `results/tables/` before visualization.
- Individual panels are saved as both PNG and SVG under `results/nature_figures/panels/`.
- Composite JCF figures are saved as 300-dpi PNG and LZW-compressed TIFF.

Run from the repository root:

```bash
PYTHONPATH=src python src/nature_figures_extended.py --bootstrap 5000 --permutations 5000
```

The cross-fitted prediction file must first be produced by `src/run_analysis.py`.
