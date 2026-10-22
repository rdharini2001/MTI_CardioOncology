# Mechanical Transition Instability in Cardio-Oncology

**Pretreatment temporal myocardial mechanics as a candidate risk-enrichment signal for cancer therapy-related cardiac dysfunction (CTRCD).**

This repository contains the complete reproducible analysis for a Journal of Cardiac Failure Brief Report evaluating **Mechanical Transition Instability (MTI)** from pretreatment Tissue Doppler Imaging (TDI) in the public `BC_cardiotox` cohort.

## Why MTI?

LVEF quantifies the net magnitude of ventricular systolic performance. MTI asks a different question: **how normally does myocardial motion resolve from systole into early relaxation before therapy begins?**

MTI is the equal-weight mean of four training-standardized features:

1. transition slope deficit;
2. transition curvature deficit;
3. early-relaxation latency; and
4. E′ attenuation.

Higher MTI indicates a slower, blunted, delayed, and lower-amplitude systolic-to-relaxation transition.

## Headline results

- TDI analysis cohort: **270 women, 27 CTRCD events**.
- MTI cross-fitted AUROC: **0.695 (95% CI 0.590–0.794)**.
- MTI + LVEF AUROC: **0.712 (95% CI 0.608–0.807)**.
- MTI adjusted for baseline LVEF: **HR 2.20 per 1 SD (95% CI 1.36–3.55)**.
- Preserved baseline LVEF ≥50%: **HR 2.22 (95% CI 1.41–3.51)**.
- Highest-risk 20% by cross-fitted MTI + LVEF contained **13/27 events (48.1%)**.
- Incremental discrimination over LVEF was not statistically established; paired uncertainty crossed zero.

## Interactive research companion

The `site/` directory contains a dependency-free interactive research website inspired by instrument-style scientific interfaces. It includes:

- an interactive pretreatment TDI waveform explorer;
- clickable MTI component anatomy;
- an MTI-decile phenotype explorer;
- a surveillance-efficiency slider based on held-out risk rankings;
- model performance and adjusted hazard views;
- methodology and reproducibility tabs; and
- an **SCI → MTI conceptual lineage** page.

Run locally:

```bash
cd site
python -m http.server 8000
```

Then open `http://localhost:8000`.

Once this repository is pushed to GitHub, it is ready for static deployment on Vercel. A `vercel.json` file is included.

## Repository structure

```text
MTI-CardioOncology/
├── data/
│   └── README.md                       # public dataset download instructions
├── legacy/
│   └── README.md                       # earlier SCI concept and provenance note
├── manuscript/
│   ├── JCF_MTI_Brief_Report.md
│   ├── JCF_MTI_Brief_Report_BLINDED.docx
│   ├── Table1_Multipanel.md
│   └── Figure_Legends.md
├── notebooks/
│   └── 01_MTI_JCF_analysis.ipynb
├── results/
│   ├── figures/
│   ├── tables/
│   ├── RESULTS_SUMMARY.md
│   └── NATURE_STYLE_ANALYSIS_NOTES.md
├── site/
│   ├── index.html
│   ├── styles.css
│   ├── app.js
│   ├── assets/study-data.json
│   └── vercel.json
├── src/
│   ├── run_analysis.py
│   ├── nature_figures_extended.py
│   ├── nature_theme.py
│   └── build_site_data.py
├── requirements.txt
└── reproduce.sh
```

## Data

The project **does not redistribute patient-level data**.

Download the public dataset from:

- **BC_cardiotox**: https://doi.org/10.6084/m9.figshare.22650748.v4
- Data descriptor: https://doi.org/10.1038/s41597-023-02419-1

Place `BC_cardiotox_clinical_and_functional_variables.csv` under `data/` before running the analysis. The public source dataset is released under CC BY 4.0; cite the original data creators and comply with their license.

## Reproduce

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

python src/run_analysis.py
python src/nature_figures_extended.py --bootstrap 5000 --permutations 5000
python src/build_site_data.py
```

The analysis uses 100 repeated stratified 5-fold cross-validation partitions. MTI standardization and model fitting occur within training folds. Figure-generation code writes the numeric tables used for each inferential visualization before plotting.

## Manuscript and figures

The JCF package is organized around the Brief Report 3-item limit:

- **Figure 1:** pretreatment TDI phenotype, MTI construction, waveform atlas, and functional inference;
- **Figure 2:** discrimination, survival association, adjusted HRs, and surveillance risk concentration;
- **Table 1:** four-panel comprehensive cohort, MTI, performance, association, and robustness table.

TIFF versions of the two primary figures are provided under `results/figures/`.

## Relationship to the earlier SCI project

The earlier project proposed **Systolic Completion Instability (SCI)** from regional LV cavity motion in apical 4-chamber cine echocardiography. The scientific premise was that abnormal temporal mechanics may precede overt ventricular dysfunction.

The EchoRisk organizers subsequently clarified that the challenge-distributed clinical dataset could not be used for this independent journal analysis. **No EchoRisk/CARDIOCARE patient-level data or results are used in the present MTI study.**

The current work preserves the mechanistic premise in a different echocardiographic representation that is fully reproducible from public data: TDI myocardial velocity rather than A4C cavity geometry.

Legacy project link: https://github.com/rdharini2001/EchoRisk

## Statistical interpretation guardrails

This is an internally validated exploratory study with 27 events. The manuscript intentionally avoids claims that MTI:

- outperforms LVEF;
- defines a clinical cutoff;
- replaces GLS or biomarkers; or
- is ready to guide surveillance without external validation.

The strongest supported interpretation is that MTI is **independently associated with CTRCD and may provide complementary baseline risk-enrichment information**.

## License

Code in this repository is released under the MIT License. The BC_cardiotox dataset is not part of this repository and remains governed by its original CC BY 4.0 terms.
