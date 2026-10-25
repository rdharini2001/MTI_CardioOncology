# MTI final research repository release

## Release goals

This release presents Mechanical Transition Instability (MTI) as an independent cardio-oncology study and synchronizes the repository with the final manuscript.

## Major changes

- Removed all prior-project lineage and legacy-study material so the repository presents MTI as an independent study.
- Replaced the older manuscript draft with the final manuscript source.
- Updated Table 1 to include the final comparator benchmark and exploratory 1-year added-value analyses.
- Updated figure legends to match the final study interpretation and code-generated figures.
- Rebuilt the root README in the same clean research-repository style as the Curriculum Learning project.
- Rebuilt the research website as a lively, responsive, dependency-free GitHub Pages companion.
- Added interactive TDI waveform, MTI component, phenotype-decile, model-comparison, Cox association, and surveillance-capacity views.
- Added preserved-LVEF and exploratory 1-year added-value sections.
- Added GitHub Pages deployment through `.github/workflows/pages.yml`.
- Added a repository consistency workflow and `tools/check_release.py`.
- Removed the Vercel-specific deployment configuration.
- Removed row-level `crossfitted_predictions.csv` from the distributable release so the repository matches the manuscript statement that patient-level data are not redistributed.
- The original BC_cardiotox CSV is not included.

## Headline values synchronized to the manuscript

- TDI cohort: 270 women, 27 CTRCD events.
- MTI AUROC: 0.695 (95% CI 0.590–0.794).
- MTI + LVEF AUROC: 0.712 (95% CI 0.608–0.807).
- MTI adjusted for LVEF: HR 2.20 (95% CI 1.36–3.55).
- Preserved baseline LVEF: HR 2.22 (95% CI 1.41–3.51).
- Highest vs lower MTI tertiles in preserved-LVEF cohort: HR 3.38 (95% CI 1.55–7.38).
- 1-year analysis: 172 women, 16 events; LVEF AUROC 0.652 → 0.698 after adding MTI.

## GitHub Pages

After pushing this release to `main`, open **Settings → Pages** and select **GitHub Actions** as the source. The included workflow deploys `docs/`.
