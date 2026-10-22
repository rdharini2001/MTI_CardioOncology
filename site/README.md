# MTI Interactive Research Companion

This directory is a dependency-free static research website for the Mechanical Transition Instability (MTI) cardio-oncology study.

## Run locally

```bash
python -m http.server 8000 --directory site
```

Open `http://localhost:8000`.

## Deploy on Vercel

After pushing the repository to GitHub, import the repository into Vercel and set the project root directory to `site/`. No environment variables or build command are required. `vercel.json` is included.

The site is a research companion, not a clinical calculator. The risk-enrichment slider uses cross-fitted study predictions and does not define a validated treatment or surveillance threshold.
