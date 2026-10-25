# MTI research website

This directory is the static GitHub Pages research companion for the Mechanical Transition Instability study.

## Local preview

From the repository root:

```bash
python -m http.server 8000 --directory docs
```

Open `http://localhost:8000`.

## Deployment

The repository includes `.github/workflows/pages.yml`, which deploys this directory with GitHub Actions. In GitHub, open **Settings → Pages** and select **GitHub Actions** as the source if it is not already selected.

The website uses aggregate study data from `docs/assets/study-data.json`. Patient-level BC_cardiotox data are not included in the site.
