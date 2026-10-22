#!/usr/bin/env bash
set -euo pipefail
python src/run_analysis.py
python src/nature_figures_extended.py --bootstrap 5000 --permutations 5000
python src/build_site_data.py
