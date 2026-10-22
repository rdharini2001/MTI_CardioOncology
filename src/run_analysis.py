from __future__ import annotations

import json, math, warnings
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.signal import savgol_filter
from scipy.stats import mannwhitneyu, spearmanr
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score, average_precision_score, roc_curve
from statsmodels.duration.hazard_regression import PHReg

warnings.filterwarnings('ignore')
SEED = 20261022
N_REPEATS = 100
N_BOOT = 2000
SMOOTH_WINDOW = 31
POLYORDER = 3

ROOT = Path(__file__).resolve().parents[1]
RAW_DEFAULT = ROOT / 'data' / 'BC_cardiotox_clinical_and_functional_variables.csv'
RESULTS = ROOT / 'results'
TABLES = RESULTS / 'tables'
TABLES.mkdir(parents=True, exist_ok=True)

COMPONENTS = [
    'transition_slope_magnitude',
    'transition_curvature',
    'early_relaxation_latency',
    'e_prime_magnitude',
]
# Higher oriented value = greater mechanical transition abnormality
ORIENTATION = {
    'transition_slope_magnitude': -1.0,
    'transition_curvature': -1.0,
    'early_relaxation_latency': +1.0,
    'e_prime_magnitude': -1.0,
}
COMPONENT_LABELS = {
    'transition_slope_magnitude': 'Transition slope deficit',
    'transition_curvature': 'Transition curvature deficit',
    'early_relaxation_latency': 'Early-relaxation latency',
    'e_prime_magnitude': "E-prime attenuation",
}


def read_data(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep=';', decimal=',')


def extract_tdi_features(wave: np.ndarray, smooth_window: int = SMOOTH_WINDOW) -> dict:
    wave = np.asarray(wave, float)
    t = np.linspace(0.0, 1.0, wave.size)
    s = savgol_filter(wave, smooth_window, POLYORDER)
    dt = t[1] - t[0]

    # S-wave peak: positive systolic velocity occurs early in the normalized cycle.
    idx_s_window = np.where((t >= 0.02) & (t <= 0.40))[0]
    i_s = idx_s_window[np.argmax(s[idx_s_window])]

    # Systolic termination: first positive-to-negative zero crossing after S peak.
    crossings = np.where((s[i_s:-1] >= 0) & (s[i_s + 1:] < 0))[0]
    if len(crossings):
        i_zero = i_s + crossings[0] + 1
    else:
        fallback = np.where((t >= 0.35) & (t <= 0.55))[0]
        i_zero = fallback[np.argmin(np.abs(s[fallback]))]

    # E-prime is the first major negative relaxation peak after systolic termination.
    lo_e = max(i_zero + 1, np.searchsorted(t, 0.45))
    hi_e = np.searchsorted(t, 0.78)
    i_e = lo_e + np.argmin(s[lo_e:hi_e + 1])

    # Post-E-prime saddle before late diastolic A-prime motion.
    hi_saddle = np.searchsorted(t, 0.86)
    i_saddle = i_e + np.argmax(s[i_e:hi_saddle + 1]) if i_e < hi_saddle else i_e

    # A-prime for descriptive landmarking only.
    lo_a = max(i_saddle + 1, np.searchsorted(t, 0.75))
    hi_a = np.searchsorted(t, 0.99)
    i_a = lo_a + np.argmin(s[lo_a:hi_a + 1]) if lo_a <= hi_a else wave.size - 2

    d1 = savgol_filter(s, smooth_window, POLYORDER, deriv=1, delta=dt)
    d2 = savgol_filter(s, smooth_window, POLYORDER, deriv=2, delta=dt)

    slope_mag = -float(np.min(d1[i_s:i_zero + 1]))
    curvature = float(np.trapezoid(np.abs(d2[i_s:i_saddle + 1]), t[i_s:i_saddle + 1]))
    latency = float(t[i_e] - t[i_zero])
    e_mag = float(abs(s[i_e]))

    # Secondary interpretable morphology measures.
    area_total = float(np.trapezoid(np.maximum(s[:i_zero + 1], 0), t[:i_zero + 1]))
    area_tail = float(np.trapezoid(np.maximum(s[i_s:i_zero + 1], 0), t[i_s:i_zero + 1]))
    systolic_persistence_fraction = area_tail / max(area_total, 1e-12)
    systolic_resolution_time = float(t[i_zero] - t[i_s])

    half_level = 0.5 * s[i_e]
    left = i_e
    while left > i_zero and s[left] <= half_level:
        left -= 1
    right = i_e
    while right < i_saddle and s[right] <= half_level:
        right += 1
    e_dispersion = float(t[right] - t[left])

    e_area = -float(np.trapezoid(np.minimum(s[i_zero:i_saddle + 1], 0), t[i_zero:i_saddle + 1]))
    e_equiv_width = e_area / max(e_mag, 1e-12)

    # Published Piñeiro-Lamas et al. comparator.
    v_e_prime = float(t[i_e] / max(e_mag, 1e-12))

    return {
        's_peak_time': float(t[i_s]),
        's_peak_velocity': float(s[i_s]),
        'systolic_zero_time': float(t[i_zero]),
        'e_prime_time': float(t[i_e]),
        'e_prime_velocity': float(s[i_e]),
        'e_prime_magnitude': e_mag,
        'a_prime_time': float(t[i_a]),
        'a_prime_velocity': float(s[i_a]),
        'transition_slope_magnitude': slope_mag,
        'transition_curvature': curvature,
        'early_relaxation_latency': latency,
        'early_relaxation_dispersion': e_dispersion,
        'early_relaxation_equiv_width': e_equiv_width,
        'systolic_persistence_fraction': systolic_persistence_fraction,
        'systolic_resolution_time': systolic_resolution_time,
        'v_e_prime': v_e_prime,
    }


def orient_and_standardize(df: pd.DataFrame, reference: pd.DataFrame | None = None) -> pd.DataFrame:
    ref = df if reference is None else reference
    out = pd.DataFrame(index=df.index)
    for c in COMPONENTS:
        mu = float(ref[c].mean())
        sd = float(ref[c].std(ddof=0)) or 1.0
        out[c] = ORIENTATION[c] * (df[c].astype(float) - mu) / sd
    return out


def make_mti(df: pd.DataFrame, reference: pd.DataFrame | None = None) -> pd.Series:
    z = orient_and_standardize(df, reference)
    return z.mean(axis=1)


def harrell_c(time, event, risk) -> float:
    time = np.asarray(time, float)
    event = np.asarray(event, int)
    risk = np.asarray(risk, float)
    total = concordant = tied = 0
    for i in range(len(time)):
        if event[i] != 1 or not np.isfinite(risk[i]):
            continue
        js = np.where((time > time[i]) & np.isfinite(risk))[0]
        if len(js) == 0:
            continue
        total += len(js)
        concordant += int(np.sum(risk[i] > risk[js]))
        tied += int(np.sum(risk[i] == risk[js]))
    return float((concordant + 0.5 * tied) / total) if total else np.nan


def one_year_labels(df: pd.DataFrame, horizon: int = 365):
    event = df['CTRCD'].to_numpy(int)
    time = df['time'].to_numpy(float)
    case = (event == 1) & (time <= horizon)
    control = time >= horizon
    evaluable = case | control
    y = case[evaluable].astype(int)
    return evaluable, y


def metrics(df: pd.DataFrame, risk: np.ndarray) -> dict:
    y = df['CTRCD'].to_numpy(int)
    time = df['time'].to_numpy(float)
    risk = np.asarray(risk, float)
    ev1, y1 = one_year_labels(df)
    return {
        'AUROC_ever': roc_auc_score(y, risk),
        'AUPRC_ever': average_precision_score(y, risk),
        'C_index': harrell_c(time, y, risk),
        'AUROC_1y': roc_auc_score(y1, risk[ev1]) if len(np.unique(y1)) == 2 else np.nan,
        'n_1y_evaluable': int(ev1.sum()),
        'events_1y': int(y1.sum()),
    }


def fit_cox_predict(train: pd.DataFrame, test: pd.DataFrame, predictors: list[str]) -> np.ndarray:
    tr = train.copy()
    te = test.copy()
    Xtr = tr[predictors].astype(float).copy()
    Xte = te[predictors].astype(float).copy()
    med = Xtr.median(axis=0)
    Xtr = Xtr.fillna(med)
    Xte = Xte.fillna(med)
    mu = Xtr.mean(axis=0)
    sd = Xtr.std(axis=0, ddof=0).replace(0, 1)
    A = (Xtr - mu) / sd
    B = (Xte - mu) / sd
    try:
        fit = PHReg(tr['time'].astype(float), A, status=tr['CTRCD'].astype(int), ties='breslow').fit(disp=0)
        beta = np.asarray(fit.params, float)
        if not np.all(np.isfinite(beta)):
            raise FloatingPointError('nonfinite Cox coefficients')
    except Exception:
        # Small-ridge regularized Cox fallback for rare fold instability.
        fit = PHReg(tr['time'].astype(float), A, status=tr['CTRCD'].astype(int), ties='breslow').fit_regularized(alpha=0.01)
        beta = np.asarray(fit.params, float)
    return np.asarray(B, float) @ beta


def repeated_crossfit(df: pd.DataFrame, repeats: int = N_REPEATS) -> pd.DataFrame:
    n = len(df)
    accum = {m: np.zeros(n, float) for m in ['MTI', 'LVEF', 'vEprime', 'MTI+LVEF', 'Clinical', 'Clinical+MTI']}
    counts = np.zeros(n, int)
    y = df['CTRCD'].to_numpy(int)

    for r in range(repeats):
        cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=SEED + r)
        for tr_idx, te_idx in cv.split(np.zeros(n), y):
            tr = df.iloc[tr_idx].copy()
            te = df.iloc[te_idx].copy()
            mti_tr = make_mti(tr, reference=tr)
            mti_te = make_mti(te, reference=tr)
            tr['MTI_fold'] = mti_tr.to_numpy()
            te['MTI_fold'] = mti_te.to_numpy()

            preds = {
                'MTI': te['MTI_fold'].to_numpy(float),
                'LVEF': -te['LVEF'].to_numpy(float),
                'vEprime': te['v_e_prime'].to_numpy(float),
                'MTI+LVEF': fit_cox_predict(tr, te, ['MTI_fold', 'LVEF']),
                'Clinical': fit_cox_predict(tr, te, ['age', 'LVEF', 'AC', 'antiHER2']),
                'Clinical+MTI': fit_cox_predict(tr, te, ['age', 'LVEF', 'AC', 'antiHER2', 'MTI_fold']),
            }
            for m, p in preds.items():
                accum[m][te_idx] += p
            counts[te_idx] += 1

    out = df[['CTRCD', 'time']].copy().reset_index(drop=True)
    for m, a in accum.items():
        out[m] = a / counts
    return out


def stratified_bootstrap_indices(y: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    y = np.asarray(y, int)
    pos = np.where(y == 1)[0]
    neg = np.where(y == 0)[0]
    return np.r_[rng.choice(pos, len(pos), replace=True), rng.choice(neg, len(neg), replace=True)]


def bootstrap_performance(df: pd.DataFrame, pred: pd.DataFrame, B: int = N_BOOT):
    rng = np.random.default_rng(SEED + 99)
    models = [c for c in pred.columns if c not in ['CTRCD', 'time']]
    store = {m: {k: [] for k in ['AUROC_ever', 'AUPRC_ever', 'C_index', 'AUROC_1y']} for m in models}
    diffs = {m: [] for m in models if m != 'LVEF'}
    y = df['CTRCD'].to_numpy(int)

    for _ in range(B):
        idx = stratified_bootstrap_indices(y, rng)
        bdf = df.iloc[idx].reset_index(drop=True)
        for m in models:
            r = pred[m].to_numpy(float)[idx]
            mm = metrics(bdf, r)
            for k in store[m]:
                store[m][k].append(mm[k])
        # primary difference in C-index vs LVEF
        c_lvef = harrell_c(bdf.time, bdf.CTRCD, pred['LVEF'].to_numpy(float)[idx])
        for m in diffs:
            c_m = harrell_c(bdf.time, bdf.CTRCD, pred[m].to_numpy(float)[idx])
            diffs[m].append(c_m - c_lvef)

    rows = []
    for m in models:
        point = metrics(df, pred[m].to_numpy(float))
        row = {'model': m, **point}
        for k, vals in store[m].items():
            vals = np.asarray(vals, float)
            row[f'{k}_CI_low'] = float(np.nanquantile(vals, .025))
            row[f'{k}_CI_high'] = float(np.nanquantile(vals, .975))
        if m != 'LVEF':
            vals = np.asarray(diffs[m], float)
            row['delta_C_vs_LVEF'] = point['C_index'] - metrics(df, pred['LVEF'].to_numpy(float))['C_index']
            row['delta_C_CI_low'] = float(np.nanquantile(vals, .025))
            row['delta_C_CI_high'] = float(np.nanquantile(vals, .975))
        rows.append(row)
    return pd.DataFrame(rows)


def cox_table(df: pd.DataFrame, predictors: list[str], label: str, complete_case: bool = True) -> pd.DataFrame:
    cols = ['time', 'CTRCD'] + predictors
    d = df[cols].dropna().copy() if complete_case else df[cols].copy()
    X = d[predictors].astype(float)
    if not complete_case:
        X = X.fillna(X.median())
    mu = X.mean(); sd = X.std(ddof=0).replace(0, 1)
    Xs = (X - mu) / sd
    fit = PHReg(d['time'].astype(float), Xs, status=d['CTRCD'].astype(int), ties='breslow').fit(disp=0)
    rows = []
    for term, b, se, p in zip(predictors, fit.params, fit.bse, fit.pvalues):
        rows.append({
            'model': label,
            'term': term,
            'HR_per_1SD': float(np.exp(b)),
            'CI_low': float(np.exp(b - 1.96 * se)),
            'CI_high': float(np.exp(b + 1.96 * se)),
            'p_value': float(p),
            'n': len(d),
            'events': int(d['CTRCD'].sum()),
        })
    return pd.DataFrame(rows)


def bootstrap_component_auc(df: pd.DataFrame, feature: str, orientation: float, B: int = 1000):
    y = df['CTRCD'].to_numpy(int)
    x = orientation * df[feature].to_numpy(float)
    point = roc_auc_score(y, x)
    rng = np.random.default_rng(SEED + 123)
    vals = []
    for _ in range(B):
        idx = stratified_bootstrap_indices(y, rng)
        vals.append(roc_auc_score(y[idx], x[idx]))
    return point, float(np.quantile(vals, .025)), float(np.quantile(vals, .975))


def smoothing_sensitivity(raw_df: pd.DataFrame, tdi_cols: list[str]) -> pd.DataFrame:
    rows = []
    y = raw_df['CTRCD'].to_numpy(int)
    for window in [21, 31, 41, 51, 61]:
        feats = pd.DataFrame([extract_tdi_features(w, smooth_window=window) for w in raw_df[tdi_cols].to_numpy(float)])
        temp = pd.concat([raw_df.reset_index(drop=True), feats], axis=1)
        score = make_mti(temp)
        score = (score - score.mean()) / score.std(ddof=0)
        fit = PHReg(temp['time'].astype(float), score.to_numpy()[:, None], status=temp['CTRCD'].astype(int)).fit(disp=0)
        rows.append({
            'savgol_window': window,
            'AUROC_ever': roc_auc_score(y, score),
            'C_index': harrell_c(temp.time, temp.CTRCD, score),
            'HR_per_1SD': float(np.exp(fit.params[0])),
            'p_value': float(fit.pvalues[0]),
        })
    return pd.DataFrame(rows)


def main():
    raw_path = RAW_DEFAULT
    if not raw_path.exists():
        # Convenience fallback when run from supplied archive extraction.
        alt = Path('/mnt/data/bc_cardiotox_raw/BC_cardiotox_clinical_and_functional_variables.csv')
        if alt.exists():
            raw_path = alt
        else:
            raise FileNotFoundError('Place BC_cardiotox_clinical_and_functional_variables.csv in data/.')

    full = read_data(raw_path)
    tdi_cols = [c for c in full.columns if c.startswith('t ')]
    full['TDI_available'] = full[tdi_cols].notna().all(axis=1).astype(int)
    tdi = full.loc[full['TDI_available'].eq(1)].copy().reset_index(drop=True)

    feature_df = pd.DataFrame([extract_tdi_features(w) for w in tdi[tdi_cols].to_numpy(float)])
    analysis = pd.concat([tdi.drop(columns=tdi_cols).reset_index(drop=True), feature_df], axis=1)
    analysis['fractional_shortening'] = (analysis['LVDd'] - analysis['LVSd']) / analysis['LVDd'] * 100.0
    analysis['MTI_raw'] = make_mti(analysis)
    analysis['MTI'] = (analysis['MTI_raw'] - analysis['MTI_raw'].mean()) / analysis['MTI_raw'].std(ddof=0)

    # Cohort summary and provenance audit.
    cohort = pd.DataFrame([
        {'cohort': 'Full public clinical cohort', 'n': len(full), 'events': int(full.CTRCD.sum()), 'event_rate': float(full.CTRCD.mean()), 'median_followup_days': float(full.time.median())},
        {'cohort': 'TDI analysis cohort', 'n': len(analysis), 'events': int(analysis.CTRCD.sum()), 'event_rate': float(analysis.CTRCD.mean()), 'median_followup_days': float(analysis.time.median())},
        {'cohort': 'TDI cohort with baseline LVEF >=50%', 'n': int((analysis.LVEF >= 50).sum()), 'events': int(analysis.loc[analysis.LVEF >= 50, 'CTRCD'].sum()), 'event_rate': float(analysis.loc[analysis.LVEF >= 50, 'CTRCD'].mean()), 'median_followup_days': float(analysis.loc[analysis.LVEF >= 50, 'time'].median())},
    ])
    cohort.to_csv(TABLES / 'cohort_summary.csv', index=False)

    # Validate E-prime landmark extraction against the published analysis.
    validation = analysis.groupby('CTRCD').agg(
        n=('CTRCD', 'size'),
        e_prime_time_mean=('e_prime_time', 'mean'),
        e_prime_time_sd=('e_prime_time', 'std'),
        e_prime_magnitude_mean=('e_prime_magnitude', 'mean'),
        e_prime_magnitude_sd=('e_prime_magnitude', 'std'),
        v_e_prime_mean=('v_e_prime', 'mean'),
        v_e_prime_sd=('v_e_prime', 'std'),
    ).reset_index()
    validation.to_csv(TABLES / 'landmark_validation.csv', index=False)

    # Component descriptive statistics + bootstrap AUC.
    comp_rows = []
    for c in COMPONENTS:
        ev = analysis.loc[analysis.CTRCD.eq(1), c]
        ne = analysis.loc[analysis.CTRCD.eq(0), c]
        auc, lo, hi = bootstrap_component_auc(analysis, c, ORIENTATION[c])
        p = mannwhitneyu(ev, ne, alternative='two-sided').pvalue
        comp_rows.append({
            'component': COMPONENT_LABELS[c],
            'raw_feature': c,
            'event_median': float(ev.median()),
            'event_Q1': float(ev.quantile(.25)),
            'event_Q3': float(ev.quantile(.75)),
            'nonevent_median': float(ne.median()),
            'nonevent_Q1': float(ne.quantile(.25)),
            'nonevent_Q3': float(ne.quantile(.75)),
            'event_higher_oriented_AUROC': auc,
            'AUROC_CI_low': lo,
            'AUROC_CI_high': hi,
            'Mann_Whitney_p': float(p),
        })
    pd.DataFrame(comp_rows).to_csv(TABLES / 'mti_components.csv', index=False)

    # Primary Cox models.
    cox_parts = []
    for predictors, label, cc in [
        (['MTI'], 'MTI unadjusted', True),
        (['LVEF'], 'LVEF unadjusted', True),
        (['v_e_prime'], 'Published vEprime comparator', True),
        (['fractional_shortening'], 'Fractional shortening comparator', True),
        (['MTI', 'LVEF'], 'MTI + LVEF', True),
        (['MTI', 'age', 'LVEF'], 'MTI + age + LVEF', True),
        (['MTI', 'age', 'LVEF', 'AC', 'antiHER2'], 'MTI + age + LVEF + treatment', True),
    ]:
        cox_parts.append(cox_table(analysis, predictors, label, complete_case=cc))
    cox = pd.concat(cox_parts, ignore_index=True)
    cox.to_csv(TABLES / 'cox_models.csv', index=False)

    # Preserved baseline LVEF sensitivity.
    preserved = analysis.loc[analysis.LVEF >= 50].copy()
    cox_pres = cox_table(preserved, ['MTI', 'LVEF'], 'Preserved baseline LVEF >=50%', complete_case=True)
    cox_pres.to_csv(TABLES / 'cox_preserved_lvef_sensitivity.csv', index=False)

    # Repeated leakage-safe cross-fitting. Reuse an exact checkpoint when present.
    pred_path = RESULTS / 'crossfitted_predictions.csv'
    if pred_path.exists():
        pred = pd.read_csv(pred_path)
        if len(pred) != len(analysis):
            pred = repeated_crossfit(analysis, repeats=N_REPEATS)
            pred.insert(0, 'row_id', np.arange(len(pred)))
            pred.to_csv(pred_path, index=False)
    else:
        pred = repeated_crossfit(analysis, repeats=N_REPEATS)
        pred.insert(0, 'row_id', np.arange(len(pred)))
        pred.to_csv(pred_path, index=False)
    perf = bootstrap_performance(analysis, pred.drop(columns='row_id'), B=N_BOOT)
    perf.to_csv(TABLES / 'crossvalidated_performance.csv', index=False)

    # Paired operating characteristics using cross-fitted Youden threshold on aggregate predictions (secondary only).
    op_rows = []
    y = analysis.CTRCD.to_numpy(int)
    for model in ['MTI', 'LVEF', 'vEprime', 'MTI+LVEF', 'Clinical', 'Clinical+MTI']:
        risk = pred[model].to_numpy(float)
        fpr, tpr, th = roc_curve(y, risk)
        j = int(np.nanargmax(tpr - fpr))
        threshold = float(th[j])
        z = risk >= threshold
        tp = int(((z == 1) & (y == 1)).sum()); tn = int(((z == 0) & (y == 0)).sum())
        fp = int(((z == 1) & (y == 0)).sum()); fn = int(((z == 0) & (y == 1)).sum())
        op_rows.append({
            'model': model, 'threshold': threshold,
            'sensitivity': tp / max(tp + fn, 1), 'specificity': tn / max(tn + fp, 1),
            'PPV': tp / max(tp + fp, 1), 'NPV': tn / max(tn + fn, 1),
            'balanced_accuracy': 0.5 * (tp / max(tp + fn, 1) + tn / max(tn + fp, 1)),
        })
    pd.DataFrame(op_rows).to_csv(TABLES / 'operating_metrics_secondary.csv', index=False)

    # MTI tertiles for KM visualization and descriptive gradient.
    q1, q2 = analysis['MTI'].quantile([1/3, 2/3]).to_numpy()
    analysis['MTI_tertile'] = pd.cut(analysis['MTI'], [-np.inf, q1, q2, np.inf], labels=['Low', 'Middle', 'High'])
    tert = analysis.groupby('MTI_tertile', observed=True).agg(n=('CTRCD','size'), events=('CTRCD','sum'), event_rate=('CTRCD','mean'), median_followup=('time','median')).reset_index()
    tert.to_csv(TABLES / 'mti_tertiles.csv', index=False)

    # Component ablation: equal-weight 3-component scores.
    ab_rows = []
    full_auc = roc_auc_score(analysis.CTRCD, analysis.MTI)
    full_c = harrell_c(analysis.time, analysis.CTRCD, analysis.MTI)
    for omit in [None] + COMPONENTS:
        keep = [c for c in COMPONENTS if c != omit]
        z = orient_and_standardize(analysis)[keep].mean(axis=1)
        z = (z - z.mean()) / z.std(ddof=0)
        ab_rows.append({
            'formulation': 'Full MTI' if omit is None else f'Omit {COMPONENT_LABELS[omit]}',
            'AUROC_ever': roc_auc_score(analysis.CTRCD, z),
            'delta_AUROC_vs_full': roc_auc_score(analysis.CTRCD, z) - full_auc,
            'C_index': harrell_c(analysis.time, analysis.CTRCD, z),
            'delta_C_vs_full': harrell_c(analysis.time, analysis.CTRCD, z) - full_c,
        })
    pd.DataFrame(ab_rows).to_csv(TABLES / 'component_ablation.csv', index=False)

    # Smoothing-window sensitivity.
    smooth = smoothing_sensitivity(tdi, tdi_cols)
    smooth.to_csv(TABLES / 'smoothing_sensitivity.csv', index=False)

    # TDI-selection audit.
    audit_rows = []
    for c in ['age', 'LVEF', 'heart_rate', 'time']:
        g1 = full.loc[full.TDI_available.eq(1), c].dropna(); g0 = full.loc[full.TDI_available.eq(0), c].dropna()
        pooled = np.sqrt((g1.var(ddof=1) + g0.var(ddof=1)) / 2)
        smd = (g1.mean() - g0.mean()) / pooled if pooled else np.nan
        audit_rows.append({'variable': c, 'TDI_mean': g1.mean(), 'No_TDI_mean': g0.mean(), 'standardized_mean_difference': smd})
    for c in ['CTRCD', 'AC', 'antiHER2', 'HTA', 'DL', 'DM', 'RTprev']:
        p1 = full.loc[full.TDI_available.eq(1), c].mean(); p0 = full.loc[full.TDI_available.eq(0), c].mean()
        pbar = (p1 + p0) / 2
        denom = np.sqrt(pbar * (1-pbar)) if 0 < pbar < 1 else np.nan
        smd = (p1 - p0) / denom if denom else np.nan
        audit_rows.append({'variable': c, 'TDI_mean': p1, 'No_TDI_mean': p0, 'standardized_mean_difference': smd})
    pd.DataFrame(audit_rows).to_csv(TABLES / 'tdi_selection_audit.csv', index=False)

    # Save derived analysis data without waveform samples for reproducibility/audit.
    analysis.to_csv(RESULTS / 'derived_analysis_dataset.csv', index=False)

    # Compact JSON for website.
    perf_index = perf.set_index('model')
    mti_cox = cox[(cox.model == 'MTI + LVEF') & (cox.term == 'MTI')].iloc[0]
    mti_unadj = cox[(cox.model == 'MTI unadjusted') & (cox.term == 'MTI')].iloc[0]
    payload = {
        'n': int(len(analysis)), 'events': int(analysis.CTRCD.sum()),
        'event_rate': float(analysis.CTRCD.mean()),
        'median_followup_days': float(analysis.time.median()),
        'median_event_time_days': float(analysis.loc[analysis.CTRCD.eq(1), 'time'].median()),
        'MTI_AUROC': float(perf_index.loc['MTI','AUROC_ever']),
        'MTI_AUROC_low': float(perf_index.loc['MTI','AUROC_ever_CI_low']),
        'MTI_AUROC_high': float(perf_index.loc['MTI','AUROC_ever_CI_high']),
        'MTI_C': float(perf_index.loc['MTI','C_index']),
        'MTI_C_low': float(perf_index.loc['MTI','C_index_CI_low']),
        'MTI_C_high': float(perf_index.loc['MTI','C_index_CI_high']),
        'LVEF_AUROC': float(perf_index.loc['LVEF','AUROC_ever']),
        'LVEF_C': float(perf_index.loc['LVEF','C_index']),
        'vE_AUROC': float(perf_index.loc['vEprime','AUROC_ever']),
        'MTI_HR_unadjusted': float(mti_unadj.HR_per_1SD),
        'MTI_HR_adjusted_LVEF': float(mti_cox.HR_per_1SD),
        'MTI_HR_adjusted_LVEF_low': float(mti_cox.CI_low),
        'MTI_HR_adjusted_LVEF_high': float(mti_cox.CI_high),
        'MTI_HR_adjusted_LVEF_p': float(mti_cox.p_value),
        'high_tertile_events': int(tert.loc[tert.MTI_tertile.astype(str).eq('High'),'events'].iloc[0]),
        'low_tertile_events': int(tert.loc[tert.MTI_tertile.astype(str).eq('Low'),'events'].iloc[0]),
    }
    (RESULTS / 'website_summary.json').write_text(json.dumps(payload, indent=2))

    print('\n=== ANALYSIS COMPLETE ===')
    print(cohort.to_string(index=False))
    print('\nCross-fitted performance:')
    print(perf[['model','AUROC_ever','AUROC_ever_CI_low','AUROC_ever_CI_high','AUPRC_ever','C_index','C_index_CI_low','C_index_CI_high','AUROC_1y']].to_string(index=False))
    print('\nCox models:')
    print(cox.to_string(index=False))

if __name__ == '__main__':
    main()
