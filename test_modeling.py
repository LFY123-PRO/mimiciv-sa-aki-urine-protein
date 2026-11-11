"""Synthetic checks only; never opens clinical data."""
from pathlib import Path
import json
import hashlib
import time
import numpy as np
import pandas as pd
from scipy.special import expit
import statsmodels.api as sm
from modeling import (PROTEIN_LEVELS, Preprocessor, ModelError, fit_comparison,
                      evaluate_comparison, bootstrap_attempt, aggregate_bootstrap,
                      fit_logistic, performance, auc)

OUT = Path(__file__).resolve().parent
rng = np.random.default_rng(998133)
n = 1800
d = pd.DataFrame({
    'age': rng.normal(65, 15, n), 'gender': rng.choice(['F', 'M'], n),
    'prior_ckd_recorded': rng.choice([False, True], n, p=[.8, .2]),
    'nonrenal_sofa': rng.integers(0, 16, n),
    'nonrenal_sofa_missing_components': rng.integers(0, 4, n),
    'max_cruo_stage_to_L': rng.choice([1, 2, 3], n),
    'log2_last_scr': rng.normal(0, 1, n), 'log2_scr_change': rng.normal(0, .5, n),
    'uo_rate_6hr_valid': rng.exponential(.5, n),
    'uo_record_hour_coverage': rng.choice(np.arange(7)/6, n),
    'rrt_state_L': rng.choice(['no_supported_evidence', 'other_evidence_no_main', 'supported_main'], n),
    'protein_band': rng.choice(PROTEIN_LEVELS, n),
    'urine_specific_gravity': rng.normal(1.025, .01, n),
    'log1p_urine_rbc_hpf': rng.normal(2, 1, n),
    'log1p_urine_wbc_hpf': rng.normal(2, 1, n)})
eta = -3 + .018*(d.age-65) + .13*d.nonrenal_sofa + .2*d.log2_last_scr + \
      .35*(d.protein_band != PROTEIN_LEVELS[0]) + .6*(d.gender == 'M')
d['death28_lower'] = rng.binomial(1, expit(eta))
d['death28_upper'] = d.death28_lower
for c in ['log2_last_scr', 'log2_scr_change', 'uo_rate_6hr_valid',
          'urine_specific_gravity', 'log1p_urine_rbc_hpf', 'log1p_urine_wbc_hpf']:
    d.loc[rng.choice(n, n//8, replace=False), c] = np.nan
checks = []
def check(name, truth, note=''):
    checks.append({'check': name, 'passed': bool(truth), 'note': note})
    if not truth:
        raise AssertionError(name)

tic = time.perf_counter()
p = fit_comparison(d)
pe = fit_comparison(d, enhanced=True)
check('locked_candidate_df', (p.fit0.diagnostics['predictor_df'], p.fit1.diagnostics['predictor_df'],
                             pe.fit0.diagnostics['predictor_df'], pe.fit1.diagnostics['predictor_df']) == (16, 21, 22, 27))
check('rank_includes_intercept', p.fit1.diagnostics['rank_including_intercept'] == 22)
check('training_calibration_expected', abs(p.train_metrics['M1_calibration_intercept']) < 1e-7 and
      abs(p.train_metrics['M1_calibration_slope']-1) < 1e-7)
check('paired_auc_direction', np.isclose(p.train_metrics['delta_auc'],
      p.train_metrics['M1_auc']-p.train_metrics['M0_auc']))
check('paired_brier_direction', np.isclose(p.train_metrics['brier_gain'],
      p.train_metrics['M0_brier']-p.train_metrics['M1_brier']))
check('training_versus_validation_equal_same_data', p.metrics_flat() == evaluate_comparison(p, d))
check('five_OR_contrasts', len(p.or_table()) == 5 and all(z['status'] == 'estimated' for z in p.or_table()))
check('six_risks_not_a_multiclass_simplex', len(p.risks) == 6 and not np.isclose(sum(r['risk'] for r in p.risks), 1))

# Independent statsmodels IRLS on the original (not scaled) matrix.
x = p.preprocess1.transform(d)
ref = sm.GLM(d.death28_lower, x, family=sm.families.Binomial()).fit(tol=1e-11, maxiter=200)
check('MLE_predictions_vs_independent_IRLS', np.max(np.abs(ref.predict()-p.fit1.predict(x))) < 1e-8)
check('MLE_covariance_vs_independent_IRLS', np.allclose(ref.cov_params(), p.fit1.covariance, atol=1e-6, rtol=1e-6))
check('AUC_tie_ranks', np.isclose(auc(np.array([0, 1, 0, 1]), np.array([.2, .2, .8, .8])), .5))

validation = d.copy()
validation.loc[:, 'log2_last_scr'] = np.nan
original_median = p.preprocess0.medians['log2_last_scr']
transformed = p.preprocess0.transform(validation)
col = p.preprocess0.names.index('log2_last_scr')
check('validation_uses_training_median', np.allclose(transformed[:, col], original_median))
check('validation_never_changes_training_fill', p.preprocess0.medians['log2_last_scr'] == original_median)
bad = d.copy()
bad.loc[0, 'protein_band'] = np.nan
try:
    evaluate_comparison(p, bad)
    unknown_stopped = False
except ModelError as exc:
    unknown_stopped = str(exc).startswith('outside_training_support:protein_band')
check('missing_protein_not_negative', unknown_stopped)

no_missing_train = d.copy()
no_missing_train['log2_last_scr'] = no_missing_train.log2_last_scr.fillna(0)
no_missing_pp = Preprocessor(protein=True).fit(no_missing_train)
new_missing_validation = no_missing_train.copy()
new_missing_validation.loc[0, 'log2_last_scr'] = np.nan
check('new_validation_missing_preserves_training_schema',
      no_missing_pp.transform(no_missing_train).shape == no_missing_pp.transform(new_missing_validation).shape and
      'log2_last_scr__missing' not in no_missing_pp.names)

all_missing_train = d.copy()
all_missing_train['log2_last_scr'] = np.nan
all_missing_pp = Preprocessor(protein=True).fit(all_missing_train)
check('all_missing_training_recorded_and_constants_removed',
      'log2_last_scr' in all_missing_pp.all_missing and
      'log2_last_scr' not in all_missing_pp.names and
      'log2_last_scr__missing' in all_missing_pp.dropped_constants)

small_train = d.loc[d.protein_band != PROTEIN_LEVELS[-1]].copy()
small_pair = fit_comparison(small_train)
try:
    evaluate_comparison(small_pair, d)
    unsupported_stopped = False
except ModelError as exc:
    unsupported_stopped = str(exc).startswith('outside_training_support:protein_band')
check('unseen_protein_category_not_zero_or_reference', unsupported_stopped)

avg = np.mean(expit(x@p.fit1.beta))
atmeans = expit(x.mean(axis=0)@p.fit1.beta)
altered = d.assign(protein_band=PROTEIN_LEVELS[0])
ref_average = float(p.fit1.predict(p.preprocess1.transform(altered)).mean())
check('standardization_uses_each_person_then_average', np.isclose(p.risks[0]['risk'], ref_average))
check('empirical_average_differs_from_at_means', abs(avg-atmeans) > 1e-4)
check('reference_RD_is_zero', p.risks[0]['risk_difference'] == 0)

bt = []
for i in range(3):
    bt.append(bootstrap_attempt(d, 'death28_lower', rng.integers(0, n, n)))
check('three_actual_synthetic_bootstraps', all(a['status'] == 'success' for a in bt))
agg = aggregate_bootstrap(p, bt, requested_attempts=3)
check('all_attempts_aggregate', agg['complete'] and agg['CI_published'])
row = next(z for z in agg['metrics'] if z['metric'] == 'delta_auc')
q = np.quantile([a['train']['delta_auc'] for a in bt], [.025, .975])
optim = np.mean([a['train']['delta_auc']-a['evaluation']['delta_auc'] for a in bt])
check('paired_location_shift_formula', np.isclose(row['CI_lower'], q[0]-optim) and
      np.isclose(row['corrected'], p.train_metrics['delta_auc']-optim))
check('RD_CI_not_optimism_shifted', np.isclose(agg['risk_differences'][1]['RD_CI_lower'],
      np.quantile([a['risks'][1]['risk_difference'] for a in bt], .025)))
check('fixed_attempts_no_success_subset_CI', not aggregate_bootstrap(p, bt+[{'status':'failure','reason':'synthetic'}],4)['CI_published'])

try:
    fit_logistic(np.column_stack([np.ones(30), np.arange(30)]), np.r_[np.zeros(15), np.ones(15)])
    separated_stopped = False
except ModelError as exc:
    separated_stopped = 'separation' in str(exc)
check('perfect_separation_stops', separated_stopped)
try:
    fit_logistic(np.column_stack([np.ones(40), np.arange(40), np.arange(40)]), np.tile([0, 1], 20))
    rank_stopped = False
except ModelError as exc:
    rank_stopped = str(exc) == 'rank_deficient_design'
check('rank_deficiency_stops_without_P_selection', rank_stopped)

report = {'data': 'synthetic only; no clinical records opened', 'status': 'PASS',
          'checks': checks, 'elapsed_seconds': time.perf_counter()-tic,
          'source_sha256': hashlib.sha256((OUT/'modeling.py').read_bytes()).hexdigest(),
          'note': 'Same-model engineering checks; not P04 or independent clinical replication.'}
(OUT/'SYNTHETIC_MODEL_QA.json').write_text(json.dumps(report, ensure_ascii=False, indent=2)+'\n', encoding='utf-8')
pd.DataFrame(checks).to_csv(OUT/'SYNTHETIC_MODEL_QA.csv', index=False)
print(json.dumps({'status': 'PASS', 'checks': len(checks), 'elapsed_seconds': report['elapsed_seconds'],
                  'source_sha256': report['source_sha256']}))
