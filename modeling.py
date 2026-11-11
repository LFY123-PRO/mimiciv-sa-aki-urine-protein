"""Locked P03.3 logistic models; this module does not read patient files."""
from dataclasses import dataclass
import warnings
import numpy as np
import pandas as pd
from scipy.optimize import brentq, linprog
from scipy.special import expit
from scipy.stats import rankdata, chi2
import statsmodels.api as sm

SEED = 20261029
PROTEIN_LEVELS = ('0_negative_or_report0', '1_trace', '2_positive_lt30',
                  '3_30_lt100', '4_100_lt300', '5_ge300')
BASE_CONTINUOUS = ('age', 'nonrenal_sofa', 'nonrenal_sofa_missing_components',
                   'log2_last_scr', 'log2_scr_change', 'uo_rate_6hr_valid',
                   'uo_record_hour_coverage')
BASE_CATEGORICAL = ('gender', 'prior_ckd_recorded', 'max_cruo_stage_to_L', 'rrt_state_L')
ENHANCED = ('urine_specific_gravity', 'log1p_urine_rbc_hpf', 'log1p_urine_wbc_hpf')
UNKNOWN = 'Unknown'


class ModelError(RuntimeError):
    pass


def numeric(s):
    v = pd.to_numeric(s, errors='coerce').astype(float)
    if np.isinf(v.to_numpy()).any():
        raise ModelError('infinite_numeric_predictor')
    return v


def canonical(s, name):
    out = s.astype('string').str.strip()
    out = out.mask(s.isna() | out.eq(''), UNKNOWN)
    if name == 'gender':
        out = out.replace({'Female': 'F', 'female': 'F', 'Male': 'M', 'male': 'M'})
    elif name == 'prior_ckd_recorded':
        out = out.str.lower().replace({'false': '0', 'f': '0', '0.0': '0',
                                      'true': '1', 't': '1', '1.0': '1',
                                      'unknown': UNKNOWN})
    elif name == 'max_cruo_stage_to_L':
        out = out.replace({'1.0': '1', '2.0': '2', '3.0': '3'})
    return out.astype(str)


def normalized_data(data):
    d = data.copy(deep=False)
    if 'rrt_state_L' not in d and 'rrt_state_l' in d:
        d = d.rename(columns={'rrt_state_l': 'rrt_state_L'})
    return d


class Preprocessor:
    """All fill values, supported levels and constant drops are learned in training."""
    def __init__(self, enhanced=False, protein=False, extra_continuous=(), extra_categorical=()):
        self.continuous = BASE_CONTINUOUS + (ENHANCED if enhanced else ()) + tuple(extra_continuous)
        self.categorical = BASE_CATEGORICAL + tuple(extra_categorical) + (('protein_band',) if protein else ())
        self.protein = protein
        self.medians = {}
        self.levels = {}
        self.references = {}

    def fit(self, data):
        d = normalized_data(data)
        needed = set(self.continuous + self.categorical)
        if not needed.issubset(d.columns):
            raise ModelError('missing_model_columns:' + ','.join(sorted(needed-set(d.columns))))
        for c in self.continuous:
            x = numeric(d[c])
            self.medians[c] = float(x.median()) if x.notna().any() else 0.0
        preferred = {'gender': 'F', 'prior_ckd_recorded': '0',
                     'rrt_state_L': 'no_supported_evidence', 'protein_band': PROTEIN_LEVELS[0]}
        for c in self.categorical:
            x = canonical(d[c], c)
            observed = set(x)
            if c == 'protein_band':
                if not observed.issubset(PROTEIN_LEVELS):
                    raise ModelError('missing_or_unsupported_protein')
                if PROTEIN_LEVELS[0] not in observed:
                    raise ModelError('protein_reference_absent_in_training')
                levels = [v for v in PROTEIN_LEVELS if v in observed]
            elif c == 'max_cruo_stage_to_L':
                levels = sorted(observed, key=lambda z: (z == UNKNOWN, z))
            else:
                ref = preferred.get(c)
                levels = ([ref] if ref in observed else []) + sorted(observed-{ref})
            self.levels[c] = tuple(levels)
            self.references[c] = levels[0]
        raw, names = self._raw_matrix(d)
        self.all_names = tuple(names)
        # All-missing columns retain their source/indicator information in the log,
        # then constant columns are omitted because the intercept already spans them.
        keep = np.ptp(raw, axis=0) > 0
        self.keep = keep
        self.names = ('Intercept',) + tuple(np.asarray(names)[keep])
        self.dropped_constants = tuple(np.asarray(names)[~keep])
        self.all_missing = tuple(c for c in self.continuous if numeric(d[c]).isna().all())
        self.fitted = True
        return self

    def _raw_matrix(self, d):
        cols, names = [], []
        for c in self.continuous:
            x = numeric(d[c])
            value = x.fillna(self.medians[c]).to_numpy()
            if c == 'age':
                value = value/10.0
            cols.extend([value, x.isna().to_numpy(dtype=float)])
            names.extend([c + ('_per10' if c == 'age' else ''), c + '__missing'])
        for c in self.categorical:
            x = canonical(d[c], c)
            unsupported = set(x)-set(self.levels[c])
            if unsupported:
                raise ModelError('outside_training_support:' + c + ':' + ','.join(sorted(unsupported)))
            for level in self.levels[c][1:]:
                cols.append(x.eq(level).to_numpy(dtype=float))
                names.append(c + '[' + level + ']')
        return np.column_stack(cols).astype(float), names

    def transform(self, data):
        if not getattr(self, 'fitted', False):
            raise ModelError('preprocessor_not_fitted')
        raw, names = self._raw_matrix(normalized_data(data))
        if tuple(names) != self.all_names:
            raise ModelError('matrix_schema_changed')
        return np.column_stack([np.ones(len(data)), raw[:, self.keep]])

    def metadata(self):
        return {'continuous_medians': self.medians, 'categorical_levels': self.levels,
                'references': self.references, 'columns': self.names,
                'dropped_constant_columns': self.dropped_constants,
                'all_missing_continuous': self.all_missing,
                'age_scale': 'per 10 years', 'imputation': 'training only',
                'categorical_missing': 'Unknown (no duplicate missing dummy)'}


def _binary_outcome(y):
    y = np.asarray(y, dtype=float)
    if y.ndim != 1 or not np.isfinite(y).all() or not np.isin(y, [0, 1]).all():
        raise ModelError('outcome_not_complete_binary')
    if np.unique(y).size != 2:
        raise ModelError('outcome_has_one_class')
    return y


def separation_lp(x, y):
    """Detect a nonzero separating direction, including quasi-separation."""
    z = (2*y-1)[:, None]*x
    # With |beta|<=1, a positive sum of nonnegative signed margins proves
    # complete or quasi separation. This is called only for suspect fits.
    result = linprog(-z.sum(axis=0), A_ub=-z, b_ub=np.zeros(len(y)),
                     bounds=[(-1.0, 1.0)]*x.shape[1], method='highs')
    if not result.success:
        raise ModelError('separation_diagnostic_failed')
    return bool(-result.fun > 1e-7)


@dataclass
class LogisticFit:
    beta: np.ndarray
    covariance: np.ndarray
    loglike: float
    names: tuple
    diagnostics: dict

    def linear_predictor(self, x):
        return np.asarray(x) @ self.beta

    def predict(self, x):
        return expit(self.linear_predictor(x))


def fit_logistic(x, y, names=None):
    """Unpenalized MLE; numerical centering/scaling is undone for reporting."""
    x = np.asarray(x, float)
    y = _binary_outcome(y)
    if x.ndim != 2 or x.shape[0] != len(y) or not np.isfinite(x).all():
        raise ModelError('invalid_design_matrix')
    if not np.all(x[:, 0] == 1):
        raise ModelError('first_column_must_be_intercept')
    mean, scale = x[:, 1:].mean(axis=0), x[:, 1:].std(axis=0)
    if (scale == 0).any():
        raise ModelError('unremoved_constant_predictor')
    xs = np.column_stack([x[:, 0], (x[:, 1:]-mean)/scale])
    rank = np.linalg.matrix_rank(xs)
    if rank != xs.shape[1]:
        raise ModelError('rank_deficient_design')
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter('always')
        try:
            fit = sm.Logit(y, xs, check_rank=False).fit(method='newton', maxiter=100,
                                                       tol=1e-9, disp=False)
        except Exception as exc:
            separated = separation_lp(xs, y)
            raise ModelError('logistic_fit_failed:' + ('separation' if separated else type(exc).__name__)) from exc
    params = np.asarray(fit.params)
    cov = np.asarray(fit.cov_params())
    converged = bool(fit.mle_retvals.get('converged', False))
    score_max = float(np.max(np.abs(xs.T@(y-expit(xs@params))))/len(y))
    suspicious = (not converged or not np.isfinite(params).all() or
                  np.max(np.abs(params)) > 25 or score_max > 1e-7 or
                  any('separation' in str(w.message).lower() for w in caught))
    lp_checked = False
    if suspicious:
        lp_checked = True
        if separation_lp(xs, y):
            raise ModelError('complete_or_quasi_separation')
        if not converged or not np.isfinite(params).all() or score_max > 1e-7:
            raise ModelError('logistic_nonconvergence')
    if not np.isfinite(cov).all() or np.linalg.eigvalsh(cov).min() <= 0:
        raise ModelError('nonpositive_or_nonfinite_covariance')
    undo = np.eye(x.shape[1])
    undo[0, 1:] = -mean/scale
    undo[np.arange(1, len(scale)+1), np.arange(1, len(scale)+1)] = 1/scale
    beta, covariance = undo@params, undo@cov@undo.T
    diag = {'rank_including_intercept': int(rank), 'columns_including_intercept': int(x.shape[1]),
            'predictor_df': int(x.shape[1]-1), 'converged': converged,
            'iterations': int(fit.mle_retvals.get('iterations', -1)),
            'scaled_score_max_per_person': score_max,
            'scaled_information_condition_number': float(np.linalg.cond(np.linalg.inv(cov))),
            'separation_LP_triggered': lp_checked, 'separation_check': 'finite MLE diagnostics; LP if suspect',
            'warnings': [type(w.message).__name__ for w in caught],
            'min_fitted_probability': float(fit.predict().min()),
            'max_fitted_probability': float(fit.predict().max())}
    return LogisticFit(beta, covariance, float(fit.llf), tuple(names or range(x.shape[1])), diag)


def auc(y, p):
    y = _binary_outcome(y)
    n1 = float(y.sum()); n0 = float(len(y)-n1)
    return float((rankdata(p, method='average')[y == 1].sum()-n1*(n1+1)/2)/(n1*n0))


def performance(y, eta):
    y = _binary_outcome(y)
    eta = np.asarray(eta, float)
    p = expit(eta)
    a = float(brentq(lambda b: (y-expit(eta+b)).sum(), -1000, 1000, xtol=1e-12))
    # Calibration slope also estimates its own intercept, unlike calibration-in-the-large.
    cal = fit_logistic(np.column_stack([np.ones(len(y)), eta]), y,
                       names=('Intercept', 'prediction_logit'))
    return {'auc': auc(y, p), 'brier': float(np.mean((y-p)**2)),
            'calibration_intercept': a, 'calibration_slope': float(cal.beta[1]),
            'calibration_slope_model_intercept': float(cal.beta[0])}


def paired_metrics(y, fit0, x0, fit1, x1):
    m0, m1 = performance(y, fit0.linear_predictor(x0)), performance(y, fit1.linear_predictor(x1))
    out = {**{'M0_'+k: v for k, v in m0.items()}, **{'M1_'+k: v for k, v in m1.items()}}
    out.update(delta_auc=m1['auc']-m0['auc'], brier_gain=m0['brier']-m1['brier'])
    return out


@dataclass
class PairResult:
    preprocess0: Preprocessor
    preprocess1: Preprocessor
    fit0: LogisticFit
    fit1: LogisticFit
    train_metrics: dict
    risks: list
    outcome: str
    n: int
    events: int

    def metrics_flat(self):
        return self.train_metrics.copy()

    def or_table(self):
        table = []
        for k in PROTEIN_LEVELS[1:]:
            name = 'protein_band['+k+']'
            if name not in self.fit1.names:
                table.append({'protein_band': k, 'status': 'NE', 'reason': 'not_present_in_training'})
                continue
            i = self.fit1.names.index(name)
            beta, se = float(self.fit1.beta[i]), float(np.sqrt(self.fit1.covariance[i, i]))
            table.append({'protein_band': k, 'status': 'estimated', 'log_OR': beta, 'SE': se,
                          'OR': float(np.exp(beta)), 'OR_lower': float(np.exp(beta-1.96*se)),
                          'OR_upper': float(np.exp(beta+1.96*se))})
        return table

    def diagnostics(self):
        df = self.fit1.diagnostics['predictor_df']-self.fit0.diagnostics['predictor_df']
        stat = 2*(self.fit1.loglike-self.fit0.loglike)
        if stat < -1e-6 or df <= 0:
            raise ModelError('invalid_nested_likelihood_comparison')
        return {'M0': self.fit0.diagnostics, 'M1': self.fit1.diagnostics,
                'preprocessing_M0': self.preprocess0.metadata(),
                'preprocessing_M1': self.preprocess1.metadata(),
                'likelihood_ratio': max(0.0, stat), 'likelihood_ratio_df': df,
                'likelihood_ratio_p': float(chi2.sf(max(0.0, stat), df)),
                'n': self.n, 'events': self.events}


def standardized_risks(preprocessor, fit, data):
    rows = []
    for level in PROTEIN_LEVELS:
        if level not in preprocessor.levels['protein_band']:
            rows.append({'protein_band': level, 'status': 'NE', 'risk': None})
            continue
        altered = data.copy(deep=False).assign(protein_band=level)
        p = fit.predict(preprocessor.transform(altered))
        rows.append({'protein_band': level, 'status': 'estimated', 'risk': float(p.mean())})
    r0 = rows[0]['risk']
    for r in rows:
        r['risk_difference'] = None if r['risk'] is None else r['risk']-r0
        r['risk_difference_pp'] = None if r['risk'] is None else 100*(r['risk']-r0)
        r['risk_ratio_to_reference'] = None if r['risk'] is None else r['risk']/r0
    return rows


def fit_comparison(data, outcome='death28_lower', enhanced=False,
                   extra_continuous=(), extra_categorical=()):
    d = normalized_data(data)
    # Caller supplies the already defined common cohort; rows are never dropped silently.
    y = _binary_outcome(d[outcome])
    pp0 = Preprocessor(enhanced, False, extra_continuous, extra_categorical).fit(d)
    pp1 = Preprocessor(enhanced, True, extra_continuous, extra_categorical).fit(d)
    x0, x1 = pp0.transform(d), pp1.transform(d)
    fit0, fit1 = fit_logistic(x0, y, pp0.names), fit_logistic(x1, y, pp1.names)
    train = paired_metrics(y, fit0, x0, fit1, x1)
    pair = PairResult(pp0, pp1, fit0, fit1, train,
                      standardized_risks(pp1, fit1, d), outcome, len(d), int(y.sum()))
    pair.diagnostics()  # Also checks the nested likelihood comparison.
    return pair


def evaluate_comparison(pair, validation):
    y = _binary_outcome(validation[pair.outcome])
    return paired_metrics(y, pair.fit0, pair.preprocess0.transform(validation),
                          pair.fit1, pair.preprocess1.transform(validation))


def bootstrap_attempt(data, outcome, indices, enhanced=False,
                       extra_continuous=(), extra_categorical=()):
    indices = np.asarray(indices)
    if len(indices) != len(data) or not np.issubdtype(indices.dtype, np.integer) or \
            np.any(indices < 0) or np.any(indices >= len(data)):
        raise ModelError('invalid_bootstrap_indices')
    try:
        train = data.iloc[indices].reset_index(drop=True)
        pair = fit_comparison(train, outcome, enhanced, extra_continuous, extra_categorical)
        evaluation = evaluate_comparison(pair, data)
        return {'status': 'success', 'train': pair.metrics_flat(), 'evaluation': evaluation,
                'risks': pair.risks, 'diagnostics': pair.diagnostics()}
    except ModelError as exc:
        return {'status': 'failure', 'reason': str(exc)}


def aggregate_bootstrap(apparent, attempts, requested_attempts=500):
    """No refilling or success-only inferential distribution if any attempt failed."""
    complete = len(attempts) == requested_attempts and all(a['status'] == 'success' for a in attempts)
    status = {'requested_attempts': requested_attempts, 'attempted': len(attempts),
              'successful': sum(a['status'] == 'success' for a in attempts),
              'complete': complete, 'CI_published': complete}
    if not complete:
        status['reason'] = 'incomplete_or_failed_fixed_budget; no success-only CI'
        return status
    point = apparent.metrics_flat() if isinstance(apparent, PairResult) else apparent
    metrics = []
    for metric, value in point.items():
        train = np.array([a['train'][metric] for a in attempts])
        evaluation = np.array([a['evaluation'][metric] for a in attempts])
        optimism = float((train-evaluation).mean())
        low, high = np.quantile(train, [.025, .975])-optimism
        metrics.append({'metric': metric, 'apparent': value, 'optimism': optimism,
                        'corrected': value-optimism, 'CI_lower': float(low), 'CI_upper': float(high),
                        'CI_method': 'training statistic percentile minus mean optimism (location shift)'})
    rd = []
    for i, level in enumerate(PROTEIN_LEVELS):
        values = [a['risks'][i]['risk_difference'] for a in attempts]
        if any(v is None for v in values):
            rd.append({'protein_band': level, 'status': 'NE', 'reason': 'unsupported_in_at_least_one_attempt'})
            continue
        low, high = np.quantile(values, [.025, .975])
        rd.append({'protein_band': level, 'status': 'estimated',
                   'RD_CI_lower': float(low), 'RD_CI_upper': float(high),
                   'RD_CI_lower_pp': 100*float(low), 'RD_CI_upper_pp': 100*float(high),
                   'CI_method': 'bootstrap refit and resampled empirical X; no optimism shift'})
    status.update(metrics=metrics, risk_differences=rd)
    return status
