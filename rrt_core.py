import warnings

import numpy as np
from scipy.stats import chi2, norm
from statsmodels.duration.hazard_regression import PHReg
from statsmodels.duration.survfunc import CumIncidenceRight


EVENT_CODES = {'administrative_7d': 0, 'RRT': 1, 'death': 2, 'ICU_exit': 3}


def _inputs(time, event_type):
    time = np.asarray(time, dtype=float)
    labels = np.asarray(event_type).astype(str)
    if time.ndim != 1 or labels.shape != time.shape or not len(time):
        raise ValueError('Nonempty equal-length time/event arrays required')
    if not np.isfinite(time).all() or (time < 0).any():
        raise ValueError('Times must be finite, nonnegative days since L')
    unknown = set(labels) - set(EVENT_CODES)
    if unknown:
        raise ValueError('Unknown event types: ' + str(sorted(unknown)))
    return time, np.array([EVENT_CODES[x] for x in labels], dtype=int)


def aj_summary(time, event_type, horizon=7.0):
    time, code = _inputs(time, event_type)
    if (time > horizon + 1e-9).any():
        raise ValueError('Inputs exceed observation horizon')
    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter('always')
        fit = CumIncidenceRight(time, code)
    ix = int(np.searchsorted(fit.times, horizon, side='right') - 1)
    result = {}
    for label, event in EVENT_CODES.items():
        if event == 0:
            continue
        if event > len(fit.cinc) or ix < 0:
            point, se = 0.0, 0.0
        else:
            point = float(fit.cinc[event - 1][ix])
            se0 = float(fit.cinc_se[event - 1][ix])
            se = se0 if np.isfinite(se0) else None
        result[label] = {'cif': point, 'se': se,
                         'count': int(np.sum(code == event))}
    early_censor = (code == 0) & (time < horizon - 1e-9)
    complete = not early_censor.any()
    identity = None
    if complete:
        errors = [abs(v['cif'] - v['count'] / len(time)) for v in result.values()]
        identity = max(errors, default=0.0)
        if identity > 1e-10:
            raise AssertionError('Common-horizon complete-followup AJ != count/N')
    return {'n': len(time), 'horizon_days': float(horizon), 'causes': result,
            'early_censor_count': int(early_censor.sum()),
            'complete_common_horizon': complete,
            'count_over_N_identity_max_error': identity,
            'warnings': sorted(set(str(w.message) for w in ws))}


def _efron_groups(time, status, X, beta):
    eta = X @ beta
    weight = np.exp(eta - np.max(eta))
    groups = []
    for event_time in np.unique(time[status == 1]):
        at_event = (time == event_time) & (status == 1)
        risk = time >= event_time
        xr, xe = X[risk], X[at_event]
        wr, we = weight[risk], weight[at_event]
        m = int(at_event.sum())
        s0, e0 = wr.sum(), we.sum()
        s1, e1 = wr @ xr, we @ xe
        s2 = np.einsum('i,ij,ik->jk', wr, xr, xr)
        e2 = np.einsum('i,ij,ik->jk', we, xe, xe)
        expected_sum = np.zeros(X.shape[1])
        information = np.zeros((X.shape[1], X.shape[1]))
        for fraction in np.arange(m) / m:
            den = s0 - fraction * e0
            if not np.isfinite(den) or den <= 0:
                raise ValueError('Nonpositive Efron denominator')
            mean = (s1 - fraction * e1) / den
            expected_sum += mean
            information += (s2 - fraction * e2) / den - np.outer(mean, mean)
        score = xe.sum(axis=0) - expected_sum
        groups.append((float(event_time), m, score, information))
    return groups


def efron_ph_score(time, status, X, beta, column_names):
    """Score for beta(t)=beta+gamma*centered(log1p(days)), Efron ties."""
    time, status = np.asarray(time, float), np.asarray(status, int)
    X, beta = np.asarray(X, float), np.asarray(beta, float)
    groups = _efron_groups(time, status, X, beta)
    if not groups:
        raise ValueError('No cause-specific events')
    g = np.log1p([v[0] for v in groups])
    g -= np.average(g, weights=[v[1] for v in groups])
    residual = np.asarray([v[2] for v in groups])
    info = np.asarray([v[3] for v in groups])
    u_beta = residual.sum(axis=0)
    u_gamma = g @ residual
    i_bb = info.sum(axis=0)
    i_bg = np.einsum('i,ijk->jk', g, info)
    i_gg = np.einsum('i,ijk->jk', g * g, info)
    if np.linalg.matrix_rank(i_bb) != X.shape[1]:
        raise ValueError('Singular PH nuisance information')
    efficient_score = u_gamma - i_bg.T @ np.linalg.solve(i_bb, u_beta)
    efficient_info = i_gg - i_bg.T @ np.linalg.solve(i_bb, i_bg)
    efficient_info = (efficient_info + efficient_info.T) / 2

    def test(indices):
        ix = np.asarray(indices, dtype=int)
        s = efficient_score[ix]
        inf = efficient_info[np.ix_(ix, ix)]
        rank = int(np.linalg.matrix_rank(inf))
        if rank != len(ix) or np.linalg.eigvalsh(inf).min() <= 0:
            return {'status': 'not_estimable', 'df': rank, 'requested_df': len(ix),
                    'chi2': None, 'p': None}
        stat = float(s @ np.linalg.solve(inf, s))
        return {'status': 'ok', 'df': len(ix), 'chi2': stat,
                'p': float(chi2.sf(stat, len(ix)))}

    return {'kind': 'Efron efficient score for log1p(days-since-L) time interactions',
            'global': test(range(X.shape[1])),
            'per_column': {name: test([j]) for j, name in enumerate(column_names)},
            'max_abs_nuisance_score': float(np.max(np.abs(u_beta))),
            'number_unique_event_times': len(groups),
            'number_tied_event_times': sum(v[1] > 1 for v in groups),
            'interpretation': 'Diagnostic only; no automatic variable selection/refit; no package Schoenfeld residuals used',
            '_efficient_score': efficient_score,
            '_efficient_info': efficient_info,
            '_nuisance_info': i_bb,
            '_nuisance_cross_info': i_bg,
            '_interaction_info': i_gg}


def fit_rrt(time, event_type, X, column_names):
    """X already preprocessed/scaled, no intercept; params/cov retain X scale."""
    time, code = _inputs(time, event_type)
    X = np.asarray(X, dtype=float)
    names = list(column_names)
    if X.ndim != 2 or X.shape[0] != len(time) or X.shape[1] != len(names):
        raise ValueError('Invalid design matrix or column names')
    if not np.isfinite(X).all() or X.shape[1] == 0:
        raise ValueError('Nonempty finite design required')
    if np.any(np.ptp(X, axis=0) == 0) or np.linalg.matrix_rank(X) != X.shape[1]:
        raise ValueError('Constant/intercept or rank-deficient Cox matrix')
    status = (code == 1).astype(int)
    if not status.any():
        raise ValueError('No RRT events')
    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter('always')
        model = PHReg(time, X, status=status, ties='efron', missing='raise')
        fit = model.fit(disp=0, maxiter=100, tol=1e-8)
    params = np.asarray(fit.params, float)
    cov = np.asarray(fit.cov_params(), float)
    information = -model.hessian(params)
    score_max = float(np.max(np.abs(model.score(params))))
    warning_text = sorted(set(str(w.message) for w in ws))
    failed_warning = any('failed to converge' in x.lower() or 'overflow' in x.lower()
                         or 'invalid value' in x.lower() for x in warning_text)
    if (not np.isfinite(params).all() or not np.isfinite(cov).all()
            or np.linalg.matrix_rank(information) != X.shape[1]
            or (np.diag(cov) <= 0).any() or score_max > 1e-5 or failed_warning):
        raise RuntimeError('Cox numerical/convergence diagnostic failed: ' + str(warning_text))
    se = np.sqrt(np.diag(cov))
    ph = efron_ph_score(time, status, X, params, names)
    group_score = ph['max_abs_nuisance_score']
    if group_score > 1e-5:
        raise AssertionError('Efron group scores do not sum to fitted score')
    effects = []
    for j, name in enumerate(names):
        effects.append({'column': name, 'log_HR_X_scale': float(params[j]),
                        'SE_log_HR_X_scale': float(se[j]),
                        'HR_X_scale': float(np.exp(params[j])),
                        'lo95_X_scale': float(np.exp(params[j] - 1.96 * se[j])),
                        'hi95_X_scale': float(np.exp(params[j] + 1.96 * se[j])),
                        'p': float(2 * norm.sf(abs(params[j] / se[j])))})
    return {'params': params, 'cov': cov, 'effects_X_scale': effects,
            'n': len(time), 'RRT_events': int(status.sum()), 'df': X.shape[1],
            'ties': 'efron', 'log_partial_likelihood': float(model.loglike(params)),
            'max_abs_score': score_max,
            'information_condition_number': float(np.linalg.cond(information)),
            'warnings': warning_text, 'ph_diagnostic': ph,
            'cause_specific_handling': 'death and ICU_exit terminate risk at their times; administrative censor at 7d',
            'scale_note': 'All numeric effects/params use provided X columns; caller must back-transform continuous HR. Protein dummy contrasts unchanged if unscaled.'}


def compare_scenarios(time_early, type_early, time_late, type_late,
                      fit_early=None, fit_late=None):
    """Arrays must align to the same ordered patients; no identity inference."""
    t1, e1 = _inputs(time_early, type_early)
    t2, e2 = _inputs(time_late, type_late)
    if len(t1) != len(t2):
        raise ValueError('Scenarios require same ordered patients')
    event_times = np.unique(np.r_[t1[e1 == 1], t2[e2 == 1]])
    count_diff, membership_diff = [], []
    for t in event_times:
        r1, r2 = t1 >= t, t2 >= t
        count_diff.append(int(r1.sum() - r2.sum()))
        membership_diff.append(int(np.sum(r1 != r2)))
    aj1, aj2 = aj_summary(t1, type_early), aj_summary(t2, type_late)
    out = {'same_ordered_patient_input_required': True,
           'observation_time_differs': int(np.sum(t1 != t2)),
           'competing_type_differs': int(np.sum(e1 != e2)),
           'RRT_union_event_times': len(event_times),
           'RRT_event_times_with_riskset_count_difference': int(np.count_nonzero(count_diff)),
           'RRT_event_times_with_riskset_membership_difference': int(np.count_nonzero(membership_diff)),
           'max_riskset_count_difference': int(max(map(abs, count_diff), default=0)),
           'AJ_RRT_CIF_early': aj1['causes']['RRT']['cif'],
           'AJ_RRT_CIF_late': aj2['causes']['RRT']['cif'],
           'AJ_RRT_CIF_late_minus_early': aj2['causes']['RRT']['cif'] - aj1['causes']['RRT']['cif']}
    if fit_early is not None and fit_late is not None:
        delta = np.asarray(fit_late['params']) - np.asarray(fit_early['params'])
        out['Cox_log_HR_X_scale_late_minus_early'] = delta.tolist()
        out['Cox_max_abs_log_HR_X_scale_difference'] = float(np.max(np.abs(delta)))
    return out
