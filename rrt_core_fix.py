import warnings

import numpy as np
from scipy.optimize import minimize
from scipy.stats import norm
from statsmodels.duration.hazard_regression import PHReg

from rrt_core import _inputs, aj_summary, compare_scenarios, efron_ph_score


def stable_efron(time, status, X, beta, need_info=True):
    """Same unpenalized Efron likelihood, direct risk sums with local shifts."""
    eta = X @ beta
    ll = 0.0
    grad = np.zeros(X.shape[1])
    info = np.zeros((X.shape[1], X.shape[1]))
    for t in np.unique(time[status == 1]):
        event = (time == t) & (status == 1)
        risk = time >= t
        xr, xe = X[risk], X[event]
        shift = float(np.max(eta[risk]))
        wr, we = np.exp(eta[risk] - shift), np.exp(eta[event] - shift)
        m = int(event.sum())
        s0, e0 = wr.sum(), we.sum()
        s1, e1 = wr @ xr, we @ xe
        if need_info:
            s2 = np.einsum('i,ij,ik->jk', wr, xr, xr)
            e2 = np.einsum('i,ij,ik->jk', we, xe, xe)
        ll += eta[event].sum() - m * shift
        grad += xe.sum(axis=0)
        for frac in np.arange(m) / m:
            den = s0 - frac * e0
            if den <= 0 or not np.isfinite(den):
                raise RuntimeError('Invalid Efron risk denominator')
            mean = (s1 - frac * e1) / den
            ll -= np.log(den)
            grad -= mean
            if need_info:
                info += (s2 - frac * e2) / den - np.outer(mean, mean)
    return ll, grad, (info + info.T) / 2


def fit_rrt(time, event_type, X, column_names):
    """Optimizer fix only; identical predictors, Efron ties, no penalty."""
    time, code = _inputs(time, event_type)
    X = np.asarray(X, dtype=float)
    names = list(column_names)
    if X.ndim != 2 or X.shape[0] != len(time) or X.shape[1] != len(names):
        raise ValueError('Invalid Cox design matrix or names')
    if not np.isfinite(X).all() or X.shape[1] == 0:
        raise ValueError('Nonempty finite Cox design required')
    if np.any(np.ptp(X, axis=0) == 0) or np.linalg.matrix_rank(X) != X.shape[1]:
        raise ValueError('Constant/intercept or rank-deficient Cox design')
    status = (code == 1).astype(int)
    if not status.any():
        raise ValueError('No RRT events')

    def objective(beta):
        ll, score, _ = stable_efron(time, status, X, beta, need_info=False)
        return -ll, -score

    with warnings.catch_warnings(record=True) as ws:
        warnings.simplefilter('always')
        opt = minimize(objective, np.zeros(X.shape[1]), method='BFGS', jac=True,
                       options={'gtol': 1e-6, 'maxiter': 500})
    beta = opt.x.copy()
    refinements = []
    for iteration in range(5):
        ll, score, info = stable_efron(time, status, X, beta)
        if not np.isfinite(score).all() or not np.isfinite(info).all():
            raise RuntimeError('Nonfinite stable Efron derivative')
        if np.max(np.abs(score)) <= 1e-7:
            break
        if np.linalg.eigvalsh(info).min() <= 0:
            raise RuntimeError('Nonpositive observed Cox information')
        step = np.linalg.solve(info, score)
        length = 1.0
        accepted = False
        for halving in range(30):
            candidate = beta + length * step
            next_ll, next_score, _ = stable_efron(time, status, X, candidate, need_info=False)
            # At numerical precision, permit tiny ll roundoff only with improved score.
            if (next_ll >= ll + 1e-4 * length * score @ step
                    or (next_ll >= ll - 1e-10 and np.max(abs(next_score)) < np.max(abs(score)))):
                beta = candidate
                refinements.append({'iteration': iteration, 'step_fraction': length,
                                    'll_increase': float(next_ll - ll),
                                    'max_abs_score': float(np.max(abs(next_score)))})
                accepted = True
                break
            length *= 0.5
        if not accepted:
            raise RuntimeError('Cox line-search refinement failed')
    ll, score, info = stable_efron(time, status, X, beta)
    eig = np.linalg.eigvalsh(info)
    if (not np.isfinite(beta).all() or not np.isfinite(info).all()
            or np.max(abs(score)) > 1e-5
            or np.linalg.matrix_rank(info) != X.shape[1] or eig.min() <= 0):
        raise RuntimeError('Stable Cox MLE derivative/information check failed')
    cov = np.linalg.inv(info)
    if not np.isfinite(cov).all() or (np.diag(cov) <= 0).any():
        raise RuntimeError('Stable Cox covariance check failed')
    # The package objective/derivatives are an independently implemented cross-check.
    model = PHReg(time, X, status=status, ties='efron', missing='raise')
    pll, ps = float(model.loglike(beta)), model.score(beta)
    phess = -model.hessian(beta)
    differences = {'ll': float(abs(pll - ll)),
                   'score_max_abs': float(np.max(abs(ps - score))),
                   'information_max_abs': float(np.max(abs(phess - info)))}
    if (differences['ll'] > 1e-7 or differences['score_max_abs'] > 1e-7
            or differences['information_max_abs'] > 1e-6):
        raise AssertionError('Stable Efron and package objective/derivatives disagree')
    se = np.sqrt(np.diag(cov))
    ph = efron_ph_score(time, status, X, beta, names)
    if ph['max_abs_nuisance_score'] > 1e-5:
        raise AssertionError('Efron PH group scores differ from fitted MLE score')
    effects = [{'column': name, 'log_HR_X_scale': float(beta[j]),
                'SE_log_HR_X_scale': float(se[j]), 'HR_X_scale': float(np.exp(beta[j])),
                'lo95_X_scale': float(np.exp(beta[j] - 1.96 * se[j])),
                'hi95_X_scale': float(np.exp(beta[j] + 1.96 * se[j])),
                'p': float(2 * norm.sf(abs(beta[j] / se[j])))} for j, name in enumerate(names)]
    return {'params': beta, 'cov': cov, 'effects_X_scale': effects,
            'n': len(time), 'RRT_events': int(status.sum()), 'df': X.shape[1],
            'ties': 'efron', 'log_partial_likelihood': float(ll),
            'max_abs_score': float(np.max(abs(score))),
            'information_condition_number': float(np.linalg.cond(info)),
            'warnings': sorted(set(str(w.message) for w in ws)),
            'ph_diagnostic': ph,
            'optimizer': {'method': 'stable Efron BFGS with line search, optional step-halved Newton refinement',
                          'initial_BFGS_success': bool(opt.success), 'initial_BFGS_message': str(opt.message),
                          'initial_BFGS_iterations': int(opt.nit), 'refinements': refinements,
                          'final_acceptance': 'independent derivative/information checks; optimizer message retained',
                          'penalty_added': False, 'predictors_removed': False,
                          'min_information_eigenvalue': float(eig.min()),
                          'package_vs_stable_differences': differences},
            'cause_specific_handling': 'death and ICU_exit terminate risk at their times; administrative censor at 7d',
            'scale_note': 'Params/HR use provided X scale; caller back-transforms continuous/dummy scaling.'}
