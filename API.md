# P03.3 logistic core

`modeling.py` does not load files or print data. Supply one already defined common analysis cohort; it never silently excludes rows. `protein_band` must be one of the six locked values and the outcome must be complete binary data.

```python
pair = fit_comparison(data, outcome='death28_lower', enhanced=False)
apparent = pair.metrics_flat()
effects = pair.or_table()
risks = pair.risks
diagnostics = pair.diagnostics()
validation = evaluate_comparison(pair, validation_data)
one = bootstrap_attempt(data, 'death28_lower', indices, enhanced=False)
summary = aggregate_bootstrap(pair, attempts, requested_attempts=500)
```

For an enhanced pair, internal labels remain `M0`/`M1`; the caller should label these E0/E1. `extra_continuous` and `extra_categorical` can append already specified sensitivity predictors. New science choices are not selected by this module.

Imputation parameters, categorical support, reference levels and constant-column drops are learned only from each training sample. Validation uses that fitted object. The categorical Unknown level is not an extra duplicated missing dummy. Protein cannot be imputed; an absent reference stops fitting and an unseen validation protein band stops evaluation.

Numerical scaling during MLE is undone in coefficients and covariance. Age is explicitly per 10 years. Column names, rank including intercept, fit convergence and diagnostic conditions are recorded. Finite MLE diagnostics are supplemented by a linear-programming separation test when a fit is suspect. There is no penalized fallback or P-based selection.

`M0_calibration_intercept` is the offset estimate with slope fixed to 1. `M0_calibration_slope` estimates its own intercept as well; the accompanying `_slope_model_intercept` is diagnostic. Paired increment directions are `delta_auc=M1-M0` and `brier_gain=M0-M1`.

Standardized risks average predicted binary death probabilities over every person's covariates, independently assigning each of the six protein bands. They are not category probabilities and need not sum to 1. RD is signed and includes probability and percentage-point forms.

The caller controls bootstrap indices and separate trial/formal seed streams. This module does not draw until a success target is reached. Only a complete, wholly successful fixed-attempt distribution produces bootstrap CIs. RD percentile intervals are computed from resampled empirical covariate distributions and receive no optimism shift. Performance intervals use the locked location-shift approximation.

Run `test_modeling.py` for 27 synthetic checks. The tests compare the scaled Newton fit against unscaled statsmodels GLM IRLS, plus schema, leakage, standardization, increment direction, bootstrap formulas and failure behavior. These are engineering checks, not independent clinical replication or P04.
