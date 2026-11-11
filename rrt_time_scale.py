import numpy as np
import pandas as pd


def _band(n):
    n = int(n)
    return '0' if n == 0 else '<10' if n < 10 else f'{n//10*10}-{n//10*10+9}'


def _truth(s):
    return s.fillna(False).astype(str).str.lower().isin(['1', 'true', 't']).to_numpy()


def read_rrt_scenarios(path):
    return pd.read_csv(path, float_precision='round_trip')


def risk_membership_check(hours, days, event_type):
    hours, days, event_type = np.asarray(hours), np.asarray(days), np.asarray(event_type)
    changed_times, max_changes, pair_changes = 0, 0, 0
    for row in np.flatnonzero(event_type == 'RRT'):
        changes = int(np.sum((hours >= hours[row]) != (days >= days[row])))
        changed_times += changes > 0
        max_changes = max(max_changes, changes)
        pair_changes += changes
    return {'RRT_times_with_risk_membership_disagreement': _band(changed_times),
            'max_patient_changes_at_one_RRT_time': _band(max_changes),
            'patient_RRT_time_pairs_disagree': _band(pair_changes),
            'all_RRT_risk_memberships_equal': changed_times == 0}


def prepare_rrt_days(frame, upper_hours, scenario_label):
    """Aligned same-patient frame must be loaded by read_rrt_scenarios."""
    if scenario_label not in ['earliest', 'latest']:
        raise ValueError('Unknown RRT date-boundary scenario')
    hours = frame.analysis_hours_after_L.to_numpy(float)
    upper = np.asarray(upper_hours, float)
    if upper.shape != hours.shape or not np.isfinite(hours).all():
        raise ValueError('RRT times/upper bounds are invalid or not aligned')
    intended = ((scenario_label == 'latest') & _truth(frame.death_is_date_interval)
                & frame.event_type.eq('death').to_numpy())
    if np.any(intended & ~np.isfinite(upper)):
        raise ValueError('Latest interval death requires finite source upper bound')
    original_is_left = hours == np.nextafter(upper, -np.inf)
    if np.any(intended & ~original_is_left):
        raise ValueError('Latest interval death is not the verified source-hour left limit; check round_trip parsing')
    days = hours / 24
    naive = days.copy()
    days[intended] = np.nextafter(upper[intended] / 24, -np.inf)
    strict_left = bool(np.all(days[intended] < upper[intended] / 24))
    if not strict_left:
        raise AssertionError('Latest interval death must be strictly before upper-day boundary')
    membership = risk_membership_check(hours, days, frame.event_type.to_numpy())
    if not membership['all_RRT_risk_memberships_equal']:
        raise AssertionError('Hours-to-days conversion changed an RRT risk membership')
    qa = {'scenario': scenario_label, 'scenario_CSV_float_parser': 'round_trip',
          'source_hour_left_limit_verified': True,
          'source_hour_left_limit_people_text': _band(np.sum(intended)),
          'days_reconstructed_people_text': _band(np.sum(intended)),
          'day_values_changed_by_reconstruction_text': _band(np.sum(days != naive)),
          'latest_interval_death_strictly_before_day_upper': strict_left,
          'source_hours_vs_final_days_RRT_risk_memberships': membership,
          'main_predictor_matrix_changed': False}
    return days, qa
