"""Run the four prespecified mortality comparisons on a controlled local CSV."""

import argparse
import csv
from hashlib import sha256
import json
import os
from pathlib import Path

for variable in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ[variable] = "1"

import numpy as np
import pandas as pd

from modeling import (
    PROTEIN_LEVELS, SEED, aggregate_bootstrap, bootstrap_attempt, fit_comparison,
)


ATTEMPTS = 500
CONFIGURATIONS = (
    ("main_lower", "death28_lower", False),
    ("main_upper", "death28_upper", False),
    ("enhanced_lower", "death28_lower", True),
    ("enhanced_upper", "death28_upper", True),
)


def load_measured(path):
    schema = Path(__file__).with_name("analysis_input_columns.csv")
    with schema.open(encoding="utf-8", newline="") as handle:
        required_header = next(csv.reader(handle))
    data = pd.read_csv(path)
    if list(data.columns) != required_header:
        raise ValueError("The input header does not match the frozen main-cohort column map.")
    for column in ("subject_id", "stay_id"):
        if data[column].isna().any() or data[column].duplicated().any():
            raise ValueError(f"{column} must be present and unique before selecting measured reports.")
    return data.loc[data.protein_band.notna()].reset_index(drop=True)


def interval_summary(attempts):
    risks = np.asarray([[category["risk"] for category in row["risks"]] for row in attempts], dtype=float)
    intervals = {}
    for index, level in enumerate(PROTEIN_LEVELS):
        estimates = risks[:, index]
        if not np.isfinite(estimates).all():
            continue
        row = {"risk_CI95": np.quantile(estimates, [.025, .975]).tolist()}
        if index and np.isfinite(risks[:, 0]).all() and (risks[:, 0] > 0).all():
            row["risk_ratio_to_reference_CI95"] = np.quantile(estimates / risks[:, 0], [.025, .975]).tolist()
        intervals[level] = row
    return intervals


def run(data):
    apparent = {
        name: fit_comparison(data, outcome, enhanced)
        for name, outcome, enhanced in CONFIGURATIONS
    }
    attempts = {name: [] for name, _, _ in CONFIGURATIONS}
    generator = np.random.default_rng(np.random.SeedSequence(SEED, spawn_key=(1,)))
    for _ in range(ATTEMPTS):
        indices = generator.integers(0, len(data), len(data))
        for name, outcome, enhanced in CONFIGURATIONS:
            result = bootstrap_attempt(data, outcome, indices, enhanced)
            attempts[name].append(result)
    comparisons = {}
    for name, _, _ in CONFIGURATIONS:
        pair = apparent[name]
        summary = aggregate_bootstrap(pair, attempts[name], ATTEMPTS)
        if summary["complete"]:
            summary["standardized_risk_and_ratio_intervals"] = interval_summary(attempts[name])
            for row in summary["metrics"]:
                if row["metric"].endswith(("calibration_intercept", "calibration_slope", "slope_model_intercept")):
                    row["CI_lower"] = row["CI_upper"] = None
                    row["CI_method"] = "Withheld because training calibration statistics are structurally degenerate"
        comparisons[name] = {
            "apparent": pair.metrics_flat(),
            "standardized_risks": pair.risks,
            "adjusted_OR": pair.or_table(),
            "n": pair.n,
            "events": pair.events,
            "bootstrap": summary,
        }
    return comparisons


def serialize(value):
    if isinstance(value, np.generic):
        return value.item()
    raise TypeError(type(value).__name__)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="Controlled P032_ANALYSIS_INPUT_WITH_EXTRA.csv")
    parser.add_argument("--output", type=Path, required=True, help="New JSON file in an authorized local directory")
    args = parser.parse_args()
    if args.output.resolve() == args.input.resolve() or args.output.exists():
        parser.error("--output must be a new file distinct from the input")
    if not args.output.parent.is_dir():
        parser.error("--output parent directory must already exist")
    result = {
        "input_sha256": sha256(args.input.read_bytes()).hexdigest(),
        "seed": SEED,
        "spawn_key": [1],
        "comparisons": run(load_measured(args.input)),
    }
    document = json.dumps(result, indent=2, default=serialize, allow_nan=False) + "\n"
    with args.output.open("x", encoding="utf-8") as handle:
        handle.write(document)


if __name__ == "__main__":
    main()
