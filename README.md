# Supplementary Code S1

**Reproduction prerequisites:** First obtain authorized access to MIMIC-IV v3.1. Run the SQL extraction scripts in the required local PostgreSQL environment, complete the local cleaning and covariate-merge steps to produce the ordered analysis-ready dataset, and then run the analysis code. The sections below specify the additional local pipeline dependencies.

**复现前提：**须先取得 MIMIC-IV 3.1 访问权限，运行 SQL 脚本提取数据，并完成清洗、协变量合并和分析输入准备后，再执行分析代码。

**Reproduction input order:** Use the original analysis-ready CSV in its original row order. The fixed bootstrap stream samples row positions. Reordering the file maps the same draws to different patients and can change the bootstrap results. The output JSON records the input SHA-256 for exact replay.

## Scope

This code accompanies a study of the additional mortality information in early routine urine protein reports among adults with early sepsis-associated acute kidney injury. The primary study population consists of patients with an interpretable urine protein report who survived and remained in the ICU at the 72-hour landmark. Follow-up for all-cause mortality ends 28 days after ICU admission. The two death columns encode the prespecified alternative classifications at the date-only endpoint boundary.

The archive contains the frozen mortality and renal replacement therapy (RRT) modeling modules, source SQL for cohort and laboratory definitions, the main analysis input header, and a portable entry point for the four mortality comparisons. MIMIC-IV patient data, row-level exports, fitted objects, stored resamples and access details are kept in the approved data environment. Access to the database requires the applicable [PhysioNet credentialing and data-use process](https://physionet.org/content/mimiciv/3.1/).

## Run the four mortality comparisons

The production environment used Python 3.12.5 and the package versions in `requirements.txt`. In a dedicated environment, install those dependencies and provide the locally prepared main-cohort file.

```sh
python -m pip install -r requirements.txt
python fit_mortality.py --input /authorized/local/P032_ANALYSIS_INPUT_WITH_EXTRA.csv --output /authorized/local/mortality_summary.json
```

The required input is the main-cohort cleaning and covariate-merge output, with the exact header in `analysis_input_columns.csv` and one row per patient and index ICU stay. Leave genuinely unavailable predictors empty. Do not use example or synthetic values in place of unavailable fields. Death indicators must be complete binary values for the eligible cohort. The adapter selects only nonmissing `protein_band` and never imputes urine protein. The six report levels and predictor columns are defined in `modeling.py`.

The script uses the production seed and direct child stream, pairs all configurations on the same 500 resamples, and does not refill failures. Preprocessing is learned within each training sample. The JSON provides apparent measures, standardized risk, odds ratios, and fixed-budget bootstrap summaries. Calibration intervals are withheld when the training statistics are structurally degenerate. An incomplete fixed attempt set yields no inferential interval based on successful attempts alone. Store the output in an authorized local directory. It contains exact aggregate sample and event counts.

## Source SQL and additional analyses

| Files | Role |
| --- | --- |
| `sql/01_primary_cohort.sql` through `sql/04_uncertain_landmark_day.sql` | The executed source for the primary and three alternate cohorts. In a credentialed PostgreSQL session, the final temporary `extraction_cohort` relation supplies the source for `P032_CLINICAL_COHORT.csv`. |
| `sql/05_extra_covariates.sql` | Additional year group, prior/current dialysis coding, and pre-RRT creatinine definitions. The production wrapper loaded local keys into the temporary `keys` relation and exported `extra_cov` as `P032_EXTRA_COVARIATES.csv`. |
| `sql/06_measurement_supplement.sql` | Measurement selection and loss definitions. The production wrapper loaded temporary `keys` and `memberships` relations. Local `measurement_supplement` and banded `measurement_summary` are the corresponding temporary outputs. |
| `sql/protein_parser.sql` | Standalone definition embedded in the applicable extraction and measurement SQL. |
| `modeling.py`, `API.md`, `test_modeling.py` | Frozen logistic methods, their interface and a synthetic engineering test. Run the test in a disposable local copy because it writes QA files. |
| `rrt_core.py`, `rrt_core_fix.py`, `rrt_time_scale.py` | The underlying RRT and competing-event routines. The final analysis used `rrt_core_fix.fit_rrt`, `rrt_core_fix.aj_summary` and the scenario time handling in `rrt_time_scale.py`. |

The SQL relies on PostgreSQL temporary relations and local `mimiciv_hosp`, `mimiciv_icu` and `mimiciv_derived` schemas. The derivation build commit and execution date for the study's upstream tables could not be established. A newly built upstream schema may change cohort membership or predictors. The source SQL documents the implemented rules. End-to-end regeneration of the cleaned CSVs, the RRT event-scenario inputs, the remaining sensitivity orchestration, and the final manuscript tables still requires the study specification and an approved local pipeline.

## Provenance and checks

`CODE_MANIFEST.json` records the checksums and scope of this archive. The frozen statistical modules and executed SQL match the recorded production hashes. The portable mortality adapter is a packaging addition. Its packaging checks covered syntax, imports and the help entry point. They did not use the source database or refit clinical models.

The documented P11 R09 replay of the preceding archive ran the four mortality comparisons on the registered ordered analysis-ready CSV, with 500 completed resamples per configuration. This documentation revision preserves the executed adapter and analysis files byte for byte.

## Release, license and citation

Version **1.0.2**, released 2026-11-11. The code is distributed under the [MIT License](LICENSE). See [CITATION.cff](CITATION.cff) for the software authors and citation metadata.

Repository: https://github.com/LFY123-PRO/mimiciv-sa-aki-urine-protein

The software license applies to the code. Access to MIMIC-IV data remains governed by PhysioNet credentialing and the applicable data-use agreement.
