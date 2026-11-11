BEGIN;
SET LOCAL statement_timeout='20min';
SET LOCAL work_mem='128MB';
CREATE TEMP TABLE keys(subject_id integer,hadm_id integer,stay_id integer,source_t0 timestamp,source_l timestamp,first_supported_tx_hour double precision);
DO $$ BEGIN
 IF EXISTS(SELECT 1 FROM keys GROUP BY stay_id HAVING count(*)<>1)
 OR EXISTS(SELECT 1 FROM keys WHERE subject_id IS NULL OR hadm_id IS NULL OR stay_id IS NULL
 OR source_t0 IS NULL OR source_l IS NULL OR source_l<=source_t0)
 THEN RAISE EXCEPTION 'QC_EXTRA_INPUT_KEYS'; END IF;
END $$;

CREATE TEMP TABLE extra_keys AS
SELECT k.*,a.admittime,p.anchor_year_group,
 CASE WHEN first_supported_tx_hour IS NOT NULL
 AND first_supported_tx_hour::text!~*'^(nan|[-+]?infinity)$' AND first_supported_tx_hour<=72
 THEN source_t0+first_supported_tx_hour*interval '1h' END AS first_tx_time
FROM keys k
JOIN mimiciv_hosp.admissions a ON a.hadm_id=k.hadm_id AND a.subject_id=k.subject_id
JOIN mimiciv_hosp.patients p ON p.subject_id=k.subject_id;
CREATE UNIQUE INDEX ON extra_keys(stay_id);
ANALYZE extra_keys;
DO $$ BEGIN
 IF (SELECT count(*) FROM extra_keys)<>(SELECT count(*) FROM keys)
 THEN RAISE EXCEPTION 'QC_EXTRA_SOURCE_KEY_MATCH'; END IF;
END $$;

CREATE TEMP TABLE extra_hadm_scope AS
SELECT k.stay_id,a.subject_id,a.hadm_id,true AS is_prior
FROM extra_keys k JOIN mimiciv_hosp.admissions a
 ON a.subject_id=k.subject_id AND a.dischtime<k.admittime
UNION ALL
SELECT stay_id,subject_id,hadm_id,false FROM extra_keys;
CREATE INDEX ON extra_hadm_scope(hadm_id,subject_id);
ANALYZE extra_hadm_scope;

CREATE TEMP TABLE extra_hadm_codes AS
SELECT d.subject_id,d.hadm_id,
 bool_or((d.icd_version=9 AND replace(trim(d.icd_code),'.','')='5856')
  OR(d.icd_version=10 AND upper(replace(trim(d.icd_code),'.',''))='N186')) AS esrd,
 bool_or((d.icd_version=9 AND upper(replace(trim(d.icd_code),'.',''))='V4511')
  OR(d.icd_version=10 AND upper(replace(trim(d.icd_code),'.',''))='Z992')) AS dependence
FROM mimiciv_hosp.diagnoses_icd d
JOIN(SELECT DISTINCT subject_id,hadm_id FROM extra_hadm_scope) h
 ON h.subject_id=d.subject_id AND h.hadm_id=d.hadm_id
WHERE(d.icd_version=9 AND upper(replace(trim(d.icd_code),'.','')) IN('5856','V4511'))
 OR(d.icd_version=10 AND upper(replace(trim(d.icd_code),'.','')) IN('N186','Z992'))
GROUP BY d.subject_id,d.hadm_id;
CREATE UNIQUE INDEX ON extra_hadm_codes(subject_id,hadm_id);

CREATE TEMP TABLE extra_diagnoses AS
SELECT h.stay_id,
 coalesce(bool_or(h.is_prior AND c.esrd),false) AS prior_esrd,
 coalesce(bool_or(h.is_prior AND c.dependence),false) AS prior_dependence,
 coalesce(bool_or(NOT h.is_prior AND c.esrd),false) AS current_esrd,
 coalesce(bool_or(NOT h.is_prior AND c.dependence),false) AS current_dependence
FROM extra_hadm_scope h LEFT JOIN extra_hadm_codes c
 ON c.subject_id=h.subject_id AND c.hadm_id=h.hadm_id
GROUP BY h.stay_id;
CREATE UNIQUE INDEX ON extra_diagnoses(stay_id);

CREATE TEMP TABLE extra_scr_raw AS
SELECT k.stay_id,k.source_t0,k.source_l,k.first_tx_time,
 l.labevent_id,l.charttime,l.storetime,l.valuenum,
 lower(trim(coalesce(l.value,''))) AS raw_value,
 lower(regexp_replace(coalesce(trim(l.valueuom),''),'\s','','g')) AS unit,
 coalesce(l.comments,'')~*'(cancel|not (done|performed)|insufficient|unsatisfactory|unable|disregard)' AS invalid,
 upper(trim(coalesce(l.comments,'')))~'^(NEG|NEGATIVE|TRACE|TR)[.!]?$' AS incompatible_exact_comment
FROM extra_keys k JOIN mimiciv_hosp.labevents l
 ON l.subject_id=k.subject_id AND(l.hadm_id=k.hadm_id OR l.hadm_id IS NULL)
 AND l.charttime>=k.source_t0 AND l.charttime<k.source_l AND l.charttime<k.first_tx_time
 AND l.storetime<=k.source_l
WHERE l.itemid=50912 AND k.first_tx_time IS NOT NULL AND k.first_tx_time>k.source_t0;
CREATE INDEX ON extra_scr_raw(stay_id,charttime);
ANALYZE extra_scr_raw;

CREATE TEMP TABLE extra_scr_parsed AS
SELECT r.*,CASE WHEN raw_value~'^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$'
 THEN raw_value::double precision END AS numeric_text
FROM extra_scr_raw r;
CREATE TEMP TABLE extra_scr_checked AS
SELECT r.*,
 CASE WHEN numeric_text IS NOT NULL THEN numeric_text WHEN raw_value='' THEN valuenum END AS scr_candidate,
 (numeric_text IS NOT NULL AND valuenum IS NOT NULL
  AND abs(numeric_text-valuenum)>1e-6*greatest(1,abs(numeric_text))) AS numeric_conflict
FROM extra_scr_parsed r;
CREATE TEMP TABLE extra_scr_last AS
SELECT DISTINCT ON(stay_id) * FROM extra_scr_checked
WHERE unit='mg/dl' AND NOT invalid AND NOT incompatible_exact_comment
 AND NOT coalesce(numeric_conflict,false)
 AND(raw_value='' OR numeric_text IS NOT NULL)
 AND scr_candidate>0 AND scr_candidate::text!~*'^(nan|[-+]?infinity)$'
 AND(valuenum IS NULL OR valuenum::text!~*'^(nan|[-+]?infinity)$')
ORDER BY stay_id,charttime DESC,storetime DESC,labevent_id DESC;
CREATE UNIQUE INDEX ON extra_scr_last(stay_id);

CREATE TEMP TABLE extra_cov AS
SELECT k.stay_id,k.anchor_year_group,
 d.prior_esrd,d.prior_dependence,d.current_esrd,d.current_dependence,
 s.scr_candidate AS pre_rrt_last_scr_mgdl,
 extract(epoch FROM(s.charttime-k.source_t0))/3600 AS pre_rrt_last_scr_hour,
 extract(epoch FROM(s.storetime-k.source_t0))/3600 AS pre_rrt_last_scr_report_hour,
 s.labevent_id AS pre_rrt_last_scr_labevent_id
FROM extra_keys k JOIN extra_diagnoses d USING(stay_id)
LEFT JOIN extra_scr_last s USING(stay_id);
CREATE UNIQUE INDEX ON extra_cov(stay_id);
DO $$ BEGIN
 IF(SELECT count(*) FROM extra_cov)<>(SELECT count(*) FROM keys)
 THEN RAISE EXCEPTION 'QC_EXTRA_OUTPUT_CARDINALITY'; END IF;
 IF EXISTS(SELECT 1 FROM extra_cov e JOIN keys k USING(stay_id)
 WHERE e.pre_rrt_last_scr_mgdl IS NOT NULL
 AND(e.pre_rrt_last_scr_hour<0 OR e.pre_rrt_last_scr_hour>=k.first_supported_tx_hour
 OR e.pre_rrt_last_scr_hour>=extract(epoch FROM(k.source_l-k.source_t0))/3600
 OR e.pre_rrt_last_scr_report_hour>extract(epoch FROM(k.source_l-k.source_t0))/3600))
 THEN RAISE EXCEPTION 'QC_EXTRA_PRE_RRT_WINDOW'; END IF;
END $$;
ANALYZE extra_cov;

COMMIT;
