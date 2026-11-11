BEGIN;
SET LOCAL statement_timeout='20min';
SET LOCAL lock_timeout='5s';
SET LOCAL work_mem='128MB';
CREATE TEMP TABLE keys(subject_id integer,hadm_id integer,stay_id integer,source_t0 timestamp,
 source_l timestamp,first_supported_tx_hour double precision);
CREATE TEMP TABLE memberships(cohort text,stay_id integer,expected_band text,expected_specimen_id bigint);
-- Native client COPY supplies local keys and membership tables before this assertion.
DO $$ BEGIN
 IF EXISTS(SELECT 1 FROM keys GROUP BY stay_id HAVING count(*)<>1)
 OR EXISTS(SELECT 1 FROM memberships GROUP BY cohort,stay_id HAVING count(*)<>1)
 THEN RAISE EXCEPTION 'QC_MEASUREMENT_INPUT_KEYS'; END IF;
END $$;
CREATE UNIQUE INDEX ON keys(stay_id);
ANALYZE keys;
CREATE TEMP TABLE p_target AS
SELECT k.*,k.source_t0 AS intime,k.source_l AS landmark_time,a.admittime,i.outtime
FROM keys k JOIN mimiciv_hosp.admissions a USING(subject_id,hadm_id)
JOIN mimiciv_icu.icustays i USING(subject_id,hadm_id,stay_id);
CREATE UNIQUE INDEX ON p_target(stay_id);
ANALYZE p_target;
DO $$ BEGIN
 IF(SELECT count(*) FROM p_target)<>(SELECT count(*) FROM keys)
 THEN RAISE EXCEPTION 'QC_MEASUREMENT_SOURCE_MATCH'; END IF;
END $$;
CREATE OR REPLACE FUNCTION pg_temp.protein_parse(v_raw text,n_raw double precision,u_raw text,c_raw text)
RETURNS TABLE(protein_band text,field_class text,comment_class text,loss_reason text)
LANGUAGE sql IMMUTABLE AS $$
WITH norm AS(
 SELECT lower(trim(coalesce(v_raw,''))) AS v,upper(trim(coalesce(c_raw,''))) AS c,
 lower(regexp_replace(coalesce(trim(u_raw),''),'\s','','g')) AS u,
 n_raw,coalesce(c_raw,'')~*'(cancel|not (done|performed)|insufficient|unsatisfactory|unable|disregard)' AS invalid
), cls AS(
 SELECT *,CASE WHEN v~'^[<>]=?\s*[0-9]+([.][0-9]+)?$' THEN 'inequality'
 WHEN v IN('neg','negative') THEN 'negative' WHEN v IN('tr','trace') THEN 'trace'
 WHEN regexp_replace(v,'\s','','g') IN('1+','2+','3+','4+','+','++','+++','++++') THEN 'explicit_grade'
 WHEN v~'^[0-9]+([.][0-9]+)?$' OR(v='' AND n_raw IS NOT NULL) THEN 'numeric'
 WHEN v='' THEN 'empty' ELSE 'unparsed_field' END AS fc,
 CASE WHEN c~'^(NEG|NEGATIVE)[.!]?$' THEN 'negative' WHEN c~'^(TRACE|TR)[.!]?$' THEN 'trace'
 WHEN c IN('1+','2+','3+','4+','+','++','+++','++++') THEN 'grade' WHEN c='' THEN 'empty' ELSE 'other' END AS cc,
 CASE WHEN v~'^[0-9]+([.][0-9]+)?$' THEN v::double precision END AS numtext,
 CASE WHEN v~'^[<>]=?\s*[0-9]+([.][0-9]+)?$' THEN regexp_replace(v,'^[<>]=?\s*','')::double precision END AS bound,
 CASE WHEN v~'^[<>]=?\s*[0-9]+([.][0-9]+)?$' THEN substring(v FROM '^([<>]=?)') END AS op
 FROM norm
), checked AS(
 SELECT *,coalesce(numtext,CASE WHEN v='' THEN n_raw END) AS n,
 CASE WHEN invalid THEN 'invalid_or_cancelled'
 WHEN u<>'mg/dl' THEN 'unit_not_supported_or_missing'
 WHEN n_raw::text~*'^(nan|[-+]?infinity)$' THEN 'nonfinite_numeric'
 WHEN numtext IS NOT NULL AND n_raw IS NOT NULL AND abs(numtext-n_raw)>1e-6*greatest(1,abs(numtext)) THEN 'numeric_conflict'
 WHEN fc='inequality' AND n_raw IS NOT NULL AND abs(bound-n_raw)>1e-6*greatest(1,abs(bound)) THEN 'inequality_limit_conflict'
 WHEN fc='inequality' AND cc IN('negative','trace') AND NOT(cc='negative' AND op='<=' AND bound=0) THEN 'interval_token_semantics_unverified'
 WHEN cc='grade' AND fc NOT IN('empty','explicit_grade') THEN 'grade_field_semantics_unverified'
 WHEN fc='trace' AND n_raw IS NOT NULL THEN 'trace_numeric_semantics_unverified'
 WHEN fc='numeric' AND cc='trace' THEN 'trace_numeric_semantics_unverified'
 WHEN fc='negative' AND n_raw IS NOT NULL AND n_raw<>0 THEN 'negative_numeric_conflict'
 WHEN (fc='trace' AND cc='negative') OR(fc='negative' AND cc='trace') THEN 'token_conflict'
 WHEN cc='negative' AND coalesce(numtext,n_raw)>0 THEN 'negative_numeric_conflict'
 WHEN fc='unparsed_field' THEN 'unparsed_field_no_numeric_override'
 WHEN fc='explicit_grade' OR(fc='empty' AND cc='grade') THEN 'explicit_plus_grade_separate'
 END AS bad
 FROM cls
), mapped AS(
 SELECT *,CASE WHEN bad IS NOT NULL THEN NULL
 WHEN fc='negative' OR(fc='empty' AND cc='negative') THEN '0_negative_or_report0'
 WHEN fc='trace' OR(fc='empty' AND cc='trace') THEN '1_trace'
 WHEN fc='inequality' AND op='<=' AND bound=0 THEN '0_negative_or_report0'
 WHEN fc='inequality' AND op IN('>','>=') AND bound>=300 THEN '5_ge300'
 WHEN fc='numeric' AND n=0 THEN '0_negative_or_report0'
 WHEN fc='numeric' AND n>0 AND n<30 THEN '2_positive_lt30'
 WHEN fc='numeric' AND n>=30 AND n<100 THEN '3_30_lt100'
 WHEN fc='numeric' AND n>=100 AND n<300 THEN '4_100_lt300'
 WHEN fc='numeric' AND n>=300 THEN '5_ge300' END AS band
 FROM checked
)
SELECT band,fc,cc,CASE WHEN band IS NOT NULL THEN NULL ELSE coalesce(bad,
 CASE WHEN fc='inequality' THEN 'inequality_not_unique_band'
 WHEN fc='empty' AND cc='other' THEN 'other_comment_without_field' ELSE 'missing_or_unparsed' END) END FROM mapped
$$;
CREATE TEMP TABLE p_lab AS
SELECT t.stay_id,t.subject_id,t.hadm_id,t.intime,t.landmark_time,t.outtime,l.labevent_id,l.specimen_id,l.charttime,l.storetime,
 x.*,(lower(trim(coalesce(l.value,'')))='' AND l.valuenum IS NULL AND x.comment_class IN('negative','trace')) AS comment_only_source,
 l.value IS NOT NULL AND trim(l.value)<>'' AND x.field_class='unparsed_field' AND l.valuenum IS NOT NULL AS unparsed_with_number
FROM p_target t JOIN mimiciv_hosp.labevents l ON l.subject_id=t.subject_id AND(l.hadm_id=t.hadm_id OR l.hadm_id IS NULL)
 AND l.charttime>=greatest(t.intime,t.admittime) AND l.charttime<t.intime+interval '48h' AND l.charttime<t.outtime
CROSS JOIN LATERAL pg_temp.protein_parse(l.value,l.valuenum::double precision,l.valueuom,l.comments) x
WHERE l.itemid=51492;
CREATE INDEX ON p_lab(stay_id,specimen_id);ANALYZE p_lab;
CREATE TEMP TABLE p_canonical AS
SELECT DISTINCT ON(stay_id,specimen_id) * FROM p_lab WHERE storetime<=landmark_time AND specimen_id IS NOT NULL
ORDER BY stay_id,specimen_id,storetime DESC,labevent_id DESC;
CREATE UNIQUE INDEX ON p_canonical(stay_id,specimen_id);
CREATE TEMP TABLE p_first_good AS
SELECT DISTINCT ON(stay_id) * FROM p_canonical WHERE protein_band IS NOT NULL
ORDER BY stay_id,charttime,specimen_id,labevent_id;
CREATE UNIQUE INDEX ON p_first_good(stay_id);
CREATE TEMP TABLE p_first_any AS
SELECT DISTINCT ON(stay_id) * FROM p_lab ORDER BY stay_id,charttime,specimen_id,storetime NULLS LAST,labevent_id;
CREATE UNIQUE INDEX ON p_first_any(stay_id);

CREATE TEMP TABLE cov_anchor AS SELECT * FROM p_target;
CREATE UNIQUE INDEX ON cov_anchor(stay_id);
SET LOCAL enable_nestloop=on;
CREATE TEMP TABLE cov_specimen_all AS
SELECT g.stay_id,l.labevent_id,l.itemid,l.charttime,l.storetime,l.valuenum,
 lower(regexp_replace(coalesce(trim(l.valueuom),''),'\s','','g')) AS unit,
 lower(trim(coalesce(l.value,''))) AS raw_value,
 coalesce(l.comments,'')~*'(cancel|not (done|performed)|insufficient|unsatisfactory|unable|disregard)' AS invalid,
 upper(trim(coalesce(l.comments,'')))~'^(NEG|NEGATIVE)[.!]?$' AS exact_negative_comment,
 upper(trim(coalesce(l.comments,'')))~'^(TRACE|TR)[.!]?$' AS exact_trace_comment
FROM p_first_good g JOIN cov_anchor t USING(stay_id)
JOIN mimiciv_hosp.labevents l ON l.subject_id=t.subject_id AND l.specimen_id=g.specimen_id
 AND (l.hadm_id=t.hadm_id OR l.hadm_id IS NULL)
WHERE l.itemid IN(51498,51493,51516);
CREATE TEMP TABLE cov_specimen_raw AS
SELECT r.* FROM cov_specimen_all r JOIN keys k USING(stay_id) WHERE r.storetime<=k.source_l;
CREATE INDEX ON cov_specimen_all(stay_id,itemid);
CREATE TEMP TABLE cov_specimen_latest AS
SELECT DISTINCT ON(stay_id,itemid) * FROM cov_specimen_raw
ORDER BY stay_id,itemid,storetime DESC,labevent_id DESC;
CREATE UNIQUE INDEX ON cov_specimen_latest(stay_id,itemid);
CREATE TEMP TABLE cov_specimen_parsed AS
SELECT r.*,
 CASE WHEN raw_value~'^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$'
 THEN raw_value::double precision END AS numeric_text
FROM cov_specimen_latest r;
CREATE TEMP TABLE cov_specimen_class AS
SELECT r.*,
 CASE WHEN numeric_text IS NOT NULL THEN numeric_text WHEN raw_value='' THEN valuenum END AS numeric_candidate,
 CASE WHEN numeric_text IS NOT NULL AND valuenum IS NOT NULL
 THEN abs(numeric_text-valuenum)>1e-6*greatest(1,abs(numeric_text)) ELSE false END AS numeric_conflict,
 raw_value<>'' AND numeric_text IS NULL AS unparsed_value,
 (unit IN('/hpf','#/hpf','cells/hpf','hpf') OR (itemid=51493 AND unit='rbc/hpf') OR (itemid=51516 AND unit='wbc/hpf')) AS hpf_unit
FROM cov_specimen_parsed r;
CREATE TEMP TABLE cov_specimen_values AS
SELECT *,CASE WHEN invalid OR numeric_conflict OR unparsed_value OR exact_trace_comment
 OR (exact_negative_comment AND numeric_candidate<>0)
 OR numeric_candidate::text~*'^(nan|[-+]?infinity)$' THEN NULL
 WHEN itemid=51498 AND numeric_candidate BETWEEN 1.000 AND 1.100 THEN numeric_candidate
 WHEN itemid IN(51493,51516) AND hpf_unit AND numeric_candidate>=0 THEN numeric_candidate END AS clean_numeric
FROM cov_specimen_class;
CREATE TEMP TABLE cov_specimen_pivot AS
SELECT stay_id,
 max(clean_numeric) FILTER(WHERE itemid=51498) AS urine_specific_gravity,
 max(clean_numeric) FILTER(WHERE itemid=51493) AS urine_rbc_hpf,
 max(clean_numeric) FILTER(WHERE itemid=51516) AS urine_wbc_hpf,
 count(*) FILTER(WHERE numeric_conflict) AS specimen_numeric_conflict_items,
 count(*) FILTER(WHERE unparsed_value AND valuenum IS NOT NULL) AS specimen_unparsed_value_with_valuenum_items,
 count(*) FILTER(WHERE itemid IN(51493,51516) AND NOT hpf_unit) AS specimen_nonhpf_unit_items,
 count(*) FILTER(WHERE invalid) AS specimen_invalid_items
FROM cov_specimen_values GROUP BY stay_id;
CREATE UNIQUE INDEX ON cov_specimen_pivot(stay_id);


CREATE TEMP TABLE extra_keys AS
SELECT k.*,CASE WHEN first_supported_tx_hour IS NOT NULL
 AND first_supported_tx_hour::text!~*'^(nan|[-+]?infinity)$' AND first_supported_tx_hour<=72
 THEN source_t0+first_supported_tx_hour*interval '1h' END AS first_tx_time
FROM keys k;
CREATE UNIQUE INDEX ON extra_keys(stay_id);
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


CREATE TEMP TABLE scuf_after_L AS
SELECT x.stay_id,min(x.tx_time) AS first_scuf_after_L FROM(
 SELECT t.stay_id,p.starttime AS tx_time FROM p_target t JOIN mimiciv_icu.procedureevents p
 ON p.stay_id=t.stay_id AND p.hadm_id=t.hadm_id AND p.subject_id=t.subject_id
 WHERE true
 AND p.itemid=225955 AND coalesce(lower(p.statusdescription),'')!~'(cancel|rewrit)'
 AND p.starttime>t.landmark_time AND p.starttime<=least(t.landmark_time+interval '7d',t.outtime)
 UNION ALL SELECT t.stay_id,r.charttime FROM p_target t JOIN mimiciv_derived.rrt r ON r.stay_id=t.stay_id
 WHERE r.dialysis_active=1 AND upper(trim(r.dialysis_type))='SCUF'
 AND r.charttime>t.landmark_time AND r.charttime<=least(t.landmark_time+interval '7d',t.outtime)
) x GROUP BY x.stay_id;
CREATE UNIQUE INDEX ON scuf_after_L(stay_id);
CREATE TEMP TABLE uo_raw6h AS
SELECT t.stay_id,sum(o.value) AS raw_observed_uo_sum_ml_6h
FROM p_target t JOIN mimiciv_icu.outputevents o
 ON o.stay_id=t.stay_id AND o.subject_id=t.subject_id AND o.hadm_id=t.hadm_id
 AND o.charttime>t.landmark_time-interval '6h' AND o.charttime<=t.landmark_time
 AND o.storetime<=t.landmark_time
WHERE o.itemid IN(226559,226560,226561,226584,226563,226564,226565,226567,226557,226558)
 AND lower(regexp_replace(coalesce(trim(o.valueuom),''),'\s','','g'))='ml'
 AND o.value>=0 AND o.value::text!~*'^(nan|[-+]?infinity)$'
GROUP BY t.stay_id;
CREATE UNIQUE INDEX ON uo_raw6h(stay_id);
CREATE TEMP TABLE measurement_supplement AS
SELECT t.stay_id,
 extract(epoch FROM(a.charttime-t.intime))/3600 AS first_any_sample_hour,
 a.specimen_id AS first_any_specimen_id,
 g.specimen_id AS selected_protein_specimen_id,g.protein_band AS selected_protein_band,
 extract(epoch FROM(g.charttime-t.intime))/3600 AS selected_protein_sample_hour,
 (g.specimen_id IS NOT NULL AND g.specimen_id<>a.specimen_id) AS first_specimen_taken_over,
 CASE WHEN g.specimen_id IS NOT NULL AND g.specimen_id<>a.specimen_id
 THEN extract(epoch FROM(g.charttime-a.charttime))/3600 END AS takeover_delay_hours,
 coalesce(z.raw_rows,0) AS protein_raw_rows,
 coalesce(z.nonnull_specimen_rows,0) AS protein_nonnull_specimen_rows,
 coalesce(z.available_rows,0) AS protein_available_rows,
 coalesce(z.storetime_null_rows,0) AS protein_storetime_null_rows,
 coalesce(z.report_after_L_rows,0) AS protein_report_after_L_rows,
 coalesce(z.interpretable_unavailable_rows,0) AS interpretable_unavailable_rows,
 extract(epoch FROM(f.first_scuf_after_L-t.landmark_time))/3600 AS first_SCUF_after_L_hours,
 u.raw_observed_uo_sum_ml_6h,
 CASE WHEN g.stay_id IS NOT NULL THEN 'measured'
 WHEN coalesce(z.raw_rows,0)=0 THEN 'no_51492_record'
 WHEN coalesce(z.nonnull_specimen_rows,0)=0 THEN 'all_specimen_id_missing'
 WHEN coalesce(z.available_rows,0)=0 THEN 'no_available_identifiable_specimen_by_L'
 ELSE 'available_canonical_none_interpretable' END AS measurement_flow_reason
FROM p_target t LEFT JOIN p_first_any a USING(stay_id)
LEFT JOIN p_first_good g USING(stay_id)
LEFT JOIN scuf_after_L f USING(stay_id) LEFT JOIN uo_raw6h u USING(stay_id)
LEFT JOIN LATERAL(
 SELECT count(*) AS raw_rows,
 count(*) FILTER(WHERE specimen_id IS NOT NULL) AS nonnull_specimen_rows,
 count(*) FILTER(WHERE specimen_id IS NOT NULL AND storetime<=landmark_time) AS available_rows,
 count(*) FILTER(WHERE storetime IS NULL) AS storetime_null_rows,
 count(*) FILTER(WHERE storetime>landmark_time) AS report_after_L_rows,
 count(*) FILTER(WHERE specimen_id IS NOT NULL AND protein_band IS NOT NULL
 AND(storetime IS NULL OR storetime>landmark_time)) AS interpretable_unavailable_rows
 FROM p_lab b WHERE b.stay_id=t.stay_id
) z ON true;
CREATE UNIQUE INDEX ON measurement_supplement(stay_id);
DO $$ BEGIN
 IF EXISTS(SELECT 1 FROM memberships m JOIN measurement_supplement s USING(stay_id)
 WHERE m.expected_band IS DISTINCT FROM s.selected_protein_band
 OR m.expected_specimen_id IS DISTINCT FROM s.selected_protein_specimen_id)
 THEN RAISE EXCEPTION 'QC_MEASUREMENT_SELECTION_CHANGED'; END IF;
 IF EXISTS(SELECT 1 FROM measurement_supplement WHERE takeover_delay_hours<0)
 THEN RAISE EXCEPTION 'QC_MEASUREMENT_TAKEOVER_TIME'; END IF;
END $$;
CREATE TEMP TABLE measurement_summary(section text,cohort text,stratum text,metric text,
 n bigint,denominator bigint,value_text text,note text);
INSERT INTO measurement_summary
SELECT 'measurement_flow',m.cohort,s.measurement_flow_reason,'people',count(*),
 (SELECT count(*) FROM memberships x WHERE x.cohort=m.cohort),NULL,
 'Mutually exclusive patient flow; no_51492_record does not prove test not ordered'
FROM memberships m JOIN measurement_supplement s USING(stay_id)
GROUP BY m.cohort,s.measurement_flow_reason;
INSERT INTO measurement_summary
SELECT 'measurement_process',m.cohort,'all',v.metric,count(*) FILTER(WHERE v.met),count(*),NULL,v.note
FROM memberships m JOIN measurement_supplement s USING(stay_id)
CROSS JOIN LATERAL(VALUES
 ('first_specimen_taken_over',coalesce(s.first_specimen_taken_over,false),'Original non-null specimen comparison; may include zero-hour specimen ties'),
 ('first_raw_specimen_id_missing_but_later_usable',s.first_any_specimen_id IS NULL AND s.selected_protein_band IS NOT NULL,'Unknown first specimen identity; separate from proven different-specimen takeover'),
 ('protein_storetime_missing',s.protein_storetime_null_rows>0,'Overlapping patient quality flags'),
 ('protein_report_after_L',s.protein_report_after_L_rows>0,'Overlapping patient quality flags'),
 ('unmeasured_with_potentially_interpretable_unavailable_record',s.selected_protein_band IS NULL AND s.interpretable_unavailable_rows>0,'Availability-associated; not necessarily the sole cause'),
 ('availability_only_at_record_level',s.measurement_flow_reason='no_available_identifiable_specimen_by_L' AND s.interpretable_unavailable_rows>0,'At least one parsed valid identifiable record; none of the identifiable records available by L')
)v(metric,met,note) GROUP BY m.cohort,v.metric,v.note;
INSERT INTO measurement_summary
SELECT 'sensitivity_input_coverage',m.cohort,'all',v.metric,count(*) FILTER(WHERE v.met),count(*),NULL,
 'Local supplementary inputs only; no new model or replacement of locked M0'
FROM memberships m JOIN measurement_supplement s USING(stay_id)
CROSS JOIN LATERAL(VALUES
 ('SCUF_first_observed_after_L',s.first_SCUF_after_L_hours IS NOT NULL),
 ('raw_recorded_UO_volume_6h_observed',s.raw_observed_uo_sum_ml_6h IS NOT NULL)
)v(metric,met) GROUP BY m.cohort,v.metric;
INSERT INTO measurement_summary
SELECT 'takeover_delay',m.cohort,'first_specimen_taken_over','hours_Q25_median_Q75',count(*),NULL,
 CASE WHEN count(*)>=10 THEN round((percentile_cont(0.25) WITHIN GROUP(ORDER BY s.takeover_delay_hours))::numeric,2)::text||'|'||
 round((percentile_cont(0.50) WITHIN GROUP(ORDER BY s.takeover_delay_hours))::numeric,2)::text||'|'||
 round((percentile_cont(0.75) WITHIN GROUP(ORDER BY s.takeover_delay_hours))::numeric,2)::text
 ELSE 'suppressed_n_lt10' END,'Hours; type-7 continuous percentiles; 2 decimals; only takeover cases'
FROM memberships m JOIN measurement_supplement s USING(stay_id)
WHERE s.first_specimen_taken_over GROUP BY m.cohort;
INSERT INTO measurement_summary
SELECT 'canonical_protein_specimens',m.cohort,coalesce(p.loss_reason,'interpretable'),
 'specimens',count(*),NULL,NULL,'Hierarchical protein parser loss reasons; stay/specimen deduplicated'
FROM memberships m JOIN p_canonical p USING(stay_id) GROUP BY m.cohort,p.loss_reason;

CREATE TEMP TABLE specimen_loss_class AS
SELECT m.cohort,g.stay_id,it.itemid,v.labevent_id,v.raw_value,v.valuenum,v.unit,
 v.clean_numeric,v.numeric_candidate,v.numeric_conflict,v.unparsed_value,v.invalid,
 v.exact_trace_comment,v.exact_negative_comment,v.hpf_unit,
 CASE WHEN v.labevent_id IS NULL THEN
 CASE WHEN EXISTS(SELECT 1 FROM cov_specimen_all z WHERE z.stay_id=g.stay_id AND z.itemid=it.itemid)
 THEN 'only_unavailable_by_L' ELSE 'no_matching_record' END
 WHEN v.clean_numeric IS NOT NULL THEN 'usable_numeric'
 WHEN v.invalid THEN 'invalid_or_cancelled'
 WHEN v.numeric_conflict THEN 'numeric_conflict'
 WHEN v.unparsed_value AND v.valuenum IS NOT NULL THEN 'nonempty_unparsed_with_valuenum'
 WHEN v.unparsed_value AND v.valuenum IS NULL THEN 'nonempty_unparsed_without_valuenum'
 WHEN v.exact_trace_comment OR(v.exact_negative_comment AND v.numeric_candidate<>0) THEN 'exact_token_semantics'
 WHEN v.numeric_candidate IS NULL THEN 'missing_numeric'
 WHEN v.numeric_candidate::text~*'^(nan|[-+]?infinity)$' THEN 'nonfinite_numeric'
 WHEN it.itemid IN(51493,51516) AND NOT v.hpf_unit THEN 'nonhpf_unit'
 ELSE 'range_or_sign_invalid' END AS missing_reason
FROM memberships m JOIN p_first_good g USING(stay_id) JOIN p_target t USING(stay_id)
CROSS JOIN(VALUES(51498),(51493),(51516)) it(itemid)
LEFT JOIN cov_specimen_values v ON v.stay_id=g.stay_id AND v.itemid=it.itemid;
INSERT INTO measurement_summary
SELECT 'specimen_item_flow',cohort,itemid::text||'|'||missing_reason,'selected_specimen_items',count(*),
 (SELECT count(*) FROM memberships x JOIN p_first_good y USING(stay_id) WHERE x.cohort=s.cohort),NULL,
 'Mutually exclusive item reasons, including usable; three independent per-item partitions'
FROM specimen_loss_class s GROUP BY cohort,itemid,missing_reason;
INSERT INTO measurement_summary
SELECT 'specimen_five_flags',m.cohort,v.itemid::text||'|'||q.metric,'selected_available_items',
 count(*) FILTER(WHERE q.met),count(*),NULL,'Overlapping quality flags; not an additive partition'
FROM memberships m JOIN cov_specimen_values v USING(stay_id)
CROSS JOIN LATERAL(VALUES
 ('numeric_conflict',v.numeric_conflict),
 ('nonempty_unparsed_with_valuenum',v.unparsed_value AND v.valuenum IS NOT NULL),
 ('nonempty_unparsed_without_valuenum',v.unparsed_value AND v.valuenum IS NULL),
 ('nonhpf_unit',v.itemid IN(51493,51516) AND NOT v.hpf_unit),
 ('invalid_or_cancelled',v.invalid)
)q(metric,met) GROUP BY m.cohort,v.itemid,q.metric;
INSERT INTO measurement_summary
SELECT 'specimen_units',m.cohort,v.itemid::text||'|'||
 CASE WHEN v.itemid=51498 AND v.unit='' THEN 'dimensionless'
 WHEN v.unit IN('/hpf','#/hpf','cells/hpf','rbc/hpf','wbc/hpf','hpf','') THEN coalesce(nullif(v.unit,''),'missing') ELSE 'other' END,
 'selected_available_items',count(*),NULL,NULL,'51498 empty unit is labelled dimensionless; no statement of assay equivalence'
FROM memberships m JOIN cov_specimen_latest v USING(stay_id) GROUP BY 2,3;

INSERT INTO measurement_summary
SELECT 'pre_RRT_SCr',m.cohort,'all',q.metric,count(*) FILTER(WHERE q.met),count(*),NULL,
 'T0-L samples before first supported treatment with report available by L; row flags overlap'
FROM memberships m JOIN extra_scr_checked v USING(stay_id)
CROSS JOIN LATERAL(VALUES
 ('nonempty_unparsed_with_valuenum',v.raw_value<>'' AND v.numeric_text IS NULL AND v.valuenum IS NOT NULL),
 ('nonempty_unparsed_without_valuenum',v.raw_value<>'' AND v.numeric_text IS NULL AND v.valuenum IS NULL),
 ('numeric_conflict',coalesce(v.numeric_conflict,false)),
 ('invalid_or_cancelled',v.invalid),('unsupported_unit',v.unit<>'mg/dl')
)q(metric,met) GROUP BY m.cohort,q.metric;
INSERT INTO measurement_summary
SELECT 'pre_RRT_SCr_patient_flow',m.cohort,q.reason,'people',count(*),NULL,NULL,
 'Mutually exclusive among first supported treatment in (T0,L]; no row may mean no eligible observed record'
FROM memberships m JOIN extra_keys k USING(stay_id)
LEFT JOIN extra_scr_last z USING(stay_id)
CROSS JOIN LATERAL(SELECT CASE WHEN z.stay_id IS NOT NULL THEN 'usable_pretherapy_last'
 WHEN EXISTS(SELECT 1 FROM extra_scr_raw w WHERE w.stay_id=k.stay_id) THEN 'observed_but_no_valid_pretherapy_SCr'
 ELSE 'no_available_pretherapy_SCr_record' END AS reason) q
WHERE k.first_tx_time>k.source_t0 AND k.first_tx_time<=k.source_l GROUP BY m.cohort,q.reason;
DO $$ BEGIN
 IF EXISTS(SELECT cohort FROM specimen_loss_class GROUP BY cohort HAVING count(*)<>
  3*(SELECT count(*) FROM memberships m JOIN p_first_good g USING(stay_id) WHERE m.cohort=specimen_loss_class.cohort))
 THEN RAISE EXCEPTION 'QC_SPECIMEN_ITEM_PARTITION'; END IF;
END $$;
COMMIT;
