BEGIN;
SET LOCAL statement_timeout='30min';
SET LOCAL lock_timeout='5s';
SET LOCAL enable_nestloop=off;
SET LOCAL work_mem='128MB';

CREATE TEMP TABLE renal_result (
 section text, window_h integer, stratum text, metric text,
 n bigint, denominator bigint, value_text text, note text
);
CREATE OR REPLACE FUNCTION pg_temp.renal_band(n bigint) RETURNS text
LANGUAGE sql IMMUTABLE AS $$
SELECT CASE WHEN n IS NULL THEN NULL WHEN n=0 THEN '0' WHEN n<10 THEN '<10'
ELSE (n/10*10)::text||'-'||(n/10*10+9)::text END $$;
CREATE TEMP TABLE renal_icu AS
SELECT i.*,a.admittime,a.dischtime,a.deathtime,a.hospital_expire_flag,
 p.gender,p.dod,p.anchor_age+extract(year FROM a.admittime)::integer-p.anchor_year AS age,
 row_number() OVER(PARTITION BY i.subject_id ORDER BY i.intime,i.stay_id) AS icu_order
FROM mimiciv_icu.icustays i
JOIN mimiciv_hosp.admissions a USING(hadm_id,subject_id)
JOIN mimiciv_hosp.patients p USING(subject_id);
CREATE INDEX ON renal_icu(stay_id);
CREATE INDEX ON renal_icu(hadm_id);
ANALYZE renal_icu;

CREATE TEMP TABLE renal_s AS
SELECT DISTINCT ON(i.stay_id) i.*,s.suspected_infection_time AS s_time,
 s.sofa_time,s.sofa_score
FROM renal_icu i JOIN mimiciv_derived.sepsis3 s USING(stay_id)
WHERE i.age>=18 AND s.sepsis3 IS TRUE
 AND s.suspected_infection_time>=i.intime-interval '24h'
 AND s.suspected_infection_time<=i.intime+interval '24h'
ORDER BY i.stay_id,s.suspected_infection_time,s.sofa_time;
CREATE INDEX ON renal_s(stay_id);
ANALYZE renal_s;

CREATE TEMP TABLE renal_k AS
SELECT s.stay_id,k.charttime,k.aki_stage,k.aki_stage_creat,k.aki_stage_uo,k.aki_stage_crrt,
 k.creat,k.creat_low_past_48hr,k.creat_low_past_7day,k.uo_rt_6hr
FROM renal_s s JOIN mimiciv_derived.kdigo_stages k USING(stay_id)
WHERE k.charttime>=s.intime-interval '24h' AND k.charttime<least(s.outtime,s.intime+interval '96h');
CREATE INDEX ON renal_k(stay_id,charttime);
ANALYZE renal_k;
CREATE TEMP TABLE renal_a AS
SELECT s.stay_id,
 min(k.charttime) FILTER(WHERE k.charttime>=s.intime AND k.charttime<=s.intime+interval '48h'
   AND greatest(coalesce(k.aki_stage_creat,0),coalesce(k.aki_stage_uo,0))>0) AS a_time,
 bool_or(k.charttime<s.intime AND greatest(coalesce(k.aki_stage_creat,0),coalesce(k.aki_stage_uo,0))>0) AS pre_aki,
 max(k.aki_stage) FILTER(WHERE k.charttime>=s.intime AND k.charttime<=s.intime+interval '48h') AS max_full_stage_48_audit_only,
 count(k.charttime) FILTER(WHERE k.charttime>=s.intime AND k.charttime<=s.intime+interval '48h') AS kdigo_rows_48,
 bool_or(k.charttime>=s.intime AND k.charttime<=s.intime+interval '48h' AND k.aki_stage_crrt>0) AS crrt_48
FROM renal_s s LEFT JOIN renal_k k USING(stay_id) GROUP BY s.stay_id;
CREATE TEMP TABLE renal_c AS
SELECT s.*,a.a_time,a.pre_aki,a.max_full_stage_48_audit_only,
 row_number() OVER(PARTITION BY s.subject_id ORDER BY s.intime,s.stay_id) AS candidate_order
FROM renal_s s JOIN renal_a a USING(stay_id) WHERE a.a_time IS NOT NULL;
CREATE UNIQUE INDEX ON renal_c(stay_id);
CREATE INDEX ON renal_c(subject_id);
CREATE INDEX ON renal_c(hadm_id);
ANALYZE renal_c;


CREATE TEMP TABLE p_prior AS
SELECT c.stay_id,
 bool_or((d.icd_version=9 AND replace(d.icd_code,'.','')='5856') OR (d.icd_version=10 AND replace(d.icd_code,'.','')='N186')) AS esrd,
 bool_or((d.icd_version=9 AND replace(d.icd_code,'.','')='V4511') OR (d.icd_version=10 AND replace(d.icd_code,'.','')='Z992')) AS dependence
FROM renal_c c JOIN mimiciv_hosp.admissions a ON a.subject_id=c.subject_id AND a.dischtime<c.admittime
JOIN mimiciv_hosp.diagnoses_icd d ON d.hadm_id=a.hadm_id GROUP BY c.stay_id;
CREATE UNIQUE INDEX ON p_prior(stay_id);
CREATE TEMP TABLE p_clinical AS
SELECT c.*,coalesce(p.esrd AND p.dependence,false) AS maintenance,
 EXISTS(SELECT 1 FROM renal_k k WHERE k.stay_id=c.stay_id AND k.charttime<c.intime AND k.charttime<c.s_time
 AND greatest(coalesce(k.aki_stage_creat,0),coalesce(k.aki_stage_uo,0))>0) AS preS_positive
FROM renal_c c LEFT JOIN p_prior p USING(stay_id) WHERE c.a_time>=c.s_time;
CREATE TEMP TABLE p_eligible AS SELECT * FROM p_clinical WHERE NOT preS_positive;
CREATE TEMP TABLE p_rank_old AS SELECT *,row_number() OVER(PARTITION BY subject_id ORDER BY intime,stay_id) AS rn FROM p_eligible;
CREATE TEMP TABLE p_rank_new AS SELECT *,row_number() OVER(PARTITION BY subject_id ORDER BY intime,stay_id) AS rn FROM p_eligible WHERE NOT maintenance;
CREATE TEMP TABLE p_first_old AS SELECT * FROM p_rank_old WHERE rn=1 AND NOT maintenance;
CREATE TEMP TABLE p_first AS SELECT * FROM p_rank_new WHERE rn=1;
CREATE UNIQUE INDEX ON p_first(stay_id); ANALYZE p_first;
CREATE TEMP TABLE p_target AS SELECT *,intime+interval '72h' AS landmark_time FROM p_first
WHERE outtime>intime+interval '72h' AND (deathtime IS NULL OR deathtime>intime+interval '72h')
AND (dod IS NULL OR deathtime IS NOT NULL OR dod>(intime+interval '72h')::date);
CREATE UNIQUE INDEX ON p_target(stay_id); ANALYZE p_target;

CREATE TEMP TABLE p_events AS
SELECT i.hadm_id,r.stay_id,r.charttime AS event_time,
 CASE WHEN r.dialysis_active=1 AND upper(trim(r.dialysis_type))='SCUF' THEN 'scuf'
 WHEN r.dialysis_active=1 AND upper(trim(r.dialysis_type)) IN ('IHD','HD','CRRT','CVVH','CVVHD','CVVHDF','CAVH','CAVHD','PERITONEAL','PD') THEN 'main'
 WHEN r.dialysis_active=1 THEN 'unknown_active' ELSE 'presence' END AS category,
 r.dialysis_active=1 AS old_main,r.dialysis_active=1 AS active_record
FROM mimiciv_derived.rrt r JOIN renal_icu i USING(stay_id)
WHERE (r.dialysis_active=1 OR r.dialysis_present=1) AND EXISTS(SELECT 1 FROM p_target t WHERE t.hadm_id=i.hadm_id)
UNION ALL SELECT p.hadm_id,p.stay_id,p.starttime,
 CASE WHEN p.itemid=224270 THEN 'catheter' WHEN p.itemid=225955 THEN 'scuf' ELSE 'main' END,p.itemid<>224270,false
FROM mimiciv_icu.procedureevents p WHERE p.itemid IN (224270,225441,225802,225803,225805,225809,225955)
AND coalesce(lower(p.statusdescription),'')!~'(cancel|rewrit)' AND EXISTS(SELECT 1 FROM p_target t WHERE t.hadm_id=p.hadm_id);
CREATE INDEX ON p_events(hadm_id,event_time); ANALYZE p_events;
CREATE TEMP TABLE p_r AS
SELECT t.stay_id,
 min(e.event_time) FILTER(WHERE e.category='main') AS first_supported_tx,
 coalesce(bool_or(e.active_record AND e.event_time=t.landmark_time),false) AS any_active_at_L,
 coalesce(bool_or(e.old_main AND e.event_time<=t.landmark_time),false) AS old_pre,
 coalesce(bool_or(e.category='main' AND e.event_time<=t.landmark_time),false) AS main_pre,
 coalesce(bool_or(e.category='scuf' AND e.event_time<=t.landmark_time),false) AS scuf_pre,
 coalesce(bool_or(e.category='unknown_active' AND e.event_time<=t.landmark_time),false) AS unknown_pre,
 coalesce(bool_or(e.category='catheter' AND e.event_time<=t.landmark_time),false) AS catheter_pre,
 coalesce(bool_or(e.category='presence' AND e.event_time<=t.landmark_time),false) AS presence_pre,
 coalesce(bool_or(e.category='main' AND e.stay_id<>t.stay_id AND e.event_time<=t.landmark_time),false) AS cross_stay_main_pre,
 min(e.event_time) FILTER(WHERE e.category='main' AND e.stay_id=t.stay_id AND e.event_time>t.landmark_time AND e.event_time<=least(t.landmark_time+interval '7 days',t.outtime) AND (t.deathtime IS NULL OR e.event_time<t.deathtime)) AS new_index_tx,
 min(e.event_time) FILTER(WHERE e.category='main' AND e.stay_id=t.stay_id AND e.event_time>t.landmark_time AND e.event_time<=least(t.landmark_time+interval '7 days',t.outtime)) AS post_tx_before_death_adjudication,
 coalesce(bool_or(e.old_main AND e.stay_id<>t.stay_id AND e.event_time>t.landmark_time),false) AS cross_stay_old_future
FROM p_target t LEFT JOIN p_events e ON e.hadm_id=t.hadm_id AND e.event_time>=t.admittime AND e.event_time<=least(t.dischtime,t.landmark_time+interval '7 days') GROUP BY t.stay_id;

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
CREATE TEMP TABLE cov_anchor AS
SELECT t.*
FROM p_target t;
CREATE UNIQUE INDEX ON cov_anchor(stay_id);
ANALYZE cov_anchor;

CREATE TEMP TABLE cov_prior AS
SELECT t.stay_id,coalesce(bool_or(d.hadm_id IS NOT NULL),false) AS prior_ckd_recorded
FROM cov_anchor t
LEFT JOIN mimiciv_hosp.admissions a ON a.subject_id=t.subject_id AND a.dischtime<t.admittime
LEFT JOIN mimiciv_hosp.diagnoses_icd d ON d.hadm_id=a.hadm_id
 AND ((d.icd_version=9 AND replace(trim(d.icd_code),'.','') LIKE '585%')
   OR (d.icd_version=10 AND upper(replace(trim(d.icd_code),'.','')) LIKE 'N18%'))
GROUP BY t.stay_id;
CREATE UNIQUE INDEX ON cov_prior(stay_id);

CREATE TEMP TABLE cov_sofa AS
SELECT s.* FROM mimiciv_derived.first_day_sofa s
JOIN cov_anchor t ON t.stay_id=s.stay_id;
CREATE UNIQUE INDEX ON cov_sofa(stay_id);

CREATE TEMP TABLE cov_stage AS
SELECT t.stay_id,
 max(greatest(coalesce(k.aki_stage_creat,0),coalesce(k.aki_stage_uo,0))) AS max_cruo_stage_to_L,
 count(k.charttime) AS stage_observed_rows,
 count(*) FILTER(WHERE k.charttime IS NOT NULL AND k.aki_stage_creat IS NULL AND k.aki_stage_uo IS NULL) AS stage_both_components_null_rows,
 count(*) FILTER(WHERE k.aki_stage_creat NOT BETWEEN 0 AND 3 OR k.aki_stage_uo NOT BETWEEN 0 AND 3) AS stage_invalid_rows
FROM cov_anchor t LEFT JOIN renal_k k ON k.stay_id=t.stay_id
 AND k.charttime>=t.intime AND k.charttime<=t.landmark_time
GROUP BY t.stay_id;
CREATE UNIQUE INDEX ON cov_stage(stay_id);

CREATE TEMP TABLE cov_scr_raw AS
SELECT t.stay_id,t.intime,t.landmark_time,l.labevent_id,l.specimen_id,l.charttime,l.storetime,l.valuenum,
 lower(regexp_replace(coalesce(trim(l.valueuom),''),'\s','','g')) AS unit,
 lower(trim(coalesce(l.value,''))) AS raw_value,
 coalesce(l.comments,'')~*'(cancel|not (done|performed)|insufficient|unsatisfactory|unable|disregard)' AS invalid,
 upper(trim(coalesce(l.comments,'')))~'^(NEG|NEGATIVE)[.!]?$' AS exact_negative_comment,
 upper(trim(coalesce(l.comments,'')))~'^(TRACE|TR)[.!]?$' AS exact_trace_comment
FROM cov_anchor t JOIN mimiciv_hosp.labevents l
 ON l.subject_id=t.subject_id AND (l.hadm_id=t.hadm_id OR l.hadm_id IS NULL)
 AND l.charttime>=greatest(t.intime,t.admittime) AND l.charttime<t.landmark_time
WHERE l.itemid=50912;
CREATE INDEX ON cov_scr_raw(stay_id,charttime);
ANALYZE cov_scr_raw;
CREATE TEMP TABLE cov_scr_parsed AS
SELECT r.*,
 CASE WHEN raw_value~'^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][+-]?[0-9]+)?$'
 THEN raw_value::double precision END AS numeric_text
FROM cov_scr_raw r;
CREATE TEMP TABLE cov_scr_class AS
SELECT r.*,
 CASE WHEN numeric_text IS NOT NULL THEN numeric_text WHEN raw_value='' THEN valuenum END AS numeric_candidate,
 CASE WHEN numeric_text IS NOT NULL AND valuenum IS NOT NULL
 THEN abs(numeric_text-valuenum)>1e-6*greatest(1,abs(numeric_text)) ELSE false END AS numeric_conflict,
 raw_value<>'' AND numeric_text IS NULL AS unparsed_value
FROM cov_scr_parsed r;
CREATE TEMP TABLE cov_scr_valid AS
SELECT *,numeric_candidate AS scr_mgdl
FROM cov_scr_class
WHERE unit='mg/dl' AND NOT invalid AND NOT numeric_conflict AND NOT unparsed_value
 AND numeric_candidate>0 AND numeric_candidate::text!~*'^(nan|[-+]?infinity)$'
 AND NOT exact_negative_comment AND NOT exact_trace_comment
 AND storetime<=landmark_time;
CREATE TEMP TABLE cov_scr_first AS
SELECT DISTINCT ON(stay_id) * FROM cov_scr_valid
WHERE charttime<intime+interval '24h'
ORDER BY stay_id,charttime,storetime,labevent_id;
CREATE UNIQUE INDEX ON cov_scr_first(stay_id);
CREATE TEMP TABLE cov_scr_last AS
SELECT DISTINCT ON(stay_id) * FROM cov_scr_valid
ORDER BY stay_id,charttime DESC,storetime DESC,labevent_id DESC;
CREATE UNIQUE INDEX ON cov_scr_last(stay_id);
CREATE TEMP TABLE cov_scr_quality AS
SELECT stay_id,count(*) AS scr_raw_rows,
 count(*) FILTER(WHERE storetime IS NULL) AS scr_storetime_null_rows,
 count(*) FILTER(WHERE storetime>landmark_time) AS scr_unavailable_by_L_rows,
 count(*) FILTER(WHERE numeric_conflict) AS scr_conflict_rows,
 count(*) FILTER(WHERE unparsed_value AND valuenum IS NOT NULL) AS scr_unparsed_value_with_valuenum_rows,
 count(*) FILTER(WHERE unit<>'mg/dl') AS scr_unsupported_unit_rows,
 count(*) FILTER(WHERE numeric_candidate<=0 OR numeric_candidate::text~*'^(nan|[-+]?infinity)$') AS scr_nonpositive_or_nonfinite_rows,
 count(*) FILTER(WHERE invalid) AS scr_invalid_rows
FROM cov_scr_class GROUP BY stay_id;
CREATE UNIQUE INDEX ON cov_scr_quality(stay_id);

CREATE TEMP TABLE cov_uo_rows AS
SELECT t.stay_id,u.charttime,u.weight,u.uo_rt_6hr,u.uo_tm_6hr
FROM cov_anchor t JOIN mimiciv_derived.kdigo_uo u ON u.stay_id=t.stay_id
 AND u.charttime>t.landmark_time-interval '1h' AND u.charttime<=t.landmark_time;
CREATE UNIQUE INDEX ON cov_uo_rows(stay_id,charttime);
CREATE TEMP TABLE cov_uo_latest AS
SELECT DISTINCT ON(stay_id) *,
 CASE WHEN weight>0 AND weight::text!~*'^(nan|[-+]?infinity)$'
 AND uo_tm_6hr>=6 AND uo_tm_6hr<12
 AND uo_rt_6hr>=0 AND uo_rt_6hr::text!~*'^(nan|[-+]?infinity)$'
 THEN uo_rt_6hr END AS uo_rate_6hr_valid
FROM cov_uo_rows ORDER BY stay_id,charttime DESC;
CREATE UNIQUE INDEX ON cov_uo_latest(stay_id);

CREATE TEMP TABLE cov_uo_output AS
SELECT t.stay_id,o.charttime,o.value,
 floor(extract(epoch FROM(t.landmark_time-o.charttime))/3600)::integer AS right_closed_hour
FROM cov_anchor t JOIN mimiciv_icu.outputevents o
 ON o.stay_id=t.stay_id AND o.subject_id=t.subject_id AND o.hadm_id=t.hadm_id
 AND o.charttime>t.landmark_time-interval '6h' AND o.charttime<=t.landmark_time
 AND o.storetime<=t.landmark_time
WHERE o.itemid IN(226559,226560,226561,226584,226563,226564,226565,226567,226557,226558)
 AND lower(regexp_replace(coalesce(trim(o.valueuom),''),'\s','','g'))='ml'
 AND o.value>=0 AND o.value::text!~*'^(nan|[-+]?infinity)$';
CREATE TEMP TABLE cov_uo_coverage AS
SELECT stay_id,count(DISTINCT right_closed_hour) AS uo_covered_hours_6h,
 count(*) AS uo_valid_record_rows_6h,
 count(*) FILTER(WHERE value=0) AS uo_explicit_zero_record_rows_6h
FROM cov_uo_output GROUP BY stay_id;
CREATE UNIQUE INDEX ON cov_uo_coverage(stay_id);

SET LOCAL enable_nestloop=on;
CREATE TEMP TABLE cov_specimen_raw AS
SELECT g.stay_id,l.labevent_id,l.itemid,l.charttime,l.storetime,l.valuenum,
 lower(regexp_replace(coalesce(trim(l.valueuom),''),'\s','','g')) AS unit,
 lower(trim(coalesce(l.value,''))) AS raw_value,
 coalesce(l.comments,'')~*'(cancel|not (done|performed)|insufficient|unsatisfactory|unable|disregard)' AS invalid,
 upper(trim(coalesce(l.comments,'')))~'^(NEG|NEGATIVE)[.!]?$' AS exact_negative_comment,
 upper(trim(coalesce(l.comments,'')))~'^(TRACE|TR)[.!]?$' AS exact_trace_comment
FROM p_first_good g JOIN cov_anchor t USING(stay_id)
JOIN mimiciv_hosp.labevents l ON l.subject_id=t.subject_id AND l.specimen_id=g.specimen_id
 AND (l.hadm_id=t.hadm_id OR l.hadm_id IS NULL) AND l.storetime<=t.landmark_time
WHERE l.itemid IN(51498,51493,51516);
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

CREATE TEMP TABLE covariates AS
SELECT t.stay_id,c.prior_ckd_recorded,
 s.respiration AS sofa_respiration,s.coagulation AS sofa_coagulation,s.liver AS sofa_liver,
 s.cardiovascular AS sofa_cardiovascular,s.cns AS sofa_cns,
 coalesce(s.respiration,0)+coalesce(s.coagulation,0)+coalesce(s.liver,0)+coalesce(s.cardiovascular,0)+coalesce(s.cns,0) AS nonrenal_sofa,
 (s.respiration IS NULL)::integer+(s.coagulation IS NULL)::integer+(s.liver IS NULL)::integer+
 (s.cardiovascular IS NULL)::integer+(s.cns IS NULL)::integer AS nonrenal_sofa_missing_components,
 (s.stay_id IS NOT NULL) AS sofa_source_row_present,
 CASE WHEN s.stay_id IS NOT NULL THEN s.subject_id IS NOT DISTINCT FROM t.subject_id
 AND s.hadm_id IS NOT DISTINCT FROM t.hadm_id END AS sofa_source_identifiers_match,
 k.max_cruo_stage_to_L,k.stage_observed_rows,k.stage_both_components_null_rows,k.stage_invalid_rows,
 f.scr_mgdl AS first_scr_0_24_mgdl,z.scr_mgdl AS last_scr_0_L_mgdl,
 ln(z.scr_mgdl)/ln(2.0) AS log2_last_scr,
 (ln(z.scr_mgdl)-ln(f.scr_mgdl))/ln(2.0) AS log2_scr_change,
 f.charttime AS first_scr_charttime,z.charttime AS last_scr_charttime,
 f.storetime AS first_scr_storetime,z.storetime AS last_scr_storetime,
 q.scr_raw_rows,q.scr_storetime_null_rows,q.scr_unavailable_by_L_rows,q.scr_conflict_rows,
 q.scr_unparsed_value_with_valuenum_rows,q.scr_unsupported_unit_rows,q.scr_nonpositive_or_nonfinite_rows,q.scr_invalid_rows,
 u.uo_rate_6hr_valid,u.charttime AS uo_rate_charttime,u.weight AS uo_rate_weight,u.uo_tm_6hr AS uo_rate_tm6,
 coalesce(h.uo_covered_hours_6h,0) AS uo_covered_hours_6h,
 coalesce(h.uo_covered_hours_6h,0)/6.0 AS uo_record_hour_coverage,
 coalesce(h.uo_valid_record_rows_6h,0) AS uo_valid_record_rows_6h,
 coalesce(h.uo_explicit_zero_record_rows_6h,0) AS uo_explicit_zero_record_rows_6h,
 CASE WHEN r.main_pre THEN 'supported_main'
 WHEN r.unknown_pre OR r.presence_pre OR r.catheter_pre OR r.scuf_pre THEN 'other_evidence_no_main'
 ELSE 'no_supported_evidence' END AS rrt_state_L,
 r.main_pre,r.unknown_pre,r.presence_pre,r.catheter_pre,r.scuf_pre,
 e.urine_specific_gravity,e.urine_rbc_hpf,e.urine_wbc_hpf,
 ln(1+e.urine_rbc_hpf) AS log1p_urine_rbc_hpf,ln(1+e.urine_wbc_hpf) AS log1p_urine_wbc_hpf,
 e.specimen_numeric_conflict_items,e.specimen_unparsed_value_with_valuenum_items,
 e.specimen_nonhpf_unit_items,e.specimen_invalid_items
FROM cov_anchor t
JOIN cov_prior c USING(stay_id)
JOIN cov_stage k USING(stay_id)
JOIN p_r r USING(stay_id)
LEFT JOIN cov_sofa s ON s.stay_id=t.stay_id
LEFT JOIN cov_scr_first f ON f.stay_id=t.stay_id
LEFT JOIN cov_scr_last z ON z.stay_id=t.stay_id
LEFT JOIN cov_scr_quality q ON q.stay_id=t.stay_id
LEFT JOIN cov_uo_latest u ON u.stay_id=t.stay_id
LEFT JOIN cov_uo_coverage h ON h.stay_id=t.stay_id
LEFT JOIN cov_specimen_pivot e ON e.stay_id=t.stay_id;
CREATE UNIQUE INDEX ON covariates(stay_id);
ANALYZE covariates;

CREATE TEMP TABLE extraction_cohort AS
SELECT t.subject_id,t.hadm_id,t.stay_id,t.age,t.gender,
 extract(epoch FROM(t.s_time-t.intime))/3600 AS s_hour,
 extract(epoch FROM(t.a_time-t.intime))/3600 AS a_hour,
 extract(epoch FROM(t.outtime-t.landmark_time))/3600 AS icu_exit_after_L_hours,
 extract(epoch FROM(t.deathtime-t.landmark_time))/3600 AS precise_death_after_L_hours,
 CASE WHEN t.deathtime IS NULL AND t.dod IS NOT NULL THEN extract(epoch FROM(t.dod::timestamp-t.landmark_time))/3600 END AS death_date_lower_after_L_hours,
 CASE WHEN t.deathtime IS NULL AND t.dod IS NOT NULL THEN extract(epoch FROM(t.dod::timestamp+interval '1 day'-t.landmark_time))/3600 END AS death_date_upper_after_L_hours,
 CASE WHEN(t.deathtime>t.landmark_time AND t.deathtime<=t.intime+interval '28 days')
 OR(t.deathtime IS NULL AND t.dod>t.landmark_time::date AND t.dod<(t.intime+interval '28 days')::date) THEN 1 ELSE 0 END AS death28_lower,
 CASE WHEN(t.deathtime>t.landmark_time AND t.deathtime<=t.intime+interval '28 days')
 OR(t.deathtime IS NULL AND t.dod>t.landmark_time::date AND t.dod<=(t.intime+interval '28 days')::date) THEN 1 ELSE 0 END AS death28_upper,
 t.hospital_expire_flag AS in_hospital_death_flag,t.intime AS source_T0,t.landmark_time AS source_L,
 (t.deathtime IS NOT NULL AND t.dod IS NOT NULL AND t.deathtime::date<>t.dod) AS death_date_conflict,
 g.protein_band,g.specimen_id AS protein_specimen_id,g.field_class AS protein_field_class,g.comment_class AS protein_comment_class,
 g.comment_only_source,extract(epoch FROM(g.charttime-t.intime))/3600 AS protein_sample_hour,
 extract(epoch FROM(g.storetime-t.intime))/3600 AS protein_report_hour,
 (g.specimen_id IS NOT NULL AND g.specimen_id<>a.specimen_id) AS first_specimen_taken_over,
 (NOT r.main_pre AND NOT r.unknown_pre) AS incident_rrt_eligible,
 extract(epoch FROM(r.new_index_tx-t.landmark_time))/3600 AS rrt_candidate_after_L_hours,
 extract(epoch FROM(r.post_tx_before_death_adjudication-t.landmark_time))/3600 AS rrt_unadjudicated_after_L_hours,
 (r.first_supported_tx=t.landmark_time) AS first_supported_tx_exactly_L,
 r.any_active_at_L,
 CASE WHEN r.first_supported_tx IS NOT NULL THEN extract(epoch FROM(r.first_supported_tx-t.intime))/3600 END AS first_supported_tx_hour,
 c.prior_ckd_recorded,
 c.sofa_respiration,
 c.sofa_coagulation,
 c.sofa_liver,
 c.sofa_cardiovascular,
 c.sofa_cns,
 c.nonrenal_sofa,
 c.nonrenal_sofa_missing_components,
 c.sofa_source_row_present,
 c.sofa_source_identifiers_match,
 c.max_cruo_stage_to_L,
 c.stage_observed_rows,
 c.stage_both_components_null_rows,
 c.stage_invalid_rows,
 c.first_scr_0_24_mgdl,
 c.last_scr_0_L_mgdl,
 c.log2_last_scr,
 c.log2_scr_change,
 c.first_scr_charttime,
 c.last_scr_charttime,
 c.first_scr_storetime,
 c.last_scr_storetime,
 c.scr_raw_rows,
 c.scr_storetime_null_rows,
 c.scr_unavailable_by_L_rows,
 c.scr_conflict_rows,
 c.scr_unparsed_value_with_valuenum_rows,
 c.scr_unsupported_unit_rows,
 c.scr_nonpositive_or_nonfinite_rows,
 c.scr_invalid_rows,
 c.uo_rate_6hr_valid,
 c.uo_rate_charttime,
 c.uo_rate_weight,
 c.uo_rate_tm6,
 c.uo_covered_hours_6h,
 c.uo_record_hour_coverage,
 c.uo_valid_record_rows_6h,
 c.uo_explicit_zero_record_rows_6h,
 c.rrt_state_L,
 c.main_pre,
 c.unknown_pre,
 c.presence_pre,
 c.catheter_pre,
 c.scuf_pre,
 c.urine_specific_gravity,
 c.urine_rbc_hpf,
 c.urine_wbc_hpf,
 c.log1p_urine_rbc_hpf,
 c.log1p_urine_wbc_hpf,
 c.specimen_numeric_conflict_items,
 c.specimen_unparsed_value_with_valuenum_items,
 c.specimen_nonhpf_unit_items,
 c.specimen_invalid_items
FROM p_target t JOIN p_r r USING(stay_id)
JOIN covariates c USING(stay_id)
LEFT JOIN p_first_good g USING(stay_id) LEFT JOIN p_first_any a USING(stay_id);
CREATE UNIQUE INDEX ON extraction_cohort(stay_id);
CREATE UNIQUE INDEX ON extraction_cohort(subject_id);
DO $$ BEGIN
 IF (SELECT count(*) FROM extraction_cohort)<>(SELECT count(*) FROM p_target) THEN RAISE EXCEPTION 'QC_COHORT_JOIN_CARDINALITY'; END IF;
 IF EXISTS(SELECT 1 FROM extraction_cohort WHERE max_cruo_stage_to_L NOT BETWEEN 1 AND 3) THEN RAISE EXCEPTION 'QC_STAGE'; END IF;
 IF EXISTS(SELECT 1 FROM extraction_cohort WHERE protein_band IS NOT NULL AND(protein_sample_hour<0 OR protein_sample_hour>=48 OR protein_report_hour>72)) THEN RAISE EXCEPTION 'QC_PROTEIN_WINDOW'; END IF;
 IF EXISTS(SELECT 1 FROM extraction_cohort WHERE uo_covered_hours_6h NOT BETWEEN 0 AND 6) THEN RAISE EXCEPTION 'QC_UO_HOURS'; END IF;
 IF EXISTS(SELECT 1 FROM extraction_cohort WHERE sofa_source_identifiers_match=false) THEN RAISE EXCEPTION 'QC_SOFA_IDENTITY'; END IF;
 IF EXISTS(SELECT 1 FROM extraction_cohort WHERE last_scr_charttime>=source_L OR last_scr_charttime<source_T0 OR last_scr_storetime>source_L OR first_scr_charttime>=source_T0+interval '24h') THEN RAISE EXCEPTION 'QC_SCR_WINDOW'; END IF;
END $$;
INSERT INTO renal_result
SELECT 'sample_flow',72,'main',v.metric,v.n,NULL,NULL,'Counts remain local exact; shared output is banded'
FROM(VALUES
 ('all_ICU',(SELECT count(*) FROM renal_icu)),('early_sepsis_adult',(SELECT count(*) FROM renal_s)),
 ('CrUO_AKI48',(SELECT count(*) FROM renal_c)),('A_ge_S',(SELECT count(*) FROM p_clinical)),
 ('no_preS',(SELECT count(*) FROM p_eligible)),('no_maintenance',(SELECT count(*) FROM p_eligible WHERE NOT maintenance)),
 ('first_clinically_qualified',(SELECT count(*) FROM p_first)),('alive_in_ICU_L',(SELECT count(*) FROM p_target)),
 ('any_protein_record',(SELECT count(DISTINCT stay_id) FROM p_lab)),('measured_six_band',(SELECT count(*) FROM p_first_good))
)v(metric,n);
INSERT INTO renal_result
SELECT 'boundaries',72,'main',v.metric,count(*) FILTER(WHERE v.met),count(*),NULL,'Observed boundary counts; nonzero not automatically a failure'
FROM extraction_cohort CROSS JOIN LATERAL(VALUES('first_supported_RRT_exactly_L',first_supported_tx_exactly_L),
 ('any_active_at_L',any_active_at_L),('death_date_conflict',death_date_conflict),
 ('enhanced_all_numeric',protein_band IS NOT NULL AND urine_specific_gravity IS NOT NULL AND urine_rbc_hpf IS NOT NULL AND urine_wbc_hpf IS NOT NULL))v(metric,met) GROUP BY v.metric;
INSERT INTO renal_result
SELECT 'protein_loss',48,coalesce(loss_reason,'interpretable'),'canonical_specimens',count(*),NULL,NULL,'Hierarchical specimen reasons, not patient losses'
FROM p_canonical GROUP BY loss_reason;
INSERT INTO renal_result
SELECT 'specimen_units',72,itemid::text||'|'||CASE WHEN unit IN('/hpf','#/hpf','cells/hpf','rbc/hpf','wbc/hpf','hpf','') THEN coalesce(nullif(unit,''),'missing') ELSE 'other' END,
 'selected_specimen_items',count(*),NULL,NULL,'Controlled unit labels only' FROM cov_specimen_latest GROUP BY 3;
COMMIT;
SELECT section,window_h,stratum,metric,pg_temp.renal_band(n) AS n_text,pg_temp.renal_band(denominator) AS denominator_text,value_text,note FROM renal_result ORDER BY section,stratum,metric;
