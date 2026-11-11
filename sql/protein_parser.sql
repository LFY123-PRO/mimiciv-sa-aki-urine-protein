
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
