
CREATE OR REPLACE FUNCTION AnnotatedVDB.generate_alt_metaseq_id(metaseqId TEXT)
       RETURNS TEXT AS $$
DECLARE
	altId TEXT;
BEGIN
	WITH _array AS (SELECT regexp_split_to_array(metaseqId, ':') AS VALUES)
	SELECT CONCAT(_array.values[1], ':'::text, _array.values[2], ':'::text, _array.values[4], ':'::text, _array.values[3]) INTO altId FROM _array;
	RETURN altId;
END;

$$ LANGUAGE plpgsql;

CREATE OR REPLACE FUNCTION AnnotatedVDB.find_variant_by_metaseq_id_variations(metaseqId TEXT, firstHitOnly BOOLEAN DEFAULT TRUE)
       RETURNS TABLE(record_primary_key CHARACTER VARYING(100)) AS $$

BEGIN
	RETURN QUERY
    	SELECT v.record_primary_key FROM AnnotatedVDB.find_variant_by_metaseq_id(metaseqId, firstHitOnly) v
        UNION ALL
        SELECT v.record_primary_key FROM AnnotatedVDB.find_variant_by_metaseq_id(AnnotatedVDB.generate_alt_metaseq_id(metaseqId), firstHitOnly) v
	LIMIT CASE WHEN firstHitOnly THEN 1 END;
END;

$$ LANGUAGE plpgsql;

CREATE OR REPLACE FUNCTION AnnotatedVDB.find_variant_by_metaseq_id(metaseqId TEXT, firstHitOnly BOOLEAN DEFAULT TRUE)
       RETURNS TABLE(record_primary_key CHARACTER VARYING(100)) AS $$
BEGIN
	RETURN QUERY
	SELECT v.record_primary_key
	FROM AnnotatedVDB.Variant v
	WHERE v.metaseq_id = metaseqId
	AND LEFT(v.metaseq_id, 50) = LEFT(metaseqId, 50)
	AND chromosome = 'chr' || split_part(metaseqId, ':', 1)::text
	LIMIT CASE WHEN firstHitOnly THEN 1 END;
END;

$$ LANGUAGE plpgsql;



