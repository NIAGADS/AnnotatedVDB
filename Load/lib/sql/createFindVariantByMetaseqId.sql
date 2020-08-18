-- finds variant by chr:pos:ref:alt id

CREATE OR REPLACE FUNCTION find_variant_by_metaseq_id(metaseqId TEXT, firstHitOnly BOOLEAN DEFAULT FALSE)
       RETURNS TABLE(record_primary_key TEXT, ref_snp_id CHARACTER VARYING, metaseq_id TEXT,
       	             has_genomicsdb_annotation BOOLEAN, is_adsp_variant BOOLEAN, bin_index LTREE) AS $$
BEGIN
	RETURN QUERY
	SELECT v.record_primary_key, v.ref_snp_id, v.metaseq_id, v.has_genomicsdb_annotation, v.is_adsp_variant,
	v.bin_index
	FROM Variant v
	WHERE v.metaseq_id = metaseqId
	AND chromosome = 'chr' || split_part(metaseqId, ':', 1)::text
	LIMIT CASE WHEN firstHitOnly THEN 1 END;
END;

$$ LANGUAGE plpgsql;


CREATE OR REPLACE FUNCTION find_variant_by_metaseq_id(metaseqId TEXT, chrm TEXT, firstHitOnly BOOLEAN DEFAULT FALSE)
       RETURNS TABLE(record_primary_key TEXT, ref_snp_id CHARACTER VARYING, metaseq_id TEXT,
       	             has_genomicsdb_annotation BOOLEAN, is_adsp_variant BOOLEAN, bin_index LTREE) AS $$
BEGIN
	RETURN QUERY
	SELECT v.record_primary_key, v.ref_snp_id, v.metaseq_id, v.has_genomicsdb_annotation, v.is_adsp_variant,
	v.bin_index
	FROM Variant v
	WHERE v.metaseq_id = metaseqId
	AND chromosome = chrm
	LIMIT CASE WHEN firstHitOnly THEN 1 END;
END;

$$ LANGUAGE plpgsql;

--DROP FUNCTION get_variant_annotation_by_metaseq_id(text,text,boolean);

CREATE OR REPLACE FUNCTION get_variant_annotation_by_metaseq_id(metaseqId TEXT, chrm TEXT, firstHitOnly BOOLEAN DEFAULT FALSE)
       RETURNS TABLE(record_primary_key TEXT, ref_snp_id CHARACTER VARYING, metaseq_id TEXT,
       	             has_genomicsdb_annotation BOOLEAN, is_adsp_variant BOOLEAN, bin_index LTREE,
adsp_most_severe_consequence JSONB, cadd_scores JSONB, allele_frequencies JSONB) AS $$

BEGIN
	RETURN QUERY
	SELECT v.record_primary_key, v.ref_snp_id, v.metaseq_id, v.has_genomicsdb_annotation, v.is_adsp_variant,
	v.bin_index, v.adsp_most_severe_consequence, v.cadd_scores, v.allele_frequencies 
	FROM Variant v
	WHERE v.metaseq_id = metaseqId
	AND chromosome = chrm
	LIMIT CASE WHEN firstHitOnly THEN 1 END;
END;

$$ LANGUAGE plpgsql;


CREATE OR REPLACE FUNCTION get_variant_annotation_by_metaseq_id(metaseqId TEXT, firstHitOnly BOOLEAN DEFAULT FALSE)
       RETURNS TABLE(record_primary_key TEXT, ref_snp_id CHARACTER VARYING, metaseq_id TEXT,
       	             has_genomicsdb_annotation BOOLEAN, is_adsp_variant BOOLEAN, bin_index LTREE,
adsp_most_severe_consequence JSONB, cadd_scores JSONB, allele_frequencies JSONB) AS $$

BEGIN
	RETURN QUERY
	SELECT v.record_primary_key, v.ref_snp_id, v.metaseq_id, v.has_genomicsdb_annotation, v.is_adsp_variant,
	v.bin_index, v.adsp_most_severe_consequence, v.cadd_scores, v.allele_frequencies 
	FROM Variant v
	WHERE v.metaseq_id = metaseqId
	AND chromosome = 'chr' || split_part(metaseqId, ':', 1)::text
	LIMIT CASE WHEN firstHitOnly THEN 1 END;
END;

$$ LANGUAGE plpgsql;
