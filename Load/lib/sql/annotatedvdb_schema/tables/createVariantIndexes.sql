

-- CREATE INDEX VARIANT_LOCATION ON Variant(CHROMOSOME, LOCATION ASC);
-- CREATE INDEX VARIANT_NESTED_LOCATION ON VARIANT(CHROMOSOME, BIN_INDEX, LOCATION ASC);

-- CREATE INDEX VARIANT_IS_ADSP ON Variant(IS_ADSP_VARIANT) WHERE IS_ADSP_VARIANT IS TRUE;
-- CREATE INDEX VARIANT_IS_ADSP_WES ON Variant(is_adsp_variant) WHERE is_adsp_variant IS TRUE AND (other_annotation->>'GenomicsDB')::jsonb @> '["ADSP_WES"]';
