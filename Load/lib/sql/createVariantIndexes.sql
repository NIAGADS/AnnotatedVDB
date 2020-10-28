-- CREATE INDEX VARIANT_RS_HASH ON Variant USING HASH(REF_SNP_ID);
-- CREATE INDEX VARIANT_UNDO ON Variant(CHROMOSOME, ROW_ALGORITHM_ID);
-- CREATE INDEX VARIANT_METASEQ_LEFT50 ON Variant(LEFT(METASEQ_ID,50));
-- CREATE INDEX VARIANT_BIN_INDEX ON Variant USING GIST(BIN_INDEX);

-- CREATE INDEX VARIANT_LOCATION ON Variant(CHROMOSOME, LOCATION ASC);
-- CREATE INDEX VARIANT_NESTED_LOCATION ON VARIANT(CHROMOSOME, BIN_INDEX, LOCATION ASC);
-- CREATE INDEX VARIANT_IS_ADSP ON Variant(IS_ADSP_VARIANT) WHERE IS_ADSP_VARIANT IS TRUE;
-- CREATE INDEX VARIANT_IS_ADSP_WES ON Variant(is_adsp_variant) WHERE is_adsp_variant IS TRUE AND (other_annotation->>'GenomicsDB')::jsonb @> '["ADSP_WES"]';



CREATE OR REPLACE FUNCTION "public"."create_variant_pk_indexes" ()  RETURNS integer
  VOLATILE
  AS $body$
DECLARE
      partition TEXT;
      index_name TEXT;
      chr TEXT;
    BEGIN

      FOR chr IN
        SELECT UNNEST(string_to_array('chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM', ' '))
      LOOP   
        partition := 'VARIANT' || '_' || chr::text;
        index_name := 'GENOMICSDB_PK' || '_' || chr::text;
        RAISE NOTICE 'Creating Index %', index_name;
        EXECUTE 'CREATE INDEX ' || index_name || ' ON ' || partition || '(record_primary_key(' || partition || '))';
        RAISE NOTICE 'A index has been created %',partition;
      END LOOP; 
      RETURN NULL;
    END;
$body$ LANGUAGE plpgsql;

SELECT create_variant_pk_indexes();

DROP FUNCTION public.create_variant_pk_indexes();
