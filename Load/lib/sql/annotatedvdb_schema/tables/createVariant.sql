-- table for variants
DROP TABLE IF EXISTS AnnotatedVDB.Variant CASCADE;

CREATE UNLOGGED TABLE AnnotatedVDB.Variant (
       CHROMOSOME           CHARACTER VARYING(10) NOT NULL,
       RECORD_PRIMARY_KEY   CHARACTER VARYING(100) NOT NULL, 
       POSITION 	    INTEGER NOT NULL,
       IS_MULTI_ALLELIC	    BOOLEAN,
       IS_ADSP_VARIANT	    BOOLEAN,
       REF_SNP_ID	    CHARACTER VARYING(25),
       METASEQ_ID	    TEXT NOT NULL,
       BIN_INDEX	    LTREE NOT NULL,
       DISPLAY_ATTRIBUTES   JSONB,
       ALLELE_FREQUENCIES      JSONB,
       CADD_SCORES	       JSONB,
       ADSP_MOST_SEVERE_CONSEQUENCE JSONB,
       ADSP_RANKED_CONSEQUENCES JSONB,
       LOSS_OF_FUNCTION		JSONB,
       VEP_OUTPUT	       JSONB,
       ADSP_QC		       JSONB,
       GWAS_FLAGS	       JSONB,
       OTHER_ANNOTATION	       JSONB,
       ROW_ALGORITHM_ID	       INTEGER NOT NULL
) PARTITION BY LIST (CHROMOSOME);

-- CREATE PARTITIONS

CREATE OR REPLACE FUNCTION "public"."create_variant_partitions" ()  RETURNS integer
  VOLATILE
  AS $body$
DECLARE
      partition TEXT;
      chr TEXT;
    BEGIN

      FOR chr IN
        SELECT UNNEST(string_to_array('chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM', ' '))
      LOOP   
        partition := 'AnnotatedVDB.VARIANT' || '_' || chr::text;
        IF NOT EXISTS(SELECT relname FROM pg_class WHERE relname=partition) THEN
           EXECUTE 'CREATE TABLE ' || partition || ' PARTITION OF AnnotatedVDB.Variant FOR VALUES IN (''' || chr || ''')';
           RAISE NOTICE 'A partition has been created %',partition;
         END IF;
      END LOOP; 
      RETURN NULL;
    END;
$body$ LANGUAGE plpgsql;

SELECT create_variant_partitions();

-- TRIGGERS

CREATE OR REPLACE FUNCTION AnnotatedVDB.set_bin_index() RETURNS TRIGGER
    LANGUAGE plpgsql
    AS $$
BEGIN
  IF NEW.bin_index IS NULL THEN
    NEW.bin_index = find_bin_index(NEW.chromosome, NEW.location, NEW.location);
  END IF;
  RETURN NEW;
END;
$$;

CREATE TRIGGER variant_set_bin_index
       AFTER INSERT ON AnnotatedVDB.Variant
       EXECUTE PROCEDURE AnnotatedVDB.set_bin_index();

------- #######

/* CREATE OR REPLACE FUNCTION AnnotatedVDB.set_record_primary_key() RETURNS TRIGGER
    LANGUAGE plpgsql
    AS $$
BEGIN
  IF NEW.record_primary_key IS NULL THEN
    NEW.record_primary_key = truncate_str(NEW.metaseq_id, 350) || COALESCE('_' || NEW.ref_snp_id, '');
  END IF;
  RETURN NEW;
END;
$$;


CREATE TRIGGER variant_set_record_primary_key
       AFTER INSERT ON AnnotatedVDB.Variant
       EXECUTE PROCEDURE AnnotatedVDB.set_record_primary_key(); */

-- INDEXES

CREATE INDEX VARIANT_RECORD_PK_HASH ON AnnotatedVDB.Variant USING HASH(RECORD_PRIMARY_KEY);
CREATE INDEX VARIANT_RS_HASH ON AnnotatedVDB.Variant USING HASH(REF_SNP_ID);
CREATE INDEX VARIANT_METASEQ_LEFT50 ON AnnotatedVDB.Variant(LEFT(METASEQ_ID,50));
CREATE INDEX VARIANT_BIN_INDEX ON AnnotatedVDB.Variant USING GIST(BIN_INDEX);
CREATE INDEX VARIANT_UNDO_INDEX ON AnnotatedVDB.Variant USING BRIN(CHROMOSOME, ROW_ALGORITHM_ID);


