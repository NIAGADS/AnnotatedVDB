DROP TABLE IF EXISTS AnnotatedVDB.:CHRM;
CREATE TABLE AnnotatedVDB.:CHRM AS (SELECT * FROM (
    WITH Duplicates AS (
        SELECT record_primary_key FROM (
            SELECT record_primary_key, count(record_primary_key) AS row_count
            FROM AnnotatedVDB.Variant 
            WHERE chromosome = (:'CHRM')::text
            GROUP BY record_primary_key
        ) c WHERE c.row_count > 1
    )
    SELECT DISTINCT v.*
    FROM AnnotatedVDB.Variant v, Duplicates d
    WHERE v.chromosome = (:'CHRM')::text
    AND v.record_primary_key = d.record_primary_key
    ) a
);

SELECT count(*) FROM AnnotatedVDB.:CHRM;

DELETE FROM AnnotatedVDB.Variant
WHERE chromosome = (:'CHRM')::text
AND record_primary_key IN (SELECT record_primary_key FROM AnnotatedVDB.:CHRM);

INSERT INTO AnnotatedVDB.Variant SELECT * FROM AnnotatedVDB.:CHRM;

-- /////////////////////////////////////////////////////////
-- irregular variants
-- /////////////////////////////////////////////////////////
DROP TABLE IF EXISTS AnnotatedVDB.:CHRM;
CREATE TABLE AnnotatedVDB.:CHRM AS (
    SELECT * FROM AnnotatedVDB.Variant
    WHERE chromosome = (:'CHRM')::text
    AND metaseq_id SIMILAR TO '%I%|%D%|%R%|%N%|%[?]%'
); 

SELECT count(*) FROM AnnotatedVDB.:CHRM;
SELECT count(*) FROM AnnotatedVDB.:CHRM WHERE ref_snp_id IS NULL;

DELETE FROM AnnotatedVDB.Variant
WHERE chromosome = (:'CHRM')::text
AND record_primary_key IN (SELECT record_primary_key FROM AnnotatedVDB.:CHRM WHERE ref_snp_id IS NULL);

DROP TABLE IF EXISTS AnnotatedVDB.:CHRM;

--psql -h $DB_HOST -U allenem -d $DB_NAME --file /home/allenem/GRCh38/project_home/AnnotatedVDB/Load/lib/sql/annotatedvdb_schema/removeDuplicates.sql  -a -v CHRM=chr5