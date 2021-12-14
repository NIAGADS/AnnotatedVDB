-- table for tracking variant loads to facilitate undos
DROP TABLE IF EXISTS AnnotatedVDB.AlgorithmInvocation;

CREATE TABLE AnnotatedVDB.AlgorithmInvocation (
       ALGORITHM_INVOCATION_ID	  SERIAL PRIMARY KEY,
       SCRIPT_NAME	          CHARACTER VARYING(50) NOT NULL,
       SCRIPT_PARAMETERS          TEXT,
       COMMIT_MODE		  BOOLEAN,
       RUN_TIME 	          TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP
);

-- INDEXES

CREATE INDEX ALG_INV_INDX01 ON AnnotatedVDB.AlgorithmInvocation(SCRIPT_NAME);
CREATE INDEX ALG_INV_INDX02 ON AnnotatedVDB.AlgorithmInvocation(RUN_TIME);
