-- table for variants
DROP TABLE IF EXISTS AlgorithmInvocation;

CREATE TABLE AlgorithmInvocation (
       ALGORITHM_INVOCATION_ID	  SERIAL PRIMARY KEY,
       SCRIPT_NAME	          CHARACTER VARYING(50) NOT NULL,
       SCRIPT_PARAMETERS          TEXT,
       COMMIT_MODE		  BOOLEAN,
       RUN_TIME 	          TIMESTAMP WITH TIME ZONE DEFAULT CURRENT_TIMESTAMP
);

-- INDEXES

CREATE INDEX ALG_INV_INDX01 ON AlgorithmInvocation(SCRIPT_NAME);
CREATE INDEX ALG_INV_INDX02 ON AlgorithmInvocation(RUN_TIME);
