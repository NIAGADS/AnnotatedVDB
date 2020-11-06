-- finds most specific inclusive bin

CREATE OR REPLACE FUNCTION get_chr_size(chr VARCHAR)
       RETURNS BIGINT AS $$
DECLARE 
	chrSize BIGINT;
BEGIN
	SELECT UPPER(location) INTO chrSize
	FROM BinIndexRef
	WHERE chromosome = chr
	AND level = 0;
	
	RETURN chrSize;
END;
$$
LANGUAGE plpgsql;


CREATE OR REPLACE FUNCTION find_bin_index(chr VARCHAR, loc_start BIGINT, loc_end BIGINT) 
        RETURNS LTREE AS $$
DECLARE 
        bin LTREE;
	chrSize BIGINT;
BEGIN
	SELECT get_chr_size(chr) INTO chrSize;

        SELECT global_bin_path INTO bin 
        FROM BinIndexRef
        WHERE chromosome = chr 
	
	-- case statements handle programmatic flanking regions that may extend past chr boundaries	
        AND location @> int8range(CASE WHEN loc_start < 1 THEN 1 
	    	     		       WHEN loc_start >= chrSize THEN chrSize - 1 
				       ELSE loc_start END, 
	CASE WHEN loc_end >= chrSize THEN chrSize - 1 ELSE loc_end END, '[]')
	
	ORDER BY nlevel(global_bin_path) DESC
	LIMIT 1;	
        
        RETURN bin;
END;
$$
LANGUAGE plpgsql


