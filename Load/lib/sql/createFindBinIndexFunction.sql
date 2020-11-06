-- finds most specific inclusive bin

CREATE OR REPLACE FUNCTION get_chr_size(chr VARCHAR)
       RETURNS BIGINT AS $$
DELCARE 
	chrSize BIGINT;
BEGIN
	SELECT UPPER(location) INTO chrSize
	FROM BinINdexRef
	WHERE chromosome = chr
	AND level = 0;
	
	RETURN chrSize;
END;
$$
LANGUAGE plpgsql

CREATE OR REPLACE FUNCTION find_bin_index(chr VARCHAR, loc_start BIGINT, loc_end BIGINT) 
        RETURNS LTREE AS $$
DECLARE 
        bin LTREE;
BEGIN
        SELECT global_bin_path INTO bin 
        FROM BinIndexRef 
        WHERE chromosome = chr 
        AND location @> int8range(loc_start, loc_end, '[]')
	ORDER BY nlevel(global_bin_path) DESC
	LIMIT 1;	
        
        RETURN bin;
END;
$$
LANGUAGE plpgsql


