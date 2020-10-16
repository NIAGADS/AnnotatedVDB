CREATE OR replace FUNCTION json_array_cast_to_text(json) RETURNS text[] AS $f$
    SELECT array_agg(x#>>'{}') FROM json_array_elements($1) t(x);
$f$ LANGUAGE SQL IMMUTABLE;
