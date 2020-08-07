-- Fast DISTINCT using the index
CREATE OR REPLACE FUNCTION fast_distinct(
   tableName varchar, fieldName varchar, sample anyelement = ''::varchar)
   -- Search a few distinct values in a possibly huge table
   -- Parameters: tableName or query expression, fieldName,
   --             sample: any value to specify result type (defaut is varchar)
   -- Author: T.Husson, 2012-09-17, distribute/use freely
   RETURNS TABLE ( result anyelement ) AS
$BODY$
BEGIN
   EXECUTE 'SELECT '||fieldName||' FROM '||tableName||' ORDER BY '||fieldName
      ||' LIMIT 1'  INTO result;
   WHILE result IS NOT NULL LOOP
      RETURN NEXT;
      EXECUTE 'SELECT '||fieldName||' FROM '||tableName
         ||' WHERE '||fieldName||' > $1 ORDER BY ' || fieldName || ' LIMIT 1'
         INTO result USING result;
   END LOOP;
END;
$BODY$ LANGUAGE plpgsql VOLATILE;

-- Disk usage
CREATE OR REPLACE VIEW diskusage AS SELECT nspname || '.' || relname AS "relation",
    pg_size_pretty(pg_relation_size(C.oid)) AS "size"
  FROM pg_class C
  LEFT JOIN pg_namespace N ON (N.oid = C.relnamespace)
  WHERE nspname NOT IN ('pg_catalog', 'information_schema')
  ORDER BY pg_relation_size(C.oid) DESC;
