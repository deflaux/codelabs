# Reshape this data so that we have one record (row) per alternate allele.
#
# This query uses BigQuery user-defined javascript functions
# https://cloud.google.com/bigquery/user-defined-functions
#
# This approach will work for all the fields.  The example below includes
# a handful of fields just for demonstration purposes.
SELECT
   reference_name, start, end, reference_bases, alternate_bases, AC_NFE
FROM js(
  [google.com:biggene:ExAC_release_0_2.variants],
  reference_name, start, end, reference_bases, alternate_bases, AC_NFE,
  '[
  {"name": "reference_name", "type": "STRING"},
  {"name": "start", "type": "INTEGER"},
  {"name": "end", "type": "INTEGER"},
  {"name": "reference_bases", "type": "STRING"},
  {"name": "alternate_bases", "type": "STRING"},
  {"name": "AC_NFE", "type": "INTEGER"}
  ]',
  "function(r, emit) {
  for(var alt = 0; alt <= r.alternate_bases.length; alt++) {
    emit({
          reference_name: r.reference_name,
          start: r.start,
          reference_bases: r.reference_bases,
          end: r.end,
          alternate_bases: r.alternate_bases[alt],
          AC_NFE: r.AC_NFE[alt]
       });
        }
    }")
