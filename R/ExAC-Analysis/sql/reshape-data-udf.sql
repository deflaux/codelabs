# Reshape this data so that we have one record (row) per alternate allele.
#
# Note that the new BigQuery feature of user-defined javascript
# functions is in limited preview.  For more info, see
# https://www.youtube.com/watch?v=GrD7ymUPt3M#t=1377
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
  for(var alt = 1; alt <= r.alternate_bases.length; alt++) {
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
