# Count ExAC variants grouped by both chromosome and number of alternate alleles.
SELECT
  reference_name,
  num_alternates,
  COUNT(num_alternates) AS num_multi_allelic_variants
FROM (
  SELECT
    reference_name,
    start,
    reference_bases,
    COUNT(alternate_bases) AS num_alternates
  FROM
    [google.com:biggene:ExAC_release_0_2.reshaped_variants]
  GROUP EACH BY
    reference_name,
    start,
    reference_bases
    )
GROUP BY
  reference_name,
  num_alternates
ORDER BY
  reference_name,
  num_alternates
