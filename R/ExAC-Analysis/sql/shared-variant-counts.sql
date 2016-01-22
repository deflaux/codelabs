# Count the number of variants found in 0.00% of alleles (due to
# rounding), 0.01% of alleles, 0.02% of alleles, etc...
SELECT
  percentage_of_alleles_with_variant,
  COUNT(1) AS number_of_variants_shared_by_this_percentage_of_alleles
FROM (
  SELECT
    reference_name,
    start,
    end,
    reference_bases,
    alternate_bases,
    ROUND(AC_Adj/AN_Adj, 3) AS percentage_of_alleles_with_variant
    # AF could also have been used in place of AC_Adj/AN_Adj
  FROM
    [google.com:biggene:ExAC_release_0_2.reshaped_variants]
  WHERE reference_name NOT IN ("X", "Y")
  OMIT RECORD IF AN_Adj < 10000)
GROUP BY
  percentage_of_alleles_with_variant
ORDER BY
  percentage_of_alleles_with_variant
