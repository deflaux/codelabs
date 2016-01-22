# Count the number of variants found in 0.00% of alleles (due to
# rounding), 0.01% of alleles, 0.02% of alleles, etc...
SELECT
  percentage_of_alleles_with_variant,
  COUNT(1) AS number_of_variants_shared_by_this_percentage_of_alleles,
  population,
  description
FROM (
  SELECT
    reference_name,
    start,
    end,
    reference_bases,
    alternate_bases,
    population,
    description,
    CASE
      WHEN population = 'AFR' THEN ROUND(AC_AFR/AN_AFR, 3)
      WHEN population = 'AMR' THEN ROUND(AC_AMR/AN_AMR, 3)
      WHEN population = 'EAS' THEN ROUND(AC_EAS/AN_EAS, 3)
      WHEN population = 'FIN' THEN ROUND(AC_FIN/AN_FIN, 3)
      WHEN population = 'NFE' THEN ROUND(AC_NFE/AN_NFE, 3)
      WHEN population = 'OTH' THEN ROUND(AC_OTH/AN_OTH, 3)
      WHEN population = 'SAS' THEN ROUND(AC_SAS/AN_SAS, 3)
      END AS percentage_of_alleles_with_variant
  FROM
    [google.com:biggene:ExAC_release_0_2.reshaped_variants] AS variants
  CROSS JOIN
    [google.com:biggene:ExAC_release_0_2.populations] AS populations
  WHERE reference_name NOT IN ("X", "Y")
  OMIT RECORD IF AN_Adj < 10000
  HAVING percentage_of_alleles_with_variant IS NOT NULL)
GROUP BY
  percentage_of_alleles_with_variant,
  population,
  description
ORDER BY
  percentage_of_alleles_with_variant,
  population,
  description
