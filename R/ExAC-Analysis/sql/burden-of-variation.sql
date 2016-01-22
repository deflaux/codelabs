# Calculate the null spectrum (burden of variation) for each gene: How many 
# individuals have 2 or more missense mutations at < 2% frequency for each gene? 
SELECT
  Gene,
  ExonicFunc,
  call.call_set_name,
  COUNT(call.call_set_name) AS number_of_rare_missense_variants
FROM
  [silver-wall-555:TuteTable.hg19] AS annots
JOIN
  EACH
  (
  SELECT
    reference_name,
    start,
    reference_bases,
    # 1,000 Genomes phase 1 variants are bi-allelic
    NTH(1, alternate_bases) WITHIN RECORD AS alternate_bases,
    call.call_set_name,
  FROM
    [genomics-public-data:1000_genomes.variants]
  WHERE
    AF <= 0.02
  OMIT
    call IF EVERY(call.genotype != 1)
    ) AS vars
ON
  vars.reference_name = annots.Chr
  AND vars.start = annots.Start
  AND vars.reference_bases = annots.Ref
  AND vars.alternate_bases = annots.Alt
WHERE
  ExonicFunc = "missense"
GROUP EACH BY
  Gene,
  ExonicFunc,
  call.call_set_name,
HAVING
  number_of_rare_missense_variants > 1
ORDER BY
  Gene,
  ExonicFunc,
  call.call_set_name
