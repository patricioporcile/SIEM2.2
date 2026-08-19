# Enrichment tables pending regeneration

The raw GO Biological Process / KEGG / ChEA enrichment files for the current
20-gene concordant set (20-GCS) are NOT included here. The versions previously
saved under these names (GO_enrichment_20GCS.csv, KEGG_enrichment_20GCS.csv,
ChEA2022_TF_enrichment_20GCS.csv) were found to be stale duplicates of an
earlier, non-replicated 63-gene analysis (GSE8059) and were removed rather
than uploaded under a misleading name.

What IS correct and included:
- enrichment_on_null_sets.csv — the empirical enrichment-on-permuted-nulls
  test (manuscript Table 3), verified consistent with the current manuscript.
- The significant terms (FDR < 0.05) from a proper Enrichr run on the 20-GCS
  are already reported in the manuscript's Table S4 (Supplementary Materials),
  and can be re-derived by re-running Enrichr (GO_Biological_Process_2021,
  KEGG_2021_Human, ChEA_2022) on the 20-gene list in
  convergent_genes_RRA_GSE50838_20genes.csv.

To regenerate: run notebook 04_enrichment_analysis.ipynb against the 20-gene
list, or query https://maayanlab.cloud/Enrichr/ manually and save the full
result tables here under the original file names.
