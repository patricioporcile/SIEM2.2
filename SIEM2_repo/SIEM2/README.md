# Null Model Specification in Cross-Lineage Transcriptomic Convergence

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

> **Null Model Specification Determines Whether Cross-Lineage Transcriptional Convergence Is Detected: a Worked Example in Neutrophils and NK Cells**
> Patricio E. Porcile, Universidad de Atacama, Chile

## Overview

This repository contains the full computational pipeline supporting the manuscript above. The manuscript's central contribution is **methodological**: multi-dataset transcriptomic meta-analyses commonly test cross-dataset gene overlap against a permutation null that shuffles each dataset's fold-change vector independently, implicitly assuming the datasets are independent under the null. When two or more inputs are legitimately correlated for reasons unrelated to the effect under study, that assumption fails and an unremarkable overlap can appear highly significant.

Using the question of whether neutrophils and NK cells show transcriptional convergence as a **worked example**, this repository reports both the naive (independent-shuffle) and the corrected (correlation-preserving) null side by side, on the same data, so the size of the specification error is directly visible: **p < 0.0001 vs. p = 0.60** for the identical 20-gene overlap (Section 3.2.1). The correction is then extended to the resulting gene set's functional enrichment (Section 3.9) and validated with an independent positive control demonstrating the pipeline detects convergence when it genuinely exists (Section 3.10).

The 20-gene fold-change-concordant set identified along the way (referred to in the manuscript as the "20-gene concordant set" or 20-GCS, **not** as a statistically validated module) is reported as a hypothesis-generating byproduct of this worked example, not as this study's main claim.

## Key Findings

- **Independent-shuffle null**: p < 0.0001 (null mean 4.8 ± 2.2) for the observed 20-gene three-way overlap — the result reported in earlier versions of this analysis, and the design most multi-dataset meta-analyses default to.
- **NK-only, correlation-preserving null** (holds the real TNF-alpha/LPS pairing fixed, permutes only NK): p = 0.60 (null mean 20.7 ± 4.3) for the identical overlap — not distinguishable from chance.
- **Enrichment-on-null-sets** (Table 7): of 21 enrichment terms tested against the same corrected null, 7 survive Benjamini-Hochberg FDR < 0.05 (a coagulation/fibrinolysis cluster, cell adhesion molecules, T-cell regulatory terms, SMAD1); RELA and the NF-kappaB KEGG pathway — the basis of the "AP-1/NF-kappaB axis" framing in earlier drafts — are the two **least**-supported associations in the table (p = 0.78 and p = 0.53).
- **Positive control** (GSE288312, IL-1beta-treated neutrophils in place of NK cells): 12 of the 224 TNF-alpha/LPS-concordant genes are also concordant with IL-1beta, against a null maximum of 4 across 10,000 permutations (p < 0.0001) — confirming the pipeline is not simply too conservative to detect anything. A sensitivity check restricting the null pool to the 8,581-gene shared background (instead of all 26,485 GSE288312-profiled genes) gives the same observed overlap (12/224) against a null mean of 0.88 ± 0.92 (max 7; p < 0.0001), showing the result is robust to this choice of background.
- The 20-gene concordant set includes two immune checkpoint genes, **CD274 (PD-L1)** and **CD40**, absent from a prior non-replicated analysis using GSE8059.

All biological interpretations of the 20-gene set are explicitly framed as **hypothesis-generating**, consistent with the permutation result in Section 3.2.1: this is a purely computational study with no experimental (wet-lab) validation.

## Datasets

| Dataset | Cell Type | Stimulus | Platform | Design | DE Method | GEO |
|---|---|---|---|---|---|---|
| **GSE50838** | NK (CD56+) | IL-2, 100 U/mL, 5d | Affymetrix HuGene 1.1 ST | 4 donors, paired (naive vs IL-2) | Paired t-test | [link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50838) |
| GSE40548 | Neutrophils | TNF-alpha, 10 ng/mL, 1h | Illumina HiSeq + SOLiD | 3 samples (1 SOLiD donor + 2 Illumina donors — mixed platform, see Limitations) | Cuffdiff | [link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40548) |
| GSE221866 | Neutrophils | LPS, 1 µg/mL, 6h | Illumina NovaSeq 6000 | 3 biological replicates | DESeq2 | [link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE221866) |
| **GSE288312** | Neutrophils | IL-1beta, 1h (positive control) | RNA-seq, raw counts | 3 treated / 6 untreated | Simplified CPM log2FC | [link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE288312) |
| GSE8059* | NK (CD56+/CD16+) | IL-2, 100 IU/mL, 0–24h | Affymetrix U133 Plus 2.0 | 4 pooled donors, no replicates | Fold-change only | [link](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE8059) |

*GSE8059 is used only as a temporal reference and for the replicated-vs-non-replicated comparison; it is not part of the primary convergence analysis. GSE288312 is used only as an independent positive control (Section 3.10); it is not part of the primary three-way convergence test.

GSE50838 GSM accessions used: naive = GSM1231002–GSM1231005, IL-2-activated = GSM1231006–GSM1231009 (8 of the 16 samples in the full series; the K562+IL-2 and R69-coculture arms were not used). GSE221866 was originally reported by Ismailova et al. (2023), *Front Immunol* 14:1123344 (doi: 10.3389/fimmu.2023.1123344). GSE288312 was reported by Thomas, Edwards & Wright (2025), *Discovery Immunology* 4(1):kyaf013 (doi: 10.1093/discim/kyaf013).

## Repository Structure

```
.
├── notebooks/
│   ├── 01_GSE50838_processing.ipynb              # Series matrix -> GPL11532/RefSeq annotation -> RMA -> paired t-test
│   ├── 02_neutrophil_parsing.ipynb               # Parse Cuffdiff (GSE40548) and DESeq2 (GSE221866) outputs
│   ├── 03_meta_analysis_RRA.ipynb                # Rank-aggregation meta-analysis + convergence filtering
│   ├── 04_enrichment_analysis.ipynb              # GO / KEGG / ChEA 2022 (Enrichr) enrichment on the 20-GCS
│   ├── 05_GSE8059_comparison.ipynb               # Replicated (GSE50838) vs non-replicated (GSE8059) comparison
│   ├── 06_GSEA_analysis.ipynb                    # GSEA preranked, KEGG 2021 gene sets
│   ├── 07_figure_generation.ipynb                # Assembles Figures 1-5 for the manuscript
│   └── 08_null_models_and_positive_control.ipynb # Independent-shuffle vs. NK-only null, enrichment-on-nulls,
│                                                    positive control (GSE288312) incl. background sensitivity check
├── data/                                         # Raw inputs as downloaded from GEO
│   ├── GSE50838/                                 # Series matrix + GPL11532 platform annotation
│   ├── GSE40548_gene_exp.diff.csv                # Cuffdiff output
│   └── GSE221866_DESeq2.xlsx                     # DESeq2 output (Ismailova et al., 2023)
├── tables/                                       # Processed outputs referenced in the manuscript
│   ├── convergent_genes_RRA_GSE50838_20genes.csv        # The 20-GCS genes (Table 1)
│   ├── convergent_genes_RRA_GSE50838_UP.csv / _DOWN.csv
│   ├── meta_analysis_RRA_GSE50838_all_genes.csv         # Full ranked gene list; source of the 8,581-gene background
│   ├── GSE50838_all_genes_DE.csv                        # Full GSE50838 differential expression results
│   ├── GSE50838_convergent_genes.csv
│   ├── comparison_GSE8059_vs_GSE50838.csv               # Underlies Figure 3
│   ├── ChEA2022_TF_enrichment_20GCS.csv                 # Table 2 in manuscript
│   ├── GO_enrichment_20GCS.csv, KEGG_enrichment_20GCS.csv  # Table 4 in manuscript
│   ├── GSEA_NK_IL2_KEGG.csv / GSEA_Neut_TNF_KEGG.csv / GSEA_Neut_LPS_KEGG.csv
│   ├── GSEA_convergent_pathways.csv                     # Table 3 in manuscript
│   ├── sensitivity_analysis.csv                         # Table 5 in manuscript (fold-change threshold sensitivity)
│   ├── permutation_null_comparison.csv                  # Independent-shuffle vs. NK-only null summary (Section 3.2.1)
│   ├── enrichment_on_null_sets.csv                      # Full 21-term result (Table 7, Section 3.9)
│   ├── candidates_224_TNF_LPS_concordant.csv            # The 224-gene TNF-alpha/LPS two-way concordant candidate set
│   ├── positive_control_GSE288312.csv                   # Positive-control result, both background variants (Section 3.10)
│   └── ranked_*.rnk                                     # GSEA preranked input files, also used by notebook 08
├── figures/                                      # Final figures as they appear in the manuscript
│   ├── Fig1_volcano_plots.png
│   ├── Fig2_heatmap_convergent.png
│   ├── Fig3_TF_enrichment.png
│   ├── Fig4_GSEA_pathways.png
│   ├── Fig5_sensitivity_analysis.png
│   └── Fig6A-D_*.png                             # Four panels of Figure 6
├── requirements.txt
├── LICENSE
└── README.md
```

## Installation & Usage

```bash
git clone https://github.com/patricioporcile/SIEM-meta-analysis.git
cd SIEM-meta-analysis
python -m venv venv && source venv/bin/activate   # or conda create -n nullmodel python=3.11
pip install -r requirements.txt
jupyter notebook notebooks/01_GSE50838_processing.ipynb
```

Notebooks are numbered to reflect a logical execution order (01 → 08), but 01+02 (per-dataset processing), 04 (enrichment), 05 (GSE8059 comparison), 06 (GSEA), and 08 (null models / positive control) can each be run independently once the outputs of 01–03 are available in `tables/`. Notebook 08 is the one that reproduces this manuscript's central methodological result and does not depend on 04–07.

## Reproducing the central result (independent-shuffle vs. NK-only null)

```python
import numpy as np, pandas as pd

def load_rnk(path):
    d = {}
    for line in open(path):
        g, v = line.strip().split('\t')
        d[g] = float(v)
    return d

tnf = load_rnk("tables/ranked_Neut_TNF.rnk")
lps = load_rnk("tables/ranked_Neut_LPS.rnk")
nk  = load_rnk("tables/ranked_NK_IL2.rnk")
background = set(pd.read_csv("tables/meta_analysis_RRA_GSE50838_all_genes.csv")["Gene"])

# Convergence criterion (manuscript Section 2.2.5): |log2FC| > 1, same direction, all three datasets.
# See notebooks/08_null_models_and_positive_control.ipynb for the full independent-shuffle vs.
# NK-only permutation comparison (10,000 iterations each).
```

The exact 20-gene concordant set (20-GCS), with per-dataset fold-changes and FDR values, is in `tables/convergent_genes_RRA_GSE50838_20genes.csv`, matching Table 1 of the manuscript. The permutation null comparison is in `tables/permutation_null_comparison.csv`; the enrichment-on-nulls result is in `tables/enrichment_on_null_sets.csv`; the positive-control result (both background variants) is in `tables/positive_control_GSE288312.csv`.

## Known Limitations (see manuscript Section 4.8 for full discussion)

- **Purely computational study.** No qRT-PCR or functional validation has been performed; all interpretations are hypothesis-generating.
- **Temporal mismatch.** Stimulation durations differ by up to 120-fold across datasets (IL-2: 5 days; TNF-alpha: 1h; LPS: 6h; IL-1beta: 1h).
- **Small n in the NK arm** (4 paired donors); post-hoc power analysis suggests genes with true |log2FC| in the 0.8–1.2 range may be missed.
- **GSE40548 is not a single-platform triplicate.** The 3 TNF-alpha samples combine one SOLiD-sequenced donor and two Illumina-sequenced donors.
- **No empirical batch-effect correction** (e.g., ComBat, ComBat-seq) was applied across the microarray/RNA-seq datasets; the rank-based meta-analysis mitigates but does not eliminate platform effects.
- **The rank-aggregation statistic used here is a simplified geometric-mean-of-ranks approach**, not the exact beta-order-statistics algorithm of Kolde et al. (2012, *Bioinformatics*); see `03_meta_analysis_RRA.ipynb` and manuscript Section 2.2.4 for the precise formula.
- ChEA 2022 transcription-factor enrichment is derived predominantly from ChIP-seq in cancer cell lines and embryonic stem cells, not primary immune cells.
- **The GSE288312 positive-control fold-change is a simplified CPM-based estimate**, not a full DESeq2/edgeR model, and its primary null pool (26,485 genes) differs from the 8,581-gene background used elsewhere; the sensitivity check in `tables/positive_control_GSE288312.csv` addresses this directly.

## Dependencies

See `requirements.txt`. Core packages: pandas, numpy, scipy, statsmodels, matplotlib, seaborn, gseapy, openpyxl, jupyter.

## Citation

If you use this pipeline, please cite the manuscript (citation details to be updated upon publication) and the original GEO dataset publications:

- Sanchez-Martínez et al. (2014), *Int J Biochem Cell Biol* 49:42-52 (GSE50838)
- Wright et al. (2013), *PLoS One* 8(3):e58498 (GSE40548)
- Ismailova et al. (2023), *Front Immunol* 14:1123344 (GSE221866)
- Thomas, Edwards & Wright (2025), *Discovery Immunology* 4(1):kyaf013 (GSE288312)
- Dybkaer et al. (2007), *BMC Genomics* 8:230 (GSE8059)

## License

MIT License — see [LICENSE](LICENSE).

## Contact

Patricio E. Porcile — patricio.porcile@uda.cl — Departamento de Química y Biología, Universidad de Atacama, Copiapó, Chile.

---
**Last updated**: July 2026 — corresponds to `SIEM_manuscript_v11.docx`, submitted for peer review.
