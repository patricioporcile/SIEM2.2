# nullconv

Null models for testing whether genes are **convergent across several expression datasets** — and a check for whether the usual
"shuffle every dataset independently" null is valid for your data.

Companion code to: Porcile PE, Hernández Lanyon P. *Null model specification determines whether cross-lineage transcriptional
convergence is detected: a worked example in neutrophils and natural killer cells* (PLOS ONE, submitted).

## The problem
A common test asks whether the number of genes that change in the same direction in all datasets exceeds chance, and builds the null by
shuffling each dataset's gene labels **independently**. That is valid only if the datasets are independent under the null. If two or more
datasets are correlated for reasons unrelated to your question (same cell type, related stimuli, shared pipeline), the independent null
is too narrow and an unremarkable overlap looks highly significant. In the manuscript example the independent null gave p < 0.0001 and the
correlation-preserving null p = 0.60 for the identical overlap; at the observed correlation (r = 0.39) the independent null falsely
rejected in >90% of simulations.

## Install
```bash
pip install -e .          # from this folder; needs numpy
pip install -e .[test]    # + pytest, pandas
```

## Quick start
```python
import nullconv as nc
# focal: log2FC of the dataset you are testing for added convergence; refs: list of log2FC vectors over the SAME genes
res = nc.convergence_test(focal, [ref1, ref2], thr=1.0)        # or rank_q=0.1 for platform-free rank calls
print(res.summary())                                           # both nulls, analytic expectation, dependence warning
nc.diagnose_dependence([ref1, ref2], focal_calls, ref_calls)   # correlation among references, expected-count inflation
nc.false_positive_rate(focal_calls, ref_calls)                 # false-positive rate of each null if focal had no signal
nc.power_curve(focal_calls, ref_calls, k_values=(0, 5, 10, 15), p_detect=1.0)
```
Calls come from `calls_from_fc(x, thr)` or `calls_from_ranks(x, q)`.

## Which null should I use? (decision guide)
1. **Define roles.** One *focal* dataset (the one whose added convergence you test) and one or more *reference* datasets.
2. **Diagnose dependence among the references** (`diagnose_dependence`): pairwise correlation over the shared genes and the ratio of expected
   counts (correlation-preserving / independent). Ratio near 1 and |r| < ~0.05 → both nulls agree; use either.
3. **If the references are dependent → use the correlation-preserving null** (permute only the focal dataset). The independent null
   answers a different question ("do the datasets overlap at all given their marginals?") and will be anti-conservative here.
4. **Report the analytic expectation** (`analytic_expectation`) next to the permutation mean: it shows that the size of the reference-concordant set sets the scale.
5. **Check power** (`power_curve`) before reading a non-significant result as absence of shared genes.
6. **Check robustness**: repeat with rank calls (`rank_q`) to remove platform-scale effects, and across thresholds.
7. **Remember the limits**: both nulls permute gene labels and ignore gene-gene correlation (co-regulated modules); this can make p-values too small.
   If the result is borderline, treat it as exploratory.

## Tests
```bash
pytest -q                                   # 12 tests; synthetic data
NULLCONV_PAPER_TABLE=path/to/meta_analysis_RRA_GSE50838_all_genes.csv pytest -q   # also reproduces the manuscript numbers
```
## License
MIT (same as the parent repository).
