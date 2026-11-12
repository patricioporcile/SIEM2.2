# nullconv

Null models for testing whether genes are **convergent across several expression datasets**, and a check for whether the usual
"shuffle every dataset independently" null answers the question you are asking.

Companion code to: Porcile PE, Hernández Lanyon P. *Independent-shuffle null models overstate gene convergence across correlated datasets:
a conditional alternative, software, and a worked example in neutrophils and natural killer cells* (PLOS ONE, submitted).

## Two nulls, two questions
| Null | Null hypothesis | Answers |
|---|---|---|
| **independent-shuffle** (as in e.g. SuperExactTest) | all datasets are mutually independent | "Do the datasets agree at all?" |
| **conditional** (permute only the focal dataset; references stay fixed gene by gene) | the focal dataset is unrelated to the references, whatever the references share | "Does the focal dataset add convergence beyond what the references already share?" |

If the reference datasets are themselves correlated (same cell type, related stimuli, shared pipeline), the independent null's hypothesis is
false, so it rejects even when the focal dataset adds nothing: for the second question those rejections are **false positives**
(for the first question they are correct detections). With one direction, the conditional test is exactly Fisher's exact (hypergeometric) test
of whether focal-called genes are over-represented among genes called in all references (`fisher_exact_conditional`); with two directions
it is the sum of the two counts. It is a conditional test, not a new principle; the package makes it, its expectation, its false-positive
behaviour and its power easy to compute.

## Install
```bash
pip install -e .          # needs numpy
pip install -e .[test]    # + pytest, pandas, scipy
```

## Quick start
```python
import nullconv as nc
res = nc.convergence_test(focal, [ref1, ref2], thr=1.0)     # log2FC vectors over the SAME genes; or rank_q=0.1 for platform-free calls
print(res.summary())                                        # both nulls, analytic expectation, dependence warning
nc.convergence_test_all_foci({"A": a, "B": b, "C": c})      # each dataset in turn as focal (symmetric questions)
nc.diagnose_dependence([ref1, ref2], focal_calls, ref_calls)  # correlation among references, expected-count ratio
nc.false_positive_rate(focal_calls, ref_calls)              # rejection rate of each null if the focal dataset added nothing
nc.power_curve(focal_calls, ref_calls, k_values=(0, 5, 10, 15), p_detect=1.0)
nc.fisher_exact_conditional(focal_calls, ref_calls, "up")   # exact hypergeometric p for one direction
```
Calls come from `calls_from_fc(x, thr)` or `calls_from_ranks(x, q)`.

## Which null should I use? (decision guide)
1. **State the question.** "Do the datasets agree at all?" → independent-shuffle null. "Does dataset F add convergence beyond datasets R1…Rk?" → conditional null with F focal.
2. **Choose the focal dataset from the question, not the data.** If the question is symmetric or no dataset is natural, run `convergence_test_all_foci` and interpret each result as
   "does this dataset add beyond the others?"; a dataset correlated with one of the others will look significant for that reason alone.
3. **Diagnose dependence among the references** (`diagnose_dependence`): correlation and ratio of expected counts. If the ratio is near 1 the two nulls agree.
4. **Report the analytic expectation** (`analytic_expectation`) next to the permutation mean.
5. **Check false-positive rate and power** (`false_positive_rate`, `power_curve`) before reading a non-significant result as absence of shared genes.
6. **Check robustness**: rank calls (`rank_q`), several thresholds, and a positive control or a semisynthetic one (plant k shared genes: `power_curve`).
7. **Limits**: both nulls permute gene labels and ignore gene-gene correlation; borderline results are exploratory. Related methods that avoid thresholds or model gene
   correlation include RRHO (Plaisier et al. 2010) and camera (Wu & Smyth 2012).

## Tests
```bash
pytest -q                                   # 14 synthetic tests, 15 with the paper table
NULLCONV_PAPER_TABLE=path/to/meta_analysis_RRA_GSE50838_all_genes.csv pytest -q   # also reproduces the manuscript numbers
```
## License
MIT (same as the parent repository).
