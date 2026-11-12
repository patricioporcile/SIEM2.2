"""Reproduces the worked example of the manuscript (NK cells + two neutrophil datasets).
Usage: python neutrophil_nk_example.py path/to/meta_analysis_RRA_GSE50838_all_genes.csv"""
import sys, pandas as pd, nullconv as nc
m = pd.read_csv(sys.argv[1])
focal, refs = m.NK_IL2_FC_50838.values, [m.Neut_TNF_FC.values, m.Neut_LPS_FC.values]
res = nc.convergence_test(focal, refs, thr=1.0, n_perm=10_000, seed=42)
print(res.summary())
f, R = nc.calls_from_fc(focal), [nc.calls_from_fc(r) for r in refs]
print("\nFalse-positive rate if the focal dataset had no convergence:", nc.false_positive_rate(f, R, seed=0))
print("Power (k shared genes, all detected):", nc.power_curve(f, R, k_values=(0, 5, 10, 15), n_sim=200, seed=0))

print("\nEach dataset as focal (others as references):")
vals = {"NK": focal, "TNF": refs[0], "LPS": refs[1]}
for k, r in nc.convergence_test_all_foci(vals, thr=1.0, n_perm=10_000, seed=42).items():
    print(f"  focal = {k}: observed {r.observed}, expected {r.expected_conditional:.1f}, p_conditional = {r.p_conditional:.4g}")
