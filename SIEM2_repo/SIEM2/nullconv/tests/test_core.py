import numpy as np, pytest, os, pandas as pd
import nullconv as nc

rng = np.random.default_rng(1)

def make(n=4000, r=0.5, seed=3):
    g = np.random.default_rng(seed)
    z1 = g.standard_normal(n); z2 = r * z1 + np.sqrt(1 - r**2) * g.standard_normal(n)
    focal = g.standard_normal(n)
    return focal * 1.2, z1 * 1.3, z2 * 1.3

def test_calls():
    x = np.array([2, -2, 0.5, np.nan]); assert list(nc.calls_from_fc(x, 1)) == [1, -1, 0, 0]
    q = nc.calls_from_ranks(np.arange(100), 0.1); assert (q == 1).sum() == 10 and (q == -1).sum() == 10

def test_input_validation():
    with pytest.raises(ValueError): nc.convergence_count(np.zeros(5), [np.zeros(4)])
    with pytest.raises(ValueError): nc.convergence_count(np.array([2, 0]), [np.array([1, 0])])
    with pytest.raises(ValueError): nc.convergence_count(np.zeros(3), [])

def test_p_never_zero():
    assert nc.empirical_p(100, np.zeros(999)) == pytest.approx(1 / 1000)

def test_analytic_matches_permutation_mean():
    f, a, b = make(); F, A, B = nc.calls_from_fc(f), nc.calls_from_fc(a), nc.calls_from_fc(b)
    null = nc.null_correlation_preserving(F, [A, B], 20000, seed=1)
    assert null.mean() == pytest.approx(nc.analytic_expectation(F, [A, B]), rel=0.03)

def test_corr_preserving_equals_explicit_permutation():
    f, a, b = make(2000); F, A, B = nc.calls_from_fc(f), nc.calls_from_fc(a), nc.calls_from_fc(b)
    fast = nc.null_correlation_preserving(F, [A, B], 6000, seed=2)
    g = np.random.default_rng(5)
    slow = np.array([nc.convergence_count(g.permutation(F), [A, B]) for _ in range(6000)])
    assert fast.mean() == pytest.approx(slow.mean(), abs=0.25) and fast.std() == pytest.approx(slow.std(), abs=0.25)

def test_independent_hypergeometric_matches_permutation():
    f, a, b = make(2000, r=0.0); F, A, B = nc.calls_from_fc(f), nc.calls_from_fc(a), nc.calls_from_fc(b)
    h = nc.null_independent(F, [A, B], 6000, seed=2)
    p = nc.null_independent(F, [A, B], 3000, seed=3, method="permutation")
    assert h.mean() == pytest.approx(p.mean(), abs=0.2) and h.std() == pytest.approx(p.std(), abs=0.2)

def test_independent_null_is_anticonservative_when_refs_correlated():
    f, a, b = make(r=0.6); F, A, B = nc.calls_from_fc(f), nc.calls_from_fc(a), nc.calls_from_fc(b)
    fp = nc.false_positive_rate(F, [A, B], n_null=4000, seed=0)
    assert fp["independent"] > 0.5 and fp["correlation_preserving"] <= 0.06

def test_both_nulls_agree_when_refs_independent():
    f, a, b = make(r=0.0); F, A, B = nc.calls_from_fc(f), nc.calls_from_fc(a), nc.calls_from_fc(b)
    fp = nc.false_positive_rate(F, [A, B], n_null=4000, seed=0)
    assert fp["independent"] < 0.09 and fp["correlation_preserving"] <= 0.06

def test_power_increases_with_k():
    f, a, b = make(r=0.4); F, A, B = nc.calls_from_fc(f), nc.calls_from_fc(a), nc.calls_from_fc(b)
    pw = nc.power_curve(F, [A, B], k_values=(0, 5, 30), n_sim=120, n_perm=400, seed=0)
    assert pw[0] <= 0.12 and pw[30] > pw[5] >= pw[0] and pw[30] > 0.6

def test_diagnose_and_warning():
    f, a, b = make(r=0.6)
    res = nc.convergence_test(f, [a, b], thr=1.0, n_perm=2000, seed=0)
    assert any("dependent" in w for w in res.warnings)
    d = nc.diagnose_dependence([a, b]); assert d["pearson"][(0, 1)] == pytest.approx(0.6, abs=0.05)
    f2, a2, b2 = make(r=0.0)
    assert not nc.convergence_test(f2, [a2, b2], n_perm=500, seed=0).warnings

def test_rank_based_calls_run():
    f, a, b = make(r=0.3)
    res = nc.convergence_test(f, [a, b], rank_q=0.1, n_perm=1000, seed=0)
    assert res.p_corr_preserving > 0.01

PAPER = os.environ.get("NULLCONV_PAPER_TABLE", "")
@pytest.mark.skipif(not os.path.exists(PAPER), reason="set NULLCONV_PAPER_TABLE to meta_analysis_RRA_GSE50838_all_genes.csv")
def test_reproduces_paper_numbers():
    m = pd.read_csv(PAPER)
    r = nc.convergence_test(m.NK_IL2_FC_50838.values, [m.Neut_TNF_FC.values, m.Neut_LPS_FC.values], n_perm=10000, seed=42)
    assert r.observed == 20 and r.p_independent < 0.001
    assert r.expected_corr_preserving == pytest.approx(20.66, abs=0.02)
    assert 0.5 < r.p_corr_preserving < 0.7
