"""Null models for testing multi-dataset gene convergence.

Question answered: given one *focal* dataset and one or more *reference* datasets,
do more genes show concordant, large changes in all of them than expected by chance?

Two null models are provided:

``independent``
    Every dataset is shuffled independently across gene labels. It is valid only if the
    datasets are independent under the null. If the reference datasets are correlated
    (e.g. two neutrophil datasets under related stimuli) it is far too permissive.

``correlation_preserving``
    The reference datasets are kept paired gene-for-gene (their real dependence is
    preserved) and only the focal dataset is permuted. It tests whether the focal dataset
    adds convergence beyond what the references already share.

Everything works on "calls" (+1 up, -1 down, 0 not called), obtained either from a
fold-change threshold (:func:`calls_from_fc`) or from within-dataset rank tails
(:func:`calls_from_ranks`, platform-free).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence

import numpy as np

__all__ = [
    "calls_from_fc", "calls_from_ranks", "convergence_count", "analytic_expectation",
    "null_independent", "null_correlation_preserving", "empirical_p",
    "diagnose_dependence", "false_positive_rate", "power_curve", "convergence_test",
    "ConvergenceResult",
]


# ----------------------------------------------------------------------------- calls
def calls_from_fc(x, thr: float = 1.0) -> np.ndarray:
    """+1 if x > thr, -1 if x < -thr, else 0 (NaN -> 0)."""
    x = np.asarray(x, dtype=float)
    return np.where(x > thr, 1, np.where(x < -thr, -1, 0)).astype(np.int8)


def calls_from_ranks(x, q: float = 0.10) -> np.ndarray:
    """+1 for the top fraction q of values, -1 for the bottom fraction q, else 0.

    Ranks are computed within the dataset, so the calls do not depend on the scale of
    fold changes (useful when mixing microarray and RNA-seq)."""
    x = np.asarray(x, dtype=float)
    if not 0 < q < 0.5:
        raise ValueError("q must be in (0, 0.5)")
    order = np.argsort(np.argsort(np.nan_to_num(x, nan=0.0)))
    r = (order + 0.5) / len(x)
    return np.where(r > 1 - q, 1, np.where(r < q, -1, 0)).astype(np.int8)


def _prep(focal, refs):
    f = np.asarray(focal)
    R = [np.asarray(r) for r in (refs if isinstance(refs, (list, tuple)) else [refs])]
    if not R:
        raise ValueError("at least one reference dataset is required")
    if any(len(r) != len(f) for r in R):
        raise ValueError("all vectors must have the same length (shared gene background)")
    if not set(np.unique(np.concatenate([f] + R))) <= {-1, 0, 1}:
        raise ValueError("inputs must be calls in {-1, 0, 1}; use calls_from_fc/calls_from_ranks")
    return f, R


def _ref_candidates(R):
    up = np.logical_and.reduce([r == 1 for r in R])
    dn = np.logical_and.reduce([r == -1 for r in R])
    return int(up.sum()), int(dn.sum()), up, dn


# ------------------------------------------------------------------------ statistics
def convergence_count(focal, refs) -> int:
    """Number of genes called in the same direction in the focal and all reference datasets."""
    f, R = _prep(focal, refs)
    _, _, up, dn = _ref_candidates(R)
    return int(((f == 1) & up).sum() + ((f == -1) & dn).sum())


def analytic_expectation(focal, refs) -> float:
    """Closed-form expected count under the correlation-preserving null:
    n_up * P(focal up) + n_dn * P(focal down), n_up/n_dn = genes concordant in all references."""
    f, R = _prep(focal, refs)
    n = len(f)
    nu, nd, _, _ = _ref_candidates(R)
    return nu * (f == 1).sum() / n + nd * (f == -1).sum() / n


def _expectation_independent(f, R):
    n = len(f)
    e = 0.0
    for s in (1, -1):
        p = (f == s).sum() / n
        for r in R:
            p *= (r == s).sum() / n
        e += n * p
    return e


def empirical_p(observed: int, null) -> float:
    """One-sided empirical p = (b + 1) / (B + 1); never exactly zero (Phipson & Smyth 2010)."""
    null = np.asarray(null)
    return float((1 + (null >= observed).sum()) / (1 + len(null)))


# ---------------------------------------------------------------------------- nulls
def null_correlation_preserving(focal, refs, n_perm: int = 10_000, seed=None) -> np.ndarray:
    """Null counts when only the focal vector is permuted (references stay paired).

    Sampled exactly from the underlying (multivariate) hypergeometric distribution, which is
    identical to permuting the focal labels but much faster."""
    f, R = _prep(focal, refs)
    rng = np.random.default_rng(seed)
    n = len(f)
    nu, nd, _, _ = _ref_candidates(R)
    a, b = int((f == 1).sum()), int((f == -1).sum())
    pool = np.array([a, b, n - a - b])
    d1 = rng.multivariate_hypergeometric(pool, nu, size=n_perm)       # focal values drawn for "up" candidates
    rem_down = pool[1] - d1[:, 1]
    d2 = rng.hypergeometric(rem_down, (n - nu) - rem_down, nd)          # "down" among those drawn for "down" candidates
    return (d1[:, 0] + d2).astype(int)


def null_independent(focal, refs, n_perm: int = 10_000, seed=None, method: str = "hypergeometric") -> np.ndarray:
    """Null counts when every dataset is shuffled independently across genes.

    ``method='hypergeometric'`` draws the count from the exact distribution of the
    intersection of independent random subsets, separately for the up and down directions (the
    tiny dependence between the two directions is ignored; verified against explicit
    permutation in the test-suite). ``method='permutation'`` shuffles explicitly."""
    f, R = _prep(focal, refs)
    rng = np.random.default_rng(seed)
    n = len(f)
    if method == "permutation":
        out = np.empty(n_perm, int)
        for i in range(n_perm):
            P = [rng.permutation(v) for v in [f] + R]
            out[i] = ((np.logical_and.reduce([p == 1 for p in P])).sum()
                      + (np.logical_and.reduce([p == -1 for p in P])).sum())
        return out
    if method != "hypergeometric":
        raise ValueError("method must be 'hypergeometric' or 'permutation'")
    out = np.zeros(n_perm, int)
    for s in (1, -1):
        x = np.full(n_perm, int((f == s).sum()))
        for r in R:
            a = int((r == s).sum())
            x = rng.hypergeometric(x, n - x, a)
        out += x
    return out


# -------------------------------------------------------------------------- diagnosis
def diagnose_dependence(ref_values: Sequence, focal_calls=None, ref_calls=None) -> Dict:
    """Quantify dependence among reference datasets.

    ``ref_values``: list of continuous vectors (e.g. log2 fold changes) over the shared genes.
    Returns pairwise Pearson/Spearman correlations. If ``focal_calls`` and ``ref_calls`` are
    given, also returns ``inflation``: expected count under the correlation-preserving null divided
    by the expected count under the independent null (>>1 means the independent null is badly
    mis-specified)."""
    V = [np.asarray(v, float) for v in ref_values]
    out = {"pearson": {}, "spearman": {}}
    for i in range(len(V)):
        for j in range(i + 1, len(V)):
            m = np.isfinite(V[i]) & np.isfinite(V[j])
            out["pearson"][(i, j)] = float(np.corrcoef(V[i][m], V[j][m])[0, 1])
            ri = np.argsort(np.argsort(V[i][m])); rj = np.argsort(np.argsort(V[j][m]))
            out["spearman"][(i, j)] = float(np.corrcoef(ri, rj)[0, 1])
    if focal_calls is not None and ref_calls is not None:
        f, R = _prep(focal_calls, ref_calls)
        e_c, e_i = analytic_expectation(f, R), _expectation_independent(f, R)
        out["expected_corr_preserving"], out["expected_independent"] = float(e_c), float(e_i)
        out["inflation"] = float(e_c / e_i) if e_i > 0 else float("inf")
    return out


def false_positive_rate(focal, refs, alpha: float = 0.05, n_null: int = 5000, seed=None) -> Dict:
    """How often would each null declare 'significant convergence' if the focal dataset carried
    NO real convergence (its labels shuffled)? For the observed dependence structure.

    If the focal labels are shuffled, the observed count is distributed exactly as the
    correlation-preserving null, so the rate for the independent null is the share of those draws
    that exceed its alpha critical value. The correlation-preserving null stays at or below alpha."""
    f, R = _prep(focal, refs)
    rng = np.random.default_rng(seed)
    cp = null_correlation_preserving(f, R, n_null, rng.integers(2**31))
    ind = null_independent(f, R, n_null, rng.integers(2**31))
    p_ind = np.array([empirical_p(o, ind) for o in cp])
    # cp null calibration: p of each draw against the same distribution
    srt = np.sort(cp)
    p_cp = (1 + (len(cp) - np.searchsorted(srt, cp, side="left"))) / (1 + len(cp))
    return {"alpha": alpha, "independent": float((p_ind <= alpha).mean()),
            "correlation_preserving": float((p_cp <= alpha).mean())}


def power_curve(focal, refs, k_values=(0, 2, 5, 10, 15), p_detect: float = 1.0, alpha: float = 0.05,
                n_sim: int = 300, n_perm: int = 1000, seed=None) -> Dict[int, float]:
    """Power of the correlation-preserving test when k truly shared genes (drawn from the
    reference-concordant candidates) are planted in the focal dataset; each planted gene is
    called in the focal dataset with probability ``p_detect`` (limited sensitivity of small designs)."""
    f, R = _prep(focal, refs)
    rng = np.random.default_rng(seed)
    nu, nd, up, dn = _ref_candidates(R)
    cand = np.where(up | dn)[0]
    sign = np.where(up, 1, -1)
    res = {}
    for k in k_values:
        hits = 0
        for _ in range(n_sim):
            g = rng.permutation(f)   # focal background with no true convergence (labels shuffled)
            if k:
                idx = rng.choice(cand, min(k, len(cand)), replace=False)
                idx = idx[rng.random(len(idx)) < p_detect]
                g[idx] = sign[idx]
            obs = convergence_count(g, R)
            null = null_correlation_preserving(g, R, n_perm, rng.integers(2**31))
            hits += empirical_p(obs, null) <= alpha
        res[int(k)] = hits / n_sim
    return res


# ------------------------------------------------------------------------ main entry
@dataclass
class ConvergenceResult:
    observed: int
    expected_corr_preserving: float
    null_independent_mean: float
    null_independent_sd: float
    p_independent: float
    null_corr_preserving_mean: float
    null_corr_preserving_sd: float
    p_corr_preserving: float
    n_genes: int
    n_candidates_up: int
    n_candidates_down: int
    ref_correlations: Dict = field(default_factory=dict)
    warnings: List[str] = field(default_factory=list)

    def summary(self) -> str:
        s = [f"observed convergent genes: {self.observed} (of {self.n_genes} shared genes)",
             f"independent-shuffle null:       {self.null_independent_mean:.1f} +/- {self.null_independent_sd:.1f}, p = {self.p_independent:.4g}",
             f"correlation-preserving null:    {self.null_corr_preserving_mean:.1f} +/- {self.null_corr_preserving_sd:.1f}"
             f" (analytic {self.expected_corr_preserving:.1f}), p = {self.p_corr_preserving:.4g}"]
        s += ["WARNING: " + w for w in self.warnings]
        return "\n".join(s)


def convergence_test(focal_values, ref_values: Sequence, thr: float = 1.0, rank_q: Optional[float] = None,
                     n_perm: int = 10_000, seed: Optional[int] = 0) -> ConvergenceResult:
    """Run both nulls on continuous inputs (e.g. log2 fold changes over the SAME genes).

    Use ``rank_q`` (e.g. 0.1) to call by within-dataset rank tails instead of the threshold ``thr``.
    The correlation-preserving result is the one to interpret; the independent result is returned
    for comparison, together with a warning when the reference datasets are dependent."""
    refs_list = list(ref_values)
    mk = (lambda v: calls_from_ranks(v, rank_q)) if rank_q else (lambda v: calls_from_fc(v, thr))
    f = mk(focal_values); R = [mk(v) for v in refs_list]
    obs = convergence_count(f, R)
    rng = np.random.default_rng(seed)
    ni = null_independent(f, R, n_perm, rng.integers(2**31))
    nc = null_correlation_preserving(f, R, n_perm, rng.integers(2**31))
    nu, nd, _, _ = _ref_candidates(R)
    dep = diagnose_dependence(refs_list, f, R)
    warns = []
    corr = [abs(v) for v in dep["pearson"].values()]
    if dep.get("inflation", 1) > 1.5 or (corr and max(corr) > 0.2):
        warns.append("reference datasets are dependent (max |r| = %.2f; expected count is %.1fx higher than under "
                     "independence). The independent-shuffle p-value is anti-conservative; interpret p_corr_preserving."
                     % (max(corr) if corr else float('nan'), dep.get("inflation", float('nan'))))
    return ConvergenceResult(obs, float(analytic_expectation(f, R)), float(ni.mean()), float(ni.std()),
                             empirical_p(obs, ni), float(nc.mean()), float(nc.std()), empirical_p(obs, nc),
                             len(f), nu, nd, dep["pearson"], warns)
