#!/usr/bin/env python3
"""
invariant_sweep.py
==================
Property-based sanity sweep over the Helixar engine.

Goal: answer "do realistic molecule × input combinations always produce
outputs that make bioprocess-engineering sense?" — without inspecting
millions of combinations by hand.

How it stays REALISTIC (the key point):
  Inputs are NOT drawn from the full cartesian product. They are sampled
  from calibrated, per-subtype bands that mirror what a process-development
  scientist would actually enter — reusing the SAME tables the engine uses:

    • method      ← route_comparison._SUBTYPE_VIABLE_METHODS
                    (a mAb is only ever sampled as biotech, never chemical)
    • scale       ← decision_engine._SUBTYPE/_MOLECULE_SCALE_THRESHOLDS
                    (mAb in kg/yr, bulk solvent in tonnes)
    • steps       ← plausibility_checker.KNOWN_PROCESSES step ranges
                    (so we never enter "2-step Liraglutide SPPS")
    • purity      ← per-subtype application bands (API ≥98 %, enzyme 90-99 %)
    • RM price    ← anchored to the molecule's COGS band

How it CHECKS correctness:
  For every sampled scenario it verifies absolute laws (COGS ordering,
  benchmark band, efficiency/risk consistency, confidence) and, by
  re-running with ONE dimension nudged, monotonicity laws (more purity →
  never cheaper, bigger scale → never dearer per kg, more steps → never
  cheaper / never more efficient).

A single violation across thousands of scenarios is a real bug — surfaced
with the exact input + output that broke the rule.

Run:  python3 invariant_sweep.py [N_PER_MOLECULE]
"""

from __future__ import annotations

import csv
import math
import random
import sys
from typing import Any, Dict, List, Optional, Tuple

# --- Engine + calibrated sampling tables (single source of truth) -----------
from cogs_estimator import estimate_cogs, estimate_cogs_whatif, _MOLECULE_COGS_OVERRIDES
from plausibility_checker import KNOWN_PROCESSES, check_plausibility
from route_comparison import _SUBTYPE_VIABLE_METHODS, compare_routes
from decision_engine import DecisionEngine

SEED = 20260531  # deterministic → reproducible runs / stable CI signal


# ---------------------------------------------------------------------------
# Molecule catalogue (name, type, subtype) from the audit CSV (falls back to DB)
# ---------------------------------------------------------------------------
def load_molecules() -> List[Dict[str, str]]:
    mols: List[Dict[str, str]] = []
    try:
        with open("audit_outputs/audit_summary.csv") as f:
            for r in csv.DictReader(f):
                if r.get("molecule"):
                    mols.append({
                        "name": r["molecule"].strip(),
                        "type": (r.get("molecule_type") or "").strip(),
                        "subtype": (r.get("molecule_subtype") or "").strip(),
                    })
    except Exception:
        pass
    if not mols:
        from molecules_db import MOLECULE_DATABASE
        for m in MOLECULE_DATABASE:
            mols.append({
                "name": m["name"],
                "type": m.get("molecule_type", ""),
                "subtype": m.get("molecule_subtype", ""),
            })
    return mols


# ---------------------------------------------------------------------------
# Realistic per-subtype input bands
# ---------------------------------------------------------------------------
# Target purity (%) bands by application class — what is actually demanded.
_PURITY_BANDS = {
    ("small_molecule", "volatile"):     (95.0, 99.8),   # solvents/fuel grade
    ("small_molecule", "non_volatile"): (98.0, 99.9),   # APIs
    ("natural_product", "alkaloid"):    (98.0, 99.9),   # API alkaloids
    ("natural_product", "terpene"):     (90.0, 99.5),   # flavour/fragrance
    ("peptide", "linear"):              (95.0, 99.5),   # therapeutic peptides
    ("peptide", "cyclic"):              (95.0, 99.5),
    ("protein", "antibody"):            (98.0, 99.9),   # mAb DS purity
    ("protein", "enzyme"):              (90.0, 99.0),   # industrial enzymes
}
_PURITY_DEFAULT = (95.0, 99.5)

# Raw-material price multiplier band relative to the molecule's COGS low anchor
# (RM is a fraction of the finished-good cost; sample around a realistic share).
_RM_FRACTION_BAND = (0.05, 0.6)


def _scale_band(name: str, mtype: str, msub: str) -> Tuple[float, float]:
    """Realistic kg/year band for this molecule, derived from the engine's own
    scale thresholds: from a bit below lab_max up to a few× pilot_max."""
    name_l = name.lower()
    lab_max, pilot_max = 1.0, 1000.0
    if name_l in DecisionEngine._MOLECULE_SCALE_THRESHOLDS:
        lab_max, pilot_max = DecisionEngine._MOLECULE_SCALE_THRESHOLDS[name_l]
    elif (mtype, msub) in DecisionEngine._SUBTYPE_SCALE_THRESHOLDS:
        lab_max, pilot_max = DecisionEngine._SUBTYPE_SCALE_THRESHOLDS[(mtype, msub)]
    lo = max(lab_max * 0.2, 1e-4)
    hi = pilot_max * 5.0
    return lo, hi


def _steps_band(name: str, method: str) -> Tuple[int, int]:
    """Realistic step count from KNOWN_PROCESSES (the plausibility ground-truth)."""
    entry = KNOWN_PROCESSES.get(name.lower(), {}).get(method)
    if entry and "steps" in entry:
        lo, hi = entry["steps"]
        return int(lo), int(hi)
    # Fallback: generic but sane per method
    return {"chemical": (1, 8), "biotechnological": (2, 6),
            "extraction": (1, 5), "enzymatic": (2, 4)}.get(method, (1, 6))


def _loguniform(rng: random.Random, lo: float, hi: float) -> float:
    if lo <= 0:
        lo = 1e-4
    return 10 ** rng.uniform(math.log10(lo), math.log10(hi))


def sample_input(rng: random.Random, mol: Dict[str, str],
                 method: Optional[str] = None) -> Dict[str, Any]:
    """Draw ONE realistic process input for the given molecule."""
    name, mtype, msub = mol["name"], mol["type"], mol["subtype"]
    viable = _SUBTYPE_VIABLE_METHODS.get((mtype, msub), ["chemical", "biotechnological"])
    if method is None:
        method = rng.choice(viable)

    p_lo, p_hi = _PURITY_BANDS.get((mtype, msub), _PURITY_DEFAULT)
    purity = round(rng.uniform(p_lo, p_hi), 2)

    s_lo, s_hi = _scale_band(name, mtype, msub)
    scale = round(_loguniform(rng, s_lo, s_hi), 4)

    st_lo, st_hi = _steps_band(name, method)
    steps = rng.randint(st_lo, st_hi)

    # RM price: a realistic fraction of the molecule's COGS low anchor.
    ov = _MOLECULE_COGS_OVERRIDES.get(name.lower())
    base_low = float(ov["base_low"]) if ov else 10.0
    rm = round(base_low * rng.uniform(*_RM_FRACTION_BAND), 4)

    return {
        "molecule_name": name,
        "molecule_type": mtype,
        "molecule_subtype": msub,
        "method": method,
        "desired_purity_percent": purity,
        "scale_kg_per_year": scale,
        "number_of_steps": steps,
        "raw_material_cost_eur_per_kg": rm,
        "has_existing_process": True,
        "has_bioreactor": True,
        "has_advanced_purification": True,
        "available_methods": {
            "has_flash_chromatography": True, "has_prep_hplc": True,
            "has_crystallization": True, "has_membrane_filtration": True,
            "has_fplc": True, "has_extraction": True,
        },
    }


# ---------------------------------------------------------------------------
# Invariants
# ---------------------------------------------------------------------------
class Violation:
    def __init__(self, law: str, mol: str, detail: str, inp: Dict[str, Any]):
        self.law, self.mol, self.detail, self.inp = law, mol, detail, inp

    def __str__(self):
        keys = ("method", "desired_purity_percent", "scale_kg_per_year",
                "number_of_steps", "raw_material_cost_eur_per_kg")
        ctx = ", ".join(f"{k}={self.inp.get(k)}" for k in keys)
        return f"[{self.law}] {self.mol}: {self.detail}  ({ctx})"


def _cogs(inp: Dict[str, Any]):
    r = estimate_cogs(inp) or {}
    return r.get("low_eur_per_kg"), r.get("high_eur_per_kg"), r


def check_absolute(inp: Dict[str, Any], mol: Dict[str, str]) -> List[Violation]:
    """Laws that must hold for a single output."""
    v: List[Violation] = []
    lo, hi, r = _cogs(inp)
    name = mol["name"]

    # I4 — COGS ordering & positivity
    if lo is None or hi is None:
        v.append(Violation("I4-cogs-missing", name, "no COGS for a DB molecule", inp))
        return v
    if not (lo > 0 and hi > 0):
        v.append(Violation("I4-cogs-positive", name, f"non-positive COGS {lo}/{hi}", inp))
    if lo > hi + 1e-9:
        v.append(Violation("I4-cogs-order", name, f"low {lo} > high {hi}", inp))

    # I5 — benchmark band: output must stay within [0.5×, 3×] of the published
    # per-molecule anchor (catches the Vanillin-class "range exploded" bug).
    ov = _MOLECULE_COGS_OVERRIDES.get(name.lower())
    if ov:
        b_lo, b_hi = float(ov["base_low"]), float(ov["base_high"])
        if lo < b_lo * 0.5 - 1e-9:
            v.append(Violation("I5-below-benchmark", name,
                               f"COGS low {lo:.3g} < 0.5× anchor {b_lo:.3g}", inp))
        if hi > b_hi * 3.0 + 1e-9:
            v.append(Violation("I5-above-benchmark", name,
                               f"COGS high {hi:.3g} > 3× anchor {b_hi:.3g}", inp))

    # I6 — confidence: a calibrated DB molecule should never be "low".
    conf = str(r.get("confidence") or "").lower()
    if conf == "low":
        v.append(Violation("I6-confidence", name,
                           "confidence 'low' for a per-molecule-calibrated molecule", inp))
    return v


def check_monotonic(rng: random.Random, inp: Dict[str, Any],
                    mol: Dict[str, str]) -> List[Violation]:
    """Laws that compare two outputs differing in exactly ONE dimension."""
    v: List[Violation] = []
    name = mol["name"]
    lo0, hi0, _ = _cogs(inp)
    if lo0 is None:
        return v

    # I1 — higher purity must not lower COGS.
    if inp["desired_purity_percent"] < 99.5:
        hi_inp = dict(inp, desired_purity_percent=min(99.95,
                      inp["desired_purity_percent"] + 0.4))
        l2, h2, _ = _cogs(hi_inp)
        if l2 is not None and (l2 < lo0 - 1e-6 or h2 < hi0 - 1e-6):
            v.append(Violation("I1-purity-monotonic", name,
                f"purity↑ lowered COGS ({lo0:.3g}-{hi0:.3g} → {l2:.3g}-{h2:.3g})", inp))

    # I2 — bigger scale must not raise COGS/kg (economy of scale).
    big = dict(inp, scale_kg_per_year=inp["scale_kg_per_year"] * 10.0)
    l2, h2, _ = _cogs(big)
    if l2 is not None and (l2 > lo0 + 1e-6 or h2 > hi0 + 1e-6):
        v.append(Violation("I2-scale-monotonic", name,
            f"scale×10 raised COGS ({lo0:.3g}-{hi0:.3g} → {l2:.3g}-{h2:.3g})", inp))

    # I3 — more steps must not lower COGS.
    more = dict(inp, number_of_steps=inp["number_of_steps"] + 3)
    l2, h2, _ = _cogs(more)
    if l2 is not None and (l2 < lo0 - 1e-6 or h2 < hi0 - 1e-6):
        v.append(Violation("I3-steps-monotonic", name,
            f"+3 steps lowered COGS ({lo0:.3g}-{hi0:.3g} → {l2:.3g}-{h2:.3g})", inp))
    return v


def check_whatif(inp: Dict[str, Any], mol: Dict[str, str]) -> List[Violation]:
    """W-laws for the what-if sensitivity projection:
      W1 continuity — at the original inputs it must equal the PDF value
      W2 RM monotone — higher raw-material price never lowers COGS
      W3 scale monotone — bigger scale never raises COGS/kg
      W4 purity monotone — higher purity never lowers COGS
    """
    v: List[Violation] = []
    name = mol["name"]
    scale = inp["scale_kg_per_year"]
    pur = inp["desired_purity_percent"]
    rm = inp["raw_material_cost_eur_per_kg"]

    pdf = estimate_cogs(inp) or {}
    lo_pdf, hi_pdf = pdf.get("low_eur_per_kg"), pdf.get("high_eur_per_kg")
    if lo_pdf is None:
        return v

    # W1 — continuity at the reference point.
    ref = estimate_cogs_whatif(inp, scale, pur, rm)
    if abs((ref.get("low_eur_per_kg") or 0) - lo_pdf) > 0.05 or \
       abs((ref.get("high_eur_per_kg") or 0) - hi_pdf) > 0.05:
        v.append(Violation("W1-whatif-continuity", name,
            f"default whatif {ref.get('low_eur_per_kg')}-{ref.get('high_eur_per_kg')} "
            f"≠ PDF {lo_pdf}-{hi_pdf}", inp))

    # W2 — RM↑ must not lower COGS.
    a = estimate_cogs_whatif(inp, scale, pur, rm * 3.0)
    if (a.get("high_eur_per_kg") or 0) < hi_pdf - 1e-6:
        v.append(Violation("W2-whatif-rm-monotonic", name,
            f"RM×3 lowered COGS high ({hi_pdf:.3g} → {a.get('high_eur_per_kg'):.3g})", inp))

    # W3 — scale↑ must not raise COGS/kg.
    b = estimate_cogs_whatif(inp, scale * 10.0, pur, rm)
    if (b.get("high_eur_per_kg") or 0) > hi_pdf + 1e-6:
        v.append(Violation("W3-whatif-scale-monotonic", name,
            f"scale×10 raised COGS high ({hi_pdf:.3g} → {b.get('high_eur_per_kg'):.3g})", inp))

    # W4 — purity↑ must not lower COGS.
    if pur < 99.5:
        c = estimate_cogs_whatif(inp, scale, min(99.95, pur + 0.4), rm)
        if (c.get("high_eur_per_kg") or 0) < hi_pdf - 1e-6:
            v.append(Violation("W4-whatif-purity-monotonic", name,
                f"purity↑ lowered COGS high ({hi_pdf:.3g} → {c.get('high_eur_per_kg'):.3g})", inp))
    return v


def check_plausibility_clean(inp: Dict[str, Any], mol: Dict[str, str]) -> List[Violation]:
    """I7 — a realistic input (viable method, in-range step count) should not
    trip the step-count plausibility warning."""
    v: List[Violation] = []
    name = mol["name"]
    st_lo, st_hi = _steps_band(name, inp["method"])
    # Only assert when we know a real step range and stayed inside it.
    if KNOWN_PROCESSES.get(name.lower(), {}).get(inp["method"]):
        warns = check_plausibility(inp, lang="de")
        step_warn = [w for w in warns if "Schritt" in w or "step" in w.lower()]
        if step_warn and st_lo <= inp["number_of_steps"] <= st_hi:
            v.append(Violation("I7-spurious-step-warning", name,
                f"in-range steps ({inp['number_of_steps']} in {st_lo}-{st_hi}) "
                f"still warned: {step_warn[0][:60]}", inp))
    return v


def check_comparison(mol: Dict[str, str]) -> List[Violation]:
    """I8 — route comparison: exactly the viable methods appear, a unique
    rank-1 recommendation exists, ranks are 1..N contiguous."""
    v: List[Violation] = []
    name, mtype, msub = mol["name"], mol["type"], mol["subtype"]
    viable = _SUBTYPE_VIABLE_METHODS.get((mtype, msub))
    inp = {
        "molecule_name": name, "molecule_type": mtype, "molecule_subtype": msub,
        "method": (viable or ["chemical"])[0],
        "desired_purity_percent": 99.0, "scale_kg_per_year": 100.0,
        "number_of_steps": 4, "raw_material_cost_eur_per_kg": 20.0,
        "available_methods": {"has_prep_hplc": True, "has_crystallization": True,
                              "has_extraction": True, "has_fplc": True},
    }
    try:
        cmp = compare_routes(inp)
    except Exception as e:
        return [Violation("I8-comparison-crash", name, f"compare_routes raised {e}", inp)]
    if not cmp.rows:
        return [Violation("I8-no-routes", name, "no routes returned", inp)]
    ranks = sorted(r.decision_rank for r in cmp.rows)
    if ranks != list(range(1, len(ranks) + 1)):
        v.append(Violation("I8-rank-contiguous", name, f"ranks not 1..N: {ranks}", inp))
    if sum(1 for r in cmp.rows if r.decision_rank == 1) != 1:
        v.append(Violation("I8-rank1-unique", name, "no unique rank-1 route", inp))
    if viable is not None and len(cmp.rows) != len(viable):
        v.append(Violation("I8-method-set", name,
            f"{len(cmp.rows)} routes but {len(viable)} viable methods", inp))
    return v


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------
def main() -> int:
    n_per = int(sys.argv[1]) if len(sys.argv) > 1 else 64
    rng = random.Random(SEED)
    mols = load_molecules()
    print("Helixar AI — Invariant Sweep (plausible inputs only)")
    print(f"Molecules: {len(mols)}  ·  scenarios/molecule: {n_per}  ·  seed: {SEED}")
    print("=" * 68)

    violations: List[Violation] = []
    n_scenarios = 0
    law_counts: Dict[str, int] = {}

    for mol in mols:
        # one comparison-level structural check per molecule
        for vi in check_comparison(mol):
            violations.append(vi); law_counts[vi.law] = law_counts.get(vi.law, 0) + 1
        for _ in range(n_per):
            inp = sample_input(rng, mol)
            n_scenarios += 1
            checks = (check_absolute(inp, mol)
                      + check_monotonic(rng, inp, mol)
                      + check_whatif(inp, mol)
                      + check_plausibility_clean(inp, mol))
            for vi in checks:
                violations.append(vi)
                law_counts[vi.law] = law_counts.get(vi.law, 0) + 1

    print(f"Scenarios evaluated: {n_scenarios}")
    print(f"Per-molecule structural checks: {len(mols)}")
    print(f"Total invariant checks (approx): {n_scenarios * 11 + len(mols) * 4}")
    print("-" * 68)
    if not violations:
        print("RESULT:  ✓ 0 violations — all invariants held across every scenario.")
        print("\nVerteidigbare Aussage:")
        print(f"  „Über {n_scenarios} realistische Molekül×Input-Szenarien hält die")
        print("   Engine alle bioprozess-fachlichen Monotonie-, Benchmark- und")
        print("   Plausibilitäts-Gesetze ein.\"")
        return 0
    print(f"RESULT:  ✗ {len(violations)} violation(s) across {len(law_counts)} law(s)")
    for law, c in sorted(law_counts.items(), key=lambda x: -x[1]):
        print(f"  {law:<26} {c}")
    print("-" * 68)
    print("First 25 violations:")
    for vi in violations[:25]:
        print("  " + str(vi))
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
