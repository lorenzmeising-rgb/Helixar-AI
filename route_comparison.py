"""
route_comparison.py
===================
Feature A — Side-by-Side Routenvergleich.

For every molecule, the engine can score all three viable production
methods in parallel — chemical synthesis, biotechnological production,
and extraction. Process Scientists need to see those side-by-side to
answer the question "which route should I pick?", which is the core
job-to-be-done for a CMC decision support tool.

Public API:
    compare_routes(process_input) -> RouteComparisonResult

The result carries one ComparisonRow per method, plus an engine
recommendation explaining which row would be chosen and why.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional, Tuple

# We import the COGS estimator at module level (cheap) and the
# DecisionEngine / make_engine lazily inside compare_routes so this
# module doesn't slow down report imports when route comparison is
# not requested.


METHODS_DE = {
    "chemical": "Chemische Synthese",
    "biotechnological": "Biotechnologische Produktion",
    "extraction": "Natürliche Extraktion",
}

METHODS_EN = {
    "chemical": "Chemical synthesis",
    "biotechnological": "Biotechnological production",
    "extraction": "Natural extraction",
}

# Per-subtype which methods are conceptually applicable. Some methods
# don't make physical sense for some subtypes (you don't "extract" a
# de novo synthetic SPPS-peptide, you don't run "chemical synthesis"
# of a full mAb in industrial practice).
_SUBTYPE_VIABLE_METHODS: Dict[Tuple[str, str], List[str]] = {
    ("small_molecule", "non_volatile"): ["chemical", "biotechnological"],
    ("small_molecule", "volatile"):     ["biotechnological", "chemical"],
    ("natural_product", "alkaloid"):    ["chemical", "extraction", "biotechnological"],
    ("natural_product", "terpene"):     ["extraction", "biotechnological", "chemical"],
    ("peptide", "linear"):              ["chemical", "biotechnological"],
    ("peptide", "cyclic"):              ["biotechnological", "chemical"],
    ("protein", "antibody"):            ["biotechnological"],   # CHO is the only realistic route
    ("protein", "enzyme"):              ["biotechnological"],
}


@dataclass
class ComparisonRow:
    """One row of the side-by-side comparison."""
    method: str
    method_label_de: str
    method_label_en: str
    cost_label: str
    risk_label: str
    efficiency_label: str
    toxicity_label: str
    cost_score: int
    risk_score: int
    efficiency_score: int
    cogs_low_eur_per_kg: Optional[float]
    cogs_high_eur_per_kg: Optional[float]
    cogs_anchor_de: Optional[str]
    cogs_anchor_en: Optional[str]
    cogs_confidence: Optional[str]
    has_concrete_steps: bool
    n_production_steps: int
    expected_yield: Optional[str]
    decision_rank: int = 0   # 1 = engine's preferred, 2 = second, …
    decision_reason_de: Optional[str] = None
    decision_reason_en: Optional[str] = None


@dataclass
class RouteComparisonResult:
    """Result of the side-by-side comparison."""
    molecule_name: str
    molecule_type: str
    molecule_subtype: str
    rows: List[ComparisonRow] = field(default_factory=list)
    recommended_method: Optional[str] = None
    recommendation_reason_de: Optional[str] = None
    recommendation_reason_en: Optional[str] = None


def _make_engine():
    """Lazy engine creation — same fallback chain as full_audit.py."""
    from decision_engine import DecisionEngine
    try:
        from database import MicrobialDB
        return DecisionEngine(db=MicrobialDB())
    except Exception:
        import pandas as pd
        df = pd.DataFrame(columns=[
            "microorganism", "strain", "compound_class",
            "compound_name", "confidence_score",
        ])
        class _DummyDB:
            def __init__(self, d): self._df = d
            def all(self): return self._df
        return DecisionEngine(db=_DummyDB(df))


def _score_for_label(label: str, metric: str) -> int:
    """Mirror of report_generator._score_for(). Keeps the comparison
    score in sync with what the executive summary shows."""
    key = str(label).lower()
    if metric == "efficiency":
        return {"low": 3, "niedrig": 3, "medium": 6, "mittel": 6,
                "high": 8, "hoch": 8, "very high": 10, "sehr hoch": 10}.get(key, 5)
    return {"low": 2, "niedrig": 2, "medium": 5, "mittel": 5,
            "high": 7, "hoch": 7, "very high": 9, "sehr hoch": 9}.get(key, 5)


def _check_concrete_hint(mname: str, method: str) -> Tuple[bool, int, Optional[str]]:
    """Return (has_steps, n_steps, yield_string) for a given molecule/method."""
    try:
        from concrete_recommendations import MOLECULE_HINTS
    except Exception:
        return False, 0, None
    hints = MOLECULE_HINTS.get(str(mname or "").lower(), {})
    hint = hints.get(method)
    if not hint:
        return False, 0, None
    steps = (hint.get("production") or {}).get("de") or hint.get("production") or []
    if isinstance(steps, dict):
        steps = steps.get("de") or steps.get("en") or []
    n = len(steps) if isinstance(steps, list) else 0
    # Extract yield range as a readable string. Chemical and extraction
    # routes store a percentage yield (yield_range_percent /
    # yield_range_percent_w_w); biotechnological routes store a titer in
    # g/L (yield_range_g_per_l). The latter was previously ignored, which
    # left the "Yield" column blank ("—") for every fermentation route
    # (50 of the 79 molecules).
    y_pct = hint.get("yield_range_percent") or hint.get("yield_range_percent_w_w")
    y_gl = hint.get("yield_range_g_per_l")
    yield_str = None
    if isinstance(y_pct, (list, tuple)) and len(y_pct) == 2:
        yield_str = f"{y_pct[0]}–{y_pct[1]} %"
    elif isinstance(y_gl, (list, tuple)) and len(y_gl) == 2:
        yield_str = f"{y_gl[0]}–{y_gl[1]} g/L"
    return n > 0, n, yield_str


def _rank_rows(rows: List[ComparisonRow], process_input: Dict[str, Any]) -> List[ComparisonRow]:
    """Sort rows by a weighted score so the engine's preferred route
    appears first. Scoring rationale:
      - lower cost is better (weight 1.0)
      - lower risk is better (weight 0.8)
      - higher efficiency is better (weight 0.7)
      - bonus for having concrete steps in MOLECULE_HINTS (+1 — we
        can give the user real reagents and conditions)
      - bonus when method matches process_input.method (user already
        considered it viable, slight default preference)
    """
    user_method = (process_input.get("method") or "").lower()

    def score(r: ComparisonRow) -> float:
        s = 10.0  # baseline
        s -= r.cost_score * 1.0
        s -= r.risk_score * 0.8
        s += r.efficiency_score * 0.7
        if r.has_concrete_steps:
            s += 1.0
        if r.method == user_method:
            s += 0.3
        return s

    ranked = sorted(rows, key=score, reverse=True)
    for i, r in enumerate(ranked, start=1):
        r.decision_rank = i
    return ranked


def _build_reason(top: ComparisonRow, lang: str) -> str:
    """Human-readable explanation of why `top` is the preferred row."""
    bits_de: List[str] = []
    bits_en: List[str] = []
    if top.cost_score <= 3:
        bits_de.append("niedrige Stückkosten")
        bits_en.append("low unit cost")
    elif top.cost_score >= 8:
        bits_de.append("trotz hoher Stückkosten")
        bits_en.append("despite high unit cost")
    if top.risk_score <= 3:
        bits_de.append("etablierter Prozess mit geringem Operations-Risiko")
        bits_en.append("established process with low operational risk")
    elif top.risk_score <= 5:
        bits_de.append("akzeptables Operations-Risiko")
        bits_en.append("acceptable operational risk")
    if top.efficiency_score >= 8:
        bits_de.append("robuste Standardroute")
        bits_en.append("robust standard route")
    if top.has_concrete_steps:
        bits_de.append("konkrete industrielle Vorlagen verfügbar")
        bits_en.append("concrete industrial precedents available")

    if not bits_de:
        bits_de = ["beste Gesamtbilanz aus Cost / Risk / Efficiency"]
        bits_en = ["best overall Cost / Risk / Efficiency balance"]

    if lang == "en":
        return "; ".join(bits_en)
    return "; ".join(bits_de)


def compare_routes(process_input: Dict[str, Any]) -> RouteComparisonResult:
    """Score all viable methods for a molecule and return them as a
    ranked side-by-side comparison.
    """
    from cogs_estimator import estimate_cogs

    mtype = str(process_input.get("molecule_type") or "").lower()
    msub = str(process_input.get("molecule_subtype") or "").lower()
    mname = str(process_input.get("molecule_name") or "")

    viable = _SUBTYPE_VIABLE_METHODS.get((mtype, msub), ["chemical", "biotechnological"])

    engine = _make_engine()
    rows: List[ComparisonRow] = []

    for method in viable:
        # Run a counterfactual: same inputs, override method
        pi_alt = dict(process_input)
        pi_alt["method"] = method
        try:
            recs = engine.analyze_process(pi_alt) or []
        except Exception:
            recs = []
        a = recs[0].get("analysis", {}) if recs else {}
        cogs = {}
        try:
            cogs = estimate_cogs(pi_alt) or {}
        except Exception:
            cogs = {}

        cost_lbl = (a.get("cost") or "medium").lower()
        risk_lbl = (a.get("risk") or "medium").lower()
        eff_lbl = (a.get("efficiency") or "medium").lower()
        tox_lbl = (a.get("toxicity") or a.get("final_toxicity") or "medium").lower()
        has_steps, n_steps, yield_str = _check_concrete_hint(mname, method)

        rows.append(ComparisonRow(
            method=method,
            method_label_de=METHODS_DE.get(method, method),
            method_label_en=METHODS_EN.get(method, method),
            cost_label=cost_lbl,
            risk_label=risk_lbl,
            efficiency_label=eff_lbl,
            toxicity_label=tox_lbl,
            cost_score=_score_for_label(cost_lbl, "cost"),
            risk_score=_score_for_label(risk_lbl, "risk"),
            efficiency_score=_score_for_label(eff_lbl, "efficiency"),
            cogs_low_eur_per_kg=cogs.get("low_eur_per_kg"),
            cogs_high_eur_per_kg=cogs.get("high_eur_per_kg"),
            cogs_anchor_de=cogs.get("anchor_de"),
            cogs_anchor_en=cogs.get("anchor_en"),
            cogs_confidence=cogs.get("confidence"),
            has_concrete_steps=has_steps,
            n_production_steps=n_steps,
            expected_yield=yield_str,
        ))

    ranked = _rank_rows(rows, process_input)
    top = ranked[0] if ranked else None
    reason_de = _build_reason(top, "de") if top else None
    reason_en = _build_reason(top, "en") if top else None

    # Add reason on every row (for the table)
    for r in ranked:
        r.decision_reason_de = (
            f"{r.method_label_de} — Rang {r.decision_rank}" if r != top
            else f"Empfohlen ({reason_de})"
        )
        r.decision_reason_en = (
            f"{r.method_label_en} — rank {r.decision_rank}" if r != top
            else f"Recommended ({reason_en})"
        )

    return RouteComparisonResult(
        molecule_name=mname,
        molecule_type=mtype,
        molecule_subtype=msub,
        rows=ranked,
        recommended_method=top.method if top else None,
        recommendation_reason_de=reason_de,
        recommendation_reason_en=reason_en,
    )
