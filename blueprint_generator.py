# ...existing code...
from typing import Dict, Any, Optional, List, Tuple
from database import MicrobialDB


def _adjust_range(r: Optional[Tuple[Optional[float], Optional[float]]], margin: float = 0.1) -> Optional[Tuple[Optional[float], Optional[float]]]:
    """Produce a conservative operating range by shrinking the reported range toward the center.
    If only one bound present (None for min/max) try to return a +/- margin around the known bound.
    margin is fraction of the span to shrink (0.1 => keep central 90%).
    """
    # If no range provided, nothing to do
    if r is None:
        return None

    # If r is not a tuple/list of two numeric-or-None values (e.g. descriptive string like "moderat"),
    # we cannot safely compute a numeric adjustment. Return None to keep blueprint non-operational and
    # let the explanation layer note missing numeric ranges.
    if not (isinstance(r, (list, tuple)) and len(r) == 2):
        return None

    lo, hi = r
    # if both None, nothing to do
    if lo is None and hi is None:
        return None

    # if one side missing, create a small window around the known value
    if lo is None and hi is not None:
        try:
            hi_f = float(hi)
        except Exception:
            return None
        span = max(1.0, abs(hi_f))  # fallback span
        new_lo = hi_f - span * (1 - margin)
        return (round(new_lo, 3), round(hi_f, 3))
    if hi is None and lo is not None:
        try:
            lo_f = float(lo)
        except Exception:
            return None
        span = max(1.0, abs(lo_f))
        new_hi = lo_f + span * (1 - margin)
        return (round(lo_f, 3), round(new_hi, 3))

    # both present: shrink toward center
    try:
        lo_f = float(lo)
        hi_f = float(hi)
    except Exception:
        return None
    if lo_f > hi_f:
        lo_f, hi_f = hi_f, lo_f
    span = hi_f - lo_f
    shrink = span * margin
    return (round(lo_f + shrink / 2, 3), round(hi_f - shrink / 2, 3))


def _conservative_yield(expected: Optional[float], confidence: float) -> Optional[float]:
    """Estimate a conservative expected yield factoring in confidence.
    Lower confidence reduces the conservative yield proportionally.
    """
    if expected is None:
        return None
    # base conservative factor: between 0.6 (low) and 0.95 (high)
    factor = 0.6 + 0.35 * confidence
    return round(float(expected) * factor, 6)


def generate_production_blueprint(
    selected_strategy: Dict[str, Any],
    db: Optional[MicrobialDB] = None,
    alternatives_count: int = 2,
    safety_margin: float = 0.15,
) -> Dict[str, Any]:
    """Generate a structured production blueprint for a chosen microbial strategy.

    Inputs:
      - selected_strategy: dict (as returned by DecisionEngine.recommend) describing the chosen row.
      - db: optional MicrobialDB to locate alternative strategies (same compound_class).
      - alternatives_count: how many lower-confidence alternatives to include.
      - safety_margin: fraction used to tighten operating ranges (0..1).

    Returns a dictionary with:
      - recommended: microorganism/strain summary
      - operating_parameters: temperature/ph/oxygen ranges (conservative)
      - expected_yield: reported and conservative estimate
      - risk_notes: short list of qualitative risk points
      - alternatives: list of alternative strategies with lower confidence
      - metadata: original decision_score, source_reference, compound_class
    """
    rec = {
        "microorganism": selected_strategy.get("microorganism"),
        "strain": selected_strategy.get("strain"),
        "compound_class": selected_strategy.get("compound_class"),
        "compound_name": selected_strategy.get("compound_name"),
    }

    temp_range = selected_strategy.get("temperature_range")
    ph_range = selected_strategy.get("ph_range")
    oxygen = selected_strategy.get("oxygen_level")

    operating = {
        "temperature_range_recommended": _adjust_range(temp_range, margin=safety_margin),
        "ph_range_recommended": _adjust_range(ph_range, margin=safety_margin),
        "oxygen_level_recommended": oxygen,  # no transformation; user must interpret textual levels
    }

    reported_yield = selected_strategy.get("expected_yield")
    confidence = float(selected_strategy.get("confidence_score", 0.0)) if selected_strategy.get("confidence_score") is not None else 0.0
    conservative_yield = _conservative_yield(reported_yield, confidence)

    expected = {
        "reported_yield": (None if reported_yield is None else float(reported_yield)),
        "conservative_yield_estimate": conservative_yield,
        "confidence_score": confidence,
    }

    # Risk notes: concise, actionable
    risk_notes: List[str] = []
    if confidence < 0.4:
        risk_notes.append("Low confidence in reported performance — expect high variability and additional validation.")
    elif confidence < 0.7:
        risk_notes.append("Moderate confidence — plan small-scale validation before scale-up.")
    else:
        risk_notes.append("High confidence — suitable for pilot studies with standard QA/QC.")

    if reported_yield is None:
        risk_notes.append("No quantitative yield reported in DB — treat yield projections as uncertain.")
    else:
        # if conservative yield much lower than reported, warn
        if conservative_yield is not None and reported_yield > 0 and conservative_yield / reported_yield < 0.8:
            risk_notes.append("Conservative yield significantly lower than reported; investigate scale-up risks.")

    if temp_range is None and ph_range is None and oxygen is None:
        risk_notes.append("No operating parameters available; recommend literature review and lab data collection.")

    # Alternatives: use DB if available
    alternatives: List[Dict[str, Any]] = []
    if db is not None and rec["compound_class"] is not None:
        try:
            candidates = db.query_by_compound_class(rec["compound_class"])
            # exclude the selected strain + microorganism if possible
            mask_exclude = (candidates["microorganism"].astype(str) == str(rec["microorganism"])) & (candidates["strain"].astype(str) == str(rec["strain"]))
            filtered = candidates[~mask_exclude].copy()
            # require lower or equal confidence than selected, prefer lower confidence alternatives
            filtered["confidence_score"] = filtered["confidence_score"].astype(float)
            filtered = filtered.sort_values(by="confidence_score", ascending=False)
            # pick alternatives with lower confidence first (but still informative)
            lower = filtered[filtered["confidence_score"] < confidence]
            if lower.empty:
                # fallback: take next-best distinct candidates even if confidence >= selected
                lower = filtered
            for _, row in lower.head(alternatives_count).iterrows():
                alt = {
                    "microorganism": row.get("microorganism"),
                    "strain": row.get("strain"),
                    "expected_yield": row.get("expected_yield"),
                    "confidence_score": float(row.get("confidence_score", 0.0)),
                    "temperature_range": row.get("temperature_range"),
                    "ph_range": row.get("ph_range"),
                    "oxygen_level": row.get("oxygen_level"),
                    "source_reference": row.get("source_reference"),
                }
                alternatives.append(alt)
        except Exception:
            # If DB access fails, leave alternatives empty but note it
            risk_notes.append("Failed to retrieve alternatives from DB.")

    blueprint = {
        "recommended": rec,
        "operating_parameters": operating,
        "expected_yield": expected,
        "risk_notes": risk_notes,
        "alternatives": alternatives,
        "metadata": {
            "decision_score": selected_strategy.get("decision_score"),
            "source_reference": selected_strategy.get("source_reference"),
        },
    }

    # If the selected strategy carries an analysis (process-optimization mode), include it
    try:
        if isinstance(selected_strategy, dict) and selected_strategy.get("analysis"):
            blueprint["analysis"] = selected_strategy.get("analysis")
    except Exception:
        pass

    # Include the original input context (constraints) if present on the selected strategy.
    try:
        input_ctx = selected_strategy.get("input_context") if isinstance(selected_strategy, dict) else None
        if input_ctx:
            blueprint["input_parameters"] = {
                "target_molecule": input_ctx.get("target_molecule"),
                "application": input_ctx.get("application"),
                "scale": input_ctx.get("scale"),
                "purity": input_ctx.get("purity"),
                "sustainability": input_ctx.get("sustainability"),
                "preferred_method": input_ctx.get("preferred_method"),
                "infrastructure": input_ctx.get("infrastructure"),
            }
    except Exception:
        # Don't fail blueprint creation for missing/invalid input_context
        pass

    return blueprint
# ...existing code...