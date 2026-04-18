# ...existing code...
from typing import Dict, Any, List


def explain_blueprint(blueprint: Dict[str, Any]) -> str:
    """Turn a production blueprint dict into a concise, professional explanation for founders.

    Tone: scientific, neutral, decision-support. Avoids lab instructions and guarantees.
    """
    rec = blueprint.get("recommended", {})
    op = blueprint.get("operating_parameters", {})
    exp = blueprint.get("expected_yield", {})
    risks: List[str] = blueprint.get("risk_notes", []) or []
    alts: List[Dict[str, Any]] = blueprint.get("alternatives", []) or []
    meta = blueprint.get("metadata", {})

    lines: List[str] = []

    # Executive summary
    name = rec.get("microorganism") or "Unknown microorganism"
    strain = rec.get("strain") or "unspecified strain"
    compound = rec.get("compound_name") or rec.get("compound_class") or "target compound"
    score = meta.get("decision_score")
    lines.append(f"Executive summary:")
    if score is not None:
        lines.append(f"- Recommended: {name} ({strain}) for production of {compound}.")
        lines.append(f"- Recommendation rationale: ranked by a composite decision score ({score:.3f}) that balances reported confidence and expected yield.")
    else:
        lines.append(f"- Recommended: {name} ({strain}) for production of {compound}.")
        lines.append(f"- Recommendation rationale: candidate selected based on available yield and confidence metadata.")

    # Input context (project requirements) — if present, mention them briefly
    input_params = blueprint.get("input_parameters") or blueprint.get("metadata", {}).get("input_parameters")
    if input_params:
        tm = input_params.get("target_molecule") or "(nicht angegeben)"
        sc = input_params.get("scale") or "(nicht angegeben)"
        sust = input_params.get("sustainability") or "(nicht angegeben)"
        lines.append("")
        lines.append("Project requirements:")
        lines.append(f"- Zielmolekül: {tm}")
        lines.append(f"- Produktionsmaßstab: {sc}")
        lines.append(f"- Nachhaltigkeit: {sust}")

    # Expected performance
    reported = exp.get("reported_yield")
    conservative = exp.get("conservative_yield_estimate")
    conf = exp.get("confidence_score", 0.0)
    perf_lines = []
    if reported is not None:
        perf_lines.append(f"reported yield ≈ {reported}")
    else:
        perf_lines.append("no quantitative reported yield")
    if conservative is not None:
        perf_lines.append(f"conservative estimate ≈ {conservative}")
    perf_lines.append(f"confidence score = {conf:.2f}")
    lines.append("")
    lines.append("Expected performance:")
    lines.append("- " + "; ".join(perf_lines) + ".")

    # Operating parameters (descriptive)
    lines.append("")
    lines.append("Available operating parameters (descriptive):")
    tr = op.get("temperature_range_recommended")
    pr = op.get("ph_range_recommended")
    ox = op.get("oxygen_level_recommended")
    if tr:
        lines.append(f"- Temperature (recommended range): {tr}")
    else:
        lines.append(f"- Temperature: not available")
    if pr:
        lines.append(f"- pH (recommended range): {pr}")
    else:
        lines.append(f"- pH: not available")
    if ox:
        lines.append(f"- Oxygen level (qualitative): {ox}")
    else:
        lines.append(f"- Oxygen level: not available")

    # Why recommended (decision logic summary)
    lines.append("")
    lines.append("Why this option was recommended:")
    lines.append("- The decision engine combines reported confidence and quantitative yield to produce a composite score; higher scores indicate better trade-offs between expected productivity and evidence quality.")
    lines.append("- A conservative yield estimate was computed to reflect uncertainty: this down-weights optimistic reports when confidence is low.")
    lines.append("- When operating parameter data are present, the engine biases toward entries with clearer provenance and measurable performance.")

    # Risks and caveats
    if risks:
        lines.append("")
        lines.append("Key risks and caveats:")
        for r in risks:
            lines.append(f"- {r}")
    else:
        lines.append("")
        lines.append("Key risks and caveats:")
        lines.append("- No specific risk notes available; verify data provenance and perform due diligence.")

    # Alternatives
    if alts:
        lines.append("")
        lines.append(f"Alternative strategies ({len(alts)}):")
        for a in alts:
            a_name = a.get("microorganism") or "Unknown"
            a_strain = a.get("strain") or "unspecified"
            a_yield = a.get("expected_yield")
            a_conf = a.get("confidence_score", 0.0)
            summary = f"{a_name} ({a_strain}) — confidence={a_conf:.2f}"
            if a_yield is not None:
                summary += f", reported_yield≈{a_yield}"
            lines.append(f"- {summary}")
    else:
        lines.append("")
        lines.append("Alternatives:")
        lines.append("- No lower-confidence alternatives identified in the database.")

    # Suggested non-operational next steps
    lines.append("")
    lines.append("Suggested next steps (decision-support, non-operational):")
    lines.append("- Review primary sources cited under metadata/source_reference and perform literature validation.")
    lines.append("- Conduct economic and regulatory assessment for the target process and market.")
    lines.append("- Prioritize analytical verification and targeted data collection to reduce key uncertainties before committing significant capital.")

    # Provenance
    src = meta.get("source_reference")
    if src:
        lines.append("")
        lines.append(f"Data provenance: {src}")

    return "\n".join(lines)
# ...existing code...