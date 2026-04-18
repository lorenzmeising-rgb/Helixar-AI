# ...existing code...
from typing import Dict, Any, Optional, List, Tuple
from datetime import datetime

try:
    from reportlab.lib.pagesizes import A4
    from reportlab.lib import colors
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.lib.units import mm
    from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak
except Exception as e:
    raise ImportError("reportlab is required for PDF export. Install with `pip install reportlab`.") from e


def _fmt_range(r: Optional[Tuple[Optional[float], Optional[float]]]) -> str:
    if r is None:
        return "n/a"
    lo, hi = r
    if lo is None and hi is None:
        return "n/a"
    return f"{'' if lo is None else lo} — {'' if hi is None else hi}"


def export_blueprint_pdf(
    blueprint: Dict[str, Any],
    explanation: str,
    path: str,
    title: str = "Production Blueprint",
    author: Optional[str] = None,
) -> str:
    """
    Export a production blueprint and its textual explanation into a clean PDF.

    Parameters
    ----------
    blueprint: dict
        Structured blueprint produced by generate_production_blueprint(...)
    explanation: str
        Plain-text explanation (e.g. from explain_blueprint(...))
    path: str
        Output PDF file path
    title: str
        Document title shown on top
    author: Optional[str]
        Optional author / organization line

    Returns
    -------
    str
        The path written to (same as input path)
    """
    doc = SimpleDocTemplate(path, pagesize=A4, leftMargin=18 * mm, rightMargin=18 * mm, topMargin=18 * mm, bottomMargin=18 * mm)
    styles = getSampleStyleSheet()
    # Tweak/add styles for clarity
    styles.add(ParagraphStyle(name="Heading1Center", parent=styles["Heading1"], alignment=1))
    styles.add(ParagraphStyle(name="Mono", parent=styles["Code"], fontSize=9))
    normal = styles["Normal"]
    heading = styles["Heading2"]

    elems: List[Any] = []

    # Header
    elems.append(Paragraph("Helixar", styles["Heading1"]))
    elems.append(Paragraph(title, styles["Heading1Center"]))
    meta_line = f"Generated: {datetime.utcnow().strftime('%Y-%m-%d %H:%M UTC')}"
    if author:
        meta_line = f"{author} — {meta_line}"
    elems.append(Paragraph(meta_line, normal))
    elems.append(Spacer(1, 6))

    # Recommended summary
    rec = blueprint.get("recommended", {})
    elems.append(Paragraph("Recommended strategy", heading))
    rec_table_data = [
        ["Microorganism", rec.get("microorganism", "n/a")],
        ["Strain", rec.get("strain", "n/a")],
        ["Compound class", rec.get("compound_class", "n/a")],
        ["Compound name", rec.get("compound_name", "n/a")],
    ]
    t = Table(rec_table_data, colWidths=(60 * mm, None))
    t.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, 0), colors.whitesmoke),
                           ("VALIGN", (0, 0), (-1, -1), "TOP"),
                           ("INNERGRID", (0, 0), (-1, -1), 0.25, colors.lightgrey),
                           ("BOX", (0, 0), (-1, -1), 0.25, colors.lightgrey)]))
    elems.append(t)
    elems.append(Spacer(1, 8))

    # Operating parameters
    elems.append(Paragraph("Operating parameters (recommended, conservative)", heading))
    op = blueprint.get("operating_parameters", {})
    op_table_data = [
        ["Temperature (recommended)", _fmt_range(op.get("temperature_range_recommended"))],
        ["pH (recommended)", _fmt_range(op.get("ph_range_recommended"))],
        ["Oxygen level (qualitative)", op.get("oxygen_level_recommended", "n/a") or "n/a"],
    ]
    t2 = Table(op_table_data, colWidths=(60 * mm, None))
    t2.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "TOP"),
                            ("INNERGRID", (0, 0), (-1, -1), 0.25, colors.lightgrey),
                            ("BOX", (0, 0), (-1, -1), 0.25, colors.lightgrey)]))
    elems.append(t2)
    elems.append(Spacer(1, 8))

    # Expected yield & confidence
    elems.append(Paragraph("Expected yield & confidence", heading))
    ey = blueprint.get("expected_yield", {})
    ey_table_data = [
        ["Reported yield", ey.get("reported_yield", "n/a")],
        ["Conservative estimate", ey.get("conservative_yield_estimate", "n/a")],
        ["Confidence score", f"{ey.get('confidence_score', 'n/a')}"],
    ]
    t3 = Table(ey_table_data, colWidths=(60 * mm, None))
    t3.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "TOP"),
                            ("INNERGRID", (0, 0), (-1, -1), 0.25, colors.lightgrey),
                            ("BOX", (0, 0), (-1, -1), 0.25, colors.lightgrey)]))
    elems.append(t3)
    elems.append(Spacer(1, 8))

    # Risk notes
    elems.append(Paragraph("Key risks and caveats", heading))
    risks: List[str] = blueprint.get("risk_notes", []) or []
    if risks:
        for r in risks:
            elems.append(Paragraph(f"• {r}", normal))
    else:
        elems.append(Paragraph("No specific risk notes recorded.", normal))
    elems.append(Spacer(1, 10))

    # Alternatives
    elems.append(Paragraph("Alternative strategies", heading))
    alts = blueprint.get("alternatives", []) or []
    if alts:
        alt_table_data = [["Microorganism", "Strain", "Reported yield", "Confidence", "Source"]]
        for a in alts:
            alt_table_data.append([
                a.get("microorganism", "n/a"),
                a.get("strain", "n/a"),
                str(a.get("expected_yield", "n/a")),
                f"{a.get('confidence_score', 'n/a')}",
                a.get("source_reference", "") or "",
            ])
        alt_tbl = Table(alt_table_data, colWidths=(45 * mm, 35 * mm, 30 * mm, 20 * mm, None))
        alt_tbl.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#F2F4F7")),
            ("TEXTCOLOR", (0, 0), (-1, 0), colors.black),
            ("ALIGN", (2, 1), (3, -1), "CENTER"),
            ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
            ("INNERGRID", (0, 0), (-1, -1), 0.25, colors.lightgrey),
            ("BOX", (0, 0), (-1, -1), 0.25, colors.lightgrey),
        ]))
        elems.append(alt_tbl)
    else:
        elems.append(Paragraph("No alternatives identified.", normal))

    elems.append(PageBreak())

    # Full explanation text (for founders)
    elems.append(Paragraph("Explanation (decision-support)", styles["Heading2"]))
    # Explanation may be long; use Paragraphs split by double newline
    for block in explanation.split("\n\n"):
        elems.append(Paragraph(block.replace("\n", "<br/>"), normal))
        elems.append(Spacer(1, 6))

    elems.append(Spacer(1, 12))

    # Footer / provenance
    meta = blueprint.get("metadata", {})
    src = meta.get("source_reference")
    decision_score = meta.get("decision_score")
    footer_lines = []
    if decision_score is not None:
        footer_lines.append(f"Decision score: {decision_score}")
    if src:
        footer_lines.append(f"Source: {src}")
    if footer_lines:
        elems.append(Spacer(1, 12))
        elems.append(Paragraph("Provenance & metadata", styles["Heading2"]))
        for fl in footer_lines:
            elems.append(Paragraph(fl, normal))

    # Build PDF
    doc.build(elems)

    return path
# ...existing code...