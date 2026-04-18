from typing import Dict, Any, List, Optional

from language_de import LABELS


from molecules_db import get_smiles_for as _get_smiles_for


def _safe_val(val: Optional[Any]) -> str:
    """Safe display for values in German.

    - None -> 'nicht verfügbar'
    - list/tuple -> join with a readable separator
    - otherwise -> str(value)
    """
    if val is None:
        return "nicht verfügbar"
    if isinstance(val, (list, tuple)):
        return " – ".join(["n/a" if v is None else str(v) for v in val])
    return str(val)


def _confidence_label(conf: Optional[Any]) -> str:
    """Map a numeric confidence to a German verbal label."""
    if conf is None:
        return "nicht verfügbar"
    try:
        c = float(conf)
    except Exception:
        return "nicht verfügbar"
    if c >= 0.75:
        return "hoch"
    if c >= 0.4:
        return "mittel"
    return "niedrig"


def generate_report(blueprint: Dict[str, Any]) -> str:
    """Generate a concise, German plain-text executive report from a blueprint.

    This is a lightweight, non-operational summary intended for quick reading
    or API responses. For printable PDF output use export_report_pdf().
    """
    rec = blueprint.get("recommended", {})
    meta = blueprint.get("metadata", {})
    ip = blueprint.get("input_parameters") or {}

    def _method_of(d: Dict[str, Any]) -> Optional[Any]:
        return d.get("method") or d.get("microorganism")

    def _process_of(d: Dict[str, Any]) -> Optional[Any]:
        return d.get("process_type") or d.get("strain")

    compound_name = _safe_val(
        rec.get("compound_name") or ip.get("target_molecule") or meta.get("compound_name") or "nicht angegeben"
    )
    method = _safe_val(_method_of(rec))
    process_type = _safe_val(_process_of(rec))

    # If an analysis is present (process-optimization mode), produce an analysis-led summary
    analysis = blueprint.get("analysis")
    if analysis:
        lines = [
            "Helixar AI",
            "Prozessanalyse Bericht",
            "",
            "Executive Summary:",
            f"Zielmolekül: {compound_name}",
            f"SMILES: {_get_smiles_for(rec.get('compound_name') or ip.get('target_molecule') or meta.get('compound_name'))}",
            f"Methode: {method}",
            f"Prozesstyp: {process_type}",
            "",
            f"Effizienz: {analysis.get('efficiency', 'n/a')}",
            f"Kosten (Level): {analysis.get('costLevel', analysis.get('cost', 'n/a'))}",
            f"Risiko: {analysis.get('risk', 'n/a')}",
            f"Toxizität: {analysis.get('final_toxicity', analysis.get('toxicity', 'n/a'))}",
            "",
            "Wesentliche Probleme:",
        ]
        for it in analysis.get("issues", []):
            lines.append(f"- {it}")
        lines.append("")
        lines.append("Vorgeschlagene Verbesserungen:")
        for imp in analysis.get("improvements", []):
            lines.append(f"- {imp}")
        # brief cost notes
        if analysis.get("costDrivers"):
            lines.append("")
            lines.append("Hauptkostentreiber:")
            for d in analysis.get("costDrivers", [])[:6]:
                lines.append(f"- {d}")
        return "\n".join(lines)

    # Legacy, short strategy summary when no analysis is present
    try:
        score_val = float(rec.get("decision_score") or meta.get("decision_score") or 0)
    except Exception:
        score_val = 0.0

    def _interpret(score: float) -> str:
        if score >= 0.75:
            return "Hohe Eignung für industrielle Umsetzung"
        if score >= 0.5:
            return "Gute Eignung bei weiterer Validierung"
        return "Begrenzte Eignung — weitere Validierung empfohlen"

    interpretation = _interpret(score_val)

    lines = [
        "Helixar AI",
        "Produktionsstrategie Bericht",
        "",
        f"Zielmolekül: {compound_name}",
        f"SMILES: {_get_smiles_for(rec.get('compound_name') or ip.get('target_molecule') or meta.get('compound_name'))}",
        f"Methode: {method}",
        f"Prozesstyp: {process_type}",
        "",
        f"Entscheidungs-Score: {round(score_val, 2):.2f}",
        interpretation,
        "",
        "Kernaussagen:",
        "- Überlegenswerte Produktionsroute",
        "- Weiterführende Validierung empfohlen",
    ]

    return "\n".join(lines)


# -----------------------------
# PDF export (reportlab-based)
# -----------------------------

def export_report_pdf(blueprint: Dict[str, Any], path: str, title: Optional[str] = None, author: Optional[str] = None) -> str:
    """Export a focused German PDF report from a blueprint.

    Produces a short, consultant-style PDF with header/footer and four pages:
    - Executive Summary
    - Process Analysis
    - Risks & Trade-offs
    - Recommended Actions

    Returns the path to the written file.
    """
    try:
        from reportlab.lib.pagesizes import A4
        from reportlab.lib.units import mm
        from reportlab.lib import colors
        from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, PageBreak
        from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    except Exception as exc:
        raise ImportError("reportlab is required for PDF export. Install it with: pip install reportlab") from exc

    styles = getSampleStyleSheet()
    styles.add(
        ParagraphStyle(
            name="ExecTitle",
            parent=styles["Heading1"],
            fontSize=18,
            leading=22,
            spaceAfter=12,
            alignment=1,
        )
    )
    styles.add(
        ParagraphStyle(
            name="HelixarTitle",
            parent=styles["Title"],
            fontSize=28,
            leading=32,
            alignment=1,
            spaceAfter=8,
        )
    )
    styles.add(
        ParagraphStyle(name="DocSubtitle", parent=styles["Heading2"], fontSize=12, leading=14, alignment=0)
    )
    styles.add(
        ParagraphStyle(name="Recommendation", parent=styles["Heading1"], fontSize=20, leading=24, alignment=0)
    )
    styles.add(ParagraphStyle(name="Score", parent=styles["Heading2"], fontSize=14, leading=18, alignment=0))
    styles.add(ParagraphStyle(name="Section", parent=styles["Heading2"], fontSize=13, leading=16, spaceAfter=6, spaceBefore=8))
    styles.add(ParagraphStyle(name="BodySmall", parent=styles["Normal"], fontSize=10, leading=14))
    # Monospace style for SMILES representation
    styles.add(ParagraphStyle(name="Mono", parent=styles["Normal"], fontName="Courier", fontSize=10, leading=12))

    normal = styles["BodySmall"]

    rec = blueprint.get("recommended", {})
    meta = blueprint.get("metadata", {})
    risks = blueprint.get("risk_notes", []) or []

    def _fmt_num(v: Optional[Any]) -> str:
        try:
            return f"{float(v):.2f}"
        except Exception:
            return _safe_val(v)

    METHOD_MAP = {
        "chemical": "Chemische Herstellung",
        "biotech": "Biotechnologische Herstellung",
        "extraction": "Natürliche Extraktion",
    }

    PROCESS_MAP = {
        "synthetic": "Synthese",
        "fermentation": "Fermentation",
        "extraction": "Extraktion",
    }

    def _translate_method_display(m: Optional[Any]) -> str:
        if m is None:
            return _safe_val(m)
        key = str(m).lower()
        if "biotech" in key or "biotechn" in key:
            return METHOD_MAP.get("biotech")
        if "chem" in key or "chemical" in key:
            return METHOD_MAP.get("chemical")
        if "extract" in key:
            return METHOD_MAP.get("extraction")
        return METHOD_MAP.get(key, _safe_val(m))

    def _translate_process_display(p: Optional[Any]) -> str:
        if p is None:
            return _safe_val(p)
        key = str(p).lower()
        if "synth" in key or "synt" in key:
            return PROCESS_MAP.get("synthetic")
        if "ferment" in key:
            return PROCESS_MAP.get("fermentation")
        if "extract" in key:
            return PROCESS_MAP.get("extraction")
        return PROCESS_MAP.get(key, _safe_val(p))

    def _explain_issue_pair(issue: str):
        s = str(issue or "").strip()
        low = s.lower()
        if "step" in low or "synth" in low:
            reason = "→ Jede zusätzliche Syntheseschritt erhöht Kosten, Durchlaufzeit und Komplexität, reduziert die Effizienz."
        elif "purif" in low or "reinig" in low or "purificat" in low:
            reason = "→ Hohe Reinigungsanforderungen treiben nachgelagerte Arbeit und Kosten; mindern Skalierbarkeit."
        elif "raw material" in low or "rohstoff" in low or "expensive" in low or "teuer" in low:
            reason = "→ Teure Rohstoffe sind ein wiederkehrender Kostentreiber und beeinflussen die Stückkosten stark."
        elif "bioreactor" in low or "bioreaktor" in low:
            reason = "→ Fehlende Bioreaktor-Infrastruktur schränkt biotechnologische Routen ein und führt zu kostspieligen Anpassungen."
        elif "waste" in low or "abfall" in low:
            reason = "→ Strenge Abfallvorschriften erhöhen Compliance-Last und Betriebskosten, besonders bei gefährlichen Abfällen."
        elif "stability" in low or "stabil" in low:
            reason = "→ Geringe Stabilität erhöht Ausfallrisiko beim Scale-up und treibt Entwicklungskosten."
        else:
            reason = "→ Dieses Thema wirkt sich negativ auf Kosten, Risiko oder Skalierbarkeit aus und erfordert Untersuchung."
        return s, reason

    def _explain_improvement(impr: str):
        s = str(impr or "").strip()
        low = s.lower()
        if "reduce number" in low or "reduce the number" in low or "reduce number of" in low:
            action = "Reduzieren der Anzahl an Syntheseschritten"
            impact = "→ Hohes Einsparpotenzial und bessere Skalierbarkeit"
            reason = "→ Weniger Schritte verringern Materialeinsatz, Prozesszeit und Nachreinigung."
        elif "crystall" in low or "crystalliz" in low or "purificat" in low:
            action = "Optimierung der Reinigung (z. B. Kristallisation)"
            impact = "→ Mittleres Einsparpotenzial bei Aufreinigungskosten"
            reason = "→ Alternative Methoden reduzieren Lösungsmittelverbrauch und erhöhen Durchsatz."
        elif "raw material" in low or "optimiz" in low or "sourcing" in low:
            action = "Optimierung der Rohstoffauswahl und Beschaffung"
            impact = "→ Mittleres Einsparpotenzial auf Stückkosten"
            reason = "→ Bessere Beschaffung reduziert wiederkehrende Materialkosten."
        elif "stabil" in low or "process stability" in low:
            action = "Verbesserung der Prozessstabilität"
            impact = "→ Mittleres Einsparpotenzial durch reduzierte Nacharbeit"
            reason = "→ Stabilere Prozesse unterstützen Reproduzierbarkeit und Scale-up."
        elif "toxic" in low or "waste" in low or "hazard" in low:
            action = "Reduzieren von gefährlichen Reagenzien und Abfallströmen"
            impact = "→ Mittleres Einsparpotenzial bei Compliance- und Entsorgungskosten"
            reason = "→ Weniger gefährliche Chemikalien senken Behandlungskosten und regulatorischen Aufwand."
        else:
            action = s or "Empfohlene Maßnahme"
            impact = "→ Erwarteter Effekt: Moderat"
            reason = "→ Diese Maßnahme adressiert identifizierte Schwächen und sollte Prozesskennzahlen verbessern."
        return action, impact, reason

    def draw_header_footer(canvas, doc):
        canvas.saveState()
        canvas.setFont("Helvetica-Bold", 11)
        canvas.setFillColor(colors.HexColor("#111111"))
        canvas.drawString(18 * mm, A4[1] - 18 * mm + 6, "Helixar AI")
        canvas.setFont("Helvetica", 8)
        canvas.setFillColor(colors.HexColor("#666666"))
        canvas.drawString(18 * mm, A4[1] - 18 * mm - 6, "KI-gestützte Entscheidungsunterstützung")
        canvas.setStrokeColor(colors.HexColor("#E6E9EE"))
        canvas.setLineWidth(0.5)
        canvas.line(18 * mm, A4[1] - 24 * mm, A4[0] - 18 * mm, A4[1] - 24 * mm)
        canvas.setFont("Helvetica", 8)
        canvas.setFillColor(colors.HexColor("#666666"))
        canvas.drawString(18 * mm, 12 * mm, "Nicht-operativer Entscheidungsreport – keine Laboranleitung")
        gen_date = __import__("datetime").datetime.utcnow().strftime('%Y-%m-%d %H:%M UTC')
        canvas.drawRightString(A4[0] - 18 * mm, 12 * mm, f"Erstellt: {gen_date}")
        canvas.restoreState()

    doc = SimpleDocTemplate(path, pagesize=A4, leftMargin=18 * mm, rightMargin=18 * mm, topMargin=36 * mm, bottomMargin=28 * mm)

    elems: List[Any] = []

    compound_title = _safe_val(rec.get('compound_name') or meta.get('compound_name') or (blueprint.get('input_parameters') or {}).get('target_molecule'))
    if title is None:
        title = f"Produktionsstrategie Bericht — {compound_title}"

    # PAGE 1 — Executive Summary
    elems.append(Paragraph("Helixar AI", styles["HelixarTitle"]))
    elems.append(Spacer(1, 6))
    elems.append(Paragraph("Prozessoptimierung — Executive Summary", styles["DocSubtitle"]))
    elems.append(Spacer(1, 10))

    method = rec.get('method') or rec.get('microorganism')
    ptype = rec.get('process_type') or rec.get('strain')

    current_date = __import__("datetime").datetime.utcnow().strftime('%Y-%m-%d')
    meta_table = [
        f"Zielmolekül: {compound_title}",
        f"Methode: {_translate_method_display(method)}",
        f"Prozesstyp: {_translate_process_display(ptype)}",
        f"Datum: {current_date}",
    ]
    for m in meta_table:
        elems.append(Paragraph(m, normal))
    elems.append(Spacer(1, 12))

    # SMILES: fetch from the small molecule database (case-insensitive). Show a monospace
    # SMILES line directly under the header for scientific credibility.
    raw_molecule = rec.get('compound_name') or meta.get('compound_name') or (blueprint.get('input_parameters') or {}).get('target_molecule')
    smiles = _get_smiles_for(raw_molecule)
    # spacing between name and SMILES
    elems.append(Spacer(1, 6))
    elems.append(Paragraph(f"SMILES: {smiles}", styles["Mono"]))
    elems.append(Spacer(1, 8))

    analysis = blueprint.get("analysis") or {}

    elems.append(Paragraph("Executive Summary", styles["Recommendation"]))
    elems.append(Spacer(1, 8))

    def perf_icon(v: str) -> str:
        key = str(v).lower()
        if key in ("high", "very high", "hoch", "sehr hoch"):
            return "🔴"
        if key in ("medium", "mittel"):
            return "⚠️"
        return "🟢"

    eff = analysis.get('efficiency', 'n/a')
    cost_val = analysis.get('cost', 'n/a')
    risk_val = analysis.get('risk', 'n/a')
    tox = analysis.get('final_toxicity', analysis.get('toxicity', 'n/a'))

    elems.append(Paragraph(f"<b>Efficiency:</b> {perf_icon(eff)} {str(eff).upper()}", styles["Score"]))
    elems.append(Paragraph(f"<b>Cost:</b> {perf_icon(cost_val)} {str(cost_val).upper()}", normal))
    elems.append(Paragraph(f"<b>Risk:</b> {perf_icon(risk_val)} {str(risk_val).upper()}", normal))
    elems.append(Paragraph(f"<b>Toxicity:</b> {perf_icon(tox)} {str(tox).upper()}", normal))
    elems.append(Spacer(1, 12))

    try:
        eff_score = float(analysis.get('efficiency_score', 0.5))
        cost_score = float(analysis.get('cost_score', 0.5))
        conf_metric = (eff_score + (1.0 - cost_score)) / 2.0
        if conf_metric >= 0.66:
            conf_label = "High"
        elif conf_metric >= 0.33:
            conf_label = "Medium"
        else:
            conf_label = "Low"
    except Exception:
        conf_label = "Medium"
    elems.append(Paragraph(f"<b>Confidence Level:</b> {conf_label}", normal))
    elems.append(Spacer(1, 12))

    elems.append(Paragraph("Key Issues", styles["Section"]))
    issues = analysis.get('issues', []) or []
    if issues:
        for it in issues[:4]:
            prob, why = _explain_issue_pair(it)
            elems.append(Paragraph(f"• {prob}", normal))
            elems.append(Paragraph(why, normal))
    else:
        elems.append(Paragraph("• Keine kritischen Probleme identifiziert; Monitoring empfohlen.", normal))
    elems.append(Spacer(1, 10))

    elems.append(Paragraph("Recommended Improvements", styles["Section"]))
    imps = analysis.get('improvements', []) or []
    if imps:
        for imp in imps[:5]:
            action, impact, reason = _explain_improvement(imp)
            elems.append(Paragraph(f"• {action}", normal))
            elems.append(Paragraph(impact, normal))
            elems.append(Paragraph(reason, normal))
    else:
        elems.append(Paragraph("• Keine spezifischen Verbesserungen identifiziert. Erwägen Sie ein detailliertes Audit.", normal))
    elems.append(PageBreak())

    # PAGE 2 — Process Analysis
    elems.append(Paragraph("Process Analysis", styles["ExecTitle"]))
    elems.append(Spacer(1, 8))

    input_params = blueprint.get('input_parameters') or {}
    has_adv_pur = bool(input_params.get('has_advanced_purification'))

    elems.append(Paragraph("Efficiency", styles["Section"]))
    elems.append(Paragraph(f"Aktuelle Effizienz: {eff}.", normal))
    elems.append(Paragraph("Primäre Treiber: Syntheseschrittanzahl und Prozesskomplexität.", normal))

    elems.append(Paragraph("Cost Structure", styles["Section"]))
    contributors = []
    if str(input_params.get('raw_material_cost') or '').lower() == 'high' or str((blueprint.get('analysis') or {}).get('cost')) == 'high':
        contributors.append('Rohstoffkosten')
    if not has_adv_pur and str((blueprint.get('input_parameters') or {}).get('desired_purity') or '').lower() in ('>99%', 'very high', 'very high'):
        contributors.append('Aufreinigung')
    if input_params.get('strict_waste_constraints'):
        contributors.append('Abfallbehandlung')
    if contributors:
        elems.append(Paragraph('Hauptkostentreiber: ' + ', '.join(contributors) + '.', normal))
    else:
        elems.append(Paragraph('Hauptkostentreiber: Aufreinigung, Rohstoffe und operative Overheads.', normal))

    elems.append(Paragraph("Cost Impact Analysis", styles["Section"]))
    cost_level = (analysis.get('costLevel') or analysis.get('cost') or 'n/a')
    elems.append(Paragraph(f"Cost Level: <b>{str(cost_level).upper()}</b>", normal))
    drivers = analysis.get('costDrivers') or []
    if drivers:
        elems.append(Paragraph("Main Cost Drivers:", normal))
        for d in drivers[:6]:
            elems.append(Paragraph(f"• {d}", normal))
    else:
        elems.append(Paragraph("Main Cost Drivers: Keine identifiziert.", normal))
    savings = analysis.get('savingsPotential') or 'low'
    elems.append(Paragraph(f"Estimated Savings Potential: <b>{str(savings).upper()}</b>", normal))

    # Downstream insights (type-specific guidance)
    downstream_insights = (analysis.get('downstream_insights') or [])
    if downstream_insights:
        elems.append(Spacer(1, 8))
        elems.append(Paragraph("Downstream Insights", styles["Section"]))
        for di in downstream_insights[:6]:
            elems.append(Paragraph(f"• {_safe_val(di)}", normal))

    elems.append(PageBreak())

    # PAGE 3 — Risks & Trade-offs
    elems.append(Paragraph("Risks & Trade-offs", styles["ExecTitle"]))
    elems.append(Spacer(1, 8))

    key_risks = (risks or []) + analysis.get('issues', [])
    seen = set()
    dedup_risks = []
    for r in key_risks:
        if r not in seen:
            dedup_risks.append(r)
            seen.add(r)
    dedup_risks = dedup_risks[:6]

    elems.append(Paragraph("Key Risks", styles["Section"]))
    if dedup_risks:
        for r in dedup_risks:
            elems.append(Paragraph(f"• {_safe_val(r)}", normal))
    else:
        elems.append(Paragraph("• Keine signifikanten Risiken identifiziert.", normal))

    elems.append(Spacer(1, 10))
    elems.append(Paragraph("Trade-offs", styles["Section"]))
    # Prefer tradeoffs generated by the analysis engine; fall back to sensible defaults.
    default_trade_lines = [
        "Kosten vs. Reinheit: Höhere Reinheit erhöht typischerweise Aufreinigungskosten.",
        "Effizienz vs. Komplexität: Mehr Syntheseschritte erhöhen Risiko und reduzieren Nettoeffizienz.",
    ]
    dynamic_tradeoffs = analysis.get("tradeoffs") if isinstance(analysis.get("tradeoffs"), (list, tuple)) else None
    for t in (dynamic_tradeoffs or default_trade_lines):
        elems.append(Paragraph(f"• {t}", normal))

    elems.append(PageBreak())

    # PAGE 4 — Recommended Actions
    elems.append(Paragraph("Recommended Actions", styles["ExecTitle"]))
    elems.append(Spacer(1, 8))
    actions = [
        "Bewerten Sie alternative Syntheserouten zur Reduktion der Schrittanzahl und Komplexität.",
        "Testen Sie verbesserte Aufreinigungsstrategien (z. B. Kristallisation) im Pilotmaßstab.",
        "Analysieren und optimieren Sie Rohstofflieferanten, um Kosten zu senken.",
        "Validieren Sie Prozessstabilität und Leistung in Pilotversuchen.",
    ]
    for idx, a in enumerate(actions, start=1):
        elems.append(Paragraph(f"{idx}. {a}", normal))

    elems.append(Spacer(1, 12))
    elems.append(Paragraph(f"Confidence Level: {conf_label}", styles["BodySmall"]))

    gen_date = __import__("datetime").datetime.utcnow().strftime('%Y-%m-%d %H:%M UTC')
    elems.append(Spacer(1, 18))
    elems.append(Paragraph(f"Generated: {gen_date}", styles["BodySmall"]))

    doc.build(elems, onFirstPage=draw_header_footer, onLaterPages=draw_header_footer)

    return path
