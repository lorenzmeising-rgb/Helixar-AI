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

def export_report_pdf(
    blueprint: Dict[str, Any],
    path: str,
    title: Optional[str] = None,
    author: Optional[str] = None,
    extras: Optional[Dict[str, Any]] = None,
    lang: Optional[str] = None,
) -> str:
    """Export a focused German PDF report from a blueprint.

    Produces a short, consultant-style PDF with header/footer:
    - Executive Summary (incl. Plausibilitätshinweise wenn vorhanden)
    - Process Analysis
    - Konkrete Empfehlungen (wenn extras['concrete_recs'] vorhanden)
    - Risks & Trade-offs
    - Empfohlene Optimierungen (aus extras['recommended_actions'] oder Fallback)

    `extras` darf folgende Schlüssel enthalten:
      - 'concrete_recs': dict mit production_route, expected_yield, downstream_steps, ...
      - 'recommended_actions': list von dicts mit title/rationale/expected_impact/...
      - 'plausibility_warnings': list von strings

    Returns the path to the written file.
    """
    extras = extras or {}
    # Resolve language: explicit param > extras['lang'] > default 'de'
    from i18n import t as _t, DEFAULT_LANGUAGE
    from engine_translations import translate_engine_text, translate_engine_value
    if lang is None:
        lang = extras.get("lang") or DEFAULT_LANGUAGE

    def L(key: str) -> str:
        return _t(key, lang=lang)

    def TE(text):
        """Translate dynamic engine output text to the active language."""
        return translate_engine_text(text, lang=lang)

    def TV(value):
        """Translate qualitative value labels (HIGH / LOW / MEDIUM)."""
        return translate_engine_value(value, lang=lang)
    try:
        from reportlab.lib.pagesizes import A4
        from reportlab.lib.units import mm
        from reportlab.lib import colors
        from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, PageBreak
        from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    except Exception as exc:
        raise ImportError("reportlab is required for PDF export. Install it with: pip install reportlab") from exc

    # Helixar brand palette — derived from logo (marine blue + turquoise)
    _BRAND_PRIMARY_HEX = "#1E3A5F"
    _BRAND_ACCENT_HEX = "#3DB8C5"
    _BRAND_TEXT_HEX = "#1A1A1A"
    _BRAND_MUTED_HEX = "#5A6B7B"

    styles = getSampleStyleSheet()
    styles.add(
        ParagraphStyle(
            name="ExecTitle",
            parent=styles["Heading1"],
            fontSize=18,
            leading=22,
            spaceAfter=12,
            alignment=1,
            textColor=colors.HexColor(_BRAND_PRIMARY_HEX),
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
            textColor=colors.HexColor(_BRAND_PRIMARY_HEX),
        )
    )
    styles.add(
        ParagraphStyle(name="DocSubtitle", parent=styles["Heading2"], fontSize=12, leading=14, alignment=0,
                       textColor=colors.HexColor(_BRAND_MUTED_HEX))
    )
    styles.add(
        ParagraphStyle(name="Recommendation", parent=styles["Heading1"], fontSize=20, leading=24, alignment=0,
                       textColor=colors.HexColor(_BRAND_PRIMARY_HEX))
    )
    styles.add(ParagraphStyle(name="Score", parent=styles["Heading2"], fontSize=14, leading=18, alignment=0,
                              textColor=colors.HexColor(_BRAND_TEXT_HEX)))
    styles.add(ParagraphStyle(name="Section", parent=styles["Heading2"], fontSize=13, leading=16, spaceAfter=6, spaceBefore=8,
                              textColor=colors.HexColor(_BRAND_PRIMARY_HEX)))
    styles.add(ParagraphStyle(name="BodySmall", parent=styles["Normal"], fontSize=10, leading=14,
                              textColor=colors.HexColor(_BRAND_TEXT_HEX)))
    # Monospace style for SMILES representation
    styles.add(ParagraphStyle(name="Mono", parent=styles["Normal"], fontName="Courier", fontSize=10, leading=12,
                              textColor=colors.HexColor(_BRAND_TEXT_HEX)))

    normal = styles["BodySmall"]

    rec = blueprint.get("recommended", {})
    meta = blueprint.get("metadata", {})
    risks = blueprint.get("risk_notes", []) or []

    def _fmt_num(v: Optional[Any]) -> str:
        try:
            return f"{float(v):.2f}"
        except Exception:
            return _safe_val(v)

    def _translate_method_display(m: Optional[Any], target_lang: str = "de") -> str:
        if m is None:
            return _safe_val(m)
        key = str(m).lower()
        if "biotech" in key or "biotechn" in key:
            return _t("pdf_method_biotech", lang=target_lang)
        if "chem" in key or "chemical" in key:
            return _t("pdf_method_chemical", lang=target_lang)
        if "extract" in key:
            return _t("pdf_method_extraction", lang=target_lang)
        return _safe_val(m)

    def _translate_process_display(p: Optional[Any], target_lang: str = "de") -> str:
        if p is None:
            return _safe_val(p)
        key = str(p).lower()
        if "synth" in key or "synt" in key:
            return _t("pdf_process_synthesis", lang=target_lang)
        if "ferment" in key:
            return _t("pdf_process_fermentation", lang=target_lang)
        if "extract" in key:
            return _t("pdf_process_extraction", lang=target_lang)
        return _safe_val(p)

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

    # Resolve the brand logo once (picked up automatically when present in repo).
    import os as _os
    _logo_dir = _os.path.dirname(_os.path.abspath(__file__))
    _logo_path_pdf = _os.path.join(_logo_dir, "helixar_logo.png")
    _has_logo_pdf = _os.path.exists(_logo_path_pdf)

    # Brand colors derived from the Helixar logo.
    BRAND_PRIMARY = colors.HexColor("#1E3A5F")  # marine blue
    BRAND_ACCENT = colors.HexColor("#3DB8C5")  # turquoise
    BRAND_MUTED = colors.HexColor("#5A6B7B")
    BRAND_DIVIDER = colors.HexColor("#DDE5EC")

    def draw_header_footer(canvas, doc):
        canvas.saveState()
        # Brand wordmark on the left
        canvas.setFont("Helvetica-Bold", 11)
        canvas.setFillColor(BRAND_PRIMARY)
        canvas.drawString(18 * mm, A4[1] - 18 * mm + 6, "Helixar")
        # Accent "AI" suffix in turquoise — matches the logo styling
        helixar_width = canvas.stringWidth("Helixar", "Helvetica-Bold", 11)
        canvas.setFillColor(BRAND_ACCENT)
        canvas.drawString(18 * mm + helixar_width + 2, A4[1] - 18 * mm + 6, "AI")
        # Subtitle under the wordmark
        canvas.setFont("Helvetica", 8)
        canvas.setFillColor(BRAND_MUTED)
        canvas.drawString(18 * mm, A4[1] - 18 * mm - 6, L("pdf_brand_subtitle"))
        # Logo on the top-right corner (small)
        if _has_logo_pdf:
            try:
                # Render the PNG at the top-right; size tuned to ~14 mm tall
                logo_h = 14 * mm
                logo_w = 14 * mm   # square-ish placement; aspect handled by reportlab
                canvas.drawImage(
                    _logo_path_pdf,
                    A4[0] - 18 * mm - logo_w,
                    A4[1] - 22 * mm,
                    width=logo_w,
                    height=logo_h,
                    preserveAspectRatio=True,
                    anchor="ne",
                    mask="auto",
                )
            except Exception:
                pass
        # Header divider
        canvas.setStrokeColor(BRAND_DIVIDER)
        canvas.setLineWidth(0.5)
        canvas.line(18 * mm, A4[1] - 24 * mm, A4[0] - 18 * mm, A4[1] - 24 * mm)
        # Footer
        canvas.setFont("Helvetica", 8)
        canvas.setFillColor(BRAND_MUTED)
        canvas.drawString(18 * mm, 12 * mm, L("pdf_footer_disclaimer"))
        gen_date = __import__("datetime").datetime.utcnow().strftime('%Y-%m-%d %H:%M UTC')
        canvas.drawRightString(A4[0] - 18 * mm, 12 * mm, f"{L('pdf_footer_generated')}: {gen_date}")
        canvas.restoreState()

    doc = SimpleDocTemplate(path, pagesize=A4, leftMargin=18 * mm, rightMargin=18 * mm, topMargin=36 * mm, bottomMargin=28 * mm)

    elems: List[Any] = []

    compound_title = _safe_val(rec.get('compound_name') or meta.get('compound_name') or (blueprint.get('input_parameters') or {}).get('target_molecule'))
    if title is None:
        title = f"Produktionsstrategie Bericht — {compound_title}"

    # PAGE 1 — Executive Summary
    elems.append(Paragraph("Helixar AI", styles["HelixarTitle"]))  # brand name — no translation
    elems.append(Spacer(1, 6))
    elems.append(Paragraph(L("pdf_subtitle"), styles["DocSubtitle"]))
    elems.append(Spacer(1, 10))

    method = rec.get('method') or rec.get('microorganism')
    ptype = rec.get('process_type') or rec.get('strain')

    current_date = __import__("datetime").datetime.utcnow().strftime('%Y-%m-%d')
    meta_table = [
        f"{L('pdf_target_molecule')}: {compound_title}",
        f"{L('pdf_method')}: {_translate_method_display(method, lang)}",
        f"{L('pdf_process_type')}: {_translate_process_display(ptype, lang)}",
        f"{L('pdf_date')}: {current_date}",
    ]
    for m in meta_table:
        elems.append(Paragraph(m, normal))
    elems.append(Spacer(1, 12))

    # Molecule representation: SMILES for small molecules, sequence for short
    # peptides, or external ID + descriptor for proteins/large biomolecules.
    # Prefer the full DB entry passed through from app.py; fall back to a
    # name-based lookup, then to the legacy SMILES string when the DB has no entry.
    raw_molecule = rec.get('compound_name') or meta.get('compound_name') or (blueprint.get('input_parameters') or {}).get('target_molecule')
    elems.append(Spacer(1, 6))
    db_entry = (blueprint.get('input_parameters') or {}).get('db_entry')
    if not db_entry:
        try:
            from molecules_db import get_entry_by_name as _get_entry
            db_entry = _get_entry(raw_molecule) if raw_molecule else None
        except Exception:
            db_entry = None
    if db_entry:
        try:
            from molecules_db import get_representation_label as _rep_label
            rep_text = _rep_label(db_entry)
        except Exception:
            rep_text = None
        if rep_text:
            # Use monospace only for actual SMILES/sequence strings; descriptive
            # text reads better in the body style.
            if rep_text.startswith("SMILES:") or rep_text.startswith("Sequenz:"):
                elems.append(Paragraph(rep_text, styles["Mono"]))
            else:
                elems.append(Paragraph(f"<b>{L('pdf_representation')}:</b> {rep_text}", normal))
        else:
            elems.append(Paragraph(f"SMILES: {_get_smiles_for(raw_molecule)}", styles["Mono"]))
    else:
        # No DB entry — show the legacy SMILES lookup result for backward compat.
        smiles = _get_smiles_for(raw_molecule)
        elems.append(Paragraph(f"SMILES: {smiles}", styles["Mono"]))
    elems.append(Spacer(1, 8))

    analysis = blueprint.get("analysis") or {}

    elems.append(Paragraph(L("pdf_executive_summary"), styles["Recommendation"]))
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

    elems.append(Paragraph(f"<b>{L('pdf_efficiency')}:</b> {perf_icon(eff)} {TV(str(eff).upper())}", styles["Score"]))
    elems.append(Paragraph(f"<b>{L('pdf_cost')}:</b> {perf_icon(cost_val)} {TV(str(cost_val).upper())}", normal))
    elems.append(Paragraph(f"<b>{L('pdf_risk')}:</b> {perf_icon(risk_val)} {TV(str(risk_val).upper())}", normal))
    elems.append(Paragraph(f"<b>{L('pdf_toxicity')}:</b> {perf_icon(tox)} {TV(str(tox).upper())}", normal))
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
    elems.append(Paragraph(f"<b>{L('pdf_confidence')}:</b> {TV(conf_label)}", normal))
    elems.append(Spacer(1, 12))

    elems.append(Paragraph(L("pdf_key_issues"), styles["Section"]))
    # Engine-produced issues are in English; translate inline.
    issues = analysis.get('issues', []) or []
    if issues:
        for it in issues[:4]:
            prob, why = _explain_issue_pair(it)
            # Translate the (possibly English) engine output before display.
            elems.append(Paragraph(f"• {TE(prob)}", normal))
            elems.append(Paragraph(TE(why), normal))
    else:
        elems.append(Paragraph(L("pdf_no_critical_issues"), normal))
    elems.append(Spacer(1, 10))

    elems.append(Paragraph(L("pdf_recommended_improvements"), styles["Section"]))
    imps = analysis.get('improvements', []) or []
    if imps:
        for imp in imps[:5]:
            action, impact, reason = _explain_improvement(imp)
            elems.append(Paragraph(f"• {TE(action)}", normal))
            elems.append(Paragraph(TE(impact), normal))
            elems.append(Paragraph(TE(reason), normal))
    else:
        elems.append(Paragraph(L("pdf_no_improvements"), normal))

    # Plausibility warnings (C8) — surface unrealistic input combinations
    plaus = extras.get("plausibility_warnings") or []
    if plaus:
        elems.append(Spacer(1, 10))
        elems.append(Paragraph(L("pdf_plausibility_section"), styles["Section"]))
        for w in plaus:
            elems.append(Paragraph(f"• {_safe_val(w)}", normal))
    elems.append(PageBreak())

    # PAGE 2 — Process Analysis
    elems.append(Paragraph(L("pdf_process_analysis"), styles["ExecTitle"]))
    elems.append(Spacer(1, 8))

    input_params = blueprint.get('input_parameters') or {}
    has_adv_pur = bool(input_params.get('has_advanced_purification'))

    # Numeric input summary — make sure the report reflects the exact values
    # the user entered (not just the qualitative buckets).
    pct = input_params.get("desired_purity_percent")
    kg = input_params.get("scale_kg_per_year")
    eur = input_params.get("raw_material_cost_eur_per_kg")
    n_supp = input_params.get("num_qualified_suppliers")
    lead = input_params.get("lead_time_weeks")
    single_region = input_params.get("single_region_concentration")
    has_any_numeric = any(v is not None for v in (pct, kg, eur, n_supp, lead, single_region))
    if has_any_numeric:
        elems.append(Paragraph(L("pdf_quantitative_inputs"), styles["Section"]))
        if pct is not None:
            elems.append(Paragraph(f"<b>{L('pdf_target_purity')}:</b> {pct:.2f} %", normal))
        if kg is not None:
            kg_disp = f"{kg:,.0f}" if kg >= 1000 else f"{kg:.3f}"
            kg_unit = "kg/year" if lang == "en" else "kg/Jahr"
            elems.append(Paragraph(f"<b>{L('pdf_target_scale')}:</b> {kg_disp} {kg_unit}", normal))
        if eur is not None:
            elems.append(Paragraph(f"<b>{L('pdf_raw_material_price')}:</b> {eur:,.2f} €/kg", normal))
        if n_supp is not None:
            elems.append(Paragraph(f"<b>{L('pdf_qualified_suppliers')}:</b> {int(n_supp)}", normal))
        if lead is not None:
            weeks_unit = "weeks" if lang == "en" else "Wochen"
            elems.append(Paragraph(f"<b>{L('pdf_lead_time')}:</b> {lead:.1f} {weeks_unit}", normal))
        if single_region is not None:
            geo = L("pdf_geo_concentrated") if single_region else L("pdf_geo_diversified")
            elems.append(Paragraph(f"<b>{L('pdf_supplier_geography')}:</b> {geo}", normal))
        # Biophysical numerics (only when peptide/protein gave values)
        bio_tagg = input_params.get("tagg_celsius")
        bio_tm = input_params.get("tm_celsius")
        bio_disulfides = input_params.get("num_disulfides")
        bio_domains = input_params.get("num_domains")
        bio_ptm = input_params.get("has_ptm")
        if any(v is not None for v in (bio_tagg, bio_tm, bio_disulfides, bio_domains, bio_ptm)):
            elems.append(Spacer(1, 4))
            elems.append(Paragraph(f"<b>{L('biophysical_inputs_header')}:</b>", normal))
            if bio_tagg is not None:
                elems.append(Paragraph(f"<b>{L('metric_tagg')}:</b> {bio_tagg:.1f} °C", normal))
            if bio_tm is not None:
                elems.append(Paragraph(f"<b>{L('metric_tm')}:</b> {bio_tm:.1f} °C", normal))
            if bio_disulfides is not None:
                elems.append(Paragraph(f"<b>{L('metric_disulfides')}:</b> {int(bio_disulfides)}", normal))
            if bio_domains is not None:
                elems.append(Paragraph(f"<b>{L('metric_domains')}:</b> {int(bio_domains)}", normal))
            if bio_ptm is not None:
                ptm_text = (L("yes") if bio_ptm else L("no"))
                elems.append(Paragraph(f"<b>{L('metric_ptm')}:</b> {ptm_text}", normal))
        elems.append(Spacer(1, 6))

    # Structure analysis (RDKit) — only when a SMILES is available and RDKit
    # successfully parses it. Falls back silently when not.
    smiles_for_pdf = input_params.get("smiles")
    rdkit_section_added = False
    if smiles_for_pdf:
        try:
            from rdkit_properties import (
                compute_properties as _rdk_compute,
                render_structure_png as _rdk_render,
                RDKIT_AVAILABLE as _RDK_AVAIL,
            )
            if _RDK_AVAIL:
                rdk_p = _rdk_compute(smiles_for_pdf)
                if rdk_p:
                    elems.append(Paragraph(L("structure_analysis_title"), styles["Section"]))
                    elems.append(Paragraph(L("structure_caption"), normal))
                    elems.append(Paragraph(f"<b>{L('structure_label_mw')}:</b> {rdk_p['molecular_weight']:.2f} g/mol", normal))
                    elems.append(Paragraph(f"<b>{L('structure_label_logp')}:</b> {rdk_p['logp']:.2f}", normal))
                    elems.append(Paragraph(f"<b>{L('structure_label_tpsa')}:</b> {rdk_p['tpsa']:.2f} Å²", normal))
                    elems.append(Paragraph(f"<b>{L('structure_label_arom_rings')}:</b> {rdk_p['num_aromatic_rings']}", normal))
                    elems.append(Paragraph(f"<b>{L('structure_label_hbd')}:</b> {rdk_p['num_h_donors']}  ·  <b>{L('structure_label_hba')}:</b> {rdk_p['num_h_acceptors']}", normal))
                    elems.append(Paragraph(f"<b>{L('structure_label_rotbonds')}:</b> {rdk_p['num_rotatable_bonds']}", normal))
                    if rdk_p['lipinski_violations'] == 0:
                        elems.append(Paragraph(f"✓ {L('structure_label_lipinski_pass')}", normal))
                    else:
                        elems.append(Paragraph(f"⚠ {L('structure_label_lipinski_fail')} ({rdk_p['lipinski_violations']})", normal))
                    # Embed the rendered structure image
                    try:
                        png = _rdk_render(smiles_for_pdf, width=350, height=250)
                        if png:
                            from reportlab.platypus import Image as _RLImage
                            from reportlab.lib.units import mm
                            import io as _io
                            img_buf = _io.BytesIO(png)
                            img = _RLImage(img_buf, width=70 * mm, height=50 * mm)
                            elems.append(Spacer(1, 4))
                            elems.append(img)
                            elems.append(Paragraph(f"<i>{L('structure_drawing_caption')}</i>", normal))
                    except Exception:
                        pass
                    elems.append(Spacer(1, 8))
                    rdkit_section_added = True
        except Exception:
            pass

    elems.append(Paragraph(L("pdf_efficiency"), styles["Section"]))
    eff_intro = "Current efficiency" if lang == "en" else "Aktuelle Effizienz"
    elems.append(Paragraph(f"{eff_intro}: {TV(str(eff))}.", normal))
    elems.append(Paragraph(L("pdf_primary_drivers"), normal))

    elems.append(Paragraph(L("pdf_cost_structure"), styles["Section"]))
    contributors = []
    if str(input_params.get('raw_material_cost') or '').lower() == 'high' or str((blueprint.get('analysis') or {}).get('cost')) == 'high':
        contributors.append('Rohstoffkosten')
    if not has_adv_pur and str((blueprint.get('input_parameters') or {}).get('desired_purity') or '').lower() in ('>99%', 'very high', 'very high'):
        contributors.append('Aufreinigung')
    if input_params.get('strict_waste_constraints'):
        contributors.append('Abfallbehandlung')
    if contributors:
        # Localize the contributor labels for English output.
        if lang == "en":
            translation = {
                "Rohstoffkosten": "raw-material cost",
                "Aufreinigung": "purification",
                "Abfallbehandlung": "waste handling",
            }
            loc_contribs = [translation.get(c, c) for c in contributors]
        else:
            loc_contribs = contributors
        prefix = L("pdf_main_cost_drivers")
        elems.append(Paragraph(prefix + ' ' + ', '.join(loc_contribs) + '.', normal))
    else:
        elems.append(Paragraph(L("pdf_main_cost_default"), normal))

    elems.append(Paragraph(L("pdf_cost_impact"), styles["Section"]))
    cost_level = (analysis.get('costLevel') or analysis.get('cost') or 'n/a')
    elems.append(Paragraph(f"{L('pdf_cost_level')}: <b>{TV(str(cost_level).upper())}</b>", normal))
    drivers = analysis.get('costDrivers') or []
    if drivers:
        elems.append(Paragraph(L("pdf_main_cost_drivers"), normal))
        for d in drivers[:6]:
            elems.append(Paragraph(f"• {TE(d)}", normal))
    else:
        elems.append(Paragraph(L("pdf_no_cost_drivers"), normal))
    savings = analysis.get('savingsPotential') or 'low'
    elems.append(Paragraph(f"{L('pdf_savings_potential')}: <b>{TV(str(savings).upper())}</b>", normal))

    # Downstream insights (type-specific guidance)
    downstream_insights = (analysis.get('downstream_insights') or [])
    if downstream_insights:
        elems.append(Spacer(1, 8))
        elems.append(Paragraph(L("pdf_downstream_insights"), styles["Section"]))
        for di in downstream_insights[:6]:
            elems.append(Paragraph(f"• {TE(_safe_val(di))}", normal))

    elems.append(PageBreak())

    # PAGE — Konkrete Empfehlungen (C9)
    concrete = extras.get("concrete_recs") or {}
    if concrete:
        elems.append(Paragraph(L("pdf_concrete_recommendations"), styles["ExecTitle"]))
        elems.append(Spacer(1, 8))
        if concrete.get("production_route"):
            elems.append(Paragraph(L("pdf_production_route"), styles["Section"]))
            elems.append(Paragraph(_safe_val(concrete.get("production_route")), normal))
        if concrete.get("expected_yield"):
            elems.append(Paragraph(L("pdf_expected_yield"), styles["Section"]))
            elems.append(Paragraph(_safe_val(concrete.get("expected_yield")), normal))
        if concrete.get("processing_time"):
            elems.append(Paragraph(L("pdf_processing_time"), styles["Section"]))
            elems.append(Paragraph(_safe_val(concrete.get("processing_time")), normal))
        if concrete.get("expected_final_purity"):
            elems.append(Paragraph(L("pdf_achievable_purity"), styles["Section"]))
            elems.append(Paragraph(_safe_val(concrete.get("expected_final_purity")), normal))
        steps = concrete.get("downstream_steps") or []
        if steps:
            elems.append(Paragraph(L("pdf_downstream_steps"), styles["Section"]))
            for i, s in enumerate(steps, 1):
                elems.append(Paragraph(f"{i}. {_safe_val(s)}", normal))
        notes = concrete.get("notes") or []
        if notes:
            elems.append(Spacer(1, 6))
            elems.append(Paragraph(L("pdf_notes"), styles["Section"]))
            for n in notes:
                elems.append(Paragraph(f"• {_safe_val(n)}", normal))
        elems.append(PageBreak())

    # PAGE 3 — Risks & Trade-offs
    elems.append(Paragraph(L("pdf_risks_tradeoffs"), styles["ExecTitle"]))
    elems.append(Spacer(1, 8))

    key_risks = (risks or []) + analysis.get('issues', [])
    seen = set()
    dedup_risks = []
    for r in key_risks:
        if r not in seen:
            dedup_risks.append(r)
            seen.add(r)
    dedup_risks = dedup_risks[:6]

    elems.append(Paragraph(L("pdf_key_risks"), styles["Section"]))
    if dedup_risks:
        for r in dedup_risks:
            elems.append(Paragraph(f"• {TE(_safe_val(r))}", normal))
    else:
        elems.append(Paragraph(L("pdf_no_significant_risks"), normal))

    elems.append(Spacer(1, 10))
    elems.append(Paragraph(L("pdf_tradeoffs"), styles["Section"]))
    # Prefer tradeoffs generated by the analysis engine; fall back to sensible defaults.
    default_trade_lines = [L("pdf_tradeoff_default_1"), L("pdf_tradeoff_default_2")]
    dynamic_tradeoffs = analysis.get("tradeoffs") if isinstance(analysis.get("tradeoffs"), (list, tuple)) else None
    if dynamic_tradeoffs:
        # Engine-produced tradeoffs are English; translate inline.
        for tl in dynamic_tradeoffs:
            elems.append(Paragraph(f"• {TE(tl)}", normal))
    else:
        for tl in default_trade_lines:
            elems.append(Paragraph(f"• {tl}", normal))

    elems.append(PageBreak())

    # PAGE 4 — Empfohlene Optimierungen (C10)
    elems.append(Paragraph(L("pdf_recommended_actions"), styles["ExecTitle"]))
    elems.append(Spacer(1, 8))
    rec_actions = extras.get("recommended_actions") or []
    effort_label_keys = {"low": "effort_low", "medium": "effort_medium", "high": "effort_high"}
    if rec_actions:
        for idx, a in enumerate(rec_actions, start=1):
            elems.append(Paragraph(f"{idx}. {_safe_val(a.get('title'))}", styles["Section"]))
            if a.get("rationale"):
                elems.append(Paragraph(f"<b>{L('pdf_rationale')}:</b> {_safe_val(a.get('rationale'))}", normal))
            # Current vs. optimized state — rendered as labelled paragraphs
            if a.get("current_state"):
                elems.append(Paragraph(
                    f"<b>{L('actions_current_state')}:</b> {_safe_val(a.get('current_state'))}",
                    normal,
                ))
            if a.get("optimized_state"):
                elems.append(Paragraph(
                    f"<b>{L('actions_optimized_state')}:</b> {_safe_val(a.get('optimized_state'))}",
                    normal,
                ))
            if a.get("expected_impact"):
                elems.append(Paragraph(f"<b>{L('pdf_impact')}:</b> {_safe_val(a.get('expected_impact'))}", normal))
            if a.get("prerequisites"):
                elems.append(Paragraph(f"<b>{L('pdf_prerequisites')}:</b> {_safe_val(a.get('prerequisites'))}", normal))
            if a.get("effort"):
                eff_key = effort_label_keys.get(str(a["effort"]).lower())
                eff_label = L(eff_key) if eff_key else _safe_val(a.get("effort"))
                elems.append(Paragraph(f"<b>{L('pdf_effort')}:</b> {eff_label}", normal))
            elems.append(Spacer(1, 6))
    else:
        # Fallback to generic actions if no rule-based actions matched.
        for idx, key in enumerate(
            ("pdf_fallback_action_1", "pdf_fallback_action_2", "pdf_fallback_action_3", "pdf_fallback_action_4"),
            start=1,
        ):
            elems.append(Paragraph(f"{idx}. {L(key)}", normal))

    elems.append(Spacer(1, 12))
    elems.append(Paragraph(f"{L('pdf_confidence')}: {TV(conf_label)}", styles["BodySmall"]))

    gen_date = __import__("datetime").datetime.utcnow().strftime('%Y-%m-%d %H:%M UTC')
    elems.append(Spacer(1, 18))
    elems.append(Paragraph(f"{L('pdf_generated')}: {gen_date}", styles["BodySmall"]))

    # ----- Penultimate page: References / Literature -----
    # Aggregate literature_source entries from concrete_recs, recommended_actions
    # and plausibility checks. Deduplicate while preserving order.
    aggregated_sources = []
    seen_src = set()

    def _add_source(src):
        if src and src not in seen_src:
            seen_src.add(src)
            aggregated_sources.append(src)

    cr = extras.get("concrete_recs") or {}
    for s in (cr.get("literature_sources") or []):
        _add_source(s)
    for action in (extras.get("recommended_actions") or []):
        _add_source(action.get("literature_source"))
    plaus = extras.get("plausibility_warnings") or []
    if isinstance(plaus, dict):
        # Newer call style: {warnings, literature_sources}
        for s in (plaus.get("literature_sources") or []):
            _add_source(s)
    plaus_sources = extras.get("plausibility_sources") or []
    for s in plaus_sources:
        _add_source(s)

    if aggregated_sources:
        elems.append(PageBreak())
        elems.append(Paragraph(L("references_title"), styles["ExecTitle"]))
        elems.append(Spacer(1, 8))
        elems.append(Paragraph(L("references_intro"), normal))
        elems.append(Spacer(1, 6))
        for idx, src in enumerate(aggregated_sources, start=1):
            elems.append(Paragraph(f"[{idx}] {src}", normal))
            elems.append(Spacer(1, 3))

    # ----- Final page: Important Notices (legal / IP / liability / sources) -----
    elems.append(PageBreak())
    elems.append(Paragraph(L("notices_title"), styles["ExecTitle"]))
    elems.append(Spacer(1, 10))
    for title_key, body_key in (
        ("notice_compliance_title", "notice_compliance_body"),
        ("notice_data_title", "notice_data_body"),
        ("notice_ip_title", "notice_ip_body"),
        ("notice_liability_title", "notice_liability_body"),
        ("notice_sources_title", "notice_sources_body"),
    ):
        elems.append(Paragraph(L(title_key), styles["Section"]))
        elems.append(Paragraph(L(body_key), normal))
        elems.append(Spacer(1, 8))

    doc.build(elems, onFirstPage=draw_header_footer, onLaterPages=draw_header_footer)

    return path
