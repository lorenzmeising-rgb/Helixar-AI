from typing import Dict, Any, List, Optional

from language_de import LABELS


from molecules_db import get_smiles_for as _get_smiles_for


# Bug C1 fix: locale-aware number formatting.
#
# The old `f"{kg:,.0f}" if kg >= 1000 else f"{kg:.3f}"` produced two
# disasters in German PDFs:
#   - "{100:.3f}" → "100.000"  (= three decimals — looks identical to
#     the German "einhunderttausend" thousands-separator format)
#   - "{5000:,.0f}" → "5,000"   (= US comma thousands-separator,
#     reads in German as 5 with decimals)
# Both make readers think the scale is off by 1000x.
#
# This helper produces canonical German / English number formats:
#   format_int_de(100)     → "100"
#   format_int_de(5000)    → "5.000"
#   format_int_de(100000)  → "100.000"
#   format_int_en(100000)  → "100,000"

def _format_int_de(v: float) -> str:
    """Format integer with German thousands separator (point)."""
    return f"{int(round(float(v))):,}".replace(",", ".")


def _format_int_en(v: float) -> str:
    """Format integer with English thousands separator (comma)."""
    return f"{int(round(float(v))):,}"


def _format_money_de(v: float) -> str:
    """Format currency value with German conventions (thousands=., decimal=,)."""
    # Python `:,.2f` produces "1,234.56" (US). Swap-convert to DE.
    raw = f"{float(v):,.2f}"
    # temporary placeholder, then swap
    return raw.replace(",", "X").replace(".", ",").replace("X", ".")


def _format_money_en(v: float) -> str:
    """Format currency value with English conventions (thousands=,, decimal=.)."""
    return f"{float(v):,.2f}"


# Bug H2 fix: ReportLab's default Helvetica only supports WinAnsi/CP1252.
# Latin Extended-A characters like 'ń' (U+0144) render as '■' (U+25A0).
# Pre-transliterate problematic characters to their nearest ASCII/Latin-1
# equivalent so author names like "Celińska" are rendered as "Celinska"
# instead of "Celi■ska".
_PDF_UNICODE_MAP = {
    # Polish / Czech / Slovak Latin Extended-A
    "ą": "a", "Ą": "A", "ć": "c", "Ć": "C", "ę": "e", "Ę": "E",
    "ł": "l", "Ł": "L", "ń": "n", "Ń": "N", "ó": "o", "Ó": "O",
    "ś": "s", "Ś": "S", "ź": "z", "Ź": "Z", "ż": "z", "Ż": "Z",
    "č": "c", "Č": "C", "ď": "d", "Ď": "D", "ě": "e", "Ě": "E",
    "ň": "n", "Ň": "N", "ř": "r", "Ř": "R", "š": "s", "Š": "S",
    "ť": "t", "Ť": "T", "ů": "u", "Ů": "U", "ž": "z", "Ž": "Z",
    # Turkish / Romanian
    "ş": "s", "Ş": "S", "ț": "t", "Ț": "T", "ğ": "g", "Ğ": "G",
    "ı": "i", "İ": "I",
    # Hungarian
    "ő": "o", "Ő": "O", "ű": "u", "Ű": "U",
    # Math / typography that may sneak in
    "≈": "~", "≤": "<=", "≥": ">=", "≠": "!=",
    "·": "-",
    # Greek letters that occasionally appear in chem nomenclature
    "α": "alpha", "β": "beta", "γ": "gamma", "δ": "delta",
    "Α": "Alpha", "Β": "Beta", "Γ": "Gamma", "Δ": "Delta",
    "μ": "u",  # often "mikro" in unit context
}


def _to_pdf_safe(text: Any) -> str:
    """Transliterate Latin-Extended-A and other non-CP1252 chars so they
    don't render as '■' in ReportLab's default Helvetica."""
    if text is None:
        return ""
    s = str(text)
    if not any(ord(c) > 127 for c in s):
        return s  # fast path
    out = []
    for c in s:
        if ord(c) <= 127:
            out.append(c)
            continue
        # already safe in CP1252?
        try:
            c.encode("cp1252")
            out.append(c)
            continue
        except UnicodeEncodeError:
            pass
        out.append(_PDF_UNICODE_MAP.get(c, "?"))
    return "".join(out)


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
        """Translate dynamic engine output text to the active language.

        Also sanitises Latin-Extended-A characters (Bug H2 fix) so author
        names like 'Celińska' no longer render as 'Celi■ska' in the PDF.
        """
        return _to_pdf_safe(translate_engine_text(text, lang=lang))

    def TV(value):
        """Translate qualitative value labels (HIGH / LOW / MEDIUM)."""
        return _to_pdf_safe(translate_engine_value(value, lang=lang))
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
    # Bug B8 fix: dedicated style for explanation/sub-lines under a bullet.
    # Rendered indented, italic, in muted grey — replaces the previous
    # "→ explanation" lines that surfaced as "? explanation" when the
    # arrow character was sanitised out for the CP1252 font.
    styles.add(ParagraphStyle(name="ExplanationIndent",
                              parent=styles["Normal"],
                              fontName="Helvetica-Oblique",
                              fontSize=9, leading=12,
                              leftIndent=14,
                              spaceAfter=2,
                              textColor=colors.HexColor("#6B7B8A")))
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
        # Bug B8 fix: removed the leading "→ " arrow from the explanation
        # strings. The arrow (U+2192) is outside CP1252 and was being
        # replaced by "?" via the PDF unicode sanitiser, producing the
        # "?-as-bullet" effect users complained about. Explanation lines
        # are now rendered in italic + indent (see the renderer below).
        if "step" in low or "synth" in low:
            reason = "Jede zusätzliche Syntheseschritt erhöht Kosten, Durchlaufzeit und Komplexität, reduziert die Effizienz."
        elif "purif" in low or "reinig" in low or "purificat" in low:
            reason = "Hohe Reinigungsanforderungen treiben nachgelagerte Arbeit und Kosten; mindern Skalierbarkeit."
        elif "raw material" in low or "rohstoff" in low or "expensive" in low or "teuer" in low:
            reason = "Teure Rohstoffe sind ein wiederkehrender Kostentreiber und beeinflussen die Stückkosten stark."
        elif "bioreactor" in low or "bioreaktor" in low:
            reason = "Fehlende Bioreaktor-Infrastruktur schränkt biotechnologische Routen ein und führt zu kostspieligen Anpassungen."
        elif "waste" in low or "abfall" in low:
            reason = "Strenge Abfallvorschriften erhöhen Compliance-Last und Betriebskosten, besonders bei gefährlichen Abfällen."
        elif "stability" in low or "stabil" in low:
            reason = "Geringe Stabilität erhöht Ausfallrisiko beim Scale-up und treibt Entwicklungskosten."
        else:
            reason = "Dieses Thema wirkt sich negativ auf Kosten, Risiko oder Skalierbarkeit aus und erfordert Untersuchung."
        return s, reason

    def _explain_improvement(impr: str):
        s = str(impr or "").strip()
        low = s.lower()
        # Bug B8 fix: arrows removed — see _explain_issue_pair above.
        if "reduce number" in low or "reduce the number" in low or "reduce number of" in low:
            action = "Reduzieren der Anzahl an Syntheseschritten"
            impact = "Hohes Einsparpotenzial und bessere Skalierbarkeit"
            reason = "Weniger Schritte verringern Materialeinsatz, Prozesszeit und Nachreinigung."
        elif "crystall" in low or "crystalliz" in low or "purificat" in low:
            action = "Optimierung der Reinigung (z. B. Kristallisation)"
            impact = "Mittleres Einsparpotenzial bei Aufreinigungskosten"
            reason = "Alternative Methoden reduzieren Lösungsmittelverbrauch und erhöhen Durchsatz."
        elif "raw material" in low or "optimiz" in low or "sourcing" in low:
            action = "Optimierung der Rohstoffauswahl und Beschaffung"
            impact = "Mittleres Einsparpotenzial auf Stückkosten"
            reason = "Bessere Beschaffung reduziert wiederkehrende Materialkosten."
        elif "stabil" in low or "process stability" in low:
            action = "Verbesserung der Prozessstabilität"
            impact = "Mittleres Einsparpotenzial durch reduzierte Nacharbeit"
            reason = "Stabilere Prozesse unterstützen Reproduzierbarkeit und Scale-up."
        elif "toxic" in low or "waste" in low or "hazard" in low:
            action = "Reduzieren von gefährlichen Reagenzien und Abfallströmen"
            impact = "Mittleres Einsparpotenzial bei Compliance- und Entsorgungskosten"
            reason = "Weniger gefährliche Chemikalien senken Behandlungskosten und regulatorischen Aufwand."
        else:
            action = s or "Empfohlene Maßnahme"
            impact = "Erwarteter Effekt: Moderat"
            reason = "Diese Maßnahme adressiert identifizierte Schwächen und sollte Prozesskennzahlen verbessern."
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

    # Color-coded severity indicator — replaces the previous emoji icons
    # (🔴 ⚠️ 🟢) which ReportLab's default Helvetica font cannot render and
    # would output as ■ squares. Pure text + color is more professional too.
    def perf_color(v: str, metric: str = "cost") -> str:
        """Pick a colour for the qualitative band.

        For cost/risk/toxicity: HIGH is bad (red), LOW is good (green).
        For efficiency: HIGH is good (green), LOW is bad (red) — inverted.
        """
        key = str(v).lower()
        is_high = key in ("high", "very high", "hoch", "sehr hoch")
        is_med = key in ("medium", "mittel")
        if metric == "efficiency":
            # efficiency: high = good
            if is_high:
                return "#1E8449"
            if is_med:
                return "#D68910"
            return "#C0392B"
        # cost / risk / toxicity: high = bad
        if is_high:
            return "#C0392B"
        if is_med:
            return "#D68910"
        return "#1E8449"

    # ---- Tier label map (replaces HOCH / MITTEL / NIEDRIG) -----------------
    # Bug B6 fix: HIGH/MEDIUM/LOW reads like a school traffic-light. CMC
    # / pharma reports use domain-specific tier names instead — these
    # carry information about WHAT the band actually means in practice.
    _TIER_LABELS_DE = {
        "efficiency": {
            "low":       "Kritisch",
            "medium":    "Optimierbar",
            "high":      "Robust",
            "very high": "Best-in-class",
        },
        "cost": {
            "low":       "Commodity-Niveau",
            "medium":    "Mid-Tier-API",
            "high":      "Specialty-Tier",
            "very high": "Biologika-Niveau",
        },
        "risk": {
            "low":       "Etabliert",
            "medium":    "Akzeptabel",
            "high":      "Überwachungsbedürftig",
            "very high": "Kritisch",
        },
        "toxicity": {
            "low":       "Standard-Handhabung",
            "medium":    "Standard-API-Containment",
            "high":      "Erhöhte Containment-Anforderung",
            "very high": "GMP-Cat-3-Containment",
        },
    }
    _TIER_LABELS_EN = {
        "efficiency": {
            "low": "Critical", "medium": "Optimisable",
            "high": "Robust", "very high": "Best-in-class",
        },
        "cost": {
            "low": "Commodity-tier", "medium": "Mid-tier API",
            "high": "Specialty-tier", "very high": "Biologic-tier",
        },
        "risk": {
            "low": "Established", "medium": "Acceptable",
            "high": "Monitoring-required", "very high": "Critical",
        },
        "toxicity": {
            "low": "Standard handling", "medium": "Standard API containment",
            "high": "Elevated containment", "very high": "GMP cat-3 containment",
        },
    }

    def _tier_label_text(v: str, metric: str) -> str:
        """Return the domain-specific tier name for a qualitative band."""
        key = str(v).lower().replace(" ", " ").strip()
        # Accept German + English values
        key_map = {
            "niedrig": "low", "mittel": "medium",
            "hoch": "high", "sehr hoch": "very high",
        }
        normalized = key_map.get(key, key)
        table = _TIER_LABELS_EN if lang == "en" else _TIER_LABELS_DE
        return table.get(metric, {}).get(normalized, str(v).title())

    def perf_value_markup(v: str, metric: str = "cost") -> str:
        """Return the tier label in bold brand-navy.

        Bug B7 fix: previously the tier word itself was coloured red /
        amber / green, which looked like a traffic-light rating on a
        single word and felt unprofessional. Now the tier word is
        rendered in the same brand-navy as other emphasised terms in the
        report, with no per-word colour coding.
        """
        tier = _tier_label_text(v, metric)
        return f'<font color="#1E3A5F"><b>{tier}</b></font>'

    # ---- Numeric score mapping (label → X/10) ----
    # Provides a concrete numerical anchor next to the qualitative bucket so
    # pharma reviewers can compare metrics on a uniform scale. Asymmetric for
    # efficiency (10 = best) vs. the other three (10 = worst), so all four
    # read consistently as "10 = ungünstig"-style is wrong — we keep
    # efficiency inverted because a high efficiency is GOOD.
    def _score_for(label: str, metric: str) -> str:
        key = str(label).lower()
        # For Effizienz: higher label = higher score (good is high)
        if metric == "efficiency":
            return {
                "low": "3", "niedrig": "3",
                "medium": "6", "mittel": "6",
                "high": "8", "hoch": "8",
                "very high": "10", "sehr hoch": "10",
            }.get(key, "—") + " / 10"
        # For Cost/Risk/Toxicity: higher label = higher score (high is BAD)
        return {
            "low": "2", "niedrig": "2",
            "medium": "5", "mittel": "5",
            "high": "7", "hoch": "7",
            "very high": "9", "sehr hoch": "9",
        }.get(key, "—") + " / 10"

    # ---- Brief reason helpers ----
    # Each returns a 3-7 word German/English phrase explaining the rating.
    # Pulls from input parameters + analysis data — same data the engine used.
    # Correct path: blueprint.get('input_parameters'), not analysis.input_context
    # (which is empty in current engine output).
    ip_ctx = blueprint.get('input_parameters') or {}
    mtype_ctx = (ip_ctx.get('molecule_type') or "").lower()
    msub_ctx = (ip_ctx.get('molecule_subtype') or "").lower()
    method_ctx = (ip_ctx.get('method') or "").lower()
    steps_ctx = ip_ctx.get('number_of_steps')
    raw_cost_ctx = (ip_ctx.get('raw_material_cost') or "").lower()
    eur_kg = ip_ctx.get('raw_material_cost_eur_per_kg')
    purity_pct = ip_ctx.get('desired_purity_percent')
    props_ctx = analysis.get('properties') or {}
    pur_diff = (props_ctx.get('purification_difficulty') or "").lower()
    complexity_ctx = (props_ctx.get('complexity') or "").lower()
    strict_waste_ctx = bool(ip_ctx.get('strict_waste_constraints'))
    n_dis = int(ip_ctx.get('num_disulfides') or 0)
    has_ptm_ctx = bool(ip_ctx.get('has_ptm'))

    def _brief_eff(v: str) -> str:
        key = str(v).lower()
        if isinstance(steps_ctx, int) and steps_ctx >= 10:
            return "lange Synthesekette" if lang == "de" else "long synthesis chain"
        if key in ("low", "niedrig"):
            return "ineffizienter Prozess-Pfad" if lang == "de" else "inefficient process path"
        if key in ("medium", "mittel"):
            return "moderate Komplexität / Schrittzahl" if lang == "de" else "moderate complexity / step count"
        if msub_ctx == "antibody":
            return "etablierte CHO-Plattform" if lang == "de" else "established CHO platform"
        if mtype_ctx == "natural_product":
            # Method-aware: chemical synthesis of natural products usually
            # competes with extraction/biotech (see Caffeine plausibility note).
            if method_ctx in ("biotechnological", "biotech"):
                return "fermentative Naturstoff-Route" if lang == "de" else "fermentative natural-product route"
            if method_ctx in ("chemical", "chem"):
                return "Extraktion oder Bio-Route meist günstiger" if lang == "de" else "extraction or bio route usually preferable"
            return "Naturstoff-Standardroute" if lang == "de" else "natural-product standard route"
        return "Standardroute, robuste Prozessparameter" if lang == "de" else "standard route, robust parameters"

    def _brief_cost(v: str) -> str:
        """Short cost rationale shown next to the cost level.

        Bug B2 fix: previously this method returned subtype-specific
        rationales (e.g. "Multistep-Synthese (5 Schritte)") regardless of
        the actual cost category, so reports like "Kosten NIEDRIG (2/10)
        — Multistep-Synthese" had an internal contradiction. Now we
        BRANCH ON THE CATEGORY FIRST and pick a rationale that is
        consistent with whether the cost is low / medium / high /
        very high.
        """
        v_norm = str(v).lower()
        is_low = v_norm in ("low", "niedrig")
        is_high = v_norm in ("high", "hoch", "very high", "very_high", "sehr hoch")

        # ---- LOW-cost rationales ----
        if is_low:
            if mtype_ctx == "small_molecule" and msub_ctx == "volatile":
                return ("Bulk-Fermentation, günstige Edukte, destillative Aufarbeitung"
                        if lang == "de"
                        else "bulk fermentation, cheap feedstock, distillation workup")
            if mtype_ctx == "protein" and msub_ctx == "enzyme":
                return ("Industrieenzym-Bulk, UF/DF-Standard-DSP"
                        if lang == "de"
                        else "industrial-enzyme bulk, UF/DF standard DSP")
            if mtype_ctx == "natural_product" and method_ctx in ("extraction", "extract"):
                return ("etablierte Extraktionsroute, Commodity-Pflanzenmaterial"
                        if lang == "de"
                        else "established extraction route, commodity plant material")
            return ("günstige Rohstoffe, Standard-Aufarbeitung, niedrige Komplexität"
                    if lang == "de"
                    else "cheap raw materials, standard workup, low complexity")

        # ---- HIGH / VERY HIGH cost rationales ----
        if is_high:
            if msub_ctx == "antibody":
                return ("Protein-A-Affinität + CHO-Zellkultur"
                        if lang == "de"
                        else "protein-A capture + CHO culture")
            if isinstance(eur_kg, (int, float)) and float(eur_kg) >= 500:
                return (f"Rohstoffe ≈ {float(eur_kg):.0f} €/kg dominiert"
                        if lang == "de"
                        else f"raw material ≈ {float(eur_kg):.0f} €/kg dominates")
            if mtype_ctx == "peptide" and method_ctx == "chemical" and isinstance(steps_ctx, int) and steps_ctx >= 15:
                return (f"SPPS-Synthese ({steps_ctx} Kupplungen), Reagenz-intensiv"
                        if lang == "de"
                        else f"SPPS synthesis ({steps_ctx} couplings), reagent-intensive")
            if mtype_ctx == "peptide" and method_ctx in ("biotechnological", "biotech"):
                return ("NRPS-Fermentation + Multi-Step-Chromatographie"
                        if lang == "de"
                        else "NRPS fermentation + multi-step chromatography")
            if mtype_ctx == "protein":
                return ("Fermentation + komplexe Chromatographie-Sequenz"
                        if lang == "de"
                        else "fermentation + complex chromatography sequence")
            if mtype_ctx == "natural_product" and method_ctx in ("extraction", "extract") \
                    and pur_diff == "very high":
                return ("komplexes Naturstoff-Gemisch, aufwändige Aufreinigung"
                        if lang == "de"
                        else "complex natural-product mixture, demanding purification")
            if isinstance(steps_ctx, int) and steps_ctx >= 5:
                return (f"Multistep-Synthese ({steps_ctx} Schritte), Reagenz- und Energie-Aufwand"
                        if lang == "de"
                        else f"multistep synthesis ({steps_ctx} steps), reagent and energy load")
            if raw_cost_ctx == "high":
                return ("teure Rohstoffe (high COGS)"
                        if lang == "de"
                        else "expensive raw materials (high COGS)")
            return ("mehrere Kostentreiber kombiniert"
                    if lang == "de"
                    else "multiple cost drivers combined")

        # ---- MEDIUM cost: no single dominator ----
        if mtype_ctx == "peptide" and method_ctx == "chemical" and isinstance(steps_ctx, int) and steps_ctx >= 15:
            return (f"SPPS-Synthese ({steps_ctx} Kupplungen), Rohstoffe moderat"
                    if lang == "de"
                    else f"SPPS synthesis ({steps_ctx} couplings), moderate feedstock cost")
        if isinstance(steps_ctx, int) and steps_ctx >= 4:
            return (f"Multistep-Prozess ({steps_ctx} Schritte), Standard-Edukte"
                    if lang == "de"
                    else f"multistep process ({steps_ctx} steps), standard feedstock")
        return ("moderate Treiber, keine Einzeldominanz"
                if lang == "de"
                else "moderate drivers, no single dominator")

    def _brief_risk(v: str) -> str:
        if msub_ctx == "antibody":
            return "Faltung, Glykosylierung, Sterilität" if lang == "de" else "folding, glycosylation, sterility"
        if mtype_ctx == "protein" and (n_dis >= 3 or has_ptm_ctx):
            return "komplexe Faltung + PTM-Kontrolle" if lang == "de" else "complex folding + PTM control"
        if isinstance(steps_ctx, int) and steps_ctx >= 6:
            return f"hohe Schrittzahl ({steps_ctx}) multipliziert Fehlerquellen" if lang == "de" else f"high step count ({steps_ctx}) compounds risk"
        if method_ctx in ("biotechnological", "biotech"):
            return "Kontaminations- und Sterilitäts-Anforderungen" if lang == "de" else "contamination / sterility requirements"
        if isinstance(purity_pct, (int, float)) and float(purity_pct) >= 99.0:
            return "sehr hohe Reinheitsziele (Out-of-Spec-Risiko)" if lang == "de" else "very high purity (out-of-spec risk)"
        if strict_waste_ctx and method_ctx in ("chemical", "chem"):
            return "strikte Abfallauflagen, halogenierte Solvents" if lang == "de" else "strict waste rules, halogenated solvents"
        if str(v).lower() in ("low", "niedrig"):
            return "beherrschbar, etablierte Parameter" if lang == "de" else "manageable, established parameters"
        return "moderate Prozess-Sensitivität" if lang == "de" else "moderate process sensitivity"

    def _brief_tox(v: str) -> str:
        # Controlled-substance check restricted to ACTUAL BtMG/DEA-scheduled
        # molecules. Earlier version flagged all alkaloids → wrongly flagged
        # caffeine, atropine, quinine etc.
        controlled = {
            "morphine", "codeine", "heroin", "cocaine", "oxycodone",
            "hydrocodone", "methadone", "buprenorphine", "fentanyl",
            "amphetamine", "methamphetamine",
        }
        mname = str(ip_ctx.get('molecule_name') or '').strip().lower()
        if mname in controlled:
            return "regulatorische Klassifizierung (BtMG-relevant)" if lang == "de" else "regulatory classification (controlled-substance)"

        if mtype_ctx == "protein":
            return "Endotoxin-Kontrolle, Glykosylierungs-Konsistenz" if lang == "de" else "endotoxin control, glycosylation consistency"
        if mtype_ctx == "peptide":
            return "Endotoxin-Spezifikation, Reinheit der Sequenz" if lang == "de" else "endotoxin specification, sequence purity"
        # Solvent containment only when the synthesis is non-trivial — for a
        # 1-step esterification (Aspirin) the reagent profile is mild and
        # "Solvent-Containment" reads alarmist.
        if (method_ctx in ("chemical", "chem") and strict_waste_ctx
                and isinstance(steps_ctx, int) and int(steps_ctx) >= 3):
            return "Solvent-Containment + Reagenz-Handhabung" if lang == "de" else "solvent containment + reagent handling"
        if method_ctx in ("chemical", "chem") and strict_waste_ctx:
            return "Standard-API-Handhabung, GMP-Reagenzqualität" if lang == "de" else "standard API handling, GMP reagent grade"
        if str(v).lower() in ("low", "niedrig"):
            return "Standard-Handhabung, keine kritischen Reagenzien" if lang == "de" else "standard handling, no critical reagents"
        return "moderate Containment-Anforderungen" if lang == "de" else "moderate containment requirements"

    eff = analysis.get('efficiency', 'n/a')
    cost_val = analysis.get('cost', 'n/a')
    risk_val = analysis.get('risk', 'n/a')
    tox = analysis.get('final_toxicity', analysis.get('toxicity', 'n/a'))

    # Render each metric as: <Label>: BEWERTUNG  (X / 10 — kurze Begründung)
    # New format (Bug B6): "{Label}: {score}/10 · {Tier-Name} — {rationale}"
    # The score stays prominent (concrete number to compare on), the tier
    # name replaces the unprofessional HOCH/MITTEL/NIEDRIG, and the
    # rationale stays in muted grey as before.
    _muted = "#6B7B8A"
    def _line(label_key: str, val: str, score: str, reason: str, metric: str) -> str:
        return (
            f"<b>{L(label_key)}:</b> {score} · {perf_value_markup(val, metric=metric)} "
            f'<font color="{_muted}" size="9">— {reason}</font>'
        )

    elems.append(Paragraph(_line("pdf_efficiency", eff, _score_for(eff, "efficiency"), _brief_eff(eff), "efficiency"), styles["Score"]))
    elems.append(Paragraph(_line("pdf_cost", cost_val, _score_for(cost_val, "cost"), _brief_cost(cost_val), "cost"), normal))
    elems.append(Paragraph(_line("pdf_risk", risk_val, _score_for(risk_val, "risk"), _brief_risk(risk_val), "risk"), normal))
    elems.append(Paragraph(_line("pdf_toxicity", tox, _score_for(tox, "toxicity"), _brief_tox(tox), "toxicity"), normal))

    # ---- COGS range estimate (€/kg) — directly below cost qualifier ----
    # Pharma reviewers expect a concrete economic anchor next to the
    # qualitative cost rating. We surface the range here in the exec
    # summary, and a fuller breakdown table further down the report.
    try:
        from cogs_estimator import estimate_cogs as _est_cogs, format_cogs_range as _fmt_cogs
        cogs_result = _est_cogs(ip_ctx) or {}
        cogs_range = _fmt_cogs(cogs_result, lang=lang)
    except Exception:
        cogs_result = None
        cogs_range = None

    if cogs_range:
        conf_word = {
            "high": ("hohe Confidence", "high confidence"),
            "medium": ("mittlere Confidence", "medium confidence"),
            "low": ("niedrige Confidence", "low confidence"),
        }.get(str((cogs_result or {}).get("confidence", "low")), ("low conf", "low conf"))
        conf_text = conf_word[0] if lang == "de" else conf_word[1]
        anchor_text = (cogs_result or {}).get("anchor_de" if lang == "de" else "anchor_en") or ""
        cogs_label = ("Geschätzte COGS" if lang == "de" else "Estimated COGS")
        elems.append(Paragraph(
            f'<b>{cogs_label}:</b> <font color="#1E3A5F"><b>{cogs_range}</b></font> '
            f'<font color="{_muted}" size="9">({conf_text} — {anchor_text})</font>',
            normal,
        ))

    elems.append(Spacer(1, 12))

    # Bug B6 fix: removed standalone "Confidence Level: Hoch/Mittel/Niedrig"
    # line. It was a derived mix-score of efficiency and cost that did NOT
    # measure recommendation reliability (despite the name suggesting so).
    # The COGS-confidence already lives in the cogs line above, so this
    # extra label was redundant and misleading.

    elems.append(Paragraph(L("pdf_key_issues"), styles["Section"]))
    # Engine-produced issues are in English; translate inline.
    # Bug N4 fix: surface only the TOP issues (1–3) here, not the full
    # list. Otherwise this section duplicates the later "Wichtige Risiken"
    # section verbatim (seen on Caffeine and other PDFs).
    issues = analysis.get('issues', []) or []
    if issues:
        for it in issues[:3]:
            prob, why = _explain_issue_pair(it)
            # Bug B8 fix: bullet for the main statement, indented italic
            # in muted grey for the explanation — replaces the
            # "?" lines (formerly "→" arrows sanitised out).
            elems.append(Paragraph(f"• {TE(prob)}", normal))
            elems.append(Paragraph(TE(why), styles["ExplanationIndent"]))
    else:
        elems.append(Paragraph(L("pdf_no_critical_issues"), normal))
    elems.append(Spacer(1, 10))

    elems.append(Paragraph(L("pdf_recommended_improvements"), styles["Section"]))
    imps = analysis.get('improvements', []) or []
    if imps:
        for imp in imps[:5]:
            action, impact, reason = _explain_improvement(imp)
            elems.append(Paragraph(f"• {TE(action)}", normal))
            elems.append(Paragraph(TE(impact), styles["ExplanationIndent"]))
            elems.append(Paragraph(TE(reason), styles["ExplanationIndent"]))
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

    # Bug F2 fix: a short "Why this route?" rationale that's always
    # rendered, independent of whether issues / actions were produced.
    # Previously the "Begründung"-Wort only appeared inside per-action
    # rationale blocks, so simple molecules (Ethanol, Aspirin, Lipase, ...)
    # had no answer to the basic question "why this method?" at all.
    try:
        _mol_type = (input_params.get("molecule_type") or "").lower()
        _mol_sub = (input_params.get("molecule_subtype") or "").lower()
        _method = (input_params.get("method") or "").lower()
        if lang == "en":
            _method_label = {
                "chemical": "chemical synthesis",
                "biotechnological": "biotechnological / fermentative production",
                "extraction": "extraction from a natural source",
            }.get(_method, _method or "the chosen production method")
            _route_explainers = {
                ("small_molecule", "volatile"):
                    "for volatile small molecules, fermentation followed by distillation is the dominant industrial route (energy-efficient, well-established)",
                ("small_molecule", "non_volatile"):
                    "for non-volatile small molecules with established commodity benchmarks (aspirin, paracetamol, ibuprofen class), multi-step chemical synthesis with crystallisation polishing is the cost-optimum",
                ("natural_product", "alkaloid"):
                    "for alkaloids the route choice depends on plant-source availability vs. synthetic accessibility — for commercially established alkaloids (caffeine, theophylline) chemical synthesis is competitive",
                ("natural_product", "terpene"):
                    "for terpenes, extraction from the natural source is standard at scale; engineered fermentation (Amyris/Evolva platforms) is the emerging alternative",
                ("peptide", "linear"):
                    "for short linear peptides (< 30 AAs), SPPS is the established industrial standard; longer peptides increasingly use recombinant expression with chemical acylation",
                ("peptide", "cyclic"):
                    "for cyclic peptides with non-proteinogenic residues (Cyclosporine class), NRPS fermentation in T. inflatum / Streptomyces is the only industrially viable route",
                ("protein", "antibody"):
                    "for mAbs the industry standard is CHO fed-batch culture with Protein-A affinity capture; the choice of host is dictated by the need for human-like glycosylation",
                ("protein", "enzyme"):
                    "for industrial enzymes (detergent/food/feed grade), submerged fermentation in established hosts (A. oryzae, B. subtilis, P. pastoris) with UF/DF polishing is the cost-optimum",
            }
            _why = _route_explainers.get((_mol_type, _mol_sub))
            _heading = "Route rationale"
        else:
            _method_label = {
                "chemical": "chemische Synthese",
                "biotechnological": "biotechnologische / fermentative Produktion",
                "extraction": "Extraktion aus natürlicher Quelle",
            }.get(_method, _method or "die gewählte Produktionsmethode")
            _route_explainers = {
                ("small_molecule", "volatile"):
                    "Bei flüchtigen Kleinmolekülen ist Fermentation mit anschließender Destillation der industrielle Standard (energieeffizient, gut etabliert).",
                ("small_molecule", "non_volatile"):
                    "Bei nicht-flüchtigen Kleinmolekülen mit etablierten Commodity-Benchmarks (Aspirin-, Paracetamol-, Ibuprofen-Klasse) ist die mehrstufige chemische Synthese mit Kristallisations-Polishing kostenoptimal.",
                ("natural_product", "alkaloid"):
                    "Bei Alkaloiden hängt die Routenwahl von Pflanzenquellen-Verfügbarkeit gegenüber synthetischer Zugänglichkeit ab — für kommerziell etablierte Alkaloide (Caffeine, Theophyllin) ist die chemische Synthese konkurrenzfähig.",
                ("natural_product", "terpene"):
                    "Bei Terpenen ist die Extraktion aus der natürlichen Quelle im großen Maßstab Standard; Engineered-Fermentation (Amyris/Evolva-Plattformen) ist die aufkommende Alternative.",
                ("peptide", "linear"):
                    "Bei kurzen linearen Peptiden (< 30 AS) ist SPPS der etablierte industrielle Standard; längere Peptide nutzen zunehmend rekombinante Expression mit chemischer Acylierung.",
                ("peptide", "cyclic"):
                    "Bei cyclischen Peptiden mit nicht-proteinogenen Resten (Cyclosporin-Klasse) ist die NRPS-Fermentation in T. inflatum / Streptomyces die einzige industriell brauchbare Route.",
                ("protein", "antibody"):
                    "Bei mAbs ist der Industriestandard CHO-Fed-Batch mit Protein-A-Affinitätscapture; die Wahl des Wirts wird durch die Notwendigkeit humanähnlicher Glykosylierung diktiert.",
                ("protein", "enzyme"):
                    "Bei industriellen Enzymen (Detergens-/Food-/Feed-Grade) ist die Submerged-Fermentation in etablierten Wirten (A. oryzae, B. subtilis, P. pastoris) mit UF/DF-Polishing kostenoptimal.",
            }
            _why = _route_explainers.get((_mol_type, _mol_sub))
            _heading = "Routenwahl-Begründung"
        if _why:
            elems.append(Paragraph(_heading, styles["Section"]))
            if lang == "en":
                elems.append(Paragraph(
                    f"The engine selected <b>{_method_label}</b> for this molecule: {_why}.",
                    normal,
                ))
            else:
                elems.append(Paragraph(
                    f"Die Engine hat <b>{_method_label}</b> für dieses Molekül gewählt: {_why}",
                    normal,
                ))
            elems.append(Spacer(1, 8))
    except Exception:
        pass

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
        # Bug C1 fix: locale-aware number formatting.
        # German uses point as thousands separator and comma as decimal;
        # English uses comma as thousands and point as decimal. Old code
        # mixed both inside the same German report (e.g. "100.000" for
        # 100, "5,000" for 5000) which made readers think the scale was
        # off by 1000x.
        _fmt_int = _format_int_en if lang == "en" else _format_int_de
        _fmt_money = _format_money_en if lang == "en" else _format_money_de
        if pct is not None:
            pct_str = f"{pct:.2f}".replace(".", ",") if lang != "en" else f"{pct:.2f}"
            elems.append(Paragraph(f"<b>{L('pdf_target_purity')}:</b> {pct_str} %", normal))
        if kg is not None:
            # Integers for whole-kg scales, fractional for sub-kg (lab)
            if float(kg) >= 1.0:
                kg_disp = _fmt_int(kg)
            else:
                # Sub-kg lab scale: show with up to 3 decimals
                kg_disp = (f"{float(kg):.3f}".replace(".", ",")
                           if lang != "en" else f"{float(kg):.3f}")
            kg_unit = "kg/year" if lang == "en" else "kg/Jahr"
            elems.append(Paragraph(f"<b>{L('pdf_target_scale')}:</b> {kg_disp} {kg_unit}", normal))
        if eur is not None:
            elems.append(Paragraph(f"<b>{L('pdf_raw_material_price')}:</b> {_fmt_money(eur)} €/kg", normal))
        if n_supp is not None:
            elems.append(Paragraph(f"<b>{L('pdf_qualified_suppliers')}:</b> {int(n_supp)}", normal))
        if lead is not None:
            weeks_unit = "weeks" if lang == "en" else "Wochen"
            lead_str = (f"{lead:.1f}".replace(".", ",") if lang != "en" else f"{lead:.1f}")
            elems.append(Paragraph(f"<b>{L('pdf_lead_time')}:</b> {lead_str} {weeks_unit}", normal))
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

    # ---- Process-flow block diagram + per-step legend ----
    # Two-row block diagram (Upstream → Downstream) gives reviewers a
    # one-glance view of the production chain. The legend immediately
    # below explains what each block means for THIS molecule class —
    # otherwise labels like "Stamm / Inokulum" stay opaque to non-experts.
    try:
        from process_flow import (
            build_process_flow_drawing as _build_flow,
            build_process_flow_legend as _build_flow_legend,
        )
        _concrete = extras.get("concrete_recs") if isinstance(extras, dict) else None
        _flow_drawing = _build_flow(input_params, _concrete, lang=lang)
        _flow_legend = _build_flow_legend(input_params, _concrete, lang=lang)
    except Exception:
        _flow_drawing = None
        _flow_legend = []

    if _flow_drawing is not None:
        flow_header = "Prozess-Flow" if lang == "de" else "Process flow"
        elems.append(Paragraph(flow_header, styles["Section"]))
        elems.append(_flow_drawing)
        flow_caption = ("Schematisch — Block-Reihenfolge spiegelt die "
                        "typische Sequenz für diese Molekül-Klasse wider. "
                        "Erklärung der einzelnen Schritte folgt unten."
                        if lang == "de" else
                        "Schematic — block order reflects the typical "
                        "sequence for this molecule class. Per-step "
                        "explanations follow below.")
        elems.append(Paragraph(f"<font size='8' color='#6B7B8A'>{flow_caption}</font>", normal))
        elems.append(Spacer(1, 8))

        # Per-step legend — addresses the "what is Stamm/Inokulum?" gap
        if _flow_legend:
            upstream_items = [x for x in _flow_legend if x.get("section") == "upstream"]
            downstream_items = [x for x in _flow_legend if x.get("section") == "downstream"]

            up_header_label = ("Upstream — Schritt-Erklärungen"
                               if lang == "de" else "Upstream — step explanations")
            dn_header_label = ("Downstream — Schritt-Erklärungen"
                               if lang == "de" else "Downstream — step explanations")

            if upstream_items:
                elems.append(Paragraph(
                    f"<font size='9' color='{'#1E3A5F'}'><b>{up_header_label}</b></font>",
                    normal,
                ))
                for i, item in enumerate(upstream_items, 1):
                    lbl = item.get("label", "")
                    expl = item.get("explanation", "")
                    elems.append(Paragraph(
                        f"<font size='9'><b>{i}. {lbl}</b> — {expl}</font>",
                        normal,
                    ))
                elems.append(Spacer(1, 4))

            if downstream_items:
                elems.append(Paragraph(
                    f"<font size='9' color='{'#3DB8C5'}'><b>{dn_header_label}</b></font>",
                    normal,
                ))
                for i, item in enumerate(downstream_items, 1):
                    lbl = item.get("label", "")
                    expl = item.get("explanation", "")
                    elems.append(Paragraph(
                        f"<font size='9'><b>{i}. {lbl}</b> — {expl}</font>",
                        normal,
                    ))
                elems.append(Spacer(1, 6))

        elems.append(Spacer(1, 10))

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
                    # Use ReportLab <super> markup so "Å²" renders correctly
                    # (the bare U+00B2 is not in Helvetica's WinAnsi range → ■).
                    elems.append(Paragraph(f"<b>{L('structure_label_tpsa')}:</b> {rdk_p['tpsa']:.2f} Å<super>2</super>", normal))
                    elems.append(Paragraph(f"<b>{L('structure_label_arom_rings')}:</b> {rdk_p['num_aromatic_rings']}", normal))
                    elems.append(Paragraph(f"<b>{L('structure_label_hbd')}:</b> {rdk_p['num_h_donors']}  ·  <b>{L('structure_label_hba')}:</b> {rdk_p['num_h_acceptors']}", normal))
                    elems.append(Paragraph(f"<b>{L('structure_label_rotbonds')}:</b> {rdk_p['num_rotatable_bonds']}", normal))
                    # Color-coded text marker for Lipinski pass/fail.
                    # ✓ renders in Helvetica but ⚠ does not (→ ■), so we
                    # use coloured bold text instead for consistency.
                    if rdk_p['lipinski_violations'] == 0:
                        elems.append(Paragraph(
                            f'<font color="#1E8449"><b>OK</b></font> &nbsp;{L("structure_label_lipinski_pass")}',
                            normal,
                        ))
                    else:
                        elems.append(Paragraph(
                            f'<font color="#C0392B"><b>!</b></font> &nbsp;{L("structure_label_lipinski_fail")} ({rdk_p["lipinski_violations"]})',
                            normal,
                        ))
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

    # Bug B5 fix: in addition to the generic input-driven contributors,
    # add SUBTYPE-SPECIFIC cost drivers so the report doesn't fall back
    # to the "Aufreinigung, Rohstoffe und operative Overheads" generic
    # string for almost every molecule.
    _subtype_drivers_de = {
        ("small_molecule", "volatile"):
            ["Destillations-Energie", "Fermentationsmedium"],
        ("small_molecule", "non_volatile"):
            ["Reagenz-Kaskade", "Kristallisations-Lösungsmittel"],
        ("natural_product", "alkaloid"):
            ["Mehrstufige Reagenzien", "Kristallisation/Umkristallisation"],
        ("natural_product", "terpene"):
            ["Lösungsmittel-Volumen (Extraktion)", "Pflanzenmaterial-Logistik"],
        ("peptide", "linear"):
            ["Fmoc-Aminosäuren", "Kupplungsreagenzien (HBTU/DIPEA)", "RP-HPLC-Lösungsmittel"],
        ("peptide", "cyclic"):
            ["Fermentationsmedium (NRPS-Wirt)", "Multi-Step-Chromatographie", "Cyclisierungs-Ausbeute"],
        ("protein", "antibody"):
            ["Protein-A-Resin", "CHO-Zellkultur-Medium", "Virusfiltration"],
        ("protein", "enzyme"):
            ["Bioreaktor-Medium", "UF/DF-Membranen", "Formulierungs-Stabilisatoren"],
    }
    _subtype_drivers_en = {
        ("small_molecule", "volatile"):
            ["distillation energy", "fermentation medium"],
        ("small_molecule", "non_volatile"):
            ["reagent cascade", "crystallisation solvents"],
        ("natural_product", "alkaloid"):
            ["multi-step reagents", "crystallisation / recrystallisation"],
        ("natural_product", "terpene"):
            ["solvent volume (extraction)", "plant-material logistics"],
        ("peptide", "linear"):
            ["Fmoc amino acids", "coupling reagents (HBTU/DIPEA)", "RP-HPLC solvents"],
        ("peptide", "cyclic"):
            ["fermentation medium (NRPS host)", "multi-step chromatography", "cyclisation yield"],
        ("protein", "antibody"):
            ["protein-A resin", "CHO cell-culture medium", "viral filtration"],
        ("protein", "enzyme"):
            ["bioreactor medium", "UF/DF membranes", "formulation stabilisers"],
    }
    _mt_key = (
        str(input_params.get("molecule_type") or "").lower(),
        str(input_params.get("molecule_subtype") or "").lower(),
    )
    _subtype_drivers = (_subtype_drivers_en if lang == "en" else _subtype_drivers_de)
    subtype_extras = _subtype_drivers.get(_mt_key, [])

    if contributors or subtype_extras:
        if lang == "en":
            translation = {
                "Rohstoffkosten": "raw-material cost",
                "Aufreinigung": "purification",
                "Abfallbehandlung": "waste handling",
            }
            loc_contribs = [translation.get(c, c) for c in contributors]
        else:
            loc_contribs = contributors
        # Merge subtype drivers in, dedup while preserving order.
        seen = set()
        merged = []
        for c in loc_contribs + subtype_extras:
            if c not in seen:
                seen.add(c)
                merged.append(c)
        prefix = L("pdf_main_cost_drivers")
        elems.append(Paragraph(prefix + ' ' + ', '.join(merged) + '.', normal))
    else:
        elems.append(Paragraph(L("pdf_main_cost_default"), normal))

    elems.append(Paragraph(L("pdf_cost_impact"), styles["Section"]))
    cost_level = (analysis.get('costLevel') or analysis.get('cost') or 'n/a')
    # Bug E5 fix: render the cost-level here with the same tier label
    # used in the executive summary (Commodity-Niveau / Mid-Tier-API /
    # Specialty-Tier / Biologika-Niveau) instead of leaking the uppercase
    # HOCH/MITTEL/NIEDRIG bucket word.
    elems.append(Paragraph(
        f"{L('pdf_cost_level')}: {perf_value_markup(str(cost_level), metric='cost')}",
        normal,
    ))

    # ---- COGS model: concrete €/kg range + breakdown table ----
    # Inserted between qualitative cost level and cost drivers so reviewers
    # can move from "wieviel kostet das eigentlich?" (range) to "warum?" (drivers)
    # in one continuous reading flow.
    try:
        from cogs_estimator import estimate_cogs as _est, format_cogs_range as _fmt
        _cogs = _est(input_params) or {}
        _rng = _fmt(_cogs, lang=lang)
    except Exception:
        _cogs = {}
        _rng = None

    if _rng:
        # Section header
        elems.append(Spacer(1, 6))
        cogs_header = "Kostenmodell (COGS-Schätzung)" if lang == "de" else "Cost model (COGS estimate)"
        elems.append(Paragraph(f"<b>{cogs_header}</b>", normal))

        # Range line — highlight in brand blue
        range_label = "Geschätzte Range" if lang == "de" else "Estimated range"
        elems.append(Paragraph(
            f"{range_label}: <font color='#1E3A5F'><b>{_rng}</b></font>",
            normal,
        ))

        # Confidence + anchor
        conf = _cogs.get("confidence", "low")
        conf_de = {"high": "Hoch", "medium": "Mittel", "low": "Niedrig"}.get(conf, conf)
        conf_en = {"high": "High", "medium": "Medium", "low": "Low"}.get(conf, conf)
        conf_text = conf_de if lang == "de" else conf_en
        anchor_text = _cogs.get("anchor_de" if lang == "de" else "anchor_en") or ""
        conf_lbl = "Confidence" if lang == "en" else "Confidence"
        elems.append(Paragraph(
            f"<font size='9'>{conf_lbl}: <b>{conf_text}</b> · {anchor_text}</font>",
            normal,
        ))

        # Breakdown table — typical share of total COGS
        breakdown = _cogs.get("breakdown") or {}
        if breakdown:
            from reportlab.platypus import Table, TableStyle
            from reportlab.lib import colors
            label_map_de = {
                "rm": "Rohstoffe",
                "dsp": "Aufarbeitung (DSP)",
                "opex": "Betriebskosten (OpEx)",
                "capex": "Anlagen-Amortisation (CapEx)",
            }
            label_map_en = {
                "rm": "Raw materials",
                "dsp": "Downstream (DSP)",
                "opex": "Operating expenses",
                "capex": "Asset amortisation (CapEx)",
            }
            labels = label_map_de if lang == "de" else label_map_en
            rows = [[
                "Posten" if lang == "de" else "Item",
                "Anteil" if lang == "de" else "Share",
            ]]
            for k in ("rm", "dsp", "opex", "capex"):
                if k in breakdown:
                    rows.append([labels[k], f"{int(round(float(breakdown[k]) * 100))} %"])
            tbl = Table(rows, colWidths=[100 * mm, 25 * mm])
            tbl.setStyle(TableStyle([
                ("FONT", (0, 0), (-1, 0), "Helvetica-Bold", 9),
                ("FONT", (0, 1), (-1, -1), "Helvetica", 9),
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.HexColor("#1E3A5F")),
                ("LINEBELOW", (0, 0), (-1, 0), 0.5, colors.HexColor("#3DB8C5")),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#F4F7FA")]),
                ("LEFTPADDING", (0, 0), (-1, -1), 6),
                ("RIGHTPADDING", (0, 0), (-1, -1), 6),
                ("TOPPADDING", (0, 0), (-1, -1), 3),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 3),
            ]))
            elems.append(Spacer(1, 4))
            elems.append(tbl)

        # Drivers list (user-facing, in target language)
        drivers_list = _cogs.get("drivers_de" if lang == "de" else "drivers_en") or []
        if drivers_list:
            drv_header = "Modifikatoren der Schätzung" if lang == "de" else "Estimate modifiers"
            elems.append(Spacer(1, 4))
            elems.append(Paragraph(f"<font size='9'><b>{drv_header}:</b></font>", normal))
            for dr in drivers_list[:6]:
                elems.append(Paragraph(f"<font size='9'>• {dr}</font>", normal))

        elems.append(Spacer(1, 8))

    drivers = analysis.get('costDrivers') or []
    if drivers:
        elems.append(Paragraph(L("pdf_main_cost_drivers"), normal))
        for d in drivers[:6]:
            elems.append(Paragraph(f"• {TE(d)}", normal))
    else:
        elems.append(Paragraph(L("pdf_no_cost_drivers"), normal))
    savings = analysis.get('savingsPotential') or 'low'
    # Bug E5 fix: render savings-potential as mixed-case German /
    # English instead of the uppercase token leak.
    _savings_de = {"low": "niedrig", "medium": "moderat", "high": "hoch", "very high": "sehr hoch"}
    _savings_en = {"low": "low", "medium": "moderate", "high": "high", "very high": "very high"}
    _savings_table = _savings_en if lang == "en" else _savings_de
    savings_label = _savings_table.get(str(savings).lower(), str(savings).lower())
    elems.append(Paragraph(
        f"{L('pdf_savings_potential')}: <b>{savings_label}</b>", normal,
    ))

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
        # Production / upstream steps (Feedback-Runde-3 — vorher fehlte das
        # komplett; User sah nur Downstream und wusste nicht, wie er das
        # Roh-Produkt überhaupt herstellt).
        prod_steps = concrete.get("production_steps") or []
        if prod_steps:
            prod_header = ("Produktions-Schritte (Upstream)" if lang == "de"
                           else "Production steps (upstream)")
            elems.append(Paragraph(prod_header, styles["Section"]))
            for i, s in enumerate(prod_steps, 1):
                elems.append(Paragraph(f"{i}. {_safe_val(s)}", normal))
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
        # Bug M3 fix: previously rendered four generic bullets ("Bewerten
        # Sie alternative Syntheserouten…", etc.). That undercut the value
        # of the modules that DID have concrete optimisations elsewhere in
        # the report. We now surface a single explicit "no recommendations
        # for this molecule yet" line that's clearly distinguishable from
        # the rule-driven block.
        elems.append(Paragraph(
            L("pdf_no_specific_optimizations"),
            normal,
        ))

    elems.append(Spacer(1, 12))
    # (B6: removed the second "Confidence Level: …" line — same reason as
    # on page 1: the value was a derived efficiency/cost mix-score that did
    # not actually measure recommendation reliability.)

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
            # Bug H2: sanitise Latin-Extended-A so author names like
            # 'Celińska' don't render as 'Celi■ska' under Helvetica.
            elems.append(Paragraph(f"[{idx}] {_to_pdf_safe(src)}", normal))
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
