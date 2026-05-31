"""
comparison_report.py
====================
Dedicated PDF report for the "Prozessvergleich" sidebar page.

Where the main recommendation PDF is focused on the single chosen
route, this report is structured around comparing all viable routes:

  1. Cover page with the molecule + input parameters
  2. Side-by-side comparison table (Cost, Risk, Efficiency, COGS,
     Yield, Steps) with the engine's ranking
  3. One detail mini-section per route (key COGS drivers, route
     anchor, source)
  4. Final analysis page with the engine's recommendation, a
     quantified rationale, and the trade-offs the user should
     consider before committing

Public entry:
    export_comparison_pdf(process_input, lang="de") -> bytes
"""

from __future__ import annotations

import io
from typing import Any, Dict, List, Optional


def _format_int_de(v: float) -> str:
    return f"{int(round(float(v))):,}".replace(",", ".")


def _format_int_en(v: float) -> str:
    return f"{int(round(float(v))):,}"


def _format_money_de(v: float) -> str:
    raw = f"{float(v):,.2f}"
    return raw.replace(",", "X").replace(".", ",").replace("X", ".")


def _format_cogs_value(x: Optional[float], lang: str) -> str:
    if x is None:
        return "—"
    if x >= 1000:
        s = f"{int(round(x)):,}"
        return s.replace(",", ".") if lang != "en" else s
    if x >= 10:
        return f"{int(round(x))}"
    s = f"{x:.1f}"
    return s.replace(".", ",") if lang != "en" else s


_NAVY = "#1E3A5F"
_ACCENT = "#3DB8C5"
_MUTED = "#5A6B7B"
_DIVIDER = "#DDE5EC"


def _make_cmp_styles():
    """Build a fresh stylesheet with the Cmp* paragraph styles.

    Self-contained so the comparison sections can be rendered both as a
    standalone PDF (export_comparison_pdf) and appended into the main
    recommendation report (report_generator.export_report_pdf) without any
    shared mutable state.
    """
    from reportlab.lib import colors
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle

    styles = getSampleStyleSheet()
    defs = [
        ("CmpTitle",    dict(parent=styles["Title"], fontSize=24, leading=28,
                             alignment=1, textColor=colors.HexColor(_NAVY), spaceAfter=8)),
        ("CmpSubtitle", dict(parent=styles["Normal"], fontSize=11, leading=14,
                             alignment=1, textColor=colors.HexColor(_MUTED), spaceAfter=12)),
        ("CmpSection",  dict(parent=styles["Heading2"], fontSize=14, leading=18,
                             textColor=colors.HexColor(_NAVY), spaceBefore=8, spaceAfter=6)),
        ("CmpBody",     dict(parent=styles["Normal"], fontSize=10, leading=14,
                             textColor=colors.HexColor("#1A1A1A"))),
        ("CmpMuted",    dict(parent=styles["Normal"], fontSize=9, leading=12,
                             textColor=colors.HexColor(_MUTED))),
        ("CmpBigStat",  dict(parent=styles["Normal"], fontSize=15, leading=18,
                             alignment=1, textColor=colors.HexColor(_NAVY))),
    ]
    for name, kw in defs:
        if name not in styles:
            styles.add(ParagraphStyle(name=name, **kw))
    return styles


def _tier_resolve(label, metric, lang):
    """Tier-label resolver — mirrors the executive-summary vocabulary."""
    tier_de = {
        "efficiency": {"low": "Kritisch", "medium": "Optimierbar",
                       "high": "Robust", "very high": "Best-in-class"},
        "cost":       {"low": "Commodity-Niveau", "medium": "Mid-Tier-API",
                       "high": "Specialty-Tier", "very high": "Biologika-Niveau"},
        "risk":       {"low": "Etabliert", "medium": "Akzeptabel",
                       "high": "Überwachungsbedürftig", "very high": "Kritisch"},
    }
    tier_en = {
        "efficiency": {"low": "Critical", "medium": "Optimisable",
                       "high": "Robust", "very high": "Best-in-class"},
        "cost":       {"low": "Commodity-tier", "medium": "Mid-tier API",
                       "high": "Specialty-tier", "very high": "Biologic-tier"},
        "risk":       {"low": "Established", "medium": "Acceptable",
                       "high": "Monitoring-required", "very high": "Critical"},
    }
    try:
        key = str(label or "medium").lower()
        key = {"niedrig": "low", "mittel": "medium", "hoch": "high",
               "sehr hoch": "very high"}.get(key, key)
        table = tier_en if lang == "en" else tier_de
        return table.get(metric, {}).get(key, str(label).title())
    except Exception:
        return str(label).title()


def build_comparison_flowables(process_input: Dict[str, Any],
                               lang: str = "de",
                               styles=None,
                               heading: Optional[str] = None) -> List[Any]:
    """Render the comparison body — table + per-route detail + final
    analysis — as a list of flowables, WITHOUT a cover page or document
    setup.

    Reusable both by the standalone comparison PDF and as an appended
    section inside the main recommendation report. Pass `heading` to emit a
    top-level section title (used when appended to the recommendation PDF).
    """
    from reportlab.lib.units import mm
    from reportlab.lib import colors
    from reportlab.platypus import Paragraph, Spacer, PageBreak, Table, TableStyle

    from report_generator import _to_pdf_safe
    from route_comparison import compare_routes, METHODS_DE, METHODS_EN

    if styles is None:
        styles = _make_cmp_styles()
    methods_table = METHODS_EN if lang == "en" else METHODS_DE

    comparison = compare_routes(process_input)
    mname = str(process_input.get("molecule_name") or "Molekül").strip()

    elems: List[Any] = []

    # When a heading is supplied (appended-into-recommendation mode) it
    # doubles as the overview title — avoids two near-identical headings
    # ("Routenvergleich" + "Vergleichsübersicht") stacked on top of each other.
    _overview_title = heading or ("Vergleichsübersicht" if lang != "en" else "Comparison Overview")

    # ----------- COMPARISON TABLE -----------
    elems.append(Paragraph(_overview_title, styles["CmpSection"]))
    elems.append(Spacer(1, 4))
    elems.append(Paragraph(
        ("Die Engine bewertet jede plausible Route nach Cost, Risk, "
         "Effizienz, COGS und Konkretheit der hinterlegten Schritte. "
         "Die empfohlene Route ist mit ★ markiert."
         if lang != "en" else
         "The engine scores every plausible route on Cost, Risk, "
         "Efficiency, COGS, and step concreteness. The recommended "
         "route is marked with ★."),
        styles["CmpMuted"],
    ))
    elems.append(Spacer(1, 12))

    if not comparison or not comparison.rows:
        elems.append(Paragraph(
            ("Keine vergleichbaren Routen für dieses Molekül identifiziert."
             if lang != "en" else
             "No comparable routes identified for this molecule."),
            styles["CmpBody"],
        ))
    else:
        header_de = ["", "Methode", "Cost", "Risk", "Effizienz", "COGS €/kg", "Yield", "Schritte"]
        header_en = ["", "Method", "Cost", "Risk", "Efficiency", "COGS €/kg", "Yield", "Steps"]
        header = header_en if lang == "en" else header_de

        # Text cells (method + tier labels) are rendered as Paragraphs so long
        # labels like "Commodity-Niveau" wrap within their column instead of
        # overrunning into the neighbouring cell. Numeric cells stay strings.
        from reportlab.lib.styles import ParagraphStyle as _PS
        _cell = _PS(name="CmpCell", parent=styles["CmpBody"], fontSize=8, leading=9.5)
        _cell_b = _PS(name="CmpCellB", parent=_cell, fontName="Helvetica-Bold")

        rows_data = [header]
        for row in comparison.rows:
            mark = "★" if row.decision_rank == 1 else f"{row.decision_rank}."
            method_lbl = row.method_label_en if lang == "en" else row.method_label_de
            _is_top = row.decision_rank == 1
            cogs_cell = (f"{_format_cogs_value(row.cogs_low_eur_per_kg, lang)}–"
                         f"{_format_cogs_value(row.cogs_high_eur_per_kg, lang)}"
                         if row.cogs_low_eur_per_kg is not None else "—")
            yield_cell = row.expected_yield or "—"
            steps_cell = (f"{row.n_production_steps}" if row.has_concrete_steps
                          else ("generic" if lang == "en" else "generisch"))
            rows_data.append([
                mark,
                Paragraph(method_lbl, _cell_b if _is_top else _cell),
                Paragraph(_tier_resolve(row.cost_label, "cost", lang), _cell),
                Paragraph(_tier_resolve(row.risk_label, "risk", lang), _cell),
                Paragraph(_tier_resolve(row.efficiency_label, "efficiency", lang), _cell),
                cogs_cell, yield_cell, steps_cell,
            ])

        # Total kept a few mm under the 174 mm usable width (A4 − 2×18 mm)
        # so the right-most "Schritte" column never clips. Cost gets 27 mm so
        # "Commodity-Niveau" / "Biologika-Niveau" fit on one line; Yield gets
        # 18 mm for the g/L titer strings ("0.3–1.5 g/L").
        col_widths = [
            9 * mm, 36 * mm, 29 * mm, 20 * mm, 20 * mm, 21 * mm, 18 * mm, 16 * mm,
        ]
        table = Table(rows_data, colWidths=col_widths, hAlign="LEFT", repeatRows=1)
        table.setStyle(TableStyle([
            ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor(_NAVY)),
            ("TEXTCOLOR",  (0, 0), (-1, 0), colors.white),
            ("FONTNAME",   (0, 0), (-1, 0), "Helvetica-Bold"),
            ("FONTSIZE",   (0, 0), (-1, 0), 8.5),
            ("BOTTOMPADDING", (0, 0), (-1, 0), 6),
            ("TOPPADDING",    (0, 0), (-1, 0), 6),
            ("FONTSIZE", (0, 1), (-1, -1), 8.5),
            ("VALIGN",   (0, 0), (-1, -1), "MIDDLE"),
            ("ALIGN",    (0, 0), (0, -1), "CENTER"),
            ("ALIGN",    (5, 1), (5, -1), "RIGHT"),
            ("ALIGN",    (6, 1), (-1, -1), "CENTER"),
            ("BOTTOMPADDING", (0, 1), (-1, -1), 5),
            ("TOPPADDING",    (0, 1), (-1, -1), 5),
            ("BACKGROUND", (0, 1), (-1, 1), colors.HexColor("#E8F0F7")),
            ("FONTNAME",   (1, 1), (1, 1), "Helvetica-Bold"),
            ("LINEBELOW", (0, 0), (-1, 0), 0.5, colors.HexColor("#B0BEC5")),
            ("LINEBELOW", (0, 1), (-1, -2), 0.25, colors.HexColor(_DIVIDER)),
        ]))
        elems.append(table)
        elems.append(Spacer(1, 12))

        elems.append(Paragraph(
            ("<i><font size='8' color='" + _MUTED + "'>★ = Engine-Empfehlung · "
             "„generisch\" = Cluster-Schätzung ohne hinterlegte Schritte · "
             "Tier-Begriffe identisch mit dem Empfehlungs-Report.</font></i>"
             if lang != "en" else
             "<i><font size='8' color='" + _MUTED + "'>★ = engine recommendation · "
             "\"generic\" = cluster estimate without stored steps · "
             "tier terms identical to the recommendation report.</font></i>"),
            styles["CmpBody"],
        ))
    elems.append(PageBreak())

    # ----------- PER-ROUTE DETAIL -----------
    if comparison and comparison.rows:
        for row in comparison.rows:
            elems.append(Paragraph(
                (f"{'★ ' if row.decision_rank == 1 else ''}Route #{row.decision_rank}: "
                 + (row.method_label_en if lang == "en" else row.method_label_de)),
                styles["CmpSection"],
            ))
            elems.append(Spacer(1, 4))

            facts_label_de = ["Cost-Niveau", "Risk-Niveau", "Effizienz", "COGS-Bereich",
                              "Erwartete Ausbeute", "Schritte", "Konfidenz"]
            facts_label_en = ["Cost level", "Risk level", "Efficiency", "COGS range",
                              "Expected yield", "Steps", "Confidence"]
            facts_label = facts_label_en if lang == "en" else facts_label_de

            cogs_low_str = _format_cogs_value(row.cogs_low_eur_per_kg, lang)
            cogs_high_str = _format_cogs_value(row.cogs_high_eur_per_kg, lang)
            cogs_text = (f"{cogs_low_str}–{cogs_high_str} €/kg"
                         if row.cogs_low_eur_per_kg is not None else "—")

            facts_values = [
                f"<b>{_tier_resolve(row.cost_label, 'cost', lang)}</b> ({row.cost_score}/10)",
                f"<b>{_tier_resolve(row.risk_label, 'risk', lang)}</b> ({row.risk_score}/10)",
                f"<b>{_tier_resolve(row.efficiency_label, 'efficiency', lang)}</b> ({row.efficiency_score}/10)",
                f"<b>{cogs_text}</b>",
                row.expected_yield or "—",
                (f"{row.n_production_steps} {'concrete steps' if lang == 'en' else 'konkrete Schritte'}"
                 if row.has_concrete_steps else
                 ("generic cluster estimate" if lang == "en" else "generische Cluster-Schätzung")),
                row.cogs_confidence or "—",
            ]

            kv_rows = []
            for lbl, val in zip(facts_label, facts_values):
                kv_rows.append([
                    Paragraph(f"<font color='{_MUTED}'>{lbl}</font>", styles["CmpBody"]),
                    Paragraph(str(val), styles["CmpBody"]),
                ])
            kv_table = Table(kv_rows, colWidths=[55 * mm, 119 * mm], hAlign="LEFT")
            kv_table.setStyle(TableStyle([
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 4),
                ("TOPPADDING",    (0, 0), (-1, -1), 4),
                ("LINEBELOW",     (0, 0), (-1, -2), 0.25, colors.HexColor(_DIVIDER)),
            ]))
            elems.append(kv_table)
            elems.append(Spacer(1, 8))

            _anchor = (row.cogs_anchor_en if lang == "en" else row.cogs_anchor_de) or ""
            if _anchor:
                elems.append(Paragraph(
                    f"<i><font size='8' color='{_MUTED}'>"
                    f"{'COGS anchor: ' if lang == 'en' else 'COGS-Anker: '}"
                    f"{_to_pdf_safe(_anchor)}</font></i>",
                    styles["CmpBody"],
                ))
            elems.append(Spacer(1, 16))

    elems.append(PageBreak())

    # ----------- FINAL ANALYSIS -----------
    elems.append(Paragraph(
        "Engine-Analyse und Empfehlung" if lang != "en" else "Engine analysis and recommendation",
        styles["CmpSection"],
    ))
    elems.append(Spacer(1, 6))

    if comparison and comparison.recommended_method:
        top_method_lbl = methods_table.get(comparison.recommended_method, comparison.recommended_method)
        elems.append(Paragraph(
            (f"Die Engine empfiehlt für <b>{_to_pdf_safe(mname)}</b> die folgende Produktionsroute:"
             if lang != "en" else
             f"For <b>{_to_pdf_safe(mname)}</b>, the engine recommends the following production route:"),
            styles["CmpBody"],
        ))
        elems.append(Spacer(1, 8))
        elems.append(Paragraph(
            f"<para alignment='center'><font size='15' color='{_NAVY}'><b>{top_method_lbl}</b></font></para>",
            styles["CmpBody"],
        ))
        elems.append(Spacer(1, 16))

        if len(comparison.rows) >= 2:
            top = comparison.rows[0]
            nxt = comparison.rows[1]
            cost_diff = nxt.cost_score - top.cost_score
            risk_diff = nxt.risk_score - top.risk_score
            eff_diff = top.efficiency_score - nxt.efficiency_score

            rationale_bits_de: List[str] = []
            rationale_bits_en: List[str] = []
            if cost_diff >= 2:
                rationale_bits_de.append(
                    f"Cost-Vorteil: {cost_diff} Bänder günstiger als {nxt.method_label_de}")
                rationale_bits_en.append(
                    f"Cost advantage: {cost_diff} bands cheaper than {nxt.method_label_en}")
            if risk_diff >= 2:
                rationale_bits_de.append(
                    f"Risk-Vorteil: {risk_diff} Bänder geringeres Operations-Risiko")
                rationale_bits_en.append(
                    f"Risk advantage: {risk_diff} bands lower operational risk")
            if eff_diff >= 2:
                rationale_bits_de.append(f"Effizienz-Vorteil: {eff_diff} Bänder robuster")
                rationale_bits_en.append(f"Efficiency advantage: {eff_diff} bands more robust")
            if top.has_concrete_steps and not nxt.has_concrete_steps:
                rationale_bits_de.append(
                    "Konkrete industrielle Vorlagen sind hinterlegt; die Alternative ist nur als Cluster-Schätzung verfügbar")
                rationale_bits_en.append(
                    "Concrete industrial precedents are available; the alternative is only a cluster estimate")

            if rationale_bits_de:
                elems.append(Paragraph(
                    "Quantifizierte Begründung:" if lang != "en" else "Quantified rationale:",
                    styles["CmpBody"],
                ))
                bits = rationale_bits_en if lang == "en" else rationale_bits_de
                for b in bits:
                    elems.append(Paragraph(f"• {b}", styles["CmpBody"]))
                elems.append(Spacer(1, 10))

        reason = (comparison.recommendation_reason_en if lang == "en"
                  else comparison.recommendation_reason_de)
        if reason:
            elems.append(Paragraph(
                "Engine-Begründung:" if lang != "en" else "Engine reasoning:",
                styles["CmpBody"],
            ))
            elems.append(Paragraph(
                f"<i><font color='{_MUTED}'>{_to_pdf_safe(reason)}</font></i>",
                styles["CmpBody"],
            ))
            elems.append(Spacer(1, 14))

        elems.append(Paragraph(
            "Zu beachten vor der Entscheidung:" if lang != "en" else
            "Trade-offs to consider before committing:",
            styles["CmpBody"],
        ))
        tradeoff_text_de = (
            "• Die COGS-Schätzungen sind Cluster-Anker mit ±30 % Unsicherheit — vor "
            "Investitions-Entscheidungen durch Pilot-Daten validieren.\n"
            "• Generische Cluster-Schätzungen (Spalte „Schritte: generisch\") sollten "
            "durch eine Literaturrecherche oder einen CDMO-Quote ergänzt werden.\n"
            "• Time-to-Pilot-Schätzungen gelten für etablierte Verfahren; neuartige "
            "Stammoptimierung kann 6–18 Monate zusätzlich kosten."
        )
        tradeoff_text_en = (
            "• COGS estimates are cluster-anchored with ±30 % uncertainty — validate "
            "with pilot data before investment decisions.\n"
            "• Generic cluster estimates (column \"Steps: generic\") should be "
            "supplemented by a literature search or a CDMO quote.\n"
            "• Time-to-pilot estimates assume an established procedure; novel "
            "strain engineering can add 6–18 months."
        )
        for line in (tradeoff_text_en if lang == "en" else tradeoff_text_de).split("\n"):
            if line.strip():
                elems.append(Paragraph(line.strip(), styles["CmpBody"]))
        elems.append(Spacer(1, 14))

    return elems


def export_comparison_pdf(process_input: Dict[str, Any],
                          lang: str = "de") -> bytes:
    """Generate the process-comparison PDF as bytes."""
    from reportlab.lib.pagesizes import A4
    from reportlab.lib.units import mm
    from reportlab.lib import colors
    from reportlab.platypus import (
        SimpleDocTemplate, Paragraph, Spacer, PageBreak, Table, TableStyle,
    )
    from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
    from reportlab.pdfgen.canvas import Canvas

    # Reuse the same unicode sanitizer + locale formatters from the main
    # report so the look-and-feel is consistent.
    from report_generator import _to_pdf_safe
    from route_comparison import compare_routes, METHODS_DE, METHODS_EN
    from cogs_estimator import estimate_cogs, format_cogs_range
    from i18n import t as _t

    _NAVY = "#1E3A5F"
    _ACCENT = "#3DB8C5"
    _MUTED = "#5A6B7B"
    _DIVIDER = "#DDE5EC"

    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(
        name="CmpTitle", parent=styles["Title"], fontSize=24, leading=28,
        alignment=1, textColor=colors.HexColor(_NAVY), spaceAfter=8,
    ))
    styles.add(ParagraphStyle(
        name="CmpSubtitle", parent=styles["Normal"], fontSize=11, leading=14,
        alignment=1, textColor=colors.HexColor(_MUTED), spaceAfter=12,
    ))
    styles.add(ParagraphStyle(
        name="CmpSection", parent=styles["Heading2"], fontSize=14, leading=18,
        textColor=colors.HexColor(_NAVY), spaceBefore=8, spaceAfter=6,
    ))
    styles.add(ParagraphStyle(
        name="CmpBody", parent=styles["Normal"], fontSize=10, leading=14,
        textColor=colors.HexColor("#1A1A1A"),
    ))
    styles.add(ParagraphStyle(
        name="CmpMuted", parent=styles["Normal"], fontSize=9, leading=12,
        textColor=colors.HexColor(_MUTED),
    ))
    styles.add(ParagraphStyle(
        name="CmpBigStat", parent=styles["Normal"], fontSize=15, leading=18,
        alignment=1, textColor=colors.HexColor(_NAVY),
    ))

    # Compute the comparison
    comparison = compare_routes(process_input)
    mname = str(process_input.get("molecule_name") or "Molekül").strip()
    mtype = str(process_input.get("molecule_type") or "").lower()
    msub = str(process_input.get("molecule_subtype") or "").lower()

    elems = []
    methods_table = METHODS_EN if lang == "en" else METHODS_DE

    # ----------- COVER PAGE -----------
    elems.append(Spacer(1, 60))
    elems.append(Paragraph("Helixar AI", styles["CmpTitle"]))
    elems.append(Paragraph(
        "Prozessvergleich" if lang != "en" else "Process Comparison",
        styles["CmpSubtitle"],
    ))
    elems.append(Spacer(1, 40))
    elems.append(Paragraph(
        f"<para alignment='center'><font size='18' color='{_NAVY}'><b>{_to_pdf_safe(mname)}</b></font></para>",
        styles["CmpBody"],
    ))
    elems.append(Spacer(1, 16))
    # Input parameters mini-block
    _kg = process_input.get("scale_kg_per_year")
    _pct = process_input.get("desired_purity_percent")
    _eur = process_input.get("raw_material_cost_eur_per_kg")
    fmt_int = _format_int_en if lang == "en" else _format_int_de
    fmt_money = (lambda v: f"{v:,.2f}") if lang == "en" else _format_money_de
    kg_unit = "kg/year" if lang == "en" else "kg/Jahr"
    pct_lbl = "Target purity" if lang == "en" else "Zielreinheit"
    scale_lbl = "Target scale" if lang == "en" else "Zielmaßstab"
    rm_lbl = "Raw-material price" if lang == "en" else "Marktpreis Edukte"

    input_lines = []
    if _pct is not None:
        pct_str = f"{_pct:.2f}".replace(".", ",") if lang != "en" else f"{_pct:.2f}"
        input_lines.append(f"<b>{pct_lbl}:</b> {pct_str} %")
    if _kg is not None:
        kg_str = fmt_int(_kg) if float(_kg) >= 1.0 else (f"{float(_kg):.3f}".replace(".", ",") if lang != "en" else f"{float(_kg):.3f}")
        input_lines.append(f"<b>{scale_lbl}:</b> {kg_str} {kg_unit}")
    if _eur is not None:
        input_lines.append(f"<b>{rm_lbl}:</b> {fmt_money(_eur)} €/kg")
    for line in input_lines:
        elems.append(Paragraph(
            f"<para alignment='center'>{line}</para>",
            styles["CmpBody"],
        ))
    elems.append(Spacer(1, 50))
    elems.append(Paragraph(
        f"<para alignment='center'><i><font color='{_MUTED}' size='9'>"
        f"{'This document compares all production routes the engine considers viable for ' if lang == 'en' else 'Dieses Dokument vergleicht alle Produktionsrouten, die die Engine für '}"
        f"{_to_pdf_safe(mname)}{'.' if lang == 'en' else ' als realistisch einstuft.'}"
        f"</font></i></para>",
        styles["CmpBody"],
    ))
    elems.append(PageBreak())

    # ----------- COMPARISON BODY (table + per-route detail + analysis) -----------
    # Reuses the shared renderer so the standalone PDF and the appended
    # recommendation-report section stay byte-for-byte identical.
    elems.extend(build_comparison_flowables(process_input, lang=lang, styles=styles))

    # Footer
    elems.append(Spacer(1, 30))
    _footer_msg = (
        'For the detailed single-route recommendation: sidebar → "Empfehlung generieren"'
        if lang == "en" else
        'Für den detaillierten Einzelrouten-Report: Sidebar → „Empfehlung generieren\"'
    )
    elems.append(Paragraph(
        f"<para alignment='center'><i><font size='8' color='{_MUTED}'>"
        f"{_footer_msg}</font></i></para>",
        styles["CmpBody"],
    ))

    # Build PDF
    buf = io.BytesIO()

    def _draw_header_footer(canvas: Canvas, doc):
        canvas.saveState()
        canvas.setFont("Helvetica-Bold", 11)
        canvas.setFillColor(colors.HexColor(_NAVY))
        canvas.drawString(18 * mm, A4[1] - 12 * mm, "Helixar")
        width_hx = canvas.stringWidth("Helixar", "Helvetica-Bold", 11)
        canvas.setFillColor(colors.HexColor(_ACCENT))
        canvas.drawString(18 * mm + width_hx + 2, A4[1] - 12 * mm, "AI")
        canvas.setFont("Helvetica", 8)
        canvas.setFillColor(colors.HexColor(_MUTED))
        canvas.drawString(
            18 * mm,
            A4[1] - 12 * mm - 9,
            "Prozessvergleich" if lang != "en" else "Process Comparison",
        )
        # Page number
        page_num = canvas.getPageNumber()
        canvas.drawRightString(
            A4[0] - 18 * mm, 12 * mm,
            f"{'Page' if lang == 'en' else 'Seite'} {page_num}",
        )
        canvas.restoreState()

    doc = SimpleDocTemplate(
        buf, pagesize=A4,
        leftMargin=18 * mm, rightMargin=18 * mm,
        topMargin=22 * mm, bottomMargin=18 * mm,
        title=f"Helixar AI — Prozessvergleich {mname}",
    )
    doc.build(elems, onFirstPage=_draw_header_footer, onLaterPages=_draw_header_footer)
    return buf.getvalue()
