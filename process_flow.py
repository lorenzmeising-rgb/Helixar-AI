"""
process_flow.py
===============
Block-diagram visualisation of the production process for the PDF report,
plus a per-step explanation legend.

The earlier version truncated long downstream labels and gave reviewers
unexplained jargon ("Stamm / Inokulum"). This rewrite:

- Limits blocks to 4 per row (with auto-wrap to a second row) → wider boxes
- Wraps block labels onto up to 2 lines instead of truncating
- Emits a structured legend (block name → one-sentence explanation) that
  the PDF caller renders directly below the diagram

Public surface:
    build_process_flow_drawing(process_input, concrete_recs, lang) -> Drawing | None
    build_process_flow_legend(process_input, concrete_recs, lang) -> List[Tuple[str, str]]
        Returns [(block_label, explanation), ...] in display order:
        upstream first, then downstream.

Design constraints (from project memory):
- Pure ReportLab Shapes — no new dependencies.
- Brand-aligned colors (marine blue + turquoise).
- Honest representation: if we don't have a specific explanation for a
  block, return a generic one rather than fabricating detail.
"""

from typing import Dict, Any, List, Optional, Tuple

# Brand palette (mirrors .streamlit/config.toml)
BRAND_PRIMARY = "#1E3A5F"
BRAND_ACCENT = "#3DB8C5"
BRAND_NEUTRAL = "#6B7B8A"


# ---------------------------------------------------------------------------
# Block-label → explanation mapping
# ---------------------------------------------------------------------------
# Keys are normalised lowercase labels. We do a substring match: an explanation
# fires if its key appears in the block label (case-insensitive). The first
# match wins, so order more-specific keys before more-generic ones.

_BLOCK_EXPLANATIONS_DE: List[Tuple[str, str]] = [
    # --- Upstream (biotech) ---
    ("stamm / inokulum",
     "Produktionsorganismus aus der Stammbank — definiertes Klon-Backup, mehrfach getestet"),
    ("inokulum",
     "Anzucht des Produktionsstamms aus der Stammbank"),
    ("vorkultur",
     "Biomassenaufbau im Shake-Flask, typ. 16–24 h Wachstum vor dem Hauptreaktor"),
    ("hauptfermentation",
     "Fed-Batch im Bioreaktor, Produktbildung unter Prozesskontrolle (pH, O2, T)"),
    ("fermentation",
     "Mikrobielle Produktion im Bioreaktor, Produkt sammelt sich extrazellulär oder im Zellpellet"),
    # --- Upstream (extraction) ---
    ("rohmaterial",
     "Pflanzliches/biologisches Ausgangsmaterial — Qualität bestimmt Ausbeute und Verunreinigungs-Profil"),
    ("vorbehandlung",
     "Mechanische und thermische Aufbereitung (Mahlung, Trocknung) für effiziente Extraktion"),
    # --- Upstream (chemical) ---
    ("schritt 1",
     "Erste Synthese-Stufe — Einführung des Kerngerüsts oder Schutzgruppen-Strategie"),
    ("schritt 2",
     "Funktionalisierung — Substitutionen, Kupplungen oder Stereoeinführung"),
    ("schritt 3",
     "Finale Transformation und Aufreinigung zum Zielmolekül"),
    ("synthese",
     "Chemische Hauptsynthese, zusammenfassender Schritt"),
    # --- Downstream (mAb) ---
    ("protein-a-capture",
     "Affinitäts-Chromatographie — > 95 % Wirkstoff-Selektivität in einem Schritt"),
    ("protein-a capture",
     "Affinitäts-Chromatographie — > 95 % Wirkstoff-Selektivität in einem Schritt"),
    ("low-ph-virusinakt",
     "Inaktivierung viraler Kontaminanten bei pH 3,5 für 60–90 min — GMP-Standard für Säuger-Zellkulturen"),
    ("low-ph virus inact",
     "Inaktivierung viraler Kontaminanten bei pH 3,5 für 60–90 min — GMP-Standard für Säuger-Zellkulturen"),
    ("virusinaktivierung",
     "Inaktivierung viraler Kontaminanten bei pH 3,5 für 60–90 min — GMP-Standard für Säuger-Zellkulturen"),
    ("niedrig-ph",
     "Inaktivierung viraler Kontaminanten bei niedrigem pH — GMP-Standard für Säuger-Zellkulturen"),
    ("kation-iex",
     "Kationenaustausch-Chromatographie — entfernt Aggregate und Wirts-Proteine über Ladungsunterschiede"),
    ("cation iex",
     "Kationenaustausch-Chromatographie — entfernt Aggregate und Wirts-Proteine über Ladungsunterschiede"),
    ("kationenaustausch",
     "Kationenaustausch-Chromatographie — entfernt Aggregate und Wirts-Proteine über Ladungsunterschiede"),
    ("anion-iex",
     "Anionenaustausch + Virusfiltration — DNA/Endotoxin-Entfernung als Polishing-Schritt"),
    ("anion iex",
     "Anionenaustausch + Virusfiltration — DNA/Endotoxin-Entfernung als Polishing-Schritt"),
    ("anionenaustausch",
     "Anionenaustausch-Chromatographie — DNA/Endotoxin-Entfernung als Polishing-Schritt"),
    ("tff + formulation",
     "Tangential-Flow-Filtration — Pufferaustausch und Konzentration auf Zielspezifikation"),
    ("tff",
     "Tangential-Flow-Filtration — Pufferaustausch und Konzentration"),
    # --- Downstream (Protein general) ---
    # NOTE: "klärung" / "klärfiltration" come BEFORE the granular keys
    # (sterilfiltration, tiefenfiltration, virusfiltration) so that a step
    # labelled "Klärung via Tiefenfiltration + 0,2 µm Sterilfiltration"
    # matches the holistic concept rather than a sub-component.
    ("zell-lyse",
     "Mechanischer oder chemischer Zellaufschluss zur Freisetzung intrazellulärer Proteine"),
    ("klärfiltration",
     "Entfernung von Zelltrümmern via Tiefenfiltration oder Zentrifugation"),
    ("klärung",
     "Erste Trennung von Zelltrümmern und Medium-Bestandteilen (Tiefenfiltration + Sterilfiltration)"),
    ("virusfiltration",
     "Mechanische Virusfiltration — finale Sicherheitsstufe vor der Formulierung"),
    ("sterilfiltration",
     "0,2 µm Sterilfiltration — entfernt bakterielle Kontamination vor weiteren Schritten"),
    ("tiefenfiltration",
     "Tiefenfiltration — Klärung des Mediums durch Entfernung feiner Zelltrümmer"),
    ("affinitäts",
     "Affinitäts-Chromatographie — selektive Bindung über spezifischen Tag oder natürlichen Liganden"),
    ("sec-polishing",
     "Größenausschluss-Chromatographie — Aggregat- und Monomeren-Trennung als Polishing"),
    ("sec polish",
     "Größenausschluss-Chromatographie — Polishing-Schritt"),
    ("uf/df",
     "Ultrafiltration/Diafiltration — Pufferaustausch und Konzentration auf Endspezifikation"),
    # --- Downstream (Peptide) ---
    ("spps-cleavage",
     "Abspaltung vom Synthese-Harz mit TFA-Cocktail, gleichzeitig Schutzgruppenentfernung"),
    ("spps cleavage",
     "Abspaltung vom Synthese-Harz mit TFA-Cocktail, gleichzeitig Schutzgruppenentfernung"),
    ("prep-hplc (rp)",
     "Präparative Reversed-Phase-HPLC — Hauptaufreinigung auf > 95 % Reinheit"),
    ("prep-hplc",
     "Präparative HPLC — hochauflösende Aufreinigung"),
    ("counter-ion-exchange",
     "Austausch von TFA gegen pharmakologisch akzeptables Gegenion (typ. Acetat)"),
    ("counter-ion exchange",
     "Austausch von TFA gegen pharmakologisch akzeptables Gegenion (typ. Acetat)"),
    ("lyophilisation",
     "Gefriertrocknung — entzieht Wasser für stabile Lager-Formulierung"),
    ("fermentations-ernte",
     "Zellernte aus dem Bioreaktor durch Zentrifugation oder Filtration"),
    ("fermentation harvest",
     "Cell harvest from bioreactor by centrifugation or filtration"),
    # --- Downstream (Natural products) ---
    ("mahlung",
     "Mechanischer Aufschluss des Rohmaterials für effiziente Extraktion"),
    ("lösungsmittel-extraktion",
     "Selektive Extraktion mit Methanol/Ethanol/CO2 je nach Polarität der Zielsubstanz"),
    ("solvent extraction",
     "Selective extraction with methanol/ethanol/CO2 depending on target polarity"),
    ("flüssig-flüssig-extraktion",
     "Überführung des Zielmoleküls in eine organische Phase (z. B. Ethylacetat)"),
    ("fraktionierung",
     "Grobe Trennung nach Polarität oder Molekulargewicht vor der Chromatographie"),
    ("säulen-chromatographie",
     "Feinaufreinigung an Silica oder umgekehrt-Phase — Trennung verbleibender Verunreinigungen"),
    ("column chromatography",
     "Feinaufreinigung an Silica oder umgekehrt-Phase — Trennung verbleibender Verunreinigungen"),
    ("aufkonzentrierung",
     "Lösungsmittelentfernung im Rotationsverdampfer — Anreicherung der Zielsubstanz"),
    ("zellabtrennung",
     "Entfernt Zellmasse aus dem Fermentationsmedium via Zentrifugation oder Filtration"),
    ("kristallisation",
     "Finale Aufreinigung zu kristalliner API-Form — typischerweise Pharma-Grade-Qualität"),
    ("umkristallisation",
     "Erneute Kristallisation aus Wasser/Ethanol für höchste Reinheit"),
    # --- Downstream (Small molecule generic) ---
    ("extraktion",
     "Trennung von Reaktionsnebenprodukten via Lösungsmittel"),
    ("eindampfen",
     "Lösungsmittelentfernung im Rotationsverdampfer"),
    ("flash-chromatographie",
     "Säulen-Chromatographie zur Trennung verbleibender Verunreinigungen"),
    ("aufarbeitung",
     "Standard-Workup (Filtration, Trocknung) nach der Synthese"),
    ("polishing",
     "Finale Aufreinigung zur Erreichung der Reinheits-Spezifikation"),
]


_BLOCK_EXPLANATIONS_EN: List[Tuple[str, str]] = [
    # --- Upstream (biotech) ---
    ("strain / inoculum",
     "Production organism from the strain bank — defined clone backup, multiply tested"),
    ("inoculum",
     "Growing the production strain from the bank"),
    ("pre-culture",
     "Biomass build-up in shake flasks, typically 16–24 h growth before the main reactor"),
    ("main fermentation",
     "Fed-batch in the bioreactor, product accumulation under process control (pH, O2, T)"),
    ("fermentation",
     "Microbial production in the bioreactor, product accumulates extracellularly or in cell pellet"),
    # --- Upstream (extraction) ---
    ("raw material",
     "Plant/biological starting material — quality determines yield and impurity profile"),
    ("pre-treatment",
     "Mechanical and thermal preparation (milling, drying) for efficient extraction"),
    # --- Upstream (chemical) ---
    ("step 1",
     "First synthesis stage — introduction of the core scaffold or protecting-group strategy"),
    ("step 2",
     "Functionalisation — substitutions, couplings or stereo introduction"),
    ("step 3",
     "Final transformation and purification to the target molecule"),
    ("synthesis",
     "Main chemical synthesis, summary step"),
    # --- Downstream (mAb) ---
    ("protein-a capture",
     "Affinity chromatography — > 95 % API selectivity in one step"),
    ("low-ph virus inact",
     "Inactivation of viral contaminants at pH 3.5 for 60–90 min — GMP standard for mammalian cells"),
    ("cation iex",
     "Cation-exchange chromatography — removes aggregates and host proteins via charge differences"),
    ("anion iex",
     "Anion-exchange + virus filtration — DNA/endotoxin removal as a polishing step"),
    ("tff + formulation",
     "Tangential-flow filtration — buffer exchange and concentration to target specification"),
    ("tff",
     "Tangential-flow filtration — buffer exchange and concentration"),
    # --- Downstream (Protein general) ---
    ("cell lysis",
     "Mechanical or chemical cell disruption to release intracellular proteins"),
    ("clarification",
     "Removal of cell debris via depth filtration or centrifugation"),
    ("affinity/iex",
     "Affinity or ion-exchange chromatography for primary capture"),
    ("sec polish",
     "Size-exclusion chromatography — aggregate/monomer separation as polishing"),
    ("uf/df",
     "Ultrafiltration/diafiltration — buffer exchange and concentration"),
    # --- Downstream (Peptide) ---
    ("spps cleavage",
     "Cleavage from synthesis resin with TFA cocktail, also removes protecting groups"),
    ("prep-hplc (rp)",
     "Preparative reversed-phase HPLC — main purification to > 95 % purity"),
    ("prep-hplc",
     "Preparative HPLC — high-resolution purification"),
    ("counter-ion exchange",
     "Exchange TFA for a pharmacologically acceptable counter-ion (usually acetate)"),
    ("lyophilisation",
     "Freeze-drying — removes water for stable storage formulation"),
    ("fermentation harvest",
     "Cell harvest from bioreactor by centrifugation or filtration"),
    # --- Downstream (Natural products) ---
    ("milling",
     "Mechanical disruption of raw material for efficient extraction"),
    ("solvent extraction",
     "Selective extraction with methanol/ethanol/CO2 depending on target polarity"),
    ("fractionation",
     "Coarse separation by polarity or molecular weight before chromatography"),
    ("column chromatography",
     "Fine purification on silica or reverse phase — removes remaining impurities"),
    ("crystallisation",
     "Final purification to crystalline API form — typically pharma-grade quality"),
    # --- Downstream (Small molecule generic) ---
    ("extraction",
     "Removal of reaction side-products via solvent"),
    ("evaporation",
     "Solvent removal in rotary evaporator"),
    ("flash chromatography",
     "Column chromatography for separating residual impurities"),
    ("workup",
     "Standard workup (filtration, drying) after synthesis"),
    ("polishing",
     "Final purification to reach the purity specification"),
]


def _find_explanation(label: str, lang: str = "de") -> Optional[str]:
    """Return a one-sentence explanation for `label` or None if no match."""
    if not label:
        return None
    table = _BLOCK_EXPLANATIONS_DE if lang == "de" else _BLOCK_EXPLANATIONS_EN
    s = str(label).strip().lower()
    for key, expl in table:
        if key in s:
            return expl
    return None


# ---------------------------------------------------------------------------
# Text wrapping for blocks
# ---------------------------------------------------------------------------

def _wrap_label(text: str, max_chars_per_line: int = 22) -> List[str]:
    """Wrap `text` onto up to 2 lines, breaking at natural break points.

    Strategy:
      1. If text fits on one line → return [text]
      2. Find a space, "/", or "+" near the midpoint and split there
      3. If still longer, split at the latest possible break before
         max_chars_per_line and ellipsis the second line if needed
    """
    if not text:
        return [""]
    s = str(text).strip()
    if len(s) <= max_chars_per_line:
        return [s]

    # Find break candidates (positions of space or punctuation we can break on)
    breakable = {" ", "/", "+", "-"}
    # Aim for break around max_chars_per_line, but allow a window
    target = max_chars_per_line
    window_lo = max(8, target - 8)
    window_hi = min(len(s) - 1, target)
    best = -1
    for i in range(window_hi, window_lo - 1, -1):
        if s[i] in breakable:
            best = i
            break
    if best == -1:
        # No good break — split at target
        best = target
    line1 = s[:best].rstrip(" /+-")
    line2 = s[best:].lstrip(" /+-")
    if len(line2) > max_chars_per_line:
        line2 = line2[: max_chars_per_line - 1] + "…"
    return [line1, line2]


# ---------------------------------------------------------------------------
# Step gathering — pulls upstream + downstream labels for the diagram
# ---------------------------------------------------------------------------

def _gather_downstream_steps(
    process_input: Dict[str, Any],
    concrete_recs: Optional[Dict[str, Any]],
    lang: str = "de",
) -> List[str]:
    """Return the downstream-step labels for the current molecule context.

    Priority:
    1. concrete_recommendations.downstream_steps (curated, most specific)
    2. Inferred from molecule_type + available methods (generic fallback)
    """
    if concrete_recs and isinstance(concrete_recs, dict):
        steps = concrete_recs.get("downstream_steps") or []
        if steps:
            return list(steps[:6])

    mtype = (process_input.get("molecule_type") or "").lower()
    msub = (process_input.get("molecule_subtype") or "").lower()
    method = (process_input.get("method") or "").lower()

    if mtype == "protein" and msub == "antibody":
        return (["Protein-A-Capture", "Low-pH-Virusinakt.", "Kation-IEX",
                 "Anion-IEX/VF", "TFF + Formulation"]
                if lang == "de" else
                ["Protein-A capture", "Low-pH virus inact.", "Cation IEX",
                 "Anion IEX/VF", "TFF + formulation"])

    if mtype == "protein":
        return (["Zell-Lyse", "Klärfiltration", "Affinitäts/IEX",
                 "SEC-Polishing", "UF/DF + Lyo"]
                if lang == "de" else
                ["Cell lysis", "Clarification", "Affinity/IEX",
                 "SEC polish", "UF/DF + lyo"])

    if mtype == "peptide":
        if method in ("biotechnological", "biotech"):
            return (["Fermentations-Ernte", "Klärung + Capture",
                     "Prep-HPLC (RP)", "Lyophilisation"]
                    if lang == "de" else
                    ["Fermentation harvest", "Clarification + capture",
                     "Prep-HPLC (RP)", "Lyophilisation"])
        return (["SPPS-Cleavage", "Prep-HPLC (RP)",
                 "Counter-Ion-Exchange", "Lyophilisation"]
                if lang == "de" else
                ["SPPS cleavage", "Prep-HPLC (RP)",
                 "Counter-ion exchange", "Lyophilisation"])

    if mtype == "natural_product":
        if method in ("extraction", "extract"):
            return (["Mahlung", "Lösungsmittel-Extraktion", "Fraktionierung",
                     "Säulen-Chromatographie", "Kristallisation"]
                    if lang == "de" else
                    ["Milling", "Solvent extraction", "Fractionation",
                     "Column chromatography", "Crystallisation"])
        return (["Fermentations-Ernte", "Extraktion",
                 "Säulen-Chromatographie", "Kristallisation"]
                if lang == "de" else
                ["Fermentation harvest", "Extraction",
                 "Column chromatography", "Crystallisation"])

    # Small molecule — derived from available methods
    methods = process_input.get("available_methods") or {}
    chain = []
    if methods.get("has_extraction"):
        chain.append("Extraktion" if lang == "de" else "Extraction")
    if methods.get("has_rotary_evaporator"):
        chain.append("Eindampfen" if lang == "de" else "Evaporation")
    if methods.get("has_crystallization"):
        chain.append("Kristallisation" if lang == "de" else "Crystallisation")
    elif methods.get("has_flash_chromatography"):
        chain.append("Flash-Chromatographie" if lang == "de" else "Flash chromatography")
    if methods.get("has_prep_hplc"):
        chain.append("prep-HPLC" if lang == "de" else "prep-HPLC")
    if not chain:
        chain = (["Aufarbeitung", "Polishing"]
                 if lang == "de" else ["Workup", "Polishing"])
    return chain


def _gather_upstream_steps(
    process_input: Dict[str, Any],
    concrete_recs: Optional[Dict[str, Any]] = None,
    lang: str = "de",
) -> List[str]:
    """Build upstream/synthesis-step labels.

    Priority (Feedback-Runde-3):
    1. concrete_recommendations.production_steps (curated, most specific)
    2. Inferred from method (generic fallback)
    """
    # Curated production steps win over generic fallback
    if concrete_recs and isinstance(concrete_recs, dict):
        prod = concrete_recs.get("production_steps") or []
        if prod:
            return list(prod[:6])

    method = (process_input.get("method") or "").lower()
    n_steps = process_input.get("number_of_steps")
    n = int(n_steps) if isinstance(n_steps, (int, float)) and n_steps else None

    if method in ("biotechnological", "biotech"):
        return (["Stamm / Inokulum", "Vorkultur", "Hauptfermentation"]
                if lang == "de" else
                ["Strain / inoculum", "Pre-culture", "Main fermentation"])

    if method in ("extraction", "extract"):
        return (["Rohmaterial", "Vorbehandlung"]
                if lang == "de" else
                ["Raw material", "Pre-treatment"])

    if n and n > 0:
        if lang == "de":
            base = [f"Schritt {i}" for i in range(1, min(n, 3) + 1)]
        else:
            base = [f"Step {i}" for i in range(1, min(n, 3) + 1)]
        if n > 3:
            base.append(f"… (+{n - 3})")
        return base

    return ["Synthese"] if lang == "de" else ["Synthesis"]


# ---------------------------------------------------------------------------
# Public: build the Drawing
# ---------------------------------------------------------------------------

def build_process_flow_drawing(
    process_input: Dict[str, Any],
    concrete_recs: Optional[Dict[str, Any]] = None,
    lang: str = "de",
    width_mm: float = 175.0,
):
    """Return a ReportLab Drawing with the block diagram, or None on failure.

    Layout: Upstream row + arrow + Downstream row.
    Max 4 blocks per row. Block text wraps to 2 lines when needed.
    """
    try:
        from reportlab.graphics.shapes import Drawing, Rect, String, PolyLine, Polygon
        from reportlab.lib.colors import HexColor
        from reportlab.lib.units import mm
    except Exception:
        return None

    upstream = _gather_upstream_steps(process_input, concrete_recs, lang)
    downstream = _gather_downstream_steps(process_input, concrete_recs, lang)
    if not upstream and not downstream:
        return None

    # Cap blocks per row at 4 for readability
    MAX_PER_ROW = 4
    upstream = upstream[:MAX_PER_ROW]
    downstream = downstream[:MAX_PER_ROW]

    # Layout constants (points; 1 pt = 1/72 inch)
    width_pt = width_mm * mm
    pad = 2 * mm
    avail_w = width_pt - 2 * pad
    block_h = 22 * mm        # taller for 2-line text
    row_gap = 12 * mm
    header_h = 7 * mm
    line_height = 11         # for fontSize=8 Bold
    total_h = header_h + block_h + row_gap + block_h + 4 * mm

    drawing = Drawing(width_pt, total_h)

    def _block_width(n: int) -> float:
        return (avail_w - (n - 1) * 4) / max(n, 1)

    def _render_row(items: List[str], y: float, fill_color: str, max_chars: int):
        n = len(items)
        bw = _block_width(n)
        for i, label in enumerate(items):
            x = pad + i * (bw + 4)
            drawing.add(Rect(x, y, bw, block_h,
                             fillColor=HexColor(fill_color),
                             strokeColor=HexColor(fill_color),
                             strokeWidth=0.5))
            # Text wrap onto up to 2 lines
            lines = _wrap_label(label, max_chars_per_line=max_chars)
            # Vertical centering: 1 line centered, 2 lines split symmetrically
            if len(lines) == 1:
                drawing.add(String(x + bw / 2, y + block_h / 2 - 3,
                                   lines[0],
                                   fontName="Helvetica-Bold", fontSize=8,
                                   fillColor=HexColor("#FFFFFF"),
                                   textAnchor="middle"))
            else:
                # 2 lines, centred around middle
                drawing.add(String(x + bw / 2, y + block_h / 2 + 2,
                                   lines[0],
                                   fontName="Helvetica-Bold", fontSize=8,
                                   fillColor=HexColor("#FFFFFF"),
                                   textAnchor="middle"))
                drawing.add(String(x + bw / 2, y + block_h / 2 - 8,
                                   lines[1],
                                   fontName="Helvetica-Bold", fontSize=8,
                                   fillColor=HexColor("#FFFFFF"),
                                   textAnchor="middle"))
            # Small connector arrow between blocks
            if i < n - 1:
                ax = x + bw
                ay = y + block_h / 2
                drawing.add(PolyLine([ax, ay, ax + 4, ay],
                                     strokeColor=HexColor(BRAND_NEUTRAL),
                                     strokeWidth=0.7))

    # Estimate per-line char capacity from block width (avg char ≈ 4.5pt at 8pt Helvetica)
    up_chars = max(12, int((_block_width(len(upstream)) - 6) / 4.5))
    dn_chars = max(12, int((_block_width(len(downstream)) - 6) / 4.5))

    # --- Upstream row ---
    y_up = total_h - header_h - block_h
    drawing.add(String(pad, total_h - 4, "1. Upstream",
                       fontName="Helvetica-Bold", fontSize=9,
                       fillColor=HexColor(BRAND_PRIMARY)))
    _render_row(upstream, y_up, BRAND_PRIMARY, up_chars)

    # --- Vertical arrow ---
    mid_x = width_pt / 2
    arrow_top = y_up
    arrow_bot = y_up - row_gap + 2 * mm
    drawing.add(PolyLine([mid_x, arrow_top, mid_x, arrow_bot + 4],
                         strokeColor=HexColor(BRAND_NEUTRAL),
                         strokeWidth=1.0))
    drawing.add(Polygon(points=[mid_x - 3, arrow_bot + 4,
                                mid_x + 3, arrow_bot + 4,
                                mid_x, arrow_bot - 1],
                        fillColor=HexColor(BRAND_NEUTRAL),
                        strokeColor=HexColor(BRAND_NEUTRAL)))

    # --- Downstream row ---
    y_dn = arrow_bot - block_h
    dn_header = ("2. Downstream / Aufarbeitung"
                 if lang == "de" else "2. Downstream / Workup")
    drawing.add(String(pad, y_dn + block_h + 2, dn_header,
                       fontName="Helvetica-Bold", fontSize=9,
                       fillColor=HexColor(BRAND_ACCENT)))
    _render_row(downstream, y_dn, BRAND_ACCENT, dn_chars)

    return drawing


# ---------------------------------------------------------------------------
# Public: build the legend
# ---------------------------------------------------------------------------

def build_process_flow_legend(
    process_input: Dict[str, Any],
    concrete_recs: Optional[Dict[str, Any]] = None,
    lang: str = "de",
) -> List[Dict[str, str]]:
    """Return a structured legend describing each block in the diagram.

    Each entry is a dict:
        { "section": "upstream" | "downstream",
          "label": <block label as shown in the diagram>,
          "explanation": <one-sentence what-this-step-does> }

    Used by the PDF report to render a description list directly below
    the diagram, addressing the "what does Stamm/Inokulum mean?" gap.
    """
    upstream = _gather_upstream_steps(process_input, concrete_recs, lang)[:4]
    downstream = _gather_downstream_steps(process_input, concrete_recs, lang)[:4]

    fallback_de = "Standardschritt für diese Prozessklasse — Details abhängig vom konkreten Equipment."
    fallback_en = "Standard step for this process class — details depend on actual equipment."

    legend: List[Dict[str, str]] = []
    for lbl in upstream:
        expl = _find_explanation(lbl, lang) or (fallback_de if lang == "de" else fallback_en)
        legend.append({"section": "upstream", "label": lbl, "explanation": expl})
    for lbl in downstream:
        expl = _find_explanation(lbl, lang) or (fallback_de if lang == "de" else fallback_en)
        legend.append({"section": "downstream", "label": lbl, "explanation": expl})
    return legend
