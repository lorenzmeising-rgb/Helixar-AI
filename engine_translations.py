"""Translation layer for dynamic engine outputs.

The decision_engine.py module produces dynamic strings (issues, improvements,
cost drivers, downstream insights, tradeoffs) hardcoded in English. This module
maps them to German on the rendering side, without touching the 2,400-line
engine itself.

Strategy:
  1. PHRASE_MAP holds the most frequently produced English phrases (~30 entries)
     and their German translations. Substring matching keeps it robust against
     dynamic `(count=N)` style insertions from the engine.
  2. VALUE_MAP translates the qualitative bucket labels (HIGH / LOW / MEDIUM)
     and a few level-style words ("high", "low", "medium") used inline in
     paragraph text.
  3. translate_engine_text() applies value substitution AND phrase substitution
     to a single string. If no rule matches, the input is returned unchanged
     (graceful fallback — better English than crash).

Public API:
  - translate_engine_text(text, lang) -> str
  - translate_engine_value(value, lang) -> str   # for HIGH/LOW/MEDIUM cells
"""

import re
from typing import Optional


# ---------------------------------------------------------------------------
# Value-level translations (qualitative buckets used as cell values)
# ---------------------------------------------------------------------------

VALUE_MAP = {
    # uppercase bucket values shown verbatim in PDF
    "HIGH": {"de": "HOCH", "en": "HIGH"},
    "LOW": {"de": "NIEDRIG", "en": "LOW"},
    "MEDIUM": {"de": "MITTEL", "en": "MEDIUM"},
    "VERY HIGH": {"de": "SEHR HOCH", "en": "VERY HIGH"},
    "VERY LOW": {"de": "SEHR NIEDRIG", "en": "VERY LOW"},
    "N/A": {"de": "n/v", "en": "N/A"},
    # lowercase variants embedded in body text
    "high": {"de": "hoch", "en": "high"},
    "low": {"de": "niedrig", "en": "low"},
    "medium": {"de": "mittel", "en": "medium"},
    "very high": {"de": "sehr hoch", "en": "very high"},
    "very low": {"de": "sehr niedrig", "en": "very low"},
    # confidence-style words
    "High": {"de": "Hoch", "en": "High"},
    "Low": {"de": "Niedrig", "en": "Low"},
    "Medium": {"de": "Mittel", "en": "Medium"},
}


def translate_engine_value(value: Optional[str], lang: str = "de") -> str:
    """Translate a single bucket value (HIGH/LOW/MEDIUM and variants)."""
    if value is None:
        return ""
    s = str(value)
    entry = VALUE_MAP.get(s) or VALUE_MAP.get(s.strip())
    if not entry:
        return s
    return entry.get(lang) or entry.get("de") or s


# ---------------------------------------------------------------------------
# Phrase-level translations
# ---------------------------------------------------------------------------

# Pairs of (regex pattern matched against the engine string, German replacement).
# The regex form lets us absorb dynamic insertions like "(count=2)" or numbers.
# Order matters: more specific patterns first.
_PHRASE_RULES = [
    # ---- Issues / observations ----
    (r"^Aromatic rings detected \(count=(\d+)\)\s*[—-]\s*supports crystallization screening.*$",
     r"Aromatische Ringe erkannt (Anzahl=\1) — unterstützt Kristallisations-Screening und potenzielle Festphasen-Trennung"),
    (r"^Aromatic rings detected \(count=(\d+)\)\s*[—-]\s*increases chemical stability.*$",
     r"Aromatische Ringe erkannt (Anzahl=\1) — erhöht chemische Stabilität und ermöglicht oft Kristallisations-Reinigung"),
    (r"^Aromatic rings detected \(count=(\d+)\)\s*[—-]\s*enable crystallization screening.*$",
     r"Aromatische Ringe erkannt (Anzahl=\1) — Kristallisations-Screening kann HPLC-Bedarf reduzieren"),

    (r"^Multiple polar groups suggest reversed-phase HPLC or polarity-based partitioning.*$",
     r"Mehrere polare Gruppen — RP-HPLC oder polaritätsbasierte Phasen-Trennung wird empfohlen"),

    (r"^Difficult purification \(HPLC/chromatography or distillation\) increases solvent, buffer and resin costs.*$",
     r"Anspruchsvolle Aufarbeitung (HPLC/Chromatographie oder Destillation) — erhöht Lösungsmittel-, Puffer- und Harzkosten und verlängert die Aufarbeitungszeit"),
    (r"^Difficult purification \(chromatography/HPLC\) increases solvent.*$",
     r"Anspruchsvolle Aufarbeitung (Chromatographie/HPLC) — erhöht Lösungsmittel-, Harz- und Personalkosten"),
    # New short / structure-context variants emitted by the engine (caught by batch test):
    (r"^Difficult purification increases chromatography/HPLC time and solvent consumption \(aromatic system, multiple polar groups\).*$",
     r"Anspruchsvolle Aufarbeitung — verlängert Chromatographie/HPLC-Zeit und erhöht Lösungsmittelverbrauch (aromatisches System, mehrere polare Gruppen)"),
    (r"^Difficult purification increases chromatography/HPLC time and solvent consumption \(hydrophobic chain, multiple polar groups\).*$",
     r"Anspruchsvolle Aufarbeitung — verlängert Chromatographie/HPLC-Zeit und erhöht Lösungsmittelverbrauch (hydrophobe Kette, mehrere polare Gruppen)"),
    (r"^Difficult purification increases chromatography/HPLC time and solvent consumption.*$",
     r"Anspruchsvolle Aufarbeitung — verlängert Chromatographie/HPLC-Zeit und erhöht den Lösungsmittelverbrauch"),
    (r"^Difficult purification expected due to: .*$",
     r"Anspruchsvolle Aufarbeitung erwartet — strukturelle Gründe identifiziert"),
    (r"^Difficult downstream separations inferred from structure.*$",
     r"Strukturbasiert anspruchsvolle Aufarbeitung erwartet"),
    (r"^Challenging purification increases downstream cost.*$",
     r"Anspruchsvolle Reinigung erhöht die Aufarbeitungskosten"),

    (r"^Biotechnological production increases contamination and aggregation risk.*$",
     r"Biotechnologische Produktion erhöht Kontaminations- und Aggregationsrisiko — erfordert sterile Bioreaktor-Operationen und validierte Medien (höhere Betriebskosten)"),
    (r"^Biotech process(es)? requires? careful control.*$",
     r"Biotechnologische Prozesse erfordern sorgfältige Kontrolle der biologischen Systeme"),
    (r"^Biotech process without bioreactor limits feasibility.*$",
     r"Biotechnologischer Prozess ohne Bioreaktor — Machbarkeit eingeschränkt"),

    (r"^Antibody downstream dominated by affinity capture.*$",
     r"Antikörper-Aufarbeitung dominiert von Affinitäts-Capture und Multi-Step-Chromatographie; Resin-Wiederverwendung und CIP-Zyklen einplanen"),
    (r"^Antibody purification is highly complex.*$",
     r"Antikörper-Aufreinigung ist hochkomplex und teuer"),
    # Antibody recovery (doi:... citations stay as-is, only the prefix is translated)
    (r"^Antibody recovery and purification \(doi:([^\)]+)\).*$",
     r"Antikörper-Gewinnung und -Aufreinigung (DOI: \1) — Standard-Sequenz: Protein-A → Low-pH-VI → IEX → AEX/VF → TFF"),
    (r"^Antibody recovery and purification.*$",
     r"Antikörper-Gewinnung und -Aufreinigung — Standard-Sequenz: Protein-A-Capture, Low-pH-Virusinaktivierung, IEX-Polishing, TFF-Konzentrierung"),

    # Peptide downstream guidance
    (r"^Peptide downstreams often require preparative HPLC \(reverse-phase\) and orthogonal polishing steps;\s*plan for solvent recycling.*$",
     r"Peptid-Aufarbeitung erfordert üblicherweise präparative HPLC (Reversed-Phase) und orthogonale Polishing-Schritte; Lösungsmittel-Recycling einplanen"),
    (r"^Peptide downstreams require preparative HPLC.*$",
     r"Peptid-Aufarbeitung erfordert präparative HPLC und orthogonale Polishing-Schritte"),

    # Natural product downstream guidance
    (r"^Natural products often require selective extraction, fractionation and targeted chromatographic separation due to complex mixtures.*$",
     r"Naturstoffe erfordern selektive Extraktion, Fraktionierung und gezielte chromatographische Trennung aufgrund komplexer Gemische"),
    (r"^Natural products require selective extraction.*$",
     r"Naturstoffe erfordern selektive Extraktion und chromatographische Trennung"),

    # Protein downstream guidance (proactive — likely to surface for enzymes/non-antibodies)
    (r"^Protein downstream typically uses ion-exchange, hydrophobic-interaction and size-exclusion chromatography.*$",
     r"Protein-Aufarbeitung nutzt typischerweise Ionenaustausch-, hydrophobe Interaktions- und Größenausschluss-Chromatographie"),
    (r"^Protein downstream requires multi-step chromatography.*$",
     r"Protein-Aufarbeitung erfordert mehrstufige Chromatographie und Polishing-Schritte"),

    (r"^Complex folding increases process development.*$",
     r"Komplexe Faltung erhöht Prozessentwicklungs- und Refolding-Kosten"),
    (r"^Complex folding increases process sensitivity.*$",
     r"Komplexe Faltung erhöht Prozessempfindlichkeit und Instabilität"),
    (r"^Complex folding increases sensitivity to process conditions.*$",
     r"Komplexe Faltung erhöht Empfindlichkeit gegenüber Prozessbedingungen und treibt Refolding-/Aufarbeitungskosten"),

    (r"^Complex molecule combined with many synthesis steps.*$",
     r"Komplexes Molekül kombiniert mit vielen Syntheseschritten — ineffizienter Prozess"),
    (r"^Complex molecule plus many steps reduces efficiency.*$",
     r"Komplexes Molekül plus viele Schritte — reduziert Effizienz und erhöht Kosten"),
    (r"^Complex molecule likely requires expensive raw materials.*$",
     r"Komplexes Molekül erfordert wahrscheinlich teure Rohstoffe"),
    (r"^Complex mixture of similar compounds increases process risk.*$",
     r"Komplexes Gemisch ähnlicher Verbindungen erhöht das Prozessrisiko"),

    (r"^Cyclic peptides are more stable but harder to purify.*$",
     r"Cyclische Peptide sind stabiler, aber schwerer zu reinigen"),

    (r"^Crystallization or solvent-based purification expected due to aromatic/rigid structure.*$",
     r"Kristallisation oder lösungsmittelbasierte Reinigung empfohlen (aromatische/rigide Struktur)"),
    (r"^Crystallization or solvent-based purification expected for non-volatile small molecules.*$",
     r"Kristallisation oder lösungsmittelbasierte Reinigung für nicht-flüchtige kleine Moleküle empfohlen"),
    (r"^Crystallization or solvent-based purification is expected.*$",
     r"Kristallisation oder lösungsmittelbasierte Reinigung erwartet"),

    (r"^Higher purity targets increase purification cost and process complexity.*$",
     r"Höhere Reinheitsziele erhöhen Aufreinigungskosten und Prozesskomplexität"),

    # ---- Raw-material / cost issues ----
    (r"^Expensive or scarce starting materials increase COGS and procurement risk.*$",
     r"Teure oder knappe Rohstoffe erhöhen die Herstellkosten (COGS) und das Beschaffungsrisiko"),
    (r"^Expensive raw materials drive up unit cost.*$",
     r"Teure Rohstoffe treiben die Stückkosten in die Höhe"),
    (r"^Raw material shortages or single-source supply.*$",
     r"Rohstoff-Engpässe oder Single-Source-Versorgung erhöhen das Versorgungsrisiko"),

    # ---- Biophysical stability / aggregation issues ----
    (r"^Low biophysical stability increases degradation during processing and reduces batch recovery.*$",
     r"Geringe biophysikalische Stabilität erhöht Degradation während des Prozesses und reduziert die Batch-Ausbeute"),
    (r"^High aggregation risk reduces yield and increases polishing and formulation costs.*$",
     r"Hohes Aggregations-Risiko reduziert die Ausbeute und erhöht Polishing- und Formulierungskosten"),
    (r"^Aggregation risk drives up downstream losses.*$",
     r"Aggregations-Risiko erhöht die Aufarbeitungs-Verluste"),
    (r"^Thermal instability requires cold chain or cryoprotectant.*$",
     r"Thermische Instabilität erfordert Kühlkette oder Cryoprotectant-Formulierung"),

    # ---- Step / process complexity cost drivers ----
    (r"^Multiple synthesis steps \((\d+)\) increase reagent, labour and unit-operation costs.*$",
     r"Mehrstufige Synthese (\1 Schritte) erhöht Reagenzien-, Personal- und Betriebskosten"),
    # Engine emits both "operational complexity" and "operational costs (reagents, labour, unit ops)":
    (r"^Multiple synthesis steps \((\d+)\) increase operational costs \(reagents, labour, unit ops\).*$",
     r"Mehrstufige Synthese (\1 Schritte) erhöht die Betriebskosten (Reagenzien, Personal, Unit-Operations)"),
    (r"^Multiple synthesis steps increase operational complexity.*$",
     r"Mehrstufige Synthese erhöht die operative Komplexität"),
    (r"^Many synthesis steps drive up reagent and labour costs.*$",
     r"Viele Syntheseschritte treiben Reagenzien- und Personalkosten in die Höhe"),
    (r"^Each additional step reduces overall yield.*$",
     r"Jeder zusätzliche Schritt reduziert den Gesamt-Yield"),
    (r"^High synthesis step count increases operational complexity, extends cycle time, and limits scalability.*$",
     r"Hohe Schrittzahl erhöht operative Komplexität, verlängert die Durchlaufzeit und limitiert die Skalierbarkeit"),

    # ---- Process simplification / scalability suggestions ----
    (r"^Simplify process design for improved scalability.*$",
     r"Prozessdesign vereinfachen für bessere Skalierbarkeit"),
    (r"^Consider process simplification.*$",
     r"Prozessvereinfachung in Betracht ziehen"),
    (r"^Reduce reagent inventory and waste streams.*$",
     r"Reagenzien-Inventar und Abfallströme reduzieren"),

    # ---- Waste / compliance cost drivers ----
    (r"^Strict waste handling and disposal requirements \(compliance cost\).*$",
     r"Strenge Abfall- und Entsorgungsauflagen (Compliance-Kosten)"),
    (r"^Waste treatment compliance adds overhead.*$",
     r"Compliance bei Abfallbehandlung verursacht zusätzlichen Aufwand"),

    # ---- Purity / trade-off statements ----
    (r"^Purity targets above standard require advanced purification assets or increase risk/cost.*$",
     r"Reinheitsziele über Standardniveau erfordern hochwertige Aufreinigungs-Assets oder erhöhen Risiko/Kosten"),
    (r"^High purity requires advanced purification.*$",
     r"Hohe Reinheit erfordert hochwertige Aufreinigung"),
    (r"^Aromatic / crystallizable structure.*HPLC not strictly required.*$",
     r"Aromatische / kristallisierbare Struktur — Kristallisations-basierte Aufarbeitung typ. ausreichend, HPLC nicht zwingend erforderlich"),
    (r"^Strict waste regulations increase process cost and complexity.*$",
     r"Strenge Abfallauflagen erhöhen Prozesskosten und -komplexität"),
    (r"^High purity significantly increases downstream complexity.*$",
     r"Hohe Reinheit erhöht die Aufarbeitungs-Komplexität deutlich"),
    (r"^Chemical synthesis includes energy-intensive steps.*$",
     r"Chemische Synthese beinhaltet energieintensive Schritte"),
    (r"^Process may conflict with waste regulations.*$",
     r"Prozess kann mit Abfallauflagen kollidieren"),
    (r"^Scaling introduces operational and mixing challenges.*$",
     r"Skalierung bringt operative und Durchmischungs-Herausforderungen mit sich"),
    (r"^Energy-intensive separation steps likely dominate cost.*$",
     r"Energieintensive Trenntechnik dominiert wahrscheinlich die Kosten"),
    (r"^Extremely challenging purification expected.*$",
     r"Extrem anspruchsvolle Aufarbeitung erwartet"),
    (r"^Low molecular stability increases process risk.*$",
     r"Niedrige molekulare Stabilität erhöht das Prozessrisiko"),
    (r"^High number of synthesis steps increases cost and complexity.*$",
     r"Hohe Anzahl an Syntheseschritten erhöht Kosten und Komplexität"),
    (r"^Protein stability and aggregation risk must be considered.*$",
     r"Protein-Stabilität und Aggregationsrisiko müssen berücksichtigt werden"),
    (r"^Very high purity targets increase downstream separation cost.*$",
     r"Sehr hohe Reinheitsziele erhöhen die Aufarbeitungskosten (längere HPLC / zusätzliche Polishing-Schritte)"),

    (r"^Chromatography likely a major cost driver.*$",
     r"Chromatographie ist wahrscheinlich ein Hauptkostentreiber"),
    (r"^Chemical synthesis may involve energy-intensive steps.*$",
     r"Chemische Synthese kann energieintensive Schritte enthalten"),
    (r"^Cost driven by standard materials and moderate downstream operations.*$",
     r"Kosten getrieben durch Standardmaterialien und moderate Aufarbeitung"),
    (r"^Cost explanation not available.*$",
     r"Kostenerklärung nicht verfügbar"),

    (r"^Amide groups present\s*[—-]\s*monitor stereochemical integrity.*$",
     r"Amid-Gruppen vorhanden — stereochemische Integrität und Hydrolysestabilität überwachen"),

    # ---- Confidence / generic notes ----
    (r"^High confidence\s*[—-]\s*suitable for pilot studies with standard QA/QC\.?$",
     r"Hohe Konfidenz — geeignet für Pilotstudien mit Standard-QA/QC"),
    (r"^Medium confidence.*$",
     r"Mittlere Konfidenz — Validierung mit Pilotdaten empfohlen"),
    (r"^Low confidence.*$",
     r"Niedrige Konfidenz — zusätzliche Literaturrecherche und Pilotdaten empfohlen"),
    (r"^No operating parameters available; recommend literature review and lab data collection\.?$",
     r"Keine Betriebsparameter verfügbar — Literaturrecherche und Labordaten-Erhebung empfohlen"),

    # ---- Improvements / actions ----
    (r"^Adjust process conditions to improve molecular stability.*$",
     r"Prozessbedingungen zur Verbesserung der molekularen Stabilität anpassen"),
    (r"^Avoid harsh pH during workup; consider protecting labile groups.*$",
     r"Aggressive pH-Werte bei der Aufarbeitung vermeiden; Schutzgruppen oder selektive Transformationen für labile Gruppen erwägen"),
    (r"^Assess convergent synthesis or protecting-group minimization.*$",
     r"Konvergente Synthese oder Minimierung von Schutzgruppen prüfen — reduziert Schrittzahl und verbessert Gesamt-Yield"),
    (r"^Consider alternative purification methods such as crystallization.*$",
     r"Alternative Aufreinigungsmethoden wie Kristallisation prüfen"),
    (r"^Consider phase-separation or flash chromatography.*$",
     r"Phasen-Trennung oder Flash-Chromatographie mit hydrophobie-orientierten Lösungsmittelsystemen prüfen; grünere Lösungsmittel zur Kostenreduktion bewerten"),
    (r"^Design multi-modal chromatography train for enzyme purification.*$",
     r"Multi-modaler Chromatographie-Zug für Enzym-Aufreinigung (Capture + Polish); Feedstock-Klarheit zur Minimierung von Säulen-Fouling bewerten"),
    (r"^Design safer reagent choices and include effluent neutralization.*$",
     r"Sicherere Reagenzien wählen und Abwasser-Neutralisation in den Prozessfluss integrieren"),
    (r"^Design solvent-based extraction and flash-chromatography workflows.*$",
     r"Lösungsmittel-basierte Extraktion und Flash-Chromatographie mit Lösungsmittel-Minimierung entwerfen"),
    (r"^Develop advanced chromatographic methods for cyclic peptides.*$",
     r"Fortgeschrittene Chromatographie-Methoden für cyclische Peptide entwickeln (Resin/Gradient/Desolvation optimieren)"),
    (r"^Develop preparative reversed-phase HPLC methods.*$",
     r"Präparative RP-HPLC-Methoden entwickeln (Gradient, Säulenbeladung, Lösungsmittel-Recovery optimieren)"),
    (r"^Develop selective extraction and fractionation.*$",
     r"Selektive Extraktion und Fraktionierung gefolgt von gezielter Chromatographie entwickeln"),
    (r"^Develop targeted purification train.*$",
     r"Zielgerichteten Aufreinigungs-Zug (Kristallisation, HPLC, Polishing) entwickeln, um Reinheit bei minimalen Kosten zu erreichen"),

    # ---- Headlines / labels ----
    (r"^Biotechnological Process$", r"Biotechnologischer Prozess"),
    (r"^Chemical Process$", r"Chemischer Prozess"),
    (r"^Additional adjustments based on recognized process class$",
     r"Zusätzliche Anpassungen basierend auf erkannter Prozessklasse"),

    # ---- Bug B1 i18n fix: improvement strings in the decision engine ----
    # These are generated when cost == "high"/"very high", risk == "high",
    # strict_waste etc.; they bypassed the existing H7 phrase rules because
    # they end with "Potential cost reduction: <impact>" / "Erwartete
    # Kosten-Reduzierung: <impact>" suffix that wasn't anchored.
    (r"^Reduce number of synthesis steps to improve efficiency and lower cost.*$",
     r"Schrittanzahl reduzieren, um Effizienz zu steigern und Kosten zu senken — Erwartete Kosten-Reduzierung: MITTEL bis HOCH"),
    (r"^Evaluate alternative suppliers for high-cost precursors and reduce protecting-group use to lower COGS.*$",
     r"Alternative Lieferanten für teure Vorprodukte prüfen und Schutzgruppen-Einsatz reduzieren, um die Herstellkosten (COGS) zu senken — Erwartete Kosten-Reduzierung: MITTEL"),
    (r"^Implement process control improvements \(robust pH/temperature control, inline monitoring\) and simplify critical unit ops to reduce failure modes.*$",
     r"Prozesskontrolle verbessern (robuste pH-/Temperatur-Regelung, Inline-Monitoring) und kritische Unit-Operations vereinfachen, um Fehlerquellen zu reduzieren — Erwartete Kosten-Reduzierung: MITTEL"),
    (r"^Switch to less toxic reagents or reduce hazardous waste streams.*$",
     r"Auf weniger toxische Reagenzien umsteigen oder gefährliche Abfallströme reduzieren — Erwartete Kosten-Reduzierung: MITTEL"),

    # ---- Bug H7 i18n fix: English engine strings flagged by deep audit ----
    (r"^Distillation likely feasible due to volatility; consider azeotrope management and reflux optimization\.?$",
     r"Destillation aufgrund der Flüchtigkeit gut machbar — Azeotrop-Management und Rückfluss-Optimierung berücksichtigen"),
    (r"^Extraction processes (?:may have|is likely to have) low yield and high variability\.?$",
     r"Extraktionsprozesse haben oft niedrige Ausbeute und hohe Variabilität"),
    (r"^Reducing steps may impact yield or product quality and requires route redesign\.?$",
     r"Schrittreduktion kann Ausbeute oder Produktqualität beeinträchtigen und erfordert ein Route-Redesign"),
    (r"^Protein downstreams require chromatography-based capture and polishing; ensure robust column regeneration and buffer systems\.?$",
     r"Protein-Aufarbeitung erfordert chromatographisches Capture und Polishing — robuste Säulen-Regeneration und Puffer-Systeme sicherstellen"),
    (r"^For large biomolecules, structural representation has limited predictive power compared to biophysical behavior\.?$",
     r"Für große Biomoleküle hat die strukturelle Darstellung begrenzte Vorhersagekraft gegenüber dem biophysikalischen Verhalten"),
    # extra defensive variants (subject-verb mismatch seen in audit: "processes is")
    (r"^Extraction processes is likely to have low yield and high variability\.?$",
     r"Extraktionsprozesse haben oft niedrige Ausbeute und hohe Variabilität"),
]


def _apply_phrase_rules(text: str) -> str:
    """Apply the first matching regex rule and return the substituted string.
    If no rule matches, return text unchanged.
    """
    for pattern, replacement in _PHRASE_RULES:
        try:
            new = re.sub(pattern, replacement, text, flags=re.IGNORECASE)
            if new != text:
                return new
        except re.error:
            continue
    return text


# ---------------------------------------------------------------------------
# Inline word translations (applied as a final pass for short words inside
# longer paragraphs that didn't match any phrase rule).
# ---------------------------------------------------------------------------

_WORD_RULES_DE = [
    # (regex, replacement). Word-boundaries protect against partial hits.
    (re.compile(r"\bHIGH\b"), "HOCH"),
    (re.compile(r"\bLOW\b"), "NIEDRIG"),
    (re.compile(r"\bMEDIUM\b"), "MITTEL"),
    (re.compile(r"\bVERY HIGH\b"), "SEHR HOCH"),
    # only translate lowercase 'high/low/medium' that follow ': ' or are wrapped — safer than a global swap
    (re.compile(r":\s+high(?=[\s.,)])"), ": hoch"),
    (re.compile(r":\s+low(?=[\s.,)])"), ": niedrig"),
    (re.compile(r":\s+medium(?=[\s.,)])"), ": mittel"),
]


def translate_engine_text(text: Optional[str], lang: str = "de") -> str:
    """Translate a dynamic engine string to the requested language.

    Strategy:
      - 'en' → return as-is (engine already produces English).
      - 'de' → try phrase-level rules first; if no full match, apply
        word-level substitutions for HIGH/LOW/MEDIUM. Anything else stays in
        English (graceful fallback — better English than a crash).
    """
    if text is None:
        return ""
    s = str(text)
    if lang != "de":
        return s
    # Phrase pass
    new = _apply_phrase_rules(s)
    # Word-level pass (always — useful for embedded HIGH/LOW etc.)
    for rx, repl in _WORD_RULES_DE:
        new = rx.sub(repl, new)
    return new
