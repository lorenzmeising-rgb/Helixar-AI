"""Lightweight i18n / language helper for Helixar AI.

Design:
  - Two supported languages: 'de' (default) and 'en'.
  - Centralized TRANSLATIONS dict — one place to add new strings.
  - t(key, lang) returns the translated string, falling back to the key itself
    when no translation is registered (so an unfinished migration won't crash).
  - get_lang() reads the current language from Streamlit session_state when
    available, else falls back to 'de'.

Usage:
    from i18n import t
    st.subheader(t("section_2_title"))

    # Or pass a language explicitly (e.g. when generating a PDF):
    title = t("section_2_title", lang="en")

Migration strategy: strings can be moved into this dict incrementally.
Anything still hard-coded in DE just stays DE — no breakage.
"""

from typing import Optional


SUPPORTED_LANGUAGES = ("de", "en")
DEFAULT_LANGUAGE = "de"


TRANSLATIONS = {
    # ====== Sprach-Toggle ======
    "lang_label": {
        "de": "Sprache",
        "en": "Language",
    },
    "lang_de": {"de": "Deutsch", "en": "German"},
    "lang_en": {"de": "Englisch", "en": "English"},

    # ====== Landing / Header ======
    "app_title": {
        "de": "Helixar · Produktions-Blueprint-Generator",
        "en": "Helixar · Production Blueprint Generator",
    },
    "app_caption": {
        "de": "Regelbasierte Entscheidungsunterstützung für die Prozess-Optimierung (nicht-operativ).",
        "en": "Rule-based decision support for process optimization (non-operational).",
    },
    "describe_process_header": {
        "de": "Beschreiben Sie Ihren Prozess",
        "en": "Describe Your Process",
    },
    "describe_process_intro": {
        "de": "Beschreiben Sie Ihren bestehenden Produktionsprozess, damit wir ihn analysieren und Verbesserungen vorschlagen können.",
        "en": "Describe your existing production process so we can analyze it and suggest improvements.",
    },

    # ====== Molekül-Block (oben) ======
    "target_molecule": {"de": "Zielmolekül", "en": "Target Molecule"},
    "molecule_type": {"de": "Molekültyp", "en": "Molecule Type"},
    "molecule_subtype": {"de": "Molekül-Untertyp (optional)", "en": "Molecule Subtype (optional)"},

    "type_small_molecule": {"de": "Kleines Molekül", "en": "Small molecule"},
    "type_peptide": {"de": "Peptid", "en": "Peptide"},
    "type_protein": {"de": "Protein", "en": "Protein"},
    "type_natural_product": {"de": "Naturstoff", "en": "Natural product"},

    "subtype_none": {"de": "Keiner", "en": "None"},
    "subtype_volatile": {"de": "Flüchtig", "en": "Volatile"},
    "subtype_non_volatile": {"de": "Nicht flüchtig", "en": "Non-volatile"},
    "subtype_linear": {"de": "Linear", "en": "Linear"},
    "subtype_cyclic": {"de": "Zyklisch", "en": "Cyclic"},
    "subtype_antibody": {"de": "Antikörper", "en": "Antibody"},
    "subtype_enzyme": {"de": "Enzym", "en": "Enzyme"},
    "subtype_terpene": {"de": "Terpen", "en": "Terpene"},
    "subtype_alkaloid": {"de": "Alkaloid", "en": "Alkaloid"},

    "molecule_type_help": {
        "de": ("Kleines Molekül: organische Verbindungen < 900 Da (z. B. Vanillin, Aspirin). "
               "Peptid: Aminosäureketten typ. < 50 Reste (z. B. Glutathion, Cyclosporin). "
               "Protein: > 50 Aminosäuren, gefaltete 3D-Struktur (z. B. Antikörper, Enzyme). "
               "Naturstoff: aus Pflanzen/Mikroben gewonnene Sekundärmetaboliten (z. B. Terpene, Alkaloide)."),
        "en": ("Small molecule: organic compounds < 900 Da (e.g. vanillin, aspirin). "
               "Peptide: amino-acid chains typically < 50 residues (e.g. glutathione, cyclosporine). "
               "Protein: > 50 amino acids, folded 3D structure (e.g. antibodies, enzymes). "
               "Natural product: secondary metabolites from plants/microbes (e.g. terpenes, alkaloids)."),
    },
    "molecule_subtype_help": {
        "de": ("Optionale Verfeinerung des Molekültyps. Beeinflusst Reinigungs- und Skalierungsempfehlungen "
               "(z. B. Antikörper benötigen Protein-A-Affinität, zyklische Peptide haben höhere Stabilität als lineare)."),
        "en": ("Optional refinement of the molecule type. Influences purification and scale-up recommendations "
               "(e.g. antibodies require protein-A affinity, cyclic peptides are more stable than linear ones)."),
    },

    # ====== Mode-Switch ======
    "no_existing_process_label": {
        "de": "Ich habe noch keinen Produktionsprozess — schlage mir einen vor",
        "en": "I do not have a production process yet — suggest one for me",
    },
    "no_existing_process_help": {
        "de": ("Wenn aktiv, beschreiben Sie nur das Zielprodukt, die Anforderungen und Ihr verfügbares Equipment. "
               "Helixar entwirft dann einen Prozess von Grund auf, statt einen bestehenden zu analysieren."),
        "en": ("If enabled, you only describe the target product, requirements, and available equipment. "
               "Helixar then drafts a process from scratch instead of analyzing an existing one."),
    },
    "greenfield_info": {
        "de": ("Greenfield-Modus: Helixar wählt eine geeignete Produktionsmethode und schätzt die Anzahl Schritte "
               "selbst, basierend auf Molekül, Anforderungen und Equipment."),
        "en": ("Greenfield mode: Helixar will pick a suitable production method and estimate the number of steps, "
               "based on the molecule, requirements, and equipment."),
    },
    "greenfield_market_info": {
        "de": ("Greenfield-Modus: Marktbedingungen (Lieferanten, Edukt-Preis, Lieferzeit) sind ohne konkrete "
               "Prozess-Route nicht bekannt — Helixar nutzt Industrie-Benchmark-Werte für deine Molekül-Klasse "
               "als Schätzung. Sobald du eine Prozess-Route festlegst, kannst du die konkreten Werte "
               "eingeben für eine präzisere Analyse."),
        "en": ("Greenfield mode: market conditions (suppliers, raw-material price, lead time) are unknown without "
               "a concrete process route — Helixar uses industry-benchmark values for your molecule class as an "
               "estimate. Once you commit to a route, you can enter the actual values for a more precise analysis."),
    },

    # ====== Section 1 — Aktueller Prozess ======
    "section_1_title": {"de": "1. Aktueller Prozess", "en": "1. Current Process"},
    "current_method_label": {"de": "Aktuelle Produktionsmethode", "en": "Current Production Method"},
    "current_method_help": {
        "de": ("chemisch: klassische organische Synthese (Reagenzien, Lösungsmittel). "
               "biotechnologisch: Fermentation oder rekombinante Expression in Mikroorganismen/Zellkultur. "
               "Extraktion: Isolation aus natürlichen Quellen (Pflanzen, Mikroben, Tiere)."),
        "en": ("chemical: classical organic synthesis (reagents, solvents). "
               "biotechnological: fermentation or recombinant expression in microorganisms/cell culture. "
               "extraction: isolation from natural sources (plants, microbes, animals)."),
    },
    "method_chemical": {"de": "chemisch", "en": "chemical"},
    "method_biotechnological": {"de": "biotechnologisch", "en": "biotechnological"},
    "method_extraction": {"de": "Extraktion", "en": "extraction"},
    "number_of_steps_label": {"de": "Anzahl Prozessschritte", "en": "Number of Process Steps"},
    "number_of_steps_help": {
        "de": ("Gesamtzahl der Prozessschritte vom Edukt/Rohmaterial zum Zielprodukt — "
               "inklusive Synthese-, Fermentations-, Extraktions- und Reinigungsschritten. "
               "Beispiele: Aspirin chemisch = ca. 2 Schritte (Acetylierung + Umkristallisation); "
               "Vanillin biotech = ca. 5 Schritte (Vorkultur + Fermentation + 3 DSP-Stufen); "
               "Trastuzumab biotech = ca. 6–8 Schritte (Fermentation + Capture + Polishing + Polishing-Sequenz)."),
        "en": ("Total number of process steps from starting material to target product — "
               "including synthesis, fermentation, extraction and purification steps. "
               "Examples: aspirin chemical = ca. 2 steps (acetylation + recrystallization); "
               "vanillin biotech = ca. 5 steps (pre-culture + fermentation + 3 DSP stages); "
               "trastuzumab biotech = ca. 6–8 steps (fermentation + capture + polishing sequence)."),
    },

    # ====== Section 2 — Anforderungen ======
    "section_2_title": {"de": "2. Anforderungen an den Prozess", "en": "2. Process Requirements"},
    "purity_label": {"de": "Zielreinheit in % (Reinheitsanforderung)", "en": "Target purity in % (Purity Requirement)"},
    "purity_help": {
        "de": ("Geforderte Reinheit des Endprodukts in Prozent. "
               "Orientierungswerte: ~95 % technische Qualität, ≥ 98 % Lebensmittel-/Kosmetikqualität, "
               "≥ 99 % pharmazeutische / API-Qualität, ≥ 99,9 % klinisch."),
        "en": ("Target purity of the final product in percent. "
               "Reference: ~95 % technical grade, ≥ 98 % food/cosmetic grade, "
               "≥ 99 % pharmaceutical / API grade, ≥ 99.9 % clinical."),
    },
    "scale_label": {"de": "Zielmaßstab in kg/Jahr (Produktionsmaßstab)", "en": "Target scale in kg/year (Production Scale)"},
    "scale_help": {
        "de": ("Geplante Produktionsmenge pro Jahr in kg. "
               "Orientierung: < 1 kg/a Labormaßstab; 1–1.000 kg/a Pilot; "
               "≥ 1.000 kg/a industrieller Maßstab; ≥ 100.000 kg/a Commodity."),
        "en": ("Planned production volume per year in kg. "
               "Reference: < 1 kg/y lab scale; 1–1,000 kg/y pilot; "
               "≥ 1,000 kg/y industrial scale; ≥ 100,000 kg/y commodity."),
    },

    # ====== Section 3 — Marktbedingungen ======
    "section_3_title": {"de": "3. Marktbedingungen (Ist-Zustand)", "en": "3. Market Conditions (Current State)"},
    "section_3_caption": {
        "de": ("Sourcing-Profil der wichtigsten Edukte. Die Kombination aus Lieferanten-Anzahl, Lieferzeit "
               "und geografischer Konzentration bestimmt das Versorgungsrisiko."),
        "en": ("Sourcing profile of the key starting materials. The combination of supplier count, lead time, "
               "and geographic concentration determines supply risk."),
    },
    "num_suppliers_label": {"de": "Anzahl qualifizierter Lieferanten", "en": "Number of qualified suppliers"},
    "num_suppliers_help": {
        "de": ("Anzahl unterschiedlicher Lieferanten, die das Edukt in geforderter Qualität liefern können. "
               "1 = Single-Source-Risiko; ≥ 4 = robuste Mehrfach-Quellen."),
        "en": ("Number of distinct suppliers able to deliver the starting material at required quality. "
               "1 = single-source risk; ≥ 4 = robust multi-source."),
    },
    "lead_time_label": {"de": "Typische Lieferzeit in Wochen", "en": "Typical lead time in weeks"},
    "lead_time_help": {
        "de": ("Durchschnittliche Lead-Time vom Bestellen bis zum Wareneingang. "
               "Orientierung: ≤ 2 Wochen Commodity / Lagerware; 2–8 Wochen Standard-Reagenzien; "
               "> 8 Wochen Spezialchemikalien oder Auftragssynthese."),
        "en": ("Average lead time from order to delivery. "
               "Reference: ≤ 2 weeks commodity / stock; 2–8 weeks standard reagents; "
               "> 8 weeks specialty chemicals or custom synthesis."),
    },
    "single_region_label": {
        "de": "Lieferanten konzentriert in einer Region (z. B. nur China oder nur EU)",
        "en": "Suppliers concentrated in one region (e.g. only China or only EU)",
    },
    "single_region_help": {
        "de": ("Markieren, wenn alle qualifizierten Lieferanten geografisch in einer Region sitzen. "
               "Dies erhöht das Geopolitik- und Logistikrisiko (z. B. Export-Restriktionen, regionale Produktionsausfälle)."),
        "en": ("Check if all qualified suppliers are located in a single region. "
               "This increases geopolitical and logistics risk (e.g. export restrictions, regional production outages)."),
    },
    "raw_cost_label": {"de": "Aktuelle Marktpreise der Rohstoffe in €/kg", "en": "Current market price of raw materials in €/kg"},
    "raw_cost_help": {
        "de": ("Durchschnittlicher Marktpreis der wichtigsten Edukte in € pro kg. "
               "Orientierung: < 10 €/kg Bulk-Chemikalien (Zucker, einfache Salze); "
               "10–500 €/kg gängige Reagenzien und Lösungsmittel; "
               "> 500 €/kg Spezialreagenzien, chirale Bausteine, seltene Enzyme."),
        "en": ("Average market price of key starting materials in € per kg. "
               "Reference: < 10 €/kg bulk chemicals (sugar, basic salts); "
               "10–500 €/kg common reagents and solvents; "
               "> 500 €/kg specialty reagents, chiral building blocks, rare enzymes."),
    },
    "strict_waste_label": {
        "de": "Strenge Auflagen für Abfall / Toxizität",
        "en": "Strict waste / toxicity constraints",
    },
    "strict_waste_help": {
        "de": ("Gibt es strenge Auflagen für den Umgang mit Abfällen / toxischen Lösungsmitteln? "
               "z. B. halogenierte Lösungsmittel, Schwermetalle, GMO-Abfälle, GMP-Pflicht zur Aufbereitung."),
        "en": ("Are there strict constraints on handling waste / toxic solvents? "
               "e.g. halogenated solvents, heavy metals, GMO waste, GMP-mandated treatment."),
    },

    # ====== Section 4 — Verfügbare Methodiken ======
    "section_4_title": {
        "de": "4. Verfügbare Methodiken & Equipment im Labor",
        "en": "4. Available Methods & Lab Equipment",
    },
    "section_4_caption": {
        "de": ("Geben Sie an, welche Geräte und Methoden grundsätzlich verfügbar sind — "
               "unabhängig davon, ob sie aktuell genutzt werden."),
        "en": ("Indicate which equipment and methods are available in principle — "
               "regardless of whether they are currently in use."),
    },
    "bioreactor_label": {
        "de": "Bioreaktor / Fermenter (≥ 1 L, mit pH/Temperatur-/O2-Kontrolle)",
        "en": "Bioreactor / fermenter (≥ 1 L, with pH/temperature/O2 control)",
    },
    "bioreactor_help": {
        "de": "Voraussetzung für skalierbare biotechnologische Produktion.",
        "en": "Prerequisite for scalable biotechnological production.",
    },
    "purification_methods_header": {
        "de": "**Reinigungs- und Aufarbeitungsmethoden:**",
        "en": "**Purification and downstream methods:**",
    },
    "method_flash_chrom_label": {
        "de": "Flash-Chromatographie (Säule, SiO2/RP)",
        "en": "Flash chromatography (column, SiO2/RP)",
    },
    "method_flash_chrom_help": {
        "de": "Standard-Reinigung für kleine Moleküle nach Synthese; Trennung über Polaritätsgradienten.",
        "en": "Standard purification for small molecules after synthesis; separation via polarity gradients.",
    },
    "method_prep_hplc_label": {"de": "Präparative HPLC", "en": "Preparative HPLC"},
    "method_prep_hplc_help": {
        "de": "Hochauflösende Reinigung für sensible/wertvolle Produkte (Peptide, API-Grade kleine Moleküle).",
        "en": "High-resolution purification for sensitive/valuable products (peptides, API-grade small molecules).",
    },
    "method_rotavap_label": {"de": "Rotationsverdampfer", "en": "Rotary evaporator"},
    "method_rotavap_help": {
        "de": "Lösungsmittelentfernung unter Vakuum bei moderater Temperatur.",
        "en": "Solvent removal under vacuum at moderate temperature.",
    },
    "method_lyo_label": {"de": "Lyophilisator (Gefriertrocknung)", "en": "Lyophilizer (freeze-drying)"},
    "method_lyo_help": {
        "de": "Schonende Trocknung für thermolabile Substanzen, Peptide, Proteine.",
        "en": "Gentle drying for thermolabile substances, peptides, proteins.",
    },
    "method_cryst_label": {"de": "Kristallisation / Umkristallisation", "en": "Crystallization / recrystallization"},
    "method_cryst_help": {
        "de": "Klassische Reinigungs- und Endformulierungsmethode für stabile kleine Moleküle.",
        "en": "Classical purification and final formulation method for stable small molecules.",
    },
    "method_membrane_label": {"de": "Membranfiltration / Ultrafiltration", "en": "Membrane filtration / ultrafiltration"},
    "method_membrane_help": {
        "de": "Größenbasierte Aufkonzentrierung und Pufferaustausch für Proteine/Peptide.",
        "en": "Size-based concentration and buffer exchange for proteins/peptides.",
    },
    "method_fplc_label": {"de": "FPLC (Protein-Chromatographie)", "en": "FPLC (protein chromatography)"},
    "method_fplc_help": {
        "de": "Affinitäts-, Ionenaustausch- oder Größenausschluss-Chromatographie für Proteine.",
        "en": "Affinity, ion-exchange, or size-exclusion chromatography for proteins.",
    },
    "method_extraction_label": {"de": "Liquid-Liquid-Extraktion / Soxhlet", "en": "Liquid-liquid extraction / Soxhlet"},
    "method_extraction_help": {
        "de": "Trennung über Verteilung zwischen zwei nicht mischbaren Phasen; auch für Naturstoffe.",
        "en": "Separation via partition between two immiscible phases; also for natural products.",
    },

    # ====== Section 5 — Biophysikalisch ======
    "section_5_title": {
        "de": "Biophysikalische Eigenschaften (nur für Peptide & Proteine)",
        "en": "Biophysical properties (only for peptides & proteins)",
    },
    "aggregation_label": {"de": "Aggregations-Risiko", "en": "Aggregation risk"},
    "aggregation_help": {
        "de": ("Neigung zur Bildung von Oligomeren oder Aggregaten unter Prozessbedingungen. "
               "low: stabil bis hohe Konzentrationen (> 50 mg/mL). "
               "medium: gelegentliche Aggregation, Formulierung beherrschbar. "
               "high: starke Aggregation, häufig bei hydrophoben Sequenzen / Antikörpern mit ungünstigem pI."),
        "en": ("Tendency to form oligomers or aggregates under process conditions. "
               "low: stable up to high concentrations (> 50 mg/mL). "
               "medium: occasional aggregation, manageable in formulation. "
               "high: strong aggregation, common in hydrophobic sequences / antibodies with unfavourable pI."),
    },
    "folding_label": {"de": "Faltungs-Komplexität", "en": "Folding complexity"},
    "folding_help": {
        "de": ("Komplexität der korrekten 3D-Faltung. "
               "low: kurze, lineare Peptide ohne Disulfide. "
               "medium: 1–2 Disulfidbrücken, einfache Domänen. "
               "high: mehrere Disulfide, Multi-Domänen, posttranslationale Modifikationen (Glykosylierung), Antikörper."),
        "en": ("Complexity of correct 3D folding. "
               "low: short, linear peptides without disulfides. "
               "medium: 1–2 disulfide bridges, simple domains. "
               "high: multiple disulfides, multi-domain, post-translational modifications (glycosylation), antibodies."),
    },
    "stability_label": {"de": "Biophysikalische Stabilität", "en": "Biophysical stability"},
    "stability_help": {
        "de": ("Stabilität gegenüber Temperatur, pH-Schwankungen, Scherung und Lagerung. "
               "low: Tm < 50 °C, schnelle Degradation, kühlpflichtig. "
               "medium: Tm 50–70 °C, lagerfähig bei 4 °C. "
               "high: Tm > 70 °C, raumtemperaturstabil, für Lyophilisation geeignet."),
        "en": ("Stability against temperature, pH fluctuations, shear, and storage. "
               "low: Tm < 50 °C, rapid degradation, requires cold chain. "
               "medium: Tm 50–70 °C, storable at 4 °C. "
               "high: Tm > 70 °C, room-temperature stable, suitable for lyophilization."),
    },

    # ====== Numeric biophysical inputs (replaces qualitative dropdowns) ======
    "tagg_label": {"de": "Aggregations-Onset Tagg in °C", "en": "Aggregation onset Tagg in °C"},
    "tagg_help": {
        "de": ("Temperatur, ab der das Protein/Peptid aggregiert (typisch aus DSF-Messung). "
               "Orientierung: ≥ 65 °C niedriges Risiko (Antikörper-Standard); "
               "50–65 °C moderates Risiko; < 50 °C hohes Risiko, Formulation-Engineering nötig."),
        "en": ("Temperature at which the protein/peptide starts to aggregate (typically from DSF). "
               "Reference: ≥ 65 °C low risk (antibody standard); "
               "50–65 °C moderate risk; < 50 °C high risk, formulation engineering required."),
    },
    "disulfides_label": {"de": "Anzahl Disulfidbrücken", "en": "Number of disulfide bridges"},
    "disulfides_help": {
        "de": ("Anzahl der intra- und intermolekularen Disulfidbrücken (S–S). "
               "Beispiele: lineare Peptide 0; Insulin 3; Antikörper (IgG) 16."),
        "en": ("Number of intra- and intermolecular disulfide bridges (S–S). "
               "Examples: linear peptides 0; insulin 3; antibodies (IgG) 16."),
    },
    "domains_label": {"de": "Anzahl Domänen", "en": "Number of domains"},
    "domains_help": {
        "de": ("Strukturell und funktionell unabhängige Untereinheiten. "
               "Beispiele: kleine Peptide / Hormone 1; einfache Enzyme 1–2; "
               "Antikörper (IgG) 12; Multi-Domain-Fusionsproteine 3+."),
        "en": ("Structurally and functionally independent subunits. "
               "Examples: small peptides / hormones 1; simple enzymes 1–2; "
               "antibodies (IgG) 12; multi-domain fusion proteins 3+."),
    },
    "ptm_label": {
        "de": "Posttranslationale Modifikationen (z. B. Glykosylierung, Phosphorylierung)",
        "en": "Post-translational modifications (e.g. glycosylation, phosphorylation)",
    },
    "ptm_help": {
        "de": ("Markieren, wenn das Molekül PTMs trägt. Erfordert eukaryotische Expression "
               "(CHO, HEK293, Pichia) — E. coli kann die meisten PTMs nicht durchführen."),
        "en": ("Check if the molecule carries PTMs. Requires eukaryotic expression "
               "(CHO, HEK293, Pichia) — E. coli cannot perform most PTMs."),
    },
    "tm_label": {"de": "Schmelzpunkt Tm in °C", "en": "Melting temperature Tm in °C"},
    "tm_help": {
        "de": ("Thermische Entfaltungstemperatur (typisch aus DSC oder DSF). "
               "Orientierung: ≥ 70 °C raumtemperaturstabil, lyophilisierbar; "
               "50–70 °C 4 °C-Lagerung; < 50 °C Tiefkühl-/Cryoprotectant nötig."),
        "en": ("Thermal unfolding temperature (typically from DSC or DSF). "
               "Reference: ≥ 70 °C room-temperature stable, lyophilizable; "
               "50–70 °C 4 °C storage; < 50 °C deep-freeze / cryoprotectant required."),
    },

    # Display labels for the report (numeric input echoes)
    "metric_tagg": {"de": "Aggregations-Onset Tagg", "en": "Aggregation onset Tagg"},
    "metric_tm": {"de": "Schmelzpunkt Tm", "en": "Melting temperature Tm"},
    "metric_disulfides": {"de": "Disulfidbrücken", "en": "Disulfide bridges"},
    "metric_domains": {"de": "Domänen", "en": "Domains"},
    "metric_ptm": {"de": "Posttranslationale Modifikationen", "en": "Post-translational modifications"},
    "biophysical_inputs_header": {"de": "Biophysikalisches Profil", "en": "Biophysical profile"},

    # ====== Buttons / Actions ======
    "analyze_button": {"de": "Prozess analysieren", "en": "Analyze process"},
    "back_to_start": {"de": "Zurück zur Startseite", "en": "Back to start"},
    "preview_button": {"de": "Vorschau anzeigen", "en": "Show preview"},
    "export_pdf_button": {"de": "Als PDF exportieren", "en": "Export as PDF"},
    "download_pdf_button": {"de": "PDF herunterladen", "en": "Download PDF"},

    # ====== Validierungs-Meldungen ======
    "warn_target_required": {
        "de": "Bitte geben Sie das Zielmolekül an (Pflichtfeld).",
        "en": "Please provide the target molecule (required).",
    },
    "warn_steps_required": {
        "de": "Bitte geben Sie die Anzahl Prozessschritte an (Pflichtfeld).",
        "en": "Please provide the number of process steps (required).",
    },
    "error_protein_chemical": {
        "de": ("Unmögliche Kombination: Proteine können nicht chemisch synthetisiert werden. "
               "Proteine werden rekombinant in lebenden Zellen exprimiert (z. B. CHO, E. coli, Hefe). "
               "Bitte wählen Sie 'biotechnologisch' als Produktionsmethode. "
               "Hinweis: Peptide (kürzere Aminosäureketten) können chemisch über SPPS synthetisiert werden."),
        "en": ("Impossible combination: proteins cannot be synthesized chemically. "
               "Proteins are expressed recombinantly in living cells (e.g. CHO, E. coli, yeast). "
               "Please choose 'biotechnological' as production method. "
               "Note: peptides (shorter amino-acid chains) can be synthesized chemically via SPPS."),
    },
    "spinner_analyzing": {
        "de": "Prozess wird analysiert …",
        "en": "Analyzing process …",
    },
    "warn_no_results": {
        "de": "Analyse konnte keine Ergebnisse erzeugen. Bitte prüfen Sie die Eingaben.",
        "en": "Analysis returned no results. Please review the inputs.",
    },

    # ====== Plausibility / Recommendations Header ======
    "plausibility_header": {
        "de": "**Plausibilitätshinweise zu Ihren Eingaben:**",
        "en": "**Plausibility notes on your inputs:**",
    },
    "plausibility_continue": {
        "de": "Die Analyse wird trotzdem fortgesetzt — bitte prüfen Sie ggf. die Eingaben.",
        "en": "The analysis will continue regardless — please review the inputs if appropriate.",
    },

    # ====== Report-Seite ======
    "report_success_title": {
        "de": "Erfolgreiche Produktionsempfehlung für {compound} generiert",
        "en": "Successfully generated production recommendation for {compound}",
    },
    "your_inputs_header": {"de": "Ihre Eingaben", "en": "Your inputs"},
    # Feature D — what-if sensitivity panel labels
    "sensitivity_section_title": {
        "de": "🎚 Was-wäre-wenn: Sensitivitäts-Analyse",
        "en": "🎚 What-if: Sensitivity analysis",
    },
    "sensitivity_section_caption": {
        "de": "Verschieben Sie Skala, Reinheit oder Rohstoffpreis und sehen Sie sofort, wie sich die COGS-Schätzung verändert. Das Modell rechnet live — kein neuer PDF-Lauf nötig.",
        "en": "Move scale, purity or raw-material price and watch the COGS estimate update live. The model recomputes on the fly — no new PDF run needed.",
    },
    "sensitivity_not_available": {
        "de": "Sensitivitäts-Modul nicht verfügbar (cogs_estimator konnte nicht geladen werden).",
        "en": "Sensitivity module unavailable (cogs_estimator could not be loaded).",
    },
    "sensitivity_scale_label": {
        "de": "Skala (log10 kg/Jahr)",
        "en": "Scale (log10 kg/year)",
    },
    "sensitivity_scale_help": {
        "de": "Slider arbeitet logarithmisch: −3 = 1 g/Jahr, 0 = 1 kg/Jahr, 3 = 1 t/Jahr, 6 = 1 000 t/Jahr.",
        "en": "Slider is log-scale: −3 = 1 g/year, 0 = 1 kg/year, 3 = 1 t/year, 6 = 1 000 t/year.",
    },
    "sensitivity_purity_label": {
        "de": "Zielreinheit (%)",
        "en": "Target purity (%)",
    },
    "sensitivity_purity_help": {
        "de": "Höhere Reinheit erhöht typ. die Aufreinigungskosten (×1,15 bei 99 %, ×1,3 bei 99,5 %, ×1,6 bei 99,9 %).",
        "en": "Higher purity typically increases purification cost (×1.15 at 99 %, ×1.3 at 99.5 %, ×1.6 at 99.9 %).",
    },
    "sensitivity_rm_label": {
        "de": "Rohstoffpreis (€/kg)",
        "en": "Raw-material price (€/kg)",
    },
    "sensitivity_rm_help": {
        "de": "Bei Rohstoff-Anteil ≥30 % der COGS schlägt eine 20 %-Preisänderung mit ≈6–10 % auf die Stückkosten durch.",
        "en": "If raw materials are ≥ 30 % of COGS, a 20 % price change translates to ≈ 6–10 % on unit cost.",
    },
    "sensitivity_original_label": {
        "de": "Original-Eingabe",
        "en": "Original input",
    },
    "sensitivity_simulated_label": {
        "de": "Simuliert",
        "en": "Simulated",
    },
    "sensitivity_chart_title": {
        "de": "COGS-Sensitivität gegenüber Skala",
        "en": "COGS sensitivity vs. scale",
    },
    "sensitivity_chart_caption": {
        "de": "Linke Y-Achse: COGS-Bandbreite (€/kg). X-Achse: Jahres-Maßstab in kg, logarithmisch.",
        "en": "Left Y-axis: COGS range (€/kg). X-axis: annual scale in kg, log scale.",
    },
    "sensitivity_reset_button": {
        "de": "🔄 Zurück zu Original-Eingabe",
        "en": "🔄 Reset to original input",
    },
    "metric_purity": {"de": "Zielreinheit", "en": "Target purity"},
    "metric_scale": {"de": "Zielmaßstab", "en": "Target scale"},
    "metric_cost": {"de": "Marktpreis Edukte", "en": "Raw-material price"},
    "metric_suppliers": {"de": "Qualifizierte Lieferanten", "en": "Qualified suppliers"},
    "metric_lead_time": {"de": "Lieferzeit", "en": "Lead time"},
    "metric_geo_concentration": {"de": "Geo-Konzentration", "en": "Geo concentration"},
    "yes": {"de": "ja", "en": "yes"},
    "no": {"de": "nein", "en": "no"},
    "weeks": {"de": "Wochen", "en": "weeks"},
    "kg_per_year": {"de": "kg/Jahr", "en": "kg/year"},

    "concrete_recs_header": {"de": "Konkrete Empfehlungen", "en": "Concrete recommendations"},
    "concrete_recs_route": {"de": "Produktionsroute", "en": "Production route"},
    "concrete_recs_yield": {"de": "Erwartete Ausbeute", "en": "Expected yield"},
    "concrete_recs_time": {"de": "Prozessdauer", "en": "Process duration"},
    "concrete_recs_purity": {"de": "Erreichbare Reinheit nach Aufarbeitung", "en": "Achievable purity after workup"},
    "concrete_recs_steps": {"de": "Empfohlene Aufarbeitungsschritte", "en": "Recommended downstream steps"},
    "concrete_recs_notes": {"de": "Hinweise", "en": "Notes"},

    "actions_header": {"de": "Empfohlene Optimierungen", "en": "Recommended optimizations"},
    "actions_caption": {
        "de": "Konkrete Maßnahmen, mit denen sich der Prozess voraussichtlich verbessern lässt.",
        "en": "Concrete measures that are likely to improve the process.",
    },
    "actions_rationale": {"de": "Begründung", "en": "Rationale"},
    "actions_impact": {"de": "Erwartete Wirkung", "en": "Expected impact"},
    "actions_prerequisites": {"de": "Voraussetzungen", "en": "Prerequisites"},
    "actions_effort": {"de": "Aufwand", "en": "Effort"},
    "actions_current_state": {"de": "Aktueller Stand", "en": "Current state"},
    "actions_optimized_state": {"de": "Optimierter Stand", "en": "Optimized state"},
    "effort_low": {"de": "niedrig", "en": "low"},
    "effort_medium": {"de": "mittel", "en": "medium"},
    "effort_high": {"de": "hoch", "en": "high"},

    # ====== PDF-Strings ======
    "pdf_title_default": {"de": "Helixar Produktionsbericht", "en": "Helixar Production Report"},
    "pdf_subtitle": {"de": "Prozessoptimierung — Zusammenfassung", "en": "Process Optimization — Executive Summary"},
    "pdf_executive_summary": {"de": "Zusammenfassung", "en": "Executive Summary"},
    "pdf_efficiency": {"de": "Effizienz", "en": "Efficiency"},
    "pdf_cost": {"de": "Kosten", "en": "Cost"},
    "pdf_risk": {"de": "Risiko", "en": "Risk"},
    # Bug M2 fix: row aggregates intrinsic toxicity AND process containment
    # / QC requirements (endotoxin spec, controlled-substance handling).
    # The previous "Toxizität" label was misread as PRODUCT toxicity for
    # Met-enkephalin / Liraglutide where the driver is actually
    # endotoxin/QC.
    "pdf_toxicity": {"de": "Stoff- & QC-Risiko", "en": "Substance & QC risk"},
    "pdf_confidence": {"de": "Confidence Level", "en": "Confidence level"},
    "pdf_key_issues": {"de": "Wichtige Themen", "en": "Key issues"},
    "pdf_recommended_improvements": {"de": "Empfohlene Verbesserungen", "en": "Recommended improvements"},
    "pdf_plausibility_section": {"de": "Plausibilitätshinweise zu den Eingaben", "en": "Plausibility notes on inputs"},
    "pdf_process_analysis": {"de": "Prozessanalyse", "en": "Process Analysis"},
    "pdf_quantitative_inputs": {"de": "Eingaben (quantitativ)", "en": "Inputs (quantitative)"},
    "pdf_target_purity": {"de": "Zielreinheit", "en": "Target purity"},
    "pdf_target_scale": {"de": "Zielmaßstab", "en": "Target scale"},
    "pdf_raw_material_price": {"de": "Marktpreis Edukte", "en": "Raw-material price"},
    "pdf_qualified_suppliers": {"de": "Qualifizierte Lieferanten", "en": "Qualified suppliers"},
    "pdf_lead_time": {"de": "Lieferzeit", "en": "Lead time"},
    "pdf_supplier_geography": {"de": "Lieferanten-Geografie", "en": "Supplier geography"},
    "pdf_geo_concentrated": {"de": "in einer Region konzentriert (Geopolitik-Risiko)", "en": "concentrated in one region (geopolitical risk)"},
    "pdf_geo_diversified": {"de": "geografisch diversifiziert", "en": "geographically diversified"},
    "pdf_concrete_recommendations": {"de": "Konkrete Empfehlungen", "en": "Concrete Recommendations"},
    "pdf_production_route": {"de": "Produktionsroute", "en": "Production route"},
    "pdf_expected_yield": {"de": "Erwartete Ausbeute", "en": "Expected yield"},
    "pdf_processing_time": {"de": "Prozessdauer", "en": "Process duration"},
    "pdf_achievable_purity": {"de": "Erreichbare Reinheit nach Aufarbeitung", "en": "Achievable purity after workup"},
    "pdf_downstream_steps": {"de": "Empfohlene Aufarbeitungsschritte", "en": "Recommended downstream steps"},
    "pdf_notes": {"de": "Hinweise", "en": "Notes"},
    "pdf_risks_tradeoffs": {"de": "Risiken & Trade-offs", "en": "Risks & Trade-offs"},
    "pdf_key_risks": {"de": "Wichtige Risiken", "en": "Key risks"},
    "pdf_tradeoffs": {"de": "Trade-offs", "en": "Trade-offs"},
    "pdf_recommended_actions": {"de": "Empfohlene Optimierungen", "en": "Recommended Optimizations"},
    "pdf_no_critical_issues": {
        "de": "• Keine kritischen Probleme identifiziert; Monitoring empfohlen.",
        "en": "• No critical issues identified; monitoring recommended.",
    },
    "pdf_no_improvements": {
        "de": "• Keine spezifischen Verbesserungen identifiziert. Erwägen Sie ein detailliertes Audit.",
        "en": "• No specific improvements identified. Consider a detailed audit.",
    },
    "pdf_primary_drivers": {
        "de": "Primäre Treiber: Syntheseschrittanzahl und Prozesskomplexität.",
        "en": "Primary drivers: number of synthesis steps and process complexity.",
    },
    "pdf_cost_structure": {"de": "Kostenstruktur", "en": "Cost Structure"},
    "pdf_cost_impact": {"de": "Kosten-Impact-Analyse", "en": "Cost Impact Analysis"},
    "pdf_main_cost_drivers": {"de": "Hauptkostentreiber:", "en": "Main cost drivers:"},
    "pdf_no_cost_drivers": {"de": "Hauptkostentreiber: Keine identifiziert.", "en": "Main cost drivers: none identified."},
    "pdf_savings_potential": {"de": "Geschätztes Einsparpotenzial", "en": "Estimated savings potential"},
    "pdf_downstream_insights": {"de": "Aufarbeitungs-Hinweise", "en": "Downstream insights"},
    "pdf_no_significant_risks": {"de": "• Keine signifikanten Risiken identifiziert.", "en": "• No significant risks identified."},
    "pdf_main_cost_default": {
        "de": "Hauptkostentreiber: Aufreinigung, Rohstoffe und operative Overheads.",
        "en": "Main cost drivers: purification, raw materials, and operational overhead.",
    },
    "pdf_cost_level": {"de": "Kosten-Level", "en": "Cost level"},
    "pdf_tradeoff_default_1": {
        "de": "Kosten vs. Reinheit: Höhere Reinheit erhöht typischerweise Aufreinigungskosten.",
        "en": "Cost vs. purity: higher purity typically increases purification costs.",
    },
    "pdf_tradeoff_default_2": {
        "de": "Effizienz vs. Komplexität: Mehr Prozessschritte erhöhen Risiko und reduzieren Nettoeffizienz.",
        "en": "Efficiency vs. complexity: more process steps increase risk and reduce net efficiency.",
    },
    # Bug M3 fix: replaced four generic fallback bullets with one honest
    # "no specific optimisations yet" line. Generic bullets undercut the
    # rule-driven optimisations elsewhere in the report.
    "pdf_no_specific_optimizations": {
        "de": "Für dieses Molekül liegen aktuell keine modul-spezifischen Optimierungsempfehlungen vor. Ein manuelles CMC-Audit auf Basis der Prozessparameter wird empfohlen.",
        "en": "No module-specific optimisation recommendations are available for this molecule yet. A manual CMC audit based on the process parameters is recommended.",
    },
    "pdf_fallback_action_1": {
        "de": "Bewerten Sie alternative Syntheserouten zur Reduktion der Schrittanzahl und Komplexität.",
        "en": "Evaluate alternative synthesis routes to reduce step count and complexity.",
    },
    "pdf_fallback_action_2": {
        "de": "Testen Sie verbesserte Aufreinigungsstrategien (z. B. Kristallisation) im Pilotmaßstab.",
        "en": "Test improved purification strategies (e.g. crystallization) at pilot scale.",
    },
    "pdf_fallback_action_3": {
        "de": "Analysieren und optimieren Sie Rohstofflieferanten, um Kosten zu senken.",
        "en": "Analyze and optimize raw-material suppliers to reduce costs.",
    },
    "pdf_fallback_action_4": {
        "de": "Validieren Sie Prozessstabilität und Leistung in Pilotversuchen.",
        "en": "Validate process stability and performance in pilot trials.",
    },
    "pdf_generated": {"de": "Erstellt", "en": "Generated"},
    "pdf_method_chemical": {"de": "Chemische Herstellung", "en": "Chemical production"},
    "pdf_method_biotech": {"de": "Biotechnologische Herstellung", "en": "Biotechnological production"},
    "pdf_method_extraction": {"de": "Natürliche Extraktion", "en": "Natural extraction"},
    "pdf_process_synthesis": {"de": "Synthese", "en": "Synthesis"},
    "pdf_process_fermentation": {"de": "Fermentation", "en": "Fermentation"},
    "pdf_process_extraction": {"de": "Extraktion", "en": "Extraction"},
    "pdf_rationale": {"de": "Begründung", "en": "Rationale"},
    "pdf_impact": {"de": "Erwartete Wirkung", "en": "Expected impact"},
    "pdf_prerequisites": {"de": "Voraussetzungen", "en": "Prerequisites"},
    "pdf_effort": {"de": "Aufwand", "en": "Effort"},
    "pdf_representation": {"de": "Repräsentation", "en": "Representation"},
    "pdf_target_molecule": {"de": "Zielmolekül", "en": "Target molecule"},
    "pdf_method": {"de": "Methode", "en": "Method"},
    "pdf_process_type": {"de": "Prozesstyp", "en": "Process type"},
    "pdf_date": {"de": "Datum", "en": "Date"},
    "pdf_brand_subtitle": {"de": "KI-gestützte Entscheidungsunterstützung", "en": "AI-assisted decision support"},
    "pdf_footer_disclaimer": {"de": "Nicht-operativer Entscheidungsreport – keine Laboranleitung", "en": "Non-operational decision report — not a lab protocol"},
    "pdf_footer_generated": {"de": "Erstellt", "en": "Generated"},

    # ====== Navigation / Landing / Feedback ======
    "nav_label": {"de": "Navigation", "en": "Navigation"},
    "nav_home": {"de": "Startseite", "en": "Home"},
    "nav_generate": {"de": "Empfehlung generieren", "en": "Generate recommendation"},
    "nav_demos": {"de": "Demos / Beispiele", "en": "Demos / Examples"},
    "nav_feedback": {"de": "Feedback geben", "en": "Give feedback"},
    "nav_settings": {"de": "Einstellungen", "en": "Settings"},

    # ====== Settings page ======
    "settings_page_title": {"de": "Einstellungen", "en": "Settings"},
    "settings_page_intro": {
        "de": "Konfigurieren Sie Sprache und weitere Optionen.",
        "en": "Configure language and other options.",
    },
    "settings_section_language": {"de": "Sprache", "en": "Language"},
    "settings_current_language": {"de": "Aktuelle Sprache", "en": "Current language"},
    "settings_current_language_help": {
        "de": "Sprache der aktuellen Sitzung — wird sofort angewandt.",
        "en": "Language for the current session — applied immediately.",
    },
    "settings_default_language": {"de": "Standardsprache", "en": "Default language"},
    "settings_default_language_help": {
        "de": "Sprache, die beim nächsten Start automatisch gewählt wird.",
        "en": "Language automatically selected on next start.",
    },
    "settings_saved": {"de": "Einstellung gespeichert.", "en": "Setting saved."},

    # ====== Homescreen hero + cards ======
    "hero_title": {
        "de": "KI-gestützte Prozess-Empfehlungen",
        "en": "AI-Powered Process Recommendations",
    },
    "hero_subtitle": {
        "de": "Fortschrittliche KI-Modelle für die Prozessentwicklung in Biotechnologie und chemischer Fertigung.",
        "en": "Advanced AI models for biotech and chemical manufacturing process development.",
    },
    "card_smart_analysis_title": {"de": "Intelligente Analyse", "en": "Smart Analysis"},
    "card_smart_analysis_body": {
        "de": "KI analysiert Ihre Prozess-Eingaben und Rahmenbedingungen.",
        "en": "AI analyzes your process inputs and constraints.",
    },
    "card_optimized_results_title": {"de": "Optimierte Ergebnisse", "en": "Optimized Results"},
    "card_optimized_results_body": {
        "de": "Datengetriebene Empfehlungen basierend auf etablierten Quellen.",
        "en": "Data-driven recommendations based on established sources.",
    },
    "card_export_share_title": {"de": "Export & Teilen", "en": "Export & Share"},
    "card_export_share_body": {
        "de": "Berichte als PDF herunterladen und mit Stakeholdern teilen.",
        "en": "Download reports as PDF and share with stakeholders.",
    },
    "start_card_title": {"de": "Neue Analyse starten", "en": "Start a new analysis"},
    "start_card_body": {
        "de": "Geben Sie Ihre Prozess-Daten ein und lassen Sie unsere KI die beste Empfehlung generieren.",
        "en": "Provide your process inputs and let our AI generate the best recommendation.",
    },
    "start_card_button": {"de": "Empfehlung generieren", "en": "Generate Recommendation"},

    # ====== Step indicator ======
    "step_inputs": {"de": "Eingaben", "en": "Inputs"},
    "step_analysis": {"de": "Analyse", "en": "Analysis"},
    "step_results": {"de": "Ergebnisse", "en": "Results"},

    # ====== Overview / Step 2 ======
    "overview_title": {"de": "Eingaben prüfen", "en": "Review inputs"},
    "overview_intro": {
        "de": "Bitte überprüfen Sie Ihre Eingaben. Klicken Sie auf 'Bericht anzeigen', um den Analysebericht zu generieren.",
        "en": "Please review your inputs. Click 'Show report' to generate the analysis report.",
    },
    "overview_back": {"de": "Zurück zu den Eingaben", "en": "Back to inputs"},
    "overview_continue": {"de": "Bericht anzeigen", "en": "Show report"},
    "overview_section_molecule": {"de": "Molekül", "en": "Molecule"},
    "overview_section_process": {"de": "Aktueller Prozess", "en": "Current process"},
    "overview_section_requirements": {"de": "Anforderungen", "en": "Requirements"},
    "overview_section_market": {"de": "Marktbedingungen", "en": "Market conditions"},
    "overview_section_methods": {"de": "Verfügbare Methodiken", "en": "Available methods"},
    "overview_section_biophys": {"de": "Biophysikalisches Profil", "en": "Biophysical profile"},

    # ====== Sensitivity / What-If analysis ======
    "sensitivity_title": {"de": "Sensitivitätsanalyse — Was wäre wenn?", "en": "Sensitivity Analysis — What If?"},
    "sensitivity_intro": {
        "de": ("Zeigt, wie sich Optimierungsempfehlungen und Warnungen ändern würden, "
               "wenn Sie einzelne Eingaben verändern. Hilft, den größten Hebel für "
               "Kostenreduktion / Prozessvereinfachung zu identifizieren."),
        "en": ("Shows how optimization recommendations and warnings would change if "
               "you altered individual inputs. Helps identify the largest lever for "
               "cost reduction / process simplification."),
    },
    "sensitivity_baseline_label": {
        "de": "Aktuell: {n_actions} Optimierung(en), {n_warnings} Warnung(en)",
        "en": "Current: {n_actions} optimization(s), {n_warnings} warning(s)",
    },
    "sensitivity_lever_current": {"de": "aktuell", "en": "current"},
    "sensitivity_insight_label": {"de": "💡 Größter Hebel", "en": "💡 Largest lever"},
    "pdf_sensitivity_title": {"de": "Sensitivitätsanalyse", "en": "Sensitivity Analysis"},
    "pdf_sensitivity_intro": {
        "de": "Auswirkung einzelner Eingaben auf die Anzahl der Empfehlungen und Warnungen.",
        "en": "Impact of individual inputs on the number of recommendations and warnings.",
    },

    # ====== Demos page ======
    "demo_page_title": {"de": "Demos / Beispielanalysen", "en": "Demos / Sample Analyses"},
    "demo_page_intro": {
        "de": ("Vorgefertigte Use-Cases zur schnellen Demonstration. Wählen Sie ein Beispiel — "
               "die Eingabewerte werden automatisch in das Formular geladen, dann genügt ein Klick "
               "auf „Prozess analysieren\"."),
        "en": ("Pre-built use cases for quick demonstration. Pick an example — the input values are "
               "loaded automatically into the form, then a single click on \"Analyze process\" runs "
               "the analysis."),
    },
    "demo_input_values": {"de": "Eingabewerte", "en": "Input values"},
    "demo_expected_outcome": {"de": "Erwartetes Ergebnis", "en": "Expected outcome"},
    "demo_load_button": {
        "de": "Diese Demo ins Formular laden",
        "en": "Load this demo into the form",
    },
    "demo_loaded_message": {
        "de": "Demo-Werte geladen — klicken Sie unten auf „Prozess analysieren\".",
        "en": "Demo values loaded — click \"Analyze process\" below.",
    },

    # Vanillin demo
    "demo_vanillin_title": {
        "de": "Vanillin via Biotechnologie (Pilot-Maßstab)",
        "en": "Vanillin via Biotechnology (pilot scale)",
    },
    "demo_vanillin_story": {
        "de": ("Ein mittelständischer Aroma-Hersteller plant, Vanillin nicht mehr aus Erdöl-Vorstufen, "
               "sondern fermentativ über engineered Hefe zu produzieren. Ziel: 100 kg/Jahr in food-grade "
               "Reinheit (98,5 %), bestehender Bioreaktor, robuste Sourcing-Lage."),
        "en": ("A mid-sized flavor manufacturer plans to switch from petroleum-based vanillin to "
               "fermentative production via engineered yeast. Target: 100 kg/year at food-grade "
               "purity (98.5 %), existing bioreactor, robust sourcing position."),
    },
    "demo_vanillin_outcome": {
        "de": ("Helixar bestätigt die Route (S. cerevisiae / A. niger auf Glucose oder Ferulasäure), "
               "schlägt 0,3–1,5 g/L Yield-Bereich vor und nennt konkrete Aufarbeitungsschritte. "
               "Wenige bis keine Plausibilitätswarnungen — Prozess ist sauber aufgesetzt."),
        "en": ("Helixar confirms the route (S. cerevisiae / A. niger on glucose or ferulic acid), "
               "suggests a 0.3–1.5 g/L yield range, and lists concrete downstream steps. "
               "Few or no plausibility warnings — the process is well configured."),
    },

    # Trastuzumab demo
    "demo_trastuzumab_title": {
        "de": "Trastuzumab — Antikörper-Produktion (industriell)",
        "en": "Trastuzumab — Antibody Production (industrial)",
    },
    "demo_trastuzumab_story": {
        "de": ("Ein Biotech-Unternehmen evaluiert die rekombinante Produktion eines HER2-Antikörpers "
               "(Trastuzumab-Biosimilar) im industriellen Maßstab. 200 kg/Jahr, 99,5 % Reinheit, "
               "vorhandene Protein-A-Affinitätschromatographie und TFF-Anlage."),
        "en": ("A biotech company evaluates recombinant production of a HER2 antibody (trastuzumab "
               "biosimilar) at industrial scale. 200 kg/year, 99.5 % purity, existing protein-A "
               "affinity chromatography and TFF setup."),
    },
    "demo_trastuzumab_outcome": {
        "de": ("Helixar erkennt korrekt, dass SMILES bei einem ~145 kDa Antikörper nicht anwendbar ist, "
               "und referenziert DrugBank:DB00072 mit dem Target HER2/ERBB2 (UniProt:P04626). "
               "Die Empfehlungen kommen aus typ-spezifischer Antikörper-Logik (Affinitäts-Capture, TFF)."),
        "en": ("Helixar correctly recognizes that SMILES is not applicable for a ~145 kDa antibody "
               "and references DrugBank:DB00072 with the target HER2/ERBB2 (UniProt:P04626). "
               "Recommendations are drawn from antibody-specific logic (affinity capture, TFF)."),
    },

    # Aspirin demo
    "demo_aspirin_title": {
        "de": "Aspirin — klassische chemische Synthese",
        "en": "Aspirin — Classical Chemical Synthesis",
    },
    "demo_aspirin_story": {
        "de": ("Ein Pharma-Generikahersteller produziert Acetylsalicylsäure (ASS) im großen Stil über "
               "die klassische 1-Schritt-Synthese aus Salicylsäure und Essigsäureanhydrid. 5.000 kg/Jahr, "
               "99 % API-grade, robuste Beschaffung, GMP-Abfallauflagen."),
        "en": ("A pharma generics manufacturer produces acetylsalicylic acid (ASA) at scale via the "
               "classical one-step synthesis from salicylic acid and acetic anhydride. 5,000 kg/year, "
               "99 % API-grade, robust sourcing, GMP waste constraints."),
    },
    "demo_aspirin_outcome": {
        "de": ("Helixar bestätigt den 1-Schritt-Pfad (Salicylsäure + Essigsäureanhydrid, kat. H2SO4), "
               "schlägt 75–92 % Yield vor und löst die Recommended Action „Lösungsmittelinventar auf "
               "grüne Alternativen prüfen\" aus (wegen GMP-Abfallauflagen + chemischer Methode)."),
        "en": ("Helixar confirms the one-step pathway (salicylic acid + acetic anhydride, cat. H2SO4), "
               "suggests 75–92 % yield, and triggers the recommended action \"Audit the solvent inventory "
               "for greener alternatives\" (due to GMP waste constraints + chemical method)."),
    },
    "back_to_home": {"de": "Zurück zur Startseite", "en": "Back to home"},
    "landing_slogan": {
        "de": "Von Molekül zu Produktionsstrategie",
        "en": "From molecule to production strategy",
    },
    "landing_call_to_action": {
        "de": "Wählen Sie im Menü 'Empfehlung generieren', um zu starten.",
        "en": "Select 'Generate recommendation' in the menu to get started.",
    },
    "feedback_page_title": {"de": "Feedback geben", "en": "Give feedback"},
    "feedback_input_label": {"de": "Ihr Feedback", "en": "Your feedback"},
    "feedback_input_placeholder": {
        "de": "Ihre Nachricht an das Team",
        "en": "Your message to the team",
    },
    "feedback_submit_button": {"de": "Feedback absenden", "en": "Submit feedback"},
    "feedback_warn_empty": {
        "de": "Bitte geben Sie zuerst Feedback ein.",
        "en": "Please enter your feedback first.",
    },
    "feedback_thanks": {
        "de": "Danke — Ihr Feedback wurde empfangen.",
        "en": "Thanks — your feedback was received.",
    },
    "feedback_page_intro": {
        "de": ("Hilf uns, Helixar AI besser zu machen. Dein Feedback geht direkt "
               "an das Helixar-Team. Pflichtfelder sind mit * markiert; alles "
               "andere ist optional."),
        "en": ("Help us make Helixar AI better. Your feedback goes directly to "
               "the Helixar team. Required fields are marked with *; everything "
               "else is optional."),
    },
    "feedback_label_name": {"de": "Name", "en": "Name"},
    "feedback_label_company": {"de": "Firma / Organisation", "en": "Company / organisation"},
    "feedback_label_email": {
        "de": "Deine E-Mail (für Rückfragen)",
        "en": "Your email (for follow-up)",
    },
    "feedback_label_category": {"de": "Kategorie *", "en": "Category *"},
    "feedback_label_text": {"de": "Dein Feedback *", "en": "Your feedback *"},
    "feedback_placeholder_name": {"de": "z. B. Anna Müller", "en": "e.g. Anna Miller"},
    "feedback_placeholder_company": {"de": "z. B. BioPharma GmbH", "en": "e.g. BioPharma Inc."},
    "feedback_placeholder_email": {"de": "name@firma.com", "en": "name@company.com"},
    "feedback_placeholder_text": {
        "de": "Beschreibe so konkret wie möglich — z. B. welches Molekül, welcher Schritt, was war erwartet vs. tatsächlich.",
        "en": "Be as specific as possible — e.g. which molecule, which step, expected vs. actual behaviour.",
    },
    "feedback_category_bug": {"de": "🐛 Bug / Fehler", "en": "🐛 Bug / error"},
    "feedback_category_feature": {"de": "✨ Feature-Wunsch", "en": "✨ Feature request"},
    "feedback_category_content": {"de": "🧬 Inhaltlich / Fachlich", "en": "🧬 Domain / scientific content"},
    "feedback_category_ux": {"de": "🎨 UX / Bedienung", "en": "🎨 UX / usability"},
    "feedback_category_other": {"de": "💬 Sonstiges", "en": "💬 Other"},
    "feedback_submit_success": {
        "de": "Danke! Dein Feedback ist beim Helixar-Team angekommen.",
        "en": "Thanks! Your feedback was delivered to the Helixar team.",
    },
    "feedback_submit_error_generic": {
        "de": "Das Senden ist fehlgeschlagen. Bitte versuche es in ein paar Minuten nochmal — oder schreib direkt an lorenzmeising@icloud.com.",
        "en": "Submission failed. Please try again in a few minutes — or write directly to lorenzmeising@icloud.com.",
    },
    "feedback_submit_error_network": {
        "de": "Keine Verbindung zum Mail-Service. Bitte prüfe deine Internetverbindung.",
        "en": "No connection to the mail service. Please check your internet connection.",
    },
    "feedback_warn_empty_text": {
        "de": "Bitte gib einen Feedback-Text ein, bevor du absendest.",
        "en": "Please enter feedback text before submitting.",
    },
    "feedback_email_invalid": {
        "de": "Die E-Mail-Adresse sieht nicht gültig aus. Du kannst das Feld auch leer lassen.",
        "en": "The email address does not look valid. You can also leave this field empty.",
    },
    "feedback_optional_hint": {
        "de": "Name, Firma und Mail sind optional — anonymes Feedback ist auch willkommen.",
        "en": "Name, company and email are optional — anonymous feedback is welcome too.",
    },
    "feedback_sending": {"de": "Wird gesendet …", "en": "Sending …"},
    "feedback_needs_activation": {
        "de": ("Einmaliger Schritt: Der Mail-Service wartet auf eine Aktivierung. "
               "Bitte prüfe das Postfach lorenzmeising@icloud.com, klicke den "
               "'Activate Form'-Link in der Mail von Formsubmit, und schicke "
               "danach diese Nachricht nochmal ab — dann kommt sie an."),
        "en": ("One-time step: the mail service is waiting for activation. "
               "Please check the lorenzmeising@icloud.com inbox, click the "
               "'Activate Form' link in the email from Formsubmit, and then "
               "resubmit this message — it will arrive then."),
    },
    "candidates_disabled": {
        "de": "Die Kandidatenauswahl wurde deaktiviert. Bitte nutzen Sie die Startseite, um eine neue Anfrage zu stellen.",
        "en": "Candidate selection has been disabled. Please use the home page to start a new request.",
    },

    # ====== Footer-Disclaimer ======
    "footer_disclaimer": {
        "de": "Hinweis: Dieses Tool bietet nicht-operatives Entscheidungs-Support. Es ersetzt keine Laborvalidierung oder regulatorische Beratung.",
        "en": "Note: This tool provides non-operational decision support. It does not replace lab validation or regulatory advice.",
    },

    # ====== Erweiterte Hinweise / Important Notices ======
    "notices_title": {"de": "Wichtige Hinweise", "en": "Important Notices"},
    "notice_compliance_title": {"de": "Nicht-GMP / Nicht-regulatorisch", "en": "Non-GMP / Non-regulatory"},
    "notice_compliance_body": {
        "de": ("Helixar AI ist ein Entscheidungsunterstützungs-Tool und kein GMP-, FDA- oder EMA-konformes "
               "System. Die Empfehlungen ersetzen weder regulatorische Beratung noch Laborvalidierung. "
               "Vor jeder Umsetzung in Produktion sind eigene Tests, Audits und Compliance-Prüfungen erforderlich."),
        "en": ("Helixar AI is a decision-support tool and not a GMP-, FDA-, or EMA-compliant system. "
               "Its recommendations do not replace regulatory advice or lab validation. "
               "Independent testing, audits, and compliance checks are required before any production use."),
    },
    "notice_data_title": {"de": "Datenschutz & lokale Verarbeitung", "en": "Data Privacy & Local Processing"},
    "notice_data_body": {
        "de": ("Eingegebene Daten (Molekülnamen, Prozessparameter, Lieferanteninformationen) werden lokal "
               "in Ihrer Sitzung verarbeitet und nicht an Dritte übertragen oder dauerhaft gespeichert. "
               "Bei Cloud-Deployment gelten die Datenschutzbestimmungen des Hosting-Providers."),
        "en": ("Inputs (molecule names, process parameters, supplier information) are processed locally "
               "within your session and are neither transmitted to third parties nor permanently stored. "
               "For cloud deployments, the hosting provider's privacy policy applies."),
    },
    "notice_ip_title": {"de": "Geistiges Eigentum (IP)", "en": "Intellectual Property (IP)"},
    "notice_ip_body": {
        "de": ("Eingegebene Molekül- und Prozessdaten verbleiben im Eigentum des Nutzers. Helixar AI "
               "beansprucht keine Rechte an Eingaben oder generierten Berichten. Empfehlungen basieren "
               "auf öffentlich zugänglichen Wissensquellen und allgemeinen Heuristiken."),
        "en": ("All input molecule and process data remain the property of the user. Helixar AI claims "
               "no rights to inputs or generated reports. Recommendations are based on publicly available "
               "knowledge sources and general heuristics."),
    },
    "notice_liability_title": {"de": "Haftungsausschluss", "en": "Liability Disclaimer"},
    "notice_liability_body": {
        "de": ("Empfehlungen sind regelbasiert und heuristisch — sie ersetzen kein Fachgutachten und "
               "keine experimentelle Validierung. Helixar AI übernimmt keine Haftung für Schäden, die "
               "aus der Anwendung der Empfehlungen ohne fachliche Prüfung entstehen."),
        "en": ("Recommendations are rule-based and heuristic — they do not replace expert review or "
               "experimental validation. Helixar AI assumes no liability for damages arising from "
               "applying the recommendations without professional review."),
    },
    "notice_sources_title": {"de": "Quellen / Literatur", "en": "Sources / Literature"},
    "notice_sources_body": {
        "de": ("Zitierte Quellen dienen der Nachvollziehbarkeit und ersetzen keine eigene "
               "systematische Literaturrecherche. Größenordnungen (Yields, Lieferzeiten, "
               "Reinheiten) sind als Plausibilitäts-Anker zu verstehen, nicht als rigorose "
               "Zitate für regulatorische Einreichungen."),
        "en": ("Cited sources are provided for traceability and do not replace systematic "
               "literature review. Reported ranges (yields, lead times, purities) are intended "
               "as plausibility anchors, not as rigorous citations for regulatory submissions."),
    },

    # ====== References section (Streamlit + PDF) ======
    "references_title": {"de": "Quellen / Literatur", "en": "References / Literature"},
    "references_intro": {
        "de": "Die folgenden Quellen lagen den oben gezeigten Empfehlungen, Plausibilitätshinweisen und Optimierungsvorschlägen zugrunde:",
        "en": "The following sources informed the recommendations, plausibility notes, and optimization actions shown above:",
    },
    "source_label_inline": {"de": "Quelle", "en": "Source"},

    # ====== Structure analysis (RDKit) ======
    "structure_analysis_title": {"de": "Strukturanalyse (RDKit)", "en": "Structure Analysis (RDKit)"},
    "structure_caption": {
        "de": "Eigenschaften berechnet aus der SMILES-Struktur — lokale Berechnung ohne externe Datenübertragung.",
        "en": "Properties computed from the SMILES structure — local computation, no external data transfer.",
    },
    "structure_label_mw": {"de": "Molekulargewicht", "en": "Molecular weight"},
    "structure_label_logp": {"de": "logP (Lipophilie)", "en": "logP (lipophilicity)"},
    "structure_label_tpsa": {"de": "TPSA (polare Oberfläche)", "en": "TPSA (polar surface area)"},
    "structure_label_arom_rings": {"de": "Aromatische Ringe", "en": "Aromatic rings"},
    "structure_label_hbd": {"de": "H-Brücken-Donoren", "en": "H-bond donors"},
    "structure_label_hba": {"de": "H-Brücken-Akzeptoren", "en": "H-bond acceptors"},
    "structure_label_rotbonds": {"de": "Rotierbare Bindungen", "en": "Rotatable bonds"},
    "structure_label_lipinski": {"de": "Lipinski-Verletzungen", "en": "Lipinski violations"},
    "structure_label_lipinski_pass": {"de": "Rule of 5: bestanden", "en": "Rule of 5: pass"},
    "structure_label_lipinski_fail": {"de": "Rule of 5: verletzt", "en": "Rule of 5: violated"},
    "structure_drawing_caption": {
        "de": "Strukturformel (gerendert via RDKit)",
        "en": "Structure (rendered via RDKit)",
    },
}


def get_lang() -> str:
    """Return the currently selected language from session_state, default to DE."""
    try:
        import streamlit as st  # local import to keep this module testable without streamlit
        lang = st.session_state.get("ui_lang") if hasattr(st, "session_state") else None
        if lang in SUPPORTED_LANGUAGES:
            return lang
    except Exception:
        pass
    return DEFAULT_LANGUAGE


def t(key: str, lang: Optional[str] = None) -> str:
    """Translate a key. Falls back to the key itself if no translation found."""
    if lang is None:
        lang = get_lang()
    if lang not in SUPPORTED_LANGUAGES:
        lang = DEFAULT_LANGUAGE
    entry = TRANSLATIONS.get(key)
    if not entry:
        return key
    return entry.get(lang) or entry.get(DEFAULT_LANGUAGE) or key
