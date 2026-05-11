"""Recommended optimization actions for an existing or proposed process.

Each action is a structured dict with:
  - title: short imperative phrase
  - rationale: why this would help (drawn from the input context)
  - expected_impact: qualitative + quantitative deltas where known
  - prerequisites: equipment / capabilities needed
  - effort: rough effort indicator (low | medium | high)
  - literature_source (optional): traceability reference for the rule

Outputs are localized in DE or EN. The `effort` field stays as a raw code
(low / medium / high) — translation to a localized label happens at the
display layer (UI / PDF) so the field can be filtered/sorted programmatically.
"""


# Literature anchors per rule — surfaced inline (per action) and aggregated
# into the references section of the PDF report.
RULE_SOURCES_BIOPHYS = {
    "low_tm_cold_chain": "Wang W., Int. J. Pharm. 1999 — Instability, stabilization, and formulation of liquid protein pharmaceuticals (doi:10.1016/S0378-5173(99)00152-0); Walsh G., Pharmaceutical Biotechnology, 2nd ed. 2018 (ISBN 978-1-119-11518-7)",
    "high_tm_lyophilizable": "Pikal MJ., Lyophilization Process Design, in Encyclopedia of Pharmaceutical Technology (doi:10.1081/E-EPT-120006023); Carpenter & Pikal, Pharm. Res. 1997 (doi:10.1023/A:1012180501551)",
    "low_tagg_formulation": "Wang W., Int. J. Pharm. 2005 — Protein aggregation and its inhibition (doi:10.1016/j.ijpharm.2004.11.014); Roberts CJ., Trends Biotechnol. 2014 (doi:10.1016/j.tibtech.2014.05.005)",
    "multi_domain_refolding": "Singh & Panda, J. Biosci. Bioeng. 2005 — Solubilization and refolding of bacterial inclusion body proteins (doi:10.1263/jbb.99.303); Walsh G., Pharmaceutical Biotechnology, 2nd ed. 2018 (ISBN 978-1-119-11518-7)",
    "ptm_eukaryotic_host": "Walsh G., Nat. Biotechnol. 2018 — Biopharmaceutical benchmarks (doi:10.1038/nbt.4313); Wurm FM., Nat. Biotechnol. 2004 — Production of recombinant protein therapeutics in CHO (doi:10.1038/nbt1026)",
}


RULE_SOURCES = {
    "switch_chem_to_biotech": "Anastas & Warner, Green Chemistry: Theory and Practice, Oxford Univ. Press 1998 (ISBN 978-0-19-850698-0); Sheldon, Green Chem. 2017 — E-factor analysis (doi:10.1039/C7GC00149E)",
    "step_consolidation": "Vollhardt & Schore, Organic Chemistry, 8th ed. 2018, Macmillan/Freeman (ISBN 978-1-319-07945-1) — telescoping and step economy; Newhouse et al., Angew. Chem. Int. Ed. 2009 (doi:10.1002/anie.200806239)",
    "spps_to_recombinant": "Bray BL., Nat. Rev. Drug Discov. 2003 — Large-scale manufacture of peptide therapeutics by chemical synthesis (doi:10.1038/nrd1133); Sanchez-Garcia et al., Microb. Cell Fact. 2016 — Recombinant pharmaceuticals from microbial cells (doi:10.1186/s12934-016-0437-3)",
    "biotech_titer_optimization": "Wurm FM., Nat. Biotechnol. 2004 — Production of recombinant protein therapeutics in CHO (doi:10.1038/nbt1026); Lim et al., Curr. Opin. Biotechnol. 2010 — High-yield expression systems for biopharmaceuticals (doi:10.1016/j.copbio.2010.07.005)",
    "purity_gap": "Walsh G., Pharmaceutical Biotechnology, 2nd ed. 2018, Wiley (ISBN 978-1-119-11518-7) — API-grade purification standards",
    "industrial_no_bioreactor": "Doran, Bioprocess Engineering Principles, 2nd ed. 2013, Academic Press (ISBN 978-0-12-220851-5) — bioreactor scale-up",
    "strict_waste_chemical": "Prat et al., Green Chem. 2016 — CHEM21 Solvent Selection Guide (doi:10.1039/C5GC01008J); Anastas & Warner, Green Chemistry: Theory and Practice, Oxford Univ. Press 1998 (ISBN 978-0-19-850698-0)",
    "single_source_dual": "ICH Q7 Good Manufacturing Practice Guide for Active Pharmaceutical Ingredients (Standard, frei verfügbar via ich.org) — supplier qualification requirements",
    "geo_diversification": "Yossi Sheffi, The New (Ab)Normal: Reshaping Business and Supply Chain Strategy Beyond Covid-19, MIT CTL Media 2020 (ISBN 978-1-73327-301-2)",
    "long_lead_time": "Christopher M., Logistics and Supply Chain Management, 5th ed. 2016, Pearson (ISBN 978-1-292-08379-7) — vendor-managed inventory and forecast collaboration",
    "raw_material_legacy": "ICH Q7 Good Manufacturing Practice Guide for APIs (Standard) — supplier qualification; Walsh G., Pharmaceutical Biotechnology, 2nd ed. 2018, Wiley (ISBN 978-1-119-11518-7)",
    "edukt_substitution": "Sheldon, Green Chem. 2017 — E-factor (doi:10.1039/C7GC00149E); Anastas & Warner, Green Chemistry: Theory and Practice, Oxford Univ. Press 1998 (ISBN 978-0-19-850698-0)",
}

from typing import Dict, Any, List


# Per-rule, per-language string templates. Each rule produces a single action;
# we keep the texts together for readability and for future translators.
ACTION_TEXTS = {
    "switch_chem_to_biotech": {
        "de": {
            "title": "Wechsel von chemischer zu biotechnologischer Synthese prüfen",
            "rationale": (
                "Für dieses Molekül existieren etablierte fermentative Routen "
                "(engineered Hefe oder E. coli), die meist mit weniger toxischen "
                "Lösungsmitteln und besserer Skalierbarkeit auskommen."
            ),
            "expected_impact": (
                "Reduktion halogenierter/organischer Lösungsmittel um typ. 60–90 %, "
                "Yield abhängig vom Stamm — bei optimierten Stämmen vergleichbar oder höher; "
                "Time-to-market 6–18 Monate für Stammoptimierung."
            ),
            "prerequisites": "Bioreaktor ≥ 1 L, Fermentations-Know-how, ggf. Stammlizenz.",
        },
        "en": {
            "title": "Evaluate switch from chemical to biotechnological synthesis",
            "rationale": (
                "For this molecule, established fermentative routes (engineered yeast "
                "or E. coli) exist that typically require fewer toxic solvents and "
                "scale better."
            ),
            "expected_impact": (
                "Reduction of halogenated/organic solvents by typically 60–90 %; "
                "yield depends on strain — comparable or higher for optimized strains; "
                "time-to-market 6–18 months for strain optimization."
            ),
            "prerequisites": "Bioreactor ≥ 1 L, fermentation know-how, possibly strain license.",
        },
    },
    "spps_to_recombinant": {
        "de": {
            "title": "Von SPPS zu rekombinanter Produktion wechseln ({steps} Kupplungen)",
            "rationale": (
                "SPPS-Schrittzahl entspricht der Sequenzlänge — sie lässt sich nicht durch "
                "Prozessoptimierung reduzieren. Bei {steps} Aminosäure-Kupplungen "
                "akkumulieren Yield-Verluste exponentiell. Rekombinante Expression umgeht "
                "diese Limitation und ist für längere Peptide (> 15 AS) kommerziell etabliert."
            ),
            "expected_impact": (
                "Gesamtausbeute steigt von typ. 30–60 % (SPPS) auf 80–95 % (rekombinant); "
                "COGS-Reduktion 50–80 %. Lösungsmittelverbrauch sinkt drastisch (kein DMF/DCM). "
                "Trade-off: Time-to-Market 12–24 Monate für Expressions-System-Etablierung."
            ),
            "prerequisites": "Fermentations-Infrastruktur (≥ 5 L Bioreaktor), Molekularbiologie-Know-how, ggf. Fusionsprotein-Strategie für schwierige Sequenzen.",
        },
        "en": {
            "title": "Switch from SPPS to recombinant production ({steps} couplings)",
            "rationale": (
                "SPPS step count equals sequence length — it cannot be reduced by process "
                "optimization. At {steps} amino-acid couplings, yield losses accumulate "
                "exponentially. Recombinant expression bypasses this and is commercially "
                "established for longer peptides (> 15 aa)."
            ),
            "expected_impact": (
                "Overall yield rises from typically 30–60 % (SPPS) to 80–95 % (recombinant); "
                "COGS reduction 50–80 %. Solvent use drops dramatically (no DMF/DCM). "
                "Trade-off: 12–24 months time-to-market for expression-system establishment."
            ),
            "prerequisites": "Fermentation infrastructure (≥ 5 L bioreactor), molecular-biology know-how, possibly fusion-protein strategy for difficult sequences.",
        },
    },
    "biotech_titer_optimization": {
        "de": {
            "title": "Fermentations-Titer und produktive Wirtschaftlichkeit optimieren",
            "rationale": (
                "Bei biotechnologischer Produktion von Peptiden/Proteinen dominieren Bioreaktor-Footprint, "
                "Medienkosten und Aufarbeitung die Herstellkosten. Höhere Titer (g/L) reduzieren alle drei "
                "proportional. Hebel: Stammoptimierung (CRISPR-Engineering, Promoter-Tuning), Fed-Batch- "
                "oder Perfusions-Strategien, Medien-DoE."
            ),
            "expected_impact": (
                "Titer-Verdopplung (typisch erreichbar in 6–12 Monaten) halbiert Bioreaktor-Stunden "
                "und Medienkosten — entspricht typischerweise 30–50 % Senkung der gesamten COGS. "
                "Skalierungs-Risiken sinken, weil weniger Equipment für gleiche Output benötigt wird."
            ),
            "prerequisites": "Fermentations-Labor (≥ 1 L), Analytik für Titer-Bestimmung (HPLC/ELISA), Stammbank/Engineering-Zugang.",
        },
        "en": {
            "title": "Optimize fermentation titer and process economics",
            "rationale": (
                "In biotech production of peptides/proteins, bioreactor footprint, media cost and "
                "downstream processing dominate COGS. Higher titers (g/L) reduce all three proportionally. "
                "Levers: strain engineering (CRISPR, promoter tuning), fed-batch or perfusion strategies, "
                "media DoE."
            ),
            "expected_impact": (
                "Titer doubling (typically achievable in 6–12 months) halves bioreactor hours and media "
                "cost — corresponds to ~30–50 % overall COGS reduction. Scale-up risks decrease because "
                "less equipment is needed for the same output."
            ),
            "prerequisites": "Fermentation lab (≥ 1 L), analytics for titer measurement (HPLC/ELISA), strain-bank / engineering access.",
        },
    },
    "step_consolidation": {
        "de": {
            "title": "Schrittzahl reduzieren ({steps} → ggf. {target_steps})",
            "rationale": (
                "Mehrstufige Synthesen verlieren typischerweise 10–25 % Yield pro Schritt. "
                "Prüfen Sie One-Pot-Reaktionen, telescoping (mehrere Stufen ohne Aufarbeitung) "
                "oder Schutzgruppen-Vermeidung."
            ),
            "expected_impact": (
                "Yield-Steigerung pro eingesparter Stufe ca. 10–20 %. "
                "Bei {steps} → {target_steps} Stufen: zusätzliche Gesamt-Yield ≈ {pct} %. "
                "Lösungsmittel- und Personalkosten sinken proportional."
            ),
            "prerequisites": "Reaktionsoptimierung im Labormaßstab, ggf. DoE.",
        },
        "en": {
            "title": "Reduce step count ({steps} → potentially {target_steps})",
            "rationale": (
                "Multi-step syntheses typically lose 10–25 % yield per step. "
                "Consider one-pot reactions, telescoping (multiple stages without workup), "
                "or protecting-group avoidance."
            ),
            "expected_impact": (
                "Per-step yield gain ~10–20 %. "
                "For {steps} → {target_steps} stages: additional overall yield ≈ {pct} %. "
                "Solvent and labor costs decrease proportionally."
            ),
            "prerequisites": "Reaction optimization at lab scale, possibly DoE.",
        },
    },
    "purity_gap": {
        "de": {
            "title": "Hochauflösende Reinigungsmethode beschaffen oder auslagern",
            "rationale": (
                "Die geforderte Reinheit {purity_str} ist mit den aktuell verfügbaren Methoden "
                "in der Regel nicht erreichbar. Prep-HPLC, FPLC oder Umkristallisation sind "
                "Standardlösungen."
            ),
            "expected_impact": (
                "Erreichen der Zielreinheit ohne diese Methoden ist unrealistisch — "
                "ohne Investition wird die spezifizierte Reinheit nicht zuverlässig erreicht."
            ),
            "prerequisites": "Investition in prep-HPLC (~50–150 k€) oder externer CDMO-Service.",
        },
        "en": {
            "title": "Procure or outsource a high-resolution purification method",
            "rationale": (
                "The required purity of {purity_str} is generally not achievable with the "
                "currently available methods. Prep-HPLC, FPLC, or recrystallization are "
                "standard solutions."
            ),
            "expected_impact": (
                "Reaching the target purity without these methods is unrealistic — "
                "without investment, the specified purity will not be reliably reached."
            ),
            "prerequisites": "Investment in prep-HPLC (~€50–150 k) or external CDMO service.",
        },
    },
    "industrial_no_bioreactor": {
        "de": {
            "title": "Bioreaktor-Kapazität sichern",
            "rationale": (
                "Industrielle biotechnologische Produktion{scale_qty} ohne Bioreaktor ist nicht möglich. "
                "Auswahl je nach Maßstab: Stainless-Steel-Fermenter (≥ 1 m³) oder "
                "Single-Use-Bioreaktor (bis ~2 m³) für flexible Produktion."
            ),
            "expected_impact": (
                "Ohne diese Investition kann die geplante Skala nicht erreicht werden. "
                "CDMO-Auslagerung als kurzfristige Alternative möglich."
            ),
            "prerequisites": "Investition oder CDMO-Vertrag.",
        },
        "en": {
            "title": "Secure bioreactor capacity",
            "rationale": (
                "Industrial biotechnological production{scale_qty} is not possible without a bioreactor. "
                "Choice by scale: stainless-steel fermenter (≥ 1 m³) or single-use "
                "bioreactor (up to ~2 m³) for flexible production."
            ),
            "expected_impact": (
                "Without this investment the planned scale cannot be reached. "
                "CDMO outsourcing is possible as a short-term alternative."
            ),
            "prerequisites": "Capital investment or CDMO contract.",
        },
    },
    "strict_waste_chemical": {
        "de": {
            "title": "Lösungsmittelinventar auf grüne Alternativen prüfen",
            "rationale": (
                "Bei strikten Abfallauflagen verursachen halogenierte und stark toxische "
                "Lösungsmittel (DCM, DMF, Pyridin) hohe Entsorgungskosten und regulatorische "
                "Hürden. CHEM21-Solvent-Selection-Guide ist ein etablierter Startpunkt."
            ),
            "expected_impact": (
                "Reduktion der Abfall-Disposal-Kosten typ. 30–60 %; "
                "verbesserte ESG- und GMP-Compliance."
            ),
            "prerequisites": "Lösungsmittel-Mapping und ggf. Reaktionsoptimierung.",
        },
        "en": {
            "title": "Audit the solvent inventory for greener alternatives",
            "rationale": (
                "Under strict waste constraints, halogenated and highly toxic solvents "
                "(DCM, DMF, pyridine) incur high disposal costs and regulatory hurdles. "
                "The CHEM21 solvent selection guide is an established starting point."
            ),
            "expected_impact": (
                "Reduction of waste-disposal costs typically 30–60 %; "
                "improved ESG and GMP compliance."
            ),
            "prerequisites": "Solvent mapping and possibly reaction optimization.",
        },
    },
    "single_source_dual": {
        "de": {
            "title": "Zweitlieferanten qualifizieren (Single-Source-Risiko reduzieren)",
            "rationale": (
                "Aktuell ist nur 1 qualifizierter Lieferant gelistet{lead_str}. "
                "Bei Ausfall (technisch, logistisch, regulatorisch) steht die Produktion "
                "still. Für GMP-Produktion ist ein qualifizierter Backup ohnehin Pflicht."
            ),
            "expected_impact": (
                "Eliminiert Single-Point-of-Failure; verbessert Verhandlungsposition "
                "bei Preisgesprächen typ. um 5–15 %; ermöglicht GMP-Pfad."
            ),
            "prerequisites": "Sourcing + Qualifizierung (Audits, Spezifikations-Tests, ggf. Re-Validierung).",
        },
        "en": {
            "title": "Qualify a second supplier (reduce single-source risk)",
            "rationale": (
                "Currently only 1 qualified supplier is listed{lead_str}. "
                "On supplier outage (technical, logistical, regulatory) production halts. "
                "For GMP production, a qualified backup is mandatory anyway."
            ),
            "expected_impact": (
                "Eliminates single point of failure; improves price negotiation position "
                "typically by 5–15 %; enables GMP path."
            ),
            "prerequisites": "Sourcing + qualification (audits, specification tests, possibly re-validation).",
        },
    },
    "geo_diversification": {
        "de": {
            "title": "Lieferanten geografisch diversifizieren",
            "rationale": (
                "Alle qualifizierten Lieferanten sitzen in einer Region. Geopolitische "
                "Ereignisse (Export-Restriktionen, Sanktionen, regionale Ausfälle) "
                "treffen die Versorgung mit voller Wucht."
            ),
            "expected_impact": (
                "Reduktion des Geopolitik-Risikos; resilientere Supply Chain. "
                "Insbesondere relevant seit Covid und aktueller China-/EU-Spannungen."
            ),
            "prerequisites": "Sourcing-Audit in alternativer Region + Qualifizierung.",
        },
        "en": {
            "title": "Diversify suppliers geographically",
            "rationale": (
                "All qualified suppliers sit in a single region. Geopolitical events "
                "(export restrictions, sanctions, regional outages) hit supply with full force."
            ),
            "expected_impact": (
                "Reduces geopolitical risk; more resilient supply chain. "
                "Especially relevant since Covid and current China/EU tensions."
            ),
            "prerequisites": "Sourcing audit in an alternative region + qualification.",
        },
    },
    "long_lead_time": {
        "de": {
            "title": "Lieferzeit reduzieren oder Sicherheitsbestand strategisch erhöhen",
            "rationale": (
                "Lieferzeit von {lead:.1f} Wochen erfordert große Sicherheitsbestände "
                "(Kapitalbindung) oder Vorab-Forecast-Vereinbarungen mit Lieferanten "
                "(Mindestabnahmen)."
            ),
            "expected_impact": (
                "Bei Vendor-Managed-Inventory (VMI) oder Konsignationslager: Reduktion "
                "der Kapitalbindung um 30–60 % bei gleicher Versorgungssicherheit."
            ),
            "prerequisites": "Verhandlung mit Lieferanten, ggf. Mindestabnahme-Commitment.",
        },
        "en": {
            "title": "Reduce lead time or strategically increase safety stock",
            "rationale": (
                "Lead time of {lead:.1f} weeks requires large safety stock (tied capital) "
                "or forward forecast agreements with suppliers (minimum take-or-pay)."
            ),
            "expected_impact": (
                "With vendor-managed inventory (VMI) or consignment stock: 30–60 % "
                "reduction in tied capital at the same supply security."
            ),
            "prerequisites": "Supplier negotiation, possibly minimum-purchase commitment.",
        },
    },
    "raw_material_legacy": {
        "de": {
            "title": "Rohstoff-Risiko durch Zweitlieferant oder Substitution mindern",
            "rationale": (
                "Niedrige Marktverfügbarkeit ist ein Single-Point-of-Failure für die "
                "Versorgungssicherheit. Ein qualifizierter Zweitlieferant oder ein "
                "alternativer Edukt-Pfad reduziert das Lieferrisiko."
            ),
            "expected_impact": "Vermeidung von Produktionsausfällen; Verhandlungsspielraum bei Preisen.",
            "prerequisites": "Sourcing-Aufwand + ggf. Re-Validierung mit zweitem Edukt.",
        },
        "en": {
            "title": "Mitigate raw-material risk via second supplier or substitution",
            "rationale": (
                "Low market availability is a single point of failure for supply security. "
                "A qualified second supplier or an alternative starting-material path reduces supply risk."
            ),
            "expected_impact": "Avoids production outages; provides negotiation room on prices.",
            "prerequisites": "Sourcing effort + possibly re-validation with second starting material.",
        },
    },
    # ----- Biophysical numeric-driven actions -----
    "low_tm_cold_chain": {
        "de": {
            "title": "Kühlkette aufbauen / Cryoprotectant-Formulierung evaluieren",
            "rationale": (
                "Mit Tm = {tm:.1f} °C ist das Protein bei Raumtemperatur thermisch instabil. "
                "Standard-Logistik bei 4 °C reicht meist nicht — Lagerung typ. < 0 °C oder "
                "Cryoprotectant (Trehalose, Sucrose, Saccharose) zur Stabilisierung notwendig."
            ),
            "expected_impact": (
                "Sicherstellung der Produkt-Integrität bis zum Patienten; "
                "Reduktion von Aktivitätsverlust und Aggregation während der Lagerung."
            ),
            "prerequisites": "Cold-Chain-Logistik (≤ –20 °C oder ≤ –70 °C); Formulation-Studien.",
        },
        "en": {
            "title": "Establish cold chain / evaluate cryoprotectant formulation",
            "rationale": (
                "With Tm = {tm:.1f} °C the protein is thermally unstable at room temperature. "
                "Standard 4 °C logistics is often insufficient — storage typically < 0 °C "
                "or cryoprotectant (trehalose, sucrose) for stabilization required."
            ),
            "expected_impact": (
                "Ensures product integrity to the patient; reduces activity loss and "
                "aggregation during storage."
            ),
            "prerequisites": "Cold-chain logistics (≤ –20 °C or ≤ –70 °C); formulation studies.",
        },
    },
    "high_tm_lyophilizable": {
        "de": {
            "title": "Lyophilisation als bevorzugte Endformulierung in Betracht ziehen",
            "rationale": (
                "Mit Tm = {tm:.1f} °C ist das Protein thermisch robust und für die Gefriertrocknung "
                "gut geeignet. Lyophilisat ermöglicht Lagerung bei Raumtemperatur und reduziert "
                "Logistikkosten gegenüber kühlpflichtigen Flüssigformulierungen."
            ),
            "expected_impact": (
                "Verlängerte Haltbarkeit (24–36 Monate bei RT statt 12–18 Monate bei 2–8 °C); "
                "deutlich geringere Lagerungs- und Versandkosten."
            ),
            "prerequisites": "Lyophilisator + Formulation-Studien zur Auswahl Cryoprotectant/Bulking-Agent.",
        },
        "en": {
            "title": "Consider lyophilization as the preferred final formulation",
            "rationale": (
                "With Tm = {tm:.1f} °C the protein is thermally robust and well-suited for "
                "freeze-drying. Lyophilisate enables room-temperature storage and reduces "
                "logistics costs vs. cold-chain liquid formulations."
            ),
            "expected_impact": (
                "Extended shelf life (24–36 months at RT vs. 12–18 months at 2–8 °C); "
                "significantly lower storage and shipping cost."
            ),
            "prerequisites": "Lyophilizer + formulation studies for cryoprotectant/bulking agent selection.",
        },
    },
    "low_tagg_formulation": {
        "de": {
            "title": "Aggregations-Risiko durch Formulation-Engineering mindern",
            "rationale": (
                "Tagg = {tagg:.1f} °C bedeutet Aggregation tritt deutlich vor klassischen "
                "Prozesstemperaturen (37 °C, 25 °C) ein. Risiko: Wirkverlust, Immunogenität. "
                "Surfactant (Polysorbat 80, Poloxamer 188), pH-Optimierung, Argininzusatz "
                "können den Onset deutlich anheben."
            ),
            "expected_impact": (
                "Tagg-Verschiebung um typ. +10–20 °C möglich; "
                "Reduktion von Aggregat-Spezies < 1 % über Lagerzeit."
            ),
            "prerequisites": "Formulation-Screening (DSF, SEC-HPLC); High-Throughput-Tools wie SoluPro.",
        },
        "en": {
            "title": "Mitigate aggregation risk via formulation engineering",
            "rationale": (
                "Tagg = {tagg:.1f} °C means aggregation onset is well below typical process "
                "temperatures (37 °C, 25 °C). Risk: activity loss, immunogenicity. "
                "Surfactants (polysorbate 80, poloxamer 188), pH optimization, arginine "
                "additives can substantially raise the onset."
            ),
            "expected_impact": (
                "Tagg shift of typically +10–20 °C achievable; aggregate species reduction "
                "below 1 % over shelf life."
            ),
            "prerequisites": "Formulation screening (DSF, SEC-HPLC); high-throughput tools such as SoluPro.",
        },
    },
    "multi_domain_refolding": {
        "de": {
            "title": "Refolding-Strategie für Multi-Domain-Protein etablieren",
            "rationale": (
                "Mit {domains} Domänen und {disulfides} Disulfidbrücken ist die korrekte Faltung "
                "nach E. coli Inclusion-Body-Expression nicht trivial. Falsch gefaltete Proteine "
                "verlieren Aktivität und können immunogene Epitope exponieren. "
                "Alternative: eukaryotische Expression (CHO/Hefe) — meist korrekt gefaltet, aber teurer."
            ),
            "expected_impact": (
                "Korrekt gefalteter Anteil bei optimiertem Refolding typ. 30–60 %; "
                "ohne Refolding-Optimierung oft < 10 %. CHO/Hefe: > 90 % korrekt gefaltet."
            ),
            "prerequisites": "Refolding-Screen (Puffer, pH, Redoxpaar GSH/GSSG, Detergent); oder eukaryot. Wirt.",
        },
        "en": {
            "title": "Establish refolding strategy for multi-domain protein",
            "rationale": (
                "With {domains} domains and {disulfides} disulfide bridges, correct folding "
                "after E. coli inclusion-body expression is non-trivial. Misfolded protein "
                "loses activity and may expose immunogenic epitopes. "
                "Alternative: eukaryotic expression (CHO/yeast) — usually folded correctly, but pricier."
            ),
            "expected_impact": (
                "Correctly folded fraction with optimized refolding typically 30–60 %; "
                "without optimization often < 10 %. CHO/yeast: > 90 % correctly folded."
            ),
            "prerequisites": "Refolding screen (buffer, pH, GSH/GSSG redox pair, detergent); or eukaryotic host.",
        },
    },
    "ptm_eukaryotic_host": {
        "de": {
            "title": "Eukaryotisches Expressionssystem für PTMs einsetzen",
            "rationale": (
                "Posttranslationale Modifikationen (Glykosylierung, korrekte Disulfidverknüpfung, "
                "Phosphorylierung) können von E. coli nicht oder nur unzureichend ausgeführt werden. "
                "Wahl: CHO-Zellen (komplexe Glykosylierung, Antikörper-Standard), HEK293 "
                "(humanähnlich), Pichia pastoris (Hefe-Glykosylierung, kostengünstiger)."
            ),
            "expected_impact": (
                "Aktive, korrekt modifizierte Proteine ab Stufe 1; "
                "vermeidet kostspielige Reformulierungen oder klinische Misserfolge."
            ),
            "prerequisites": "Zellbank-Aufbau in eukaryot. System; entsprechende Bioreaktor-Infrastruktur.",
        },
        "en": {
            "title": "Use eukaryotic expression system for PTMs",
            "rationale": (
                "Post-translational modifications (glycosylation, correct disulfide formation, "
                "phosphorylation) cannot be performed by E. coli or only insufficiently. "
                "Options: CHO cells (complex glycosylation, antibody standard), HEK293 "
                "(human-like), Pichia pastoris (yeast glycosylation, more cost-efficient)."
            ),
            "expected_impact": (
                "Active, correctly modified protein from step 1; "
                "avoids costly reformulations or clinical failures."
            ),
            "prerequisites": "Cell-bank development in eukaryotic system; matching bioreactor infrastructure.",
        },
    },
    "edukt_substitution": {
        "de": {
            "title": "Edukt-Substitution oder Recyclingsstrategie evaluieren",
            "rationale": (
                "Bei Edukt-Kosten > 500 €/kg{cost_qty} dominieren die Rohstoffe die Herstellkosten. "
                "Optionen: günstigere Vorstufen synthetisieren, Lösungsmittel- und "
                "Reagenzien-Recycling, Optimierung der Stoichiometrie."
            ),
            "expected_impact": (
                "Senkung der Herstellkosten um 20–50 % bei substantieller Edukt-Substitution."
            ),
            "prerequisites": "Methoden-Re-Entwicklung, ggf. Stoffbilanz-Analyse.",
        },
        "en": {
            "title": "Evaluate starting-material substitution or recycling strategy",
            "rationale": (
                "At starting-material costs > €500/kg{cost_qty}, raw materials dominate "
                "manufacturing cost. Options: synthesize cheaper precursors, solvent and "
                "reagent recycling, stoichiometry optimization."
            ),
            "expected_impact": "Manufacturing cost reduction of 20–50 % with substantial substitution.",
            "prerequisites": "Method redevelopment, possibly mass-balance analysis.",
        },
    },
}


def _make(rule_key: str, lang: str, effort: str,
          current_state: str = "", optimized_state: str = "",
          **fmt) -> Dict[str, Any]:
    """Format a localized action dict from ACTION_TEXTS.

    Optional `current_state` and `optimized_state` are pre-formatted strings
    (caller-side, lang-specific) describing the as-is situation and the
    post-optimization state — they replace what the sensitivity panel did before.
    """
    bundle = ACTION_TEXTS.get(rule_key, {}).get(lang) or ACTION_TEXTS.get(rule_key, {}).get("de") or {}
    out: Dict[str, Any] = {"effort": effort}
    for field in ("title", "rationale", "expected_impact", "prerequisites"):
        template = bundle.get(field, "")
        try:
            out[field] = template.format(**fmt) if fmt else template
        except Exception:
            out[field] = template
    if current_state:
        out["current_state"] = current_state
    if optimized_state:
        out["optimized_state"] = optimized_state
    # Sources may live in either RULE_SOURCES or RULE_SOURCES_BIOPHYS
    src = RULE_SOURCES.get(rule_key) or RULE_SOURCES_BIOPHYS.get(rule_key)
    if src:
        out["literature_source"] = src
    return out


# --- Biophysical numeric-driven rules ---

def _rule_biophys_low_tm(p: Dict[str, Any], lang: str) -> List[Dict[str, Any]]:
    """Cold-chain or cryoprotectant action when Tm is low."""
    tm = p.get("tm_celsius")
    mtype = (p.get("molecule_type") or "").lower()
    if mtype not in ("peptide", "protein"):
        return []
    if not isinstance(tm, (int, float)):
        return []
    if float(tm) >= 50.0:
        return []
    if lang == "en":
        cur = f"Tm = {tm:.1f} °C — protein thermally unstable at room temperature"
        opt = "With ≤ –20 °C deep-freeze logistics or cryoprotectant formulation (trehalose, sucrose): full activity preserved through shelf life"
    else:
        cur = f"Tm = {tm:.1f} °C — Protein bei Raumtemperatur thermisch instabil"
        opt = "Mit ≤ –20 °C Tiefkühl-Logistik oder Cryoprotectant-Formulierung (Trehalose, Saccharose): volle Aktivität über die Lagerzeit erhalten"
    return [_make("low_tm_cold_chain", lang, effort="medium",
                  current_state=cur, optimized_state=opt, tm=float(tm))]


def _rule_biophys_high_tm(p: Dict[str, Any], lang: str) -> List[Dict[str, Any]]:
    """Lyophilization opportunity when Tm is high."""
    tm = p.get("tm_celsius")
    mtype = (p.get("molecule_type") or "").lower()
    if mtype not in ("peptide", "protein"):
        return []
    if not isinstance(tm, (int, float)):
        return []
    if float(tm) < 70.0:
        return []
    methods = p.get("available_methods") or {}
    if methods.get("has_lyophilizer"):
        return []
    if lang == "en":
        cur = f"Tm = {tm:.1f} °C — thermally robust, but only stored as cold-chain liquid (no lyophilizer in-house)"
        opt = "With lyophilization: 24–36 months room-temperature shelf life (vs. 12–18 months at 2–8 °C), substantial logistics-cost reduction"
    else:
        cur = f"Tm = {tm:.1f} °C — thermisch robust, aber nur als kühlpflichtige Flüssigformulierung gelagert (kein Lyophilisator vorhanden)"
        opt = "Mit Lyophilisation: 24–36 Monate Raumtemperatur-Haltbarkeit (statt 12–18 Monate bei 2–8 °C), deutliche Logistikkosten-Senkung"
    return [_make("high_tm_lyophilizable", lang, effort="medium",
                  current_state=cur, optimized_state=opt, tm=float(tm))]


def _rule_biophys_low_tagg(p: Dict[str, Any], lang: str) -> List[Dict[str, Any]]:
    """Formulation engineering when Tagg is low."""
    tagg = p.get("tagg_celsius")
    mtype = (p.get("molecule_type") or "").lower()
    if mtype not in ("peptide", "protein"):
        return []
    if not isinstance(tagg, (int, float)):
        return []
    if float(tagg) >= 55.0:
        return []
    if lang == "en":
        cur = f"Tagg = {tagg:.1f} °C — aggregation onset below typical process temperatures (risk of activity loss, immunogenicity)"
        opt = "With formulation engineering (polysorbate 80, poloxamer 188, arginine, pH optimization): Tagg shift of typically +10–20 °C, aggregates < 1 %"
    else:
        cur = f"Tagg = {tagg:.1f} °C — Aggregations-Onset unter typischen Prozesstemperaturen (Risiko: Aktivitätsverlust, Immunogenität)"
        opt = "Mit Formulation-Engineering (Polysorbat 80, Poloxamer 188, Arginin, pH-Optimierung): Tagg-Verschiebung um typ. +10–20 °C, Aggregate < 1 %"
    return [_make("low_tagg_formulation", lang, effort="high",
                  current_state=cur, optimized_state=opt, tagg=float(tagg))]


def _rule_biophys_multi_domain(p: Dict[str, Any], lang: str) -> List[Dict[str, Any]]:
    """Refolding strategy when many domains and/or disulfides."""
    domains = p.get("num_domains")
    disulfides = p.get("num_disulfides")
    mtype = (p.get("molecule_type") or "").lower()
    if mtype != "protein":
        return []
    if not (isinstance(domains, int) and isinstance(disulfides, int)):
        return []
    if domains < 3 and disulfides < 4:
        return []
    if lang == "en":
        cur = f"{domains} domains, {disulfides} disulfides — high folding complexity, E. coli inclusion bodies typically < 10 % correctly folded"
        opt = "With optimized refolding screen (GSH/GSSG redox, detergent, pH): 30–60 % correctly folded; alternative eukaryotic host (CHO/yeast): > 90 %"
    else:
        cur = f"{domains} Domänen, {disulfides} Disulfide — hohe Faltungs-Komplexität, E. coli Inclusion Bodies typisch < 10 % korrekt gefaltet"
        opt = "Mit optimiertem Refolding-Screen (GSH/GSSG-Redox, Detergent, pH): 30–60 % korrekt gefaltet; alternativ eukaryot. Wirt (CHO/Hefe): > 90 %"
    return [_make("multi_domain_refolding", lang, effort="high",
                  current_state=cur, optimized_state=opt,
                  domains=domains, disulfides=disulfides)]


def _rule_biophys_ptm(p: Dict[str, Any], lang: str) -> List[Dict[str, Any]]:
    """Eukaryotic host required when PTMs are needed."""
    ptm = p.get("has_ptm")
    mtype = (p.get("molecule_type") or "").lower()
    if mtype != "protein" or not ptm:
        return []
    if lang == "en":
        cur = "Post-translational modifications required, but no eukaryotic expression system specified (E. coli cannot perform PTMs)"
        opt = "With CHO cells (complex glycosylation, antibody standard), HEK293 (human-like) or Pichia pastoris (cheaper yeast glycosylation): correctly modified protein from step 1"
    else:
        cur = "Posttranslationale Modifikationen erforderlich, aber kein eukaryotisches Expressionssystem vorgesehen (E. coli kann keine PTMs)"
        opt = "Mit CHO-Zellen (komplexe Glykosylierung, Antikörper-Standard), HEK293 (humanähnlich) oder Pichia pastoris (günstigere Hefe-Glykosylierung): korrekt modifiziertes Protein ab Stufe 1"
    return [_make("ptm_eukaryotic_host", lang, effort="high",
                  current_state=cur, optimized_state=opt)]


# ---------------------------------------------------------------------------
# Rule implementations
# ---------------------------------------------------------------------------

def _rule_switch_chem_to_biotech(p: Dict[str, Any], lang: str) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    if (p.get("method") or "").lower() != "chemical":
        return out
    mtype = (p.get("molecule_type") or "").lower()
    msub = (p.get("molecule_subtype") or "").lower()
    if mtype == "natural_product" or msub == "terpene" or (p.get("molecule_name") or "").lower() in {
        "vanillin", "linalool", "citral", "limonene", "geraniol", "ethanol", "glutathione"
    }:
        mname = p.get("molecule_name") or ""
        if lang == "en":
            cur = f"Chemical synthesis for {mname} — typically requires toxic/halogenated solvents, lower scale-up flexibility"
            opt = "Fermentative route via engineered yeast / E. coli: 60–90 % less halogenated solvent use, comparable or higher yield with optimized strain, 6–18 months time-to-market for strain optimization"
        else:
            cur = f"Chemische Synthese für {mname} — erfordert typ. toxische/halogenierte Lösungsmittel, geringere Skalierbarkeit"
            opt = "Fermentative Route via engineered Hefe / E. coli: 60–90 % weniger halogenierte Lösungsmittel, bei optimierten Stämmen vergleichbarer oder höherer Yield, Time-to-Market 6–18 Monate für Stammoptimierung"
        out.append(_make("switch_chem_to_biotech", lang, effort="high",
                         current_state=cur, optimized_state=opt))
    return out


def _rule_step_consolidation(p: Dict[str, Any], lang: str) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    steps = p.get("number_of_steps")
    if not p.get("has_existing_process") or not steps:
        return out

    # Bug D: SPPS (Solid-Phase Peptide Synthesis) detection.
    # When a peptide is produced via chemical synthesis with many steps, the
    # step count = number of amino-acid couplings = sequence length, which is
    # a hard biological constraint. Suggesting "reduce 30 → 28 steps" is
    # nonsensical there. The right advice is to switch to recombinant
    # production or evaluate a shorter analogue, not to telescope steps.
    mtype = (p.get("molecule_type") or "").lower()
    method = (p.get("method") or "").lower()
    is_spps = (mtype == "peptide" and method == "chemical" and int(steps) >= 15)
    if is_spps:
        # Emit a route-switch suggestion instead of a step-count consolidation.
        if lang == "en":
            cur = f"SPPS with {int(steps)} couplings — step count is fixed by the peptide sequence; each coupling loses 1–5 % yield (overall yield typically 30–60 %)"
            opt = "Recombinant production (E. coli / yeast fusion expression with protease cleavage): 80–95 % overall yield, COGS reduction of 50–80 % for sequences > 15 aa, time-to-market 12–24 months"
        else:
            cur = f"SPPS mit {int(steps)} Kupplungen — Schrittzahl ist durch die Peptidsequenz vorgegeben; jede Kupplung verliert 1–5 % Yield (Gesamtausbeute typ. 30–60 %)"
            opt = "Rekombinante Produktion (E. coli / Hefe-Fusionsexpression mit Protease-Spaltung): 80–95 % Gesamtausbeute, COGS-Senkung von 50–80 % für Sequenzen > 15 AS, Time-to-Market 12–24 Monate"
        out.append(_make("spps_to_recombinant", lang, effort="high",
                         current_state=cur, optimized_state=opt,
                         steps=int(steps)))
        return out

    if int(steps) >= 4:
        target = max(int(steps) - 2, 2)
        gain_pct = round((1.15 ** 2 - 1) * 100)
        if lang == "en":
            cur = f"{int(steps)} synthesis steps — each step loses typically 10–25 % yield, drives reagent/labour cost"
            opt = f"With {target} steps via telescoping / one-pot / protecting-group avoidance: ≈ {gain_pct} % additional overall yield, solvent and labour costs drop proportionally"
        else:
            cur = f"{int(steps)} Syntheseschritte — jeder Schritt verliert typ. 10–25 % Yield, treibt Reagenzien-/Personalkosten"
            opt = f"Mit {target} Schritten via Telescoping / One-Pot / Vermeidung von Schutzgruppen: ≈ {gain_pct} % zusätzlicher Gesamt-Yield, Lösungsmittel- und Personalkosten sinken proportional"
        out.append(_make("step_consolidation", lang, effort="medium",
                         current_state=cur, optimized_state=opt,
                         steps=int(steps), target_steps=target, pct=gain_pct))
    return out


def _rule_purity_gap(p: Dict[str, Any], lang: str) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    pct = p.get("desired_purity_percent")
    methods = p.get("available_methods") or {}
    if isinstance(pct, (int, float)):
        high_purity = pct >= 99.0
        purity_str = f"{pct:.2f} %"
    else:
        high_purity = (p.get("desired_purity") or "").lower() in (">99%", "very high")
        purity_str = "> 99 %"
    if high_purity and not (
        methods.get("has_prep_hplc") or methods.get("has_fplc") or methods.get("has_crystallization")
    ):
        if lang == "en":
            cur = f"Target {purity_str}, no high-resolution method available (no prep-HPLC, FPLC or crystallization)"
            opt = "With prep-HPLC (~€50–150 k investment) or crystallization: target purity reliably reachable, GMP-pathway possible"
        else:
            cur = f"Zielreinheit {purity_str}, keine hochauflösende Methode verfügbar (keine prep-HPLC, FPLC oder Kristallisation)"
            opt = "Mit prep-HPLC (~50–150 k€ Investition) oder Kristallisation: Zielreinheit zuverlässig erreichbar, GMP-Pfad möglich"
        out.append(_make("purity_gap", lang, effort="high",
                         current_state=cur, optimized_state=opt,
                         purity_str=purity_str))
    return out


def _rule_industrial_no_bioreactor(p: Dict[str, Any], lang: str) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    scale = (p.get("scale") or "").lower()
    method = (p.get("method") or "").lower()
    if scale == "industrial" and method == "biotechnological" and not p.get("has_bioreactor"):
        kg = p.get("scale_kg_per_year")
        unit = "kg/year" if lang == "en" else "kg/Jahr"
        scale_qty = f" ({kg:.0f} {unit})" if isinstance(kg, (int, float)) else ""
        kg_str = f"{kg:.0f} {unit}" if isinstance(kg, (int, float)) else (
            "industrial scale" if lang == "en" else "industrieller Maßstab")
        if lang == "en":
            cur = f"Industrial biotechnological production at {kg_str} without bioreactor in-house"
            opt = "With ≥ 1 m³ stainless-steel fermenter or single-use bioreactor: planned scale reachable, GMP-validatable; alternatively CDMO outsourcing as bridge"
        else:
            cur = f"Industrielle biotechnologische Produktion bei {kg_str} ohne eigenen Bioreaktor"
            opt = "Mit ≥ 1 m³ Stainless-Steel-Fermenter oder Single-Use-Bioreaktor: geplante Skala erreichbar, GMP-validierbar; alternativ CDMO-Auslagerung als Brücke"
        out.append(_make("industrial_no_bioreactor", lang, effort="high",
                         current_state=cur, optimized_state=opt,
                         scale_qty=scale_qty))
    return out


def _rule_strict_waste(p: Dict[str, Any], lang: str) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    if not p.get("strict_waste_constraints"):
        return out
    if (p.get("method") or "").lower() == "chemical":
        if lang == "en":
            cur = "Strict waste constraints + chemical synthesis — halogenated / highly toxic solvents (DCM, DMF, pyridine) drive disposal cost and regulatory burden"
            opt = "With CHEM21 solvent selection (green replacements, recycling): 30–60 % lower waste-disposal cost, improved ESG and GMP compliance"
        else:
            cur = "Strenge Abfallauflagen + chemische Synthese — halogenierte / stark toxische Lösungsmittel (DCM, DMF, Pyridin) treiben Entsorgungs- und Compliance-Kosten"
            opt = "Mit CHEM21-Lösungsmittel-Auswahl (grüne Alternativen, Recycling): 30–60 % geringere Entsorgungskosten, verbesserte ESG- und GMP-Compliance"
        out.append(_make("strict_waste_chemical", lang, effort="medium",
                         current_state=cur, optimized_state=opt))
    return out


def _rule_raw_materials(p: Dict[str, Any], lang: str) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    avail = (p.get("raw_material_availability") or "").lower()
    cost = (p.get("raw_material_cost") or "").lower()
    n_supp = p.get("num_qualified_suppliers")
    lead = p.get("lead_time_weeks")
    single_region = bool(p.get("single_region_concentration"))

    # Single-source qualification — most concrete and high-leverage action
    if isinstance(n_supp, (int, float)) and int(n_supp) <= 1:
        if isinstance(lead, (int, float)):
            lead_str = (
                f", Lieferzeit {lead:.1f} Wochen" if lang == "de"
                else f", lead time {lead:.1f} weeks"
            )
        else:
            lead_str = ""
        if lang == "en":
            cur = f"1 qualified supplier listed{lead_str} — single point of failure"
            opt = "With 2–3 qualified suppliers: GMP audit-trail enabled, 5–15 % better negotiation position, no production halt on supplier outage"
        else:
            cur = f"1 qualifizierter Lieferant gelistet{lead_str} — Single-Point-of-Failure"
            opt = "Mit 2–3 qualifizierten Lieferanten: GMP-Audit-Trail ermöglicht, 5–15 % bessere Verhandlungsposition, keine Produktionsunterbrechung bei Lieferantenausfall"
        out.append(_make("single_source_dual", lang, effort="medium",
                         current_state=cur, optimized_state=opt,
                         lead_str=lead_str))

    # Geographic diversification when single-region
    if single_region:
        if lang == "en":
            cur = "All qualified suppliers in a single region (geopolitical / logistics single point of failure)"
            opt = "Second supplier in alternative region: protection against export restrictions, sanctions, regional outages — resilient supply chain"
        else:
            cur = "Alle qualifizierten Lieferanten in einer Region (Geopolitik-/Logistik-Single-Point-of-Failure)"
            opt = "Zweiter Lieferant in alternativer Region: Schutz vor Export-Restriktionen, Sanktionen, regionalen Ausfällen — resilientere Supply Chain"
        out.append(_make("geo_diversification", lang, effort="medium",
                         current_state=cur, optimized_state=opt))

    # Long lead time — different action than supplier diversification
    if isinstance(lead, (int, float)) and float(lead) >= 8.0 and not (
        isinstance(n_supp, (int, float)) and int(n_supp) <= 1
    ):
        if lang == "en":
            cur = f"Lead time {lead:.1f} weeks — requires large safety stock (capital lock-up) or forward-pay commitments"
            opt = "With VMI (Vendor-Managed Inventory) or consignment stock: 30–60 % less capital tied up at the same supply security"
        else:
            cur = f"Lieferzeit {lead:.1f} Wochen — erfordert hohe Sicherheitsbestände (Kapitalbindung) oder Vorab-Forecast-Commitments"
            opt = "Mit Vendor-Managed-Inventory (VMI) oder Konsignationslager: 30–60 % weniger Kapitalbindung bei gleicher Versorgungssicherheit"
        out.append(_make("long_lead_time", lang, effort="low",
                         current_state=cur, optimized_state=opt,
                         lead=float(lead)))

    # Legacy fallback when availability=low but no specific signal triggered
    if avail == "low" and not (isinstance(n_supp, (int, float)) and int(n_supp) <= 1) and not single_region:
        if lang == "en":
            cur = "Low market availability of raw materials — supply-security single point of failure"
            opt = "With qualified second supplier or alternative starting-material path: protection against production outages, negotiation leverage on prices"
        else:
            cur = "Niedrige Marktverfügbarkeit der Rohstoffe — Single-Point-of-Failure für die Versorgungssicherheit"
            opt = "Mit qualifiziertem Zweitlieferant oder alternativem Edukt-Pfad: Schutz vor Produktionsausfällen, Verhandlungsspielraum bei Preisen"
        out.append(_make("raw_material_legacy", lang, effort="medium",
                         current_state=cur, optimized_state=opt))

    if cost == "high":
        eur_per_kg = p.get("raw_material_cost_eur_per_kg")
        if isinstance(eur_per_kg, (int, float)):
            cost_qty = (
                f" (aktuell {eur_per_kg:.2f} €/kg)" if lang == "de"
                else f" (currently {eur_per_kg:.2f} €/kg)"
            )
            if lang == "en":
                cur = f"Starting-material price {eur_per_kg:.2f} €/kg — dominates manufacturing cost"
                opt = "With cheaper precursor synthesis or solvent / reagent recycling: 20–50 % manufacturing-cost reduction at substantial substitution"
            else:
                cur = f"Edukt-Preis {eur_per_kg:.2f} €/kg — dominiert die Herstellkosten"
                opt = "Mit günstigerer Vorstufen-Synthese oder Lösungsmittel-/Reagenzien-Recycling: 20–50 % Senkung der Herstellkosten bei substantieller Substitution"
        else:
            cost_qty = ""
            cur = ""
            opt = ""
        out.append(_make("edukt_substitution", lang, effort="medium",
                         current_state=cur, optimized_state=opt,
                         cost_qty=cost_qty))

    return out


def _rule_biotech_titer_optimization(p: Dict[str, Any], lang: str) -> List[Dict[str, Any]]:
    """Bug F fallback: biomolecule production with biotech method has implicit cost
    pressure (bioreactor, media, downstream). When NO other rule produces an
    action but cost signals are present, suggest titer/yield optimization —
    the most universally applicable improvement for fermentation-based COGS.

    Triggers when:
    - molecule is peptide or protein
    - method is biotechnological
    - scale is at least pilot (i.e. not pure lab-scale R&D)
    - at least one of: raw_cost ≥ medium, purity ≥ 99 %, scale = industrial,
      multi-domain protein, or has_ptm
    """
    out: List[Dict[str, Any]] = []
    mtype = (p.get("molecule_type") or "").lower()
    if mtype not in ("peptide", "protein"):
        return out
    method = (p.get("method") or "").lower()
    if method not in ("biotechnological", "biotech"):
        return out
    scale = (p.get("scale") or "").lower()
    if scale not in ("pilot", "industrial"):
        return out
    raw_cost = (p.get("raw_material_cost") or "").lower()
    desired_purity = (p.get("desired_purity") or "").lower()
    has_ptm = bool(p.get("has_ptm"))
    num_domains = int(p.get("num_domains") or 0)
    eur_per_kg = p.get("raw_material_cost_eur_per_kg")
    cost_signal = (
        raw_cost in ("medium", "high")
        or desired_purity in (">99%", "very high")
        or scale == "industrial"
        or has_ptm
        or num_domains >= 2
        or (isinstance(eur_per_kg, (int, float)) and eur_per_kg >= 50.0)
    )
    if not cost_signal:
        return out

    if lang == "en":
        cur = "Biotech production at pilot/industrial scale — bioreactor footprint, media and downstream dominate COGS; baseline titer not optimised"
        opt = "Strain engineering + fed-batch/perfusion + media DoE → typically 2× titer in 6–12 months, ~30–50 % COGS reduction"
    else:
        cur = "Biotechnologische Produktion im Pilot-/Industriemaßstab — Bioreaktor-Footprint, Medien und Aufarbeitung dominieren die Herstellkosten; Basis-Titer nicht optimiert"
        opt = "Stamm-Engineering + Fed-Batch/Perfusion + Medien-DoE → typ. Titer-Verdopplung in 6–12 Monaten, ~30–50 % Senkung der COGS"
    out.append(_make("biotech_titer_optimization", lang, effort="high",
                     current_state=cur, optimized_state=opt))
    return out


def build_recommended_actions(process_input: Dict[str, Any], lang: str = "de") -> List[Dict[str, Any]]:
    """Apply all rules and return the deduplicated list of recommended actions."""
    if not isinstance(process_input, dict):
        return []
    if lang not in ("de", "en"):
        lang = "de"
    actions: List[Dict[str, Any]] = []
    for rule in (
        _rule_switch_chem_to_biotech,
        _rule_step_consolidation,
        _rule_purity_gap,
        _rule_industrial_no_bioreactor,
        _rule_strict_waste,
        _rule_raw_materials,
        # Biophysical numeric-driven rules
        _rule_biophys_low_tm,
        _rule_biophys_high_tm,
        _rule_biophys_low_tagg,
        _rule_biophys_multi_domain,
        _rule_biophys_ptm,
        # Bug F: catch-all for biomolecules with cost pressure
        _rule_biotech_titer_optimization,
    ):
        try:
            actions.extend(rule(process_input, lang))
        except Exception:
            continue
    # Deduplicate by title
    seen = set()
    unique: List[Dict[str, Any]] = []
    for a in actions:
        t_ = a.get("title")
        if t_ in seen:
            continue
        seen.add(t_)
        unique.append(a)
    return unique
