"""Plausibility checks for user-provided process descriptions.

Produces a list of advisory messages (warnings, not hard errors) that flag
implausible inputs based on a small curated knowledge base of known molecules
and general process heuristics.

The checks are intentionally conservative: they raise a warning only when the
deviation from known practice is large. The user can always override and
proceed — this layer only nudges them to double-check.

All output messages are returned in the requested language (de | en).
"""

from typing import Dict, Any, List


# ---------------------------------------------------------------------------
# Curated knowledge base
# ---------------------------------------------------------------------------

# Per-molecule notes are language-tagged. Step ranges are language-independent.
KNOWN_PROCESSES: Dict[str, Dict[str, Dict[str, Any]]] = {
    "vanillin": {
        "chemical": {"steps": (1, 3),
                     "note": {"de": "Guaiacol → Vanillin via Reimer-Tiemann oder Glyoxylsäure-Route.",
                              "en": "Guaiacol → vanillin via Reimer-Tiemann or glyoxylic-acid route."},
                     "literature_source": "Fache et al., Acta Chim. Slov. 2000 (vanillin synthesis routes overview); Borregaard / Rhodia industrial process"},
        "biotechnological": {"steps": (2, 5),
                             "note": {"de": "Glucose oder Ferulasäure → Vanillin via S. cerevisiae oder A. niger.",
                                      "en": "Glucose or ferulic acid → vanillin via S. cerevisiae or A. niger."},
                             "literature_source": "Hansen et al., Appl. Environ. Microbiol. 2009 (doi:10.1128/AEM.02074-08)"},
        "extraction": {"steps": (1, 2),
                       "note": {"de": "Aus Vanilleschoten — sehr niedrige Ausbeute (~2 % w/w).",
                                "en": "From vanilla pods — very low yield (~2 % w/w)."},
                       "literature_source": "Sinha et al., Crit. Rev. Food Sci. Nutr. 2008 (doi:10.1080/10408390701764235)"},
    },
    "aspirin": {
        "chemical": {"steps": (1, 2),
                     "note": {"de": "Salicylsäure + Essigsäureanhydrid → Aspirin (1-Schritt-Synthese).",
                              "en": "Salicylic acid + acetic anhydride → aspirin (one-step synthesis)."},
                     "literature_source": "Vollhardt & Schore, Organic Chemistry, 8th ed. 2018, Macmillan/Freeman (ISBN 978-1-319-07945-1)"},
    },
    "ibuprofen": {
        "chemical": {"steps": (3, 6),
                     "note": {"de": "BHC- oder Boots-Verfahren, typisch 3–6 Schritte.",
                              "en": "BHC or Boots process, typically 3–6 steps."},
                     "literature_source": "BHC Boots-Hoechst-Celanese process — US Patent 4,981,995 (Elango & Davenport 1991); Presidential Green Chemistry Award 1997"},
    },
    "ethanol": {
        "biotechnological": {"steps": (1, 2),
                             "note": {"de": "Hefe-Fermentation, Standard-Bioprozess.",
                                      "en": "Yeast fermentation, standard bioprocess."},
                             "literature_source": "Doran, Bioprocess Engineering Principles, 2nd ed. 2013, Academic Press (ISBN 978-0-12-220851-5)"},
        "chemical": {"steps": (1, 2),
                     "note": {"de": "Ethylen-Hydratation.", "en": "Ethylene hydration."},
                     "literature_source": "Ullmann's Encyclopedia of Industrial Chemistry — Ethanol (Wiley-VCH, ISBN 978-3-527-30385-4; Online: doi:10.1002/14356007.a09_587.pub2)"},
    },
    "penicillin": {
        "biotechnological": {"steps": (1, 3),
                             "note": {"de": "Fermentation mit P. chrysogenum, anschließende Aufarbeitung.",
                                      "en": "Fermentation with P. chrysogenum, followed by downstream workup."},
                             "literature_source": "Walsh G., Pharmaceutical Biotechnology, 2nd ed. 2018, Wiley (ISBN 978-1-119-11518-7); Elander, Appl. Microbiol. Biotechnol. 2003 (doi:10.1007/s00253-003-1274-y)"},
    },
    "linalool": {
        "extraction": {"steps": (1, 2),
                       "note": {"de": "Aus ätherischen Ölen (Lavendel, Koriander).",
                                "en": "From essential oils (lavender, coriander)."},
                       "literature_source": "Lapczynski et al., Food Chem. Toxicol. 2008 — linalool review (doi:10.1016/j.fct.2008.06.025)"},
        "biotechnological": {"steps": (2, 4),
                             "note": {"de": "Engineered S. cerevisiae oder E. coli Stämme.",
                                      "en": "Engineered S. cerevisiae or E. coli strains."},
                             "literature_source": "Amiri et al., Curr. Pharm. Biotechnol. 2014 — microbial terpene production review (doi:10.2174/138920101504140825195713)"},
    },
    "citral": {
        "extraction": {"steps": (1, 2),
                       "note": {"de": "Aus Lemongrass-Öl.", "en": "From lemongrass oil."},
                       "literature_source": "Schaneberg & Khan, J. Agric. Food Chem. 2002 (doi:10.1021/jf015818g)"},
        "chemical": {"steps": (3, 6),
                     "note": {"de": "Klassische Synthese aus Isobuten und Formaldehyd.",
                              "en": "Classical synthesis from isobutene and formaldehyde."},
                     "literature_source": "BASF Citral Process — Ullmann's Encyclopedia of Industrial Chemistry (Wiley-VCH, ISBN 978-3-527-30385-4)"},
    },
    "limonene": {
        "extraction": {"steps": (1, 2),
                       "note": {"de": "Aus Citrusschalen — Massenrohstoff.",
                                "en": "From citrus peels — bulk feedstock."},
                       "literature_source": "Ciriminna et al., Chem. Eur. J. 2014 — limonene from citrus processing (doi:10.1002/chem.201400290)"},
    },
    "cyclosporine": {
        "biotechnological": {"steps": (2, 4),
                             "note": {"de": "Fermentation mit Tolypocladium inflatum, gefolgt von chromatographischer Reinigung.",
                                      "en": "Fermentation with Tolypocladium inflatum, followed by chromatographic purification."},
                             "literature_source": "Borel et al., Agents Actions 1976 — discovery (doi:10.1007/BF01972257); Survase et al., Crit. Rev. Biotechnol. 2011 — production review (doi:10.3109/07388551.2010.497461)"},
    },
    "glutathione": {
        "biotechnological": {"steps": (2, 4),
                             "note": {"de": "Fermentative Produktion in Hefe.",
                                      "en": "Fermentative production in yeast."},
                             "literature_source": "Li et al., Appl. Microbiol. Biotechnol. 2004 — Glutathione fermentation (doi:10.1007/s00253-004-1559-9)"},
        "chemical": {"steps": (3, 8),
                     "note": {"de": "SPPS möglich, aber für 3-AS-Peptid selten kosteneffizient.",
                              "en": "SPPS feasible, but rarely cost-efficient for a 3-AA peptide."},
                     "literature_source": "Chan & White, Fmoc Solid Phase Peptide Synthesis, Oxford Univ. Press 2000 (ISBN 978-0-19-963724-9)"},
    },
    "insulin": {
        "biotechnological": {"steps": (3, 6),
                             "note": {"de": "Rekombinant in E. coli oder Hefe + komplexe Aufarbeitung.",
                                      "en": "Recombinant in E. coli or yeast + complex downstream processing."},
                             "literature_source": "Walsh G., Nat. Biotechnol. 2005 — insulin production review (doi:10.1038/nbt0405-401); Walsh G., Pharmaceutical Biotechnology, 2nd ed. 2018, Wiley (ISBN 978-1-119-11518-7)"},
    },
    # ---------------- Extension Set ----------------
    "citric acid": {
        "biotechnological": {"steps": (2, 4),
                             "note": {"de": "A. niger Submersfermentation, klassische Calciumcitrat-Fällung.",
                                      "en": "A. niger submerged fermentation, classical calcium-citrate precipitation."},
                             "literature_source": "Berovic & Legisa, Biotechnol. Annu. Rev. 2007 (doi:10.1016/S1387-2656(07)13011-8)"},
    },
    "lactic acid": {
        "biotechnological": {"steps": (2, 4),
                             "note": {"de": "Lactobacillus / B. coagulans Fermentation; Aufarbeitung kritisch (chiral).",
                                      "en": "Lactobacillus / B. coagulans fermentation; downstream is critical (chiral purity)."},
                             "literature_source": "Abdel-Rahman et al., Biotechnol. Adv. 2013 (doi:10.1016/j.biotechadv.2013.04.002)"},
    },
    "succinic acid": {
        "biotechnological": {"steps": (2, 4),
                             "note": {"de": "Anaerobe Fermentation mit CO2-Fixierung; Direktkristallisation möglich.",
                                      "en": "Anaerobic fermentation with CO2 fixation; direct crystallization feasible."},
                             "literature_source": "Cok et al., Biofuels Bioprod. Bioref. 2014 (doi:10.1002/bbb.1427)"},
    },
    "itaconic acid": {
        "biotechnological": {"steps": (2, 3),
                             "note": {"de": "Aspergillus terreus Fermentation; relativ einfache Aufarbeitung.",
                                      "en": "Aspergillus terreus fermentation; relatively simple downstream."},
                             "literature_source": "Klement & Büchs, Bioresour. Technol. 2013 (doi:10.1016/j.biortech.2013.02.054)"},
    },
    "paracetamol": {
        "chemical": {"steps": (1, 2),
                     "note": {"de": "Acetylierung von p-Aminophenol mit Essigsäureanhydrid (1-Schritt-Synthese).",
                              "en": "Acetylation of p-aminophenol with acetic anhydride (one-step synthesis)."},
                     "literature_source": "Joncour et al., Green Chem. 2014 (doi:10.1039/C4GC00603H); Ellis F., Paracetamol: A Curriculum Resource, RSC 2002 (ISBN 978-0-85404-365-1)"},
    },
    "glucoamylase": {
        "biotechnological": {"steps": (2, 4),
                             "note": {"de": "A. niger Fermentation + UF-Aufkonzentrierung; meist als Multikomponenten-Formulierung verkauft.",
                                      "en": "A. niger fermentation + UF concentration; typically sold as multi-component formulation."},
                             "literature_source": "Norouzian et al., Biotechnol. Adv. 2006 (doi:10.1016/j.biotechadv.2005.06.003)"},
    },
    "glucose isomerase": {
        "biotechnological": {"steps": (3, 5),
                             "note": {"de": "Streptomyces-Fermentation + Immobilisierung — Immobilisierung ist der Kostentreiber, nicht die Fermentation.",
                                      "en": "Streptomyces fermentation + immobilization — immobilization is the cost driver, not fermentation."},
                             "literature_source": "Bhosale et al., Microbiol. Rev. 1996 (doi:10.1128/mr.60.2.280-300.1996)"},
    },
    "asparaginase": {
        "biotechnological": {"steps": (4, 6),
                             "note": {"de": "E. coli Fermentation, mehrstufige Reinigung (IEX + SEC); klinisch / Pharma-Grade.",
                                      "en": "E. coli fermentation, multi-step purification (IEX + SEC); clinical / pharma grade."},
                             "literature_source": "Kotzia & Labrou, J. Biotechnol. 2007 (doi:10.1016/j.jbiotec.2007.07.939)"},
    },
    "pembrolizumab": {
        "biotechnological": {"steps": (5, 8),
                             "note": {"de": "CHO-Zellkultur + Standard-Antikörper-Aufarbeitung (Protein A → IEX → VRF → TFF).",
                                      "en": "CHO cell culture + standard antibody downstream (Protein A → IEX → VRF → TFF)."},
                             "literature_source": "Liu et al., mAbs 2010 (doi:10.4161/mabs.2.5.12796); Walsh G., Nat. Biotechnol. 2018 (doi:10.1038/nbt.4313)"},
    },
    "β-carotene": {
        "biotechnological": {"steps": (3, 5),
                             "note": {"de": "Blakeslea trispora oder engineered Yarrowia + Lipid-Extraktion + Kristallisation.",
                                      "en": "Blakeslea trispora or engineered Yarrowia + lipid extraction + crystallization."},
                             "literature_source": "Mata-Gómez et al., Microb. Cell Fact. 2014 (doi:10.1186/1475-2859-13-12)"},
    },
    "astaxanthin": {
        "biotechnological": {"steps": (3, 5),
                             "note": {"de": "Algenkultur (Haematococcus) oder Phaffia-Fermentation; Zellaufschluss bei Algen anspruchsvoll.",
                                      "en": "Algal cultivation (Haematococcus) or Phaffia fermentation; algal cell lysis is the bottleneck."},
                             "literature_source": "Shah et al., Front. Plant Sci. 2016 (doi:10.3389/fpls.2016.00531)"},
    },
    "caffeine": {
        "extraction": {"steps": (2, 4),
                       "note": {"de": "Aus Kaffeebohnen-Entkoffeinierung als Nebenprodukt; CO2-überkritisch oder klassisch wässrig.",
                                "en": "By-product of coffee-bean decaffeination; supercritical CO2 or classical aqueous."},
                       "literature_source": "Heilmann W., Lebensmittelchemie 2001 (DLG-Verlag); Ramalakshmi & Raghavan, Crit. Rev. Food Sci. Nutr. 1999 (doi:10.1080/10408699991279231)"},
    },
    "morphine": {
        "extraction": {"steps": (3, 5),
                       "note": {"de": "Aus Mohnstroh oder Latex; mehrstufige Extraktion und pH-Cycling, dann Salzbildung als API.",
                                "en": "From poppy straw or latex; multi-step extraction and pH cycling, then salt formation as API."},
                       "literature_source": "Kutchan, The Alkaloids 1998 (doi:10.1016/S0099-9598(08)60052-6); Beaudoin & Facchini, Planta 2014 (doi:10.1007/s00425-014-2056-8)"},
    },
}


# ---------------------------------------------------------------------------
# Localized message templates
# ---------------------------------------------------------------------------

MSG_BIOPHYS = {
    "tm_below_tagg": {
        "de": ("Inkonsistent: Schmelzpunkt Tm = {tm:.1f} °C liegt unter dem Aggregations-Onset "
               "Tagg = {tagg:.1f} °C. Üblicherweise gilt Tm ≥ Tagg, weil Aggregation typisch "
               "der Entfaltung folgt. Bitte Werte prüfen."),
        "en": ("Inconsistent: melting temperature Tm = {tm:.1f} °C is below aggregation onset "
               "Tagg = {tagg:.1f} °C. Usually Tm ≥ Tagg because aggregation typically follows "
               "unfolding. Please review the values."),
    },
    "antibody_no_disulfides": {
        "de": ("Antikörper mit 0 Disulfidbrücken angegeben — IgG-Antikörper haben standardmäßig "
               "16 Disulfidbrücken (12 intrachain + 4 interchain). Bitte Wert prüfen."),
        "en": ("Antibody specified with 0 disulfide bridges — IgG antibodies have 16 disulfide "
               "bridges by default (12 intrachain + 4 interchain). Please verify."),
    },
    "antibody_too_few_domains": {
        "de": ("Antikörper mit nur {domains} Domäne(n) angegeben — IgG-Antikörper haben "
               "standardmäßig 12 Domänen (4 schwere × 2 Untereinheiten + 2 leichte × 2). "
               "Möglicherweise Antikörper-Fragment (Fab, scFv) gemeint?"),
        "en": ("Antibody specified with only {domains} domain(s) — IgG antibodies have "
               "12 domains by default (4 heavy × 2 subunits + 2 light × 2). "
               "Did you mean an antibody fragment (Fab, scFv)?"),
    },
    "very_low_tm": {
        "de": ("Sehr niedrige Schmelzpunkt Tm = {tm:.1f} °C — Tiefkühlpflicht (≤ –20 °C) und "
               "Cryoprotectant-Formulierung sind erforderlich. Lyophilisation ohne Schutzstoffe ist nicht möglich."),
        "en": ("Very low melting temperature Tm = {tm:.1f} °C — deep-freeze storage (≤ –20 °C) "
               "and cryoprotectant formulation required. Lyophilization without protectants is not feasible."),
    },
    "ptm_chemical": {
        "de": ("Posttranslationale Modifikationen wurden markiert, aber chemische Synthese ist "
               "ausgewählt. Glykosylierung und ähnliche PTMs sind chemisch nicht oder nur "
               "extrem aufwändig durchführbar. Bitte biotechnologische Methode wählen."),
        "en": ("Post-translational modifications are checked, but chemical synthesis is selected. "
               "Glycosylation and similar PTMs are not feasible chemically (or only with extreme "
               "effort). Please switch to a biotechnological method."),
    },
}


MSG = {
    "terpene_chemical": {
        "de": ("Terpene werden in der industriellen Praxis fast immer durch Extraktion oder "
               "biotechnologische Produktion (engineered Hefen/E. coli) hergestellt; "
               "rein chemische Synthese ist meist unwirtschaftlich."),
        "en": ("Terpenes are almost always produced by extraction or biotechnological routes "
               "(engineered yeast/E. coli) at industrial scale; purely chemical synthesis is "
               "usually uneconomical."),
    },
    "natural_product_industrial_chemical": {
        "de": ("Industrielle Chemo-Synthese von Naturstoffen ist selten konkurrenzfähig zur Extraktion "
               "oder biotechnologischen Produktion — bitte Wirtschaftlichkeit prüfen."),
        "en": ("Industrial chemical synthesis of natural products is rarely competitive with extraction "
               "or biotechnological routes — please review the business case."),
    },
    "step_too_low": {
        "de": "Sie haben {steps} Schritt(e) angegeben für {name} via {method}. Realistisch sind {lo}–{hi} Schritte. {note}",
        "en": "You entered {steps} step(s) for {name} via {method}. Realistic range: {lo}–{hi} steps. {note}",
    },
    "step_too_high": {
        "de": "Sie haben {steps} Schritte angegeben für {name} via {method}. Üblich sind {lo}–{hi} Schritte — sind Aufarbeitungs-/Schutzgruppenschritte mitgezählt?",
        "en": "You entered {steps} steps for {name} via {method}. Typical range: {lo}–{hi} steps — are workup/protecting-group steps counted?",
    },
    "industrial_no_bioreactor": {
        "de": ("Industrielle biotechnologische Produktion{scale_qty} ohne Bioreaktor ist nicht realistisch — "
               "ein Fermenter ≥ 1 L (für ≥ 1.000 kg/Jahr typ. ≥ 1 m³) ist Voraussetzung."),
        "en": ("Industrial biotechnological production{scale_qty} without a bioreactor is not realistic — "
               "a fermenter ≥ 1 L (for ≥ 1,000 kg/year typically ≥ 1 m³) is required."),
    },
    "high_purity_no_method": {
        "de": ("Zielreinheit {purity_str} verlangt typisch prep-HPLC, FPLC oder Umkristallisation — "
               "keine dieser Methoden ist als verfügbar markiert."),
        "en": ("Target purity {purity_str} typically requires prep-HPLC, FPLC, or recrystallization — "
               "none of these methods is marked as available."),
    },
    "single_source": {
        "de": ("Single-Source-Risiko: nur 1 qualifizierter Lieferant. Bei Lieferantenausfall steht die gesamte "
               "Produktion still. Qualifizierung eines Zweitlieferanten ist Pflicht für GMP-Produktion und "
               "stark empfohlen darunter."),
        "en": ("Single-source risk: only 1 qualified supplier. A supplier outage halts the entire production. "
               "Qualifying a second supplier is mandatory for GMP production and strongly recommended otherwise."),
    },
    "long_lead": {
        "de": ("Sehr lange Lieferzeit ({lead:.1f} Wochen): erfordert hohe Sicherheitsbestände und macht "
               "Bedarfsplanung anfällig für Forecast-Fehler. Lange Lead-Times deuten oft auf Auftragssynthese "
               "hin — Preisrisiko prüfen."),
        "en": ("Very long lead time ({lead:.1f} weeks): requires high safety stock and makes demand planning "
               "prone to forecast errors. Long lead times often indicate custom synthesis — review price risk."),
    },
    "single_region_few_suppliers": {
        "de": ("Geografische Konzentration kombiniert mit nur 1–2 Lieferanten = kumuliertes Versorgungsrisiko "
               "(Geopolitik, regionale Ausfälle, Logistik). Diversifikation in eine zweite Region wird empfohlen."),
        "en": ("Geographic concentration combined with only 1–2 suppliers = compounded supply risk "
               "(geopolitics, regional outages, logistics). Diversification into a second region is recommended."),
    },
    "industrial_long_lead": {
        "de": ("Industrieller Maßstab mit {lead:.1f} Wochen Lieferzeit ist operational schwierig: erfordert "
               "≥ 3 Monate Sicherheitsbestand oder strategische Vorab-Vereinbarungen mit Lieferanten."),
        "en": ("Industrial scale with {lead:.1f} weeks lead time is operationally difficult: requires "
               "≥ 3 months safety stock or strategic forward agreements with suppliers."),
    },
    "api_purity_single_source": {
        "de": ("Pharma-/API-grade Reinheit aus Single-Source-Lieferung ist regulatorisch schwierig: GMP-Audit-Trail "
               "typischerweise mit mind. einem qualifizierten Backup-Lieferanten gefordert."),
        "en": ("Pharma/API-grade purity with a single-source supply is regulatorily challenging: GMP audit trails "
               "typically require at least one qualified backup supplier."),
    },
}


def _m(key: str, lang: str, **kwargs) -> str:
    bundle = MSG.get(key) or MSG_BIOPHYS.get(key) or {}
    template = bundle.get(lang) or bundle.get("de") or key
    if kwargs:
        try:
            return template.format(**kwargs)
        except Exception:
            return template
    return template


# --- Biophysical numeric plausibility checks ---

def _biophysical_warnings(p: Dict[str, Any], lang: str) -> List[str]:
    warnings: List[str] = []
    mtype = (p.get("molecule_type") or "").lower()
    msub = (p.get("molecule_subtype") or "").lower()
    if mtype not in ("peptide", "protein"):
        return warnings

    # Bug M4 fix: cyclic peptides (Cyclosporine class) have no tertiary
    # fold, so Tm and Tagg are conceptually not meaningful. Skip Tm/Tagg
    # plausibility checks for them — otherwise placeholder values lead
    # to self-contradictory warnings ("Tm < Tagg" for inputs the user
    # never explicitly meant to compare).
    if mtype == "peptide" and msub == "cyclic":
        # Only run the disulfide / PTM checks below, not the Tm/Tagg.
        return warnings

    tm = p.get("tm_celsius")
    tagg = p.get("tagg_celsius")
    disulfides = p.get("num_disulfides")
    domains = p.get("num_domains")
    has_ptm = bool(p.get("has_ptm"))
    method = (p.get("method") or "").lower()

    # 1. Tm below Tagg is unphysical
    if isinstance(tm, (int, float)) and isinstance(tagg, (int, float)):
        if float(tm) < float(tagg) - 2.0:  # 2 °C tolerance for measurement noise
            warnings.append(_m("tm_below_tagg", lang, tm=float(tm), tagg=float(tagg)))

    # 2. Antibody with 0 disulfides
    if msub == "antibody" and isinstance(disulfides, int) and disulfides == 0:
        warnings.append(_m("antibody_no_disulfides", lang))

    # 3. Antibody with too few domains
    if msub == "antibody" and isinstance(domains, int) and 0 < domains < 4:
        warnings.append(_m("antibody_too_few_domains", lang, domains=domains))

    # 4. Very low Tm
    if isinstance(tm, (int, float)) and float(tm) < 35.0:
        warnings.append(_m("very_low_tm", lang, tm=float(tm)))

    # 5. PTMs marked, but chemical method selected — combinatorially infeasible
    if has_ptm and method == "chemical":
        warnings.append(_m("ptm_chemical", lang))

    return warnings


# ---------------------------------------------------------------------------
# Rule implementations
# ---------------------------------------------------------------------------

def _categorical_warnings(p: Dict[str, Any], lang: str) -> List[str]:
    warnings: List[str] = []
    mtype = (p.get("molecule_type") or "").lower()
    msub = (p.get("molecule_subtype") or "").lower()
    method = (p.get("method") or "").lower()

    if msub == "terpene" and method == "chemical":
        warnings.append(_m("terpene_chemical", lang))
    if mtype == "natural_product" and method == "chemical" and (p.get("scale") or "").lower() == "industrial":
        warnings.append(_m("natural_product_industrial_chemical", lang))
    return warnings


def _step_count_warnings(p: Dict[str, Any], lang: str, sources: List[str]) -> List[str]:
    warnings: List[str] = []
    if not p.get("has_existing_process"):
        return warnings
    steps = p.get("number_of_steps")
    if steps is None:
        return warnings
    name = (p.get("molecule_name") or "").strip().lower()
    method = (p.get("method") or "").lower()
    entry = KNOWN_PROCESSES.get(name, {}).get(method)
    if not entry:
        return warnings
    lo, hi = entry["steps"]
    note = entry["note"].get(lang) or entry["note"].get("de") or ""
    pretty_name = name.title()
    triggered = False
    if int(steps) < lo:
        warnings.append(_m("step_too_low", lang, steps=int(steps), name=pretty_name, method=method, lo=lo, hi=hi, note=note))
        triggered = True
    elif int(steps) > hi + 2:
        warnings.append(_m("step_too_high", lang, steps=int(steps), name=pretty_name, method=method, lo=lo, hi=hi))
        triggered = True
    if triggered and entry.get("literature_source"):
        sources.append(entry["literature_source"])
    return warnings


def _scale_equipment_warnings(p: Dict[str, Any], lang: str) -> List[str]:
    warnings: List[str] = []
    scale = (p.get("scale") or "").lower()
    method = (p.get("method") or "").lower()
    if scale == "industrial" and method == "biotechnological" and not p.get("has_bioreactor"):
        kg = p.get("scale_kg_per_year")
        unit = "kg/year" if lang == "en" else "kg/Jahr"
        scale_qty = f" ({kg:.0f} {unit})" if isinstance(kg, (int, float)) else ""
        warnings.append(_m("industrial_no_bioreactor", lang, scale_qty=scale_qty))
    return warnings


def _purity_method_warnings(p: Dict[str, Any], lang: str) -> List[str]:
    warnings: List[str] = []
    pct = p.get("desired_purity_percent")
    methods = p.get("available_methods") or {}
    high_purity = False
    if isinstance(pct, (int, float)):
        high_purity = pct >= 99.0
        purity_str = f"{pct:.2f} %"
    else:
        high_purity = (p.get("desired_purity") or "").lower() in (">99%", "very high")
        purity_str = "> 99 %"
    if high_purity:
        # Bug H6 fix: the "needs prep-HPLC / FPLC / recrystallization" warning
        # was previously fired pauschal. For some subtypes the de-facto
        # purification path is something else entirely:
        #   - small_molecule/volatile  → distillation
        #   - protein/antibody         → Protein-A + IEX + virus filtration
        #   - protein/enzyme           → UF/DF + IEX
        # Suppress the warning for those subtypes when the realistic
        # purification stack is available.
        mtype = (p.get("molecule_type") or "").lower()
        msub = (p.get("molecule_subtype") or "").lower()

        # Volatiles: distillation is the standard path → suppress
        if msub == "volatile" or mtype == "small_molecule" and msub == "volatile":
            return warnings

        # Antibodies: Protein-A + IEX + virus filtration is the standard.
        # If FPLC or membrane filtration is available (proxy for that stack),
        # or has_advanced_purification flag is set, suppress.
        if msub == "antibody":
            if (methods.get("has_fplc")
                    or methods.get("has_membrane_filtration")
                    or p.get("has_advanced_purification")):
                return warnings

        # Enzymes: UF/DF + IEX standard. Same proxy as antibody.
        if msub == "enzyme":
            if (methods.get("has_membrane_filtration")
                    or methods.get("has_fplc")
                    or p.get("has_advanced_purification")):
                return warnings

        if not (methods.get("has_prep_hplc") or methods.get("has_fplc") or methods.get("has_crystallization")):
            warnings.append(_m("high_purity_no_method", lang, purity_str=purity_str))
    return warnings


def _sourcing_warnings(p: Dict[str, Any], lang: str) -> List[str]:
    warnings: List[str] = []
    n = p.get("num_qualified_suppliers")
    lead = p.get("lead_time_weeks")
    single = bool(p.get("single_region_concentration"))
    scale = (p.get("scale") or "").lower()
    purity = (p.get("desired_purity") or "").lower()
    pct = p.get("desired_purity_percent")
    api_purity = (isinstance(pct, (int, float)) and pct >= 99.0) or purity in (">99%", "very high")

    if isinstance(n, (int, float)) and int(n) <= 1:
        warnings.append(_m("single_source", lang))
    if isinstance(lead, (int, float)) and float(lead) >= 12.0:
        warnings.append(_m("long_lead", lang, lead=float(lead)))
    if single and isinstance(n, (int, float)) and int(n) <= 2:
        warnings.append(_m("single_region_few_suppliers", lang))
    if scale == "industrial" and isinstance(lead, (int, float)) and float(lead) >= 8.0:
        warnings.append(_m("industrial_long_lead", lang, lead=float(lead)))
    if api_purity and isinstance(n, (int, float)) and int(n) <= 1:
        warnings.append(_m("api_purity_single_source", lang))
    return warnings


def check_plausibility(process_input: Dict[str, Any], lang: str = "de") -> List[str]:
    """Return a list of plausibility warnings for the given process input.

    Backward-compatible wrapper — returns only warnings, drops sources.
    Use check_plausibility_with_sources() to get warnings + literature refs.
    """
    return check_plausibility_with_sources(process_input, lang).get("warnings", [])


def check_plausibility_with_sources(process_input: Dict[str, Any], lang: str = "de") -> Dict[str, Any]:
    """Return both warnings and literature sources for the given process input.

    Output:
      {
        "warnings": [str, ...],
        "literature_sources": [str, ...]   # deduplicated, only sources whose
                                            # rule actually fired for this input
      }
    """
    if not isinstance(process_input, dict):
        return {"warnings": [], "literature_sources": []}
    if lang not in ("de", "en"):
        lang = "de"
    warnings: List[str] = []
    sources: List[str] = []
    warnings.extend(_categorical_warnings(process_input, lang))
    warnings.extend(_step_count_warnings(process_input, lang, sources))
    warnings.extend(_scale_equipment_warnings(process_input, lang))
    warnings.extend(_purity_method_warnings(process_input, lang))
    warnings.extend(_sourcing_warnings(process_input, lang))
    warnings.extend(_biophysical_warnings(process_input, lang))
    # Dedupe sources while preserving order
    seen = set()
    deduped = []
    for s in sources:
        if s not in seen:
            seen.add(s)
            deduped.append(s)
    return {"warnings": warnings, "literature_sources": deduped}
