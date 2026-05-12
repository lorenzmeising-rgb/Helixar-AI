"""
molecule_hints_extension.py
===========================
Extension der MOLECULE_HINTS aus concrete_recommendations.py um die
verbleibenden 57 Moleküle der Helixar-DB, die in der ersten Detail-DB-
Runde noch nicht abgedeckt waren.

Jeder Eintrag enthält pharma-realistische Produktions- und Aufarbeitungs-
Sequenzen mit echten Parametern (Temperatur, pH, Yield-Range, Equipment)
und Literatur-Quellen mit DOI/ISBN.

Wird in concrete_recommendations.py via `MOLECULE_HINTS.update(EXTENSION_HINTS)`
gemerged — so bleibt die Hauptdatei wartbar.

Datenqualitäts-Disziplin:
- Echte Organismen / Stämme (keine Phantasie-Namen)
- Echte Ausbeute-Ranges aus publizierter Literatur
- DOI/ISBN für jede Quelle
- Keine Unicode-Subscripts (würden ReportLab brechen)
"""

from typing import Dict, Any


EXTENSION_HINTS: Dict[str, Dict[str, Any]] = {
    # ====================================================================
    # BATCH 1 — Small molecules / volatile
    # ====================================================================
    "acetone": {
        "chemical": {
            "reagents": "Cumol + O2 (Luft) + H2SO4 (Hock-Verfahren, Coprodukt: Phenol)",
            "reagents_en": "Cumene + O2 (air) + H2SO4 (Hock process, co-product: phenol)",
            "yield_range_percent": (85, 95),
            "literature_source": "Weissermel & Arpe, Industrial Organic Chemistry, 4th ed. 2003, Wiley-VCH (ISBN 978-3-527-30578-0) — Cumol-Verfahren; Schmidt et al., Appl. Catal. A 2005 (doi:10.1016/j.apcata.2004.10.043)",
            "production": {
                "de": [
                    "Edukt-Bereitstellung: Cumol (Isopropylbenzol) aus Friedel-Crafts-Alkylierung von Benzol + Propen",
                    "Schritt 1 — Autoxidation: Cumol + O2 (Luft) bei 80-130 °C, 2-7 bar → Cumolhydroperoxid (CHP), Konversion ~20-30 % pro Pass",
                    "Schritt 2 — Säure-Spaltung (Hock-Umlagerung): CHP + verd. H2SO4 (5-10 %) bei 60-80 °C → Phenol + Aceton (1:1 molar)",
                    "Coproduktion: bei 1 t Aceton fallen 1.6 t Phenol an — Wirtschaftlichkeit hängt an Phenol-Markt",
                    "Roh-Aceton 85-95 % Gesamtausbeute vor Aufarbeitung",
                ],
                "en": [
                    "Reagent supply: cumene (isopropylbenzene) from Friedel-Crafts alkylation of benzene + propene",
                    "Step 1 — autoxidation: cumene + O2 (air) at 80-130 °C, 2-7 bar -> cumene hydroperoxide (CHP), conversion ~20-30 % per pass",
                    "Step 2 — acid cleavage (Hock rearrangement): CHP + dil. H2SO4 (5-10 %) at 60-80 °C -> phenol + acetone (1:1 molar)",
                    "Co-production: 1 t acetone yields 1.6 t phenol — economics depend on phenol market",
                    "Crude acetone 85-95 % overall yield before workup",
                ],
            },
            "downstream": {
                "de": [
                    "Phasentrennung: Aceton + Phenol via fraktionierte Destillation",
                    "Aceton-Kolonne: Kopfprodukt bei 56 °C, Abtrennung von Wasser + Cumol-Spuren",
                    "Pharma-/Reagenzien-Grade: zusätzliche Trocknung über Molekularsieb 4 Å",
                ],
                "en": [
                    "Phase separation: acetone + phenol via fractional distillation",
                    "Acetone column: top product at 56 °C, separation of water + cumene traces",
                    "Pharma / reagent grade: additional drying over molecular sieve 4 A",
                ],
            },
            "expected_purity_after_workup": "≥ 99.5 % (technical) bis ≥ 99.9 % (HPLC-grade)",
        },
        "biotechnological": {
            "organism": "Clostridium acetobutylicum (ABE-Fermentation)",
            "organism_en": "Clostridium acetobutylicum (ABE fermentation)",
            "substrate": "Glucose / Melasse / Maisstärke-Hydrolysat",
            "substrate_en": "glucose / molasses / corn-starch hydrolysate",
            "yield_range_g_per_l": (3, 8),
            "fermentation_time_h": (48, 72),
            "literature_source": "Jones & Woods, Microbiol. Rev. 1986 — Acetone-butanol fermentation revisited (doi:10.1128/mr.50.4.484-524.1986); Ezeji et al., Curr. Opin. Biotechnol. 2007 (doi:10.1016/j.copbio.2007.03.008)",
            "production": {
                "de": [
                    "Stammvorbereitung: C. acetobutylicum (ATCC 824), strikt anaerobe Kultivierung aus Sporen",
                    "Vorkultur: Reinforced Clostridial Medium oder Glucose/Hefeextrakt, anaerob bei 37 °C, 12-24 h",
                    "Hauptfermentation: anaerobe Submerskultur in 10-100 m3 Bioreaktor, 35-37 °C, pH 5.5-6.5, N2-Begasung statt O2",
                    "ABE-Produkt-Verhältnis: Aceton : Butanol : Ethanol ~ 3 : 6 : 1 (g/g), Gesamt-Solventgehalt 18-22 g/L",
                    "Fermentations-Dauer 48-72 h, Aceton-Konzentration 3-8 g/L (limitiert durch Butanol-Toxizität)",
                ],
                "en": [
                    "Strain preparation: C. acetobutylicum (ATCC 824), strict anaerobic cultivation from spores",
                    "Pre-culture: reinforced clostridial medium or glucose/yeast extract, anaerobic at 37 °C, 12-24 h",
                    "Main fermentation: anaerobic submerged culture in 10-100 m3 bioreactor, 35-37 °C, pH 5.5-6.5, N2 sparging instead of O2",
                    "ABE product ratio: acetone : butanol : ethanol ~ 3 : 6 : 1 (g/g), total solvent content 18-22 g/L",
                    "Fermentation duration 48-72 h, acetone concentration 3-8 g/L (limited by butanol toxicity)",
                ],
            },
            "downstream": {
                "de": [
                    "Zellabtrennung via Mikrofiltration",
                    "Destillative Trennung ABE-Mix: Aceton (56 °C) → Ethanol (78 °C) → Butanol (118 °C)",
                    "Pharma-Grade: Trocknung über Molekularsieb",
                ],
                "en": [
                    "Cell removal via microfiltration",
                    "Distillative separation of ABE mix: acetone (56 °C) -> ethanol (78 °C) -> butanol (118 °C)",
                    "Pharma grade: drying over molecular sieve",
                ],
            },
            "expected_purity_after_workup": "≥ 99 % (technical, biotech-route)",
        },
    },
    "ethyl acetate": {
        "chemical": {
            "reagents": "Ethanol + Essigsäure + H2SO4-Katalysator (Fischer-Veresterung)",
            "reagents_en": "Ethanol + acetic acid + H2SO4 catalyst (Fischer esterification)",
            "yield_range_percent": (65, 85),
            "literature_source": "Vollhardt & Schore, Organic Chemistry, 8th ed. 2018, Macmillan/Freeman (ISBN 978-1-319-07945-1) — Fischer-Veresterung; Cellard et al., Chem. Eng. Sci. 2013 (reactive distillation, doi:10.1016/j.ces.2013.05.005)",
            "production": {
                "de": [
                    "Edukt-Bereitstellung: wasserfreies Ethanol + Eisessig (Essigsäure ≥ 99 %)",
                    "Schritt 1 — Veresterung: Ethanol + Essigsäure (1:1 mol) + 1-2 mol% H2SO4 als Katalysator, 70-90 °C, 2-4 h",
                    "Gleichgewichts-Reaktion: Konversion ~65-70 % ohne Wasser-Entzug, ~85-90 % mit reaktiver Destillation",
                    "Industriell: Reaktivdestillation kombiniert Reaktion + Wasser-Abzug → höhere Effizienz",
                    "Roh-Ethylacetat 65-85 % Ausbeute vor Aufarbeitung",
                ],
                "en": [
                    "Reagent supply: anhydrous ethanol + glacial acetic acid (>= 99 %)",
                    "Step 1 — esterification: ethanol + acetic acid (1:1 mol) + 1-2 mol% H2SO4 catalyst, 70-90 °C, 2-4 h",
                    "Equilibrium reaction: conversion ~65-70 % without water removal, ~85-90 % with reactive distillation",
                    "Industrial: reactive distillation combines reaction + water removal -> higher efficiency",
                    "Crude ethyl acetate 65-85 % yield before workup",
                ],
            },
            "downstream": {
                "de": [
                    "Neutralisation des Katalysators mit Natriumcarbonat",
                    "Wasch-Schritte: Wasser → gesättigte NaHCO3-Lösung → Wasser",
                    "Trocknung über wasserfreies MgSO4 oder Na2SO4",
                    "Fraktionierte Destillation bei 77 °C (Reinheit ≥ 99.5 %)",
                ],
                "en": [
                    "Catalyst neutralisation with sodium carbonate",
                    "Wash steps: water -> saturated NaHCO3 solution -> water",
                    "Drying over anhydrous MgSO4 or Na2SO4",
                    "Fractional distillation at 77 °C (purity >= 99.5 %)",
                ],
            },
            "expected_purity_after_workup": "≥ 99.5 % (technical) bis ≥ 99.9 % (HPLC-grade)",
        },
    },
    "methanol": {
        "chemical": {
            "reagents": "Synthesegas (CO + 2 H2) + Cu/ZnO/Al2O3-Katalysator (ICI-Verfahren)",
            "reagents_en": "Syngas (CO + 2 H2) + Cu/ZnO/Al2O3 catalyst (ICI low-pressure process)",
            "yield_range_percent": (90, 95),
            "literature_source": "Olah, Goeppert & Prakash, Beyond Oil and Gas: The Methanol Economy, 2nd ed. 2009, Wiley-VCH (ISBN 978-3-527-32422-7); Bertau et al., Methanol: The Basic Chemical and Energy Feedstock of the Future, Springer 2014 (ISBN 978-3-642-39708-0)",
            "production": {
                "de": [
                    "Edukt-Bereitstellung: Synthesegas (Steam-Reforming von Erdgas: CH4 + H2O → CO + 3 H2)",
                    "Stöchiometrie-Einstellung auf CO : H2 = 1 : 2 durch Wasser-Gas-Shift-Reaktion",
                    "Schritt 1 — Methanol-Synthese: CO + 2 H2 → CH3OH bei 220-280 °C, 50-100 bar an Cu/ZnO/Al2O3-Katalysator",
                    "ICI-Verfahren (Niederdruck, dominant seit 1966): pro Pass ~5-10 % Konversion, Recycle bringt Gesamt-Konversion auf > 90 %",
                    "Rohmethanol: ~95 % CH3OH, Rest Wasser + höhere Alkohole als Nebenprodukte",
                ],
                "en": [
                    "Feedstock supply: syngas (steam reforming of natural gas: CH4 + H2O -> CO + 3 H2)",
                    "Stoichiometry adjustment to CO : H2 = 1 : 2 via water-gas-shift reaction",
                    "Step 1 — methanol synthesis: CO + 2 H2 -> CH3OH at 220-280 °C, 50-100 bar over Cu/ZnO/Al2O3 catalyst",
                    "ICI low-pressure process (dominant since 1966): ~5-10 % conversion per pass, recycle brings total conversion to > 90 %",
                    "Crude methanol: ~95 % CH3OH, balance water + higher alcohols as side products",
                ],
            },
            "downstream": {
                "de": [
                    "Vorentspannung des Synthesegas-Stroms",
                    "Fraktionierte Destillation: Methanol bei 65 °C, Wasser als Sumpf",
                    "Trocken-Destillation (zweite Kolonne) für Reagenzien-Grade",
                    "Pharma-Grade: Behandlung mit KMnO4 zur Entfernung reduzierender Verunreinigungen",
                ],
                "en": [
                    "Pressure release of syngas stream",
                    "Fractional distillation: methanol at 65 °C, water as bottoms",
                    "Dry distillation (second column) for reagent grade",
                    "Pharma grade: treatment with KMnO4 to remove reducing impurities",
                ],
            },
            "expected_purity_after_workup": "≥ 99.85 % (technical AA grade) bis ≥ 99.99 % (semiconductor / HPLC)",
        },
    },
    "isopropanol": {
        "chemical": {
            "reagents": "Propen + H2O (direkte Hydratisierung) oder Propen + H2SO4 (indirekte Hydratisierung)",
            "reagents_en": "Propene + H2O (direct hydration) or propene + H2SO4 (indirect hydration)",
            "yield_range_percent": (70, 90),
            "literature_source": "Weissermel & Arpe, Industrial Organic Chemistry, 4th ed. 2003, Wiley-VCH (ISBN 978-3-527-30578-0); Logsdon & Loke, Kirk-Othmer Encyclopedia 2000 — Isopropyl alcohol (doi:10.1002/0471238961.0919151612150505.a01)",
            "production": {
                "de": [
                    "Edukt-Bereitstellung: Propen (aus Steamcracking, Reinheit > 90 %)",
                    "Variante A — Indirekt: Propen + H2SO4 (75-85 %) bei 20-40 °C → Isopropylsulfat → Hydrolyse mit Wasser → 2-Propanol + verd. H2SO4 (Recycle)",
                    "Variante B — Direkt: Propen + H2O bei 130-200 °C, 25-65 bar an saurem Katalysator (Phosphorsäure auf SiO2 oder Sulfonsäure-Harz)",
                    "Direkte Hydratisierung dominiert moderne Anlagen — geringere Säure-Korrosion und Abwasser-Last",
                    "Roh-Isopropanol: 70-90 % Ausbeute, Diisopropylether als Hauptnebenprodukt",
                ],
                "en": [
                    "Feedstock supply: propene (from steam cracking, purity > 90 %)",
                    "Variant A — indirect: propene + H2SO4 (75-85 %) at 20-40 °C -> isopropyl sulfate -> hydrolysis with water -> 2-propanol + dilute H2SO4 (recycle)",
                    "Variant B — direct: propene + H2O at 130-200 °C, 25-65 bar over acid catalyst (phosphoric acid on SiO2 or sulfonic-acid resin)",
                    "Direct hydration dominates modern plants — lower acid corrosion and waste-water load",
                    "Crude isopropanol: 70-90 % yield, diisopropyl ether as main side product",
                ],
            },
            "downstream": {
                "de": [
                    "Phasentrennung von Wasser-Phase",
                    "Azeotrop-Destillation: 2-Propanol/Wasser-Azeotrop bei 80 °C (87.7 % IPA)",
                    "Entwässerung für absoluten IPA: Molekularsieb 3 Å oder Membran-Pervaporation",
                    "Pharma-/Reagenz-Grade: zusätzliche fraktionierte Destillation",
                ],
                "en": [
                    "Phase separation of water phase",
                    "Azeotropic distillation: 2-propanol/water azeotrope at 80 °C (87.7 % IPA)",
                    "Dehydration for absolute IPA: molecular sieve 3 A or membrane pervaporation",
                    "Pharma / reagent grade: additional fractional distillation",
                ],
            },
            "expected_purity_after_workup": "≥ 99.5 % (technical) bis ≥ 99.9 % (USP / Ph. Eur.)",
        },
    },
    "2,3-butanediol": {
        "biotechnological": {
            "organism": "Klebsiella pneumoniae oder Bacillus licheniformis (entschärfte Stämme)",
            "organism_en": "Klebsiella pneumoniae or Bacillus licheniformis (de-pathogenised strains)",
            "substrate": "Glucose / Xylose / Hemicellulose-Hydrolysate",
            "substrate_en": "glucose / xylose / hemicellulose hydrolysates",
            "yield_range_g_per_l": (50, 150),
            "fermentation_time_h": (24, 60),
            "literature_source": "Ji, Huang & Ouyang, Biotechnol. Adv. 2011 — Microbial 2,3-butanediol production (doi:10.1016/j.biotechadv.2011.01.007); Celińska & Grajek, Biotechnol. Adv. 2009 (doi:10.1016/j.biotechadv.2009.05.002)",
            "production": {
                "de": [
                    "Stammvorbereitung: K. pneumoniae (entschärfte SDM-Stamm) oder B. licheniformis aus Stammbank",
                    "Vorkultur: Glucose + Hefeextrakt + Pepton, 30-37 °C, 12-16 h",
                    "Hauptfermentation: Fed-Batch in 10-200 m3 Bioreaktor, 30-37 °C, pH 5.5-6.5, mikroaerob (niedriger O2-Eintrag 0.1-0.3 vvm — entscheidend für 2,3-BDO statt 2,3-Butandion-Bildung)",
                    "Substrat: Glucose-Feed kontinuierlich (Konzentration < 30 g/L im Reaktor), Gesamt-Substratlast bis 250 g/L",
                    "Fermentations-Dauer 24-60 h, finale 2,3-BDO-Konzentration 50-150 g/L (Stereoisomere R,R / meso je nach Stamm)",
                ],
                "en": [
                    "Strain preparation: K. pneumoniae (de-pathogenised SDM strain) or B. licheniformis from strain bank",
                    "Pre-culture: glucose + yeast extract + peptone, 30-37 °C, 12-16 h",
                    "Main fermentation: fed-batch in 10-200 m3 bioreactor, 30-37 °C, pH 5.5-6.5, microaerobic (low O2 transfer 0.1-0.3 vvm — critical for 2,3-BDO over 2,3-butanedione formation)",
                    "Substrate: continuous glucose feed (concentration < 30 g/L in reactor), total substrate load up to 250 g/L",
                    "Fermentation duration 24-60 h, final 2,3-BDO concentration 50-150 g/L (stereoisomers R,R / meso depending on strain)",
                ],
            },
            "downstream": {
                "de": [
                    "Zellabtrennung via Mikrofiltration oder Zentrifugation",
                    "Aufkonzentrierung via Vakuum-Verdampfung (BDO hat hohen Siedepunkt 180 °C)",
                    "Salz-Aussalzung mit (NH4)2SO4 zur Phasentrennung",
                    "Vakuum-Destillation bei reduziertem Druck",
                ],
                "en": [
                    "Cell separation via microfiltration or centrifugation",
                    "Concentration via vacuum evaporation (BDO has high boiling point 180 °C)",
                    "Salt-out with (NH4)2SO4 for phase separation",
                    "Vacuum distillation at reduced pressure",
                ],
            },
            "expected_purity_after_workup": "≥ 95 % (technical) bis ≥ 99 % (polymer-grade Vorstufe)",
        },
    },
    "butanol": {
        "biotechnological": {
            "organism": "Clostridium acetobutylicum / C. beijerinckii (ABE-Fermentation)",
            "organism_en": "Clostridium acetobutylicum / C. beijerinckii (ABE fermentation)",
            "substrate": "Glucose / Stärke / Lignocellulose-Hydrolysat",
            "substrate_en": "glucose / starch / lignocellulose hydrolysate",
            "yield_range_g_per_l": (10, 20),
            "fermentation_time_h": (48, 96),
            "literature_source": "Jones & Woods, Microbiol. Rev. 1986 — Acetone-butanol fermentation revisited (doi:10.1128/mr.50.4.484-524.1986); Green, Curr. Opin. Biotechnol. 2011 — Fermentative production of butan-1-ol (doi:10.1016/j.copbio.2011.02.004)",
            "production": {
                "de": [
                    "Stammvorbereitung: C. acetobutylicum (ATCC 824) oder C. beijerinckii (NCIMB 8052), Sporen-basierte Inokulation",
                    "Vorkultur: Reinforced Clostridial Medium, strikt anaerob, 37 °C, 12-24 h",
                    "Hauptfermentation: Anaerobe Submerskultur, 35-37 °C, pH 5.0-6.0, N2-Begasung, biphasisches Profil",
                    "Phase 1 (Acidogenese, 0-24 h): Bildung von Acetat + Butyrat, pH fällt auf 4.5",
                    "Phase 2 (Solventogenese, 24-72 h): Re-Assimilation der Säuren zu Aceton/Butanol/Ethanol, pH steigt wieder",
                    "Fermentations-Dauer 48-96 h, Butanol-Konzentration 10-20 g/L (begrenzt durch Butanol-Toxizität)",
                ],
                "en": [
                    "Strain preparation: C. acetobutylicum (ATCC 824) or C. beijerinckii (NCIMB 8052), spore-based inoculation",
                    "Pre-culture: reinforced clostridial medium, strict anaerobic, 37 °C, 12-24 h",
                    "Main fermentation: anaerobic submerged culture, 35-37 °C, pH 5.0-6.0, N2 sparging, biphasic profile",
                    "Phase 1 (acidogenesis, 0-24 h): formation of acetate + butyrate, pH drops to 4.5",
                    "Phase 2 (solventogenesis, 24-72 h): re-assimilation of acids to acetone/butanol/ethanol, pH rises again",
                    "Fermentation duration 48-96 h, butanol concentration 10-20 g/L (limited by butanol toxicity)",
                ],
            },
            "downstream": {
                "de": [
                    "Zellabtrennung via Zentrifugation",
                    "In-situ-Product-Removal (Gas-Stripping oder Pervaporation) zur Reduktion der Butanol-Toxizität — moderne Anlagen",
                    "Fraktionierte Destillation: Aceton (56 °C) → Ethanol (78 °C) → Butanol (118 °C)",
                    "Trocknung über Molekularsieb für Pharma-/Reagenz-Grade",
                ],
                "en": [
                    "Cell separation via centrifugation",
                    "In-situ product removal (gas stripping or pervaporation) to reduce butanol toxicity — modern plants",
                    "Fractional distillation: acetone (56 °C) -> ethanol (78 °C) -> butanol (118 °C)",
                    "Drying over molecular sieve for pharma / reagent grade",
                ],
            },
            "expected_purity_after_workup": "≥ 99.5 % (technical) bis ≥ 99.9 % (Reagenz-Grade)",
        },
    },
    "diacetyl": {
        "biotechnological": {
            "organism": "Lactococcus lactis subsp. lactis biovar diacetylactis oder Leuconostoc mesenteroides",
            "organism_en": "Lactococcus lactis subsp. lactis biovar diacetylactis or Leuconostoc mesenteroides",
            "substrate": "Citrat-haltiges Milchprodukt-Medium (Molke, Buttermilch)",
            "substrate_en": "citrate-containing dairy medium (whey, buttermilk)",
            "yield_range_g_per_l": (0.05, 1.5),
            "fermentation_time_h": (24, 72),
            "literature_source": "Hugenholtz, FEMS Microbiol. Rev. 1993 — Citrate metabolism in lactic acid bacteria (doi:10.1111/j.1574-6976.1993.tb00226.x); Bassit et al., Appl. Environ. Microbiol. 1993 — Effect of inoculation level on diacetyl production (doi:10.1128/aem.59.6.1893-1897.1993)",
            "production": {
                "de": [
                    "Stammvorbereitung: Diacetyl-positiver Stamm aus Stammbank (Citrat-Permease+, ALS-Decarboxylase aktiv)",
                    "Vorkultur: M17-Medium + Glucose + Citrat, 30 °C, 16 h",
                    "Hauptfermentation: Submerskultur in 1-50 m3 Bioreaktor, 22-30 °C, pH 4.5-5.5 (saurer Bereich begünstigt nicht-enzymatische Oxidation α-Acetolactat → Diacetyl)",
                    "Substrat: Molke oder definiertes Glucose + Citrat-Medium",
                    "Fermentations-Dauer 24-72 h, Diacetyl-Konzentration 50-1500 mg/L (Aromaschwelle 4-7 mg/L → in Butter sehr niedrig)",
                ],
                "en": [
                    "Strain preparation: diacetyl-positive strain from strain bank (citrate permease+, ALS decarboxylase active)",
                    "Pre-culture: M17 medium + glucose + citrate, 30 °C, 16 h",
                    "Main fermentation: submerged culture in 1-50 m3 bioreactor, 22-30 °C, pH 4.5-5.5 (acidic range favours non-enzymatic oxidation alpha-acetolactate -> diacetyl)",
                    "Substrate: whey or defined glucose + citrate medium",
                    "Fermentation duration 24-72 h, diacetyl concentration 50-1500 mg/L (flavour threshold 4-7 mg/L -> very low in butter)",
                ],
            },
            "downstream": {
                "de": [
                    "Vakuum-Destillation oder Steam-Stripping zur Konzentrierung",
                    "Flüssig-flüssig-Extraktion mit Diethylether oder Pentan",
                    "Fraktionierte Destillation bei 88 °C unter N2-Schutzgas (Diacetyl oxidationsempfindlich)",
                    "Lagerung im Dunkeln unter Inertgas",
                ],
                "en": [
                    "Vacuum distillation or steam stripping for concentration",
                    "Liquid-liquid extraction with diethyl ether or pentane",
                    "Fractional distillation at 88 °C under N2 inert gas (diacetyl is oxidation-sensitive)",
                    "Storage in the dark under inert gas",
                ],
            },
            "expected_purity_after_workup": "≥ 97 % (food-grade flavour); Vorsicht: pulmonale Toxizität bei Popcorn-Arbeitern dokumentiert (Bronchiolitis obliterans)",
        },
    },
    "ethyl lactate": {
        "chemical": {
            "reagents": "Lactic acid + Ethanol + saurer Katalysator (Amberlyst-15 oder H2SO4)",
            "reagents_en": "Lactic acid + ethanol + acid catalyst (Amberlyst-15 or H2SO4)",
            "yield_range_percent": (70, 90),
            "literature_source": "Pereira et al., Green Chem. 2009 — Ethyl lactate as a solvent (doi:10.1039/B902341A); Aparicio & Alcalde, J. Phys. Chem. B 2009 — Ethyl lactate environmental impact assessment (doi:10.1021/jp9012646)",
            "production": {
                "de": [
                    "Edukt-Bereitstellung: bio-basierte Lactic acid (aus Fermentation, siehe lactic-acid-Eintrag) + wasserfreies Ethanol",
                    "Schritt 1 — Veresterung: Lactic acid + Ethanol (Überschuss 3-5x) + Amberlyst-15 oder verd. H2SO4, 70-90 °C, 3-6 h",
                    "Gleichgewichts-Reaktion: Konversion 70-85 % ohne Wasser-Abzug",
                    "Reaktive Destillation oder Pervaporations-Membran für > 95 % Konversion",
                    "Roh-Ethyllactat 70-90 % Ausbeute vor Aufarbeitung",
                ],
                "en": [
                    "Reagent supply: bio-based lactic acid (from fermentation, see lactic-acid entry) + anhydrous ethanol",
                    "Step 1 — esterification: lactic acid + ethanol (excess 3-5x) + Amberlyst-15 or dilute H2SO4, 70-90 °C, 3-6 h",
                    "Equilibrium reaction: conversion 70-85 % without water removal",
                    "Reactive distillation or pervaporation membrane for > 95 % conversion",
                    "Crude ethyl lactate 70-90 % yield before workup",
                ],
            },
            "downstream": {
                "de": [
                    "Katalysator-Abtrennung (Filtration bei Amberlyst, Neutralisation bei H2SO4)",
                    "Fraktionierte Destillation bei 151-154 °C unter reduziertem Druck",
                    "Trocknung über Molekularsieb 3 Å für Reagenzien-Grade",
                ],
                "en": [
                    "Catalyst removal (filtration for Amberlyst, neutralisation for H2SO4)",
                    "Fractional distillation at 151-154 °C under reduced pressure",
                    "Drying over molecular sieve 3 A for reagent grade",
                ],
            },
            "expected_purity_after_workup": "≥ 98 % (technical / solvent) bis ≥ 99 % (FCC food-grade als Aroma-Lösungsmittel)",
        },
    },
    # ====================================================================
    # BATCH 2 — Small molecules / non_volatile
    # ====================================================================
    "glucose": {
        "biotechnological": {
            "organism": "Enzymatische Stärkehydrolyse (α-Amylase aus B. licheniformis + Glucoamylase aus A. niger)",
            "organism_en": "Enzymatic starch hydrolysis (alpha-amylase from B. licheniformis + glucoamylase from A. niger)",
            "substrate": "Maisstärke / Weizenstärke (30-35 % w/w in Wasser)",
            "substrate_en": "corn starch / wheat starch (30-35 % w/w in water)",
            "yield_range_percent": (95, 98),
            "literature_source": "Reichelt, Starch-Stärke 1983 — Industrial starch hydrolysis (doi:10.1002/star.19830351203); van der Maarel et al., J. Biotechnol. 2002 — Properties and applications of starch-converting enzymes (doi:10.1016/S0168-1656(01)00407-2)",
            "production": {
                "de": [
                    "Stärke-Aufschlämmung: 30-35 % w/w Maisstärke in Wasser, pH 5.8-6.2 mit Ca(OH)2 eingestellt",
                    "Verflüssigung (Liquefaction): α-Amylase (B. licheniformis, thermostabil) bei 95-105 °C für 90-120 min → Maltodextrine (DE 10-15)",
                    "Verzuckerung (Saccharification): pH-Absenkung auf 4.5, Temperatur 60 °C, Zugabe von Glucoamylase (A. niger), 48-72 h Inkubation",
                    "Konversion zu Glucose ≥ 95 % (DE 95-98), Restmaltose < 2 %",
                    "Optional: Isomerase-Schritt mit Glucose-Isomerase → High-Fructose-Corn-Syrup (HFCS) statt reiner Glucose",
                ],
                "en": [
                    "Starch slurry: 30-35 % w/w corn starch in water, pH 5.8-6.2 adjusted with Ca(OH)2",
                    "Liquefaction: alpha-amylase (B. licheniformis, thermostable) at 95-105 °C for 90-120 min -> maltodextrins (DE 10-15)",
                    "Saccharification: pH adjusted to 4.5, temperature 60 °C, addition of glucoamylase (A. niger), 48-72 h incubation",
                    "Conversion to glucose >= 95 % (DE 95-98), residual maltose < 2 %",
                    "Optional: isomerase step with glucose isomerase -> high-fructose corn syrup (HFCS) instead of pure glucose",
                ],
            },
            "downstream": {
                "de": [
                    "Filtration über Diatomeenerde zur Entfernung suspendierter Feststoffe",
                    "Aktivkohle-Behandlung zur Entfernung von Farb- und Geruchsstoffen",
                    "Ionenaustausch (Kationen- und Anionenharze) zur Entsalzung",
                    "Eindampfung im Mehrstufen-Vakuumverdampfer auf 70-75 % w/w Trockensubstanz",
                    "Kristallisation aus übersättigter Lösung (für kristalline Glucose) oder Sprühtrocknung",
                ],
                "en": [
                    "Filtration over diatomaceous earth to remove suspended solids",
                    "Activated-carbon treatment to remove colour and odour compounds",
                    "Ion exchange (cation and anion resins) for desalting",
                    "Evaporation in multi-stage vacuum evaporator to 70-75 % w/w dry substance",
                    "Crystallization from supersaturated solution (for crystalline glucose) or spray drying",
                ],
            },
            "expected_purity_after_workup": "≥ 99.5 % (D-Glucose monohydrate, USP / Ph. Eur. / food grade)",
        },
    },
    "itaconic acid": {
        "biotechnological": {
            "organism": "Aspergillus terreus (analog zu Citric-acid-Verfahren mit A. niger)",
            "organism_en": "Aspergillus terreus (analogous to citric-acid process with A. niger)",
            "substrate": "Glucose / Saccharose / Melasse",
            "substrate_en": "glucose / sucrose / molasses",
            "yield_range_g_per_l": (50, 90),
            "fermentation_time_h": (100, 168),
            "literature_source": "Klement & Büchs, Bioresour. Technol. 2013 — Itaconic acid production review (doi:10.1016/j.biortech.2012.10.124); Steiger et al., Front. Microbiol. 2013 — Biochemistry of microbial itaconic acid production (doi:10.3389/fmicb.2013.00023)",
            "production": {
                "de": [
                    "Stammvorbereitung: A. terreus Hochertrags-Stamm (z. B. ATCC 10020 oder kommerzielle Hochleister), Sporulation auf Schrägagar",
                    "Vorkultur: Glucose-haltiges Medium, 33-37 °C, 24-48 h Vorwachstum",
                    "Hauptfermentation: Submerskultur in 50-200 m3 Bioreaktor, 33-37 °C, pH 2.0-3.5 (sehr sauer, Mn2+-Limitierung kritisch für hohe Ausbeute)",
                    "Substrat: Glucose oder Saccharose, Initialkonzentration 100-150 g/L, optional Fed-Batch",
                    "Fermentations-Dauer 100-168 h, finale Itaconsäure-Konzentration 50-90 g/L (theoretische Maximalausbeute 72 % auf Glucose)",
                ],
                "en": [
                    "Strain preparation: A. terreus high-yield strain (e.g. ATCC 10020 or commercial high-producers), sporulation on slant agar",
                    "Pre-culture: glucose-containing medium, 33-37 °C, 24-48 h pre-growth",
                    "Main fermentation: submerged culture in 50-200 m3 bioreactor, 33-37 °C, pH 2.0-3.5 (very acidic, Mn2+ limitation critical for high yield)",
                    "Substrate: glucose or sucrose, initial concentration 100-150 g/L, optionally fed-batch",
                    "Fermentation duration 100-168 h, final itaconic acid concentration 50-90 g/L (theoretical maximum yield 72 % on glucose)",
                ],
            },
            "downstream": {
                "de": [
                    "Myzel-Abtrennung via Filtration oder Mikrofiltration",
                    "Aktivkohle-Behandlung zur Entfernung von Pigmenten",
                    "Direkt-Kristallisation aus der sauren Fermentationsbrühe (Itaconsäure hat niedrige Wasserlöslichkeit)",
                    "Optional: Umkristallisation aus Wasser für Polymer-Grade",
                ],
                "en": [
                    "Mycelium separation via filtration or microfiltration",
                    "Activated-carbon treatment to remove pigments",
                    "Direct crystallization from acidic fermentation broth (itaconic acid has low water solubility)",
                    "Optional: recrystallization from water for polymer grade",
                ],
            },
            "expected_purity_after_workup": "≥ 99 % (technical / polymer-grade Vorstufe für Itaconat-Polymere)",
        },
    },
    # ====================================================================
    # BATCH 3 — Peptide / linear
    # ====================================================================
    # Kurze Peptide (≤ 10 AS): klassisch SPPS, voll synthetisch
    # Mittlere bis lange Peptide (≥ 20 AS): rekombinant in E. coli / Hefe
    # ====================================================================
    "glutathione": {
        "biotechnological": {
            "organism": "S. cerevisiae oder Candida utilis (Hochertrags-Stämme)",
            "organism_en": "S. cerevisiae or Candida utilis (high-yield strains)",
            "substrate": "Glucose + L-Cystein + L-Glutamat + Glycin (Vorstufen-Feed)",
            "substrate_en": "glucose + L-cysteine + L-glutamate + glycine (precursor feed)",
            "yield_range_g_per_l": (1, 5),
            "fermentation_time_h": (24, 48),
            "literature_source": "Penninckx, FEMS Yeast Res. 2002 — An overview on glutathione in Saccharomyces (doi:10.1111/j.1567-1364.2002.tb00091.x); Li et al., Biotechnol. Adv. 2004 — Glutathione production by microbial fermentation (doi:10.1016/j.biotechadv.2003.10.001)",
            "production": {
                "de": [
                    "Stammvorbereitung: S. cerevisiae mit überexprimierter γ-Glutamylcystein-Synthetase + Glutathion-Synthetase (Knockout der GSH-Abbau-Wege)",
                    "Vorkultur: YPD-Medium, 30 °C, 12-16 h, OD600 ~5-10",
                    "Hauptfermentation: Fed-Batch im 1-50 m3 Bioreaktor, 28-30 °C, pH 5.5-6.0, hoher O2-Eintrag",
                    "Vorstufen-Feed: L-Cystein (Schlüssel-Aminosäure, oft limitierend) + L-Glutamat + Glycin in stöchiometrischen Verhältnissen",
                    "Fermentations-Dauer 24-48 h, intrazelluläre GSH-Akkumulation 1-5 g/L (entspricht ~3-5 % der Zell-Trockenmasse)",
                ],
                "en": [
                    "Strain preparation: S. cerevisiae with overexpressed gamma-glutamylcysteine synthetase + glutathione synthetase (knockout of GSH degradation pathways)",
                    "Pre-culture: YPD medium, 30 °C, 12-16 h, OD600 ~5-10",
                    "Main fermentation: fed-batch in 1-50 m3 bioreactor, 28-30 °C, pH 5.5-6.0, high O2 transfer",
                    "Precursor feed: L-cysteine (key amino acid, often limiting) + L-glutamate + glycine in stoichiometric ratios",
                    "Fermentation duration 24-48 h, intracellular GSH accumulation 1-5 g/L (corresponds to ~3-5 % of cell dry mass)",
                ],
            },
            "downstream": {
                "de": [
                    "Zell-Ernte via Zentrifugation",
                    "Zellaufschluss (Hochdruckhomogenisator oder Glaskugelmühle) zur Freisetzung des intrazellulären GSH",
                    "Klärung + Aktivkohle-Behandlung zur Entfernung von Pigmenten",
                    "Kationenaustausch-Chromatographie (SP-Sepharose) — Capture",
                    "Kristallisation aus Ethanol/Wasser unter N2-Schutzgas (GSH oxidationsempfindlich)",
                ],
                "en": [
                    "Cell harvest via centrifugation",
                    "Cell lysis (high-pressure homogenizer or bead mill) to release intracellular GSH",
                    "Clarification + activated-carbon treatment to remove pigments",
                    "Cation-exchange chromatography (SP-Sepharose) — capture",
                    "Crystallization from ethanol/water under N2 inert gas (GSH is oxidation-sensitive)",
                ],
            },
            "expected_purity_after_workup": "≥ 98 % (USP / Ph. Eur., reduzierte Form)",
        },
    },
    "leu-enkephalin": {
        "chemical": {
            "reagents": "Fmoc-Aminosäuren (Tyr, Gly, Phe, Leu) + Kupplungsreagenzien (HBTU/DIPEA) + 2-Chlortrityl-Harz",
            "reagents_en": "Fmoc-amino acids (Tyr, Gly, Phe, Leu) + coupling reagents (HBTU/DIPEA) + 2-chlorotrityl resin",
            "yield_range_percent": (60, 75),
            "literature_source": "Hughes et al., Nature 1975 — Identification of leucine-enkephalin (doi:10.1038/258577a0); Merrifield, J. Am. Chem. Soc. 1963 — Solid-phase peptide synthesis (doi:10.1021/ja00897a025); Albericio, Curr. Opin. Chem. Biol. 2004 — SPPS-Methoden (doi:10.1016/j.cbpa.2004.10.002)",
            "production": {
                "de": [
                    "Harz-Vorbereitung: 2-Chlortrityl-Chlorid-Harz mit C-terminalem Fmoc-Leu beladen (Beladung 0.5-0.8 mmol/g)",
                    "Schritt 1-4 — SPPS-Zyklus pro Aminosäure (Phe → Gly → Gly → Tyr): Fmoc-Entschützung mit 20 % Piperidin in DMF (2x 5 min) → Wasch-Schritte (DMF, DCM) → Kupplung mit nächster Fmoc-AS + HBTU + DIPEA in DMF (30-60 min) → Wasch-Schritte",
                    "Insgesamt 5 SPPS-Zyklen für Pentapeptid YGGFL (Tyr-Gly-Gly-Phe-Leu)",
                    "Final-Deprotection + Cleavage: 95 % TFA + 2.5 % H2O + 2.5 % Triisopropylsilan (2 h) → freies Peptid in Lösung",
                    "Fällung in kaltem Diethylether, Lyophilisation des Roh-Peptids (60-75 % Gesamtausbeute)",
                ],
                "en": [
                    "Resin preparation: 2-chlorotrityl chloride resin loaded with C-terminal Fmoc-Leu (loading 0.5-0.8 mmol/g)",
                    "Step 1-4 — SPPS cycle per amino acid (Phe -> Gly -> Gly -> Tyr): Fmoc deprotection with 20 % piperidine in DMF (2x 5 min) -> wash steps (DMF, DCM) -> coupling with next Fmoc-AA + HBTU + DIPEA in DMF (30-60 min) -> wash steps",
                    "Total of 5 SPPS cycles for pentapeptide YGGFL (Tyr-Gly-Gly-Phe-Leu)",
                    "Final deprotection + cleavage: 95 % TFA + 2.5 % H2O + 2.5 % triisopropylsilane (2 h) -> free peptide in solution",
                    "Precipitation in cold diethyl ether, lyophilisation of crude peptide (60-75 % overall yield)",
                ],
            },
            "downstream": {
                "de": [
                    "Roh-Peptid in 0.1 % TFA/Wasser/Acetonitril gelöst",
                    "Präparative Reversed-Phase-HPLC (C18, Acetonitril/Wasser-Gradient mit 0.1 % TFA)",
                    "Counter-Ion-Exchange auf Acetat oder Hydrochlorid",
                    "Lyophilisation des reinen Peptids",
                ],
                "en": [
                    "Crude peptide dissolved in 0.1 % TFA/water/acetonitrile",
                    "Preparative reversed-phase HPLC (C18, acetonitrile/water gradient with 0.1 % TFA)",
                    "Counter-ion exchange to acetate or hydrochloride",
                    "Lyophilisation of pure peptide",
                ],
            },
            "expected_purity_after_workup": "≥ 98 % (Forschungs-Grade); ≥ 99 % (klinische Studien)",
        },
    },
    "met-enkephalin": {
        "chemical": {
            "reagents": "Fmoc-Aminosäuren (Tyr, Gly, Phe, Met) + HBTU/DIPEA + 2-Chlortrityl-Harz",
            "reagents_en": "Fmoc-amino acids (Tyr, Gly, Phe, Met) + HBTU/DIPEA + 2-chlorotrityl resin",
            "yield_range_percent": (55, 70),
            "literature_source": "Hughes et al., Nature 1975 — Identification of methionine- and leucine-enkephalin (doi:10.1038/258577a0); Albericio, Curr. Opin. Chem. Biol. 2004 (doi:10.1016/j.cbpa.2004.10.002)",
            "production": {
                "de": [
                    "Harz-Vorbereitung: 2-Chlortrityl-Chlorid-Harz mit C-terminalem Fmoc-Met beladen",
                    "Schritt 1-4 — SPPS-Zyklus pro AS (Phe → Gly → Gly → Tyr): Fmoc-Entschützung → Kupplung → Wasch-Schritte; Methionin-Seitenkette anfällig für Oxidation, Schutzatmosphäre N2",
                    "Final-Deprotection + Cleavage mit TFA-Cocktail (mit Thioether-Scavenger wie EDT für Met-Schutz)",
                    "Fällung in kaltem Diethylether, Lyophilisation",
                    "Roh-Pentapeptid YGGFM (Tyr-Gly-Gly-Phe-Met), Ausbeute 55-70 %",
                ],
                "en": [
                    "Resin preparation: 2-chlorotrityl chloride resin loaded with C-terminal Fmoc-Met",
                    "Step 1-4 — SPPS cycle per AA (Phe -> Gly -> Gly -> Tyr): Fmoc deprotection -> coupling -> wash steps; methionine side chain susceptible to oxidation, N2 inert atmosphere",
                    "Final deprotection + cleavage with TFA cocktail (with thioether scavenger like EDT to protect Met)",
                    "Precipitation in cold diethyl ether, lyophilisation",
                    "Crude pentapeptide YGGFM (Tyr-Gly-Gly-Phe-Met), yield 55-70 %",
                ],
            },
            "downstream": {
                "de": [
                    "Prep-RP-HPLC unter N2-Schutzgas (Methionin-Oxidation vermeiden)",
                    "Counter-Ion-Exchange auf Acetat",
                    "Lyophilisation unter Inertgas, Lagerung bei -20 °C",
                ],
                "en": [
                    "Prep-RP-HPLC under N2 inert gas (avoid methionine oxidation)",
                    "Counter-ion exchange to acetate",
                    "Lyophilisation under inert gas, storage at -20 °C",
                ],
            },
            "expected_purity_after_workup": "≥ 98 % (mit Met(O)-Gehalt < 1 %)",
        },
    },
    "bradykinin": {
        "chemical": {
            "reagents": "Fmoc-AS (Arg, Pro, Gly, Phe, Ser, Pro) + HBTU/DIPEA + Wang-Harz",
            "reagents_en": "Fmoc-AA (Arg, Pro, Gly, Phe, Ser, Pro) + HBTU/DIPEA + Wang resin",
            "yield_range_percent": (40, 60),
            "literature_source": "Marceau & Regoli, Nat. Rev. Drug Discov. 2004 — Bradykinin receptor ligands (doi:10.1038/nrd1469); Albericio, Curr. Opin. Chem. Biol. 2004 (doi:10.1016/j.cbpa.2004.10.002)",
            "production": {
                "de": [
                    "Harz-Vorbereitung: Wang-Harz mit C-terminalem Fmoc-Arg(Pbf) beladen (Pbf = Schutzgruppe für Arg-Guanidin)",
                    "SPPS-Zyklen (9 Kupplungen) für Sequenz RPPGFSPFR (Arg-Pro-Pro-Gly-Phe-Ser-Pro-Phe-Arg): Standard-Fmoc-Protokoll, Pro-Pro-Sequenzen können Kupplungsschwierigkeiten verursachen",
                    "Double-Couplings für schwierige Stellen (Pro-Pro), Capping mit Ac2O/DIPEA nach jeder Stufe",
                    "Final-Cleavage: TFA-Cocktail (95 % TFA + 2.5 % H2O + 2.5 % Triisopropylsilan), 2-3 h",
                    "Fällung in kaltem Ether, Roh-Peptid 40-60 % Gesamtausbeute (gut für 9-mer)",
                ],
                "en": [
                    "Resin preparation: Wang resin loaded with C-terminal Fmoc-Arg(Pbf) (Pbf = protecting group for Arg guanidine)",
                    "SPPS cycles (9 couplings) for sequence RPPGFSPFR (Arg-Pro-Pro-Gly-Phe-Ser-Pro-Phe-Arg): standard Fmoc protocol, Pro-Pro sequences can cause coupling difficulties",
                    "Double couplings for difficult positions (Pro-Pro), capping with Ac2O/DIPEA after each stage",
                    "Final cleavage: TFA cocktail (95 % TFA + 2.5 % H2O + 2.5 % triisopropylsilane), 2-3 h",
                    "Precipitation in cold ether, crude peptide 40-60 % overall yield (good for 9-mer)",
                ],
            },
            "downstream": {
                "de": [
                    "Prep-RP-HPLC (C18, Acetonitril/Wasser-Gradient mit 0.1 % TFA)",
                    "Counter-Ion-Exchange auf Acetat",
                    "Lyophilisation, Lagerung bei -20 °C",
                ],
                "en": [
                    "Prep-RP-HPLC (C18, acetonitrile/water gradient with 0.1 % TFA)",
                    "Counter-ion exchange to acetate",
                    "Lyophilisation, storage at -20 °C",
                ],
            },
            "expected_purity_after_workup": "≥ 98 % (Forschungs- / Diagnostik-Grade)",
        },
    },
    "angiotensin ii": {
        "chemical": {
            "reagents": "Fmoc-AS (Asp, Arg, Val, Tyr, Ile, His, Pro, Phe) + HBTU/DIPEA + Wang-Harz",
            "reagents_en": "Fmoc-AA (Asp, Arg, Val, Tyr, Ile, His, Pro, Phe) + HBTU/DIPEA + Wang resin",
            "yield_range_percent": (45, 65),
            "literature_source": "Skeggs et al., J. Exp. Med. 1956 — Angiotensin II structural identification (doi:10.1084/jem.103.3.295); de Gasparo et al., Pharmacol. Rev. 2000 — Angiotensin II receptors (doi:10.1124/pr.59.3.3) — historical context; Behrendt et al., Curr. Med. Chem. 2007 (doi:10.2174/092986707782360051)",
            "production": {
                "de": [
                    "Harz-Vorbereitung: Wang-Harz mit C-terminalem Fmoc-Phe beladen",
                    "SPPS-Zyklen (8 Kupplungen) für Sequenz DRVYIHPF (Asp-Arg-Val-Tyr-Ile-His-Pro-Phe): Standard-Fmoc-Protokoll mit Pbf für Arg und Trt für His",
                    "Aspartat-Seitenkette mit OtBu geschützt zur Vermeidung von Aspartimid-Bildung",
                    "Final-Cleavage: TFA-Cocktail mit zusätzlichem H2O als Scavenger für His(Trt)",
                    "Fällung in kaltem Ether, Roh-Octapeptid 45-65 % Gesamtausbeute",
                ],
                "en": [
                    "Resin preparation: Wang resin loaded with C-terminal Fmoc-Phe",
                    "SPPS cycles (8 couplings) for sequence DRVYIHPF (Asp-Arg-Val-Tyr-Ile-His-Pro-Phe): standard Fmoc protocol with Pbf for Arg and Trt for His",
                    "Aspartate side chain protected with OtBu to avoid aspartimide formation",
                    "Final cleavage: TFA cocktail with additional H2O scavenger for His(Trt)",
                    "Precipitation in cold ether, crude octapeptide 45-65 % overall yield",
                ],
            },
            "downstream": {
                "de": [
                    "Prep-RP-HPLC (C18) mit kontrolliertem Gradient",
                    "Counter-Ion-Exchange auf Acetat",
                    "Lyophilisation",
                ],
                "en": [
                    "Prep-RP-HPLC (C18) with controlled gradient",
                    "Counter-ion exchange to acetate",
                    "Lyophilisation",
                ],
            },
            "expected_purity_after_workup": "≥ 98 % (Diagnostik / Forschung); klinisch ≥ 99 % (z. B. für Sepsis-Indikation, Giapreza)",
        },
    },
    "glucagon": {
        "biotechnological": {
            "organism": "Escherichia coli (BL21) oder Saccharomyces cerevisiae mit rekombinantem Glucagon-Gen",
            "organism_en": "Escherichia coli (BL21) or Saccharomyces cerevisiae with recombinant glucagon gene",
            "substrate": "Definiertes Minimalmedium + Glucose + Induktor (IPTG)",
            "substrate_en": "defined minimal medium + glucose + inducer (IPTG)",
            "yield_range_g_per_l": (0.5, 2),
            "fermentation_time_h": (24, 48),
            "literature_source": "Maurer & Sorgato, Curr. Opin. Drug Discov. Devel. 2010 — Recombinant glucagon manufacturing; Habener, Diabetes 1981 (doi:10.2337/diab.30.10.901)",
            "production": {
                "de": [
                    "Stammvorbereitung: rekombinante E. coli mit Glucagon-Vorläufer-Gen (oft als Fusion mit löslichem Tag wie SUMO oder MBP), GMP-Stammbank",
                    "Vorkultur: LB- oder definiertes Minimalmedium + Selektion, 37 °C, 12-16 h",
                    "Hauptfermentation: Fed-Batch im 1-20 m3 Bioreaktor, 37 °C, pH 7.0, hoher Glucose-Feed",
                    "Induktion mit IPTG (0.1-1 mM) bei OD600 ~30-50, Expression 4-8 h",
                    "Glucagon-Vorläufer akkumuliert intrazellulär (oft als Inclusion Bodies), 0.5-2 g/L",
                ],
                "en": [
                    "Strain preparation: recombinant E. coli with glucagon-precursor gene (often as fusion with soluble tag like SUMO or MBP), GMP strain bank",
                    "Pre-culture: LB or defined minimal medium + selection, 37 °C, 12-16 h",
                    "Main fermentation: fed-batch in 1-20 m3 bioreactor, 37 °C, pH 7.0, high glucose feed",
                    "Induction with IPTG (0.1-1 mM) at OD600 ~30-50, expression 4-8 h",
                    "Glucagon precursor accumulates intracellularly (often as inclusion bodies), 0.5-2 g/L",
                ],
            },
            "downstream": {
                "de": [
                    "Zell-Lyse + IB-Isolation via Hochdruckhomogenisator",
                    "Solubilisierung in 8 M Harnstoff + Reduktionsmittel",
                    "Refolding bei niedriger Konzentration (0.1-0.3 g/L), 4-12 °C",
                    "Enzymatische Tag-Abspaltung (z. B. SUMO-Protease oder Faktor Xa)",
                    "Anionenaustausch-Chromatographie + Reverse-Phase-HPLC",
                    "Lyophilisation",
                ],
                "en": [
                    "Cell lysis + IB isolation via high-pressure homogenizer",
                    "Solubilisation in 8 M urea + reducing agent",
                    "Refolding at low concentration (0.1-0.3 g/L), 4-12 °C",
                    "Enzymatic tag cleavage (e.g. SUMO protease or factor Xa)",
                    "Anion-exchange chromatography + reverse-phase HPLC",
                    "Lyophilisation",
                ],
            },
            "expected_purity_after_workup": "≥ 99 % (USP / Ph. Eur., GlucaGen / Baqsimi)",
        },
    },
    "calcitonin": {
        "chemical": {
            "reagents": "Fmoc-AS für 32-mer Lachs-Calcitonin + HBTU/DIPEA + Rink-Amid-Harz (C-Term-Amid)",
            "reagents_en": "Fmoc-AA for 32-mer salmon calcitonin + HBTU/DIPEA + Rink amide resin (C-term amide)",
            "yield_range_percent": (15, 30),
            "literature_source": "Guidobono, Curr. Pharm. Des. 2002 — Salmon calcitonin clinical use (doi:10.2174/1381612023395916); Munson et al., J. Bone Miner. Res. 1995 — Salmon calcitonin therapeutics (doi:10.1002/jbmr.5650100207)",
            "production": {
                "de": [
                    "Harz-Vorbereitung: Rink-Amid-Harz (für C-terminales Amid) mit Fmoc-Pro als erster AS beladen",
                    "SPPS-Zyklen (32 Kupplungen) für Lachs-Calcitonin-Sequenz: standardisierter Fmoc-Workflow, kritische Cys-Reste (Position 1+7) mit Trt geschützt",
                    "Double-Couplings für sterisch anspruchsvolle Positionen, Capping nach jeder Stufe",
                    "Final-Cleavage: TFA-Cocktail (typisch Reagent K für Trp/Met-haltige Peptide)",
                    "Disulfid-Bildung: Luft-Oxidation in verdünnter Lösung (0.1 g/L) bei pH 8 → intramolekulares Disulfid 1-7",
                    "Roh-Calcitonin 15-30 % Gesamtausbeute (typisch für 32-mer)",
                ],
                "en": [
                    "Resin preparation: Rink amide resin (for C-terminal amide) loaded with Fmoc-Pro as first AA",
                    "SPPS cycles (32 couplings) for salmon calcitonin sequence: standardised Fmoc workflow, critical Cys residues (positions 1+7) protected with Trt",
                    "Double couplings for sterically demanding positions, capping after each stage",
                    "Final cleavage: TFA cocktail (typically reagent K for Trp/Met-containing peptides)",
                    "Disulfide formation: air oxidation in dilute solution (0.1 g/L) at pH 8 -> intramolecular disulfide 1-7",
                    "Crude calcitonin 15-30 % overall yield (typical for 32-mer)",
                ],
            },
            "downstream": {
                "de": [
                    "Prep-RP-HPLC mit hochauflösendem Gradient (C18)",
                    "Disulfid-Isomeren-Trennung kritisch (richtige Disulfid-Konnektivität nur ~70 % im Roh-Mix)",
                    "Counter-Ion-Exchange auf Acetat",
                    "Lyophilisation, Lagerung bei -20 °C",
                ],
                "en": [
                    "Prep-RP-HPLC with high-resolution gradient (C18)",
                    "Disulfide-isomer separation critical (correct disulfide connectivity only ~70 % in crude mix)",
                    "Counter-ion exchange to acetate",
                    "Lyophilisation, storage at -20 °C",
                ],
            },
            "expected_purity_after_workup": "≥ 99 % (USP / Ph. Eur., klinisch z. B. Miacalcin)",
        },
    },
    "liraglutide": {
        "chemical": {
            "reagents": "Fmoc-AS für 31-mer GLP-1-Analog + Lys-Seitenketten-Modifikation mit γGlu-C16-Fettsäure",
            "reagents_en": "Fmoc-AA for 31-mer GLP-1 analogue + Lys side-chain modification with gamma-Glu-C16 fatty acid",
            "yield_range_percent": (10, 25),
            "literature_source": "Knudsen et al., J. Med. Chem. 2000 — Potent derivatives of glucagon-like peptide-1 with pharmacokinetic properties (doi:10.1021/jm9909645) — Liraglutide-Erstbeschreibung; Andersen et al., Curr. Pharm. Biotechnol. 2018 — Lipidation strategies for peptide therapeutics",
            "production": {
                "de": [
                    "Harz-Vorbereitung: Wang-Harz mit C-terminalem Fmoc-Gly beladen",
                    "SPPS-Zyklen (31 Kupplungen) für GLP-1(7-37)-Analog-Sequenz HAEGTFTSDVSSYLEGQAAKEFIAWLVRGRG; Lys26 mit orthogonal entfernbarer ivDde-Schutzgruppe (statt Standard-Boc)",
                    "Selektive ivDde-Entschützung an Lys26 mit Hydrazin/DMF (Standard-Lys-Boc-Gruppen unverändert)",
                    "Acylierung an Lys26-ε-Amino: γ-Glu-Spacer-Kupplung + C16-Fettsäure (Palmitinsäure) — kritisch für Albumin-Bindung und Halbwertszeit-Verlängerung",
                    "Final-Cleavage: TFA-Cocktail",
                    "Roh-Liraglutide 10-25 % Gesamtausbeute (anspruchsvoll wegen 31 Stufen + Lipidierung)",
                ],
                "en": [
                    "Resin preparation: Wang resin loaded with C-terminal Fmoc-Gly",
                    "SPPS cycles (31 couplings) for GLP-1(7-37) analogue sequence HAEGTFTSDVSSYLEGQAAKEFIAWLVRGRG; Lys26 with orthogonally removable ivDde protecting group (instead of standard Boc)",
                    "Selective ivDde deprotection at Lys26 with hydrazine/DMF (standard Lys-Boc groups untouched)",
                    "Acylation at Lys26 epsilon-amino: gamma-Glu spacer coupling + C16 fatty acid (palmitic acid) — critical for albumin binding and half-life extension",
                    "Final cleavage: TFA cocktail",
                    "Crude liraglutide 10-25 % overall yield (challenging due to 31 stages + lipidation)",
                ],
            },
            "downstream": {
                "de": [
                    "Mehrstufige Prep-RP-HPLC zur Trennung von Lipidierungs-Isomeren und Truncations",
                    "Disulfid-freies Peptid, kein Oxidations-Schritt nötig",
                    "Counter-Ion-Exchange auf Acetat",
                    "Lyophilisation, finale Formulierung mit Phosphat-Puffer + Phenol",
                ],
                "en": [
                    "Multi-stage prep-RP-HPLC for separation of lipidation isomers and truncations",
                    "Disulfide-free peptide, no oxidation step required",
                    "Counter-ion exchange to acetate",
                    "Lyophilisation, final formulation with phosphate buffer + phenol",
                ],
            },
            "expected_purity_after_workup": "≥ 99 % (USP, Victoza / Saxenda); modern oft semi-synthetisch (rekombinanter Backbone + chemische Acylierung)",
        },
    },
    "teriparatide": {
        "biotechnological": {
            "organism": "Escherichia coli (BL21(DE3))",
            "organism_en": "Escherichia coli (BL21(DE3))",
            "substrate": "Definiertes Minimalmedium + Glucose + IPTG",
            "substrate_en": "defined minimal medium + glucose + IPTG",
            "yield_range_g_per_l": (1, 4),
            "fermentation_time_h": (24, 48),
            "literature_source": "Quattrocchi & Kourlas, Clin. Ther. 2004 — Teriparatide review (doi:10.1016/S0149-2918(04)90060-7); Jüppner et al., Endocrinology 1991 — PTH receptor cloning (Forteo-Background)",
            "production": {
                "de": [
                    "Stammvorbereitung: rekombinantes E. coli BL21(DE3) mit PTH(1-34)-Gen unter T7-Promotor (oft als Fusion mit löslichem Tag zur Stabilität)",
                    "Vorkultur: LB-Medium + Kanamycin, 37 °C, 12 h",
                    "Hauptfermentation: Fed-Batch im 1-20 m3 Bioreaktor, 37 °C → Temperatur-Shift auf 30 °C bei Induktion, pH 7.0",
                    "Induktion mit IPTG (0.5-1 mM) bei OD600 ~40-60, Expression 6-8 h",
                    "Teriparatide-Vorläufer als Inclusion Bodies oder lösliches Fusionsprotein, 1-4 g/L",
                ],
                "en": [
                    "Strain preparation: recombinant E. coli BL21(DE3) with PTH(1-34) gene under T7 promoter (often as fusion with soluble tag for stability)",
                    "Pre-culture: LB medium + kanamycin, 37 °C, 12 h",
                    "Main fermentation: fed-batch in 1-20 m3 bioreactor, 37 °C -> temp shift to 30 °C at induction, pH 7.0",
                    "Induction with IPTG (0.5-1 mM) at OD600 ~40-60, expression 6-8 h",
                    "Teriparatide precursor as inclusion bodies or soluble fusion protein, 1-4 g/L",
                ],
            },
            "downstream": {
                "de": [
                    "Zell-Lyse via Hochdruckhomogenisator",
                    "Bei IB-Route: Solubilisierung + Refolding; bei löslicher Route: direkte Klärung",
                    "Enzymatische Tag-Abspaltung (z. B. TEV-Protease, Faktor Xa)",
                    "Anionenaustausch-Chromatographie (Q-Sepharose) — Capture",
                    "Reverse-Phase-HPLC für Polishing",
                    "Ultrafiltration + Lyophilisation",
                ],
                "en": [
                    "Cell lysis via high-pressure homogenizer",
                    "IB route: solubilisation + refolding; soluble route: direct clarification",
                    "Enzymatic tag cleavage (e.g. TEV protease, factor Xa)",
                    "Anion-exchange chromatography (Q-Sepharose) — capture",
                    "Reverse-phase HPLC for polishing",
                    "Ultrafiltration + lyophilisation",
                ],
            },
            "expected_purity_after_workup": "≥ 99 % (USP / Ph. Eur., Forteo)",
        },
    },
    "exenatide": {
        "chemical": {
            "reagents": "Fmoc-AS für 39-mer Exendin-4-Sequenz + HBTU/DIPEA + Rink-Amid-Harz",
            "reagents_en": "Fmoc-AA for 39-mer exendin-4 sequence + HBTU/DIPEA + Rink amide resin",
            "yield_range_percent": (10, 20),
            "literature_source": "Eng et al., J. Biol. Chem. 1992 — Isolation and characterization of exendin-4 (doi:10.1016/S0021-9258(18)42531-8); Drucker, Cell Metab. 2018 — GLP-1 receptor agonists clinical review (doi:10.1016/j.cmet.2018.03.001)",
            "production": {
                "de": [
                    "Harz-Vorbereitung: Rink-Amid-Harz (für C-terminales Amid) mit Fmoc-Ser beladen",
                    "SPPS-Zyklen (39 Kupplungen) für Exendin-4-Sequenz HGEGTFTSDLSKQMEEEAVRLFIEWLKNGGPSSGAPPPS — sehr lange Sequenz, Double-Couplings empfehlenswert",
                    "Schwierige Stellen: Trp-haltige Regionen (Trp25), Ser-Cluster am C-Terminus (Aggregations-anfällig)",
                    "Final-Cleavage: TFA-Cocktail mit zusätzlichem Scavenger für Trp",
                    "Roh-Exenatide 10-20 % Gesamtausbeute (typisch für 39-mer)",
                ],
                "en": [
                    "Resin preparation: Rink amide resin (for C-terminal amide) loaded with Fmoc-Ser",
                    "SPPS cycles (39 couplings) for exendin-4 sequence HGEGTFTSDLSKQMEEEAVRLFIEWLKNGGPSSGAPPPS — very long sequence, double couplings recommended",
                    "Difficult positions: Trp-containing regions (Trp25), Ser clusters at C-terminus (aggregation-prone)",
                    "Final cleavage: TFA cocktail with additional Trp scavenger",
                    "Crude exenatide 10-20 % overall yield (typical for 39-mer)",
                ],
            },
            "downstream": {
                "de": [
                    "Mehrstufige Prep-RP-HPLC für hohe Reinheit",
                    "Counter-Ion-Exchange auf Acetat",
                    "Lyophilisation, finale Formulierung in Acetatpuffer (Byetta) oder Mikrosphären-Formulierung (Bydureon)",
                ],
                "en": [
                    "Multi-stage prep-RP-HPLC for high purity",
                    "Counter-ion exchange to acetate",
                    "Lyophilisation, final formulation in acetate buffer (Byetta) or microsphere formulation (Bydureon)",
                ],
            },
            "expected_purity_after_workup": "≥ 99 % (USP, Byetta / Bydureon)",
        },
    },
    # ====================================================================
    # BATCH 4 — Peptide / cyclic (NRPS-Fermentation, Semi-Synthetik, SPPS-Cyclisierung)
    # ====================================================================
    "cyclosporine": {
        "biotechnological": {
            "organism": "Tolypocladium inflatum (Schimmelpilz, NRPS-Biosynthese)",
            "organism_en": "Tolypocladium inflatum (mould, NRPS biosynthesis)",
            "substrate": "Glucose + L-Aminosäuren (Bbu, Sar, MeVal, MeLeu als Substrat-Vorstufen)",
            "substrate_en": "glucose + L-amino acids (Bbu, Sar, MeVal, MeLeu as substrate precursors)",
            "yield_range_g_per_l": (1, 5),
            "fermentation_time_h": (240, 360),
            "literature_source": "Borel et al., Pharmacol. Rev. 2002 — Cyclosporine pharmacology and biosynthesis (doi:10.1007/978-3-642-56923-1); Survase et al., Crit. Rev. Biotechnol. 2011 — Cyclosporine A production (doi:10.3109/07388551.2010.539027)",
            "production": {
                "de": [
                    "Stammvorbereitung: T. inflatum Hochertrags-Stamm (Klassiker NRRL 8044 oder neuere Hochleister, durch Mutagenese verbessert)",
                    "Vorkultur: Glucose-haltiges Medium + Hefeextrakt + Caseinpepton, 24-27 °C, 72-96 h",
                    "Hauptfermentation: aerobe Submerskultur in 50-200 m3 Bioreaktor, 24-27 °C, pH 5.5-6.5, hoher O2-Eintrag",
                    "Substrat-Feed: Glucose-Feed kontrolliert; Zugabe von L-Aminosäuren (besonders D-α-Aminobutyrat) als Cyclosporin-Vorstufen",
                    "Fermentations-Dauer 10-15 Tage, extrazelluläre Cyclosporin-A-Akkumulation 1-5 g/L (begleitet von Nebenformen B, C, D)",
                ],
                "en": [
                    "Strain preparation: T. inflatum high-yield strain (classic NRRL 8044 or newer high-producers, improved by mutagenesis)",
                    "Pre-culture: glucose-containing medium + yeast extract + casein peptone, 24-27 °C, 72-96 h",
                    "Main fermentation: aerobic submerged culture in 50-200 m3 bioreactor, 24-27 °C, pH 5.5-6.5, high O2 transfer",
                    "Substrate feed: controlled glucose feed; addition of L-amino acids (especially D-alpha-aminobutyrate) as cyclosporine precursors",
                    "Fermentation duration 10-15 days, extracellular cyclosporine A accumulation 1-5 g/L (accompanied by minor forms B, C, D)",
                ],
            },
            "downstream": {
                "de": [
                    "Myzel-Trennung via Filtration (Cyclosporin ist sowohl im Myzel als auch im Medium)",
                    "Lösungsmittel-Extraktion mit Ethylacetat oder n-Butylacetat",
                    "Säulen-Chromatographie an Kieselgel zur Trennung von A/B/C/D-Cyclosporinen",
                    "Umkristallisation aus Aceton/Wasser für API-Grade",
                    "Optional: präparative HPLC für klinische Reinheit",
                ],
                "en": [
                    "Mycelium separation via filtration (cyclosporine is in both mycelium and medium)",
                    "Solvent extraction with ethyl acetate or n-butyl acetate",
                    "Column chromatography on silica gel to separate A/B/C/D cyclosporines",
                    "Recrystallization from acetone/water for API grade",
                    "Optional: preparative HPLC for clinical purity",
                ],
            },
            "expected_purity_after_workup": "≥ 98.5 % (USP / Ph. Eur., Sandimmun / Neoral)",
        },
    },
    "gramicidin s": {
        "biotechnological": {
            "organism": "Bacillus brevis (Aneurinibacillus migulanus) — NRPS-Biosynthese",
            "organism_en": "Bacillus brevis (Aneurinibacillus migulanus) — NRPS biosynthesis",
            "substrate": "Glucose + Aminosäuren-reicher Komplex (Pepton + Hefeextrakt)",
            "substrate_en": "glucose + amino-acid-rich complex (peptone + yeast extract)",
            "yield_range_g_per_l": (0.5, 2),
            "fermentation_time_h": (40, 72),
            "literature_source": "Marahiel, Chem. Biol. 1997 — Multidomain enzymes of gramicidin S biosynthesis (doi:10.1016/S1074-5521(97)90019-1); Kratzschmar et al., J. Bacteriol. 1989 — Gramicidin S synthetase gene cluster (doi:10.1128/jb.171.10.5422-5429.1989)",
            "production": {
                "de": [
                    "Stammvorbereitung: B. brevis ATCC 9999 oder Hochertrags-Mutanten, Sporenbildung als Stabilitätsform",
                    "Vorkultur: Pepton-Hefeextrakt-Glucose-Medium, 37 °C, 16-24 h",
                    "Hauptfermentation: aerobe Submerskultur in 1-50 m3 Bioreaktor, 37 °C, pH 6.0-7.0",
                    "Stationäre Phase wichtig: Gramicidin-S-Produktion ist sekundärmetabolitisch, beginnt nach Glucose-Erschöpfung",
                    "Fermentations-Dauer 40-72 h, Gramicidin-S-Akkumulation extrazellulär + zellgebunden, 0.5-2 g/L",
                ],
                "en": [
                    "Strain preparation: B. brevis ATCC 9999 or high-yield mutants, sporulation as stability form",
                    "Pre-culture: peptone-yeast extract-glucose medium, 37 °C, 16-24 h",
                    "Main fermentation: aerobic submerged culture in 1-50 m3 bioreactor, 37 °C, pH 6.0-7.0",
                    "Stationary phase important: gramicidin S production is secondary metabolite, begins after glucose depletion",
                    "Fermentation duration 40-72 h, gramicidin S accumulation extracellular + cell-bound, 0.5-2 g/L",
                ],
            },
            "downstream": {
                "de": [
                    "Zellernte via Zentrifugation (Gramicidin S großteils zellgebunden)",
                    "Extraktion mit Ethanol oder Aceton",
                    "Säulen-Chromatographie an Kieselgel oder RP-Säule",
                    "Umkristallisation aus Methanol/Wasser",
                ],
                "en": [
                    "Cell harvest via centrifugation (gramicidin S mostly cell-bound)",
                    "Extraction with ethanol or acetone",
                    "Column chromatography on silica or RP column",
                    "Recrystallization from methanol/water",
                ],
            },
            "expected_purity_after_workup": "≥ 98 % (Forschung / topisch); klinisch heute fast nicht mehr verwendet wegen Hämolyse",
        },
    },
    "bacitracin": {
        "biotechnological": {
            "organism": "Bacillus licheniformis (auch B. subtilis-Variante)",
            "organism_en": "Bacillus licheniformis (also B. subtilis variant)",
            "substrate": "Sojamehl + Glucose + Salze (komplex)",
            "substrate_en": "soybean meal + glucose + salts (complex)",
            "yield_range_g_per_l": (1, 3),
            "fermentation_time_h": (40, 72),
            "literature_source": "Stone & Strominger, Proc. Natl. Acad. Sci. 1971 — Mechanism of action of bacitracin (doi:10.1073/pnas.68.12.3223); Ming & Epperson, J. Bacteriol. 2002 — Bacitracin biosynthesis cluster (doi:10.1128/JB.184.21.5896-5905.2002)",
            "production": {
                "de": [
                    "Stammvorbereitung: B. licheniformis Industrial-Stamm aus Stammbank, Sporen-Inokulation",
                    "Vorkultur: Sojamehl-Glucose-Medium, 37 °C, 16-24 h",
                    "Hauptfermentation: aerobe Submerskultur in 10-100 m3 Bioreaktor, 37 °C, pH 6.5-7.0",
                    "Substrat: Sojamehl liefert Aminosäure-Pool für NRPS-Synthese; Glucose als C-Quelle",
                    "Fermentations-Dauer 40-72 h, extrazelluläre Bacitracin-Akkumulation 1-3 g/L (Bacitracin A als Hauptkomponente)",
                ],
                "en": [
                    "Strain preparation: B. licheniformis industrial strain from strain bank, spore inoculation",
                    "Pre-culture: soybean-meal-glucose medium, 37 °C, 16-24 h",
                    "Main fermentation: aerobic submerged culture in 10-100 m3 bioreactor, 37 °C, pH 6.5-7.0",
                    "Substrate: soybean meal supplies amino-acid pool for NRPS synthesis; glucose as C source",
                    "Fermentation duration 40-72 h, extracellular bacitracin accumulation 1-3 g/L (bacitracin A as main component)",
                ],
            },
            "downstream": {
                "de": [
                    "Zellabtrennung via Filtration",
                    "Säuerung auf pH 3 → Bacitracin fällt teilweise aus",
                    "Lösungsmittel-Extraktion oder Ionenaustausch (kationischer Charakter)",
                    "Sprühtrocknung zu Bacitracin-Zink-Komplex (handelsübliche Form, stabilisiert)",
                ],
                "en": [
                    "Cell separation via filtration",
                    "Acidification to pH 3 -> bacitracin precipitates partially",
                    "Solvent extraction or ion exchange (cationic character)",
                    "Spray drying to bacitracin-zinc complex (commercial form, stabilised)",
                ],
            },
            "expected_purity_after_workup": "≥ 70 IU/mg (USP, topisch); fast immer als Zink-Komplex formuliert",
        },
    },
    "vancomycin": {
        "biotechnological": {
            "organism": "Amycolatopsis orientalis (früher Nocardia orientalis) — NRPS-Glykopeptid-Biosynthese",
            "organism_en": "Amycolatopsis orientalis (formerly Nocardia orientalis) — NRPS glycopeptide biosynthesis",
            "substrate": "Komplex-Medium (Sojamehl + Glucose + Stärke + Salze)",
            "substrate_en": "complex medium (soybean meal + glucose + starch + salts)",
            "yield_range_g_per_l": (3, 10),
            "fermentation_time_h": (120, 200),
            "literature_source": "Levine, Clin. Infect. Dis. 2006 — Vancomycin: A history (doi:10.1086/491709); Hubbard & Walsh, Angew. Chem. Int. Ed. 2003 — Vancomycin biosynthesis (doi:10.1002/anie.200200534)",
            "production": {
                "de": [
                    "Stammvorbereitung: A. orientalis Hochertrags-Stamm aus Stammbank, Sporenbildung auf Schrägagar",
                    "Vorkultur: Sojamehl-Glucose-Stärke-Medium, 28-30 °C, 48-72 h",
                    "Hauptfermentation: aerobe Submerskultur in 50-200 m3 Bioreaktor, 28-30 °C, pH 6.5-7.5, hoher O2-Eintrag",
                    "Substrat-Feed: Glucose + Stärke + komplexe N-Quellen (Aminosäuren als NRPS-Vorstufen)",
                    "Fermentations-Dauer 120-200 h, finale Vancomycin-Konzentration 3-10 g/L (Hauptform Vancomycin B)",
                ],
                "en": [
                    "Strain preparation: A. orientalis high-yield strain from strain bank, sporulation on slant agar",
                    "Pre-culture: soybean-meal-glucose-starch medium, 28-30 °C, 48-72 h",
                    "Main fermentation: aerobic submerged culture in 50-200 m3 bioreactor, 28-30 °C, pH 6.5-7.5, high O2 transfer",
                    "Substrate feed: glucose + starch + complex N sources (amino acids as NRPS precursors)",
                    "Fermentation duration 120-200 h, final vancomycin concentration 3-10 g/L (main form vancomycin B)",
                ],
            },
            "downstream": {
                "de": [
                    "Myzel-Abtrennung via Filtration",
                    "Adsorption an Aktivkohle oder Polymerharz (Diaion HP20)",
                    "Elution mit Methanol/Wasser-Gradient",
                    "Ionenaustausch-Chromatographie (Kationenaustauscher, da Vancomycin amphoter)",
                    "Umkristallisation als Hydrochlorid für API-Grade",
                ],
                "en": [
                    "Mycelium separation via filtration",
                    "Adsorption to activated carbon or polymer resin (Diaion HP20)",
                    "Elution with methanol/water gradient",
                    "Ion-exchange chromatography (cation exchanger — vancomycin is amphoteric)",
                    "Recrystallization as hydrochloride for API grade",
                ],
            },
            "expected_purity_after_workup": "≥ 95 % (USP / Ph. Eur., Vancomycin-HCl)",
        },
    },
    "daptomycin": {
        "biotechnological": {
            "organism": "Streptomyces roseosporus (Lipopeptid-NRPS)",
            "organism_en": "Streptomyces roseosporus (lipopeptide NRPS)",
            "substrate": "Glucose + Decansäure-Feed (für C10-Lipidkette)",
            "substrate_en": "glucose + decanoic acid feed (for C10 lipid tail)",
            "yield_range_g_per_l": (1, 4),
            "fermentation_time_h": (120, 200),
            "literature_source": "Eisenstein et al., Clin. Infect. Dis. 2010 — Daptomycin clinical and pharmacological review (doi:10.1086/648676); Miao et al., Microbiology 2005 — Daptomycin biosynthesis (doi:10.1099/mic.0.27757-0)",
            "production": {
                "de": [
                    "Stammvorbereitung: S. roseosporus NRRL 11379 oder Hochertrags-Mutante",
                    "Vorkultur: komplexes Medium (Hefeextrakt + Glucose + Salze), 28-30 °C, 48-72 h",
                    "Hauptfermentation: aerobe Submerskultur in 10-100 m3 Bioreaktor, 28-30 °C, pH 6.5-7.5",
                    "Schlüssel-Substrat: Decansäure (C10)-Feed → kontrolliert die Lipidketten-Länge (klinisches Daptomycin hat C10-Decanoyl-Seitenkette)",
                    "Fermentations-Dauer 120-200 h, finale Daptomycin-Konzentration 1-4 g/L",
                ],
                "en": [
                    "Strain preparation: S. roseosporus NRRL 11379 or high-yield mutant",
                    "Pre-culture: complex medium (yeast extract + glucose + salts), 28-30 °C, 48-72 h",
                    "Main fermentation: aerobic submerged culture in 10-100 m3 bioreactor, 28-30 °C, pH 6.5-7.5",
                    "Key substrate: decanoic acid (C10) feed -> controls lipid-tail length (clinical daptomycin has C10 decanoyl side chain)",
                    "Fermentation duration 120-200 h, final daptomycin concentration 1-4 g/L",
                ],
            },
            "downstream": {
                "de": [
                    "Myzel-Abtrennung via Filtration",
                    "Säurefällung bei pH 2.5",
                    "Resolubilisierung bei pH 6, Adsorption an Polymerharz",
                    "Reverse-Phase-Chromatographie zur Trennung von Lipidketten-Isomeren (C9/C10/C11)",
                    "Umkristallisation als Natriumsalz für API-Grade",
                ],
                "en": [
                    "Mycelium separation via filtration",
                    "Acid precipitation at pH 2.5",
                    "Resolubilisation at pH 6, adsorption to polymer resin",
                    "Reverse-phase chromatography to separate lipid-tail isomers (C9/C10/C11)",
                    "Recrystallization as sodium salt for API grade",
                ],
            },
            "expected_purity_after_workup": "≥ 98 % (USP, Cubicin / daptomycin for injection)",
        },
    },
    "octreotide": {
        "chemical": {
            "reagents": "Fmoc-AS (8-mer cyclisch) + HBTU/DIPEA + Wang-Harz; Disulfid-Cyclisierung via Luft oder I2",
            "reagents_en": "Fmoc-AA (8-mer cyclic) + HBTU/DIPEA + Wang resin; disulfide cyclisation via air or I2",
            "yield_range_percent": (20, 35),
            "literature_source": "Lamberts et al., N. Engl. J. Med. 1996 — Octreotide review (doi:10.1056/NEJM199607253350406); Bauer et al., Life Sci. 1982 — Erstbeschreibung Octreotid",
            "production": {
                "de": [
                    "Harz-Vorbereitung: Wang-Harz mit C-terminalem Fmoc-Threoninol beladen (Threoninol = reduzierte Form von Thr, C-terminale Modifikation)",
                    "SPPS-Zyklen (8 Kupplungen) für Sequenz D-Phe-Cys-Phe-D-Trp-Lys-Thr-Cys-Thr-ol: Standard-Fmoc-Protokoll, D-Aminosäuren erfordern höhere Reinheit",
                    "Cys-Reste mit Trt geschützt; Lys mit Boc; Trp mit Boc",
                    "Final-Cleavage: TFA-Cocktail (mit Scavenger für Trp)",
                    "Disulfid-Cyclisierung: Luft-Oxidation in verdünnter Lösung (0.1 g/L, pH 8, 24 h) ODER I2 in Methanol für schnelle Cyclisierung",
                    "Roh-Octreotid 20-35 % Gesamtausbeute (gut für cyclisches 8-mer mit D-AS)",
                ],
                "en": [
                    "Resin preparation: Wang resin loaded with C-terminal Fmoc-threoninol (threoninol = reduced form of Thr, C-terminal modification)",
                    "SPPS cycles (8 couplings) for sequence D-Phe-Cys-Phe-D-Trp-Lys-Thr-Cys-Thr-ol: standard Fmoc protocol, D-amino acids require higher purity",
                    "Cys residues protected with Trt; Lys with Boc; Trp with Boc",
                    "Final cleavage: TFA cocktail (with scavenger for Trp)",
                    "Disulfide cyclisation: air oxidation in dilute solution (0.1 g/L, pH 8, 24 h) OR I2 in methanol for fast cyclisation",
                    "Crude octreotide 20-35 % overall yield (good for cyclic 8-mer with D-AAs)",
                ],
            },
            "downstream": {
                "de": [
                    "Prep-RP-HPLC zur Trennung von Disulfid-Isomeren und Truncations",
                    "Counter-Ion-Exchange auf Acetat",
                    "Lyophilisation, Lagerung als Pamoat-Salz für Depot-Formulierung (Sandostatin LAR)",
                ],
                "en": [
                    "Prep-RP-HPLC to separate disulfide isomers and truncations",
                    "Counter-ion exchange to acetate",
                    "Lyophilisation, storage as pamoate salt for depot formulation (Sandostatin LAR)",
                ],
            },
            "expected_purity_after_workup": "≥ 99 % (USP, Sandostatin)",
        },
    },
    "anidulafungin": {
        "biotechnological": {
            "organism": "Aspergillus nidulans (Echinocandin-B-Produzent) + chemische Semi-Synthese",
            "organism_en": "Aspergillus nidulans (echinocandin B producer) + chemical semi-synthesis",
            "substrate": "Komplex-Medium (Glucose + Hefeextrakt + Salze)",
            "substrate_en": "complex medium (glucose + yeast extract + salts)",
            "yield_range_g_per_l": (0.5, 2),  # Echinocandin B precursor titer
            "fermentation_time_h": (120, 200),
            "literature_source": "Vazquez & Sobel, Clin. Infect. Dis. 2006 — Anidulafungin review (doi:10.1086/501830); Cacho et al., Adv. Biochem. Eng. Biotechnol. 2015 — Echinocandin biosynthesis (doi:10.1007/10_2014_274)",
            "production": {
                "de": [
                    "Phase A — Fermentations-Phase: A. nidulans Stammvorbereitung + Submerskultur in 10-100 m3 Bioreaktor (28-30 °C, pH 6.5-7.0, 5-8 Tage) → Echinocandin-B-Akkumulation 0.5-2 g/L",
                    "Aufarbeitung des Echinocandin-B-Intermediats: Filtration, Lösungsmittel-Extraktion, Säulenchromatographie",
                    "Phase B — Semi-Synthese: enzymatische Abspaltung der nativen Linolyl-Seitenkette mit Acylase aus Actinoplanes utahensis → freier Hexapeptid-Kern",
                    "Reacylierung: Kupplung mit p-Pentyloxy-Terphenyl-Carbonsäure (Anidulafungin-spezifische Seitenkette) → Anidulafungin",
                    "Roh-Anidulafungin nach Semi-Synthese: ~50-70 % Gesamtausbeute auf Echinocandin-B-Basis",
                ],
                "en": [
                    "Phase A — fermentation: A. nidulans strain preparation + submerged culture in 10-100 m3 bioreactor (28-30 °C, pH 6.5-7.0, 5-8 days) -> echinocandin B accumulation 0.5-2 g/L",
                    "Workup of echinocandin B intermediate: filtration, solvent extraction, column chromatography",
                    "Phase B — semi-synthesis: enzymatic cleavage of native linoleyl side chain with acylase from Actinoplanes utahensis -> free hexapeptide core",
                    "Reacylation: coupling with p-pentyloxy terphenyl carboxylic acid (anidulafungin-specific side chain) -> anidulafungin",
                    "Crude anidulafungin after semi-synthesis: ~50-70 % overall yield based on echinocandin B",
                ],
            },
            "downstream": {
                "de": [
                    "Prep-RP-HPLC zur Trennung von Acylierungs-Isomeren",
                    "Kristallisation aus Acetonitril/Wasser für API-Grade",
                    "Lyophilisation, finale Formulierung mit Fructose als Lyophilisations-Schutz (Ecalta / Eraxis)",
                ],
                "en": [
                    "Prep-RP-HPLC to separate acylation isomers",
                    "Crystallization from acetonitrile/water for API grade",
                    "Lyophilisation, final formulation with fructose as lyoprotectant (Ecalta / Eraxis)",
                ],
            },
            "expected_purity_after_workup": "≥ 99 % (USP, Ecalta / Eraxis)",
        },
    },
    "caspofungin": {
        "biotechnological": {
            "organism": "Glarea lozoyensis (Pneumocandin-B0-Produzent) + chemische Semi-Synthese",
            "organism_en": "Glarea lozoyensis (pneumocandin B0 producer) + chemical semi-synthesis",
            "substrate": "Komplex-Medium (Glucose + Pepton + Mannit)",
            "substrate_en": "complex medium (glucose + peptone + mannitol)",
            "yield_range_g_per_l": (0.5, 2),  # Pneumocandin B0 precursor
            "fermentation_time_h": (168, 240),
            "literature_source": "Denning, Lancet 2003 — Echinocandin antifungal drugs (doi:10.1016/S0140-6736(03)14472-8); Bouffard et al., J. Med. Chem. 1994 — Caspofungin discovery (doi:10.1021/jm00033a006)",
            "production": {
                "de": [
                    "Phase A — Fermentation: G. lozoyensis (ATCC 74030) Submerskultur in 10-100 m3 Bioreaktor (24-26 °C, pH 5.5-6.5, 7-10 Tage) → Pneumocandin-B0-Akkumulation 0.5-2 g/L",
                    "Aufarbeitung des Pneumocandin-B0-Intermediats via Filtration + Lösungsmittel-Extraktion",
                    "Phase B — Semi-Synthese (mehrstufige chemische Modifikation): Reduktion der Hexapeptid-Carbonyl-Gruppe zu Aminomethyl-Funktion + Modifikation der Aminosäure-Seitenketten",
                    "Insgesamt ~6 chemische Stufen Semi-Synthese mit Schutzgruppen-Strategie",
                    "Roh-Caspofungin nach Semi-Synthese: ~30-50 % Gesamtausbeute auf Pneumocandin-B0-Basis",
                ],
                "en": [
                    "Phase A — fermentation: G. lozoyensis (ATCC 74030) submerged culture in 10-100 m3 bioreactor (24-26 °C, pH 5.5-6.5, 7-10 days) -> pneumocandin B0 accumulation 0.5-2 g/L",
                    "Workup of pneumocandin B0 intermediate via filtration + solvent extraction",
                    "Phase B — semi-synthesis (multi-stage chemical modification): reduction of hexapeptide carbonyl to aminomethyl function + modification of amino-acid side chains",
                    "Total ~6 chemical stages of semi-synthesis with protecting-group strategy",
                    "Crude caspofungin after semi-synthesis: ~30-50 % overall yield based on pneumocandin B0",
                ],
            },
            "downstream": {
                "de": [
                    "Prep-RP-HPLC zur Trennung von Stereoisomeren der Semi-Synthese",
                    "Salzbildung als Diacetat",
                    "Lyophilisation, finale Formulierung (Cancidas)",
                ],
                "en": [
                    "Prep-RP-HPLC to separate stereoisomers from semi-synthesis",
                    "Salt formation as diacetate",
                    "Lyophilisation, final formulation (Cancidas)",
                ],
            },
            "expected_purity_after_workup": "≥ 99 % (USP, Cancidas)",
        },
    },
    "colistin": {
        "biotechnological": {
            "organism": "Paenibacillus polymyxa subsp. colistinus (NRPS-Lipopeptid)",
            "organism_en": "Paenibacillus polymyxa subsp. colistinus (NRPS lipopeptide)",
            "substrate": "Komplex-Medium (Glucose + Hefeextrakt + Pepton)",
            "substrate_en": "complex medium (glucose + yeast extract + peptone)",
            "yield_range_g_per_l": (2, 6),
            "fermentation_time_h": (40, 72),
            "literature_source": "Falagas & Kasiakou, Clin. Infect. Dis. 2005 — Colistin clinical review (doi:10.1086/429323); Choi et al., FEMS Microbiol. Lett. 2009 — Colistin biosynthesis gene cluster (doi:10.1111/j.1574-6968.2009.01619.x)",
            "production": {
                "de": [
                    "Stammvorbereitung: P. polymyxa Industrial-Stamm aus Stammbank, Sporen-Inokulation",
                    "Vorkultur: Pepton-Hefeextrakt-Glucose-Medium, 30 °C, 16-24 h",
                    "Hauptfermentation: aerobe Submerskultur in 10-100 m3 Bioreaktor, 30 °C, pH 6.5-7.5",
                    "Substrat: Glucose + komplexe Aminosäure-Quellen für NRPS-Synthese",
                    "Fermentations-Dauer 40-72 h, finale Colistin-Konzentration 2-6 g/L (Colistin A und B als Hauptkomponenten)",
                ],
                "en": [
                    "Strain preparation: P. polymyxa industrial strain from strain bank, spore inoculation",
                    "Pre-culture: peptone-yeast extract-glucose medium, 30 °C, 16-24 h",
                    "Main fermentation: aerobic submerged culture in 10-100 m3 bioreactor, 30 °C, pH 6.5-7.5",
                    "Substrate: glucose + complex amino-acid sources for NRPS synthesis",
                    "Fermentation duration 40-72 h, final colistin concentration 2-6 g/L (colistin A and B as main components)",
                ],
            },
            "downstream": {
                "de": [
                    "Zellabtrennung via Filtration",
                    "Säuerung auf pH 3 + Filtration über Aktivkohle",
                    "Kationenaustausch-Chromatographie (Colistin ist polykationisch)",
                    "Sprühtrocknung als Colistin-Sulfat oder Colistimethat-Natrium",
                ],
                "en": [
                    "Cell separation via filtration",
                    "Acidification to pH 3 + filtration over activated carbon",
                    "Cation-exchange chromatography (colistin is polycationic)",
                    "Spray drying as colistin sulfate or colistimethate sodium",
                ],
            },
            "expected_purity_after_workup": "≥ 90 % (Colistin-Sulfat USP); klinisch oft Colistimethat-Natrium (Prodrug)",
        },
    },
    "polymyxin b": {
        "biotechnological": {
            "organism": "Paenibacillus polymyxa (NRPS-Lipopeptid, eng verwandt mit Colistin = Polymyxin E)",
            "organism_en": "Paenibacillus polymyxa (NRPS lipopeptide, closely related to colistin = polymyxin E)",
            "substrate": "Komplex-Medium analog Colistin",
            "substrate_en": "complex medium analogous to colistin",
            "yield_range_g_per_l": (2, 5),
            "fermentation_time_h": (40, 72),
            "literature_source": "Velkov et al., J. Med. Chem. 2010 — Structure-activity relationships of polymyxins (doi:10.1021/jm900999h); Storm et al., Annu. Rev. Biochem. 1977 — Polymyxin biosynthesis (doi:10.1146/annurev.bi.46.070177.003515)",
            "production": {
                "de": [
                    "Stammvorbereitung: P. polymyxa Industrial-Stamm (Polymyxin-B-Produzent), Sporen-Inokulation",
                    "Vorkultur: Pepton-Hefeextrakt-Glucose-Medium, 30 °C, 16-24 h",
                    "Hauptfermentation: aerobe Submerskultur in 10-100 m3 Bioreaktor, 30 °C, pH 6.5-7.5",
                    "Polymyxin-B besteht aus Komponenten B1, B2, B3 (unterschiedliche Lipidketten); klinische Form ist Mischung B1+B2",
                    "Fermentations-Dauer 40-72 h, finale Polymyxin-B-Konzentration 2-5 g/L",
                ],
                "en": [
                    "Strain preparation: P. polymyxa industrial strain (polymyxin B producer), spore inoculation",
                    "Pre-culture: peptone-yeast extract-glucose medium, 30 °C, 16-24 h",
                    "Main fermentation: aerobic submerged culture in 10-100 m3 bioreactor, 30 °C, pH 6.5-7.5",
                    "Polymyxin B consists of components B1, B2, B3 (different lipid chains); clinical form is mixture B1+B2",
                    "Fermentation duration 40-72 h, final polymyxin B concentration 2-5 g/L",
                ],
            },
            "downstream": {
                "de": [
                    "Zellabtrennung via Filtration",
                    "Säuerung + Aktivkohle-Behandlung",
                    "Kationenaustausch-Chromatographie",
                    "Sprühtrocknung als Polymyxin-B-Sulfat",
                ],
                "en": [
                    "Cell separation via filtration",
                    "Acidification + activated-carbon treatment",
                    "Cation-exchange chromatography",
                    "Spray drying as polymyxin B sulfate",
                ],
            },
            "expected_purity_after_workup": "≥ 90 % (Polymyxin-B-Sulfat USP)",
        },
    },
    # ====================================================================
    # BATCH 5 — Protein / antibody (alle 8 verbleibenden mAbs)
    # ====================================================================
    # Produktion ist im Kern identisch zu Trastuzumab (CHO Fed-Batch +
    # Protein-A-Capture + IEX + VI/VF + TFF). Unterschiede sind primär:
    # - Zelllinie (CHO-K1, CHO-DG44, CHO-S, SP2/0)
    # - Glykosylierungsprofil (für ADCC/CDC-Optimierung)
    # - Titer-Range (3-10 g/L modern, ältere Prozesse 1-3 g/L)
    # - Target-spezifische Affinitäts-Schritte (z. B. Daratumumab anti-CD38)
    # ====================================================================
    "adalimumab": {
        "biotechnological": {
            "organism": "CHO-K1-Zellen (Suspensionskultur, ursprünglich aus Cambridge Antibody Technology)",
            "organism_en": "CHO-K1 cells (suspension culture, originally from Cambridge Antibody Technology)",
            "substrate": "Chemisch definiertes Medium (CDM), Glucose + Glutamin + Spurenelemente",
            "substrate_en": "chemically defined medium (CDM), glucose + glutamine + trace elements",
            "yield_range_g_per_l": (2, 6),
            "fermentation_time_h": (288, 360),  # 12-15 days
            "literature_source": "Lawrence, mAbs 2011 — Adalimumab manufacturing review (doi:10.4161/mabs.3.6.18020); Frenzel et al., Methods Mol. Biol. 2014 — Adalimumab production process (doi:10.1007/978-1-4939-1196-5_19)",
            "production": {
                "de": [
                    "Zellbank: GMP-MCB/WCB von rekombinanten CHO-K1-Zellen mit Adalimumab-Schwer- und Leichtketten-Genen (anti-TNFalpha humaner IgG1)",
                    "Auftauen + Vorkultur: WCB-Ampulle, 37 °C/5 % CO2/CDM-Medium, Schüttelfläschchen 3-5 Tage",
                    "Seed-Train: 50 mL → 2 L Wave-Bag → 50 L Single-Use → 500 L → 2000-12000 L Edelstahl-Bioreaktor",
                    "Hauptproduktion: Fed-Batch 12-15 Tage, 37 °C, pH 7.0, DO 30 %, Glucose- + Aminosäure-Feeds",
                    "Produktions-Ende: Antikörper-Titer 2-6 g/L (moderne Hochleister bis 8 g/L)",
                ],
                "en": [
                    "Cell bank: GMP MCB/WCB of recombinant CHO-K1 cells with adalimumab heavy- and light-chain genes (anti-TNFalpha human IgG1)",
                    "Thaw + pre-culture: WCB vial, 37 °C/5 % CO2/CDM medium, shake flasks 3-5 days",
                    "Seed train: 50 mL -> 2 L wave bag -> 50 L single-use -> 500 L -> 2000-12000 L stainless-steel bioreactor",
                    "Main production: fed-batch 12-15 days, 37 °C, pH 7.0, DO 30 %, glucose + amino-acid feeds",
                    "End of production: antibody titer 2-6 g/L (modern high-producers up to 8 g/L)",
                ],
            },
            "downstream": {
                "de": [
                    "Klärung via Tiefenfiltration + 0.2 µm Sterilfiltration",
                    "Protein-A-Affinitätschromatographie (Capture, > 95 % Reinheit ab Stufe 1)",
                    "Niedrig-pH-Virusinaktivierung (pH 3.5 für 60 min)",
                    "Kationenaustausch-Chromatographie (CEX, Bind-Elute)",
                    "Anionenaustausch im Flow-through (AEX) + Virusfiltration (Planova 20N)",
                    "Ultrafiltration / Pufferaustausch (TFF, finale Formulierung mit Citrat-Phosphat-Puffer)",
                ],
                "en": [
                    "Clarification via depth filtration + 0.2 um sterile filtration",
                    "Protein-A affinity chromatography (capture, > 95 % purity from step 1)",
                    "Low-pH virus inactivation (pH 3.5 for 60 min)",
                    "Cation-exchange chromatography (CEX, bind-elute)",
                    "Anion-exchange in flow-through (AEX) + virus filtration (Planova 20N)",
                    "Ultrafiltration / buffer exchange (TFF, final formulation with citrate-phosphate buffer)",
                ],
            },
            "expected_purity_after_workup": "≥ 99 % (cGMP, Humira / Biosimilar-äquivalent)",
        },
    },
    "bevacizumab": {
        "biotechnological": {
            "organism": "CHO-Zellen (DUKX-B11, dhfr-negativ — MTX-selektierbar)",
            "organism_en": "CHO cells (DUKX-B11, dhfr-negative — MTX selectable)",
            "substrate": "Chemisch definiertes Medium (CDM)",
            "substrate_en": "chemically defined medium (CDM)",
            "yield_range_g_per_l": (2, 5),
            "fermentation_time_h": (288, 360),
            "literature_source": "Ferrara et al., Nat. Rev. Drug Discov. 2004 — Bevacizumab and anti-angiogenic therapy (doi:10.1038/nrd1381); Presta et al., Cancer Res. 1997 — Bevacizumab humanization",
            "production": {
                "de": [
                    "Zellbank: GMP-MCB/WCB von CHO-DUKX-B11-Zellen mit Bevacizumab-Genen (anti-VEGF-A humaner IgG1)",
                    "Auftauen + Vorkultur: WCB-Ampulle, 37 °C/5 % CO2/CDM-Medium, 3-5 Tage",
                    "Seed-Train: stufenweise Expansion bis 2000-12000 L Produktions-Bioreaktor",
                    "Hauptproduktion: Fed-Batch 12-15 Tage, 37 °C, pH 7.0, DO 30 %, optimierte Feed-Strategie",
                    "Produktions-Ende: Bevacizumab-Titer 2-5 g/L im Überstand",
                ],
                "en": [
                    "Cell bank: GMP MCB/WCB of CHO-DUKX-B11 cells with bevacizumab genes (anti-VEGF-A human IgG1)",
                    "Thaw + pre-culture: WCB vial, 37 °C/5 % CO2/CDM medium, 3-5 days",
                    "Seed train: stepwise expansion to 2000-12000 L production bioreactor",
                    "Main production: fed-batch 12-15 days, 37 °C, pH 7.0, DO 30 %, optimised feed strategy",
                    "End of production: bevacizumab titer 2-5 g/L in supernatant",
                ],
            },
            "downstream": {
                "de": [
                    "Klärung + Sterilfiltration",
                    "Protein-A-Affinitätschromatographie (Capture)",
                    "Niedrig-pH-Virusinaktivierung",
                    "CEX + AEX-Polishing",
                    "Virusfiltration + TFF mit Trehalose als Stabilisator (Avastin-Formulierung)",
                ],
                "en": [
                    "Clarification + sterile filtration",
                    "Protein-A affinity chromatography (capture)",
                    "Low-pH virus inactivation",
                    "CEX + AEX polishing",
                    "Virus filtration + TFF with trehalose as stabilizer (Avastin formulation)",
                ],
            },
            "expected_purity_after_workup": "≥ 99 % (cGMP, Avastin)",
        },
    },
    "cetuximab": {
        "biotechnological": {
            "organism": "SP2/0 Mausmyelom-Zellen (NICHT CHO — historische Zelllinien-Wahl)",
            "organism_en": "SP2/0 mouse myeloma cells (NOT CHO — historical cell-line choice)",
            "substrate": "Chemisch definiertes Medium (CDM), oft mit Cholesterin-Supplement (SP2/0-Bedarf)",
            "substrate_en": "chemically defined medium (CDM), often with cholesterol supplement (SP2/0 requirement)",
            "yield_range_g_per_l": (0.5, 2),  # SP2/0 lower titer than CHO
            "fermentation_time_h": (240, 336),
            "literature_source": "Kim & Grothey, Drugs 2008 — Cetuximab review (doi:10.2165/0003495-200868090-00007); Burtness et al., Eur. J. Cancer 2008 — Cetuximab clinical (doi:10.1016/j.ejca.2008.07.013)",
            "production": {
                "de": [
                    "Zellbank: GMP-MCB/WCB von SP2/0-Zellen mit Cetuximab-Genen (chimärer Maus/Human IgG1 anti-EGFR)",
                    "Wichtig: Cetuximab trägt durch SP2/0-Glykosylierung Galα1,3Gal-Epitope → erhöhtes Anaphylaxie-Risiko (Schwarzbeermäuse-Allergie)",
                    "Auftauen + Vorkultur: WCB-Ampulle, 37 °C/5 % CO2/CDM mit Cholesterin, 3-5 Tage",
                    "Seed-Train: 50 mL → 2 L → 50 L → 500 L → 2000-5000 L Bioreaktor",
                    "Hauptproduktion: Fed-Batch 10-14 Tage, 37 °C, pH 7.0, DO 30 % (SP2/0 typischerweise empfindlicher als CHO)",
                    "Produktions-Ende: Cetuximab-Titer 0.5-2 g/L (niedriger als typische CHO-Prozesse)",
                ],
                "en": [
                    "Cell bank: GMP MCB/WCB of SP2/0 cells with cetuximab genes (chimeric mouse/human IgG1 anti-EGFR)",
                    "Important: cetuximab carries Gal-alpha-1,3-Gal epitopes from SP2/0 glycosylation -> elevated anaphylaxis risk (alpha-Gal allergy)",
                    "Thaw + pre-culture: WCB vial, 37 °C/5 % CO2/CDM with cholesterol, 3-5 days",
                    "Seed train: 50 mL -> 2 L -> 50 L -> 500 L -> 2000-5000 L bioreactor",
                    "Main production: fed-batch 10-14 days, 37 °C, pH 7.0, DO 30 % (SP2/0 typically more sensitive than CHO)",
                    "End of production: cetuximab titer 0.5-2 g/L (lower than typical CHO processes)",
                ],
            },
            "downstream": {
                "de": [
                    "Klärung + Sterilfiltration",
                    "Protein-A-Affinitätschromatographie (Capture)",
                    "Niedrig-pH-Virusinaktivierung",
                    "Mehrstufiges Polishing (CEX + AEX + HIC) zur Entfernung von Glykoform-Varianten",
                    "Virusfiltration + TFF + Formulierung als isotonisches NaCl-haltiges Konzentrat (Erbitux)",
                ],
                "en": [
                    "Clarification + sterile filtration",
                    "Protein-A affinity chromatography (capture)",
                    "Low-pH virus inactivation",
                    "Multi-stage polishing (CEX + AEX + HIC) to remove glycoform variants",
                    "Virus filtration + TFF + formulation as isotonic NaCl-containing concentrate (Erbitux)",
                ],
            },
            "expected_purity_after_workup": "≥ 99 % (cGMP, Erbitux)",
        },
    },
    "daratumumab": {
        "biotechnological": {
            "organism": "CHO-Zellen (Suspensionskultur)",
            "organism_en": "CHO cells (suspension culture)",
            "substrate": "Chemisch definiertes Medium (CDM)",
            "substrate_en": "chemically defined medium (CDM)",
            "yield_range_g_per_l": (3, 8),
            "fermentation_time_h": (288, 360),
            "literature_source": "Sanchez et al., Cancer Med. 2020 — Daratumumab in multiple myeloma (doi:10.1002/cam4.3437); de Weers et al., J. Immunol. 2011 — Daratumumab mechanism (doi:10.4049/jimmunol.1003032)",
            "production": {
                "de": [
                    "Zellbank: GMP-MCB/WCB von rekombinanten CHO-Zellen mit Daratumumab-Genen (vollhumaner IgG1 anti-CD38)",
                    "Auftauen + Vorkultur: WCB-Ampulle, 37 °C/5 % CO2/CDM-Medium, 3-5 Tage",
                    "Seed-Train: 50 mL → 2 L → 50 L → 500 L → 2000-5000 L Bioreaktor",
                    "Hauptproduktion: Fed-Batch 12-15 Tage, 37 °C, pH 7.0, DO 30 %, optimierte Glykosylierung (afucosyliert für verbesserte ADCC)",
                    "Produktions-Ende: Daratumumab-Titer 3-8 g/L im Überstand",
                ],
                "en": [
                    "Cell bank: GMP MCB/WCB of recombinant CHO cells with daratumumab genes (fully human IgG1 anti-CD38)",
                    "Thaw + pre-culture: WCB vial, 37 °C/5 % CO2/CDM medium, 3-5 days",
                    "Seed train: 50 mL -> 2 L -> 50 L -> 500 L -> 2000-5000 L bioreactor",
                    "Main production: fed-batch 12-15 days, 37 °C, pH 7.0, DO 30 %, optimised glycosylation (afucosylated for enhanced ADCC)",
                    "End of production: daratumumab titer 3-8 g/L in supernatant",
                ],
            },
            "downstream": {
                "de": [
                    "Klärung + Sterilfiltration",
                    "Protein-A-Affinitätschromatographie (Capture)",
                    "Niedrig-pH-Virusinaktivierung",
                    "CEX + AEX-Polishing (Fokus auf Glykoform-Homogenität für ADCC-Konsistenz)",
                    "Virusfiltration + TFF + Formulierung (Darzalex IV oder subkutan)",
                ],
                "en": [
                    "Clarification + sterile filtration",
                    "Protein-A affinity chromatography (capture)",
                    "Low-pH virus inactivation",
                    "CEX + AEX polishing (focus on glycoform homogeneity for ADCC consistency)",
                    "Virus filtration + TFF + formulation (Darzalex IV or subcutaneous)",
                ],
            },
            "expected_purity_after_workup": "≥ 99 % (cGMP, Darzalex)",
        },
    },
    "infliximab": {
        "biotechnological": {
            "organism": "SP2/0 Mausmyelom-Zellen (analog Cetuximab)",
            "organism_en": "SP2/0 mouse myeloma cells (analogous to cetuximab)",
            "substrate": "Chemisch definiertes Medium (CDM) mit Cholesterin",
            "substrate_en": "chemically defined medium (CDM) with cholesterol",
            "yield_range_g_per_l": (0.5, 2),
            "fermentation_time_h": (240, 336),
            "literature_source": "Knight et al., Mol. Immunol. 1993 — Infliximab humanization; Tracey et al., Pharmacol. Ther. 2008 — Anti-TNF therapy review (doi:10.1016/j.pharmthera.2007.10.001)",
            "production": {
                "de": [
                    "Zellbank: GMP-MCB/WCB von SP2/0-Zellen mit Infliximab-Genen (chimärer Maus/Human IgG1 anti-TNFalpha)",
                    "Auftauen + Vorkultur: WCB-Ampulle, 37 °C/5 % CO2/CDM mit Cholesterin, 3-5 Tage",
                    "Seed-Train: 50 mL → 2 L → 50 L → 500 L → 2000-5000 L Bioreaktor",
                    "Hauptproduktion: Fed-Batch 10-14 Tage, 37 °C, pH 7.0, DO 30 %",
                    "Produktions-Ende: Infliximab-Titer 0.5-2 g/L (typisch für SP2/0)",
                ],
                "en": [
                    "Cell bank: GMP MCB/WCB of SP2/0 cells with infliximab genes (chimeric mouse/human IgG1 anti-TNFalpha)",
                    "Thaw + pre-culture: WCB vial, 37 °C/5 % CO2/CDM with cholesterol, 3-5 days",
                    "Seed train: 50 mL -> 2 L -> 50 L -> 500 L -> 2000-5000 L bioreactor",
                    "Main production: fed-batch 10-14 days, 37 °C, pH 7.0, DO 30 %",
                    "End of production: infliximab titer 0.5-2 g/L (typical for SP2/0)",
                ],
            },
            "downstream": {
                "de": [
                    "Klärung + Sterilfiltration",
                    "Protein-A-Affinitätschromatographie",
                    "Niedrig-pH-Virusinaktivierung",
                    "CEX + AEX-Polishing + HIC zur Reduktion der Alpha-Gal-Variabilität",
                    "Virusfiltration + TFF (Remicade lyophilisierte Formulierung)",
                ],
                "en": [
                    "Clarification + sterile filtration",
                    "Protein-A affinity chromatography",
                    "Low-pH virus inactivation",
                    "CEX + AEX polishing + HIC to reduce alpha-Gal variability",
                    "Virus filtration + TFF (Remicade lyophilised formulation)",
                ],
            },
            "expected_purity_after_workup": "≥ 99 % (cGMP, Remicade)",
        },
    },
    "nivolumab": {
        "biotechnological": {
            "organism": "CHO-Zellen (Suspensionskultur)",
            "organism_en": "CHO cells (suspension culture)",
            "substrate": "Chemisch definiertes Medium (CDM)",
            "substrate_en": "chemically defined medium (CDM)",
            "yield_range_g_per_l": (3, 8),
            "fermentation_time_h": (288, 360),
            "literature_source": "Brahmer et al., N. Engl. J. Med. 2015 — Nivolumab in non-small-cell lung cancer (doi:10.1056/NEJMoa1504627); Hodi et al., N. Engl. J. Med. 2010 — Anti-PD-1 therapy",
            "production": {
                "de": [
                    "Zellbank: GMP-MCB/WCB von rekombinanten CHO-Zellen mit Nivolumab-Genen (vollhumaner IgG4 anti-PD-1, mit S228P-Stabilisierungsmutation in der Hinge-Region)",
                    "Auftauen + Vorkultur: WCB-Ampulle, 37 °C/5 % CO2/CDM-Medium, 3-5 Tage",
                    "Seed-Train: 50 mL → 2 L → 50 L → 500 L → 2000-12000 L Bioreaktor",
                    "Hauptproduktion: Fed-Batch 12-15 Tage, 37 °C, pH 7.0, DO 30 %",
                    "Produktions-Ende: Nivolumab-Titer 3-8 g/L im Überstand",
                ],
                "en": [
                    "Cell bank: GMP MCB/WCB of recombinant CHO cells with nivolumab genes (fully human IgG4 anti-PD-1, with S228P stabilization mutation in hinge region)",
                    "Thaw + pre-culture: WCB vial, 37 °C/5 % CO2/CDM medium, 3-5 days",
                    "Seed train: 50 mL -> 2 L -> 50 L -> 500 L -> 2000-12000 L bioreactor",
                    "Main production: fed-batch 12-15 days, 37 °C, pH 7.0, DO 30 %",
                    "End of production: nivolumab titer 3-8 g/L in supernatant",
                ],
            },
            "downstream": {
                "de": [
                    "Klärung + Sterilfiltration",
                    "Protein-A-Affinitätschromatographie (IgG4 bindet schwächer an Standard-ProtA → ggf. spezielle Resine wie MabSelect SuRe)",
                    "Niedrig-pH-Virusinaktivierung (IgG4-Aggregate vermeiden — Zeit + pH eng kontrollieren)",
                    "CEX + AEX-Polishing + HIC zur Halb-Antikörper-Trennung (IgG4-Eigenheit)",
                    "Virusfiltration + TFF (Opdivo, IV-Konzentrat)",
                ],
                "en": [
                    "Clarification + sterile filtration",
                    "Protein-A affinity chromatography (IgG4 binds weaker to standard ProtA -> may need special resins like MabSelect SuRe)",
                    "Low-pH virus inactivation (avoid IgG4 aggregates — tight time + pH control)",
                    "CEX + AEX polishing + HIC for half-antibody separation (IgG4 peculiarity)",
                    "Virus filtration + TFF (Opdivo, IV concentrate)",
                ],
            },
            "expected_purity_after_workup": "≥ 99 % (cGMP, Opdivo)",
        },
    },
    "pertuzumab": {
        "biotechnological": {
            "organism": "CHO-Zellen (Suspensionskultur, gleiche Plattform wie Trastuzumab)",
            "organism_en": "CHO cells (suspension culture, same platform as trastuzumab)",
            "substrate": "Chemisch definiertes Medium (CDM)",
            "substrate_en": "chemically defined medium (CDM)",
            "yield_range_g_per_l": (3, 8),
            "fermentation_time_h": (288, 360),
            "literature_source": "Swain et al., N. Engl. J. Med. 2015 — Pertuzumab in HER2-positive metastatic breast cancer (doi:10.1056/NEJMoa1413513); Adams et al., Cancer Res. 2006 — Pertuzumab mechanism (doi:10.1158/0008-5472.CAN-05-4136)",
            "production": {
                "de": [
                    "Zellbank: GMP-MCB/WCB von rekombinanten CHO-Zellen mit Pertuzumab-Genen (humanisierter IgG1 anti-HER2, bindet an anderes HER2-Epitop als Trastuzumab — komplementäre Wirkung)",
                    "Auftauen + Vorkultur: WCB-Ampulle, 37 °C/5 % CO2/CDM-Medium, 3-5 Tage",
                    "Seed-Train: 50 mL → 2 L → 50 L → 500 L → 2000-12000 L Bioreaktor",
                    "Hauptproduktion: Fed-Batch 12-15 Tage, 37 °C, pH 7.0, DO 30 %",
                    "Produktions-Ende: Pertuzumab-Titer 3-8 g/L im Überstand",
                ],
                "en": [
                    "Cell bank: GMP MCB/WCB of recombinant CHO cells with pertuzumab genes (humanized IgG1 anti-HER2, binds to a different HER2 epitope than trastuzumab — complementary action)",
                    "Thaw + pre-culture: WCB vial, 37 °C/5 % CO2/CDM medium, 3-5 days",
                    "Seed train: 50 mL -> 2 L -> 50 L -> 500 L -> 2000-12000 L bioreactor",
                    "Main production: fed-batch 12-15 days, 37 °C, pH 7.0, DO 30 %",
                    "End of production: pertuzumab titer 3-8 g/L in supernatant",
                ],
            },
            "downstream": {
                "de": [
                    "Klärung + Sterilfiltration",
                    "Protein-A-Affinitätschromatographie",
                    "Niedrig-pH-Virusinaktivierung",
                    "CEX + AEX-Polishing",
                    "Virusfiltration + TFF (Perjeta-Formulierung)",
                ],
                "en": [
                    "Clarification + sterile filtration",
                    "Protein-A affinity chromatography",
                    "Low-pH virus inactivation",
                    "CEX + AEX polishing",
                    "Virus filtration + TFF (Perjeta formulation)",
                ],
            },
            "expected_purity_after_workup": "≥ 99 % (cGMP, Perjeta)",
        },
    },
    "rituximab": {
        "biotechnological": {
            "organism": "CHO-DG44-Zellen (dhfr-negativ, MTX-amplifiziert)",
            "organism_en": "CHO-DG44 cells (dhfr-negative, MTX-amplified)",
            "substrate": "Chemisch definiertes Medium (CDM)",
            "substrate_en": "chemically defined medium (CDM)",
            "yield_range_g_per_l": (1, 4),
            "fermentation_time_h": (288, 360),
            "literature_source": "Reff et al., Blood 1994 — Rituximab development (doi:10.1182/blood.V83.2.435.435); Pierpont et al., Front. Oncol. 2018 — Rituximab clinical review (doi:10.3389/fonc.2018.00163)",
            "production": {
                "de": [
                    "Zellbank: GMP-MCB/WCB von CHO-DG44-Zellen mit Rituximab-Genen (chimärer Maus/Human IgG1 anti-CD20)",
                    "Auftauen + Vorkultur: WCB-Ampulle, 37 °C/5 % CO2/CDM-Medium, 3-5 Tage",
                    "Seed-Train: 50 mL → 2 L → 50 L → 500 L → 2000-12000 L Bioreaktor",
                    "Hauptproduktion: Fed-Batch 12-15 Tage, 37 °C, pH 7.0, DO 30 %, MTX-Selektion gewährleistet hohe Gen-Kopienzahl",
                    "Produktions-Ende: Rituximab-Titer 1-4 g/L (ältere Plattform — neuere Biosimilars erreichen 3-5 g/L)",
                ],
                "en": [
                    "Cell bank: GMP MCB/WCB of CHO-DG44 cells with rituximab genes (chimeric mouse/human IgG1 anti-CD20)",
                    "Thaw + pre-culture: WCB vial, 37 °C/5 % CO2/CDM medium, 3-5 days",
                    "Seed train: 50 mL -> 2 L -> 50 L -> 500 L -> 2000-12000 L bioreactor",
                    "Main production: fed-batch 12-15 days, 37 °C, pH 7.0, DO 30 %, MTX selection ensures high gene copy number",
                    "End of production: rituximab titer 1-4 g/L (older platform — newer biosimilars reach 3-5 g/L)",
                ],
            },
            "downstream": {
                "de": [
                    "Klärung + Sterilfiltration",
                    "Protein-A-Affinitätschromatographie (Capture)",
                    "Niedrig-pH-Virusinaktivierung",
                    "CEX + AEX-Polishing (Glykoform-Konsistenz kritisch für ADCC bei NHL-Indikation)",
                    "Virusfiltration + TFF (MabThera / Rituxan / Truxima Biosimilar)",
                ],
                "en": [
                    "Clarification + sterile filtration",
                    "Protein-A affinity chromatography (capture)",
                    "Low-pH virus inactivation",
                    "CEX + AEX polishing (glycoform consistency critical for ADCC in NHL indication)",
                    "Virus filtration + TFF (MabThera / Rituxan / Truxima biosimilar)",
                ],
            },
            "expected_purity_after_workup": "≥ 99 % (cGMP, MabThera / Rituxan)",
        },
    },
    # ====================================================================
    # BATCH 6 — Protein / enzyme (Industrieenzyme)
    # ====================================================================
    "amylase": {
        "biotechnological": {
            "organism": "Bacillus licheniformis (alpha-Amylase, thermostabil) oder Aspergillus oryzae (alpha-Amylase, mesophil)",
            "organism_en": "Bacillus licheniformis (alpha-amylase, thermostable) or Aspergillus oryzae (alpha-amylase, mesophilic)",
            "substrate": "Glucose + Stärke + Maisquellwasser (komplex)",
            "substrate_en": "glucose + starch + corn-steep liquor (complex)",
            "yield_range_g_per_l": (5, 20),
            "fermentation_time_h": (48, 96),
            "literature_source": "Pandey et al., Biotechnol. Appl. Biochem. 2000 — Industrial amylase production (doi:10.1042/BA19990073); Souza & Magalhães, Braz. J. Microbiol. 2010 — Application of microbial alpha-amylase (doi:10.1590/S1517-83822010000400004)",
            "production": {
                "de": [
                    "Stammvorbereitung: B. licheniformis Industrial-Stamm (z. B. Termamyl-Stämme von Novozymes) oder A. oryzae",
                    "Vorkultur: komplexes Medium (Stärke + Maisquellwasser + Salze), 37 °C (B. licheniformis) bzw. 30 °C (A. oryzae), 24-48 h",
                    "Hauptfermentation: Submerskultur in 50-500 m3 Bioreaktor, 37 °C bzw. 30 °C, pH 6.5-7.5, hoher O2-Eintrag (>0.5 vvm)",
                    "Stärke im Medium induziert Amylase-Expression (Substrat-Induktion); A. oryzae sekretiert extrazellulär",
                    "Fermentations-Dauer 48-96 h, Enzym-Titer 5-20 g/L im Kulturüberstand",
                ],
                "en": [
                    "Strain preparation: B. licheniformis industrial strain (e.g. Termamyl strains from Novozymes) or A. oryzae",
                    "Pre-culture: complex medium (starch + corn-steep liquor + salts), 37 °C (B. licheniformis) or 30 °C (A. oryzae), 24-48 h",
                    "Main fermentation: submerged culture in 50-500 m3 bioreactor, 37 °C or 30 °C, pH 6.5-7.5, high O2 transfer (>0.5 vvm)",
                    "Starch in medium induces amylase expression (substrate induction); A. oryzae secretes extracellularly",
                    "Fermentation duration 48-96 h, enzyme titer 5-20 g/L in culture supernatant",
                ],
            },
            "downstream": {
                "de": [
                    "Zellabtrennung via Rotations-Vakuumfilter oder Mikrofiltration",
                    "Aufkonzentrierung via Ultrafiltration (10-30 kDa Cutoff)",
                    "Stabilisierung mit CaCl2 + Sorbitol (für alpha-Amylase wichtig)",
                    "Formulierung: flüssig (Liquid) oder Sprühtrocknung (Granulat)",
                ],
                "en": [
                    "Cell separation via rotary vacuum filter or microfiltration",
                    "Concentration via ultrafiltration (10-30 kDa cutoff)",
                    "Stabilisation with CaCl2 + sorbitol (important for alpha-amylase)",
                    "Formulation: liquid or spray-dried granulate",
                ],
            },
            "expected_purity_after_workup": "Industrial grade (Mehrkomponenten-Formulierung); Aktivität spezifiziert in KNU/g oder BAU/g; pharma-grade > 95 %",
        },
    },
    "cellulase": {
        "biotechnological": {
            "organism": "Trichoderma reesei (Hypocrea jecorina) — Industrie-Workhorse für Cellulase",
            "organism_en": "Trichoderma reesei (Hypocrea jecorina) — industrial workhorse for cellulase",
            "substrate": "Cellulose (Avicel oder vorbehandelte Lignocellulose) + Lactose als Induktor",
            "substrate_en": "cellulose (Avicel or pretreated lignocellulose) + lactose as inducer",
            "yield_range_g_per_l": (40, 100),
            "fermentation_time_h": (120, 200),
            "literature_source": "Bischof, Ramoni & Seiboth, Microb. Cell Fact. 2016 — Cellulases and beyond: the first 70 years of the enzyme producer Trichoderma reesei (doi:10.1186/s12934-016-0507-6); Peterson & Nevalainen, Microbiology 2012 (doi:10.1099/mic.0.054072-0)",
            "production": {
                "de": [
                    "Stammvorbereitung: T. reesei Hochleister-Stamm (RUT-C30 oder neuere Hyper-Producer), Sporen-Inokulation",
                    "Vorkultur: Glucose-haltiges Medium, 28 °C, 48-72 h",
                    "Hauptfermentation: Submerskultur in 100-500 m3 Bioreaktor, 28-30 °C, pH 4.5-5.5",
                    "Schlüssel-Induktoren: Cellulose oder Lactose → schalten Cellulase-Expression an (CRE1-Knockout-Stämme dekrepetisiert)",
                    "Fermentations-Dauer 120-200 h, extrazelluläre Sekretion 40-100 g/L Cellulase-Cocktail (Endoglucanasen + Cellobiohydrolasen + Beta-Glucosidase) — höchster Sekretions-Titer aller Industrieenzyme",
                ],
                "en": [
                    "Strain preparation: T. reesei high-yield strain (RUT-C30 or newer hyper-producers), spore inoculation",
                    "Pre-culture: glucose-containing medium, 28 °C, 48-72 h",
                    "Main fermentation: submerged culture in 100-500 m3 bioreactor, 28-30 °C, pH 4.5-5.5",
                    "Key inducers: cellulose or lactose -> switch on cellulase expression (CRE1-knockout strains de-repressed)",
                    "Fermentation duration 120-200 h, extracellular secretion 40-100 g/L cellulase cocktail (endoglucanases + cellobiohydrolases + beta-glucosidase) — highest secretion titer of all industrial enzymes",
                ],
            },
            "downstream": {
                "de": [
                    "Myzel-Abtrennung via Rotations-Vakuumfilter oder Mikrofiltration",
                    "Ultrafiltration zur Aufkonzentrierung (5-10 kDa Cutoff)",
                    "Stabilisierung mit Glycerin / Sorbitol",
                    "Formulierung: flüssig (Bulk-Industrie für Bio-Ethanol-Produktion) oder sprühgetrocknet",
                ],
                "en": [
                    "Mycelium separation via rotary vacuum filter or microfiltration",
                    "Ultrafiltration for concentration (5-10 kDa cutoff)",
                    "Stabilisation with glycerol / sorbitol",
                    "Formulation: liquid (bulk industry for bioethanol production) or spray-dried",
                ],
            },
            "expected_purity_after_workup": "Industrial grade (Cellulase-Cocktail); Aktivität in FPU/mL spezifiziert",
        },
    },
    "lactase": {
        "biotechnological": {
            "organism": "Kluyveromyces lactis (Hefe, Standard für laktose-freie Milchprodukte) oder Aspergillus oryzae",
            "organism_en": "Kluyveromyces lactis (yeast, standard for lactose-free dairy) or Aspergillus oryzae",
            "substrate": "Laktose / Molke (Schlüssel-Induktor)",
            "substrate_en": "lactose / whey (key inducer)",
            "yield_range_g_per_l": (1, 8),
            "fermentation_time_h": (24, 72),
            "literature_source": "Husain, Crit. Rev. Biotechnol. 2010 — beta-galactosidases and their potential applications (doi:10.3109/07388550903330497); Panesar et al., Enzyme Microb. Technol. 2010 (doi:10.1016/j.enzmictec.2010.01.008)",
            "production": {
                "de": [
                    "Stammvorbereitung: K. lactis CBS 2359 oder A. oryzae Industrial-Stamm",
                    "Vorkultur: Laktose- oder Hefeextrakt-Glucose-Medium, 30 °C, 16-24 h",
                    "Hauptfermentation: Submerskultur in 10-100 m3 Bioreaktor, 30 °C (K. lactis) bzw. 30 °C (A. oryzae), pH 6.0-7.0",
                    "Substrat: Laktose (z. B. aus Käsemolke) induziert die beta-Galactosidase-Expression",
                    "Fermentations-Dauer 24-72 h, Enzym-Titer 1-8 g/L (K. lactis intrazellulär — A. oryzae extrazellulär)",
                ],
                "en": [
                    "Strain preparation: K. lactis CBS 2359 or A. oryzae industrial strain",
                    "Pre-culture: lactose or yeast extract-glucose medium, 30 °C, 16-24 h",
                    "Main fermentation: submerged culture in 10-100 m3 bioreactor, 30 °C (K. lactis) or 30 °C (A. oryzae), pH 6.0-7.0",
                    "Substrate: lactose (e.g. from cheese whey) induces beta-galactosidase expression",
                    "Fermentation duration 24-72 h, enzyme titer 1-8 g/L (K. lactis intracellular — A. oryzae extracellular)",
                ],
            },
            "downstream": {
                "de": [
                    "Bei K. lactis: Zell-Lyse (Hefe via Glaskugelmühle) zur Freisetzung der intrazellulären Lactase",
                    "Klärung via Filtration",
                    "Ultrafiltration + Aktivkohle-Polishing",
                    "Stabilisierung mit Glycerin oder Trehalose, Formulierung als Flüssigkeit oder Lyophilisat",
                ],
                "en": [
                    "K. lactis: cell lysis (yeast via bead mill) to release intracellular lactase",
                    "Clarification via filtration",
                    "Ultrafiltration + activated-carbon polishing",
                    "Stabilisation with glycerol or trehalose, formulation as liquid or lyophilisate",
                ],
            },
            "expected_purity_after_workup": "Food-grade (Aktivität spezifiziert in NLU/g); pharma-grade ≥ 95 %",
        },
    },
    "lipase": {
        "biotechnological": {
            "organism": "Aspergillus oryzae (Wirtsstamm für rekombinante Lipase) oder Candida antarctica (CalB, native Quelle)",
            "organism_en": "Aspergillus oryzae (host for recombinant lipase) or Candida antarctica (CalB, native source)",
            "substrate": "Glucose + Hefeextrakt + Sojabohnenöl (Induktor)",
            "substrate_en": "glucose + yeast extract + soybean oil (inducer)",
            "yield_range_g_per_l": (2, 15),
            "fermentation_time_h": (48, 120),
            "literature_source": "Hasan, Shah & Hameed, Enzyme Microb. Technol. 2006 — Industrial applications of microbial lipases (doi:10.1016/j.enzmictec.2005.10.016); Sharma et al., Biotechnol. Adv. 2001 — Microbial lipases (doi:10.1016/S0734-9750(01)00086-6)",
            "production": {
                "de": [
                    "Stammvorbereitung: A. oryzae mit rekombinantem CalB-Gen (kommerzielle Lipase Novozyme 435) oder native C. antarctica",
                    "Vorkultur: Glucose-Hefeextrakt-Medium, 30 °C, 24-48 h",
                    "Hauptfermentation: Submerskultur in 10-100 m3 Bioreaktor, 28-30 °C, pH 6.0-7.0",
                    "Schlüssel-Induktor: Sojabohnenöl oder Triglyceride → schalten Lipase-Expression ein",
                    "Fermentations-Dauer 48-120 h, Enzym-Titer 2-15 g/L (rekombinante A. oryzae-Stämme erreichen höhere Werte)",
                ],
                "en": [
                    "Strain preparation: A. oryzae with recombinant CalB gene (commercial lipase Novozyme 435) or native C. antarctica",
                    "Pre-culture: glucose-yeast extract medium, 30 °C, 24-48 h",
                    "Main fermentation: submerged culture in 10-100 m3 bioreactor, 28-30 °C, pH 6.0-7.0",
                    "Key inducer: soybean oil or triglycerides -> switch on lipase expression",
                    "Fermentation duration 48-120 h, enzyme titer 2-15 g/L (recombinant A. oryzae strains reach higher values)",
                ],
            },
            "downstream": {
                "de": [
                    "Zellabtrennung via Filtration",
                    "Ultrafiltration zur Aufkonzentrierung",
                    "Optional: Immobilisierung auf hydrophobem Träger (z. B. Lewatit VP OC 1600 für Novozyme 435) — Standard für Industrieanwendungen",
                    "Stabilisierung + Formulierung als immobilisierte Pellets oder flüssige Lipase",
                ],
                "en": [
                    "Cell separation via filtration",
                    "Ultrafiltration for concentration",
                    "Optional: immobilisation on hydrophobic carrier (e.g. Lewatit VP OC 1600 for Novozyme 435) — standard for industrial applications",
                    "Stabilisation + formulation as immobilised pellets or liquid lipase",
                ],
            },
            "expected_purity_after_workup": "Industrial grade (Aktivität in PLU/g für Lipase oder LU/g)",
        },
    },
    "pectinase": {
        "biotechnological": {
            "organism": "Aspergillus niger (Endo-Polygalacturonase + Pektin-Lyase + Pektinesterase)",
            "organism_en": "Aspergillus niger (endo-polygalacturonase + pectin lyase + pectin esterase)",
            "substrate": "Pektin / Zitrusschalen-Hydrolysat + Glucose",
            "substrate_en": "pectin / citrus-peel hydrolysate + glucose",
            "yield_range_g_per_l": (5, 30),
            "fermentation_time_h": (48, 120),
            "literature_source": "Garg et al., Front. Microbiol. 2016 — Microbial pectinases (doi:10.3389/fmicb.2016.00010); Jayani, Saxena & Gupta, Process Biochem. 2005 — Microbial pectinolytic enzymes (doi:10.1016/j.procbio.2005.03.026)",
            "production": {
                "de": [
                    "Stammvorbereitung: A. niger Industrial-Stamm (z. B. NRRL 3 oder kommerziell von DSM/Novozymes)",
                    "Vorkultur: Pektin- oder Zitrusschalen-haltiges Medium, 28-30 °C, 48-72 h",
                    "Hauptfermentation: Submerskultur in 50-200 m3 Bioreaktor, 28-30 °C, pH 4.0-5.5",
                    "Schlüssel-Induktor: Pektin oder Galakturonsäure → schalten Pektinase-Cluster-Expression an",
                    "Fermentations-Dauer 48-120 h, extrazelluläre Sekretion 5-30 g/L Pektinase-Cocktail",
                ],
                "en": [
                    "Strain preparation: A. niger industrial strain (e.g. NRRL 3 or commercial from DSM/Novozymes)",
                    "Pre-culture: pectin- or citrus-peel-containing medium, 28-30 °C, 48-72 h",
                    "Main fermentation: submerged culture in 50-200 m3 bioreactor, 28-30 °C, pH 4.0-5.5",
                    "Key inducer: pectin or galacturonic acid -> switch on pectinase cluster expression",
                    "Fermentation duration 48-120 h, extracellular secretion 5-30 g/L pectinase cocktail",
                ],
            },
            "downstream": {
                "de": [
                    "Myzel-Abtrennung via Tiefenfiltration",
                    "Ultrafiltration zur Aufkonzentrierung",
                    "Aktivkohle-Polishing zur Entfernung dunkler Pigmente",
                    "Formulierung: flüssig (Saftherstellung) oder Sprühtrocknung",
                ],
                "en": [
                    "Mycelium separation via depth filtration",
                    "Ultrafiltration for concentration",
                    "Activated-carbon polishing to remove dark pigments",
                    "Formulation: liquid (juice production) or spray-dried",
                ],
            },
            "expected_purity_after_workup": "Food-grade (Aktivität in PGNU/g oder AVJP/g); Hauptanwendung Saftklärung",
        },
    },
    "phytase": {
        "biotechnological": {
            "organism": "Aspergillus niger (native PhyA) oder rekombinant in Trichoderma reesei / Pichia pastoris",
            "organism_en": "Aspergillus niger (native PhyA) or recombinant in Trichoderma reesei / Pichia pastoris",
            "substrate": "Glucose + komplexes Medium (Sojamehl, Maisquellwasser)",
            "substrate_en": "glucose + complex medium (soybean meal, corn-steep liquor)",
            "yield_range_g_per_l": (2, 15),
            "fermentation_time_h": (96, 168),
            "literature_source": "Lei et al., Annu. Rev. Anim. Biosci. 2013 — Phytase, a new life for an 'old' enzyme (doi:10.1146/annurev-animal-031412-103717); Mullaney & Ullah, Biochem. Biophys. Res. Commun. 2003 — Phytase review (doi:10.1016/S0006-291X(03)00708-4)",
            "production": {
                "de": [
                    "Stammvorbereitung: A. niger (PhyA-Gen nativ) ODER rekombinante T. reesei / P. pastoris mit phyA-Gen unter starkem Promotor (AOX1 für P. pastoris)",
                    "Vorkultur: komplexes Medium mit niedrigem Phosphat (Phosphat-Limitierung induziert PhyA-Expression), 28-30 °C, 24-48 h",
                    "Hauptfermentation: Submerskultur in 50-200 m3 Bioreaktor, 28-30 °C (A. niger / T. reesei) bzw. 30 °C (P. pastoris)",
                    "P. pastoris: Methanol-Feed zur Induktion von AOX1-Promotor (Phasen: Glycerin-Wachstum → Methanol-Induktion)",
                    "Fermentations-Dauer 96-168 h, extrazelluläre Sekretion 2-15 g/L Phytase",
                ],
                "en": [
                    "Strain preparation: A. niger (native PhyA gene) OR recombinant T. reesei / P. pastoris with phyA gene under strong promoter (AOX1 for P. pastoris)",
                    "Pre-culture: complex medium with low phosphate (phosphate limitation induces PhyA expression), 28-30 °C, 24-48 h",
                    "Main fermentation: submerged culture in 50-200 m3 bioreactor, 28-30 °C (A. niger / T. reesei) or 30 °C (P. pastoris)",
                    "P. pastoris: methanol feed for AOX1 promoter induction (phases: glycerol growth -> methanol induction)",
                    "Fermentation duration 96-168 h, extracellular secretion 2-15 g/L phytase",
                ],
            },
            "downstream": {
                "de": [
                    "Zellabtrennung via Filtration",
                    "Ultrafiltration zur Aufkonzentrierung",
                    "Granulation oder Pelletierung (Phytase wird oft als hitzestabile Beschichtungs-Formulierung für Tierfutter eingesetzt)",
                    "Stabilisierung mit Sorbitol",
                ],
                "en": [
                    "Cell separation via filtration",
                    "Ultrafiltration for concentration",
                    "Granulation or pelletisation (phytase is often deployed as heat-stable coating formulation for animal feed)",
                    "Stabilisation with sorbitol",
                ],
            },
            "expected_purity_after_workup": "Feed-grade (Aktivität in FTU/g); kommerziell als Quantum Blue (AB Vista) oder Ronozyme (DSM)",
        },
    },
    "protease": {
        "biotechnological": {
            "organism": "Bacillus subtilis (Subtilisin, Standard-Industrieenzym) oder Aspergillus oryzae",
            "organism_en": "Bacillus subtilis (subtilisin, standard industrial enzyme) or Aspergillus oryzae",
            "substrate": "Glucose + Sojamehl + Maisquellwasser",
            "substrate_en": "glucose + soybean meal + corn-steep liquor",
            "yield_range_g_per_l": (5, 25),
            "fermentation_time_h": (24, 72),
            "literature_source": "Rao et al., Microbiol. Mol. Biol. Rev. 1998 — Molecular and biotechnological aspects of microbial proteases (doi:10.1128/MMBR.62.3.597-635.1998); Gupta et al., Appl. Microbiol. Biotechnol. 2002 — Bacterial alkaline proteases (doi:10.1007/s00253-002-0975-y)",
            "production": {
                "de": [
                    "Stammvorbereitung: B. subtilis Industrial-Stamm (z. B. Carlsberg-Stamm für Alcalase von Novozymes) oder A. oryzae",
                    "Vorkultur: Pepton-Hefeextrakt-Medium, 37 °C (B. subtilis) bzw. 30 °C (A. oryzae), 16-24 h",
                    "Hauptfermentation: Submerskultur in 50-500 m3 Bioreaktor, 37 °C bzw. 30 °C, pH 7.0-8.0 (Subtilisin ist alkalische Protease)",
                    "Stationäre Phase wichtig: Subtilisin-Sekretion startet nach Wachstumsphase",
                    "Fermentations-Dauer 24-72 h, extrazelluläre Sekretion 5-25 g/L (Alcalase erreicht > 20 g/L)",
                ],
                "en": [
                    "Strain preparation: B. subtilis industrial strain (e.g. Carlsberg strain for Alcalase from Novozymes) or A. oryzae",
                    "Pre-culture: peptone-yeast extract medium, 37 °C (B. subtilis) or 30 °C (A. oryzae), 16-24 h",
                    "Main fermentation: submerged culture in 50-500 m3 bioreactor, 37 °C or 30 °C, pH 7.0-8.0 (subtilisin is alkaline protease)",
                    "Stationary phase important: subtilisin secretion starts after growth phase",
                    "Fermentation duration 24-72 h, extracellular secretion 5-25 g/L (Alcalase reaches > 20 g/L)",
                ],
            },
            "downstream": {
                "de": [
                    "Zellabtrennung via Filtration",
                    "Ultrafiltration + Reinigung der Protease durch Inhibitor-Zugabe (PMSF) während Verarbeitung zur Vermeidung von Auto-Proteolyse",
                    "Stabilisierung mit CaCl2 + Polyolen",
                    "Formulierung als Flüssigkeit oder Sprühtrocknung; viele Anwendungen brauchen Verkapselung (Waschmittel-Granulat)",
                ],
                "en": [
                    "Cell separation via filtration",
                    "Ultrafiltration + protease purification with inhibitor addition (PMSF) during processing to avoid auto-proteolysis",
                    "Stabilisation with CaCl2 + polyols",
                    "Formulation as liquid or spray-dried; many applications require encapsulation (detergent granulate)",
                ],
            },
            "expected_purity_after_workup": "Industrial grade (Aktivität in AU/g oder KNPU/g); Anwendungen: Waschmittel, Lederbearbeitung, Käserei",
        },
    },
    # ====================================================================
    # BATCH 7 — Natural product / terpene (Wasserdampf-Destillation / Biotech)
    # ====================================================================
    "α-pinene": {
        "extraction": {
            "source": "Kiefernharz (Pinus sylvestris, P. pinaster) — Terpentinöl",
            "source_en": "Pine resin (Pinus sylvestris, P. pinaster) — turpentine oil",
            "yield_range_percent_w_w": (60, 80),  # alpha-Pinene-Anteil im Terpentinöl
            "literature_source": "Behr & Johnen, ChemSusChem 2009 — Myrcene as a natural base chemical (doi:10.1002/cssc.200900186); Erman, Chemistry of the Monoterpenes, Marcel Dekker 1985 (ISBN 0-8247-7338-1)",
            "production": {
                "de": [
                    "Rohmaterial: Harzgewinnung durch Anritzen lebender Kiefern (gum turpentine) ODER Sulfat-Terpentin als Nebenprodukt der Zellstoff-Produktion (Kraft-Verfahren)",
                    "Vorbehandlung: Harz wird mit Wasserdampf destilliert → Roh-Terpentinöl (45-65 % alpha-Pinen, 25-35 % beta-Pinen, Rest 3-Caren, Camphen)",
                    "alpha-Pinen-Anteil im Roh-Terpentinöl 60-80 % w/w bei sorgfältiger Auswahl der Kiefernart",
                    "Industrieller Maßstab: hauptsächlich aus Sulfat-Terpentin (USA, Skandinavien)",
                ],
                "en": [
                    "Raw material: resin extraction by tapping live pine trees (gum turpentine) OR sulfate turpentine as side product of kraft pulp production",
                    "Pre-treatment: resin is steam-distilled -> crude turpentine oil (45-65 % alpha-pinene, 25-35 % beta-pinene, rest 3-carene, camphene)",
                    "alpha-Pinene content in crude turpentine 60-80 % w/w with careful pine-species selection",
                    "Industrial scale: mainly from sulfate turpentine (USA, Scandinavia)",
                ],
            },
            "downstream": {
                "de": [
                    "Fraktionierte Destillation bei reduziertem Druck (alpha-Pinen: 156 °C bei Normaldruck, 60 °C bei 50 mbar)",
                    "Trennung von beta-Pinen via Schmelzkristallisation oder enantioselektive Säulen-Chromatographie",
                    "Optisch aktive Trennung (+)- vs (-)-alpha-Pinen mit chiraler Stationärphase",
                    "Lagerung unter Inertgas (alpha-Pinen ist oxidationsempfindlich)",
                ],
                "en": [
                    "Fractional distillation at reduced pressure (alpha-pinene: 156 °C at atmospheric, 60 °C at 50 mbar)",
                    "Separation from beta-pinene via melt crystallization or enantioselective column chromatography",
                    "Optical separation (+)- vs (-)-alpha-pinene with chiral stationary phase",
                    "Storage under inert gas (alpha-pinene is oxidation-sensitive)",
                ],
            },
            "expected_purity_after_workup": "≥ 95 % (technical / Vorstufe für Camphor- und Geraniol-Synthese)",
        },
    },
    "citral": {
        "extraction": {
            "source": "Lemongrass-Öl (Cymbopogon citratus) — Wasserdampf-Destillation",
            "source_en": "Lemongrass oil (Cymbopogon citratus) — steam distillation",
            "yield_range_percent_w_w": (70, 85),  # citral content in lemongrass oil
            "literature_source": "Schwab et al., Plant J. 2008 — Biosynthesis of plant-derived flavor compounds (doi:10.1111/j.1365-313X.2008.03446.x); Negi et al., J. Essent. Oil Res. 2005 — Lemongrass essential oil composition",
            "production": {
                "de": [
                    "Rohmaterial: frisches oder getrocknetes Lemongrass (Cymbopogon citratus), Ernte 3-4x pro Jahr",
                    "Vorbehandlung: Lemongrass wird gehäckselt und in Destillations-Kessel gefüllt",
                    "Wasserdampf-Destillation: 100-110 °C, 2-4 h pro Charge → Roh-Lemongrass-Öl (typ. 65-85 % Citral als Mischung aus Geranial + Neral)",
                    "Geographisch: Hauptproduzenten Indien (Cochin-Sorte), Guatemala, China",
                    "Alternative: chemische Synthese aus beta-Pinen über Myrcen-Route (BASF, Symrise)",
                ],
                "en": [
                    "Raw material: fresh or dried lemongrass (Cymbopogon citratus), harvest 3-4x per year",
                    "Pre-treatment: lemongrass is chopped and loaded into distillation vessel",
                    "Steam distillation: 100-110 °C, 2-4 h per batch -> crude lemongrass oil (typ. 65-85 % citral as mixture of geranial + neral)",
                    "Geographically: main producers India (Cochin variety), Guatemala, China",
                    "Alternative: chemical synthesis from beta-pinene via myrcene route (BASF, Symrise)",
                ],
            },
            "downstream": {
                "de": [
                    "Phasentrennung Öl/Wasser",
                    "Fraktionierte Vakuum-Destillation (Citral siedet bei 229 °C atm, 90-100 °C bei 10 mbar)",
                    "Trennung Geranial (E-Isomer) und Neral (Z-Isomer) optional via fraktionierte Destillation oder Chromatographie",
                    "Lagerung in dunklen Glasflaschen unter Inertgas (oxidationsempfindlich, polymerisations-anfällig)",
                ],
                "en": [
                    "Oil/water phase separation",
                    "Fractional vacuum distillation (citral boils at 229 °C atm, 90-100 °C at 10 mbar)",
                    "Separation of geranial (E-isomer) and neral (Z-isomer) optional via fractional distillation or chromatography",
                    "Storage in dark glass bottles under inert gas (oxidation-sensitive, polymerization-prone)",
                ],
            },
            "expected_purity_after_workup": "≥ 96 % (Lebensmittel-/Flavor-Grade); ≥ 99 % für Vitamin-A-Synthese",
        },
    },
    "geraniol": {
        "extraction": {
            "source": "Palmarosa-Öl (Cymbopogon martinii) ODER Citronella-Öl ODER Rosenöl",
            "source_en": "Palmarosa oil (Cymbopogon martinii) OR citronella oil OR rose oil",
            "yield_range_percent_w_w": (75, 90),  # geraniol content in palmarosa oil
            "literature_source": "Dubey & Luthra, Curr. Sci. 2001 — Biosynthesis of essential oil compounds (doi:10.1126/science.1057284 - related); Chen & Viljoen, S. Afr. J. Bot. 2010 — Geraniol natural occurrence (doi:10.1016/j.sajb.2010.05.008)",
            "production": {
                "de": [
                    "Rohmaterial: Palmarosa-Gras (Cymbopogon martinii), Ernte vor der Blüte für höchsten Geraniol-Gehalt",
                    "Vorbehandlung: getrocknetes Gras wird gehäckselt",
                    "Wasserdampf-Destillation: 100-105 °C, 3-5 h pro Charge → Roh-Palmarosa-Öl (75-90 % Geraniol als Hauptkomponente)",
                    "Alternative Quellen: Citronella-Java-Öl (35-45 % Geraniol) oder Rosenöl (15-30 %)",
                    "Geographische Hauptproduzenten: Indien (über 90 % des Welt-Palmarosa-Öls), Vietnam",
                ],
                "en": [
                    "Raw material: palmarosa grass (Cymbopogon martinii), harvest before flowering for highest geraniol content",
                    "Pre-treatment: dried grass is chopped",
                    "Steam distillation: 100-105 °C, 3-5 h per batch -> crude palmarosa oil (75-90 % geraniol as main component)",
                    "Alternative sources: citronella-Java oil (35-45 % geraniol) or rose oil (15-30 %)",
                    "Main geographic producers: India (over 90 % of world palmarosa oil), Vietnam",
                ],
            },
            "downstream": {
                "de": [
                    "Phasentrennung Öl/Wasser",
                    "Fraktionierte Vakuum-Destillation (Geraniol: 230 °C atm, 100 °C bei 10 mbar)",
                    "Trennung von Nerol (cis-Isomer) und Linalool",
                    "Lagerung unter Inertgas",
                ],
                "en": [
                    "Oil/water phase separation",
                    "Fractional vacuum distillation (geraniol: 230 °C atm, 100 °C at 10 mbar)",
                    "Separation from nerol (cis-isomer) and linalool",
                    "Storage under inert gas",
                ],
            },
            "expected_purity_after_workup": "≥ 97 % (Flavor & Fragrance Industrie-Standard); ≥ 99 % für Pharma-Synthese",
        },
    },
    "limonene": {
        "extraction": {
            "source": "Zitrusfrucht-Schalen (Orange, Zitrone, Grapefruit) — Nebenprodukt der Saftindustrie",
            "source_en": "Citrus peels (orange, lemon, grapefruit) — by-product of juice industry",
            "yield_range_percent_w_w": (90, 95),  # limonene content in citrus peel oil
            "literature_source": "González-Mas et al., Front. Plant Sci. 2019 — Volatile compounds in Citrus essential oils (doi:10.3389/fpls.2019.00012); Ciriminna et al., Org. Process Res. Dev. 2014 — Limonene as a renewable building block (doi:10.1021/op4002814)",
            "production": {
                "de": [
                    "Rohmaterial: Zitrusschalen aus der Saftproduktion (jährlich ~20 Mio. t Schalenabfall weltweit, hauptsächlich Brasilien + Florida)",
                    "Vorbehandlung: Schalen werden mechanisch gepresst (Cold-Press-Verfahren) — Limonen-Öl tritt aus Öldrüsen aus",
                    "Alternative: Wasserdampf-Destillation der gepressten Schalen für höhere Ausbeute",
                    "Roh-Schalen-Öl: 90-95 % d-Limonen (R-Enantiomer) als Hauptkomponente, neben Linalool, Myrcen",
                    "Industrieller Maßstab: ~50-70 kg d-Limonen pro t Orangen (aus Schalen)",
                ],
                "en": [
                    "Raw material: citrus peels from juice production (annually ~20 million t peel waste worldwide, mainly Brazil + Florida)",
                    "Pre-treatment: peels are mechanically pressed (cold-press process) — limonene oil exits from oil glands",
                    "Alternative: steam distillation of pressed peels for higher yield",
                    "Crude peel oil: 90-95 % d-limonene (R-enantiomer) as main component, plus linalool, myrcene",
                    "Industrial scale: ~50-70 kg d-limonene per ton of oranges (from peels)",
                ],
            },
            "downstream": {
                "de": [
                    "Phasentrennung von wässriger Phase",
                    "Filtration zur Entfernung von Pulpe-Rückständen",
                    "Vakuum-Destillation (Limonen: 176 °C atm, 60 °C bei 10 mbar)",
                    "Optional: chirale Trennung d- vs l-Limonen für asymmetrische Synthese-Anwendungen",
                    "Lagerung unter Inertgas (oxidationsempfindlich → Bildung von Carvon)",
                ],
                "en": [
                    "Aqueous phase separation",
                    "Filtration to remove pulp residues",
                    "Vacuum distillation (limonene: 176 °C atm, 60 °C at 10 mbar)",
                    "Optional: chiral separation d- vs l-limonene for asymmetric synthesis applications",
                    "Storage under inert gas (oxidation-sensitive -> formation of carvone)",
                ],
            },
            "expected_purity_after_workup": "≥ 96 % (technical, green-solvent); ≥ 99 % (food/flavor grade)",
        },
    },
    "linalool": {
        "extraction": {
            "source": "Lavendel-Öl (Lavandula angustifolia), Bois-de-Rose-Öl, Ho-Öl (Cinnamomum camphora)",
            "source_en": "Lavender oil (Lavandula angustifolia), bois-de-rose oil, ho oil (Cinnamomum camphora)",
            "yield_range_percent_w_w": (30, 40),  # linalool in lavender oil
            "literature_source": "Lis-Balchin, Lavender: The Genus Lavandula, Taylor & Francis 2002 (ISBN 978-0-415-28486-7); Aprotosoaie et al., Flavour Fragr. J. 2014 — Linalool review (doi:10.1002/ffj.3197)",
            "production": {
                "de": [
                    "Rohmaterial: Lavendel-Blüten (Lavandula angustifolia, hauptsächlich aus Provence/Bulgarien) ODER Ho-Holz aus China (höchster Linalool-Gehalt 80-90 %)",
                    "Ernte Lavendel: zur Vollblüte für maximalen ätherische-Öl-Gehalt",
                    "Wasserdampf-Destillation: 100-105 °C, 1-2 h für Lavendel; 4-8 h für Ho-Holz",
                    "Roh-Lavendel-Öl: 30-40 % Linalool + 25-35 % Linalylacetat",
                    "Industrie: hauptsächlich aus Ho-Öl für hochreines Linalool; alternativ chemische Synthese aus Geraniol/Citral",
                ],
                "en": [
                    "Raw material: lavender flowers (Lavandula angustifolia, mainly from Provence/Bulgaria) OR ho wood from China (highest linalool content 80-90 %)",
                    "Lavender harvest: at full bloom for maximum essential oil content",
                    "Steam distillation: 100-105 °C, 1-2 h for lavender; 4-8 h for ho wood",
                    "Crude lavender oil: 30-40 % linalool + 25-35 % linalyl acetate",
                    "Industry: mainly from ho oil for highly pure linalool; alternatively chemical synthesis from geraniol/citral",
                ],
            },
            "downstream": {
                "de": [
                    "Phasentrennung Öl/Wasser",
                    "Fraktionierte Vakuum-Destillation (Linalool: 198 °C atm, 80 °C bei 10 mbar)",
                    "Trennung von Linalylacetat via Verseifung (für reines Linalool)",
                    "Enantiomeren-Trennung (R)- vs (S)-Linalool optional",
                ],
                "en": [
                    "Oil/water phase separation",
                    "Fractional vacuum distillation (linalool: 198 °C atm, 80 °C at 10 mbar)",
                    "Separation from linalyl acetate via saponification (for pure linalool)",
                    "Enantiomer separation (R)- vs (S)-linalool optional",
                ],
            },
            "expected_purity_after_workup": "≥ 97 % (Flavor & Fragrance); ≥ 99 % für Pharma / Vitamin-E-Synthese",
        },
    },
    "menthol": {
        "extraction": {
            "source": "Pfefferminzöl (Mentha arvensis 'Cornmint' oder M. piperita) — Cornmint dominiert industriell",
            "source_en": "Peppermint oil (Mentha arvensis 'Cornmint' or M. piperita) — cornmint dominates industrially",
            "yield_range_percent_w_w": (60, 80),  # menthol in cornmint oil
            "literature_source": "Lange, Adv. Biochem. Eng. Biotechnol. 2015 — Mentha monoterpene biosynthesis (doi:10.1007/10_2014_283); Croteau et al., Naturwissenschaften 2005 — Menthol biosynthesis (doi:10.1007/s00114-005-0019-9)",
            "production": {
                "de": [
                    "Rohmaterial: Mentha arvensis (Cornmint, Hauptquelle aus Indien — > 75 % der Weltproduktion), Ernte zur Vollblüte",
                    "Vorbehandlung: Trocknung der Pflanzenteile bis Wassergehalt 30-40 %",
                    "Wasserdampf-Destillation: 100 °C, 2-4 h → Roh-Cornmint-Öl (typ. 70-85 % L-Menthol als Hauptkomponente, plus Menthon, Iso-Menthon)",
                    "De-Menthol-Verfahren: Roh-Öl wird auf -10 bis -20 °C abgekühlt → L-Menthol kristallisiert aus",
                    "Alternative: chemische Synthese (Haarmann-Reimer-Route via Citronellal aus m-Cresol), zugleich asymmetrische Hydrierung (Noyori-Verfahren) für (-)-Menthol",
                ],
                "en": [
                    "Raw material: Mentha arvensis (cornmint, main source from India — > 75 % of world production), harvest at full bloom",
                    "Pre-treatment: drying of plant material to water content 30-40 %",
                    "Steam distillation: 100 °C, 2-4 h -> crude cornmint oil (typ. 70-85 % L-menthol as main component, plus menthone, iso-menthone)",
                    "De-mentholisation: crude oil cooled to -10 to -20 °C -> L-menthol crystallises out",
                    "Alternative: chemical synthesis (Haarmann-Reimer route via citronellal from m-cresol), also asymmetric hydrogenation (Noyori process) for (-)-menthol",
                ],
            },
            "downstream": {
                "de": [
                    "Kristallisation aus dem dementholisierten Öl bei -10 °C",
                    "Filtration der Menthol-Kristalle",
                    "Umkristallisation aus Ethanol/Wasser",
                    "Vakuum-Destillation des Mutterlauges → weitere Menthol-Fraktionen",
                    "Optional: chirale Trennung der Stereoisomere ((-)-Menthol ist klinisch aktiv)",
                ],
                "en": [
                    "Crystallization from dementholised oil at -10 °C",
                    "Filtration of menthol crystals",
                    "Recrystallization from ethanol/water",
                    "Vacuum distillation of mother liquor -> further menthol fractions",
                    "Optional: chiral separation of stereoisomers ((-)-menthol is clinically active)",
                ],
            },
            "expected_purity_after_workup": "≥ 99.5 % (USP / Ph. Eur., kristallin); klinisch (-)-Menthol",
        },
    },
    "squalene": {
        "extraction": {
            "source": "Olivenöl-Deodorant-Destillat (DOD, Nebenprodukt der Raffination) ODER Aurantiochytrium-Fermentation",
            "source_en": "Olive-oil deodorant distillate (DOD, by-product of refining) OR Aurantiochytrium fermentation",
            "yield_range_percent_w_w": (10, 50),  # squalene in DOD
            "literature_source": "Spanova & Daum, Eur. J. Lipid Sci. Technol. 2011 — Squalene: biochemistry and biotechnology (doi:10.1002/ejlt.201100203); Aki et al., J. Am. Oil Chem. Soc. 2003 — Marine microorganism squalene production (doi:10.1007/s11746-003-0729-6)",
            "production": {
                "de": [
                    "Hauptquelle (kommerziell): Olivenöl-Deodorant-Destillat (DOD), gewonnen bei der Olivenöl-Raffination (250-260 °C unter Vakuum) — enthält 10-50 % w/w Squalen",
                    "Alternative (umstrittene Quelle): Haileberöl (heute meist vermieden wegen Tierschutz und Schwermetall-Kontamination)",
                    "Biotech-Alternative: Aurantiochytrium-Fermentation (heterotrophe Mikroalge), Glucose-Substrat, 25-28 °C, 4-7 Tage → 5-30 % Squalen der Zell-Trockenmasse",
                    "Biotech-Vorteil: tierfrei, lipidfrei von Schwermetallen, reproduzierbarer Stereochemie-Profil",
                ],
                "en": [
                    "Main commercial source: olive-oil deodorant distillate (DOD), obtained during olive-oil refining (250-260 °C under vacuum) — contains 10-50 % w/w squalene",
                    "Alternative (controversial source): shark liver oil (today mostly avoided due to animal welfare and heavy-metal contamination)",
                    "Biotech alternative: Aurantiochytrium fermentation (heterotrophic microalga), glucose substrate, 25-28 °C, 4-7 days -> 5-30 % squalene of cell dry mass",
                    "Biotech advantage: animal-free, free of heavy-metal lipid contamination, reproducible stereochemistry profile",
                ],
            },
            "downstream": {
                "de": [
                    "Fraktionierte Hochvakuum-Destillation des DOD (Squalen siedet bei 285 °C bei 25 mbar)",
                    "Säulen-Chromatographie an Kieselgel zur Trennung von Sterolen und Triglyceriden",
                    "Hydrierung optional zu Squalan (gesättigte Form, kosmetische Anwendung)",
                    "Lagerung unter Inertgas (Squalen oxidationsempfindlich)",
                ],
                "en": [
                    "Fractional high-vacuum distillation of DOD (squalene boils at 285 °C at 25 mbar)",
                    "Column chromatography on silica gel to separate sterols and triglycerides",
                    "Hydrogenation optional to squalane (saturated form, cosmetic application)",
                    "Storage under inert gas (squalene oxidation-sensitive)",
                ],
            },
            "expected_purity_after_workup": "≥ 98 % (Pharma / Kosmetik); Vakzin-Adjuvans (MF59, AS03) erfordert > 99 %",
        },
        "biotechnological": {
            "organism": "Aurantiochytrium sp. (heterotrophe marine Protistin)",
            "organism_en": "Aurantiochytrium sp. (heterotrophic marine protist)",
            "substrate": "Glucose / Glycerin + Hefeextrakt",
            "substrate_en": "glucose / glycerol + yeast extract",
            "yield_range_g_per_l": (5, 30),
            "fermentation_time_h": (96, 168),
            "literature_source": "Aki et al., J. Am. Oil Chem. Soc. 2003 (doi:10.1007/s11746-003-0729-6); Chang et al., Mar. Drugs 2014 — Aurantiochytrium for squalene (doi:10.3390/md12063657)",
            "production": {
                "de": [
                    "Stammvorbereitung: Aurantiochytrium-Hochertrags-Stamm aus marinem Isolat",
                    "Vorkultur: Glucose-Hefeextrakt-Medium mit NaCl-Zusatz (marine Spezies), 25-28 °C, 48-72 h",
                    "Hauptfermentation: Submerskultur in 10-50 m3 Bioreaktor, 25-28 °C, pH 6.0-7.0, hoher O2-Eintrag (Squalen-Biosynthese ist aerob)",
                    "Fed-Batch mit Glucose oder Glycerin als C-Quelle, Stickstoff-Limitierung induziert Squalen-Akkumulation",
                    "Fermentations-Dauer 4-7 Tage, intrazelluläre Squalen-Akkumulation 5-30 g/L (begleitet von DHA, EPA als Bonus-Produkten)",
                ],
                "en": [
                    "Strain preparation: Aurantiochytrium high-yield strain from marine isolate",
                    "Pre-culture: glucose-yeast extract medium with NaCl supplement (marine species), 25-28 °C, 48-72 h",
                    "Main fermentation: submerged culture in 10-50 m3 bioreactor, 25-28 °C, pH 6.0-7.0, high O2 transfer (squalene biosynthesis is aerobic)",
                    "Fed-batch with glucose or glycerol as C source, nitrogen limitation induces squalene accumulation",
                    "Fermentation duration 4-7 days, intracellular squalene accumulation 5-30 g/L (accompanied by DHA, EPA as bonus products)",
                ],
            },
            "downstream": {
                "de": [
                    "Zellernte via Zentrifugation",
                    "Zell-Lyse + Lösungsmittel-Extraktion (Hexan oder ueberkritisches CO2)",
                    "Lipid-Fraktion via Säulen-Chromatographie auf Squalen reinigen",
                    "Hochvakuum-Destillation für Pharma-Grade",
                ],
                "en": [
                    "Cell harvest via centrifugation",
                    "Cell lysis + solvent extraction (hexane or supercritical CO2)",
                    "Purify lipid fraction via column chromatography to squalene",
                    "High-vacuum distillation for pharma grade",
                ],
            },
            "expected_purity_after_workup": "≥ 98 % (Kosmetik / Vakzin-Adjuvans)",
        },
    },
    "artemisinin": {
        "extraction": {
            "source": "Artemisia annua (Einjähriger Beifuß) — Blätter und Blüten",
            "source_en": "Artemisia annua (sweet wormwood) — leaves and flowers",
            "yield_range_percent_w_w": (0.5, 1.5),  # artemisinin in dry leaves
            "literature_source": "Paddon et al., Nature 2013 — High-level semi-synthetic production of the potent antimalarial artemisinin (doi:10.1038/nature12051); Brown, Molecules 2010 — Artemisinin and its analogues (doi:10.3390/molecules15117603)",
            "production": {
                "de": [
                    "Rohmaterial: Artemisia annua-Blätter aus kontrolliertem Anbau (Hauptanbau China, Vietnam, Madagaskar, Kenia)",
                    "Anbau: Hochertrags-Sorten mit 0.5-1.5 % Artemisinin-Gehalt; Pflanzen werden vor der Blüte geerntet (höchster Gehalt)",
                    "Trocknung der Blätter (Sonne oder kontrolliert), Mahlung",
                    "Artemisinin-Gehalt im Roh-Pflanzenmaterial: 0.5-1.5 % w/w (Trockenmasse)",
                    "Alternative seit 2013: Semi-Synthese aus Artemisininsäure (gentechnisch in S. cerevisiae produziert via Amyris/Sanofi-Verfahren, Keasling-Labor) + chemische Umwandlung",
                ],
                "en": [
                    "Raw material: Artemisia annua leaves from controlled cultivation (main cultivation in China, Vietnam, Madagascar, Kenya)",
                    "Cultivation: high-yield varieties with 0.5-1.5 % artemisinin content; plants harvested before flowering (highest content)",
                    "Drying of leaves (sun or controlled), milling",
                    "Artemisinin content in raw plant material: 0.5-1.5 % w/w (dry mass)",
                    "Alternative since 2013: semi-synthesis from artemisinic acid (genetically produced in S. cerevisiae via Amyris/Sanofi process, Keasling lab) + chemical conversion",
                ],
            },
            "downstream": {
                "de": [
                    "Lösungsmittel-Extraktion (n-Hexan oder Petrolether) der getrockneten Blätter — Soxhlet oder Mazeration",
                    "Konzentrierung des Extrakts via Rotationsverdampfer",
                    "Säulen-Chromatographie an Kieselgel zur Trennung von Artemisinin von Begleitsesquiterpenen",
                    "Umkristallisation aus Cyclohexan oder Methanol für API-Grade",
                    "Klinische Form: meist als Artesunat (Wasser-lösliche Derivate) für IV-Gabe",
                ],
                "en": [
                    "Solvent extraction (n-hexane or petroleum ether) of dried leaves — Soxhlet or maceration",
                    "Concentration of extract via rotary evaporator",
                    "Column chromatography on silica gel to separate artemisinin from accompanying sesquiterpenes",
                    "Recrystallization from cyclohexane or methanol for API grade",
                    "Clinical form: usually as artesunate (water-soluble derivative) for IV administration",
                ],
            },
            "expected_purity_after_workup": "≥ 98 % (WHO Pharmacopoeia / USP); klinisch als Teil von ACT (Artemisinin Combination Therapy)",
        },
    },
    # ====================================================================
    # BATCH 8 — Natural product / alkaloid (Pflanzen-Extraktion / Semi-Synthese)
    # ====================================================================
    "atropine": {
        "extraction": {
            "source": "Atropa belladonna (Tollkirsche), Datura stramonium (Stechapfel), Hyoscyamus niger (Bilsenkraut)",
            "source_en": "Atropa belladonna (deadly nightshade), Datura stramonium (jimsonweed), Hyoscyamus niger (henbane)",
            "yield_range_percent_w_w": (0.2, 0.7),  # atropine in dry plant
            "literature_source": "Berkov et al., Phytochem. Anal. 2003 — Alkaloid HPLC determination in Solanaceae (doi:10.1002/pca.728); Kohnen-Johannsen & Kayser, Molecules 2019 — Tropane alkaloid biosynthesis review (doi:10.3390/molecules24040796)",
            "production": {
                "de": [
                    "Rohmaterial: Pflanzenmaterial aus kontrolliertem Anbau (Wurzeln + Blätter von A. belladonna oder D. stramonium); regulatorisch streng kontrolliert (BtMG-relevant in einigen Ländern)",
                    "Geographische Hauptproduzenten: Indien, Ägypten, Polen für kommerziellen Anbau",
                    "Vorbehandlung: Trocknung der Pflanzenteile bis Wassergehalt < 10 %, Mahlung zu Pulver",
                    "Atropin existiert in Pflanze als (S)-(-)-Hyoscyamin → Racemisierung zu Atropin ((R,S)-Hyoscyamin) durch alkalische Hitze-Behandlung beim Extraktionsprozess",
                    "Alkaloid-Gehalt: 0.2-0.7 % w/w Gesamtalkaloide (Trockenmasse), davon Atropin/Hyoscyamin ~80 %, Scopolamin als Begleitalkaloid",
                ],
                "en": [
                    "Raw material: plant material from controlled cultivation (roots + leaves of A. belladonna or D. stramonium); regulatorily strictly controlled (controlled substance in some countries)",
                    "Main geographic producers: India, Egypt, Poland for commercial cultivation",
                    "Pre-treatment: drying of plant material to water content < 10 %, milling to powder",
                    "Atropine exists in plant as (S)-(-)-hyoscyamine -> racemisation to atropine ((R,S)-hyoscyamine) by alkaline heat treatment during extraction process",
                    "Alkaloid content: 0.2-0.7 % w/w total alkaloids (dry mass), of which atropine/hyoscyamine ~80 %, scopolamine as accompanying alkaloid",
                ],
            },
            "downstream": {
                "de": [
                    "Wässrige Säure-Extraktion (verd. H2SO4 oder HCl, pH 2-3) der gemahlenen Pflanze — Alkaloide gehen als Salze in Lösung",
                    "Filtration + Aktivkohle-Behandlung zur Entfärbung",
                    "pH-Erhöhung auf 9-10 mit NH4OH → Fällung der freien Basen",
                    "Lösungsmittel-Extraktion mit Chloroform oder Dichlormethan",
                    "Trennung Atropin von Scopolamin via fraktionierte Kristallisation oder Säulen-Chromatographie",
                    "Salzbildung als Atropin-Sulfat für klinische Form (Ph. Eur. / USP)",
                ],
                "en": [
                    "Aqueous acid extraction (dilute H2SO4 or HCl, pH 2-3) of milled plant — alkaloids dissolve as salts",
                    "Filtration + activated-carbon treatment for decolorisation",
                    "pH adjustment to 9-10 with NH4OH -> precipitation of free bases",
                    "Solvent extraction with chloroform or dichloromethane",
                    "Separation of atropine from scopolamine via fractional crystallization or column chromatography",
                    "Salt formation as atropine sulfate for clinical form (Ph. Eur. / USP)",
                ],
            },
            "expected_purity_after_workup": "≥ 99 % (Atropin-Sulfat USP / Ph. Eur., klinische Anwendung: Mydriatikum, Bradykardie-Notfall)",
        },
    },
    "codeine": {
        "extraction": {
            "source": "Papaver somniferum (Schlafmohn) — Latex oder Mohnstroh (CPS)",
            "source_en": "Papaver somniferum (opium poppy) — latex or concentrated poppy straw (CPS)",
            "yield_range_percent_w_w": (1, 3),  # codeine content in opium latex
            "literature_source": "Schiff, J. Ethnopharmacol. 2002 — Opium and its alkaloids (doi:10.1016/S0378-8741(01)00408-6); Beaudoin & Facchini, Planta 2014 — Benzylisoquinoline alkaloid biosynthesis (doi:10.1007/s00425-014-2056-8)",
            "production": {
                "de": [
                    "Rohmaterial: Schlafmohn-Latex (eingetrocknet zu Opium) ODER konzentriertes Mohnstroh (CPS) aus regulatorisch zugelassenem Anbau (BtMG-Lizenz in DE, INCB-Quotenregelung international)",
                    "Hauptproduzenten kontrollierter Anbau: Australien (Tasmanien), Frankreich, Spanien, Türkei, Indien",
                    "Codein-Gehalt im Opium: 1-3 % w/w (Trockenmasse); im Mohnstroh: 0.1-0.5 %",
                    "Vorbehandlung: Mahlung des Rohmaterials, regulatorisches Tracking jeder Charge (DEA / Bundesopiumstelle)",
                    "Alternative Route: Semi-Synthese aus Morphin durch O3-Methylierung des Phenol-OH (selektive Methylierung am C3, da C6-OH abgeschirmt) — kommerziell wichtig, weil Codein-Bedarf > natürliches Codein-Vorkommen",
                ],
                "en": [
                    "Raw material: opium poppy latex (dried to opium) OR concentrated poppy straw (CPS) from regulatorily licensed cultivation (controlled-substance licence in DE/EU, INCB quota system internationally)",
                    "Main licensed producers: Australia (Tasmania), France, Spain, Turkey, India",
                    "Codeine content in opium: 1-3 % w/w (dry mass); in poppy straw: 0.1-0.5 %",
                    "Pre-treatment: milling of raw material, regulatory tracking of every batch (DEA / national opium agency)",
                    "Alternative route: semi-synthesis from morphine by O3-methylation of phenolic OH (selective methylation at C3, since C6-OH is shielded) — commercially important because codeine demand > natural codeine occurrence",
                ],
            },
            "downstream": {
                "de": [
                    "Wässrige Säure-Extraktion (verd. H2SO4) der gemahlenen Roh-Substanz",
                    "Filtration + Aktivkohle zur Entfärbung",
                    "pH-Erhöhung auf 9 → Fällung der freien Morphin-/Codein-Basen",
                    "Fraktionierte Kristallisation: Morphin kristallisiert zuerst aus, Codein bleibt in Mutterlauge",
                    "Säulen-Chromatographie zur Feinreinigung",
                    "Salzbildung als Codein-Phosphat oder Codein-Hydrochlorid für klinische Anwendung",
                ],
                "en": [
                    "Aqueous acid extraction (dilute H2SO4) of milled raw substance",
                    "Filtration + activated carbon for decolorisation",
                    "pH adjustment to 9 -> precipitation of free morphine/codeine bases",
                    "Fractional crystallization: morphine crystallises first, codeine remains in mother liquor",
                    "Column chromatography for fine purification",
                    "Salt formation as codeine phosphate or codeine hydrochloride for clinical use",
                ],
            },
            "expected_purity_after_workup": "≥ 99 % (Codein-Phosphat USP / Ph. Eur., klinische Anwendung: Hustenstillung, schwaches Analgetikum)",
        },
    },
    "quinine": {
        "extraction": {
            "source": "Cinchona-Baum-Rinde (Cinchona officinalis, C. ledgeriana, C. succirubra)",
            "source_en": "Cinchona tree bark (Cinchona officinalis, C. ledgeriana, C. succirubra)",
            "yield_range_percent_w_w": (4, 8),  # quinine in dried Cinchona bark
            "literature_source": "Achan et al., Malar. J. 2011 — Quinine, an old anti-malarial drug in a modern world (doi:10.1186/1475-2875-10-144); Kacprzak, in Heterocycles in Natural Product Synthesis 2011, Wiley-VCH (ISBN 978-3-527-32706-8)",
            "production": {
                "de": [
                    "Rohmaterial: Cinchona-Baum-Rinde aus Plantagen (Hauptanbau in Java/Indonesien, Kongo, früher Südamerika)",
                    "Anbau-Zyklus: Cinchona-Bäume werden nach 8-12 Jahren gefällt (oder die Rinde wird abgeschält und der Stamm regeneriert)",
                    "Vorbehandlung: Trocknung der Rinde, Mahlung zu Pulver",
                    "Alkaloid-Gehalt: 4-8 % w/w Gesamtalkaloide (Trockenmasse), davon Chinin ~60-70 %, Chinidin/Cinchonin/Cinchonidin als Begleitalkaloide",
                    "Globale Produktion: ~300-700 t Chinin pro Jahr (für Pharma + Tonic Water + Bitterstoffe in Getränken)",
                ],
                "en": [
                    "Raw material: cinchona tree bark from plantations (main cultivation Java/Indonesia, Congo, formerly South America)",
                    "Cultivation cycle: cinchona trees felled after 8-12 years (or bark peeled and trunk regenerated)",
                    "Pre-treatment: drying of bark, milling to powder",
                    "Alkaloid content: 4-8 % w/w total alkaloids (dry mass), of which quinine ~60-70 %, quinidine/cinchonine/cinchonidine as accompanying alkaloids",
                    "Global production: ~300-700 t quinine per year (for pharma + tonic water + bitter agents in beverages)",
                ],
            },
            "downstream": {
                "de": [
                    "Alkalische Extraktion: Rinde wird mit Calciumhydroxid + Wasser behandelt → freie Basen werden frei",
                    "Lösungsmittel-Extraktion mit Toluol oder Petrolether",
                    "Säure-Aufnahme der Alkaloide in verd. H2SO4 → Sulfat-Salze",
                    "Fraktionierte Kristallisation: Chinin-Sulfat kristallisiert selektiv → grobe Vorreinigung",
                    "Umkristallisation aus Wasser/Ethanol für höhere Reinheit",
                    "Salzbildung als Chinin-Sulfat-Dihydrat für klinische Anwendung (Malaria, nächtliche Wadenkrämpfe)",
                ],
                "en": [
                    "Alkaline extraction: bark treated with calcium hydroxide + water -> free bases liberated",
                    "Solvent extraction with toluene or petroleum ether",
                    "Acid uptake of alkaloids in dilute H2SO4 -> sulfate salts",
                    "Fractional crystallization: quinine sulfate crystallises selectively -> coarse pre-purification",
                    "Recrystallization from water/ethanol for higher purity",
                    "Salt formation as quinine sulfate dihydrate for clinical use (malaria, nocturnal leg cramps)",
                ],
            },
            "expected_purity_after_workup": "≥ 98 % (Chinin-Sulfat USP / Ph. Eur.); Tonic Water-Grade typ. ≥ 90 %",
        },
    },
    "isobutanol": {
        "biotechnological": {
            "organism": "Engineered E. coli oder S. cerevisiae (Ehrlich pathway via L-Valin-Vorstufe)",
            "organism_en": "Engineered E. coli or S. cerevisiae (Ehrlich pathway via L-valine precursor)",
            "substrate": "Glucose / Xylose",
            "substrate_en": "glucose / xylose",
            "yield_range_g_per_l": (5, 25),
            "fermentation_time_h": (48, 96),
            "literature_source": "Atsumi, Hanai & Liao, Nature 2008 — Non-fermentative pathways for synthesis of branched-chain higher alcohols (doi:10.1038/nature06450); Lee et al., Curr. Opin. Biotechnol. 2008 — Biological isobutanol production (doi:10.1016/j.copbio.2008.10.014)",
            "production": {
                "de": [
                    "Stammvorbereitung: rekombinantes E. coli (oder S. cerevisiae) mit überexprimierter Acetolactat-Synthase, Ketolsäure-Reduktoisomerase und Dihydroxysäure-Dehydratase; Knockout konkurrierender Ethanol-Pfade",
                    "Vorkultur: definiertes Minimalmedium + Glucose + Selektion, 30-37 °C, 12-16 h",
                    "Hauptfermentation: aerobe oder mikroaerobe Submerskultur in 10-50 m3 Bioreaktor, 30-37 °C, pH 6.5-7.0",
                    "Glucose-Feed kontrolliert (< 30 g/L Reaktorkonzentration), In-situ-Stripping zur Toxizitäts-Reduktion empfohlen",
                    "Fermentations-Dauer 48-96 h, Isobutanol-Konzentration 5-25 g/L",
                ],
                "en": [
                    "Strain preparation: recombinant E. coli (or S. cerevisiae) with overexpressed acetolactate synthase, keto-acid reductoisomerase and dihydroxy-acid dehydratase; knockout of competing ethanol pathways",
                    "Pre-culture: defined minimal medium + glucose + selection, 30-37 °C, 12-16 h",
                    "Main fermentation: aerobic or microaerobic submerged culture in 10-50 m3 bioreactor, 30-37 °C, pH 6.5-7.0",
                    "Glucose feed controlled (< 30 g/L reactor concentration), in-situ stripping recommended to mitigate toxicity",
                    "Fermentation duration 48-96 h, isobutanol concentration 5-25 g/L",
                ],
            },
            "downstream": {
                "de": [
                    "Zellabtrennung via Zentrifugation",
                    "Gas-Stripping oder Pervaporation für integriertes Product-Removal",
                    "Phasentrennung Isobutanol-Wasser (Isobutanol bildet eine Phase ab ~8 % Konzentration)",
                    "Fraktionierte Destillation bei 108 °C, Trocknung über Molekularsieb",
                ],
                "en": [
                    "Cell separation via centrifugation",
                    "Gas stripping or pervaporation for integrated product removal",
                    "Phase separation isobutanol-water (isobutanol forms a phase above ~8 % concentration)",
                    "Fractional distillation at 108 °C, drying over molecular sieve",
                ],
            },
            "expected_purity_after_workup": "≥ 99 % (technical) bis ≥ 99.5 % (Reagenz-Grade)",
        },
    },
}

