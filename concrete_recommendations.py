"""Concrete, parameterised recommendations for downstream reporting.

Where the engine produces qualitative outputs ("purification difficulty: high"),
this module turns them into actionable, numeric recommendations:
methods + solvents + expected yields + processing times.

All textual fields are localized — the build_concrete_recommendations()
function takes a `lang` parameter ('de' | 'en') and returns the appropriate
language. Numeric ranges (yield, time) are language-independent.
"""

from typing import Dict, Any, List, Optional, Tuple


# Per-known-molecule starting points. Keys 'production_route' and 'downstream'
# are language-tagged so reports come out cohesive. The 'source' field is a
# free-text reference for traceability — surfaced inline and on the references
# page of the PDF report.
MOLECULE_HINTS: Dict[str, Dict[str, Any]] = {
    "vanillin": {
        "biotechnological": {
            "organism": "S. cerevisiae oder A. niger (engineered)",
            "organism_en": "S. cerevisiae or A. niger (engineered)",
            "substrate": "Glucose oder Ferulasäure",
            "substrate_en": "glucose or ferulic acid",
            "yield_range_g_per_l": (0.3, 1.5),
            "fermentation_time_h": (48, 96),
            "literature_source": "Hansen et al., Appl. Environ. Microbiol. 2009 — De novo biosynthesis of vanillin in S. cerevisiae (doi:10.1128/AEM.02074-08); Walton et al., Phytochemistry 2003 (review, doi:10.1016/S0031-9422(03)00149-3)",
            "downstream": {
                "de": [
                    "Zellabtrennung via Zentrifugation oder Mikrofiltration",
                    "Flüssig-flüssig-Extraktion mit Ethylacetat",
                    "Aufkonzentrierung im Rotationsverdampfer (40 °C, ~100 mbar)",
                    "Umkristallisation aus Wasser/Ethanol für API-grade Reinheit",
                ],
                "en": [
                    "Cell removal via centrifugation or microfiltration",
                    "Liquid-liquid extraction with ethyl acetate",
                    "Concentration via rotary evaporation (40 °C, ~100 mbar)",
                    "Recrystallization from water/ethanol for API-grade purity",
                ],
            },
            "expected_purity_after_workup": "95–99 %",
        },
        "chemical": {
            "reagents": "Guaiacol, Glyoxylsäure, NaOH, anschließend Oxidation",
            "reagents_en": "Guaiacol, glyoxylic acid, NaOH, followed by oxidation",
            "yield_range_percent": (60, 80),
            "literature_source": "Borregaard / Rhodia industrial process (Glyoxylsäure-Route); Fache et al., Acta Chim. Slov. 2000 (overview of vanillin synthesis routes)",
            "downstream": {
                "de": [
                    "Säure-/Base-Extraktion zur Abtrennung unreagierter Edukte",
                    "Umkristallisation aus heißem Wasser",
                ],
                "en": [
                    "Acid/base extraction to remove unreacted starting materials",
                    "Recrystallization from hot water",
                ],
            },
        },
        "extraction": {
            "source": "Vanilla planifolia Schoten (kuriert)",
            "source_en": "Vanilla planifolia pods (cured)",
            "yield_range_percent_w_w": (1.5, 2.5),
            "literature_source": "Sinha et al., Crit. Rev. Food Sci. Nutr. 2008 — Vanilla as a flavour ingredient (doi:10.1080/10408390701764235)",
            "downstream": {
                "de": [
                    "Perkolation mit Ethanol/Wasser (60:40)",
                    "Konzentration im Rotationsverdampfer",
                ],
                "en": [
                    "Percolation with ethanol/water (60:40)",
                    "Concentration via rotary evaporation",
                ],
            },
        },
    },
    "aspirin": {
        "chemical": {
            "reagents": "Salicylsäure + Essigsäureanhydrid, kat. H2SO4",
            "reagents_en": "Salicylic acid + acetic anhydride, cat. H2SO4",
            "yield_range_percent": (75, 92),
            "literature_source": "Vollhardt & Schore, Organic Chemistry, 8th ed. 2018, Macmillan/Freeman (ISBN 978-1-319-07945-1) — Acetylation of salicylic acid (standard textbook synthesis)",
            "downstream": {
                "de": [
                    "Fällung in Eiswasser",
                    "Filtration und Waschen mit kaltem Wasser",
                    "Umkristallisation aus Ethanol/Wasser",
                ],
                "en": [
                    "Precipitation in ice water",
                    "Filtration and washing with cold water",
                    "Recrystallization from ethanol/water",
                ],
            },
            "expected_purity_after_workup": "≥ 99 %",
        },
    },
    "ibuprofen": {
        "chemical": {
            "reagents": "Isobutylbenzol → Friedel-Crafts → BHC- oder Boots-Verfahren",
            "reagents_en": "Isobutylbenzene → Friedel-Crafts → BHC or Boots process",
            "yield_range_percent": (50, 75),
            "literature_source": "BHC Boots-Hoechst-Celanese process — US Patent 4,981,995 (Elango & Davenport 1991); Presidential Green Chemistry Award 1997; Cann & Connelly, Real-World Cases in Green Chemistry, ACS 2000 (ISBN 978-0-8412-3733-8)",
            "downstream": {
                "de": [
                    "Säulenchromatographie (SiO2, Hexan/EtOAc-Gradient)",
                    "Umkristallisation aus n-Hexan",
                ],
                "en": [
                    "Column chromatography (SiO2, hexane/EtOAc gradient)",
                    "Recrystallization from n-hexane",
                ],
            },
            "expected_purity_after_workup": "98–99 %",
        },
    },
    "ethanol": {
        "biotechnological": {
            "organism": "S. cerevisiae",
            "organism_en": "S. cerevisiae",
            "substrate": "Glucose oder Stärke-Hydrolysat",
            "substrate_en": "glucose or starch hydrolysate",
            "yield_range_g_per_l": (80, 130),
            "fermentation_time_h": (48, 72),
            "literature_source": "Doran, Bioprocess Engineering Principles, 2nd ed. 2013, Academic Press (ISBN 978-0-12-220851-5) — yeast ethanol fermentation; Jacques et al., The Alcohol Textbook, 5th ed. 2003, Nottingham Univ. Press (ISBN 978-1-89745-541-6)",
            "downstream": {
                "de": [
                    "Destillation (klassisch oder Rektifikationskolonne)",
                    "Molekularsieb-Trocknung für absoluten Ethanol",
                ],
                "en": [
                    "Distillation (classical or rectification column)",
                    "Molecular-sieve drying for absolute ethanol",
                ],
            },
            "expected_purity_after_workup": "95–99,5 %",
        },
    },
    "penicillin": {
        "biotechnological": {
            "organism": "P. chrysogenum",
            "organism_en": "P. chrysogenum",
            "substrate": "Lactose / Maisquellwasser + Phenylessigsäure",
            "substrate_en": "lactose / corn-steep liquor + phenylacetic acid",
            "yield_range_g_per_l": (5, 50),
            "fermentation_time_h": (120, 200),
            "literature_source": "Walsh G., Pharmaceutical Biotechnology, 2nd ed. 2018, Wiley (ISBN 978-1-119-11518-7) — penicillin fermentation; Elander, Appl. Microbiol. Biotechnol. 2003 (industrial production review, doi:10.1007/s00253-003-1274-y)",
            "downstream": {
                "de": [
                    "Filtration des Myzels",
                    "Liquid-liquid Extraktion bei pH 2 (n-Butylacetat)",
                    "Rückextraktion bei pH 7,5",
                    "Kristallisation als Kaliumsalz",
                ],
                "en": [
                    "Filtration of mycelium",
                    "Liquid-liquid extraction at pH 2 (n-butyl acetate)",
                    "Back-extraction at pH 7.5",
                    "Crystallization as potassium salt",
                ],
            },
            "expected_purity_after_workup": "≥ 98 %",
        },
    },
    # ---------------- Extension Set: 12 high-volume biotech molecules ----------------
    "citric acid": {
        "biotechnological": {
            "organism": "Aspergillus niger",
            "organism_en": "Aspergillus niger",
            "substrate": "Saccharose / Melasse",
            "substrate_en": "sucrose / molasses",
            "yield_range_g_per_l": (150, 200),
            "fermentation_time_h": (120, 200),
            "literature_source": "Berovic & Legisa, Biotechnol. Annu. Rev. 2007 — Citric acid production (doi:10.1016/S1387-2656(07)13011-8); Show et al., Bioengineered 2015 (doi:10.4161/21655979.2014.978001)",
            "downstream": {
                "de": [
                    "Myzel-Abtrennung via Filtration",
                    "Fällung als Calciumcitrat mit Ca(OH)₂",
                    "Zersetzung mit H2SO4 → Calciumsulfat-Niederschlag + freie Citronensäure",
                    "Aktivkohle-Behandlung, Eindampfung, Kristallisation",
                ],
                "en": [
                    "Mycelium separation via filtration",
                    "Precipitation as calcium citrate with Ca(OH)₂",
                    "Decomposition with H2SO4 → calcium sulfate solid + free citric acid",
                    "Activated-carbon polish, evaporation, crystallization",
                ],
            },
            "expected_purity_after_workup": "≥ 99,5 % (food grade)",
        },
    },
    "lactic acid": {
        "biotechnological": {
            "organism": "Lactobacillus delbrueckii / Bacillus coagulans",
            "organism_en": "Lactobacillus delbrueckii / Bacillus coagulans",
            "substrate": "Glucose / Stärke-Hydrolysat",
            "substrate_en": "glucose / starch hydrolysate",
            "yield_range_g_per_l": (90, 180),
            "fermentation_time_h": (24, 72),
            "literature_source": "Abdel-Rahman et al., Biotechnol. Adv. 2013 — Lactic acid production review (doi:10.1016/j.biotechadv.2013.04.002); Datta & Henry, J. Chem. Technol. Biotechnol. 2006 (doi:10.1002/jctb.1486)",
            "downstream": {
                "de": [
                    "Zellabtrennung via Mikrofiltration",
                    "Fällung als Calciumlactat (klassisch) oder Membran-/Elektrodialyse",
                    "Säuerung mit H2SO4, Aktivkohle-Reinigung",
                    "Veresterung/Hydrolyse zur Reinheitssteigerung (PLA-Grade)",
                ],
                "en": [
                    "Cell separation via microfiltration",
                    "Precipitation as calcium lactate (classical) or membrane / electrodialysis",
                    "Acidification with H2SO4, activated-carbon cleanup",
                    "Esterification / hydrolysis for higher purity (PLA grade)",
                ],
            },
            "expected_purity_after_workup": "≥ 88 % (technical) bis ≥ 99,5 % (PLA grade)",
        },
    },
    "succinic acid": {
        "biotechnological": {
            "organism": "Actinobacillus succinogenes / Basfia succiniciproducens",
            "organism_en": "Actinobacillus succinogenes / Basfia succiniciproducens",
            "substrate": "Glucose + CO2",
            "substrate_en": "glucose + CO2",
            "yield_range_g_per_l": (50, 100),
            "fermentation_time_h": (40, 80),
            "literature_source": "Cok et al., Biofuels Bioprod. Bioref. 2014 — Succinic acid biorefinery (doi:10.1002/bbb.1427); Beauprez et al., Process Biochem. 2010 (doi:10.1016/j.procbio.2010.03.035)",
            "downstream": {
                "de": [
                    "Zellabtrennung (Filtration / Zentrifugation)",
                    "Aufkonzentrierung über Vakuumeindampfung oder Membranfiltration",
                    "Direkt-Kristallisation aus saurer Brühe",
                    "Optional: Umkristallisation aus Wasser",
                ],
                "en": [
                    "Cell separation (filtration / centrifugation)",
                    "Concentration via vacuum evaporation or membrane filtration",
                    "Direct crystallization from acidic broth",
                    "Optional: recrystallization from water",
                ],
            },
            "expected_purity_after_workup": "≥ 99 %",
        },
    },
    "paracetamol": {
        "chemical": {
            "reagents": "p-Aminophenol + Essigsäureanhydrid (Acetylierung)",
            "reagents_en": "p-aminophenol + acetic anhydride (acetylation)",
            "yield_range_percent": (85, 95),
            "literature_source": "Joncour et al., Green Chem. 2014 — Paracetamol synthesis review (doi:10.1039/C4GC00603H); Ellis, Paracetamol: A Curriculum Resource, RSC 2002 (ISBN 978-0-85404-365-1)",
            "downstream": {
                "de": [
                    "Fällung in kaltem Wasser",
                    "Filtration und Waschen",
                    "Umkristallisation aus Wasser/Ethanol",
                    "Aktivkohle-Behandlung für API-Grade",
                ],
                "en": [
                    "Precipitation in cold water",
                    "Filtration and washing",
                    "Recrystallization from water/ethanol",
                    "Activated-carbon polish for API grade",
                ],
            },
            "expected_purity_after_workup": "≥ 99,5 % (Ph. Eur. / USP)",
        },
    },
    "glucoamylase": {
        "biotechnological": {
            "organism": "Aspergillus niger",
            "organism_en": "Aspergillus niger",
            "substrate": "Komplexes Medium (Stärke/Glucose, Maisquellwasser)",
            "substrate_en": "complex medium (starch/glucose, corn-steep liquor)",
            "yield_range_g_per_l": (5, 20),
            "fermentation_time_h": (120, 168),
            "literature_source": "Norouzian et al., Biotechnol. Adv. 2006 — Fungal glucoamylase production (doi:10.1016/j.biotechadv.2005.06.003); Kumar & Satyanarayana, Crit. Rev. Biotechnol. 2009 (doi:10.1080/07388550802479237)",
            "downstream": {
                "de": [
                    "Myzel-Abtrennung (Tiefenfiltration)",
                    "Aufkonzentrierung via Ultrafiltration",
                    "Optional: Säulenchromatographie für Hochreinheit",
                    "Formulierung (Stabilisatoren, Konservierungsmittel) für Endprodukt",
                ],
                "en": [
                    "Mycelium separation (depth filtration)",
                    "Concentration via ultrafiltration",
                    "Optional: column chromatography for high purity",
                    "Formulation (stabilizers, preservatives) for final product",
                ],
            },
            "expected_purity_after_workup": "Industrial grade (Mehrkomponenten-Formulierung); pharma grade > 95 %",
        },
    },
    "glucose isomerase": {
        "biotechnological": {
            "organism": "Streptomyces rubiginosus / S. murinus",
            "organism_en": "Streptomyces rubiginosus / S. murinus",
            "substrate": "Komplexes Medium",
            "substrate_en": "complex medium",
            "yield_range_g_per_l": (1, 5),
            "fermentation_time_h": (48, 96),
            "literature_source": "Bhosale et al., Microbiol. Rev. 1996 — Glucose isomerase comprehensive review (doi:10.1128/mr.60.2.280-300.1996); Crabb & Mitchinson, Trends Biotechnol. 1997 (doi:10.1016/S0167-7799(97)01076-5)",
            "downstream": {
                "de": [
                    "Zellaufschluss (mechanisch / chemisch)",
                    "Klärung via Zentrifugation und Filtration",
                    "Immobilisierung auf Träger (z. B. quervernetzten Aggregaten oder DEAE-Cellulose) — kritischer Schritt für HFCS",
                    "Verpackung in Festbettreaktoren für kontinuierliche Isomerisierung",
                ],
                "en": [
                    "Cell lysis (mechanical / chemical)",
                    "Clarification via centrifugation and filtration",
                    "Immobilization on carrier (e.g. cross-linked aggregates or DEAE-cellulose) — critical step for HFCS",
                    "Packing into fixed-bed reactors for continuous isomerization",
                ],
            },
            "expected_purity_after_workup": "Immobilized industrial form (nicht als isolierte freie Enzymlösung verkauft)",
        },
    },
    "asparaginase": {
        "biotechnological": {
            "organism": "E. coli (Stamm K-12)",
            "organism_en": "E. coli (strain K-12)",
            "substrate": "Komplexes Medium",
            "substrate_en": "complex medium",
            "yield_range_g_per_l": (0.5, 2),
            "fermentation_time_h": (24, 48),
            "literature_source": "Kotzia & Labrou, J. Biotechnol. 2007 — Microbial L-asparaginase production (doi:10.1016/j.jbiotec.2007.07.939); Cachumba et al., Braz. J. Microbiol. 2016 (doi:10.1016/j.bjm.2016.10.004)",
            "downstream": {
                "de": [
                    "Zellaufschluss (Hochdruckhomogenisator)",
                    "Klärung + Ammonsulfat-Fällung",
                    "Anionenaustausch-Chromatographie (DEAE oder Q-Sepharose)",
                    "Größenausschluss-Chromatographie als Polishing-Schritt",
                    "Lyophilisation der Endformulierung",
                ],
                "en": [
                    "Cell lysis (high-pressure homogenizer)",
                    "Clarification + ammonium sulfate precipitation",
                    "Anion-exchange chromatography (DEAE or Q-Sepharose)",
                    "Size-exclusion chromatography as polishing step",
                    "Lyophilization of final formulation",
                ],
            },
            "expected_purity_after_workup": "≥ 99 % (Pharma / klinisch)",
        },
    },
    "pembrolizumab": {
        "biotechnological": {
            "organism": "CHO-Zellen (Suspensionskultur)",
            "organism_en": "CHO cells (suspension culture)",
            "substrate": "Chemisch definiertes Medium",
            "substrate_en": "chemically defined medium",
            "yield_range_g_per_l": (3, 8),
            "fermentation_time_h": (288, 360),  # 12-15 days
            "literature_source": "Liu et al., mAbs 2010 — Antibody recovery and purification (doi:10.4161/mabs.2.5.12796); Walsh G., Nat. Biotechnol. 2018 — Biopharmaceutical benchmarks (doi:10.1038/nbt.4313)",
            "downstream": {
                "de": [
                    "Klärung via Tiefenfiltration",
                    "Protein-A-Affinitätschromatographie (Capture)",
                    "Niedrig-pH-Virusinaktivierung",
                    "Ionenaustausch-Polishing (CEX + AEX im Flow-through)",
                    "Virusfiltration (20 nm)",
                    "Ultrafiltration / Pufferaustausch (TFF)",
                ],
                "en": [
                    "Clarification via depth filtration",
                    "Protein-A affinity chromatography (capture)",
                    "Low-pH virus inactivation",
                    "Ion-exchange polishing (CEX + AEX in flow-through)",
                    "Virus filtration (20 nm)",
                    "Ultrafiltration / buffer exchange (TFF)",
                ],
            },
            "expected_purity_after_workup": "≥ 99 % (cGMP / klinisch)",
        },
    },
    "trastuzumab": {
        "biotechnological": {
            "organism": "CHO-Zellen (Suspensionskultur)",
            "organism_en": "CHO cells (suspension culture)",
            "substrate": "Chemisch definiertes Medium (Glucose, Glutamin, Spurenelemente, Antifoam)",
            "substrate_en": "chemically defined medium (glucose, glutamine, trace elements, antifoam)",
            "yield_range_g_per_l": (3, 10),
            "fermentation_time_h": (288, 360),  # 12-15 days fed-batch
            "literature_source": "Liu et al., mAbs 2010 — Antibody recovery and purification (doi:10.4161/mabs.2.5.12796); Walsh G., Nat. Biotechnol. 2018 — Biopharmaceutical benchmarks (doi:10.1038/nbt.4313); Wurm FM., Nat. Biotechnol. 2004 — Production of recombinant protein therapeutics in CHO (doi:10.1038/nbt1026)",
            "downstream": {
                "de": [
                    "Klärung via Tiefenfiltration + 0,2 µm Sterilfiltration",
                    "Protein-A-Affinitätschromatographie (Capture, > 95 % Reinheit ab Stufe 1)",
                    "Niedrig-pH-Virusinaktivierung (pH 3,5 für 60 min)",
                    "Kationenaustausch-Chromatographie (CEX, Bind-Elute)",
                    "Anionenaustausch im Flow-through (AEX) + Virusfiltration (Planova 20N)",
                    "Ultrafiltration / Pufferaustausch (TFF, finale Formulierung)",
                ],
                "en": [
                    "Clarification via depth filtration + 0.2 µm sterile filtration",
                    "Protein-A affinity chromatography (capture, > 95 % purity from step 1)",
                    "Low-pH virus inactivation (pH 3.5 for 60 min)",
                    "Cation-exchange chromatography (CEX, bind-elute)",
                    "Anion-exchange in flow-through (AEX) + virus filtration (Planova 20N)",
                    "Ultrafiltration / buffer exchange (TFF, final formulation)",
                ],
            },
            "expected_purity_after_workup": "≥ 99 % (cGMP / klinisch, biosimilar-equivalent)",
        },
    },
    "β-carotene": {
        "biotechnological": {
            "organism": "Blakeslea trispora oder engineered Yarrowia lipolytica",
            "organism_en": "Blakeslea trispora or engineered Yarrowia lipolytica",
            "substrate": "Glucose / pflanzliche Öle",
            "substrate_en": "glucose / vegetable oils",
            "yield_range_g_per_l": (1, 5),
            "fermentation_time_h": (96, 168),
            "literature_source": "Mata-Gómez et al., Microb. Cell Fact. 2014 — Microbial carotenoid production (doi:10.1186/1475-2859-13-12); Lopes et al., World J. Microbiol. Biotechnol. 2019 — Blakeslea review (doi:10.1007/s11274-019-2647-4)",
            "downstream": {
                "de": [
                    "Zellaufschluss + Lipid-Extraktion mit Hexan / Pflanzenöl",
                    "Saponifikation zur Entfernung von Fettsäuren",
                    "Kristallisation aus Aceton oder Ethanol",
                    "Filtration und Trocknung; Lichtschutz erforderlich",
                ],
                "en": [
                    "Cell lysis + lipid extraction with hexane / vegetable oil",
                    "Saponification to remove fatty acids",
                    "Crystallization from acetone or ethanol",
                    "Filtration and drying; light protection required",
                ],
            },
            "expected_purity_after_workup": "≥ 96 % (food/feed grade)",
        },
    },
    "astaxanthin": {
        "biotechnological": {
            "organism": "Haematococcus pluvialis (Algen) oder Phaffia rhodozyma",
            "organism_en": "Haematococcus pluvialis (algae) or Phaffia rhodozyma",
            "substrate": "Photoautotroph (Algen) oder Glucose (Phaffia)",
            "substrate_en": "photoautotrophic (algae) or glucose (Phaffia)",
            "yield_range_percent_w_w": (1.5, 5.0),  # % of dry biomass
            "fermentation_time_h": (240, 480),  # algae 10-20 days
            "literature_source": "Shah et al., Front. Plant Sci. 2016 — Astaxanthin from Haematococcus (doi:10.3389/fpls.2016.00531); Rodríguez-Sáiz et al., Appl. Microbiol. Biotechnol. 2010 — Phaffia review (doi:10.1007/s00253-010-2520-8)",
            "downstream": {
                "de": [
                    "Zellaufschluss (mechanisch — Cyste/Zellwand bei Haematococcus)",
                    "Lösungsmittel-Extraktion (CO2-überkritisch oder Ethanol/Aceton)",
                    "Säulen-Chromatographie zur Trennung von Stereoisomeren",
                    "Verkapselung (Mikrokapseln) zur Stabilisierung",
                ],
                "en": [
                    "Cell lysis (mechanical — Haematococcus cyst wall is robust)",
                    "Solvent extraction (supercritical CO2 or ethanol/acetone)",
                    "Column chromatography for stereoisomer separation",
                    "Encapsulation (microcapsules) for stabilization",
                ],
            },
            "expected_purity_after_workup": "≥ 90 % (Aquakultur), ≥ 97 % (Nahrungsergänzung)",
        },
    },
    "caffeine": {
        "extraction": {
            "source": "Tee-/Kaffeebohnen-Entkoffeinierung (Nebenprodukt)",
            "source_en": "tea / coffee-bean decaffeination (by-product)",
            "yield_range_percent_w_w": (1.0, 4.0),  # of dry coffee/tea
            "literature_source": "Heilmann W., Decaffeination Methods, Lebensmittelchemie 2001 (DLG-Verlag); Ramalakshmi & Raghavan, Crit. Rev. Food Sci. Nutr. 1999 (doi:10.1080/10408699991279231)",
            "downstream": {
                "de": [
                    "Überkritische CO2-Extraktion oder Wasser/Lösungsmittel-Extraktion",
                    "Aktivkohle-Adsorption + Desorption",
                    "Eindampfung und Kristallisation",
                    "Umkristallisation für USP/Ph.Eur.-Reinheit",
                ],
                "en": [
                    "Supercritical CO2 extraction or water / solvent extraction",
                    "Activated-carbon adsorption + desorption",
                    "Evaporation and crystallization",
                    "Recrystallization for USP / Ph.Eur. grade",
                ],
            },
            "expected_purity_after_workup": "≥ 99 %",
        },
    },
    "morphine": {
        "extraction": {
            "source": "Papaver somniferum (Schlafmohn) — Latex oder konzentriertes Mohnstroh",
            "source_en": "Papaver somniferum (opium poppy) — latex or concentrated poppy straw",
            "yield_range_percent_w_w": (8, 14),  # of dry latex
            "literature_source": "Kutchan, The Alkaloids 1998 — Morphine biosynthesis and extraction (doi:10.1016/S0099-9598(08)60052-6); Beaudoin & Facchini, Planta 2014 — poppy alkaloid review (doi:10.1007/s00425-014-2056-8)",
            "downstream": {
                "de": [
                    "Wässrige Extraktion bei pH ≈ 4 (Säuerung)",
                    "Filtration und Aktivkohle-Reinigung",
                    "pH-Erhöhung auf ~9 → Fällung der freien Base",
                    "Umkristallisation aus Ethanol/Wasser",
                    "Optional: Salzbildung als Sulfat / Hydrochlorid für API-Grade",
                ],
                "en": [
                    "Aqueous extraction at pH ≈ 4 (acidification)",
                    "Filtration and activated-carbon cleanup",
                    "pH adjustment to ~9 → precipitation of free base",
                    "Recrystallization from ethanol/water",
                    "Optional: salt formation as sulfate / hydrochloride for API grade",
                ],
            },
            "expected_purity_after_workup": "≥ 99 % (USP / Ph. Eur.)",
        },
    },
}


# Generic purification recommendations by molecule type, language-tagged.
TYPE_DOWNSTREAM_HINTS: Dict[str, Dict[str, Any]] = {
    "small_molecule": {
        "literature_source": "Vollhardt & Schore, Organic Chemistry, 8th ed. 2018, Macmillan/Freeman (ISBN 978-1-319-07945-1) — standard workup procedures for small molecules; Clayden et al., Organic Chemistry, 2nd ed. 2012, Oxford Univ. Press (ISBN 978-0-19-927029-3)",
        "default_steps": {
            "de": [
                "Lösungsmittelentfernung (Rotationsverdampfer, 30–50 °C, Vakuum)",
                "Flash-Chromatographie auf SiO2 mit Polaritätsgradient",
                "Optional: Umkristallisation für ≥ 99 % Reinheit",
            ],
            "en": [
                "Solvent removal (rotary evaporator, 30–50 °C, vacuum)",
                "Flash chromatography on SiO2 with polarity gradient",
                "Optional: recrystallization for ≥ 99 % purity",
            ],
        },
        "expected_yield_loss_percent": (10, 25),
    },
    "peptide": {
        "literature_source": "Chan & White, Fmoc Solid Phase Peptide Synthesis: A Practical Approach, Oxford Univ. Press 2000 (ISBN 978-0-19-963724-9); Behrendt et al., J. Pept. Sci. 2016 (review, doi:10.1002/psc.2836)",
        "default_steps": {
            "de": [
                "TFA-Cleavage und Fällung in kaltem Diethylether (bei SPPS)",
                "Präparative RP-HPLC (C18, Wasser/Acetonitril + 0,1 % TFA)",
                "Lyophilisation zur Endformulierung",
            ],
            "en": [
                "TFA cleavage and precipitation in cold diethyl ether (for SPPS)",
                "Preparative RP-HPLC (C18, water/acetonitrile + 0.1 % TFA)",
                "Lyophilization for final formulation",
            ],
        },
        "expected_yield_loss_percent": (20, 50),
    },
    "protein": {
        "literature_source": "Walsh G., Pharmaceutical Biotechnology, 2nd ed. 2018, Wiley (ISBN 978-1-119-11518-7) — antibody/protein downstream processing; Liu et al., mAbs 2010 — Recovery and purification of monoclonal antibodies (doi:10.4161/mabs.2.5.12796)",
        "default_steps": {
            "de": [
                "Zellaufschluss (mechanisch oder Lyse-Puffer)",
                "Klärung via Zentrifugation und Tiefenfiltration",
                "Affinitätschromatographie (Protein A für Antikörper, His-Tag/Ni-NTA für Tag-Proteine)",
                "Polishing via Ionenaustausch- oder Größenausschluss-Chromatographie",
                "Ultrafiltration / Pufferaustausch (TFF)",
            ],
            "en": [
                "Cell lysis (mechanical or lysis buffer)",
                "Clarification via centrifugation and depth filtration",
                "Affinity chromatography (protein A for antibodies, His-tag/Ni-NTA for tagged proteins)",
                "Polishing via ion-exchange or size-exclusion chromatography",
                "Ultrafiltration / buffer exchange (TFF)",
            ],
        },
        "expected_yield_loss_percent": (40, 70),
    },
    "natural_product": {
        "literature_source": "Sarker & Nahar (eds.), Natural Products Isolation, 3rd ed., Methods in Molecular Biology vol. 864, Springer 2012 (ISBN 978-1-61779-623-4, doi:10.1007/978-1-61779-624-1)",
        "default_steps": {
            "de": [
                "Extraktion (Soxhlet oder Mazeration mit Ethanol/Hexan)",
                "Filtration und Konzentration",
                "Säulenchromatographie zur Auftrennung der Sekundärmetaboliten",
            ],
            "en": [
                "Extraction (Soxhlet or maceration with ethanol/hexane)",
                "Filtration and concentration",
                "Column chromatography to separate secondary metabolites",
            ],
        },
        "expected_yield_loss_percent": (30, 60),
    },
}


def _format_range_yield_per_l(rng: Tuple[float, float], lang: str) -> str:
    if lang == "en":
        return f"{rng[0]}–{rng[1]} g/L fermentation broth"
    return f"{rng[0]}–{rng[1]} g/L Fermentationsbrühe"


def _format_range_overall_pct(rng: Tuple[float, float], lang: str) -> str:
    if lang == "en":
        return f"{rng[0]}–{rng[1]} % (overall, all steps)"
    return f"{rng[0]}–{rng[1]} % (über alle Stufen)"


def _format_range_w_w(rng: Tuple[float, float], lang: str) -> str:
    if lang == "en":
        return f"{rng[0]}–{rng[1]} % (w/w of feedstock)"
    return f"{rng[0]}–{rng[1]} % (w/w bezogen auf Rohstoff)"


def _format_range_hours(rng: Tuple[float, float]) -> str:
    return f"{rng[0]}–{rng[1]} h"


def _make_route(hint: Dict[str, Any], lang: str) -> Optional[str]:
    if "organism" in hint:
        organism = hint.get("organism_en") if lang == "en" else hint.get("organism")
        substrate = hint.get("substrate_en") if lang == "en" else hint.get("substrate")
        if substrate:
            connector = " on " if lang == "en" else " auf "
            return f"{organism}{connector}{substrate}"
        return organism
    if "reagents" in hint:
        return hint.get("reagents_en") if lang == "en" else hint.get("reagents")
    if "source" in hint:
        return hint.get("source_en") if lang == "en" else hint.get("source")
    return None


def build_concrete_recommendations(process_input: Dict[str, Any], lang: str = "de") -> Dict[str, Any]:
    """Build a structured dict of concrete recommendations.

    Output keys (text fields are localized to `lang`):
      - target_molecule
      - production_route
      - expected_yield (string with units and range)
      - processing_time (string with range, when known)
      - downstream_steps (list of concrete purification operations)
      - expected_final_purity (string)
      - notes (free-text caveats)
    """
    if lang not in ("de", "en"):
        lang = "de"
    p = process_input or {}
    name = (p.get("molecule_name") or "").strip().lower()
    method = (p.get("method") or "").lower()
    mtype = (p.get("molecule_type") or "small_molecule").lower()

    out: Dict[str, Any] = {
        "target_molecule": p.get("molecule_name") or ("unknown" if lang == "en" else "unbekannt"),
        "production_route": None,
        "expected_yield": None,
        "processing_time": None,
        "downstream_steps": [],
        "expected_final_purity": None,
        "notes": [],
        "literature_sources": [],
    }

    hint = MOLECULE_HINTS.get(name, {}).get(method)
    if hint:
        out["production_route"] = _make_route(hint, lang)
        if hint.get("literature_source"):
            out["literature_sources"].append(hint["literature_source"])

        if "yield_range_g_per_l" in hint:
            out["expected_yield"] = _format_range_yield_per_l(hint["yield_range_g_per_l"], lang)
        elif "yield_range_percent" in hint:
            out["expected_yield"] = _format_range_overall_pct(hint["yield_range_percent"], lang)
        elif "yield_range_percent_w_w" in hint:
            out["expected_yield"] = _format_range_w_w(hint["yield_range_percent_w_w"], lang)

        if "fermentation_time_h" in hint:
            out["processing_time"] = _format_range_hours(hint["fermentation_time_h"])

        downstream = hint.get("downstream") or {}
        out["downstream_steps"] = list(downstream.get(lang) or downstream.get("de") or [])
        out["expected_final_purity"] = hint.get("expected_purity_after_workup")

    # Fallback: type-generic downstream if molecule unknown
    if not out["downstream_steps"]:
        type_hint = TYPE_DOWNSTREAM_HINTS.get(mtype) or TYPE_DOWNSTREAM_HINTS["small_molecule"]
        steps_loc = type_hint["default_steps"]
        out["downstream_steps"] = list(steps_loc.get(lang) or steps_loc.get("de") or [])
        if type_hint.get("literature_source"):
            out["literature_sources"].append(type_hint["literature_source"])
        loss_lo, loss_hi = type_hint["expected_yield_loss_percent"]
        # Localized, human-readable type name (replaces the bare 'protein' / 'peptide' code)
        type_name_map = {
            "de": {
                "small_molecule": "kleine Moleküle",
                "peptide": "Peptide",
                "protein": "Proteine",
                "natural_product": "Naturstoffe",
            },
            "en": {
                "small_molecule": "small molecules",
                "peptide": "peptides",
                "protein": "proteins",
                "natural_product": "natural products",
            },
        }
        type_name = type_name_map.get(lang, type_name_map["de"]).get(mtype, mtype)
        if lang == "en":
            note = f"Typical workup yield loss for {type_name}: {loss_lo}–{loss_hi} %."
        else:
            note = f"Typischer Yield-Verlust durch Aufarbeitung für {type_name}: {loss_lo}–{loss_hi} %."
        out["notes"].append(note)

    # Filter downstream steps by available methods (when user provided them)
    available = p.get("available_methods") or {}
    if available:
        usable, unavailable = _filter_steps_by_capability(out["downstream_steps"], available, lang)
        out["downstream_steps"] = usable
        if unavailable:
            if lang == "en":
                prefix = "The following recommended steps require equipment not marked as available: "
            else:
                prefix = "Folgende empfohlene Schritte erfordern Equipment, das nicht als verfügbar markiert wurde: "
            out["notes"].append(prefix + "; ".join(unavailable))

    return out


# Keyword sets for capability detection — covers DE and EN spellings.
_CAPABILITY_KEYWORDS = [
    ("has_prep_hplc", ("hplc", "präparative", "preparative", "rp-hplc")),
    ("has_flash_chromatography", ("flash", "siliciumdioxid", "sio₂", "säulenchromatograph", "column chromatograph")),
    ("has_lyophilizer", ("lyophilisat", "gefriertrock", "lyophiliz", "freeze-dry")),
    ("has_rotary_evaporator", ("rotationsverdamp", "rotavap", "rotary evapor")),
    ("has_crystallization", ("kristallisation", "umkristallisation", "crystalliz", "recrystalliz")),
    ("has_membrane_filtration", ("ultrafiltration", "tff", "tangentialfluss", "membranfiltration", "membrane filtration", "tangential")),
    ("has_fplc", ("fplc", "affinitätschromatograph", "ionenaustausch", "größenausschluss", "affinity chromatograph", "ion-exchange", "size-exclusion")),
    ("has_extraction", ("extraktion", "soxhlet", "perkolation", "extraction", "percolation")),
]


def _filter_steps_by_capability(steps: List[str], available: Dict[str, bool], lang: str) -> Tuple[List[str], List[str]]:
    usable: List[str] = []
    unavailable: List[str] = []
    for step in steps:
        s_low = step.lower()
        required_caps = [cap for cap, kws in _CAPABILITY_KEYWORDS if any(kw in s_low for kw in kws)]
        if not required_caps:
            usable.append(step)
            continue
        if all(available.get(cap) for cap in required_caps):
            usable.append(step)
        else:
            unavailable.append(step)
    return usable, unavailable
