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
            "literature_source": "Hansen et al., Appl. Environ. Microbiol. 2009 — De novo biosynthesis of vanillin in S. cerevisiae (doi:10.1128/aem.02681-08); Walton et al., Phytochemistry 2003 (review, doi:10.1016/S0031-9422(03)00149-3)",
            "production": {
                "de": [
                    "Stammvorbereitung: rekombinante S. cerevisiae mit eingebrachter Ferulasäure-Decarboxylase- und Vanillin-Synthese-Pathway (qPCR-Verifizierung der Genetik)",
                    "Vorkultur: 30 °C, YPD-Medium (10 g/L Hefeextrakt, 20 g/L Pepton, 20 g/L Glucose), 16-24 h, Ziel-OD600 ~3",
                    "Seed-Train: schrittweise Hochskalierung 1 L → 50 L → 500 L Bioreaktor",
                    "Hauptfermentation: Fed-Batch im Edelstahl-Bioreaktor (1-10 m3), pH 6.5, 28-30 °C, 0.5-2 vvm Belüftung, Ferulasäure-Feed (0.5-2 g/L/h)",
                    "Fermentations-Dauer 48-96 h, Vanillin-Titer am Ende 0.3-1.5 g/L im Medium",
                ],
                "en": [
                    "Strain preparation: recombinant S. cerevisiae carrying ferulate decarboxylase and vanillin pathway (qPCR genotype check)",
                    "Pre-culture: 30 °C, YPD medium (10 g/L yeast extract, 20 g/L peptone, 20 g/L glucose), 16-24 h, target OD600 ~3",
                    "Seed train: stepwise scale-up 1 L → 50 L → 500 L bioreactor",
                    "Main fermentation: fed-batch in stainless-steel bioreactor (1-10 m3), pH 6.5, 28-30 °C, 0.5-2 vvm aeration, ferulic-acid feed (0.5-2 g/L/h)",
                    "Fermentation duration 48-96 h, final vanillin titer 0.3-1.5 g/L in medium",
                ],
            },
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
            "production": {
                "de": [
                    "Edukt-Bereitstellung: Guaiacol (Brenzcatechinmonomethylether) + Glyoxylsäure (40 % wässr.) + NaOH",
                    "Schritt 1 — Kondensation: Guaiacol + Glyoxylsäure unter basischen Bedingungen (NaOH, 15-25 °C, 2-4 h) → Mandelsäure-Zwischenprodukt",
                    "Schritt 2 — Oxidation: Mandelsäure-Intermediat + Luft/Cu-Katalysator oder NaOCl (40-60 °C, 2-3 h) → Vanillinsäure-Übergang",
                    "Schritt 3 — Decarboxylierung: thermisch (140-160 °C) oder Säure-katalysiert → Vanillin",
                    "Roh-Vanillin im Reaktionsgemisch: 60-80 % theoretische Ausbeute, vor Aufarbeitung",
                ],
                "en": [
                    "Reagent supply: guaiacol (catechol monomethyl ether) + glyoxylic acid (40 % aq.) + NaOH",
                    "Step 1 — condensation: guaiacol + glyoxylic acid under basic conditions (NaOH, 15-25 °C, 2-4 h) → mandelic-acid intermediate",
                    "Step 2 — oxidation: mandelic intermediate + air/Cu catalyst or NaOCl (40-60 °C, 2-3 h) → vanillic-acid transition",
                    "Step 3 — decarboxylation: thermal (140-160 °C) or acid-catalysed → vanillin",
                    "Crude vanillin in reaction mixture: 60-80 % theoretical yield before workup",
                ],
            },
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
            "literature_source": "Sinha et al., Crit. Rev. Food Sci. Nutr. 2008 — Vanilla as a flavour ingredient (doi:10.1080/09687630701539350)",
            "production": {
                "de": [
                    "Rohmaterial: kurierten Vanilleschoten (Vanilla planifolia), Trockenmasse, mind. 6-9 Monate fermentiert",
                    "Vorbehandlung: Mahlen oder Schneiden in 2-5 mm Stücke; optional CO2-überkritische Vor-Entfettung",
                    "Vanillin-Gehalt im Rohmaterial: 1.5-2.5 % w/w (Trockenmasse)",
                ],
                "en": [
                    "Raw material: cured vanilla pods (Vanilla planifolia), dry mass, minimum 6-9 months cured",
                    "Pre-treatment: milling or cutting into 2-5 mm pieces; optional supercritical CO2 de-fatting",
                    "Vanillin content in raw material: 1.5-2.5 % w/w (dry mass)",
                ],
            },
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
            "production": {
                "de": [
                    "Edukt-Bereitstellung: Salicylsäure (Pharma-Grade) + Essigsäureanhydrid (≥ 99 %)",
                    "Schritt 1 — Acetylierung: Salicylsäure + Essigsäureanhydrid mit kat. H2SO4 (1-2 mol%), 70-90 °C, 1-2 h",
                    "Reaktionsabbruch: Zugabe kalten Wassers zur Hydrolyse des überschüssigen Anhydrids",
                    "Roh-Acetylsalicylsäure fällt aus dem Reaktionsgemisch aus (75-92 % Ausbeute, vor Aufarbeitung)",
                ],
                "en": [
                    "Reagent supply: salicylic acid (pharma grade) + acetic anhydride (>= 99 %)",
                    "Step 1 — acetylation: salicylic acid + acetic anhydride with cat. H2SO4 (1-2 mol%), 70-90 °C, 1-2 h",
                    "Quench: addition of cold water to hydrolyse excess anhydride",
                    "Crude acetylsalicylic acid precipitates from the reaction mixture (75-92 % yield before workup)",
                ],
            },
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
            "production": {
                "de": [
                    "Edukt-Bereitstellung: Isobutylbenzol (IBB) + Essigsäureanhydrid + HF (Boots) bzw. wasserfreie HF (BHC)",
                    "Schritt 1 — Friedel-Crafts-Acylierung: IBB + Essigsäureanhydrid → 4'-Isobutylacetophenon (Pd-katalysiert bei BHC, AlCl3 bei Boots, ~40-80 °C)",
                    "Schritt 2 — Hydrierung: Acetophenon-Intermediat → Alkohol (Raney-Ni oder Pd/C, 60-100 °C, 3-10 bar H2)",
                    "Schritt 3 — Carbonylierung: Alkohol + CO + Pd-Katalysator → Ibuprofen (BHC, 3-Stufen-Prozess); klassisch Boots benötigt 6 Stufen",
                    "Roh-Ibuprofen 50-75 % Gesamtausbeute, vor Aufarbeitung",
                ],
                "en": [
                    "Reagent supply: isobutylbenzene (IBB) + acetic anhydride + HF (Boots) or anhydrous HF (BHC)",
                    "Step 1 — Friedel-Crafts acylation: IBB + acetic anhydride -> 4'-isobutylacetophenone (Pd-catalysed in BHC, AlCl3 in Boots, ~40-80 °C)",
                    "Step 2 — hydrogenation: acetophenone intermediate -> alcohol (Raney-Ni or Pd/C, 60-100 °C, 3-10 bar H2)",
                    "Step 3 — carbonylation: alcohol + CO + Pd catalyst -> ibuprofen (BHC, 3-step process); classical Boots needs 6 steps",
                    "Crude ibuprofen 50-75 % overall yield, before workup",
                ],
            },
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
            "production": {
                "de": [
                    "Rohstoff-Vorbereitung: Stärke-Hydrolyse (alpha-Amylase + Glucoamylase, 60-95 °C) zu fermentierbarer Glucose",
                    "Hefe-Stamm-Anzucht: S. cerevisiae (Bäcker- oder Brennereihefe), Vorkultur 30 °C, YPD-Medium, 12-16 h",
                    "Hauptfermentation: Großmaßstab-Fermenter (100-1000 m3), 28-32 °C, pH 4-5, anaerob",
                    "Fermentations-Dauer 48-72 h, Endkonzentration 8-15 % v/v Ethanol (80-130 g/L)",
                    "Maische enthält Ethanol + Reste an Zucker, Hefe-Zellen, Nebenprodukten (Fuselöle, Acetaldehyd)",
                ],
                "en": [
                    "Feedstock preparation: starch hydrolysis (alpha-amylase + glucoamylase, 60-95 °C) to fermentable glucose",
                    "Yeast strain prep: S. cerevisiae (baker's or distiller's yeast), pre-culture 30 °C, YPD medium, 12-16 h",
                    "Main fermentation: large-scale fermenter (100-1000 m3), 28-32 °C, pH 4-5, anaerobic",
                    "Fermentation duration 48-72 h, final concentration 8-15 % v/v ethanol (80-130 g/L)",
                    "Mash contains ethanol + residual sugar, yeast cells, side products (fusel oils, acetaldehyde)",
                ],
            },
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
            "production": {
                "de": [
                    "Stammvorbereitung: Hochertrags-Stamm von P. chrysogenum (Wisconsin-Linie oder neuere Hochleister)",
                    "Vorkultur: Schüttelkolben oder kleiner Bioreaktor, 25-27 °C, 24-48 h",
                    "Seed-Train: 100 L → 10 m3 → finaler Produktions-Fermenter 50-200 m3",
                    "Hauptfermentation: Fed-Batch, 25-27 °C, pH 6.5-7.0, hohe Belüftung (1 vvm) + Rühren (150-300 rpm), Lactose-Feed + Phenylessigsäure-Precursor-Feed",
                    "Fermentations-Dauer 120-200 h, finale Titer 5-50 g/L Penicillin G",
                ],
                "en": [
                    "Strain preparation: high-yield P. chrysogenum strain (Wisconsin lineage or newer high-producers)",
                    "Pre-culture: shake flask or small bioreactor, 25-27 °C, 24-48 h",
                    "Seed train: 100 L -> 10 m3 -> final production fermenter 50-200 m3",
                    "Main fermentation: fed-batch, 25-27 °C, pH 6.5-7.0, high aeration (1 vvm) + agitation (150-300 rpm), lactose feed + phenylacetic acid precursor feed",
                    "Fermentation duration 120-200 h, final titer 5-50 g/L penicillin G",
                ],
            },
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
            "literature_source": "Berovic & Legisa, Biotechnol. Annu. Rev. 2007 — Citric acid production (doi:10.1016/S1387-2656(07)13011-8); Show et al., Bioengineered 2015",
            "production": {
                "de": [
                    "Stammvorbereitung: A. niger Hochertrags-Stamm aus Stammbank, Sporulation auf Schrägagar",
                    "Vorkultur: Sporensuspension in flüssiges Medium, 28-30 °C, 24-48 h Vorwachstum",
                    "Hauptfermentation: Submerskultur in 100-500 m3 Bioreaktor, 28-30 °C, pH 1.6-2.2 (sehr sauer für Citrat-Akkumulation!), hoher Sauerstoffeintrag (>0.5 vvm)",
                    "Substrat: Melasse (oder Saccharose-Sirup), Fe-Limitierung kritisch (< 0.5 ppm) für hohe Ausbeute",
                    "Fermentations-Dauer 120-200 h, finale Citronensäure-Konzentration 150-200 g/L",
                ],
                "en": [
                    "Strain preparation: A. niger high-yield strain from strain bank, sporulation on slant agar",
                    "Pre-culture: spore suspension into liquid medium, 28-30 °C, 24-48 h pre-growth",
                    "Main fermentation: submerged culture in 100-500 m3 bioreactor, 28-30 °C, pH 1.6-2.2 (very acidic for citrate accumulation!), high oxygen transfer (>0.5 vvm)",
                    "Substrate: molasses (or sucrose syrup), Fe limitation critical (< 0.5 ppm) for high yield",
                    "Fermentation duration 120-200 h, final citric-acid concentration 150-200 g/L",
                ],
            },
            "downstream": {
                "de": [
                    "Myzel-Abtrennung via Filtration",
                    "Fällung als Calciumcitrat mit Ca(OH)2",
                    "Zersetzung mit H2SO4 → Calciumsulfat-Niederschlag + freie Citronensäure",
                    "Aktivkohle-Behandlung, Eindampfung, Kristallisation",
                ],
                "en": [
                    "Mycelium separation via filtration",
                    "Precipitation as calcium citrate with Ca(OH)2",
                    "Decomposition with H2SO4 -> calcium sulfate solid + free citric acid",
                    "Activated-carbon polish, evaporation, crystallization",
                ],
            },
            "expected_purity_after_workup": ">= 99,5 % (food grade)",
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
            "production": {
                "de": [
                    "Stammvorbereitung: Lactobacillus delbrueckii (L-Milchsäure) oder Bacillus coagulans (thermophil, höhere Temperatur)",
                    "Vorkultur: MRS-Medium oder Glucose/Hefeextrakt, 37-50 °C (stammabhängig), 12-24 h",
                    "Hauptfermentation: Batch oder Fed-Batch in 10-500 m3 Bioreaktor, 37-50 °C, pH-Kontrolle bei 5.5-6.5 mittels Ca(OH)2- oder NH4OH-Zugabe (Milchsäure wird sofort neutralisiert!)",
                    "Substrat: Glucose oder Stärke-Hydrolysat (100-200 g/L Startkonzentration)",
                    "Fermentations-Dauer 24-72 h, Endkonzentration 90-180 g/L Lactat-Salz",
                ],
                "en": [
                    "Strain preparation: Lactobacillus delbrueckii (L-lactic acid) or Bacillus coagulans (thermophilic, higher temperature)",
                    "Pre-culture: MRS medium or glucose/yeast extract, 37-50 °C (strain-dependent), 12-24 h",
                    "Main fermentation: batch or fed-batch in 10-500 m3 bioreactor, 37-50 °C, pH control at 5.5-6.5 via Ca(OH)2 or NH4OH addition (lactic acid is immediately neutralised!)",
                    "Substrate: glucose or starch hydrolysate (100-200 g/L starting concentration)",
                    "Fermentation duration 24-72 h, final concentration 90-180 g/L lactate salt",
                ],
            },
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
            "production": {
                "de": [
                    "Stammvorbereitung: A. succinogenes oder B. succiniciproducens (gentechnisch optimiert für hohe Ausbeute)",
                    "Vorkultur: Komplexes Medium (Hefeextrakt + Glucose + MgCO3 als CO2-Quelle), 37 °C, 16-24 h, anaerob/mikroaerob",
                    "Hauptfermentation: Fed-Batch in 10-100 m3 Bioreaktor, 37 °C, pH 6.5-7.0 (NH4OH-Titration), CO2-Begasung kontinuierlich",
                    "Glucose-Feed-Kontrolle: 30-60 g/L im Reaktor, gesamte Substrat-Belastung 100-150 g/L",
                    "Fermentations-Dauer 40-80 h, finale Bernsteinsäure-Konzentration 50-100 g/L",
                ],
                "en": [
                    "Strain preparation: A. succinogenes or B. succiniciproducens (genetically optimised for high yield)",
                    "Pre-culture: complex medium (yeast extract + glucose + MgCO3 as CO2 source), 37 °C, 16-24 h, anaerobic/microaerobic",
                    "Main fermentation: fed-batch in 10-100 m3 bioreactor, 37 °C, pH 6.5-7.0 (NH4OH titration), continuous CO2 sparging",
                    "Glucose feed control: 30-60 g/L in reactor, total substrate load 100-150 g/L",
                    "Fermentation duration 40-80 h, final succinic-acid concentration 50-100 g/L",
                ],
            },
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
            "literature_source": "Joncour et al., Green Chem. 2014 — Paracetamol synthesis review (doi:10.1039/c4gc00166d); Ellis, Paracetamol: A Curriculum Resource, RSC 2002 (ISBN 978-0-85404-365-1)",
            "production": {
                "de": [
                    "Edukt-Bereitstellung: p-Aminophenol (Pharma-Grade, > 99 %) + Essigsäureanhydrid (≥ 99 %)",
                    "Schritt 1 — Acetylierung: p-Aminophenol + Essigsäureanhydrid in wässriger Suspension, 80-90 °C, 30-60 min (selektive N-Acetylierung, OH-Gruppe bleibt frei)",
                    "Reaktion in Wasser ohne organische Lösungsmittel (green-chemistry-konform)",
                    "Roh-Paracetamol fällt direkt als Feststoff aus (85-95 % Ausbeute)",
                ],
                "en": [
                    "Reagent supply: p-aminophenol (pharma grade, > 99 %) + acetic anhydride (>= 99 %)",
                    "Step 1 — acetylation: p-aminophenol + acetic anhydride in aqueous suspension, 80-90 °C, 30-60 min (selective N-acetylation, OH group remains free)",
                    "Reaction in water without organic solvents (green-chemistry compliant)",
                    "Crude paracetamol precipitates directly as solid (85-95 % yield)",
                ],
            },
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
            "literature_source": "Norouzian et al., Biotechnol. Adv. 2006 — Fungal glucoamylase production (doi:10.1016/j.biotechadv.2005.06.003); Kumar & Satyanarayana, Crit. Rev. Biotechnol. 2009",
            "production": {
                "de": [
                    "Stammvorbereitung: A. niger Industrial-Stamm (z. B. NRRL 3122 oder kommerzielle Hochleister von Novozymes/DSM)",
                    "Vorkultur: Schüttelkolben oder kleiner Bioreaktor, 28-30 °C, 48-72 h auf Stärke- oder Glucose-Medium",
                    "Hauptfermentation: Submerskultur in 100-500 m3 Bioreaktor, 28-30 °C, pH 4.5-5.5, hoher Sauerstoffeintrag (>0.5 vvm)",
                    "Substrat: Stärke-Hydrolysat oder Glucose + Maisquellwasser als Stickstoffquelle",
                    "Fermentations-Dauer 120-168 h, Enzym-Titer 5-20 g/L im Kulturüberstand (extrazelluläre Sekretion)",
                ],
                "en": [
                    "Strain preparation: A. niger industrial strain (e.g. NRRL 3122 or commercial high-producers from Novozymes/DSM)",
                    "Pre-culture: shake flask or small bioreactor, 28-30 °C, 48-72 h on starch or glucose medium",
                    "Main fermentation: submerged culture in 100-500 m3 bioreactor, 28-30 °C, pH 4.5-5.5, high oxygen transfer (>0.5 vvm)",
                    "Substrate: starch hydrolysate or glucose + corn-steep liquor as nitrogen source",
                    "Fermentation duration 120-168 h, enzyme titer 5-20 g/L in culture supernatant (extracellular secretion)",
                ],
            },
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
            "literature_source": "Bhosale et al., Microbiol. Rev. 1996 — Glucose isomerase comprehensive review (doi:10.1128/mr.60.2.280-300.1996); Crabb & Mitchinson, Trends Biotechnol. 1997",
            "production": {
                "de": [
                    "Stammvorbereitung: Streptomyces rubiginosus oder S. murinus, Industrial-Stamm aus Stammbank",
                    "Vorkultur: Komplexes Medium (Hefeextrakt + Sojamehl + Glucose + Mg2+), 28-30 °C, 24-48 h",
                    "Hauptfermentation: Submerskultur in 50-200 m3 Bioreaktor, 28-30 °C, pH 6.8-7.2, hoher O2-Eintrag",
                    "Induktion: Xylose oder Xylan im Medium aktiviert Enzym-Expression (intrazellulär)",
                    "Fermentations-Dauer 48-96 h, Enzym-Aktivität 1-5 g/L (intrazellulär gespeichert)",
                ],
                "en": [
                    "Strain preparation: Streptomyces rubiginosus or S. murinus, industrial strain from strain bank",
                    "Pre-culture: complex medium (yeast extract + soy meal + glucose + Mg2+), 28-30 °C, 24-48 h",
                    "Main fermentation: submerged culture in 50-200 m3 bioreactor, 28-30 °C, pH 6.8-7.2, high O2 transfer",
                    "Induction: xylose or xylan in medium activates enzyme expression (intracellular)",
                    "Fermentation duration 48-96 h, enzyme activity 1-5 g/L (stored intracellularly)",
                ],
            },
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
    # ---- Recombinant therapeutic proteins (non-mAb) ----
    "insulin": {
        "biotechnological": {
            "organism": "E. coli (Refolding) oder Pichia pastoris (sekretorisch)",
            "organism_en": "E. coli (refolding) or Pichia pastoris (secretory)",
            "substrate": "Glucose-Medium + Induktor (IPTG für E. coli, Methanol für P. pastoris)",
            "substrate_en": "glucose medium + inducer (IPTG for E. coli, methanol for P. pastoris)",
            "yield_range_g_per_l": (0.5, 3.0),
            "fermentation_time_h": (48, 96),
            "literature_source": "Walsh G., Eur. J. Pharm. Biopharm. 2005 — Therapeutic insulins and their large-scale manufacture (doi:10.1007/s00253-004-1809-x); Baeshen et al., J. Microbiol. Biotechnol. 2014 — Cell factories for insulin production (doi:10.1007/s00253-004-1809-x)",
            "production": {
                "de": [
                    "Stammvorbereitung: rekombinante E. coli K12 mit Insulin-Vorläufer-Gen (Proinsulin oder mini-Proinsulin), GMP-Stammbank",
                    "Vorkultur: definiertes Minimalmedium + Selektionsantibiotikum, 37 °C, 12-16 h",
                    "Hauptfermentation: Fed-Batch im 1-50 m3 Bioreaktor, 37 °C, pH 7.0, hoher Glucose-Feed, Induktion mit IPTG nach Erreichen hoher Zelldichte",
                    "Expression als Einschlusskörper (E. coli) — typische Akkumulation 20-30 % des Gesamtproteins",
                    "Zell-Ernte nach 48-96 h via Hochleistungszentrifugation oder Tiefenfiltration",
                ],
                "en": [
                    "Strain preparation: recombinant E. coli K12 with insulin-precursor gene (proinsulin or mini-proinsulin), GMP strain bank",
                    "Pre-culture: defined minimal medium + selection antibiotic, 37 °C, 12-16 h",
                    "Main fermentation: fed-batch in 1-50 m3 bioreactor, 37 °C, pH 7.0, high glucose feed, induction with IPTG after high cell density",
                    "Expression as inclusion bodies (E. coli) — typical accumulation 20-30 % of total protein",
                    "Cell harvest after 48-96 h via high-performance centrifugation or depth filtration",
                ],
            },
            "downstream": {
                "de": [
                    "Zell-Lyse + Inclusion-Body-Isolation via Hochdruckhomogenisator",
                    "Solubilisierung der Inclusion Bodies (8 M Harnstoff oder 6 M Guanidinium-HCl)",
                    "Refolding bei niedriger Proteinkonzentration (0.1-0.5 g/L), Redox-Puffer (GSH/GSSG), pH 10-11, 4-12 °C, 12-24 h",
                    "Enzymatische Spaltung des Fusions-Tags (Trypsin + Carboxypeptidase B) → reifes Insulin",
                    "Reverse-Phase-HPLC (Pharma-Grade C18) — Polishing-Schritt",
                    "Ultrafiltration + Lyophilisation oder Direkt-Formulierung",
                ],
                "en": [
                    "Cell lysis + inclusion-body isolation via high-pressure homogenizer",
                    "Solubilisation of inclusion bodies (8 M urea or 6 M guanidinium-HCl)",
                    "Refolding at low protein concentration (0.1-0.5 g/L), redox buffer (GSH/GSSG), pH 10-11, 4-12 °C, 12-24 h",
                    "Enzymatic cleavage of fusion tag (trypsin + carboxypeptidase B) → mature insulin",
                    "Reverse-phase HPLC (pharma-grade C18) — polishing step",
                    "Ultrafiltration + lyophilisation or direct formulation",
                ],
            },
            "expected_purity_after_workup": "≥ 99,5 % (Ph. Eur. / USP human insulin)",
        },
    },
    "erythropoietin": {
        "biotechnological": {
            "organism": "CHO-Zellen (Suspensionskultur)",
            "organism_en": "CHO cells (suspension culture)",
            "substrate": "Chemisch definiertes Medium (CDM)",
            "substrate_en": "chemically defined medium (CDM)",
            "yield_range_g_per_l": (0.1, 0.5),
            "fermentation_time_h": (240, 336),  # 10-14 days
            "literature_source": "Jelkmann W., Eur. J. Haematol. 2007 — Recombinant EPO (doi:10.1093/ndt/gfm392); Walsh G., Nat. Biotechnol. 2018 — Biopharmaceutical benchmarks (doi:10.1038/nbt.4305)",
            "production": {
                "de": [
                    "Zellbank: rekombinante CHO-Zellen mit EPO-Gen + glykosylierungs-relevanten Genen, GMP-MCB/WCB",
                    "Auftauen + Vorkultur: WCB-Ampulle, 37 °C/5 % CO2/CDM, Wachstum in Schüttelfläschchen 3-5 Tage",
                    "Seed-Train: 50 mL → 2 L → 50 L → 500 L → 2000 L Bioreaktor",
                    "Hauptproduktion: Fed-Batch im Produktions-Bioreaktor, 37 °C, pH 7.0, DO 30 %, Glucose- + Glutamin-Feed über 10-14 Tage",
                    "Kritisch: hochkonsistente Sialylierung der N-Glykane → Halbwertszeit + Bioaktivität (Darbepoetin alfa hat zusätzliche Glykosylierungs-Sites für längere Wirkdauer)",
                    "Produktions-Dauer 10-14 Tage, EPO-Titer 100-500 mg/L im Überstand",
                ],
                "en": [
                    "Cell bank: recombinant CHO cells with EPO gene + glycosylation-relevant genes, GMP MCB/WCB",
                    "Thaw + pre-culture: WCB vial, 37 °C/5 % CO2/CDM, growth in shake flasks 3-5 days",
                    "Seed train: 50 mL -> 2 L -> 50 L -> 500 L -> 2000 L bioreactor",
                    "Main production: fed-batch in production bioreactor, 37 °C, pH 7.0, DO 30 %, glucose + glutamine feed over 10-14 days",
                    "Critical: highly consistent sialylation of N-glycans -> half-life + bioactivity (darbepoetin alfa has extra glycosylation sites for prolonged action)",
                    "Production duration 10-14 days, EPO titer 100-500 mg/L in supernatant",
                ],
            },
            "downstream": {
                "de": [
                    "Klärung via Tiefenfiltration + 0,2 µm Sterilfiltration (kein Protein-A, da keine Fc-Region)",
                    "Anionenaustausch-Chromatographie (Q-Sepharose) — Capture-Schritt",
                    "Reverse-Phase-HPLC zur Trennung von Glykoformen — kritisch für Sialylierungs-Konsistenz",
                    "Hydroxyapatit-Chromatographie — Polishing",
                    "Ultrafiltration / Diafiltration zur Endkonzentration",
                    "Optional: Größenausschluss-Chromatographie als finales Polishing",
                ],
                "en": [
                    "Clarification via depth filtration + 0.2 µm sterile filtration (no Protein A — no Fc region)",
                    "Anion-exchange chromatography (Q-Sepharose) — capture step",
                    "Reverse-phase HPLC for glycoform separation — critical for sialylation consistency",
                    "Hydroxyapatite chromatography — polishing",
                    "Ultrafiltration / diafiltration to final concentration",
                    "Optional: size-exclusion chromatography as final polishing",
                ],
            },
            "expected_purity_after_workup": "≥ 99 % (Ph. Eur. / USP); Glykoform-Profil definierte Spezifikation",
        },
    },
    "filgrastim": {
        "biotechnological": {
            "organism": "E. coli (BL21 oder K12)",
            "organism_en": "E. coli (BL21 or K12)",
            "substrate": "Definiertes Minimalmedium + Glucose",
            "substrate_en": "defined minimal medium + glucose",
            "yield_range_g_per_l": (1, 5),
            "fermentation_time_h": (24, 48),
            "literature_source": "Dasari et al., Process Biochem. 2008 — Optimization of GCSF production in E. coli (doi:10.1016/j.procbio.2008.02.012); Welte K., Cell Death Differ. 2014 (G-CSF review, doi:10.14260/jemds/2015/1115)",
            "production": {
                "de": [
                    "Stammvorbereitung: rekombinante E. coli BL21(DE3) mit G-CSF-Gen unter T7-Promotor",
                    "Vorkultur: LB- oder definiertes Minimalmedium + Kanamycin, 37 °C, 12 h",
                    "Hauptfermentation: Fed-Batch im 1-20 m3 Bioreaktor, 37 °C → Temperatur-Shift auf 30 °C bei Induktion, pH 7.0, hoher O2-Eintrag",
                    "Induktion mit IPTG (0.1-1 mM) bei OD600 ~30-50, Expression 4-6 h",
                    "Filgrastim akkumuliert als Einschlusskörper (kein N-terminales Methionin im klassischen Filgrastim — wichtig für Aktivität)",
                    "Ernte nach 24-48 h, Titer 1-5 g/L im IB-Anteil",
                ],
                "en": [
                    "Strain preparation: recombinant E. coli BL21(DE3) with G-CSF gene under T7 promoter",
                    "Pre-culture: LB or defined minimal medium + kanamycin, 37 °C, 12 h",
                    "Main fermentation: fed-batch in 1-20 m3 bioreactor, 37 °C -> temp shift to 30 °C at induction, pH 7.0, high O2 transfer",
                    "Induction with IPTG (0.1-1 mM) at OD600 ~30-50, expression 4-6 h",
                    "Filgrastim accumulates as inclusion bodies (no N-terminal methionine in classical filgrastim — important for activity)",
                    "Harvest after 24-48 h, titer 1-5 g/L in IB fraction",
                ],
            },
            "downstream": {
                "de": [
                    "Zell-Lyse + Inclusion-Body-Isolation (Hochdruckhomogenisator)",
                    "Solubilisierung in 8 M Harnstoff + Reduktionsmittel (DTT)",
                    "Refolding bei niedriger Konzentration (0.1-0.3 g/L), 4-12 °C, 12-24 h",
                    "Kationenaustausch-Chromatographie (SP-Sepharose) — Capture",
                    "Reverse-Phase-HPLC — Polishing",
                    "Ultrafiltration / Diafiltration in Formulierungs-Puffer",
                ],
                "en": [
                    "Cell lysis + inclusion-body isolation (high-pressure homogenizer)",
                    "Solubilisation in 8 M urea + reducing agent (DTT)",
                    "Refolding at low concentration (0.1-0.3 g/L), 4-12 °C, 12-24 h",
                    "Cation-exchange chromatography (SP-Sepharose) — capture",
                    "Reverse-phase HPLC — polishing",
                    "Ultrafiltration / diafiltration into formulation buffer",
                ],
            },
            "expected_purity_after_workup": "≥ 99 % (Ph. Eur. / USP); monomere Form ≥ 98 %",
        },
    },
    "etanercept": {
        "biotechnological": {
            "organism": "CHO-Zellen (Suspensionskultur)",
            "organism_en": "CHO cells (suspension culture)",
            "substrate": "Chemisch definiertes Medium (CDM)",
            "substrate_en": "chemically defined medium (CDM)",
            "yield_range_g_per_l": (1, 4),
            "fermentation_time_h": (288, 360),  # 12-15 days
            "literature_source": "Mohler et al., J. Immunol. 1993 — TNFR-Fc-Fusionsprotein (Enbrel-Erstbeschreibung); Walsh G., Nat. Biotechnol. 2018 — Biopharmaceutical benchmarks (doi:10.1038/nbt.4305)",
            "production": {
                "de": [
                    "Zellbank: rekombinante CHO-Zellen mit Etanercept-Gen (TNFR2 extrazelluläre Domäne + IgG1-Fc-Fusion), GMP-MCB/WCB",
                    "Auftauen + Vorkultur: WCB-Ampulle, 37 °C/5 % CO2/CDM, Schüttelfläschchen 3-5 Tage",
                    "Seed-Train: 50 mL → 2 L → 50 L → 500 L → 2000-5000 L Edelstahl-Bioreaktor",
                    "Hauptproduktion: Fed-Batch 12-15 Tage, 37 °C, pH 7.0, DO 30 %, Glucose- + Glutamin-Feed",
                    "Glykosylierung der TNFR-Domäne ist kritisch für Bindungsaffinität — CHO-typische N-Glykane sind notwendig",
                    "Produktions-Ende: Titer 1-4 g/L im Überstand",
                ],
                "en": [
                    "Cell bank: recombinant CHO cells with etanercept gene (TNFR2 extracellular domain + IgG1-Fc fusion), GMP MCB/WCB",
                    "Thaw + pre-culture: WCB vial, 37 °C/5 % CO2/CDM, shake flasks 3-5 days",
                    "Seed train: 50 mL -> 2 L -> 50 L -> 500 L -> 2000-5000 L stainless-steel bioreactor",
                    "Main production: fed-batch 12-15 days, 37 °C, pH 7.0, DO 30 %, glucose + glutamine feed",
                    "Glycosylation of TNFR domain critical for binding affinity — CHO-typical N-glycans necessary",
                    "End of production: titer 1-4 g/L in supernatant",
                ],
            },
            "downstream": {
                "de": [
                    "Klärung via Tiefenfiltration + 0,2 µm Sterilfiltration",
                    "Protein-A-Affinitätschromatographie (Capture via Fc-Region; > 95 % Reinheit ab Stufe 1)",
                    "Niedrig-pH-Virusinaktivierung (pH 3,5 für 60 min)",
                    "Hydrophobe Interaktions-Chromatographie (HIC) — Dimer/Aggregat-Trennung",
                    "Anionenaustausch im Flow-through (AEX) + Virusfiltration",
                    "Ultrafiltration / Pufferaustausch (TFF, finale Formulierung)",
                ],
                "en": [
                    "Clarification via depth filtration + 0.2 µm sterile filtration",
                    "Protein-A affinity chromatography (capture via Fc region; > 95 % purity from step 1)",
                    "Low-pH virus inactivation (pH 3.5 for 60 min)",
                    "Hydrophobic-interaction chromatography (HIC) — dimer/aggregate separation",
                    "Anion-exchange in flow-through (AEX) + virus filtration",
                    "Ultrafiltration / buffer exchange (TFF, final formulation)",
                ],
            },
            "expected_purity_after_workup": "≥ 99 % (cGMP); Dimer/Aggregat ≤ 1 %",
        },
    },
    "somatropin": {
        "biotechnological": {
            "organism": "E. coli (periplasmatisch oder zytoplasmatisch)",
            "organism_en": "E. coli (periplasmic or cytoplasmic)",
            "substrate": "Definiertes Minimalmedium + Glucose",
            "substrate_en": "defined minimal medium + glucose",
            "yield_range_g_per_l": (0.5, 3),
            "fermentation_time_h": (24, 48),
            "literature_source": "Crommelin et al., Pharmaceutical Biotechnology, 5th ed. 2019 (ISBN 978-3-030-00710-2) — recombinant human growth hormone manufacturing; Singh & Panda, J. Biosci. Bioeng. 2005",
            "production": {
                "de": [
                    "Stammvorbereitung: rekombinante E. coli mit hGH-Gen + Signal-Sequenz (periplasmatische Sekretion) ODER ohne (zytoplasmatisch, Inclusion Bodies)",
                    "Vorkultur: definiertes Minimalmedium + Selektion, 37 °C, 12 h",
                    "Hauptfermentation: Fed-Batch im 1-50 m3 Bioreaktor, 37 °C, pH 7.0, hoher Glucose-Feed",
                    "Induktion mit IPTG oder Temperatur-Shift, Expression 4-8 h",
                    "Bei periplasmatischer Strategie: hGH wird sekretiert + bildet korrekte Disulfid-Brücken in situ (kein Refolding nötig)",
                    "Bei zytoplasmatischer Strategie: Inclusion Bodies, Refolding erforderlich",
                    "Ernte nach 24-48 h, Titer 0.5-3 g/L",
                ],
                "en": [
                    "Strain preparation: recombinant E. coli with hGH gene + signal sequence (periplasmic secretion) OR without (cytoplasmic, inclusion bodies)",
                    "Pre-culture: defined minimal medium + selection, 37 °C, 12 h",
                    "Main fermentation: fed-batch in 1-50 m3 bioreactor, 37 °C, pH 7.0, high glucose feed",
                    "Induction with IPTG or temperature shift, expression 4-8 h",
                    "Periplasmic strategy: hGH is secreted + forms correct disulfide bridges in situ (no refolding needed)",
                    "Cytoplasmic strategy: inclusion bodies, refolding required",
                    "Harvest after 24-48 h, titer 0.5-3 g/L",
                ],
            },
            "downstream": {
                "de": [
                    "Periplasmatische Extraktion (osmotischer Schock) oder Zell-Lyse + IB-Isolation (zytoplasmatisch)",
                    "Bei IB-Route: Solubilisierung + Refolding in Harnstoff/Redox-Puffer",
                    "Anionenaustausch-Chromatographie (Q-Sepharose) — Capture",
                    "Hydrophobe Interaktions-Chromatographie (HIC) — Polishing",
                    "Größenausschluss-Chromatographie als finales Polishing",
                    "Ultrafiltration / Lyophilisation",
                ],
                "en": [
                    "Periplasmic extraction (osmotic shock) or cell lysis + IB isolation (cytoplasmic)",
                    "IB route: solubilisation + refolding in urea/redox buffer",
                    "Anion-exchange chromatography (Q-Sepharose) — capture",
                    "Hydrophobic-interaction chromatography (HIC) — polishing",
                    "Size-exclusion chromatography as final polishing",
                    "Ultrafiltration / lyophilisation",
                ],
            },
            "expected_purity_after_workup": "≥ 99 % (Ph. Eur. / USP)",
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
            "literature_source": "Kotzia & Labrou, J. Biotechnol. 2007 — Microbial L-asparaginase production; Cachumba et al., Braz. J. Microbiol. 2016 (doi:10.1016/j.bjm.2016.10.004)",
            "production": {
                "de": [
                    "Stammvorbereitung: E. coli K-12 mit rekombinantem L-Asparaginase-II-Gen (typ. ansB), aus GMP-Stammbank",
                    "Vorkultur: LB- oder definiertes Minimalmedium, 37 °C, 12-16 h, mit Selektions-Antibiotikum (z. B. Kanamycin)",
                    "Hauptfermentation: Fed-Batch im 1-50 m3 Bioreaktor, 37 °C, pH 7.0, hoher O2-Eintrag, Glucose-Feed",
                    "Induktion: IPTG oder Temperatur-Shift (induzierbare Promotoren) für 4-8 h zur Enzym-Expression",
                    "Fermentations-Dauer 24-48 h, Asparaginase-Akkumulation periplasmatisch oder intrazellulär, 0.5-2 g/L",
                ],
                "en": [
                    "Strain preparation: E. coli K-12 with recombinant L-asparaginase II gene (typ. ansB), from GMP strain bank",
                    "Pre-culture: LB or defined minimal medium, 37 °C, 12-16 h, with selection antibiotic (e.g. kanamycin)",
                    "Main fermentation: fed-batch in 1-50 m3 bioreactor, 37 °C, pH 7.0, high O2 transfer, glucose feed",
                    "Induction: IPTG or temperature shift (inducible promoters) for 4-8 h for enzyme expression",
                    "Fermentation duration 24-48 h, asparaginase accumulation periplasmic or intracellular, 0.5-2 g/L",
                ],
            },
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
            "literature_source": "Liu et al., mAbs 2010 — Antibody recovery and purification (doi:10.4161/mabs.2.5.12645); Walsh G., Nat. Biotechnol. 2018 — Biopharmaceutical benchmarks (doi:10.1038/nbt.4305)",
            "production": {
                "de": [
                    "Zellbank: Master- + Working Cell Bank (MCB/WCB) von rekombinanten CHO-Zellen mit anti-PD1-Antikörper-Gen, GMP-konform charakterisiert",
                    "Auftauen + Vorkultur: WCB-Ampulle auftauen, 37 °C/5 % CO2/serumfreies CDM-Medium, Wachstum in Schüttelfläschchen 3-5 Tage",
                    "Seed-Train: schrittweise Expansion 50 mL → 2 L → 50 L → 500 L → 2000-5000 L Edelstahl-Bioreaktor",
                    "Hauptproduktion: Fed-Batch im Produktions-Bioreaktor, 37 °C, pH 7.0, DO > 30 %, Glucose- und Glutamin-Feed-Strategie",
                    "Produktions-Dauer 12-15 Tage, Antikörper-Titer 3-8 g/L im Kulturüberstand",
                ],
                "en": [
                    "Cell bank: master + working cell bank (MCB/WCB) of recombinant CHO cells with anti-PD1 antibody gene, GMP-characterised",
                    "Thaw + pre-culture: thaw WCB vial, 37 °C/5 % CO2/serum-free CDM medium, growth in shake flasks 3-5 days",
                    "Seed train: stepwise expansion 50 mL -> 2 L -> 50 L -> 500 L -> 2000-5000 L stainless-steel bioreactor",
                    "Main production: fed-batch in production bioreactor, 37 °C, pH 7.0, DO > 30 %, glucose and glutamine feed strategy",
                    "Production duration 12-15 days, antibody titer 3-8 g/L in culture supernatant",
                ],
            },
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
            "literature_source": "Liu et al., mAbs 2010 — Antibody recovery and purification (doi:10.4161/mabs.2.5.12645); Walsh G., Nat. Biotechnol. 2018 — Biopharmaceutical benchmarks (doi:10.1038/nbt.4305); Wurm FM., Nat. Biotechnol. 2004 — Production of recombinant protein therapeutics in CHO (doi:10.4161/mabs.2.5.12645)",
            "production": {
                "de": [
                    "Zellbank: GMP-MCB/WCB von rekombinanten CHO-Zellen mit Trastuzumab-Schwer- und -Leichtketten-Genen, Glykosylierungs-Konsistenz validiert",
                    "Auftauen + Vorkultur: WCB-Ampulle, 37 °C/5 % CO2/CDM-Medium, Schüttelfläschchen 3-5 Tage bis ~2-5 × 10^6 Zellen/mL",
                    "Seed-Train: 50 mL → 2 L Wave-Bag → 50 L Single-Use-Bioreaktor → 500 L → 2000-5000 L Edelstahl-Bioreaktor",
                    "Hauptproduktion: Fed-Batch 12-15 Tage, 37 °C → Temp-Shift auf 32-34 °C nach Tag 4-5 für Produktivitäts-Boost, pH 7.0, DO 30 %, Glucose- + Aminosäure-Feeds",
                    "Produktions-Ende: Antikörper-Titer 3-10 g/L im Überstand, Zell-Viabilität noch > 70 %",
                ],
                "en": [
                    "Cell bank: GMP MCB/WCB of recombinant CHO cells with trastuzumab heavy- and light-chain genes, glycosylation consistency validated",
                    "Thaw + pre-culture: WCB vial, 37 °C/5 % CO2/CDM medium, shake flasks 3-5 days to ~2-5 x 10^6 cells/mL",
                    "Seed train: 50 mL -> 2 L wave bag -> 50 L single-use bioreactor -> 500 L -> 2000-5000 L stainless-steel bioreactor",
                    "Main production: fed-batch 12-15 days, 37 °C -> temp shift to 32-34 °C after day 4-5 for productivity boost, pH 7.0, DO 30 %, glucose + amino-acid feeds",
                    "End of production: antibody titer 3-10 g/L in supernatant, cell viability still > 70 %",
                ],
            },
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
            "literature_source": "Mata-Gómez et al., Microb. Cell Fact. 2014 — Microbial carotenoid production (doi:10.1201/b17587-13); Lopes et al., World J. Microbiol. Biotechnol. 2019 — Blakeslea review (doi:10.1201/b17587-13)",
            "production": {
                "de": [
                    "Stammvorbereitung: Blakeslea trispora (Pilz, Co-Kultur + - Stämme nötig) ODER engineered Yarrowia lipolytica (Hefe, Solo-Kultur)",
                    "Vorkultur: Komplexes Medium (Glucose + pflanzliche Öle als C-Quelle), 26-28 °C, 48-72 h",
                    "Hauptfermentation: Submerskultur in 50-200 m3 Bioreaktor, 26-28 °C (Blakeslea) bzw. 28-30 °C (Yarrowia), pH 6.0-6.5",
                    "Substrat-Feed: Glucose + Sojaöl als C-Quelle, Beta-Ionon-Stimulanz (nur Blakeslea) zur Carotinoid-Produktion",
                    "Fermentations-Dauer 96-168 h, Beta-Carotin-Akkumulation intrazellulär 1-5 g/L (Lipidkörperchen)",
                ],
                "en": [
                    "Strain preparation: Blakeslea trispora (fungus, requires + and - mating strains in co-culture) OR engineered Yarrowia lipolytica (yeast, single culture)",
                    "Pre-culture: complex medium (glucose + vegetable oils as C source), 26-28 °C, 48-72 h",
                    "Main fermentation: submerged culture in 50-200 m3 bioreactor, 26-28 °C (Blakeslea) or 28-30 °C (Yarrowia), pH 6.0-6.5",
                    "Substrate feed: glucose + soy oil as C source, beta-ionone stimulant (Blakeslea only) for carotenoid production",
                    "Fermentation duration 96-168 h, beta-carotene accumulation intracellular 1-5 g/L (lipid bodies)",
                ],
            },
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
            "literature_source": "Shah et al., Front. Plant Sci. 2016 — Astaxanthin from Haematococcus (doi:10.3389/fpls.2016.00531); Rodríguez-Sáiz et al., Appl. Microbiol. Biotechnol. 2010 — Phaffia review (doi:10.3389/fpls.2016.00531)",
            "production": {
                "de": [
                    "Stammvorbereitung: Haematococcus pluvialis (Grünalge, photoautotroph) ODER Phaffia rhodozyma (Hefe, heterotroph)",
                    "Vorkultur Haematococcus: Photobioreaktor mit Modifiziertem BBM-Medium, 22-25 °C, niedrige Lichtintensität (40 µmol Photonen/m2/s), grünes Wachstumsstadium",
                    "Vorkultur Phaffia: Glucose-Hefeextrakt-Medium, 20-22 °C, 48-72 h",
                    "Stress-Induktion (Haematococcus): Wechsel zu hoher Lichtintensität (300-500 µmol/m2/s), Stickstoff-Limitierung → Aplanosporen-Bildung + Astaxanthin-Akkumulation",
                    "Produktions-Dauer 10-20 Tage (Algen) bzw. 96-168 h (Phaffia), Astaxanthin-Gehalt 1.5-5 % der Trockenmasse",
                ],
                "en": [
                    "Strain preparation: Haematococcus pluvialis (green alga, photoautotrophic) OR Phaffia rhodozyma (yeast, heterotrophic)",
                    "Pre-culture Haematococcus: photobioreactor with modified BBM medium, 22-25 °C, low light intensity (40 umol photons/m2/s), green growth stage",
                    "Pre-culture Phaffia: glucose-yeast extract medium, 20-22 °C, 48-72 h",
                    "Stress induction (Haematococcus): switch to high light intensity (300-500 umol/m2/s), nitrogen limitation -> aplanospore formation + astaxanthin accumulation",
                    "Production duration 10-20 days (algae) or 96-168 h (Phaffia), astaxanthin content 1.5-5 % of dry biomass",
                ],
            },
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
            "production": {
                "de": [
                    "Rohmaterial: rohe Kaffeebohnen (Robusta bevorzugt, höherer Koffein-Gehalt: 2-4 % w/w) oder Schwarztee-Reste",
                    "Vorbehandlung: Befeuchtung der Bohnen mit Wasserdampf (öffnet Zellstruktur) oder direkte Verarbeitung der Tee-Aufgüsse",
                    "Koffein-Gehalt im Rohmaterial: Kaffee 1.0-2.5 %, Robusta-Kaffee 2.0-4.0 %, Tee 2.0-4.0 % (Trockenmasse)",
                    "Hauptquelle in industrieller Praxis: Entkoffeinierungs-Anlagen (Mass-Balance) als Nebenprodukt",
                ],
                "en": [
                    "Raw material: raw coffee beans (Robusta preferred, higher caffeine content: 2-4 % w/w) or black tea waste",
                    "Pre-treatment: steam moistening of beans (opens cell structure) or direct processing of tea infusions",
                    "Caffeine content in raw material: coffee 1.0-2.5 %, Robusta coffee 2.0-4.0 %, tea 2.0-4.0 % (dry mass)",
                    "Main source in industrial practice: decaffeination plants (mass balance) as by-product",
                ],
            },
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
            "literature_source": "Kutchan, The Alkaloids 1998 — Morphine biosynthesis and extraction; Beaudoin & Facchini, Planta 2014 — poppy alkaloid review (doi:10.1007/s00425-014-2056-8)",
            "production": {
                "de": [
                    "Rohmaterial: Papaver somniferum-Latex (eingetrocknet zu Opium) ODER konzentriertes Mohnstroh (CPS) aus zugelassenem Anbau (BtMG-Lizenz!)",
                    "Anbau-Logistik: Vertragsanbau unter Aufsicht der zuständigen Bundesopiumstelle (Deutschland) oder DEA/INCB (international)",
                    "Vorbehandlung Mohnstroh: Mahlung auf 2-5 mm; Latex: Trocknung und Mahlung",
                    "Alkaloid-Gehalt im Rohmaterial: Latex 8-14 % Morphin (Trockenmasse); Mohnstroh 0.3-1.5 %",
                    "Zusatzliche regulatorische Anforderungen: GMP-Tracking + DEA/BtMG-Mengenbilanzierung für jeden Schritt",
                ],
                "en": [
                    "Raw material: Papaver somniferum latex (dried to opium) OR concentrated poppy straw (CPS) from licensed cultivation (controlled-substance licence required!)",
                    "Cultivation logistics: contract farming under supervision of national opium agency (BfArM in Germany) or DEA/INCB (international)",
                    "Pre-treatment poppy straw: milling to 2-5 mm; latex: drying and milling",
                    "Alkaloid content in raw material: latex 8-14 % morphine (dry mass); poppy straw 0.3-1.5 %",
                    "Additional regulatory requirements: GMP tracking + DEA/controlled-substance mass balance for every step",
                ],
            },
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


# Merge in the extended hints for the remaining 57 molecules of the DB.
# Kept in a separate module so this main file stays readable.
try:
    from molecule_hints_extension import EXTENSION_HINTS as _EXT_HINTS
    MOLECULE_HINTS.update(_EXT_HINTS)
except Exception:
    pass


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
        "literature_source": "Walsh G., Pharmaceutical Biotechnology, 2nd ed. 2018, Wiley (ISBN 978-1-119-11518-7) — antibody/protein downstream processing; Liu et al., mAbs 2010 — Recovery and purification of monoclonal antibodies (doi:10.4161/mabs.2.5.12645)",
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
      - production_steps (list of concrete upstream / synthesis operations)
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
        "production_steps": [],
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

        # Upstream / production steps (Feedback-Runde-3: war vorher nur als
        # Freitext im "production_route" — jetzt zusätzlich als strukturierte
        # Schritte-Liste für die PDF-Sektion + Process-Flow-Diagramm.
        production = hint.get("production") or {}
        if production:
            out["production_steps"] = list(production.get(lang) or production.get("de") or [])

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
