"""
industry_benchmarks.py
======================
Feature C — Industry benchmark values per subtype and method.

For each (molecule_type, molecule_subtype, method) the engine can compare
the selected route's Yield / COGS / Time-to-Pilot against industry-typical
ranges from published CMC literature. The benchmark line appears inline in
the executive summary so reviewers see "Yield 65 % · industry best-in-class
78 %" at a glance — turning descriptive numbers into comparative decision
support.

Lookup priority used by ``get_benchmarks``:

    1. Per-molecule COGS override (pulled lazily from
       ``cogs_estimator._MOLECULE_COGS_OVERRIDES`` so we stay in sync with
       the COGS estimator — see P1 blindspot fix).
    2. 3D cluster ``(type, subtype, method)``.
    3. 2D cluster ``(type, subtype)`` (legacy fallback so old callers that
       don't pass ``method`` still work).

Sources cited per cluster are real, well-known references that the
DOI/Crossref pipeline can be validated against.
"""

from __future__ import annotations

from typing import Dict, Optional, Tuple


# ---------------------------------------------------------------------------
# Method-normalization (mirrors cogs_estimator._normalize_method intentionally)
# ---------------------------------------------------------------------------
def _normalize_method(method: Optional[str]) -> str:
    m = str(method or "").lower().strip()
    if m in ("chem", "chemical"):
        return "chemical"
    if m in ("bio", "biotech", "biotechnological"):
        return "biotechnological"
    if m in ("extract", "extraction"):
        return "extraction"
    if m in ("enz", "enzymatic", "enzyme"):
        return "enzymatic"
    return m


# ---------------------------------------------------------------------------
# 3D cluster (type, subtype, method) — primary table
# ---------------------------------------------------------------------------
# Each record:
#   yield_range_pct:        (low, high, best_in_class)  in %
#   cogs_typical_eur_kg:    (low, high, best_in_class)  in €/kg
#   ttp_months:             (low, high, best_in_class)  in months
#   source_de / source_en
#
# Sources:
#   - Federsel HJ., Nat. Rev. Drug Discov. 2005 (doi:10.1038/nrd1799)
#   - Walsh G., Nat. Biotechnol. 2018 (doi:10.1038/nbt.4313)
#   - Kelley B., MAbs 2009 (doi:10.4161/mabs.1.5.9448)
#   - Bray BL., Nat. Rev. Drug Discov. 2003 (doi:10.1038/nrd1133)
#   - Marahiel MA., Chem. Rev. 1997 (doi:10.1021/cr960029e)
#   - Chapman et al., Catalysts 2018 (doi:10.3390/catal8060238)
#   - Maier, Ullmann's Encyclopedia 2011 (doi:10.1002/14356007.a07_513.pub2)
#   - Bauer, Garbe, Surburg, Common Fragrance & Flavor Materials 5th ed.,
#     Wiley-VCH 2006 (ISBN 978-3-527-31150-9)
#   - Newman & Cragg, J. Nat. Prod. 2020 (doi:10.1021/acs.jnatprod.9b01285)
#   - Anastas & Warner, Green Chemistry, Wiley 1998 (ISBN 978-0-19-850698-0)

_BENCHMARKS_BY_METHOD: Dict[Tuple[str, str, str], Dict[str, object]] = {

    # ----- small_molecule / volatile -----
    ("small_molecule", "volatile", "chemical"): {
        "yield_range_pct":     (80, 95, 97),
        "cogs_typical_eur_kg": (0.3, 3.0, 0.5),
        "ttp_months":          (6, 12, 6),
        "source_de": "Welt-Commodity-Solvents (Cumol-/Oxo-/Reforming-Prozesse) — Anastas & Warner Green Chem (Wiley 1998)",
        "source_en": "World commodity solvents (cumene/oxo/reforming) — Anastas & Warner Green Chemistry (Wiley 1998)",
    },
    ("small_molecule", "volatile", "biotechnological"): {
        "yield_range_pct":     (70, 88, 93),
        "cogs_typical_eur_kg": (0.5, 4.0, 0.8),
        "ttp_months":          (12, 24, 12),
        "source_de": "Bioethanol / Fermentationsalkohole — Industriedaten, RFA/Renewable Fuels Association 2022",
        "source_en": "Bioethanol / fermentative alcohols — industry data, RFA/Renewable Fuels Association 2022",
    },

    # ----- small_molecule / non_volatile -----
    ("small_molecule", "non_volatile", "chemical"): {
        "yield_range_pct":     (60, 85, 92),
        "cogs_typical_eur_kg": (3, 30, 5),
        "ttp_months":          (9, 18, 9),
        "source_de": "Generika-API-Industrie — Federsel HJ., Nat. Rev. Drug Discov. 2005 (doi:10.1038/nrd1799)",
        "source_en": "Generic API industry — Federsel HJ., Nat. Rev. Drug Discov. 2005 (doi:10.1038/nrd1799)",
    },
    ("small_molecule", "non_volatile", "biotechnological"): {
        "yield_range_pct":     (50, 80, 88),
        "cogs_typical_eur_kg": (1, 10, 1.5),
        "ttp_months":          (12, 24, 12),
        "source_de": "Bulk-Fermentation organischer Säuren (Citrat/Lactat/Succinat) — Industriedaten",
        "source_en": "Bulk fermentation of organic acids (citric/lactic/succinic) — industry data",
    },
    ("small_molecule", "non_volatile", "enzymatic"): {
        "yield_range_pct":     (75, 92, 96),
        "cogs_typical_eur_kg": (0.3, 5, 0.5),
        "ttp_months":          (6, 12, 6),
        "source_de": "Enzymatische Bulk-Hydrolyse (Stärke/Cellulose) — Chapman et al., Catalysts 2018 (doi:10.3390/catal8060238)",
        "source_en": "Enzymatic bulk hydrolysis (starch/cellulose) — Chapman et al., Catalysts 2018 (doi:10.3390/catal8060238)",
    },

    # ----- natural_product / alkaloid -----
    ("natural_product", "alkaloid", "extraction"): {
        "yield_range_pct":     (40, 70, 80),
        "cogs_typical_eur_kg": (40, 800, 30),
        "ttp_months":          (12, 24, 12),
        "source_de": "Alkaloid-Extraktion (Cinchona/Papaver) — Maier, Ullmann's Encyclopedia 2011 (doi:10.1002/14356007.a07_513.pub2)",
        "source_en": "Alkaloid extraction (cinchona/papaver) — Maier, Ullmann's Encyclopedia 2011 (doi:10.1002/14356007.a07_513.pub2)",
    },
    ("natural_product", "alkaloid", "chemical"): {
        "yield_range_pct":     (45, 75, 85),
        "cogs_typical_eur_kg": (8, 50, 6),
        "ttp_months":          (9, 18, 9),
        "source_de": "Synthese-Alkaloide Bulk (Coffein-Traube-Synthese) — Maier, Ullmann's Encyclopedia 2011 (doi:10.1002/14356007.a07_513.pub2)",
        "source_en": "Synthetic alkaloid bulk (caffeine Traube synthesis) — Maier, Ullmann's Encyclopedia 2011 (doi:10.1002/14356007.a07_513.pub2)",
    },
    ("natural_product", "alkaloid", "biotechnological"): {
        "yield_range_pct":     (15, 45, 60),
        "cogs_typical_eur_kg": (100, 1000, 50),
        "ttp_months":          (18, 36, 18),
        "source_de": "Engineered Mikroben für Alkaloide (Antheia/Scopolamin in Hefe) — Galanie et al., Science 2015 (doi:10.1126/science.aac9373)",
        "source_en": "Engineered microbes for alkaloids (Antheia/scopolamine in yeast) — Galanie et al., Science 2015 (doi:10.1126/science.aac9373)",
    },

    # ----- natural_product / terpene -----
    ("natural_product", "terpene", "extraction"): {
        "yield_range_pct":     (40, 75, 85),
        "cogs_typical_eur_kg": (15, 500, 12),
        "ttp_months":          (9, 18, 9),
        "source_de": "Ätherische Öle / Pflanzen-Extraktion — Bauer, Garbe, Surburg, Common Fragrance & Flavor Materials 5th ed., Wiley-VCH 2006 (ISBN 978-3-527-31150-9)",
        "source_en": "Essential oils / plant extraction — Bauer, Garbe, Surburg, Common Fragrance & Flavor Materials 5th ed., Wiley-VCH 2006 (ISBN 978-3-527-31150-9)",
    },
    ("natural_product", "terpene", "chemical"): {
        "yield_range_pct":     (60, 85, 92),
        "cogs_typical_eur_kg": (5, 80, 8),
        "ttp_months":          (6, 12, 6),
        "source_de": "Semi-synthetische Terpenoide (BASF/Symrise Citral-Route) — Bauer, Garbe, Surburg, Wiley-VCH 2006 (ISBN 978-3-527-31150-9)",
        "source_en": "Semi-synthetic terpenoids (BASF/Symrise citral route) — Bauer, Garbe, Surburg, Wiley-VCH 2006 (ISBN 978-3-527-31150-9)",
    },
    ("natural_product", "terpene", "biotechnological"): {
        "yield_range_pct":     (30, 60, 75),
        "cogs_typical_eur_kg": (50, 3000, 40),
        "ttp_months":          (18, 30, 18),
        "source_de": "Engineered Hefe / Mikroalgen (Amyris-Artemisinin, Haematococcus-Astaxanthin) — Paddon & Keasling, Nat. Rev. Microbiol. 2014 (doi:10.1038/nrmicro3240)",
        "source_en": "Engineered yeast / microalgae (Amyris artemisinin, Haematococcus astaxanthin) — Paddon & Keasling, Nat. Rev. Microbiol. 2014 (doi:10.1038/nrmicro3240)",
    },

    # ----- peptide / linear -----
    ("peptide", "linear", "chemical"): {
        "yield_range_pct":     (30, 60, 75),
        "cogs_typical_eur_kg": (500, 50000, 300),
        "ttp_months":          (12, 24, 12),
        "source_de": "SPPS-Industriestandard (Bachem/PolyPeptide) — Bray BL., Nat. Rev. Drug Discov. 2003 (doi:10.1038/nrd1133)",
        "source_en": "SPPS industry standard (Bachem/PolyPeptide) — Bray BL., Nat. Rev. Drug Discov. 2003 (doi:10.1038/nrd1133)",
    },
    ("peptide", "linear", "biotechnological"): {
        "yield_range_pct":     (40, 70, 80),
        "cogs_typical_eur_kg": (100, 5000, 80),
        "ttp_months":          (18, 30, 18),
        "source_de": "Rekombinante Peptid-Expression (Insulin in E. coli/Hefe) — Walsh G., Nat. Biotechnol. 2018 (doi:10.1038/nbt.4313)",
        "source_en": "Recombinant peptide expression (insulin in E. coli/yeast) — Walsh G., Nat. Biotechnol. 2018 (doi:10.1038/nbt.4313)",
    },

    # ----- peptide / cyclic -----
    ("peptide", "cyclic", "chemical"): {
        "yield_range_pct":     (20, 45, 55),
        "cogs_typical_eur_kg": (3000, 100000, 2000),
        "ttp_months":          (18, 30, 18),
        "source_de": "SPPS + Makrozyklisierung — Bray BL., Nat. Rev. Drug Discov. 2003 (doi:10.1038/nrd1133)",
        "source_en": "SPPS + macrocyclization — Bray BL., Nat. Rev. Drug Discov. 2003 (doi:10.1038/nrd1133)",
    },
    ("peptide", "cyclic", "biotechnological"): {
        "yield_range_pct":     (20, 50, 60),
        "cogs_typical_eur_kg": (1000, 50000, 800),
        "ttp_months":          (18, 30, 18),
        "source_de": "NRPS-Fermentation (Cyclosporin/Vancomycin-Klasse) — Marahiel et al., Chem. Rev. 1997 (doi:10.1021/cr960029e)",
        "source_en": "NRPS fermentation (cyclosporine/vancomycin class) — Marahiel et al., Chem. Rev. 1997 (doi:10.1021/cr960029e)",
    },

    # ----- protein / antibody -----
    ("protein", "antibody", "biotechnological"): {
        # mAb fed-batch titers: industry currently 5–10 g/L (best 12+).
        # Yield here in % of theoretical max — best processes hit 80–90 %.
        "yield_range_pct":     (50, 80, 90),
        "cogs_typical_eur_kg": (50000, 300000, 30000),
        "ttp_months":          (24, 36, 24),
        "source_de": "mAb CHO-Fed-Batch Industriestandard — Walsh G., Nat. Biotechnol. 2018 (doi:10.1038/nbt.4313); Kelley B., MAbs 2009 (doi:10.4161/mabs.1.5.9448)",
        "source_en": "mAb CHO fed-batch industry standard — Walsh G., Nat. Biotechnol. 2018 (doi:10.1038/nbt.4313); Kelley B., MAbs 2009 (doi:10.4161/mabs.1.5.9448)",
    },

    # ----- protein / enzyme -----
    ("protein", "enzyme", "biotechnological"): {
        # Industrial enzymes; therapeutic recombinant enzymes get per-molecule overrides.
        "yield_range_pct":     (70, 90, 95),
        "cogs_typical_eur_kg": (10, 200, 10),
        "ttp_months":          (9, 18, 9),
        "source_de": "Industrieenzym-Standard (Novozymes/DSM) — Chapman et al., Catalysts 2018 (doi:10.3390/catal8060238)",
        "source_en": "Industrial enzyme standard (Novozymes/DSM) — Chapman et al., Catalysts 2018 (doi:10.3390/catal8060238)",
    },
}


# ---------------------------------------------------------------------------
# 2D cluster (type, subtype) — legacy fallback for callers that don't pass
# ``method``. Numbers chosen as a "method-agnostic" middle ground.
# ---------------------------------------------------------------------------
_BENCHMARKS: Dict[Tuple[str, str], Dict[str, object]] = {
    ("small_molecule", "volatile"): {
        "yield_range_pct":      (75, 92, 95),
        "cogs_typical_eur_kg":  (0.5, 4.0, 0.8),
        "ttp_months":           (6, 18, 6),
        "source_de": "Welt-Commodity-Solvents — Anastas & Warner Green Chem (Wiley 1998)",
        "source_en": "World commodity solvents — Anastas & Warner Green Chemistry (Wiley 1998)",
    },
    ("small_molecule", "non_volatile"): {
        "yield_range_pct":      (60, 85, 92),
        "cogs_typical_eur_kg":  (3, 50, 5),
        "ttp_months":           (9, 18, 9),
        "source_de": "Generika-API-Industrie — Federsel HJ., Nat. Rev. Drug Discov. 2005 (doi:10.1038/nrd1799)",
        "source_en": "Generic API industry — Federsel HJ., Nat. Rev. Drug Discov. 2005 (doi:10.1038/nrd1799)",
    },
    ("natural_product", "alkaloid"): {
        "yield_range_pct":      (40, 70, 80),
        "cogs_typical_eur_kg":  (15, 500, 15),
        "ttp_months":           (12, 24, 12),
        "source_de": "Alkaloid-Industrie — Maier, Ullmann's Encyclopedia 2011 (doi:10.1002/14356007.a07_513.pub2)",
        "source_en": "Alkaloid industry — Maier, Ullmann's Encyclopedia 2011 (doi:10.1002/14356007.a07_513.pub2)",
    },
    ("natural_product", "terpene"): {
        "yield_range_pct":      (50, 80, 88),
        "cogs_typical_eur_kg":  (10, 300, 15),
        "ttp_months":           (9, 18, 9),
        "source_de": "Welt-Terpenoid-Industrie — Bauer, Garbe, Surburg, Common Fragrance & Flavor Materials 5th ed., Wiley-VCH 2006 (ISBN 978-3-527-31150-9)",
        "source_en": "Global terpenoid industry — Bauer, Garbe, Surburg, Common Fragrance & Flavor Materials 5th ed., Wiley-VCH 2006 (ISBN 978-3-527-31150-9)",
    },
    ("peptide", "linear"): {
        "yield_range_pct":      (30, 65, 75),
        "cogs_typical_eur_kg":  (200, 30000, 200),
        "ttp_months":           (12, 24, 12),
        "source_de": "SPPS/Rekombinante Peptide — Bray BL., Nat. Rev. Drug Discov. 2003 (doi:10.1038/nrd1133)",
        "source_en": "SPPS/recombinant peptides — Bray BL., Nat. Rev. Drug Discov. 2003 (doi:10.1038/nrd1133)",
    },
    ("peptide", "cyclic"): {
        "yield_range_pct":      (20, 50, 60),
        "cogs_typical_eur_kg":  (2000, 80000, 1500),
        "ttp_months":           (18, 30, 18),
        "source_de": "NRPS/SPPS-Industriestandard — Marahiel et al., Chem. Rev. 1997 (doi:10.1021/cr960029e)",
        "source_en": "NRPS/SPPS industry standard — Marahiel et al., Chem. Rev. 1997 (doi:10.1021/cr960029e)",
    },
    ("protein", "antibody"): {
        "yield_range_pct":      (50, 80, 90),
        "cogs_typical_eur_kg":  (50000, 300000, 30000),
        "ttp_months":           (24, 36, 24),
        "source_de": "mAb CHO-Fed-Batch — Walsh G., Nat. Biotechnol. 2018 (doi:10.1038/nbt.4313); Kelley B., MAbs 2009 (doi:10.4161/mabs.1.5.9448)",
        "source_en": "mAb CHO fed-batch — Walsh G., Nat. Biotechnol. 2018 (doi:10.1038/nbt.4313); Kelley B., MAbs 2009 (doi:10.4161/mabs.1.5.9448)",
    },
    ("protein", "enzyme"): {
        "yield_range_pct":      (70, 90, 95),
        "cogs_typical_eur_kg":  (10, 200, 20),
        "ttp_months":           (9, 18, 9),
        "source_de": "Industrieenzym-Standard — Chapman et al., Catalysts 2018 (doi:10.3390/catal8060238)",
        "source_en": "Industrial enzyme standard — Chapman et al., Catalysts 2018 (doi:10.3390/catal8060238)",
    },
}


# ---------------------------------------------------------------------------
# Per-molecule yield/ttp overrides where the cluster default is meaningfully
# off for a specific molecule. COGS overrides are pulled lazily from
# cogs_estimator._MOLECULE_COGS_OVERRIDES (single source of truth — see P1).
# ---------------------------------------------------------------------------
_MOLECULE_BENCHMARK_OVERRIDES: Dict[str, Dict[str, object]] = {
    # Mature commodity → very high yield, very fast time-to-pilot.
    "vanillin": {
        "yield_range_pct": (82, 92, 95),
        "ttp_months":      (3, 9, 3),
        "source_de": "Borregaard/Solvay Synthese-Prozess (Guaiacol/Glyoxylsäure-Route) — Hocking 1997, J. Chem. Educ. (doi:10.1021/ed074p1055)",
        "source_en": "Borregaard/Solvay synthesis process (guaiacol/glyoxylic acid route) — Hocking 1997, J. Chem. Educ. (doi:10.1021/ed074p1055)",
    },
    "aspirin": {
        "yield_range_pct": (85, 95, 97),
        "ttp_months":      (3, 6, 3),
    },
    "paracetamol": {
        "yield_range_pct": (80, 92, 95),
        "ttp_months":      (3, 6, 3),
    },
    "ibuprofen": {
        "yield_range_pct": (80, 90, 95),
        "ttp_months":      (3, 6, 3),
    },
    # Top-onkologische mAbs: Scale-up zur Phase III/Markt dauert lang.
    "trastuzumab":   {"yield_range_pct": (55, 78, 88), "ttp_months": (24, 36, 24)},
    "adalimumab":    {"yield_range_pct": (60, 80, 90), "ttp_months": (24, 36, 24)},
    "pembrolizumab": {"yield_range_pct": (50, 75, 88), "ttp_months": (30, 42, 30)},
    "nivolumab":     {"yield_range_pct": (50, 75, 88), "ttp_months": (30, 42, 30)},
    # Insulin: industriereife rekombinante Expression seit 1980er.
    "insulin":       {"yield_range_pct": (60, 80, 88), "ttp_months": (12, 18, 12)},
    # GLP-1-Analoga (SPPS + Lipidierung): Ausbeute durch lange Sequenz limitiert.
    "liraglutide":   {"yield_range_pct": (25, 45, 55), "ttp_months": (15, 24, 15)},
    "exenatide":     {"yield_range_pct": (25, 45, 55), "ttp_months": (15, 24, 15)},
    # Cyclische Peptide via NRPS: Sehr lange Etablierung des Produktionsstamms.
    "cyclosporine":  {"yield_range_pct": (25, 50, 60), "ttp_months": (18, 30, 18)},
    "vancomycin":    {"yield_range_pct": (30, 55, 65), "ttp_months": (15, 24, 15)},
    "daptomycin":    {"yield_range_pct": (20, 45, 55), "ttp_months": (18, 30, 18)},
    # Industrieenzyme: Hochskaliert, sehr schnell.
    "amylase":       {"yield_range_pct": (75, 92, 96), "ttp_months": (6, 12, 6)},
    "lipase":        {"yield_range_pct": (72, 90, 94), "ttp_months": (6, 12, 6)},
    "protease":      {"yield_range_pct": (75, 92, 96), "ttp_months": (6, 12, 6)},
    "cellulase":     {"yield_range_pct": (70, 88, 93), "ttp_months": (6, 12, 6)},
    # Rekombinante therapeutische Enzyme (CHO): Scale-up vergleichbar mAb.
    "erythropoietin": {"yield_range_pct": (50, 75, 85), "ttp_months": (24, 36, 24)},
    "filgrastim":    {"yield_range_pct": (55, 78, 88), "ttp_months": (18, 30, 18)},
    "somatropin":    {"yield_range_pct": (55, 78, 88), "ttp_months": (18, 30, 18)},
    "asparaginase":  {"yield_range_pct": (60, 82, 90), "ttp_months": (15, 24, 15)},
    # Naturstoffe mit etablierter Extraktion.
    "morphine":      {"yield_range_pct": (50, 75, 85), "ttp_months": (9, 18, 9)},
    "codeine":       {"yield_range_pct": (50, 75, 85), "ttp_months": (9, 18, 9)},
    "quinine":       {"yield_range_pct": (45, 70, 80), "ttp_months": (12, 18, 12)},
    "caffeine":      {"yield_range_pct": (75, 90, 95), "ttp_months": (6, 12, 6)},
    # Engineered-Microbial Terpenoide (lange Stamm-Entwicklung).
    "artemisinin":   {"yield_range_pct": (30, 55, 70), "ttp_months": (18, 30, 18)},
    "astaxanthin":   {"yield_range_pct": (25, 50, 65), "ttp_months": (18, 30, 18)},
    # Bulk-Lösungsmittel-Commodities.
    "ethanol":   {"yield_range_pct": (85, 95, 97), "ttp_months": (3, 9, 3)},
    "methanol":  {"yield_range_pct": (90, 97, 98), "ttp_months": (3, 6, 3)},
    "acetone":   {"yield_range_pct": (88, 96, 98), "ttp_months": (3, 6, 3)},
    "glucose":   {"yield_range_pct": (90, 97, 98), "ttp_months": (3, 6, 3)},
}


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------
def get_benchmarks(molecule_type: str,
                   molecule_subtype: str,
                   method: Optional[str] = None,
                   molecule_name: Optional[str] = None) -> Optional[Dict[str, object]]:
    """Return the benchmark record for a (type, subtype, method) combination,
    with optional per-molecule overrides applied on top.

    Lookup order:
        1. 3D cluster (type, subtype, method) if method given
        2. 2D cluster (type, subtype) — legacy fallback
        3. Per-molecule override (cogs from cogs_estimator, yield/ttp from
           this module's _MOLECULE_BENCHMARK_OVERRIDES) — applied on top of
           the cluster record.

    Returns None if no benchmark is known at any level.
    """
    t = str(molecule_type or "").lower().strip()
    s = str(molecule_subtype or "").lower().strip()
    m = _normalize_method(method)
    name = str(molecule_name or "").lower().strip()

    # 1. 3D cluster
    rec: Optional[Dict[str, object]] = None
    if m:
        rec = _BENCHMARKS_BY_METHOD.get((t, s, m))
    # 2. 2D fallback
    if rec is None:
        rec = _BENCHMARKS.get((t, s))
    if rec is None:
        return None

    # Copy so we don't mutate the module-level dict when applying overrides.
    rec = dict(rec)

    # 3. Per-molecule yield/ttp override.
    if name and name in _MOLECULE_BENCHMARK_OVERRIDES:
        ov = _MOLECULE_BENCHMARK_OVERRIDES[name]
        for k in ("yield_range_pct", "ttp_months", "source_de", "source_en"):
            if k in ov:
                rec[k] = ov[k]

    # 4. Per-molecule COGS override — lazy import to avoid circular dependency
    #    and keep cogs_estimator._MOLECULE_COGS_OVERRIDES as single source of
    #    truth for COGS numbers.
    if name:
        try:
            from cogs_estimator import _MOLECULE_COGS_OVERRIDES  # type: ignore
            cogs_ov = _MOLECULE_COGS_OVERRIDES.get(name)
            if cogs_ov:
                lo = float(cogs_ov["base_low"])
                hi = float(cogs_ov["base_high"])
                # Best-in-class ≈ optimized industrial process near the low end.
                rec["cogs_typical_eur_kg"] = (lo, hi, round(lo * 0.85, 2))
        except Exception:
            # Defensive: never let benchmark lookup blow up the report.
            pass

    return rec


def fmt_yield_comparison(observed_pct_str: Optional[str],
                         bench: Dict[str, object],
                         lang: str) -> Optional[str]:
    """Build a "observed vs industry-best" comparison string for the yield
    line. observed_pct_str is the engine's "60–75 %" style string."""
    if not bench:
        return None
    bic = bench.get("yield_range_pct")
    if not bic or len(bic) < 3:
        return None
    low, high, best = bic
    if lang == "en":
        return f"industry typical {low}–{high} %, best-in-class {best} %"
    return f"Industrie typ. {low}–{high} %, best-in-class {best} %"


def fmt_cogs_comparison(low_obs, high_obs, bench: Dict[str, object], lang: str) -> Optional[str]:
    """Build a "observed vs industry-typical" string for the COGS line."""
    if not bench:
        return None
    bic = bench.get("cogs_typical_eur_kg")
    if not bic or len(bic) < 3:
        return None
    low_b, high_b, best_b = bic

    def _fmt(x):
        if x >= 1000:
            s = f"{int(round(x)):,}".replace(",", ".") if lang != "en" else f"{int(round(x)):,}"
            return s
        if x >= 10:
            return f"{int(round(x))}"
        return f"{x:.1f}".replace(".", ",") if lang != "en" else f"{x:.1f}"

    if lang == "en":
        return f"industry typical {_fmt(low_b)}–{_fmt(high_b)} €/kg, best-in-class {_fmt(best_b)} €/kg"
    return f"Industrie typ. {_fmt(low_b)}–{_fmt(high_b)} €/kg, best-in-class {_fmt(best_b)} €/kg"


def fmt_ttp_comparison(bench: Dict[str, object], lang: str) -> Optional[str]:
    """Time-to-Pilot industry comparison."""
    if not bench:
        return None
    bic = bench.get("ttp_months")
    if not bic or len(bic) < 3:
        return None
    low, high, best = bic
    if lang == "en":
        return f"industry typical {low}–{high} months, best-in-class {best} months"
    return f"Industrie typ. {low}–{high} Monate, best-in-class {best} Monate"
