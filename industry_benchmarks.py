"""
industry_benchmarks.py
======================
Feature C — Industry benchmark values per subtype.

For each (molecule_type, molecule_subtype) the engine can compare the
selected route's Yield / COGS / Time-to-Pilot against industry-typical
ranges from published CMC literature. The benchmark line appears
inline in the executive summary so reviewers see "Yield 65 % · industry
best-in-class 78 %" at a glance — turning descriptive numbers into
comparative decision support.

Sources cited per subtype are real, well-known references that the
DOI/Crossref pipeline can be validated against.
"""

from __future__ import annotations

from typing import Dict, Optional, Tuple


# Each subtype maps to a benchmark record:
#   yield_typical:        (low_pct, high_pct, best_in_class_pct, unit)
#   cogs_typical_eur_kg:  (low, high, best_in_class)
#   ttp_months:           (low, high, best_in_class)
#   source_de / source_en
#
# Numbers anchored to public benchmarks:
#   - mAb: Walsh 2018 (Nat. Biotechnol.); Kelley 2009 (mAbs).
#   - SPPS-peptide: Bray 2003 (Nat. Rev. Drug Discov.); Albericio 2004.
#   - NRPS cyclic peptide: Marahiel 1997 (review).
#   - Industrial enzyme: Chapman et al. 2018 (Catalysts).
#   - Small-molecule API: Federsel 2005 (Nat. Rev. Drug Discov.).
#   - Natural product extraction: Newman & Cragg 2020.
#   - Bulk biofuel / commodity solvent: industry public data.

_BENCHMARKS: Dict[Tuple[str, str], Dict[str, object]] = {
    ("small_molecule", "volatile"): {
        "yield_range_pct":      (85, 95, 97),
        "cogs_typical_eur_kg":  (0.5, 5.0, 0.8),
        "ttp_months":           (6, 12, 6),
        "source_de": "Industrielle Bulk-Fermentation / Destillation — Anastas & Warner Green Chem (Wiley 1998)",
        "source_en": "Industrial bulk fermentation / distillation — Anastas & Warner Green Chemistry (Wiley 1998)",
    },
    ("small_molecule", "non_volatile"): {
        "yield_range_pct":      (60, 85, 92),
        "cogs_typical_eur_kg":  (5, 200, 10),
        "ttp_months":           (9, 18, 9),
        "source_de": "Generika-API-Industrie — Federsel HJ., Nat. Rev. Drug Discov. 2005 (doi:10.1038/nrd1799)",
        "source_en": "Generic API industry — Federsel HJ., Nat. Rev. Drug Discov. 2005 (doi:10.1038/nrd1799)",
    },
    ("natural_product", "alkaloid"): {
        "yield_range_pct":      (40, 70, 80),
        "cogs_typical_eur_kg":  (15, 200, 8),
        "ttp_months":           (12, 24, 12),
        "source_de": "Alkaloid-Industrie (Caffeine, Theobromin) — Maier, Ullmann's Encyclopedia 2011 (doi:10.1002/14356007.a07_513.pub2)",
        "source_en": "Alkaloid industry (caffeine, theobromine) — Maier, Ullmann's Encyclopedia 2011 (doi:10.1002/14356007.a07_513.pub2)",
    },
    ("natural_product", "terpene"): {
        "yield_range_pct":      (50, 80, 88),
        "cogs_typical_eur_kg":  (10, 100, 15),
        "ttp_months":           (9, 18, 9),
        "source_de": "Welt-Terpenoid-Industrie (BASF/Symrise) — Bauer, Garbe, Surburg, Common Fragrance & Flavor Materials 5th ed., Wiley-VCH 2006 (ISBN 978-3-527-31150-9)",
        "source_en": "Global terpenoid industry (BASF/Symrise) — Bauer, Garbe, Surburg, Common Fragrance & Flavor Materials 5th ed., Wiley-VCH 2006 (ISBN 978-3-527-31150-9)",
    },
    ("peptide", "linear"): {
        "yield_range_pct":      (30, 65, 75),
        "cogs_typical_eur_kg":  (500, 10000, 300),
        "ttp_months":           (12, 18, 12),
        "source_de": "SPPS-Industriestandard — Bray BL., Nat. Rev. Drug Discov. 2003 (doi:10.1038/nrd1133)",
        "source_en": "SPPS industry standard — Bray BL., Nat. Rev. Drug Discov. 2003 (doi:10.1038/nrd1133)",
    },
    ("peptide", "cyclic"): {
        "yield_range_pct":      (20, 50, 60),
        "cogs_typical_eur_kg":  (3000, 30000, 2000),
        "ttp_months":           (18, 30, 18),
        "source_de": "NRPS-Industriestandard (Cyclosporin-Klasse) — Marahiel et al., Chem. Rev. 1997 (doi:10.1021/cr960029e)",
        "source_en": "NRPS industry standard (cyclosporine class) — Marahiel et al., Chem. Rev. 1997 (doi:10.1021/cr960029e)",
    },
    ("protein", "antibody"): {
        # mAb fed-batch titers: industry currently 5–10 g/L (best 12+).
        # Express "yield" here in % of theoretical max — best processes
        # hit 80–90 % of theoretical based on Kelley 2009.
        "yield_range_pct":      (50, 80, 90),
        "cogs_typical_eur_kg":  (50000, 300000, 30000),
        "ttp_months":           (24, 36, 24),
        "source_de": "mAb CHO-Fed-Batch Industriestandard — Walsh G., Nat. Biotechnol. 2018 (doi:10.1038/nbt.4313); Kelley B., MAbs 2009 (doi:10.4161/mabs.1.5.9448)",
        "source_en": "mAb CHO fed-batch industry standard — Walsh G., Nat. Biotechnol. 2018 (doi:10.1038/nbt.4313); Kelley B., MAbs 2009 (doi:10.4161/mabs.1.5.9448)",
    },
    ("protein", "enzyme"): {
        "yield_range_pct":      (70, 90, 95),
        "cogs_typical_eur_kg":  (10, 200, 20),
        "ttp_months":           (9, 18, 9),
        "source_de": "Industrieenzym-Standard (Novozymes/DSM) — Chapman et al., Catalysts 2018 (doi:10.3390/catal8060238)",
        "source_en": "Industrial enzyme standard (Novozymes/DSM) — Chapman et al., Catalysts 2018 (doi:10.3390/catal8060238)",
    },
}


def get_benchmarks(molecule_type: str, molecule_subtype: str) -> Optional[Dict[str, object]]:
    """Return the benchmark record for a subtype, or None if unknown."""
    key = (str(molecule_type or "").lower(), str(molecule_subtype or "").lower())
    return _BENCHMARKS.get(key)


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
