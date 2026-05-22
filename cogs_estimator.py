"""
cogs_estimator.py
=================
Domain-grounded COGS (Cost of Goods Sold) range estimator in €/kg.

Why: Pharma/biotech reviewers expect a concrete economic anchor next to
the qualitative cost rating. "Cost: HIGH" alone is unactionable;
"€4.500–6.200 /kg with 60 % from Protein-A capture" is.

How: Rule-based per molecule_type × molecule_subtype × method, anchored
to published price benchmarks for each class. The result is always a
RANGE (low/high) with a confidence label, because point estimates would
falsely imply precision the model doesn't have.

Public surface:
    estimate_cogs(process_input) -> dict with:
        low_eur_per_kg, high_eur_per_kg : float | None
        breakdown                       : dict[str, float]   (fractions, sum=1)
        confidence                      : "low" | "medium" | "high"
        notes_de, notes_en              : str
        anchor                          : short label of the reference class
        method_visible_drivers          : list[str]          (key drivers in human form)

Design constraints (from project memory):
- Regelbasiert, kein LLM — every output is traceable to an explicit rule.
- Honest representation: if we lack a benchmark for a combination,
  return None for the range and surface confidence=low.
- Numbers are anchored to publicly cited price benchmarks (typical CDMO
  pricing, peer-reviewed CMC literature). Surfaced in PDF references.
"""

from typing import Dict, Any, Optional, Tuple, List


# ----------------------------------------------------------------------
# Benchmark anchors per (type, subtype, method).
#
# Each anchor encodes:
#   base_low / base_high  : €/kg at "default" inputs (pharma-grade,
#                            mid-scale, no unusual factors)
#   breakdown             : typical fractional split into
#                            rm (raw materials), dsp (downstream
#                            processing), opex (operations), capex (asset
#                            amortisation)
#   anchor_label_de/en    : short reference describing the comparator
#   source                : public reference for traceability
#
# Sources of these numbers (kept short here, expanded in literature_sources
# of estimate_cogs result):
# - small_molecule commodity: Anastas & Warner Green Chem; Sheldon E-factor
# - small_molecule pharma:    Federsel HJ., Nat Rev Drug Discov 2005
#                              (doi:10.1038/nrd1799)
# - natural_product extract:  Newman & Cragg, J Nat Prod 2020 (doi:10.1021/acs.jnatprod.9b01285)
# - peptide SPPS:             Bray BL., Nat Rev Drug Discov 2003 (doi:10.1038/nrd1133)
# - peptide recombinant:      Sanchez-Garcia, Microb Cell Fact 2016 (doi:10.1186/s12934-016-0437-3)
# - antibody:                 Walsh G., Nat Biotechnol 2018 (doi:10.1038/nbt.4313);
#                              Kelley B., MAbs 2009 (doi:10.4161/mabs.1.5.9448)
# - enzyme industrial:        Chapman et al., Catalysts 2018 (doi:10.3390/catal8060238)
# ----------------------------------------------------------------------

_ANCHORS: Dict[Tuple[str, str, str], Dict[str, Any]] = {
    # --- Small molecules ---
    ("small_molecule", "volatile", "biotechnological"): {
        "base_low": 0.5, "base_high": 2.5,
        "breakdown": {"rm": 0.45, "dsp": 0.20, "opex": 0.30, "capex": 0.05},
        "anchor_de": "Bioethanol-Industrie (Bulk-Fermentation)",
        "anchor_en": "industrial bioethanol fermentation",
    },
    ("small_molecule", "volatile", "chemical"): {
        "base_low": 1.0, "base_high": 5.0,
        "breakdown": {"rm": 0.55, "dsp": 0.15, "opex": 0.25, "capex": 0.05},
        "anchor_de": "Solvens-Commodity (petrochemisch)",
        "anchor_en": "commodity solvent (petrochemical)",
    },
    ("small_molecule", "non_volatile", "chemical"): {
        # API range: aspirin (~10 €/kg) to mid-tier oncology API (~5000 €/kg)
        "base_low": 5.0, "base_high": 80.0,
        "breakdown": {"rm": 0.40, "dsp": 0.30, "opex": 0.25, "capex": 0.05},
        "anchor_de": "Commodity-API (z. B. Aspirin, Paracetamol, Ibuprofen)",
        "anchor_en": "commodity API (e.g. aspirin, paracetamol, ibuprofen)",
    },
    ("small_molecule", "non_volatile", "biotechnological"): {
        "base_low": 30.0, "base_high": 250.0,
        "breakdown": {"rm": 0.25, "dsp": 0.35, "opex": 0.35, "capex": 0.05},
        "anchor_de": "biotech Small-Molecule (Vanillin-Klasse)",
        "anchor_en": "biotech small molecule (vanillin class)",
    },

    # --- Natural products ---
    ("natural_product", "alkaloid", "chemical"): {
        "base_low": 50.0, "base_high": 600.0,
        "breakdown": {"rm": 0.30, "dsp": 0.40, "opex": 0.25, "capex": 0.05},
        "anchor_de": "Alkaloid-Synthese (Caffeine-Range bis Codein)",
        "anchor_en": "alkaloid synthesis (caffeine to codeine range)",
    },
    ("natural_product", "alkaloid", "extraction"): {
        "base_low": 100.0, "base_high": 2000.0,
        "breakdown": {"rm": 0.20, "dsp": 0.50, "opex": 0.25, "capex": 0.05},
        "anchor_de": "Alkaloid-Extraktion (Pflanzenmaterial-abhängig)",
        "anchor_en": "alkaloid extraction (plant-material dependent)",
    },
    ("natural_product", "alkaloid", "biotechnological"): {
        "base_low": 200.0, "base_high": 1500.0,
        "breakdown": {"rm": 0.20, "dsp": 0.40, "opex": 0.35, "capex": 0.05},
        "anchor_de": "fermentative Alkaloid-Produktion (Stamm-Engineering)",
        "anchor_en": "fermentative alkaloid production (strain-engineered)",
    },
    ("natural_product", "terpene", "biotechnological"): {
        "base_low": 50.0, "base_high": 800.0,
        "breakdown": {"rm": 0.20, "dsp": 0.35, "opex": 0.40, "capex": 0.05},
        "anchor_de": "Terpenoid-Fermentation (Astaxanthin, β-Carotin)",
        "anchor_en": "terpenoid fermentation (astaxanthin, β-carotene)",
    },
    ("natural_product", "terpene", "extraction"): {
        "base_low": 30.0, "base_high": 400.0,
        "breakdown": {"rm": 0.25, "dsp": 0.45, "opex": 0.25, "capex": 0.05},
        "anchor_de": "Terpenoid-Extraktion (Pflanze/Algen)",
        "anchor_en": "terpenoid extraction (plant/algae)",
    },

    # --- Peptides ---
    ("peptide", "linear", "chemical"): {
        # SPPS: depends heavily on length (couplings)
        "base_low": 500.0, "base_high": 8000.0,
        "breakdown": {"rm": 0.50, "dsp": 0.30, "opex": 0.15, "capex": 0.05},
        "anchor_de": "SPPS-Peptid-Synthese (Sequenzlängen-abhängig)",
        "anchor_en": "SPPS peptide synthesis (sequence-length dependent)",
    },
    ("peptide", "linear", "biotechnological"): {
        "base_low": 100.0, "base_high": 1500.0,
        "breakdown": {"rm": 0.20, "dsp": 0.40, "opex": 0.35, "capex": 0.05},
        "anchor_de": "rekombinante Peptid-Produktion (E. coli / Hefe)",
        "anchor_en": "recombinant peptide production (E. coli / yeast)",
    },
    ("peptide", "cyclic", "biotechnological"): {
        # NRPS fermentation (cyclosporine class)
        "base_low": 800.0, "base_high": 5000.0,
        "breakdown": {"rm": 0.20, "dsp": 0.45, "opex": 0.30, "capex": 0.05},
        "anchor_de": "NRPS-Fermentation (Cyclosporin-Klasse)",
        "anchor_en": "NRPS fermentation (cyclosporine class)",
    },
    ("peptide", "cyclic", "chemical"): {
        # Cyclisation as final step
        "base_low": 1500.0, "base_high": 15000.0,
        "breakdown": {"rm": 0.45, "dsp": 0.35, "opex": 0.15, "capex": 0.05},
        "anchor_de": "cyclische Peptid-Synthese (SPPS + Cyclisierung)",
        "anchor_en": "cyclic peptide synthesis (SPPS + cyclisation)",
    },

    # --- Proteins ---
    ("protein", "antibody", "biotechnological"): {
        # mAb COGS: published CMC literature puts commercial mAb COGS at
        # $50-200/g = €50k-200k/kg for established processes. High-titer
        # next-gen processes (5+ g/L) push the low end toward €30k/kg.
        # Source: Kelley B., MAbs 2009; Klutz et al., Biotechnol Prog 2015.
        "base_low": 30000.0, "base_high": 150000.0,
        "breakdown": {"rm": 0.15, "dsp": 0.50, "opex": 0.25, "capex": 0.10},
        "anchor_de": "monoklonaler Antikörper (CHO-Fed-Batch, Protein-A-Standard-DSP)",
        "anchor_en": "monoclonal antibody (CHO fed-batch, protein-A standard DSP)",
    },
    ("protein", "enzyme", "biotechnological"): {
        # Industrial enzymes are CHEAP — Novozymes/DSM lipase, amylase
        # typically 5-50 €/kg bulk
        "base_low": 5.0, "base_high": 80.0,
        "breakdown": {"rm": 0.30, "dsp": 0.20, "opex": 0.40, "capex": 0.10},
        "anchor_de": "industrielles Enzym (Detergens/Food-Grade, Bulk)",
        "anchor_en": "industrial enzyme (detergent/food-grade, bulk)",
    },
}

# Method synonyms — engine accepts several variants
_METHOD_NORMALIZE = {
    "chem": "chemical", "chemical synthesis": "chemical",
    "biotech": "biotechnological", "biotechnological synthesis": "biotechnological",
    "extract": "extraction", "extraction-based": "extraction",
}


def _normalize_method(m: str) -> str:
    s = str(m or "").lower().strip()
    return _METHOD_NORMALIZE.get(s, s)


# ----------------------------------------------------------------------
# Adjustment factors — applied on top of anchor
# ----------------------------------------------------------------------

def _purity_multiplier(purity_percent: Optional[float]) -> float:
    """Higher purity targets multiply COGS (more polishing steps)."""
    if purity_percent is None:
        return 1.0
    p = float(purity_percent)
    if p >= 99.9:
        return 1.6   # pharma-grade, e.g. Trastuzumab
    if p >= 99.5:
        return 1.3
    if p >= 99.0:
        return 1.15
    if p >= 97.0:
        return 1.0
    return 0.9       # tech-grade, lower polishing


def _scale_multiplier(scale_kg_per_year: Optional[float]) -> Tuple[float, float]:
    """Scale economies — larger scale = lower €/kg. Returns (low_mult, high_mult)."""
    if scale_kg_per_year is None:
        return 1.0, 1.0
    kg = float(scale_kg_per_year)
    # Apply roughly log-scale economy: each 10× scale-up halves per-kg cost
    # within a reasonable range. Bounded so we don't run away.
    if kg <= 1.0:          # lab/early R&D
        return 2.5, 4.0
    if kg <= 10.0:
        return 1.6, 2.5
    if kg <= 100.0:        # late R&D / clinical
        return 1.2, 1.8
    if kg <= 1000.0:       # pilot / late clinical
        return 1.0, 1.3
    if kg <= 10000.0:      # commercial scale
        return 0.7, 1.0
    return 0.5, 0.8        # ultra-large commodity


def _steps_multiplier(steps: Optional[int]) -> float:
    """Each additional step compounds costs. Smooth function around base=3."""
    if steps is None:
        return 1.0
    n = int(steps)
    if n <= 1:
        return 0.7
    if n <= 3:
        return 1.0
    if n <= 5:
        return 1.2
    if n <= 10:
        return 1.5
    # SPPS-like — every additional coupling ~5 % extra; cap at 3.0x
    return min(3.0, 1.5 + 0.05 * (n - 10))


def _raw_material_multiplier(eur_per_kg: Optional[float],
                             base_low: float,
                             rm_fraction: float) -> float:
    """User-entered RM price overrides the embedded benchmark when given.

    The RM share of COGS is `rm_fraction`. If the user's RM price is
    much higher than the implicit one (=base_low * rm_fraction), the
    total goes up proportionally.
    """
    if eur_per_kg is None:
        return 1.0
    user_rm = float(eur_per_kg)
    if user_rm <= 0:
        return 1.0
    # Implicit RM cost embedded in the low anchor
    implicit_rm = base_low * rm_fraction
    if implicit_rm <= 0:
        return 1.0
    # Cap effect to avoid runaway numbers from outlier inputs
    ratio = max(0.3, min(5.0, user_rm / implicit_rm))
    # Apply only to the RM portion — i.e. weighted by rm_fraction
    return 1.0 + (ratio - 1.0) * rm_fraction


# ----------------------------------------------------------------------
# Public entry point
# ----------------------------------------------------------------------

# Bug H8 fix: per-molecule overrides for cases where the (type, subtype, method)
# tuple gives a misleading range. Sources cited in production.
#
# Format: molecule_name (lowercase, normalised) -> (low €/kg, high €/kg,
#         anchor_de, anchor_en, breakdown_dict_or_None)
_MOLECULE_COGS_OVERRIDES: Dict[str, Dict[str, Any]] = {
    # Caffeine: USP/Ph.Eur. bulk-synthetic caffeine from China/India is
    # commercially priced at 8–25 €/kg (FOB list prices 2024–2025). The
    # generic alkaloid_chemical anchor (50–600 €/kg) is set for rare
    # alkaloids like Morphine/Codeine and over-estimates Caffeine by
    # 1–2 orders of magnitude.
    # Source: Maier, Ullmann's Encyclopedia of Industrial Chemistry, 7th
    # ed., Wiley-VCH 2011 — Purines and Caffeine
    # (doi:10.1002/14356007.a07_513.pub2).
    "caffeine": {
        "base_low": 8.0, "base_high": 30.0,
        "breakdown": {"rm": 0.50, "dsp": 0.25, "opex": 0.20, "capex": 0.05},
        "anchor_de": "Synthese-Caffeine (USP/Ph.Eur. Bulk, kommerziell etabliert)",
        "anchor_en": "synthetic caffeine (USP/Ph.Eur. bulk, commercially established)",
    },
    # Linalool: world-commodity (~30.000 t/a) priced at 10–30 €/kg in bulk
    # (BASF/Symrise Citral-based semi-synthetic route). The generic
    # terpene anchor (extraction or biotech) over-estimates by 5–50×.
    # Source: BASF/Symrise public price benchmarks 2023; Bauer, Garbe,
    # Surburg, Common Fragrance and Flavor Materials, 5th ed., Wiley-VCH
    # 2006 (ISBN 978-3-527-31150-9).
    "linalool": {
        "base_low": 10.0, "base_high": 50.0,
        "breakdown": {"rm": 0.55, "dsp": 0.15, "opex": 0.25, "capex": 0.05},
        "anchor_de": "Linalool Welt-Commodity (Citral-basiert semi-synthetisch)",
        "anchor_en": "linalool world commodity (citral-based semi-synthetic)",
    },
    # Liraglutide: full SPPS for a 30-mer + γ-Glu-C16 lipidation is at
    # the high end of peptide COGS, around 50–500 k€/kg API for classic
    # SPPS and 10–50 k€/kg for semi-synthetic (recombinant backbone +
    # chemical acylation). The default peptide/linear/chemical anchor
    # (500–8000 €/kg) is correct for short peptides but 1–2 orders of
    # magnitude too low here.
    # Source: Knudsen et al., J. Med. Chem. 2000 (doi:10.1021/jm9909645);
    # Bray BL., Nat. Rev. Drug Discov. 2003 (doi:10.1038/nrd1133).
    "liraglutide": {
        "base_low": 50000.0, "base_high": 500000.0,
        "breakdown": {"rm": 0.40, "dsp": 0.40, "opex": 0.15, "capex": 0.05},
        "anchor_de": "GLP-1-Analogon (30-mer + C16-Lipidierung, SPPS)",
        "anchor_en": "GLP-1 analogue (30-mer + C16 lipidation, SPPS)",
    },
}


def estimate_cogs(process_input: Dict[str, Any]) -> Dict[str, Any]:
    """Estimate COGS in €/kg for the given process inputs.

    Returns a dict structured for the PDF/UI consumer. See module docstring.
    """
    if not isinstance(process_input, dict):
        return _unknown_result("Invalid input")

    mtype = str(process_input.get("molecule_type") or "").lower().strip()
    msub = str(process_input.get("molecule_subtype") or "").lower().strip()
    method = _normalize_method(process_input.get("method"))
    mname = str(process_input.get("molecule_name") or "").lower().strip()

    # Per-molecule override takes precedence over the generic cluster
    # (Bug H8 fix: avoids misleading ranges for outlier molecules).
    anchor = _MOLECULE_COGS_OVERRIDES.get(mname)
    is_override = anchor is not None
    if anchor is None:
        anchor = _ANCHORS.get((mtype, msub, method))
    confidence = "high"
    fallback_note_de: Optional[str] = None
    fallback_note_en: Optional[str] = None

    if anchor is None:
        # Try anchor with method swap (chemical→biotechnological etc.)
        # so we still produce something useful — flag low confidence.
        for alt_method in ("biotechnological", "chemical", "extraction"):
            anchor = _ANCHORS.get((mtype, msub, alt_method))
            if anchor is not None:
                confidence = "low"
                fallback_note_de = (
                    f"Keine exakte Benchmark für Methode '{method}' verfügbar — "
                    f"Schätzung basiert auf '{alt_method}'-Referenz, große Unsicherheit."
                )
                fallback_note_en = (
                    f"No exact benchmark for method '{method}' — estimate based "
                    f"on '{alt_method}' reference, large uncertainty."
                )
                break

    if anchor is None:
        return _unknown_result(
            f"No benchmark for {mtype}/{msub}/{method}",
        )

    base_low = float(anchor["base_low"])
    base_high = float(anchor["base_high"])
    breakdown = dict(anchor["breakdown"])
    rm_share = float(breakdown.get("rm", 0.4))

    # ----- Apply modifiers -----
    # Bug B4 fix: per-molecule overrides (Caffeine, Linalool, Liraglutide)
    # are calibrated for a representative real-world process and already
    # implicitly include "typical" purity, step-count and RM-price
    # assumptions. Stacking the full multiplier chain on top of an
    # override pushes the range out of the empirically validated band
    # (Caffeine 8–30 €/kg → 32–180 €/kg). When the override is active,
    # damp the purity / steps / raw-material multipliers strongly. Scale
    # economy stays in effect (it's a real, monotonic cost driver).
    pm = _purity_multiplier(process_input.get("desired_purity_percent"))
    sm_low, sm_high = _scale_multiplier(process_input.get("scale_kg_per_year"))
    stm = _steps_multiplier(process_input.get("number_of_steps"))
    rmm = _raw_material_multiplier(
        process_input.get("raw_material_cost_eur_per_kg"),
        base_low, rm_share,
    )
    if is_override:
        # Compress everything except scale toward 1.0 by ~75 %.
        pm = 1.0 + 0.25 * (pm - 1.0)
        stm = 1.0 + 0.25 * (stm - 1.0)
        rmm = 1.0 + 0.25 * (rmm - 1.0)

    low = base_low * pm * sm_low * stm * rmm
    high = base_high * pm * sm_high * stm * rmm

    # Cap compounding — prevent runaway upper bounds when multiple
    # multipliers stack (e.g. 30-step SPPS + €3000/kg RM + 99.5 % purity
    # would otherwise reach 25–30× base which exceeds real-world
    # observations). Override anchors get a tighter cap.
    _MAX_COMPOUND_LOW = 1.5 if is_override else 4.0
    _MAX_COMPOUND_HIGH = 2.0 if is_override else 6.0
    low = min(low, base_low * _MAX_COMPOUND_LOW)
    high = min(high, base_high * _MAX_COMPOUND_HIGH)
    # Ensure low <= high after capping
    if low > high:
        low, high = high, low

    # Confidence demoted to medium if many adjustments compound
    extreme_modifiers = [
        pm > 1.4,
        max(sm_low, sm_high) > 3.0 or min(sm_low, sm_high) < 0.6,
        stm > 1.8,
        rmm > 1.5 or rmm < 0.7,
    ]
    if confidence == "high" and sum(1 for x in extreme_modifiers if x) >= 2:
        confidence = "medium"

    # Build human-readable driver list (for PDF/UI)
    drivers_de: List[str] = []
    drivers_en: List[str] = []
    pct_in = process_input.get("desired_purity_percent")
    if isinstance(pct_in, (int, float)) and float(pct_in) >= 99.0:
        drivers_de.append(f"Reinheit {float(pct_in):.1f} % × {pm:.2f}")
        drivers_en.append(f"purity {float(pct_in):.1f} % × {pm:.2f}")
    kg_in = process_input.get("scale_kg_per_year")
    if isinstance(kg_in, (int, float)):
        drivers_de.append(f"Maßstab {float(kg_in):.0f} kg/Jahr × {sm_low:.2f}–{sm_high:.2f}")
        drivers_en.append(f"scale {float(kg_in):.0f} kg/year × {sm_low:.2f}–{sm_high:.2f}")
    steps_in = process_input.get("number_of_steps")
    if isinstance(steps_in, int) and stm != 1.0:
        drivers_de.append(f"{steps_in} Schritte × {stm:.2f}")
        drivers_en.append(f"{steps_in} steps × {stm:.2f}")
    eur_in = process_input.get("raw_material_cost_eur_per_kg")
    if isinstance(eur_in, (int, float)) and rmm != 1.0:
        drivers_de.append(f"Rohstoff {float(eur_in):.0f} €/kg × {rmm:.2f}")
        drivers_en.append(f"raw material {float(eur_in):.0f} €/kg × {rmm:.2f}")

    notes_de = f"Referenz: {anchor['anchor_de']}."
    notes_en = f"Reference: {anchor['anchor_en']}."
    if fallback_note_de:
        notes_de += " " + fallback_note_de
        notes_en += " " + fallback_note_en  # type: ignore

    return {
        "low_eur_per_kg": round(low, 1),
        "high_eur_per_kg": round(high, 1),
        "breakdown": breakdown,
        "confidence": confidence,
        "notes_de": notes_de,
        "notes_en": notes_en,
        "anchor_de": anchor["anchor_de"],
        "anchor_en": anchor["anchor_en"],
        "drivers_de": drivers_de,
        "drivers_en": drivers_en,
        "applied_multipliers": {
            "purity": round(pm, 2),
            "scale_low": round(sm_low, 2),
            "scale_high": round(sm_high, 2),
            "steps": round(stm, 2),
            "raw_material": round(rmm, 2),
        },
    }


def _unknown_result(reason: str) -> Dict[str, Any]:
    """Honest empty-state when we have no benchmark."""
    return {
        "low_eur_per_kg": None,
        "high_eur_per_kg": None,
        "breakdown": {},
        "confidence": "low",
        "notes_de": f"Keine COGS-Schätzung verfügbar ({reason}).",
        "notes_en": f"No COGS estimate available ({reason}).",
        "anchor_de": None,
        "anchor_en": None,
        "drivers_de": [],
        "drivers_en": [],
        "applied_multipliers": {},
    }


def format_cogs_range(result: Dict[str, Any], lang: str = "de") -> Optional[str]:
    """Format the range as 'low – high €/kg' or None if unavailable."""
    lo = result.get("low_eur_per_kg")
    hi = result.get("high_eur_per_kg")
    if lo is None or hi is None:
        return None
    # Round to 2 sig figs for human readability
    def _fmt(x: float) -> str:
        if x >= 1000:
            return f"{x:,.0f}".replace(",", ".")
        if x >= 10:
            return f"{x:.0f}"
        return f"{x:.1f}"
    return f"{_fmt(lo)} – {_fmt(hi)} €/kg"
