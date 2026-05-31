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

# Per-molecule COGS overrides — calibrated against published industry
# benchmarks and commercial bulk prices. Avoid the generic-cluster trap
# where one anchor (e.g. "Commodity-API: Aspirin–Ibuprofen") spans
# 16×, smearing the COGS-range over an unhelpful span. The override
# replaces the cluster anchor with a per-molecule one (1.5–5× span);
# the multiplier-damping logic further constrains modifier stacking.
#
# Format per entry:
#   base_low / base_high : €/kg at typical (mid-scale, mid-purity) inputs
#   breakdown            : RM / DSP / OPEX / CAPEX fractions (must sum to 1.0)
#   anchor_de / anchor_en: short reference label shown in the PDF
#
# Price sources (compact):
#   - Bulk solvents: ICIS price reports; Sasol/BASF/Dow public lists 2023-2025
#   - Commodity APIs: Federsel HJ., Nat. Rev. Drug Discov. 2005 (doi:10.1038/nrd1799)
#   - Biotech APIs: Walsh G., Nat. Biotechnol. 2018 (doi:10.1038/nbt.4313)
#   - mAbs: Kelley B., MAbs 2009 (doi:10.4161/mabs.1.5.9448)
#   - SPPS-peptides: Bray BL., Nat. Rev. Drug Discov. 2003 (doi:10.1038/nrd1133)
#   - NRPS peptides: Marahiel MA., Chem. Rev. 1997 (doi:10.1021/cr960029e)
#   - Industrial enzymes: Chapman et al., Catalysts 2018 (doi:10.3390/catal8060238)
#   - Alkaloids: Maier, Ullmann's Encyclopedia 2011 (doi:10.1002/14356007.a07_513.pub2)
#   - Terpenoids: Bauer, Garbe, Surburg, Common Fragrance & Flavor Materials,
#                 5th ed., Wiley-VCH 2006 (ISBN 978-3-527-31150-9)

_MOLECULE_COGS_OVERRIDES: Dict[str, Dict[str, Any]] = {
    # ---------- small_molecule / volatile (bulk solvents) ----------
    "ethanol": {
        "base_low": 0.5, "base_high": 1.5,
        "breakdown": {"rm": 0.45, "dsp": 0.20, "opex": 0.30, "capex": 0.05},
        "anchor_de": "Bioethanol-Bulk (industrielle Fermentation + Destillation)",
        "anchor_en": "industrial bioethanol (fermentation + distillation)",
    },
    "acetone": {
        "base_low": 0.7, "base_high": 1.8,
        "breakdown": {"rm": 0.55, "dsp": 0.15, "opex": 0.25, "capex": 0.05},
        "anchor_de": "Aceton Welt-Commodity (Cumol-Prozess)",
        "anchor_en": "acetone world commodity (cumene process)",
    },
    "ethyl acetate": {
        "base_low": 0.8, "base_high": 2.0,
        "breakdown": {"rm": 0.55, "dsp": 0.15, "opex": 0.25, "capex": 0.05},
        "anchor_de": "Ethylacetat Welt-Commodity (Veresterung)",
        "anchor_en": "ethyl acetate world commodity (esterification)",
    },
    "methanol": {
        "base_low": 0.3, "base_high": 0.8,
        "breakdown": {"rm": 0.50, "dsp": 0.10, "opex": 0.35, "capex": 0.05},
        "anchor_de": "Methanol Welt-Commodity (Erdgas-Reformierung)",
        "anchor_en": "methanol world commodity (natural-gas reforming)",
    },
    "isopropanol": {
        "base_low": 0.7, "base_high": 1.6,
        "breakdown": {"rm": 0.55, "dsp": 0.15, "opex": 0.25, "capex": 0.05},
        "anchor_de": "Isopropanol Welt-Commodity (Propylen-Hydratation)",
        "anchor_en": "isopropanol world commodity (propylene hydration)",
    },
    "butanol": {
        "base_low": 1.0, "base_high": 2.5,
        "breakdown": {"rm": 0.55, "dsp": 0.15, "opex": 0.25, "capex": 0.05},
        "anchor_de": "n-Butanol Welt-Commodity (Oxo-Synthese)",
        "anchor_en": "n-butanol world commodity (oxo synthesis)",
    },
    "isobutanol": {
        "base_low": 1.2, "base_high": 3.0,
        "breakdown": {"rm": 0.55, "dsp": 0.15, "opex": 0.25, "capex": 0.05},
        "anchor_de": "Isobutanol Welt-Commodity (Oxo-Synthese)",
        "anchor_en": "isobutanol world commodity (oxo synthesis)",
    },
    "2,3-butanediol": {
        "base_low": 2.0, "base_high": 5.0,
        "breakdown": {"rm": 0.40, "dsp": 0.25, "opex": 0.30, "capex": 0.05},
        "anchor_de": "2,3-Butandiol (fermentativ, mittlere Skala)",
        "anchor_en": "2,3-butanediol (fermentative, mid scale)",
    },
    "ethyl lactate": {
        "base_low": 2.0, "base_high": 5.0,
        "breakdown": {"rm": 0.55, "dsp": 0.20, "opex": 0.20, "capex": 0.05},
        "anchor_de": "Ethyllactat (grünes Lösungsmittel, Veresterung)",
        "anchor_en": "ethyl lactate (green solvent, esterification)",
    },
    "diacetyl": {
        "base_low": 5.0, "base_high": 15.0,
        "breakdown": {"rm": 0.45, "dsp": 0.30, "opex": 0.20, "capex": 0.05},
        "anchor_de": "Diacetyl (Food-Flavor, kleinerer Markt)",
        "anchor_en": "diacetyl (food flavour, smaller market)",
    },

    # ---------- small_molecule / non_volatile (APIs + bulk chemicals) ----------
    "vanillin": {
        "base_low": 10.0, "base_high": 30.0,
        "breakdown": {"rm": 0.45, "dsp": 0.25, "opex": 0.25, "capex": 0.05},
        "anchor_de": "Synthese-Vanillin Bulk (Guaiacol/Glyoxylsäure, Borregaard/Solvay)",
        "anchor_en": "synthetic vanillin bulk (guaiacol/glyoxylic acid, Borregaard/Solvay)",
    },
    "ibuprofen": {
        "base_low": 10.0, "base_high": 25.0,
        "breakdown": {"rm": 0.40, "dsp": 0.30, "opex": 0.25, "capex": 0.05},
        "anchor_de": "Ibuprofen-API (BHC/Boots-Verfahren, Generika-Standard)",
        "anchor_en": "ibuprofen API (BHC/Boots process, generic standard)",
    },
    "glucose": {
        "base_low": 0.3, "base_high": 1.0,
        "breakdown": {"rm": 0.60, "dsp": 0.15, "opex": 0.20, "capex": 0.05},
        "anchor_de": "Glucose Welt-Commodity (Mais-/Weizenstärke-Hydrolyse)",
        "anchor_en": "glucose world commodity (corn/wheat starch hydrolysis)",
    },
    "aspirin": {
        "base_low": 5.0, "base_high": 15.0,
        "breakdown": {"rm": 0.45, "dsp": 0.25, "opex": 0.25, "capex": 0.05},
        "anchor_de": "Aspirin-API (Salicylsäure-Acetylierung, Generika-Standard)",
        "anchor_en": "aspirin API (salicylic acid acetylation, generic standard)",
    },
    "citric acid": {
        "base_low": 1.0, "base_high": 3.0,
        "breakdown": {"rm": 0.45, "dsp": 0.25, "opex": 0.25, "capex": 0.05},
        "anchor_de": "Citronensäure Welt-Commodity (Aspergillus niger Fermentation)",
        "anchor_en": "citric acid world commodity (Aspergillus niger fermentation)",
    },
    "lactic acid": {
        "base_low": 1.5, "base_high": 4.0,
        "breakdown": {"rm": 0.45, "dsp": 0.25, "opex": 0.25, "capex": 0.05},
        "anchor_de": "Milchsäure Bulk (Lactobacillus-Fermentation)",
        "anchor_en": "lactic acid bulk (Lactobacillus fermentation)",
    },
    "succinic acid": {
        "base_low": 2.0, "base_high": 6.0,
        "breakdown": {"rm": 0.45, "dsp": 0.25, "opex": 0.25, "capex": 0.05},
        "anchor_de": "Bernsteinsäure (Biobasierte Plattform, Reverdia/Myriant)",
        "anchor_en": "succinic acid (bio-based platform, Reverdia/Myriant)",
    },
    "itaconic acid": {
        "base_low": 4.0, "base_high": 12.0,
        "breakdown": {"rm": 0.40, "dsp": 0.30, "opex": 0.25, "capex": 0.05},
        "anchor_de": "Itakonsäure (Aspergillus terreus Fermentation)",
        "anchor_en": "itaconic acid (Aspergillus terreus fermentation)",
    },
    "paracetamol": {
        "base_low": 5.0, "base_high": 12.0,
        "breakdown": {"rm": 0.45, "dsp": 0.25, "opex": 0.25, "capex": 0.05},
        "anchor_de": "Paracetamol-API (p-Aminophenol-Acetylierung, Generika)",
        "anchor_en": "paracetamol API (p-aminophenol acetylation, generic)",
    },

    # ---------- natural_product / alkaloid ----------
    "caffeine": {
        "base_low": 8.0, "base_high": 30.0,
        "breakdown": {"rm": 0.50, "dsp": 0.25, "opex": 0.20, "capex": 0.05},
        "anchor_de": "Synthese-Caffeine (USP/Ph.Eur. Bulk, Traube-Synthese)",
        "anchor_en": "synthetic caffeine (USP/Ph.Eur. bulk, Traube synthesis)",
    },
    "morphine": {
        "base_low": 200.0, "base_high": 800.0,
        "breakdown": {"rm": 0.30, "dsp": 0.40, "opex": 0.25, "capex": 0.05},
        "anchor_de": "Morphin-API (P. somniferum-Extraktion, BtMG-/DEA-lizenziert)",
        "anchor_en": "morphine API (P. somniferum extraction, controlled-substance)",
    },
    "codeine": {
        "base_low": 150.0, "base_high": 600.0,
        "breakdown": {"rm": 0.30, "dsp": 0.40, "opex": 0.25, "capex": 0.05},
        "anchor_de": "Codein-API (Morphin-Methylierung oder direkte Extraktion)",
        "anchor_en": "codeine API (morphine methylation or direct extraction)",
    },
    "atropine": {
        "base_low": 500.0, "base_high": 2500.0,
        "breakdown": {"rm": 0.30, "dsp": 0.40, "opex": 0.25, "capex": 0.05},
        "anchor_de": "Atropin-API (Datura/Atropa-Extraktion oder Semisynthese)",
        "anchor_en": "atropine API (Datura/Atropa extraction or semi-synthesis)",
    },
    "quinine": {
        "base_low": 40.0, "base_high": 200.0,
        "breakdown": {"rm": 0.40, "dsp": 0.35, "opex": 0.20, "capex": 0.05},
        "anchor_de": "Chinin-API (Cinchona-Extraktion, Bulk)",
        "anchor_en": "quinine API (Cinchona extraction, bulk)",
    },

    # ---------- natural_product / terpene ----------
    "linalool": {
        "base_low": 10.0, "base_high": 50.0,
        "breakdown": {"rm": 0.55, "dsp": 0.15, "opex": 0.25, "capex": 0.05},
        "anchor_de": "Linalool Welt-Commodity (Citral-basiert semi-synthetisch, BASF/Symrise)",
        "anchor_en": "linalool world commodity (citral-based semi-synthetic, BASF/Symrise)",
    },
    "citral": {
        "base_low": 12.0, "base_high": 40.0,
        "breakdown": {"rm": 0.55, "dsp": 0.15, "opex": 0.25, "capex": 0.05},
        "anchor_de": "Citral Welt-Commodity (BASF/Symrise, semi-synthetisch)",
        "anchor_en": "citral world commodity (BASF/Symrise, semi-synthetic)",
    },
    "limonene": {
        "base_low": 5.0, "base_high": 25.0,
        "breakdown": {"rm": 0.50, "dsp": 0.20, "opex": 0.25, "capex": 0.05},
        "anchor_de": "Limonen (Zitrusöl-Nebenprodukt, F&F-Bulk)",
        "anchor_en": "limonene (citrus oil by-product, F&F bulk)",
    },
    "menthol": {
        "base_low": 15.0, "base_high": 50.0,
        "breakdown": {"rm": 0.50, "dsp": 0.20, "opex": 0.25, "capex": 0.05},
        "anchor_de": "Menthol (BASF/Symrise/Takasago, Citronellal-Route)",
        "anchor_en": "menthol (BASF/Symrise/Takasago, citronellal route)",
    },
    "geraniol": {
        "base_low": 20.0, "base_high": 60.0,
        "breakdown": {"rm": 0.55, "dsp": 0.20, "opex": 0.20, "capex": 0.05},
        "anchor_de": "Geraniol (semi-synthetisch, Citronellol-Oxidation)",
        "anchor_en": "geraniol (semi-synthetic, citronellol oxidation)",
    },
    "α-pinene": {
        "base_low": 4.0, "base_high": 15.0,
        "breakdown": {"rm": 0.55, "dsp": 0.20, "opex": 0.20, "capex": 0.05},
        "anchor_de": "α-Pinen (Terpentinöl-Nebenprodukt der Zellstoff-Industrie)",
        "anchor_en": "α-pinene (turpentine by-product of the pulp industry)",
    },
    "β-carotene": {
        "base_low": 200.0, "base_high": 1500.0,
        "breakdown": {"rm": 0.35, "dsp": 0.35, "opex": 0.25, "capex": 0.05},
        "anchor_de": "β-Carotin (DSM/BASF synthetisch oder D. salina-Extraktion)",
        "anchor_en": "β-carotene (DSM/BASF synthetic or D. salina extraction)",
    },
    "astaxanthin": {
        "base_low": 500.0, "base_high": 3000.0,
        "breakdown": {"rm": 0.30, "dsp": 0.45, "opex": 0.20, "capex": 0.05},
        "anchor_de": "Astaxanthin (H. pluvialis Mikroalgen oder P. rhodozyma-Fermentation)",
        "anchor_en": "astaxanthin (H. pluvialis microalgae or P. rhodozyma fermentation)",
    },
    "squalene": {
        "base_low": 30.0, "base_high": 200.0,
        "breakdown": {"rm": 0.45, "dsp": 0.30, "opex": 0.20, "capex": 0.05},
        "anchor_de": "Squalen (Amaranth-Öl oder Hai-Leber-historisch; heute mikrobiell)",
        "anchor_en": "squalene (amaranth oil or shark-liver-historic; now microbial)",
    },
    "artemisinin": {
        "base_low": 200.0, "base_high": 800.0,
        "breakdown": {"rm": 0.35, "dsp": 0.40, "opex": 0.20, "capex": 0.05},
        "anchor_de": "Artemisinin (Amyris/Sanofi semi-synthetisch via Hefe + A. annua)",
        "anchor_en": "artemisinin (Amyris/Sanofi semi-synthetic via yeast + A. annua)",
    },

    # ---------- peptide / linear (SPPS or recombinant) ----------
    "glutathione": {
        "base_low": 100.0, "base_high": 500.0,
        "breakdown": {"rm": 0.40, "dsp": 0.35, "opex": 0.20, "capex": 0.05},
        "anchor_de": "Glutathion Bulk (Hefe-Fermentation, Kohjin/Tanabe)",
        "anchor_en": "glutathione bulk (yeast fermentation, Kohjin/Tanabe)",
    },
    "leu-enkephalin": {
        "base_low": 5000.0, "base_high": 30000.0,
        "breakdown": {"rm": 0.50, "dsp": 0.30, "opex": 0.15, "capex": 0.05},
        "anchor_de": "Leu-Enkephalin (5-mer SPPS, Forschungs-/Spezialitäten-Grade)",
        "anchor_en": "leu-enkephalin (5-mer SPPS, research / specialty grade)",
    },
    "met-enkephalin": {
        "base_low": 5000.0, "base_high": 30000.0,
        "breakdown": {"rm": 0.50, "dsp": 0.30, "opex": 0.15, "capex": 0.05},
        "anchor_de": "Met-Enkephalin (5-mer SPPS, Forschungs-/Spezialitäten-Grade)",
        "anchor_en": "met-enkephalin (5-mer SPPS, research / specialty grade)",
    },
    "bradykinin": {
        "base_low": 8000.0, "base_high": 50000.0,
        "breakdown": {"rm": 0.50, "dsp": 0.30, "opex": 0.15, "capex": 0.05},
        "anchor_de": "Bradykinin (9-mer SPPS, Forschungs-/Spezialitäten-Grade)",
        "anchor_en": "bradykinin (9-mer SPPS, research / specialty grade)",
    },
    "angiotensin ii": {
        "base_low": 10000.0, "base_high": 60000.0,
        "breakdown": {"rm": 0.50, "dsp": 0.30, "opex": 0.15, "capex": 0.05},
        "anchor_de": "Angiotensin II (8-mer SPPS, Giapreza-API ähnlich)",
        "anchor_en": "angiotensin II (8-mer SPPS, Giapreza-like API)",
    },
    "glucagon": {
        "base_low": 50000.0, "base_high": 300000.0,
        "breakdown": {"rm": 0.45, "dsp": 0.35, "opex": 0.15, "capex": 0.05},
        "anchor_de": "Glucagon-API (29-mer SPPS oder rekombinant, Lilly/Novo)",
        "anchor_en": "glucagon API (29-mer SPPS or recombinant, Lilly/Novo)",
    },
    "calcitonin": {
        "base_low": 100000.0, "base_high": 500000.0,
        "breakdown": {"rm": 0.45, "dsp": 0.40, "opex": 0.10, "capex": 0.05},
        "anchor_de": "Lachs-Calcitonin (32-mer SPPS, Miacalcic-API)",
        "anchor_en": "salmon calcitonin (32-mer SPPS, Miacalcic API)",
    },
    "liraglutide": {
        "base_low": 50000.0, "base_high": 500000.0,
        "breakdown": {"rm": 0.40, "dsp": 0.40, "opex": 0.15, "capex": 0.05},
        "anchor_de": "Liraglutide (GLP-1-Analog, SPPS 30-mer + C16-Lipidierung)",
        "anchor_en": "liraglutide (GLP-1 analogue, SPPS 30-mer + C16 lipidation)",
    },
    "teriparatide": {
        "base_low": 200000.0, "base_high": 800000.0,
        "breakdown": {"rm": 0.45, "dsp": 0.40, "opex": 0.10, "capex": 0.05},
        "anchor_de": "Teriparatide (PTH 1–34, SPPS oder rekombinant, Forteo-API)",
        "anchor_en": "teriparatide (PTH 1–34, SPPS or recombinant, Forteo API)",
    },
    "exenatide": {
        "base_low": 50000.0, "base_high": 300000.0,
        "breakdown": {"rm": 0.45, "dsp": 0.40, "opex": 0.10, "capex": 0.05},
        "anchor_de": "Exenatide (GLP-1-Analog, SPPS 39-mer, Byetta-API)",
        "anchor_en": "exenatide (GLP-1 analogue, SPPS 39-mer, Byetta API)",
    },
    "insulin": {
        "base_low": 500.0, "base_high": 3000.0,
        "breakdown": {"rm": 0.30, "dsp": 0.45, "opex": 0.20, "capex": 0.05},
        "anchor_de": "Insulin-API (rekombinant E. coli/Hefe, Eli Lilly/Novo Nordisk)",
        "anchor_en": "insulin API (recombinant E. coli/yeast, Lilly/Novo Nordisk)",
    },

    # ---------- peptide / cyclic (NRPS or SPPS+cyclisation) ----------
    "cyclosporine": {
        "base_low": 3000.0, "base_high": 30000.0,
        "breakdown": {"rm": 0.20, "dsp": 0.50, "opex": 0.25, "capex": 0.05},
        "anchor_de": "Ciclosporin (T. inflatum NRPS-Fermentation, Sandimmun-API)",
        "anchor_en": "cyclosporine (T. inflatum NRPS fermentation, Sandimmune API)",
    },
    "gramicidin s": {
        "base_low": 5000.0, "base_high": 25000.0,
        "breakdown": {"rm": 0.25, "dsp": 0.50, "opex": 0.20, "capex": 0.05},
        "anchor_de": "Gramicidin S (Bacillus brevis Fermentation)",
        "anchor_en": "gramicidin S (Bacillus brevis fermentation)",
    },
    "bacitracin": {
        "base_low": 1000.0, "base_high": 5000.0,
        "breakdown": {"rm": 0.25, "dsp": 0.50, "opex": 0.20, "capex": 0.05},
        "anchor_de": "Bacitracin (Bacillus licheniformis-Fermentation, Veterinär-Bulk)",
        "anchor_en": "bacitracin (Bacillus licheniformis fermentation, veterinary bulk)",
    },
    "vancomycin": {
        "base_low": 1000.0, "base_high": 8000.0,
        "breakdown": {"rm": 0.25, "dsp": 0.50, "opex": 0.20, "capex": 0.05},
        "anchor_de": "Vancomycin-API (Amycolatopsis orientalis Fermentation)",
        "anchor_en": "vancomycin API (Amycolatopsis orientalis fermentation)",
    },
    "daptomycin": {
        "base_low": 50000.0, "base_high": 250000.0,
        "breakdown": {"rm": 0.20, "dsp": 0.55, "opex": 0.20, "capex": 0.05},
        "anchor_de": "Daptomycin-API (Streptomyces roseosporus, Cubicin)",
        "anchor_en": "daptomycin API (Streptomyces roseosporus, Cubicin)",
    },
    "octreotide": {
        "base_low": 100000.0, "base_high": 500000.0,
        "breakdown": {"rm": 0.50, "dsp": 0.35, "opex": 0.10, "capex": 0.05},
        "anchor_de": "Octreotide (8-mer cyclisches SPPS-Peptid, Sandostatin-API)",
        "anchor_en": "octreotide (8-mer cyclic SPPS peptide, Sandostatin API)",
    },
    "polymyxin b": {
        "base_low": 1500.0, "base_high": 8000.0,
        "breakdown": {"rm": 0.25, "dsp": 0.50, "opex": 0.20, "capex": 0.05},
        "anchor_de": "Polymyxin B (B. polymyxa Fermentation, Last-Resort-Antibiotikum)",
        "anchor_en": "polymyxin B (B. polymyxa fermentation, last-resort antibiotic)",
    },
    "colistin": {
        "base_low": 1500.0, "base_high": 8000.0,
        "breakdown": {"rm": 0.25, "dsp": 0.50, "opex": 0.20, "capex": 0.05},
        "anchor_de": "Colistin (Polymyxin E, P. colistinus Fermentation)",
        "anchor_en": "colistin (polymyxin E, P. colistinus fermentation)",
    },
    "caspofungin": {
        "base_low": 30000.0, "base_high": 150000.0,
        "breakdown": {"rm": 0.25, "dsp": 0.55, "opex": 0.15, "capex": 0.05},
        "anchor_de": "Caspofungin (Echinocandin, Cancidas-API)",
        "anchor_en": "caspofungin (echinocandin, Cancidas API)",
    },
    "anidulafungin": {
        "base_low": 30000.0, "base_high": 150000.0,
        "breakdown": {"rm": 0.25, "dsp": 0.55, "opex": 0.15, "capex": 0.05},
        "anchor_de": "Anidulafungin (Echinocandin, Eraxis-API)",
        "anchor_en": "anidulafungin (echinocandin, Eraxis API)",
    },

    # ---------- protein / antibody (mAbs, CHO-Fed-Batch) ----------
    "adalimumab": {
        "base_low": 50000.0, "base_high": 300000.0,
        "breakdown": {"rm": 0.15, "dsp": 0.50, "opex": 0.25, "capex": 0.10},
        "anchor_de": "Adalimumab (anti-TNFα IgG1, Humira-API, CHO-Fed-Batch)",
        "anchor_en": "adalimumab (anti-TNFα IgG1, Humira API, CHO fed-batch)",
    },
    "trastuzumab": {
        "base_low": 50000.0, "base_high": 300000.0,
        "breakdown": {"rm": 0.15, "dsp": 0.50, "opex": 0.25, "capex": 0.10},
        "anchor_de": "Trastuzumab (anti-HER2 IgG1, Herceptin-API, CHO-Fed-Batch)",
        "anchor_en": "trastuzumab (anti-HER2 IgG1, Herceptin API, CHO fed-batch)",
    },
    "rituximab": {
        "base_low": 50000.0, "base_high": 300000.0,
        "breakdown": {"rm": 0.15, "dsp": 0.50, "opex": 0.25, "capex": 0.10},
        "anchor_de": "Rituximab (anti-CD20 IgG1, MabThera-API, CHO-Fed-Batch)",
        "anchor_en": "rituximab (anti-CD20 IgG1, MabThera API, CHO fed-batch)",
    },
    "bevacizumab": {
        "base_low": 50000.0, "base_high": 250000.0,
        "breakdown": {"rm": 0.15, "dsp": 0.50, "opex": 0.25, "capex": 0.10},
        "anchor_de": "Bevacizumab (anti-VEGF IgG1, Avastin-API, CHO-Fed-Batch)",
        "anchor_en": "bevacizumab (anti-VEGF IgG1, Avastin API, CHO fed-batch)",
    },
    "infliximab": {
        "base_low": 50000.0, "base_high": 300000.0,
        "breakdown": {"rm": 0.15, "dsp": 0.50, "opex": 0.25, "capex": 0.10},
        "anchor_de": "Infliximab (chimerischer anti-TNFα, Remicade-API, SP2/0)",
        "anchor_en": "infliximab (chimeric anti-TNFα, Remicade API, SP2/0)",
    },
    "pembrolizumab": {
        "base_low": 100000.0, "base_high": 500000.0,
        "breakdown": {"rm": 0.15, "dsp": 0.50, "opex": 0.25, "capex": 0.10},
        "anchor_de": "Pembrolizumab (anti-PD-1 IgG4, Keytruda-API, Top-Onkologie)",
        "anchor_en": "pembrolizumab (anti-PD-1 IgG4, Keytruda API, top oncology)",
    },
    "nivolumab": {
        "base_low": 100000.0, "base_high": 500000.0,
        "breakdown": {"rm": 0.15, "dsp": 0.50, "opex": 0.25, "capex": 0.10},
        "anchor_de": "Nivolumab (anti-PD-1 IgG4, Opdivo-API, Top-Onkologie)",
        "anchor_en": "nivolumab (anti-PD-1 IgG4, Opdivo API, top oncology)",
    },
    "cetuximab": {
        "base_low": 100000.0, "base_high": 400000.0,
        "breakdown": {"rm": 0.15, "dsp": 0.50, "opex": 0.25, "capex": 0.10},
        "anchor_de": "Cetuximab (chimerischer anti-EGFR, Erbitux-API, SP2/0)",
        "anchor_en": "cetuximab (chimeric anti-EGFR, Erbitux API, SP2/0)",
    },
    "pertuzumab": {
        "base_low": 100000.0, "base_high": 400000.0,
        "breakdown": {"rm": 0.15, "dsp": 0.50, "opex": 0.25, "capex": 0.10},
        "anchor_de": "Pertuzumab (anti-HER2-Dimerisierung, Perjeta-API, CHO)",
        "anchor_en": "pertuzumab (anti-HER2 dimerisation, Perjeta API, CHO)",
    },
    "daratumumab": {
        "base_low": 100000.0, "base_high": 500000.0,
        "breakdown": {"rm": 0.15, "dsp": 0.50, "opex": 0.25, "capex": 0.10},
        "anchor_de": "Daratumumab (anti-CD38, Darzalex-API, CHO-Fed-Batch)",
        "anchor_en": "daratumumab (anti-CD38, Darzalex API, CHO fed-batch)",
    },
    "etanercept": {
        # Fusion protein, formally protein/antibody in this DB
        "base_low": 80000.0, "base_high": 400000.0,
        "breakdown": {"rm": 0.15, "dsp": 0.55, "opex": 0.20, "capex": 0.10},
        "anchor_de": "Etanercept (TNFR-Fc-Fusion, Enbrel-API, CHO)",
        "anchor_en": "etanercept (TNFR-Fc fusion, Enbrel API, CHO)",
    },

    # ---------- protein / enzyme (industrial enzymes + recombinant therapeutics) ----------
    "amylase": {
        "base_low": 10.0, "base_high": 50.0,
        "breakdown": {"rm": 0.30, "dsp": 0.20, "opex": 0.40, "capex": 0.10},
        "anchor_de": "α-Amylase Industrieenzym (Novozymes Termamyl, Bacillus)",
        "anchor_en": "α-amylase industrial enzyme (Novozymes Termamyl, Bacillus)",
    },
    "lipase": {
        "base_low": 30.0, "base_high": 150.0,
        "breakdown": {"rm": 0.30, "dsp": 0.25, "opex": 0.35, "capex": 0.10},
        "anchor_de": "Lipase Industrieenzym (Novozymes Lipozyme/CalB, A. oryzae)",
        "anchor_en": "lipase industrial enzyme (Novozymes Lipozyme/CalB, A. oryzae)",
    },
    "protease": {
        "base_low": 20.0, "base_high": 100.0,
        "breakdown": {"rm": 0.30, "dsp": 0.25, "opex": 0.35, "capex": 0.10},
        "anchor_de": "Protease Industrieenzym (Novozymes/DSM, Detergens-Grade)",
        "anchor_en": "protease industrial enzyme (Novozymes/DSM, detergent grade)",
    },
    "cellulase": {
        "base_low": 20.0, "base_high": 100.0,
        "breakdown": {"rm": 0.30, "dsp": 0.25, "opex": 0.35, "capex": 0.10},
        "anchor_de": "Cellulase (Novozymes Cellic, T. reesei)",
        "anchor_en": "cellulase (Novozymes Cellic, T. reesei)",
    },
    "lactase": {
        "base_low": 50.0, "base_high": 200.0,
        "breakdown": {"rm": 0.30, "dsp": 0.30, "opex": 0.30, "capex": 0.10},
        "anchor_de": "Lactase (β-Galactosidase, K. lactis-Fermentation, Food-Grade)",
        "anchor_en": "lactase (β-galactosidase, K. lactis fermentation, food grade)",
    },
    "glucoamylase": {
        "base_low": 15.0, "base_high": 80.0,
        "breakdown": {"rm": 0.30, "dsp": 0.25, "opex": 0.35, "capex": 0.10},
        "anchor_de": "Glucoamylase (A. niger-Fermentation, Stärke-Industrie)",
        "anchor_en": "glucoamylase (A. niger fermentation, starch industry)",
    },
    "phytase": {
        "base_low": 20.0, "base_high": 100.0,
        "breakdown": {"rm": 0.30, "dsp": 0.25, "opex": 0.35, "capex": 0.10},
        "anchor_de": "Phytase (Tier-Feed-Additiv, A. niger / P. pastoris)",
        "anchor_en": "phytase (animal feed additive, A. niger / P. pastoris)",
    },
    "glucose isomerase": {
        "base_low": 50.0, "base_high": 200.0,
        "breakdown": {"rm": 0.25, "dsp": 0.30, "opex": 0.35, "capex": 0.10},
        "anchor_de": "Glucose-Isomerase (immobilisiert, HFCS-Produktion)",
        "anchor_en": "glucose isomerase (immobilised, HFCS production)",
    },
    "pectinase": {
        "base_low": 30.0, "base_high": 150.0,
        "breakdown": {"rm": 0.30, "dsp": 0.30, "opex": 0.30, "capex": 0.10},
        "anchor_de": "Pektinase (A. niger, Saft-/Wein-Industrie)",
        "anchor_en": "pectinase (A. niger, juice/wine industry)",
    },
    "asparaginase": {
        "base_low": 100.0, "base_high": 500.0,
        "breakdown": {"rm": 0.30, "dsp": 0.40, "opex": 0.20, "capex": 0.10},
        "anchor_de": "Asparaginase (E. coli/E. carotovora, ALL-Pharmazeutikum)",
        "anchor_en": "asparaginase (E. coli/E. carotovora, ALL pharmaceutical)",
    },
    "erythropoietin": {
        # Recombinant therapeutic — priced like a biologic, not like an
        # industrial enzyme. Classified as "enzyme" in the DB but behaves
        # as a biologic for COGS purposes.
        "base_low": 100000.0, "base_high": 500000.0,
        "breakdown": {"rm": 0.15, "dsp": 0.50, "opex": 0.25, "capex": 0.10},
        "anchor_de": "Erythropoietin (rekombinant CHO, Epogen/Eprex-API)",
        "anchor_en": "erythropoietin (recombinant CHO, Epogen/Eprex API)",
    },
    "filgrastim": {
        "base_low": 50000.0, "base_high": 250000.0,
        "breakdown": {"rm": 0.15, "dsp": 0.55, "opex": 0.20, "capex": 0.10},
        "anchor_de": "Filgrastim (G-CSF, rekombinant E. coli, Neupogen-API)",
        "anchor_en": "filgrastim (G-CSF, recombinant E. coli, Neupogen API)",
    },
    "somatropin": {
        "base_low": 30000.0, "base_high": 150000.0,
        "breakdown": {"rm": 0.15, "dsp": 0.55, "opex": 0.20, "capex": 0.10},
        "anchor_de": "Somatropin (hGH, rekombinant E. coli, Genotropin-API)",
        "anchor_en": "somatropin (hGH, recombinant E. coli, Genotropin API)",
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
    cluster_fallback = False  # P4 fix: track when we're using the generic cluster
    if anchor is None:
        anchor = _ANCHORS.get((mtype, msub, method))
        if anchor is not None:
            cluster_fallback = True
    confidence = "high"
    fallback_note_de: Optional[str] = None
    fallback_note_en: Optional[str] = None

    # P4 fix: if molecule not in the calibrated 79er-list, demote confidence
    # to "medium" and surface a banner-style note. This prevents the
    # Vanillin-class bug (wide cluster range) from being read by pilots as
    # a high-confidence per-molecule answer.
    if cluster_fallback and mname:
        confidence = "medium"
        fallback_note_de = (
            f"Cluster-Schätzung für '{mname}' — kein per-Molekül-Anker "
            f"hinterlegt. Spanne ist branchen-typisch für {mtype}/{msub}/"
            f"{method}, nicht molekül-spezifisch kalibriert."
        )
        fallback_note_en = (
            f"Cluster estimate for '{mname}' — no per-molecule anchor "
            f"on file. Range is industry-typical for {mtype}/{msub}/"
            f"{method}, not molecule-specific calibration."
        )

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
    # multipliers stack. Override anchors get a TIGHT cap because they
    # are explicit per-molecule calibrations against published price
    # benchmarks (Bug C3 fix: previous cap of 2.0× pushed Caffeine to
    # 60 €/kg vs target 30; tightened to 1.3× = 30 vs target 30).
    _MAX_COMPOUND_LOW = 1.2 if is_override else 4.0
    _MAX_COMPOUND_HIGH = 1.3 if is_override else 6.0
    low = min(low, base_low * _MAX_COMPOUND_LOW)
    high = min(high, base_high * _MAX_COMPOUND_HIGH)
    # Also floor at base × 0.7 — overrides shouldn't drop too far below
    # the empirical range even at lab scale.
    if is_override:
        low = max(low, base_low * 0.7)
        high = max(high, base_high * 0.7)
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
    # Bug W5 fix: previously only fired at ≥ 99 % purity, which caused
    # PDFs at 95 % (e.g. Linalool) to silently drop the purity modifier
    # from the list — readers couldn't tell whether the modifier had
    # applied or not. Now ALWAYS list the purity modifier so the chain
    # is fully transparent (× 1.00 is shown explicitly at standard
    # purity).
    if isinstance(pct_in, (int, float)):
        _pct_de = f"{float(pct_in):.1f}".replace(".", ",")
        _pm_de = f"{pm:.2f}".replace(".", ",")
        drivers_de.append(f"Reinheit {_pct_de} % × {_pm_de}")
        drivers_en.append(f"purity {float(pct_in):.1f} % × {pm:.2f}")
    kg_in = process_input.get("scale_kg_per_year")
    if isinstance(kg_in, (int, float)):
        # Bug E6 fix: locale-aware thousands formatting so the German
        # report doesn't show "5000 kg/Jahr" next to "5.000 kg/Jahr"
        # earlier in the same document.
        _kg_int = int(round(float(kg_in)))
        kg_de = f"{_kg_int:,}".replace(",", ".")
        kg_en = f"{_kg_int:,}"
        # Comma-based DE decimal for multipliers
        sm_low_de = f"{sm_low:.2f}".replace(".", ",")
        sm_high_de = f"{sm_high:.2f}".replace(".", ",")
        drivers_de.append(f"Maßstab {kg_de} kg/Jahr × {sm_low_de}–{sm_high_de}")
        drivers_en.append(f"scale {kg_en} kg/year × {sm_low:.2f}–{sm_high:.2f}")
    steps_in = process_input.get("number_of_steps")
    if isinstance(steps_in, int) and stm != 1.0:
        _stm_de = f"{stm:.2f}".replace(".", ",")
        drivers_de.append(f"{steps_in} Schritte × {_stm_de}")
        drivers_en.append(f"{steps_in} steps × {stm:.2f}")
    eur_in = process_input.get("raw_material_cost_eur_per_kg")
    if isinstance(eur_in, (int, float)) and rmm != 1.0:
        _eur_int = int(round(float(eur_in)))
        _eur_de = f"{_eur_int:,}".replace(",", ".")
        _eur_en = f"{_eur_int:,}"
        _rmm_de = f"{rmm:.2f}".replace(".", ",")
        drivers_de.append(f"Rohstoff {_eur_de} €/kg × {_rmm_de}")
        drivers_en.append(f"raw material {_eur_en} €/kg × {rmm:.2f}")

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


def estimate_cogs_whatif(base_input: Dict[str, Any],
                         new_scale_kg_per_year: Optional[float] = None,
                         new_purity_percent: Optional[float] = None,
                         new_rm_eur_per_kg: Optional[float] = None) -> Dict[str, Any]:
    """What-if COGS projection for the sensitivity sliders.

    Why this exists separately from estimate_cogs():
      estimate_cogs() is the calibrated POINT estimate shown in the PDF. For
      per-molecule overrides it damps the modifiers to 25 % and hard-caps the
      result at ~1.2×/1.3× of the anchor — correct for the headline number,
      but it makes every slider move plateau instantly (RM 5→500 €/kg changed
      nothing). This function instead projects the COGS RELATIVE to the PDF
      value, applying the *real, smooth, uncapped* economic elasticities so
      the sliders respond the way a process scientist expects.

    Construction (the key property the user asked for):
        COGS_whatif = COGS_pdf(base_input) × f_scale × f_purity × f_rm
      where every factor f == 1.0 exactly at the molecule's original input.
      → Sliders at default reproduce the PDF number byte-for-byte; only a
        deliberate change moves the result, and the move is monotone &
        economically grounded:
          • RM    — linear in the RM cost share of COGS
          • scale — smooth log-degressive economy of scale (~−26 %/decade,
                    calibrated to the engine's own scale bands)
          • purity— progressive "log-impurity" polishing cost (each extra
                    nine of purity adds a constant increment, weighted by the
                    DSP share)

    Returns a dict mirroring estimate_cogs() plus a `whatif_factors` block.
    Falls back to estimate_cogs() output verbatim if no benchmark exists.
    """
    import math as _m
    base = estimate_cogs(base_input) or {}
    low0 = base.get("low_eur_per_kg")
    high0 = base.get("high_eur_per_kg")
    if low0 is None or high0 is None:
        return base  # honest unknown — nothing to project from

    breakdown = base.get("breakdown") or {}
    rm_share = float(breakdown.get("rm", 0.4))
    dsp_share = float(breakdown.get("dsp", 0.25))

    orig_scale = base_input.get("scale_kg_per_year")
    orig_purity = base_input.get("desired_purity_percent")
    orig_rm = base_input.get("raw_material_cost_eur_per_kg")

    # ----- f_scale: smooth log-degressive economy of scale -----
    f_scale = 1.0
    if (new_scale_kg_per_year is not None and orig_scale is not None
            and float(orig_scale) > 0 and float(new_scale_kg_per_year) > 0):
        # k≈0.13 → each 10× scale-up lowers €/kg by ~26 %, matching the
        # engine's stepwise _scale_multiplier bands. Bounded gently so an
        # extreme slider can't run to absurd values.
        ratio = float(orig_scale) / float(new_scale_kg_per_year)
        f_scale = ratio ** 0.13
        f_scale = max(0.3, min(3.0, f_scale))

    # ----- f_purity: progressive log-impurity polishing cost -----
    def _puri_factor(p: float) -> float:
        impurity = max((100.0 - float(p)) / 100.0, 1e-4)  # cap at 99.99 %
        return -_m.log10(impurity)                         # 99 %→2, 99.9 %→3
    f_purity = 1.0
    if new_purity_percent is not None and orig_purity is not None:
        delta = _puri_factor(new_purity_percent) - _puri_factor(orig_purity)
        f_purity = 1.0 + dsp_share * delta
        f_purity = max(0.5, f_purity)  # purity can't drive cost below half

    # ----- f_rm: linear in the RM cost share -----
    f_rm = 1.0
    if (new_rm_eur_per_kg is not None and orig_rm is not None
            and float(orig_rm) > 0 and float(new_rm_eur_per_kg) > 0):
        f_rm = (1.0 - rm_share) + rm_share * (float(new_rm_eur_per_kg) / float(orig_rm))
        f_rm = max(0.4, f_rm)  # RM share can't take total below a floor

    f_total = f_scale * f_purity * f_rm
    low = round(low0 * f_total, 2)
    high = round(high0 * f_total, 2)
    if low > high:
        low, high = high, low

    out = dict(base)
    out["low_eur_per_kg"] = low
    out["high_eur_per_kg"] = high
    out["whatif_factors"] = {
        "scale": round(f_scale, 3),
        "purity": round(f_purity, 3),
        "raw_material": round(f_rm, 3),
        "total": round(f_total, 3),
    }
    out["whatif_reference"] = {"low": low0, "high": high0}
    return out


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
    """Format the range as 'low – high €/kg' or None if unavailable.

    Bug E6 fix: consistent locale-aware number formatting. German uses
    point as thousands separator and comma as decimal; English uses
    comma as thousands and point as decimal.
    """
    lo = result.get("low_eur_per_kg")
    hi = result.get("high_eur_per_kg")
    if lo is None or hi is None:
        return None
    def _fmt_de(x: float) -> str:
        if x >= 1000:
            return f"{x:,.0f}".replace(",", ".")
        if x >= 10:
            return f"{x:.0f}"
        return f"{x:.1f}".replace(".", ",")
    def _fmt_en(x: float) -> str:
        if x >= 1000:
            return f"{x:,.0f}"
        if x >= 10:
            return f"{x:.0f}"
        return f"{x:.1f}"
    _fmt = _fmt_en if lang == "en" else _fmt_de
    return f"{_fmt(lo)} – {_fmt(hi)} €/kg"
