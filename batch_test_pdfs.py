"""
batch_test_pdfs.py
==================
Generate PDF reports for a curated sample of 12 molecules covering all
8 type/subtype combinations in the molecule DB.

Output:
- test_outputs/<NN>_<type>_<name>.pdf            — one PDF per molecule
- test_outputs/test_summary.csv                  — flat per-molecule summary
- test_outputs/test_lint.csv                     — lint findings (auto-detected issues)
- test_outputs/test_engine_output.json           — full engine + plausibility dump for review

Usage:
    cd /Users/lorenzmeising/Helixar-AI
    python3 batch_test_pdfs.py

This is a developer tool, not part of the user-facing app. It is excluded
from the Streamlit deployment via .gitignore (or by simply not being
imported from app.py).
"""

import csv
import json
import os
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# Ensure we run from the repo root no matter where invoked from
REPO_ROOT = Path(__file__).resolve().parent
os.chdir(REPO_ROOT)

# ---------- Imports from the Helixar codebase ----------
from molecules_db import MOLECULE_DATABASE, get_entry_by_name
from decision_engine import DecisionEngine
from blueprint_generator import generate_production_blueprint
from report_generator import generate_report
from explanation_layer import explain_blueprint
from plausibility_checker import check_plausibility_with_sources
from concrete_recommendations import build_concrete_recommendations
from recommended_actions import build_recommended_actions
from pdf_exporter import export_report_pdf
from database import MicrobialDB

OUTPUT_DIR = REPO_ROOT / "test_outputs"
OUTPUT_DIR.mkdir(exist_ok=True)


# ---------- Label-mapping helpers (mirroring app.py) ----------

def purity_label_from_percent(pct: float) -> str:
    if pct >= 99.0:
        return ">99%"
    if pct >= 97.0:
        return "high"
    return "standard"


def scale_label_from_kg_per_year(kg: float) -> str:
    if kg >= 1000.0:
        return "industrial"
    if kg >= 1.0:
        return "pilot"
    return "lab"


def cost_label_from_eur_per_kg(eur: float) -> str:
    if eur > 500.0:
        return "high"
    if eur >= 10.0:
        return "medium"
    return "low"


def availability_label(num_sup: int, lead_w: float, single_region: bool) -> str:
    if num_sup <= 1 or lead_w >= 8.0:
        return "low"
    if num_sup >= 4 and lead_w <= 2.0 and not single_region:
        return "high"
    if single_region and num_sup <= 2:
        return "low"
    return "medium"


def aggregation_label_from_tagg(tagg_c: float) -> str:
    if tagg_c >= 65.0:
        return "low"
    if tagg_c >= 50.0:
        return "medium"
    return "high"


def folding_label_from_structure(disulfides: int, domains: int, ptm: bool) -> str:
    if domains >= 3 or disulfides >= 3 or (domains >= 2 and disulfides >= 2):
        return "high"
    if domains >= 1 and (disulfides >= 1 or ptm):
        return "medium"
    return "low"


def stability_label_from_tm(tm_c: float) -> str:
    if tm_c >= 70.0:
        return "high"
    if tm_c >= 50.0:
        return "medium"
    return "low"


# ---------- Curated sample: 12 molecules, all 8 type/subtype combos ----------
# Each entry is a dict with realistic, domain-appropriate inputs.
# Inputs reflect typical production parameters for that class:
#   - small_molecule/non_volatile: chemical synthesis, mid-large scale
#   - small_molecule/volatile:     fermentation, large scale
#   - natural_product/alkaloid:    extraction + bio, low-mid scale
#   - natural_product/terpene:     bio, mid scale
#   - peptide/linear:              SPPS or recombinant, low scale, very expensive
#   - peptide/cyclic:              fermentation (NRPS) or SPPS, low scale
#   - protein/antibody:            CHO cells, biotech, mid scale, premium price
#   - protein/enzyme:              microbial fermentation, large scale, cheap

SAMPLE: List[Dict[str, Any]] = [
    # ---- small_molecule / non_volatile ----
    {
        "name": "Aspirin",
        "method": "chemical",
        "number_of_steps": 1,
        "desired_purity_percent": 99.0,
        "scale_kg_per_year": 5000.0,
        "num_qualified_suppliers": 5,
        "lead_time_weeks": 2.0,
        "single_region_concentration": False,
        "raw_material_cost_eur_per_kg": 8.0,
        "strict_waste_constraints": True,
        "has_bioreactor": False,
        "purification_methods": {
            "has_flash_chromatography": True, "has_prep_hplc": False,
            "has_rotary_evaporator": True, "has_lyophilizer": False,
            "has_crystallization": True, "has_membrane_filtration": False,
            "has_fplc": False, "has_extraction": True,
        },
    },
    {
        "name": "Ibuprofen",
        "method": "chemical",
        "number_of_steps": 4,
        "desired_purity_percent": 99.5,
        "scale_kg_per_year": 2000.0,
        "num_qualified_suppliers": 3,
        "lead_time_weeks": 4.0,
        "single_region_concentration": False,
        "raw_material_cost_eur_per_kg": 25.0,
        "strict_waste_constraints": True,
        "has_bioreactor": False,
        "purification_methods": {
            "has_flash_chromatography": True, "has_prep_hplc": False,
            "has_rotary_evaporator": True, "has_lyophilizer": False,
            "has_crystallization": True, "has_membrane_filtration": False,
            "has_fplc": False, "has_extraction": True,
        },
    },
    {
        "name": "Vanillin",
        "method": "biotechnological",
        "number_of_steps": 3,
        "desired_purity_percent": 98.5,
        "scale_kg_per_year": 100.0,
        "num_qualified_suppliers": 3,
        "lead_time_weeks": 4.0,
        "single_region_concentration": False,
        "raw_material_cost_eur_per_kg": 50.0,
        "strict_waste_constraints": False,
        "has_bioreactor": True,
        "purification_methods": {
            "has_flash_chromatography": False, "has_prep_hplc": False,
            "has_rotary_evaporator": True, "has_lyophilizer": False,
            "has_crystallization": True, "has_membrane_filtration": False,
            "has_fplc": False, "has_extraction": True,
        },
    },
    # ---- small_molecule / volatile ----
    {
        "name": "Ethanol",
        "method": "biotechnological",
        "number_of_steps": 2,
        "desired_purity_percent": 99.5,
        "scale_kg_per_year": 100000.0,
        "num_qualified_suppliers": 8,
        "lead_time_weeks": 1.0,
        "single_region_concentration": False,
        "raw_material_cost_eur_per_kg": 1.0,
        "strict_waste_constraints": False,
        "has_bioreactor": True,
        "purification_methods": {
            "has_flash_chromatography": False, "has_prep_hplc": False,
            "has_rotary_evaporator": False, "has_lyophilizer": False,
            "has_crystallization": False, "has_membrane_filtration": False,
            "has_fplc": False, "has_extraction": False,
        },
    },
    # ---- natural_product / alkaloid ----
    {
        "name": "Caffeine",
        "method": "chemical",
        "number_of_steps": 3,
        "desired_purity_percent": 99.0,
        "scale_kg_per_year": 1000.0,
        "num_qualified_suppliers": 4,
        "lead_time_weeks": 3.0,
        "single_region_concentration": False,
        "raw_material_cost_eur_per_kg": 30.0,
        "strict_waste_constraints": False,
        "has_bioreactor": False,
        "purification_methods": {
            "has_flash_chromatography": False, "has_prep_hplc": False,
            "has_rotary_evaporator": True, "has_lyophilizer": False,
            "has_crystallization": True, "has_membrane_filtration": False,
            "has_fplc": False, "has_extraction": True,
        },
    },
    {
        "name": "Morphine",
        "method": "chemical",
        "number_of_steps": 5,
        "desired_purity_percent": 99.5,
        "scale_kg_per_year": 50.0,
        "num_qualified_suppliers": 2,
        "lead_time_weeks": 6.0,
        "single_region_concentration": True,
        "raw_material_cost_eur_per_kg": 200.0,
        "strict_waste_constraints": True,
        "has_bioreactor": False,
        "purification_methods": {
            "has_flash_chromatography": True, "has_prep_hplc": True,
            "has_rotary_evaporator": True, "has_lyophilizer": False,
            "has_crystallization": True, "has_membrane_filtration": False,
            "has_fplc": False, "has_extraction": True,
        },
    },
    # ---- natural_product / terpene ----
    {
        "name": "Astaxanthin",
        "method": "biotechnological",
        "number_of_steps": 4,
        "desired_purity_percent": 95.0,
        "scale_kg_per_year": 20.0,
        "num_qualified_suppliers": 2,
        "lead_time_weeks": 6.0,
        "single_region_concentration": False,
        "raw_material_cost_eur_per_kg": 800.0,
        "strict_waste_constraints": False,
        "has_bioreactor": True,
        "purification_methods": {
            "has_flash_chromatography": True, "has_prep_hplc": False,
            "has_rotary_evaporator": True, "has_lyophilizer": False,
            "has_crystallization": True, "has_membrane_filtration": False,
            "has_fplc": False, "has_extraction": True,
        },
    },
    # ---- peptide / linear ----
    {
        "name": "Glutathione",
        "method": "biotechnological",
        "number_of_steps": 3,
        "desired_purity_percent": 98.0,
        "scale_kg_per_year": 500.0,
        "num_qualified_suppliers": 4,
        "lead_time_weeks": 3.0,
        "single_region_concentration": False,
        "raw_material_cost_eur_per_kg": 150.0,
        "strict_waste_constraints": False,
        "has_bioreactor": True,
        "purification_methods": {
            "has_flash_chromatography": False, "has_prep_hplc": False,
            "has_rotary_evaporator": False, "has_lyophilizer": True,
            "has_crystallization": True, "has_membrane_filtration": True,
            "has_fplc": True, "has_extraction": False,
        },
        "tagg_celsius": 60.0,
        "tm_celsius": 55.0,
        "num_disulfides": 0,
        "num_domains": 1,
        "has_ptm": False,
    },
    {
        "name": "Liraglutide",
        "method": "chemical",  # SPPS is a chemical synthesis
        "number_of_steps": 30,
        "desired_purity_percent": 99.5,
        "scale_kg_per_year": 10.0,
        "num_qualified_suppliers": 2,
        "lead_time_weeks": 8.0,
        "single_region_concentration": False,
        "raw_material_cost_eur_per_kg": 3000.0,
        "strict_waste_constraints": True,
        "has_bioreactor": False,
        "purification_methods": {
            "has_flash_chromatography": False, "has_prep_hplc": True,
            "has_rotary_evaporator": False, "has_lyophilizer": True,
            "has_crystallization": False, "has_membrane_filtration": True,
            "has_fplc": True, "has_extraction": False,
        },
        "tagg_celsius": 55.0,
        "tm_celsius": 60.0,
        "num_disulfides": 0,
        "num_domains": 1,
        "has_ptm": True,  # lipidation on Lys26
    },
    # ---- peptide / cyclic ----
    {
        "name": "Cyclosporine",
        "method": "biotechnological",
        "number_of_steps": 5,
        "desired_purity_percent": 99.5,
        "scale_kg_per_year": 30.0,
        "num_qualified_suppliers": 2,
        "lead_time_weeks": 8.0,
        "single_region_concentration": True,
        "raw_material_cost_eur_per_kg": 2000.0,
        "strict_waste_constraints": True,
        "has_bioreactor": True,
        "purification_methods": {
            "has_flash_chromatography": True, "has_prep_hplc": True,
            "has_rotary_evaporator": True, "has_lyophilizer": True,
            "has_crystallization": True, "has_membrane_filtration": True,
            "has_fplc": False, "has_extraction": True,
        },
        "tagg_celsius": 70.0,
        "tm_celsius": 80.0,
        "num_disulfides": 0,
        "num_domains": 1,
        "has_ptm": False,
    },
    # ---- protein / antibody ----
    {
        "name": "Trastuzumab",
        "method": "biotechnological",
        "number_of_steps": 5,
        "desired_purity_percent": 99.5,
        "scale_kg_per_year": 200.0,
        "num_qualified_suppliers": 2,
        "lead_time_weeks": 8.0,
        "single_region_concentration": False,
        "raw_material_cost_eur_per_kg": 5000.0,
        "strict_waste_constraints": True,
        "has_bioreactor": True,
        "purification_methods": {
            "has_flash_chromatography": False, "has_prep_hplc": False,
            "has_rotary_evaporator": False, "has_lyophilizer": True,
            "has_crystallization": False, "has_membrane_filtration": True,
            "has_fplc": True, "has_extraction": False,
        },
        "tagg_celsius": 68.0,
        "tm_celsius": 75.0,
        "num_disulfides": 16,
        "num_domains": 12,
        "has_ptm": True,
    },
    # ---- protein / enzyme ----
    {
        "name": "Lipase",
        "method": "biotechnological",
        "number_of_steps": 3,
        "desired_purity_percent": 95.0,
        "scale_kg_per_year": 50000.0,
        "num_qualified_suppliers": 5,
        "lead_time_weeks": 3.0,
        "single_region_concentration": False,
        "raw_material_cost_eur_per_kg": 15.0,
        "strict_waste_constraints": False,
        "has_bioreactor": True,
        "purification_methods": {
            "has_flash_chromatography": False, "has_prep_hplc": False,
            "has_rotary_evaporator": False, "has_lyophilizer": False,
            "has_crystallization": False, "has_membrane_filtration": True,
            "has_fplc": False, "has_extraction": False,
        },
        "tagg_celsius": 55.0,
        "tm_celsius": 60.0,
        "num_disulfides": 2,
        "num_domains": 1,
        "has_ptm": False,
    },
]


def build_process_input(entry: Dict[str, Any]) -> Dict[str, Any]:
    """Convert a SAMPLE entry into the processInput dict the engine expects."""
    db_entry = get_entry_by_name(entry["name"])
    if not db_entry:
        raise ValueError(f"Molecule {entry['name']!r} not found in molecules_db")

    pi = {
        "molecule_name": entry["name"],
        "has_existing_process": True,
        "method": entry["method"],
        "number_of_steps": entry["number_of_steps"],
        "desired_purity": purity_label_from_percent(entry["desired_purity_percent"]),
        "desired_purity_percent": entry["desired_purity_percent"],
        "molecule_type": db_entry["molecule_type"],
        "molecule_subtype": db_entry["molecule_subtype"],
        "smiles": db_entry.get("smiles"),
        "db_entry": db_entry,
        "scale": scale_label_from_kg_per_year(entry["scale_kg_per_year"]),
        "scale_kg_per_year": entry["scale_kg_per_year"],
        "raw_material_availability": availability_label(
            entry["num_qualified_suppliers"],
            entry["lead_time_weeks"],
            entry["single_region_concentration"],
        ),
        "num_qualified_suppliers": entry["num_qualified_suppliers"],
        "lead_time_weeks": entry["lead_time_weeks"],
        "single_region_concentration": entry["single_region_concentration"],
        "raw_material_cost": cost_label_from_eur_per_kg(entry["raw_material_cost_eur_per_kg"]),
        "raw_material_cost_eur_per_kg": entry["raw_material_cost_eur_per_kg"],
        "strict_waste_constraints": entry["strict_waste_constraints"],
        "has_advanced_purification": any(entry["purification_methods"].values()),
        "has_bioreactor": entry["has_bioreactor"],
        "available_methods": entry["purification_methods"],
    }

    # Biophysical only for peptide/protein
    if db_entry["molecule_type"] in ("peptide", "protein"):
        tagg = entry.get("tagg_celsius")
        tm = entry.get("tm_celsius")
        ds = int(entry.get("num_disulfides") or 0)
        dom = int(entry.get("num_domains") or 1)
        ptm = bool(entry.get("has_ptm"))
        pi["aggregation_risk"] = aggregation_label_from_tagg(tagg) if tagg else "medium"
        pi["folding_complexity"] = folding_label_from_structure(ds, dom, ptm)
        pi["biophysical_stability"] = stability_label_from_tm(tm) if tm else "medium"
        pi["tagg_celsius"] = tagg
        pi["tm_celsius"] = tm
        pi["num_disulfides"] = ds
        pi["num_domains"] = dom
        pi["has_ptm"] = ptm

    return pi


# ---------- Lint patterns ----------

# Known English phrases that should appear translated in the German report.
# If they show up in the DE PDF text, that's a translation gap.
ENGLISH_LINT_PATTERNS = [
    r"\bHigh synthesis step count\b",
    r"\bDifficult purification\b",
    r"\bExpensive starting materials\b",
    r"\bLow biophysical stability\b",
    r"\bHigh aggregation risk\b",
    r"\bMultiple synthesis steps\b",
    r"\bSimplify process design\b",
    r"\bStrict waste handling\b",
    r"\bPurity targets above standard\b",
    r"\bChemical synthesis includes energy-intensive steps\b",
    r"\bHigh purity significantly increases downstream complexity\b",
    r"\b(Aromatic|crystallizable) structure\b",
    # Common engine-output strings that should be translated:
    r"\bRecommended Actions\b",
    r"\bConcrete Recommendations\b",
    r"\bIssues\b",  # Should be "Wichtige Themen" or similar in DE
]

# Unicode chars that the ReportLab default font cannot render.
# We strip them at PDF-generation time; if any leak through, that's a bug.
BAD_UNICODE_PATTERNS = [
    "■",   # ■ — replacement glyph for unsupported chars
    "₂",   # ₂
    "₃",   # ₃
    "₄",   # ₄
    "²",   # ²
    "³",   # ³
]


def lint_pdf_text(pdf_text: str, molecule_name: str, lang: str = "de") -> List[Dict[str, Any]]:
    """Scan extracted PDF text for known issues. Returns a list of findings."""
    findings = []

    if lang == "de":
        for pat in ENGLISH_LINT_PATTERNS:
            matches = re.findall(pat, pdf_text, flags=re.IGNORECASE)
            if matches:
                findings.append({
                    "molecule": molecule_name,
                    "severity": "warning",
                    "kind": "untranslated_english",
                    "detail": f"Pattern '{pat}' found {len(matches)}x",
                })

    for bad in BAD_UNICODE_PATTERNS:
        if bad in pdf_text:
            findings.append({
                "molecule": molecule_name,
                "severity": "error",
                "kind": "bad_unicode",
                "detail": f"Char {bad!r} (U+{ord(bad):04X}) found in PDF text",
            })

    # Structural checks
    if len(pdf_text) < 500:
        findings.append({
            "molecule": molecule_name,
            "severity": "error",
            "kind": "pdf_too_short",
            "detail": f"PDF text only {len(pdf_text)} chars — likely empty/broken",
        })

    return findings


def extract_pdf_text(pdf_path: Path) -> str:
    """Extract text from a PDF for linting. Uses PyPDF2 if available."""
    try:
        from pypdf import PdfReader  # pypdf is the modern fork
    except Exception:
        try:
            from PyPDF2 import PdfReader  # legacy
        except Exception:
            return ""  # no PDF parser available — skip text lint
    try:
        reader = PdfReader(str(pdf_path))
        return "\n".join((page.extract_text() or "") for page in reader.pages)
    except Exception as e:
        return f"__EXTRACT_FAILED__ {e}"


# ---------- Main batch runner ----------

def run_one(entry: Dict[str, Any], idx: int, engine: DecisionEngine, lang: str = "de") -> Dict[str, Any]:
    """Run engine + PDF generation for one molecule. Returns a summary row."""
    name = entry["name"]
    print(f"\n[{idx+1:02d}/12] {name} ...")

    process_input = build_process_input(entry)

    # 1. Plausibility
    plaus = check_plausibility_with_sources(process_input, lang=lang)

    # 2. Engine
    try:
        recs = engine.analyze_process(process_input)
    except Exception as e:
        print(f"    engine failed: {e}")
        recs = []

    selected = recs[0] if recs else {}

    # 3. Blueprint + report text
    bp = None
    report_text = ""
    try:
        bp = generate_production_blueprint(selected, db=None, alternatives_count=2, safety_margin=0.15)
        try:
            report_text = generate_report(bp)
        except Exception:
            report_text = explain_blueprint(bp) if bp else ""
    except Exception as e:
        print(f"    blueprint failed: {e}")

    # 4. Concrete recs + recommended actions
    try:
        crecs = build_concrete_recommendations(process_input, lang=lang)
    except Exception:
        crecs = None
    try:
        ractions = build_recommended_actions(process_input, lang=lang)
    except Exception:
        ractions = []

    # 5. PDF
    pdf_path = OUTPUT_DIR / f"{idx+1:02d}_{process_input['molecule_type']}_{name.replace(' ', '_').replace('α', 'a').replace('β', 'b')}.pdf"
    try:
        pdf_extras = {
            "concrete_recs": crecs,
            "recommended_actions": ractions,
            "plausibility_warnings": plaus.get("warnings", []),
            "plausibility_sources": plaus.get("literature_sources", []),
        }
        pdf_bytes = export_report_pdf(
            report_text, bp, None,
            title=f"Produktionsbericht — {name}",
            extras=pdf_extras,
            lang=lang,
        )
        if pdf_bytes:
            pdf_path.write_bytes(pdf_bytes)
            pdf_size = len(pdf_bytes)
        else:
            pdf_size = 0
            pdf_path = None
    except Exception as e:
        print(f"    PDF export failed: {e}")
        pdf_path = None
        pdf_size = 0

    # 6. Lint
    lint_findings = []
    pdf_text = ""
    if pdf_path:
        pdf_text = extract_pdf_text(pdf_path)
        lint_findings = lint_pdf_text(pdf_text, name, lang=lang)

    # 7. Build summary row.
    # Engine analysis fields live under selected["analysis"], not at top-level.
    analysis = (selected or {}).get("analysis") or {}
    props = analysis.get("properties") or {}
    row = {
        "idx": idx + 1,
        "molecule": name,
        "molecule_type": process_input["molecule_type"],
        "molecule_subtype": process_input["molecule_subtype"],
        "method": process_input["method"],
        "scale": process_input["scale"],
        "purity_target": process_input["desired_purity_percent"],
        "engine_cost": analysis.get("cost"),
        "engine_complexity": props.get("complexity"),
        "engine_risk": analysis.get("risk"),
        "engine_efficiency": analysis.get("efficiency"),
        "engine_toxicity": analysis.get("final_toxicity") or analysis.get("toxicity"),
        "engine_process_class": analysis.get("process_class"),
        "n_issues": len(analysis.get("issues") or []),
        "n_improvements": len(analysis.get("improvements") or []),
        "n_plausibility_warnings": len(plaus.get("warnings") or []),
        "n_recommended_actions": len(ractions or []),
        "n_concrete_recs": len(crecs.get("recommendations", []) if isinstance(crecs, dict) else (crecs or [])),
        "n_lint_findings": len(lint_findings),
        "pdf_size_bytes": pdf_size,
        "pdf_path": str(pdf_path.relative_to(REPO_ROOT)) if pdf_path else "",
    }

    print(f"    cost={row['engine_cost']} complexity={row['engine_complexity']} "
          f"risk={row['engine_risk']} eff={row['engine_efficiency']} "
          f"plaus={row['n_plausibility_warnings']} actions={row['n_recommended_actions']} "
          f"lint={row['n_lint_findings']} pdf={pdf_size//1024}KB")

    # Return both summary row and detailed engine output for the JSON dump
    return {
        "summary": row,
        "engine_output": selected,
        "plausibility": plaus,
        "recommended_actions": ractions,
        "concrete_recs": crecs,
        "lint_findings": lint_findings,
        "pdf_text_preview": pdf_text[:2000] if pdf_text else "",
    }


def main():
    print(f"Helixar AI — Batch test runner")
    print(f"Output dir: {OUTPUT_DIR}")
    print(f"Molecules:  {len(SAMPLE)}")

    engine = DecisionEngine(db=None, db_path=None) if False else None
    try:
        engine = DecisionEngine(db=MicrobialDB())
    except Exception:
        # Build with empty in-memory DB if MicrobialDB() needs args
        import pandas as pd
        empty_df = pd.DataFrame(columns=["microorganism", "strain", "compound_class", "compound_name", "confidence_score"])
        class _DummyDB:
            def __init__(self, df): self._df = df
            def all(self): return self._df
        engine = DecisionEngine(db=_DummyDB(empty_df))

    summary_rows = []
    all_lint = []
    full_outputs = {}

    for i, entry in enumerate(SAMPLE):
        result = run_one(entry, i, engine, lang="de")
        summary_rows.append(result["summary"])
        all_lint.extend(result["lint_findings"])
        full_outputs[entry["name"]] = {
            "engine_output": result["engine_output"],
            "plausibility": result["plausibility"],
            "recommended_actions": result["recommended_actions"],
            "concrete_recs": result["concrete_recs"],
            "lint_findings": result["lint_findings"],
            "pdf_text_preview": result["pdf_text_preview"],
        }

    # Write summary CSV
    summary_csv = OUTPUT_DIR / "test_summary.csv"
    with summary_csv.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(summary_rows[0].keys()))
        w.writeheader()
        w.writerows(summary_rows)

    # Write lint CSV
    lint_csv = OUTPUT_DIR / "test_lint.csv"
    with lint_csv.open("w", newline="", encoding="utf-8") as f:
        if all_lint:
            w = csv.DictWriter(f, fieldnames=list(all_lint[0].keys()))
            w.writeheader()
            w.writerows(all_lint)
        else:
            f.write("# No lint findings — all PDFs clean.\n")

    # Write full engine output JSON (for my own deep review later)
    json_path = OUTPUT_DIR / "test_engine_output.json"
    with json_path.open("w", encoding="utf-8") as f:
        # Pre-serialize: walk and convert non-JSON-able objects to strings
        def _serialize(obj):
            if isinstance(obj, dict):
                return {str(k): _serialize(v) for k, v in obj.items()}
            if isinstance(obj, (list, tuple)):
                return [_serialize(x) for x in obj]
            if isinstance(obj, (str, int, float, bool)) or obj is None:
                return obj
            return repr(obj)
        json.dump(_serialize(full_outputs), f, ensure_ascii=False, indent=2)

    # Print summary
    print(f"\n{'='*60}")
    print(f"DONE. Outputs in {OUTPUT_DIR.relative_to(REPO_ROOT)}/:")
    print(f"  - {summary_csv.name}: {len(summary_rows)} rows")
    print(f"  - {lint_csv.name}: {len(all_lint)} findings")
    print(f"  - {json_path.name}: full engine dump")
    print(f"  - {len(SAMPLE)} PDF files")

    if all_lint:
        print(f"\n⚠️  {len(all_lint)} lint findings:")
        for f in all_lint[:10]:
            print(f"   [{f['severity'].upper():>7}] {f['molecule']}: {f['detail']}")
        if len(all_lint) > 10:
            print(f"   ... + {len(all_lint)-10} more (see {lint_csv.name})")
    else:
        print("\n✅ Lint clean — no auto-detected issues.")


if __name__ == "__main__":
    main()
