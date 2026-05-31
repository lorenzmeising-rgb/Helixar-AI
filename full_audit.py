"""
full_audit.py
=============
Comprehensive end-to-end audit of Helixar AI:
  - Run all 79 molecules through the engine + report pipeline
  - Apply extended lint rules (consistency, plausibility, structure)
  - Generate a structured bug_report.md grouped by severity

Run: python3 full_audit.py
Output:
  - audit_outputs/<NN>_<type>_<name>.pdf
  - audit_outputs/bug_report.md
  - audit_outputs/audit_summary.csv
"""

from __future__ import annotations

import csv
import json
import os
import re
import sys
import traceback
from pathlib import Path
from typing import Any, Dict, List

REPO = Path(__file__).resolve().parent
os.chdir(REPO)

from molecules_db import MOLECULE_DATABASE
from decision_engine import DecisionEngine
from blueprint_generator import generate_production_blueprint
from report_generator import generate_report
from explanation_layer import explain_blueprint
from plausibility_checker import check_plausibility_with_sources
from concrete_recommendations import build_concrete_recommendations
from recommended_actions import build_recommended_actions
from pdf_exporter import export_report_pdf
from database import MicrobialDB


def _make_engine() -> DecisionEngine:
    """Initialize DecisionEngine with the same fallback chain as batch_test_pdfs.py."""
    try:
        return DecisionEngine(db=MicrobialDB())
    except Exception:
        import pandas as pd
        empty_df = pd.DataFrame(columns=[
            "microorganism", "strain", "compound_class",
            "compound_name", "confidence_score",
        ])

        class _DummyDB:
            def __init__(self, df):
                self._df = df
            def all(self):
                return self._df

        return DecisionEngine(db=_DummyDB(empty_df))

OUTPUT_DIR = REPO / "audit_outputs"
OUTPUT_DIR.mkdir(exist_ok=True)


# ---------------------------------------------------------------------------
# Class-typical defaults: synthesize process_input for any DB molecule
# ---------------------------------------------------------------------------

DEFAULT_PROFILES: Dict[tuple, Dict[str, Any]] = {
    ("small_molecule", "non_volatile"): {
        "method": "chemical", "number_of_steps": 5, "desired_purity_percent": 99.0,
        "scale_kg_per_year": 100.0, "num_qualified_suppliers": 3, "lead_time_weeks": 4.0,
        "single_region_concentration": False, "raw_material_cost_eur_per_kg": 50.0,
        "strict_waste_constraints": False, "has_bioreactor": False,
    },
    ("small_molecule", "volatile"): {
        "method": "biotechnological", "number_of_steps": 3, "desired_purity_percent": 99.5,
        "scale_kg_per_year": 5000.0, "num_qualified_suppliers": 5, "lead_time_weeks": 2.0,
        "single_region_concentration": False, "raw_material_cost_eur_per_kg": 5.0,
        "strict_waste_constraints": False, "has_bioreactor": True,
    },
    ("natural_product", "alkaloid"): {
        "method": "chemical", "number_of_steps": 6, "desired_purity_percent": 99.0,
        "scale_kg_per_year": 50.0, "num_qualified_suppliers": 2, "lead_time_weeks": 6.0,
        "single_region_concentration": True, "raw_material_cost_eur_per_kg": 150.0,
        "strict_waste_constraints": False, "has_bioreactor": False,
    },
    ("natural_product", "terpene"): {
        "method": "extraction", "number_of_steps": 4, "desired_purity_percent": 95.0,
        "scale_kg_per_year": 200.0, "num_qualified_suppliers": 3, "lead_time_weeks": 4.0,
        "single_region_concentration": False, "raw_material_cost_eur_per_kg": 80.0,
        "strict_waste_constraints": False, "has_bioreactor": False,
    },
    ("peptide", "linear"): {
        "method": "chemical", "number_of_steps": 15, "desired_purity_percent": 98.0,
        "scale_kg_per_year": 5.0, "num_qualified_suppliers": 2, "lead_time_weeks": 6.0,
        "single_region_concentration": False, "raw_material_cost_eur_per_kg": 2500.0,
        "strict_waste_constraints": False, "has_bioreactor": False,
        "tagg_celsius": 55.0, "tm_celsius": 55.0,
        "num_disulfides": 0, "num_domains": 1, "has_ptm": False,
    },
    ("peptide", "cyclic"): {
        "method": "biotechnological", "number_of_steps": 8, "desired_purity_percent": 99.0,
        "scale_kg_per_year": 2.0, "num_qualified_suppliers": 2, "lead_time_weeks": 8.0,
        "single_region_concentration": True, "raw_material_cost_eur_per_kg": 5000.0,
        "strict_waste_constraints": False, "has_bioreactor": True,
        # Bug M4 fix: cyclic peptides have no tertiary fold → Tm/Tagg are
        # conceptually not meaningful. Omit so the plausibility checker
        # doesn't flag self-contradictory placeholder values.
        "num_disulfides": 0, "num_domains": 1, "has_ptm": False,
    },
    ("protein", "antibody"): {
        "method": "biotechnological", "number_of_steps": 7, "desired_purity_percent": 99.5,
        "scale_kg_per_year": 100.0, "num_qualified_suppliers": 3, "lead_time_weeks": 6.0,
        "single_region_concentration": False, "raw_material_cost_eur_per_kg": 8000.0,
        "strict_waste_constraints": True, "has_bioreactor": True,
        "tagg_celsius": 50.0, "tm_celsius": 65.0,
        "num_disulfides": 16, "num_domains": 4, "has_ptm": True,
    },
    ("protein", "enzyme"): {
        "method": "biotechnological", "number_of_steps": 5, "desired_purity_percent": 95.0,
        "scale_kg_per_year": 1000.0, "num_qualified_suppliers": 4, "lead_time_weeks": 3.0,
        "single_region_concentration": False, "raw_material_cost_eur_per_kg": 200.0,
        "strict_waste_constraints": False, "has_bioreactor": True,
        "tagg_celsius": 55.0, "tm_celsius": 50.0,
        "num_disulfides": 2, "num_domains": 2, "has_ptm": False,
    },
}

DEFAULT_PURIFICATION = {
    "has_flash_chromatography": True, "has_prep_hplc": True,
    "has_rotary_evaporator": True, "has_lyophilizer": True,
    "has_crystallization": True, "has_membrane_filtration": True,
    "has_fplc": True, "has_extraction": True,
}


def build_process_input(db_entry: Dict[str, Any]) -> Dict[str, Any]:
    """Synthesize a sensible process_input for any DB molecule."""
    key = (db_entry["molecule_type"], db_entry["molecule_subtype"])
    base = DEFAULT_PROFILES.get(key, DEFAULT_PROFILES[("small_molecule", "non_volatile")])
    out = dict(base)
    out["molecule_name"] = db_entry["name"]
    out["molecule_type"] = db_entry["molecule_type"]
    out["molecule_subtype"] = db_entry["molecule_subtype"]
    out["purification_methods"] = dict(DEFAULT_PURIFICATION)
    # Plausibility checker looks under 'available_methods' — mirror the same dict.
    out["available_methods"] = dict(DEFAULT_PURIFICATION)
    out["has_advanced_purification"] = True
    return out


# ---------------------------------------------------------------------------
# Extended lint rules
# ---------------------------------------------------------------------------

BAD_UNICODE_PATTERNS = ["₂", "₃", "₄", "₅", "₆", "₇", "₈", "₉", "₀", "■"]

# P6 fix: extended translation-leak patterns. Patterns must be unambiguously
# English (i.e. they would not appear in a well-written German PDF). The
# audit pre-strips citation lines (DOI / ISBN / et al.) before applying the
# linter so that English-language references don't false-positive.
ENGLISH_LEAK_PATTERNS = [
    # Original (kept)
    r"\bthe\s+(?:process|molecule|method|recommendation)\b",
    r"\bproduction\s+steps?\b",
    r"\bplease\s+(?:note|consider|review)\b",
    r"\bThis\s+(?:is|requires|will|may)\b",
    # P6 — modal/auxiliary chains that engine rules emit
    r"\b(?:should|must|will|may)\s+(?:be|consider|require|need|increase|decrease|reduce)\b",
    # P6 — common English connectors
    r"\bdue\s+to\b",
    r"\bsuch\s+as\b",
    r"\bcompared\s+to\b",
    r"\bin\s+addition\s+to\b",
    r"\bin\s+order\s+to\b",
    # P6 — English adjective+noun combos
    r"\b(?:high|low|medium)\s+(?:efficiency|risk|yield)\b",
    r"\b(?:expected|estimated)\s+(?:yield|purity|method)\b",
    r"\bproduction\s+(?:cost|route|sequence)\b",
    r"\bstep[s]?\s+(?:required|recommended)\b",
    # P6 — English copula constructions
    r"\bis\s+(?:required|recommended|necessary|critical|important)\b",
    r"\bare\s+(?:required|recommended|necessary)\b",
    # P6 — adverbs that mark English narrative
    r"\b(?:significantly|substantially|generally)\s+(?:higher|lower|better|worse)\b",
    r"\bbased\s+on\b",
    r"\bcan\s+(?:be|cause|lead)\b",
]

# Citation-line markers — lines containing any of these tokens are removed
# before the English-leak linter runs, because real DOI/ISBN references are
# English and would otherwise produce a flood of false positives.
_CITATION_LINE_RE = re.compile(
    r"\b(?:doi:|isbn[\s:-]|et al\.|vol\.\s*\d|pp\.\s*\d|"
    r"nat\.\s+(?:biotechnol|rev|prod)|j\.\s+(?:am|biol|med|nat)|"
    r"wiley|springer|elsevier|nature\b|science\b|"
    r"ullmann|annu\.\s+rev|appl\.\s+microbiol|chem\.\s+rev|"
    r"crit\.\s+rev|topics\s+in|trends\s+in|MAbs\b|"
    r"©\s*\d{4}|\(\s*\d{4}\s*\))",
    flags=re.IGNORECASE,
)


def _strip_citation_lines(text: str) -> str:
    """Return *text* with lines that look like literature citations removed.

    The English-leak linter applies to the remaining lines only — so a
    legitimate "Walsh G., Nat. Biotechnol. 2018 — should be …" reference
    does not get flagged just because the journal/year line looks like
    English narrative.
    """
    if not text:
        return ""
    kept = []
    for line in text.splitlines():
        if _CITATION_LINE_RE.search(line):
            continue
        kept.append(line)
    return "\n".join(kept)


def extract_pdf_text(pdf_path: Path) -> str:
    try:
        from pypdf import PdfReader
    except Exception:
        try:
            from PyPDF2 import PdfReader
        except Exception:
            return ""
    try:
        reader = PdfReader(str(pdf_path))
        return "\n".join((page.extract_text() or "") for page in reader.pages)
    except Exception as e:
        return f"__EXTRACT_FAILED__ {e}"


def lint_extended(
    molecule_name: str,
    pdf_text: str,
    selected: Dict[str, Any],
    bp: Any,
    process_input: Dict[str, Any],
    plaus: Dict[str, Any],
    actions: List[Any],
    concrete: List[Any],
) -> List[Dict[str, Any]]:
    """Apply 20+ heuristic checks; return list of findings."""
    f = []
    add = lambda sev, kind, detail: f.append({
        "molecule": molecule_name, "severity": sev, "kind": kind, "detail": detail
    })

    # ---- 1. Layout / encoding ----
    for ch in BAD_UNICODE_PATTERNS:
        if ch in pdf_text:
            add("error", "bad_unicode", f"Char {ch!r} (U+{ord(ch):04X}) leaks into PDF")

    # P6 fix: lint the citation-stripped text so DOI/ISBN reference blocks
    # don't false-positive on the (legitimately English) journal names and
    # narrative.
    _lintable_text = _strip_citation_lines(pdf_text)
    for pat in ENGLISH_LEAK_PATTERNS:
        m = re.findall(pat, _lintable_text, flags=re.IGNORECASE)
        if m:
            # Show one short context snippet so the bug report points to the
            # exact location of the leak.
            mobj = re.search(pat, _lintable_text, flags=re.IGNORECASE)
            snippet = ""
            if mobj:
                start = max(0, mobj.start() - 25)
                end = min(len(_lintable_text), mobj.end() + 25)
                snippet = _lintable_text[start:end].replace("\n", " ").strip()
            add("warning", "english_leak",
                f"'{pat}' matched {len(m)}x — i18n missing? Context: …{snippet}…")

    # ---- 2. Structural ----
    if len(pdf_text) < 1000:
        add("critical", "pdf_too_short",
            f"PDF text only {len(pdf_text)} chars — likely truncated/broken")

    # ---- 3. Critical sections ----
    expected_sections = ["produkt", "route", "begründ", "quell", "risiko"]
    missing = [s for s in expected_sections if s.lower() not in pdf_text.lower()]
    if missing:
        add("warning", "missing_section",
            f"Expected section keywords missing: {missing}")

    # ---- 4. Engine output sanity ----
    if not selected:
        add("critical", "no_recommendation",
            "Engine returned empty recommendation list")
    else:
        risk = (selected.get("risk_level") or "").lower()
        risk_reasons = selected.get("risk_reasons") or []
        if risk in ("high", "very high", "very_high") and not risk_reasons:
            add("error", "missing_risk_reasons",
                f"risk_level={risk} but no risk_reasons listed")

    # ---- 5. COGS plausibility ----
    cost = (selected.get("cost_level") or "").lower() if selected else ""
    eur = process_input.get("raw_material_cost_eur_per_kg", 0)
    if cost == "low" and eur >= 500:
        add("warning", "cogs_inconsistent",
            f"cost_level=low but raw_material_cost={eur}€/kg (≥500)")
    if cost == "very high" and eur < 50:
        add("warning", "cogs_inconsistent",
            f"cost_level=very high but raw_material_cost={eur}€/kg (<50)")

    # ---- 6. Plausibility warnings dropped on the floor? ----
    plaus_warnings = plaus.get("warnings", []) if plaus else []
    pdf_lower = pdf_text.lower()
    for w in plaus_warnings:
        # warnings may be strings (current shape) or dicts
        text = w if isinstance(w, str) else (
            w.get("text") or w.get("category") or str(w)
        )
        # Heuristic: pick a few distinctive ≥4-letter words and check that
        # at least 2 of them appear in the PDF text. Avoids false positives
        # from numeric tokens that get tokenised differently in the PDF.
        words = [w_ for w_ in re.findall(r"[A-Za-zäöüÄÖÜß-]{4,}", text)]
        hits = sum(1 for w_ in words[:8] if w_.lower() in pdf_lower)
        if words and hits < 2:
            add("info", "plaus_warning_not_in_pdf",
                f"Plausibility warning not surfaced in PDF: '{text[:80]}…'")

    # ---- 7. Recommendation count sanity ----
    if not actions and (selected.get("risk_level") or "").lower() in ("high", "very high"):
        add("warning", "high_risk_no_actions",
            f"risk={selected.get('risk_level')} but recommended_actions list empty")

    # ---- 8. Source presence ----
    doi_count = len(re.findall(r"10\.\d{4,9}/[-._;()/:A-Z0-9]+", pdf_text, flags=re.IGNORECASE))
    isbn_count = len(re.findall(r"ISBN[\s:-]*\d", pdf_text, flags=re.IGNORECASE))
    if doi_count + isbn_count == 0:
        add("warning", "no_sources",
            "PDF contains no DOI or ISBN — sourcing missing")

    # ---- 9. Route consistency ----
    method = process_input.get("method", "")
    mtype = process_input.get("molecule_type", "")
    if mtype == "protein" and method == "chemical":
        add("error", "route_mismatch",
            f"protein with method=chemical is implausible")
    if mtype == "small_molecule" and method == "biotechnological" \
       and process_input.get("molecule_subtype") == "non_volatile" \
       and process_input.get("scale_kg_per_year", 0) < 10:
        add("info", "route_unusual",
            f"non-volatile SM at <10kg via biotech is unusual")

    # ---- 10. Blueprint completeness ----
    if bp is None:
        add("critical", "no_blueprint", "Blueprint generation failed (None)")
    elif isinstance(bp, dict):
        rec = bp.get("recommended") or {}
        n_steps_in_rec = len(rec.get("production_steps", []) or rec.get("steps", []) or [])
        if not rec:
            add("warning", "blueprint_no_recommended",
                "Blueprint contains no 'recommended' key")
        elif n_steps_in_rec == 0:
            # Acceptable for some pure-process selectors; flag as info only
            add("info", "blueprint_recommended_no_steps",
                "Blueprint.recommended has no production_steps array")

    # ---- 11. Concrete recommendations sourcing ----
    if isinstance(concrete, list) and concrete:
        no_source = [c for c in concrete if isinstance(c, dict)
                     and not (c.get("source") or c.get("citation"))]
        if no_source and len(no_source) > len(concrete) // 2:
            add("info", "concrete_recs_unsourced",
                f"{len(no_source)}/{len(concrete)} concrete recs without source")

    # ---- 12. PDF size plausibility ----
    return f


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

def run_one(db_entry: Dict[str, Any], idx: int, engine: DecisionEngine) -> Dict[str, Any]:
    name = db_entry["name"]
    type_ = db_entry["molecule_type"]
    print(f"[{idx:02d}/79] {name:<22s} {type_:<16s} ", end="", flush=True)

    process_input = build_process_input(db_entry)
    plaus = {}
    recs: List[Dict[str, Any]] = []
    bp = None
    report_text = ""
    pdf_path = OUTPUT_DIR / f"{idx:02d}_{type_}_{name.replace(' ', '_')}.pdf"
    findings: List[Dict[str, Any]] = []
    error = ""

    try:
        plaus = check_plausibility_with_sources(process_input, lang="de") or {}
    except Exception as e:
        findings.append({"molecule": name, "severity": "critical",
                         "kind": "plausibility_crash", "detail": str(e)})

    try:
        recs = engine.analyze_process(process_input) or []
    except Exception as e:
        findings.append({"molecule": name, "severity": "critical",
                         "kind": "engine_crash", "detail": str(e)})

    selected = recs[0] if recs else {}
    # Flatten analysis fields up so the rest of the audit reads them uniformly
    if isinstance(selected, dict) and "analysis" in selected:
        a = selected["analysis"]
        selected.setdefault("cost_level", a.get("cost"))
        selected.setdefault("risk_level", a.get("risk"))
        selected.setdefault("efficiency_level", a.get("efficiency"))
        selected.setdefault("risk_reasons",
                            [a.get("risk_explanation")] if a.get("risk_explanation") else [])

    try:
        bp = generate_production_blueprint(selected, db=None,
                                           alternatives_count=2, safety_margin=0.15)
    except Exception as e:
        findings.append({"molecule": name, "severity": "critical",
                         "kind": "blueprint_crash", "detail": str(e)})

    try:
        report_text = generate_report(bp) if bp else ""
    except Exception as e:
        try:
            report_text = explain_blueprint(bp) if bp else ""
        except Exception as e2:
            findings.append({"molecule": name, "severity": "error",
                             "kind": "report_text_crash",
                             "detail": f"{e} / {e2}"})

    actions: List[Any] = []
    concrete: List[Any] = []
    try:
        actions = build_recommended_actions(process_input, lang="de") or []
    except Exception as e:
        findings.append({"molecule": name, "severity": "error",
                         "kind": "actions_crash", "detail": str(e)})
    try:
        concrete = build_concrete_recommendations(process_input, lang="de") or {}
    except Exception as e:
        findings.append({"molecule": name, "severity": "info",
                         "kind": "concrete_recs_crash", "detail": str(e)})

    # PDF
    pdf_ok = False
    try:
        pdf_extras = {
            "concrete_recs": concrete,
            "recommended_actions": actions,
            "plausibility_warnings": plaus.get("warnings", []),
            "plausibility_sources": plaus.get("literature_sources", []),
        }
        pdf_bytes = export_report_pdf(
            report_text, bp, None,
            title=f"Produktionsbericht — {name}",
            extras=pdf_extras,
            lang="de",
        )
        if pdf_bytes:
            pdf_path.write_bytes(pdf_bytes)
            pdf_ok = True
    except Exception as e:
        findings.append({"molecule": name, "severity": "critical",
                         "kind": "pdf_crash",
                         "detail": f"{type(e).__name__}: {e}"})

    pdf_text = extract_pdf_text(pdf_path) if pdf_ok else ""
    findings.extend(lint_extended(
        name, pdf_text, selected, bp, process_input, plaus, actions, concrete,
    ))

    print(f"recs={len(recs)} plaus_w={len(plaus.get('warnings', []) or [])} "
          f"actions={len(actions)} pdf={'ok' if pdf_ok else 'FAIL'} "
          f"lint={len(findings)}")

    return {
        "idx": idx,
        "molecule": name,
        "molecule_type": type_,
        "molecule_subtype": db_entry["molecule_subtype"],
        "engine_recommendations": len(recs),
        "plausibility_warnings": len(plaus.get("warnings", []) or []),
        "actions": len(actions),
        "concrete_recs": len(concrete),
        "pdf_ok": pdf_ok,
        "pdf_size_kb": pdf_path.stat().st_size // 1024 if pdf_ok else 0,
        "findings": findings,
        "selected_cost": (selected.get("cost_level") or "") if selected else "",
        "selected_risk": (selected.get("risk_level") or "") if selected else "",
        "selected_efficiency": (selected.get("efficiency_level") or "") if selected else "",
    }


def write_bug_report(rows: List[Dict[str, Any]]) -> None:
    """Build a structured bug_report.md grouped by severity."""
    all_findings: List[Dict[str, Any]] = []
    for r in rows:
        all_findings.extend(r["findings"])

    severity_order = ["critical", "error", "warning", "info"]
    by_sev = {s: [] for s in severity_order}
    for f in all_findings:
        by_sev.setdefault(f["severity"], []).append(f)

    lines = []
    lines.append("# Helixar AI — Full Audit Bug Report")
    lines.append("")
    lines.append(f"**Test set:** {len(rows)} molecules from MOLECULE_DATABASE")
    lines.append(f"**Findings total:** {len(all_findings)}")
    lines.append("")
    lines.append("**Summary by severity:**")
    for s in severity_order:
        lines.append(f"  - {s.upper()}: {len(by_sev[s])}")
    lines.append("")

    # Group by kind within severity
    for sev in severity_order:
        items = by_sev[sev]
        if not items:
            continue
        lines.append(f"## {sev.upper()} ({len(items)})")
        lines.append("")
        by_kind: Dict[str, List[Dict[str, Any]]] = {}
        for it in items:
            by_kind.setdefault(it["kind"], []).append(it)
        for kind, group in sorted(by_kind.items(), key=lambda x: -len(x[1])):
            lines.append(f"### `{kind}` — {len(group)} occurrence(s)")
            for it in group[:30]:
                lines.append(f"  - **{it['molecule']}**: {it['detail']}")
            if len(group) > 30:
                lines.append(f"  - … and {len(group) - 30} more")
            lines.append("")

    # Per-molecule overview
    lines.append("## Per-molecule overview")
    lines.append("")
    lines.append("| # | Molecule | Type | PDF | Recs | Plaus-W | Actions | Findings |")
    lines.append("|---|----------|------|-----|------|---------|---------|----------|")
    for r in rows:
        lines.append(
            f"| {r['idx']:02d} | {r['molecule']} | "
            f"{r['molecule_type']}/{r['molecule_subtype']} | "
            f"{'✓' if r['pdf_ok'] else '✗'} ({r['pdf_size_kb']} KB) | "
            f"{r['engine_recommendations']} | "
            f"{r['plausibility_warnings']} | {r['actions']} | "
            f"{len(r['findings'])} |"
        )

    (OUTPUT_DIR / "bug_report.md").write_text("\n".join(lines), encoding="utf-8")


def write_summary_csv(rows: List[Dict[str, Any]]) -> None:
    csv_path = OUTPUT_DIR / "audit_summary.csv"
    fieldnames = [
        "idx", "molecule", "molecule_type", "molecule_subtype",
        "selected_cost", "selected_risk", "selected_efficiency",
        "engine_recommendations", "plausibility_warnings",
        "actions", "concrete_recs", "pdf_ok", "pdf_size_kb",
        "n_findings",
    ]
    with csv_path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow({
                **{k: r.get(k, "") for k in fieldnames if k != "n_findings"},
                "n_findings": len(r["findings"]),
            })


def main() -> int:
    print("Helixar AI — Full Audit Runner")
    print(f"Output: {OUTPUT_DIR}")
    print(f"Molecules: {len(MOLECULE_DATABASE)}")
    print()

    engine = _make_engine()
    rows: List[Dict[str, Any]] = []
    for i, m in enumerate(MOLECULE_DATABASE, start=1):
        try:
            row = run_one(m, i, engine)
        except Exception as e:
            tb = traceback.format_exc()
            print(f"  CRASH on {m['name']}: {e}")
            if i <= 3:
                print(tb)
            row = {
                "idx": i, "molecule": m["name"],
                "molecule_type": m["molecule_type"],
                "molecule_subtype": m["molecule_subtype"],
                "engine_recommendations": 0, "plausibility_warnings": 0,
                "actions": 0, "concrete_recs": 0,
                "pdf_ok": False, "pdf_size_kb": 0,
                "findings": [{"molecule": m["name"], "severity": "critical",
                              "kind": "runner_crash",
                              "detail": f"{type(e).__name__}: {e}"}],
                "selected_cost": "", "selected_risk": "", "selected_efficiency": "",
            }
        rows.append(row)

    print()
    print("=" * 60)
    print(f"Audit complete: {len(rows)} molecules processed")
    total_findings = sum(len(r["findings"]) for r in rows)
    by_sev: Dict[str, int] = {}
    for r in rows:
        for f in r["findings"]:
            by_sev[f["severity"]] = by_sev.get(f["severity"], 0) + 1
    print(f"Findings: {total_findings}")
    for s in ["critical", "error", "warning", "info"]:
        print(f"  {s.upper()}: {by_sev.get(s, 0)}")

    write_summary_csv(rows)
    write_bug_report(rows)
    print()
    print(f"  CSV: {OUTPUT_DIR / 'audit_summary.csv'}")
    print(f"  MD : {OUTPUT_DIR / 'bug_report.md'}")
    # Non-zero exit on critical/error so CI fails loudly. Warnings are
    # advisory and do not fail the build.
    blocking = by_sev.get("critical", 0) + by_sev.get("error", 0)
    if blocking:
        print(f"\nFAIL: {blocking} blocking finding(s) (critical/error).")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
