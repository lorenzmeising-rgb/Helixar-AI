# ...existing code...
import os
from typing import Optional, List, Dict, Any

import streamlit as st
import pandas as pd

from database import MicrobialDB
from decision_engine import DecisionEngine
from blueprint_generator import generate_production_blueprint
from explanation_layer import explain_blueprint
from report_generator import generate_report
import feedback as feedback_module
import tempfile
import base64
from molecules_db import get_entry_by_name, MOLECULE_DATABASE, get_representation_label
from plausibility_checker import check_plausibility, check_plausibility_with_sources
from concrete_recommendations import build_concrete_recommendations
from recommended_actions import build_recommended_actions
from i18n import t, SUPPORTED_LANGUAGES, DEFAULT_LANGUAGE
from demo_scenarios import DEMOS, get_demo_by_id

# Locate a simple CSV DB in the project folder (optional). If not present, start empty.
BASE_DIR = os.path.dirname(__file__)
DB_CSV = os.path.join(BASE_DIR, "db.csv")


@st.cache_resource
def load_db(path: str) -> MicrobialDB:
    if os.path.exists(path):
        try:
            db = MicrobialDB.from_csv(path)
            try:
                # Debug: show loaded DB path and a few compounds
                print(f"[DEBUG] Loaded DB from: {path}")
                all_compounds = db.all().get("compound_name")
                if all_compounds is not None:
                    uniq = sorted([str(x) for x in pd.Series(all_compounds).dropna().unique()])
                    print(f"[DEBUG] Compounds in DB (sample): {uniq[:10]}")
            except Exception:
                pass
            return db
        except Exception:
            return MicrobialDB()
    return MicrobialDB()


# initialize DB and engine
db = load_db(DB_CSV)
engine = DecisionEngine(db=db)

# Resolve the logo path once — used as page icon, sidebar header, landing
# page, and PDF header. If the file is missing, all places fall back gracefully.
_LOGO_PATH = os.path.join(BASE_DIR, "helixar_logo.png")
_logo_for_tab = _LOGO_PATH if os.path.exists(_LOGO_PATH) else None

st.set_page_config(
    page_title="Helixar · Produktions-Blueprint",
    page_icon=_logo_for_tab,
    layout="centered",
)

# --- Global custom CSS — dark biotech theme matching the Helixar brand ---
st.markdown(
    """
    <style>
      /* ===== Brand palette ===== */
      :root {
        --hx-bg: #0D1B2A;
        --hx-bg-alt: #14253A;
        --hx-card: #1B2B40;
        --hx-card-hover: #243650;
        --hx-border: #2A3F5C;
        --hx-text: #E8EEF4;
        --hx-text-muted: #8AA0B8;
        --hx-primary: #1E3A5F;
        --hx-accent: #3DB8C5;
        --hx-accent-soft: rgba(61, 184, 197, 0.12);
      }
      /* ===== Typography ===== */
      h1, h2, h3 {
        color: var(--hx-text);
        font-weight: 600;
        letter-spacing: -0.01em;
      }
      h1 { font-size: 2.0rem; }
      h2 { font-size: 1.5rem; }
      h3 { font-size: 1.15rem; }
      p, li, .stMarkdown { color: var(--hx-text); }
      /* Section dividers — subtle */
      hr { border: none; height: 1px; background-color: var(--hx-border); margin: 1.5em 0; }
      /* ===== Sidebar ===== */
      [data-testid="stSidebar"] {
        background-color: var(--hx-bg-alt);
        border-right: 1px solid var(--hx-border);
      }
      [data-testid="stSidebar"] .stRadio > label { color: var(--hx-text); }
      [data-testid="stSidebar"] [role="radiogroup"] label {
        padding: 6px 8px;
        border-radius: 6px;
        margin-bottom: 2px;
        transition: background-color 0.15s ease;
      }
      [data-testid="stSidebar"] [role="radiogroup"] label:hover {
        background-color: var(--hx-accent-soft);
      }
      /* ===== Buttons ===== */
      .stButton > button {
        border-radius: 8px;
        font-weight: 500;
        border: 1px solid var(--hx-border);
        background-color: var(--hx-card);
        color: var(--hx-text);
        transition: all 0.15s ease;
      }
      .stButton > button:hover {
        background-color: var(--hx-accent);
        border-color: var(--hx-accent);
        color: #0D1B2A;
      }
      .stButton > button[kind="primary"] {
        background-color: var(--hx-accent);
        border-color: var(--hx-accent);
        color: #0D1B2A;
        font-weight: 600;
      }
      .stButton > button[kind="primary"]:hover {
        background-color: #5BCAD7;
      }
      /* ===== Metric cards ===== */
      [data-testid="stMetricValue"] {
        color: var(--hx-accent);
        font-weight: 600;
      }
      [data-testid="stMetricLabel"] {
        color: var(--hx-text-muted);
        font-size: 0.85rem;
      }
      /* ===== Inputs ===== */
      .stTextInput > div > div > input,
      .stNumberInput > div > div > input,
      .stSelectbox > div > div,
      .stTextArea textarea {
        background-color: var(--hx-card);
        border: 1px solid var(--hx-border);
        color: var(--hx-text);
        border-radius: 6px;
      }
      /* ===== Code blocks ===== */
      code {
        color: var(--hx-accent);
        background-color: var(--hx-bg-alt);
        padding: 1px 6px;
        border-radius: 4px;
        font-size: 0.9em;
      }
      /* ===== Captions ===== */
      [data-testid="stCaptionContainer"] { color: var(--hx-text-muted); }
      .stMarkdown small { color: var(--hx-text-muted); }
      /* ===== Cards (custom) ===== */
      .hx-card {
        background-color: var(--hx-card);
        border: 1px solid var(--hx-border);
        border-radius: 12px;
        padding: 22px 20px;
        margin-bottom: 12px;
        transition: border-color 0.15s ease;
      }
      .hx-card:hover { border-color: var(--hx-accent); }
      .hx-card-icon {
        font-size: 28px;
        margin-bottom: 10px;
        color: var(--hx-accent);
      }
      .hx-card-title {
        font-size: 1.05rem;
        font-weight: 600;
        color: var(--hx-text);
        margin-bottom: 6px;
      }
      .hx-card-body {
        font-size: 0.9rem;
        color: var(--hx-text-muted);
        line-height: 1.45;
      }
      /* ===== Hero section ===== */
      .hx-hero-title {
        font-size: 2.4rem;
        font-weight: 700;
        color: var(--hx-text);
        line-height: 1.15;
        margin-bottom: 14px;
      }
      .hx-hero-subtitle {
        font-size: 1.05rem;
        color: var(--hx-text-muted);
        line-height: 1.55;
        max-width: 480px;
      }
      /* ===== Step indicator ===== */
      .hx-steps {
        display: flex;
        align-items: center;
        gap: 0;
        margin: 8px 0 28px 0;
      }
      .hx-step {
        display: flex;
        flex-direction: column;
        align-items: center;
        flex: 0 0 auto;
        min-width: 110px;
      }
      .hx-step-circle {
        width: 32px;
        height: 32px;
        border-radius: 50%;
        background-color: var(--hx-card);
        border: 2px solid var(--hx-border);
        color: var(--hx-text-muted);
        display: flex;
        align-items: center;
        justify-content: center;
        font-weight: 600;
        font-size: 0.9rem;
      }
      .hx-step.active .hx-step-circle {
        background-color: var(--hx-accent);
        border-color: var(--hx-accent);
        color: #0D1B2A;
      }
      .hx-step.done .hx-step-circle {
        background-color: var(--hx-primary);
        border-color: var(--hx-primary);
        color: var(--hx-text);
      }
      .hx-step-label {
        margin-top: 6px;
        font-size: 0.85rem;
        color: var(--hx-text-muted);
      }
      .hx-step.active .hx-step-label, .hx-step.done .hx-step-label {
        color: var(--hx-text);
      }
      .hx-step-line {
        flex: 1 1 auto;
        height: 2px;
        background-color: var(--hx-border);
        margin: 0 -8px;
        margin-top: -22px;
        z-index: -1;
      }
      .hx-step-line.done { background-color: var(--hx-primary); }
      /* ===== Demo card grid ===== */
      .hx-demo-card {
        background-color: var(--hx-card);
        border: 1px solid var(--hx-border);
        border-radius: 12px;
        padding: 18px;
        height: 100%;
        transition: all 0.15s ease;
      }
      .hx-demo-card:hover {
        border-color: var(--hx-accent);
        background-color: var(--hx-card-hover);
      }
      .hx-demo-tag {
        display: inline-block;
        background-color: var(--hx-accent-soft);
        color: var(--hx-accent);
        padding: 2px 10px;
        border-radius: 12px;
        font-size: 0.75rem;
        font-weight: 500;
        margin-bottom: 10px;
      }
      /* Hide Streamlit's default header & footer for cleaner look */
      header[data-testid="stHeader"] { background-color: transparent; }
      footer { visibility: hidden; }
      /* Tighten main block padding */
      .block-container { padding-top: 2rem; padding-bottom: 3rem; }
    </style>
    """,
    unsafe_allow_html=True,
)


# --- Session state initialization --------------------------------------------
if "ui_lang" not in st.session_state:
    st.session_state.ui_lang = DEFAULT_LANGUAGE
if "page" not in st.session_state:
    # Landing page on app load; separate from the recommendation form/page
    st.session_state.page = "landing"

# --- Language switch is now on the Settings page; no global top-bar switch.
def _render_lang_switch():
    """Kept for backward compatibility. Top-bar language switch removed —
    language is now configured on the Settings page (Einstellungen)."""
    return
if "recommendations" not in st.session_state:
    st.session_state.recommendations = []
if "selected_strategy" not in st.session_state:
    st.session_state.selected_strategy = None
if "blueprint" not in st.session_state:
    st.session_state.blueprint = None
if "german_report" not in st.session_state:
    st.session_state.german_report = None
if "input_context" not in st.session_state:
    st.session_state.input_context = {}


def show_start_page():
    """Render the start page where users enter their target and request recommendations.

    This function only writes to session_state when the user explicitly clicks the
    "Empfehlungen generieren" button.
    """
    # Step indicator: Step 1 (Inputs) is active here
    render_step_indicator(active_step=1)

    st.title(t("app_title"))
    st.caption(t("app_caption"))

    # Single unified flow: the recommendation report. The former
    # "Prozessvergleich" mode is no longer a separate output — the route
    # comparison can now be appended directly to the recommendation PDF via
    # a checkbox next to the export buttons (avoids a redundant cover page
    # and a thin standalone document).
    # ----- Recommendation flow -----
    st.subheader(t("describe_process_header"))
    st.markdown(t("describe_process_intro"))
    # Show autofill notice (if autofill recently occurred)
    if st.session_state.get("proc_autofill_msg"):
        st.info(st.session_state.proc_autofill_msg)
    st.markdown("---")

    # 1. Target molecule
    # Prepare autofill control flags in session_state
    if "proc_autofill_locked" not in st.session_state:
        st.session_state.proc_autofill_locked = False
    if "proc_last_autofill_name" not in st.session_state:
        st.session_state.proc_last_autofill_name = ""
    if "proc_smiles" not in st.session_state:
        st.session_state.proc_smiles = None
    if "proc_db_entry" not in st.session_state:
        st.session_state.proc_db_entry = None
    if "proc_autofill_msg" not in st.session_state:
        st.session_state.proc_autofill_msg = ""

    # Feedback-Runde-3: realistische Industrie-Defaults pro (type, subtype).
    # Werte angelehnt an die Anchors in cogs_estimator.py + publizierte
    # Pharma-/Bulk-Benchmarks. Form: (purity_%, scale_kg_per_year,
    # raw_cost_eur_per_kg, num_steps).
    TYPE_DEFAULTS = {
        ("small_molecule", "non_volatile"): (99.0, 1000.0, 20.0, 3),     # API-Bulk (Aspirin-Klasse)
        ("small_molecule", "volatile"):     (99.5, 50000.0, 1.0, 2),     # Ethanol-Commodity
        ("natural_product", "alkaloid"):    (99.0, 100.0, 100.0, 4),     # Morphine/Codeine-Klasse
        ("natural_product", "terpene"):     (95.0, 100.0, 80.0, 4),      # Carotinoide
        ("peptide", "linear"):              (99.5, 10.0, 1000.0, 5),     # GLP-1-Analoga, Insulin
        ("peptide", "cyclic"):              (99.5, 5.0, 2000.0, 5),      # Cyclosporin-Klasse
        ("protein", "antibody"):            (99.5, 100.0, 5000.0, 5),    # mAbs
        ("protein", "enzyme"):              (95.0, 10000.0, 15.0, 3),    # Industrieenzyme
    }

    def _apply_type_defaults():
        """Setzt realistische Industrie-Defaults für andere Felder, sobald
        der Type/Subtype manuell geändert wurde. Wir setzen NUR Felder,
        die der User noch nicht selbst angefasst hat (Streamlit-State-
        Heuristik: Feld noch im 'initial value'-Zustand).
        Wichtig: Callback läuft vor dem Rerun, darf session_state freien
        Widget-Keys zuweisen."""
        mtype = st.session_state.get("proc_molecule_type_v2")
        msub = st.session_state.get("proc_molecule_subtype_v2")
        defaults = TYPE_DEFAULTS.get((mtype, msub)) or TYPE_DEFAULTS.get((mtype, "non_volatile"))
        if not defaults:
            return
        purity, scale, raw_cost, steps = defaults
        # Immer überschreiben — der User hat ja gerade aktiv Type/Subtype
        # gewechselt, also will er typischerweise neue passende Defaults sehen.
        st.session_state["proc_purity_percent"] = float(purity)
        st.session_state["proc_scale_kg_per_year"] = float(scale)
        st.session_state["proc_raw_cost_eur_per_kg"] = float(raw_cost)
        st.session_state["proc_number_of_steps"] = int(steps)

    def _lock_autofill():
        st.session_state.proc_autofill_locked = True
        # User has manually overridden molecule type/subtype.
        # The autofilled SMILES + DB entry no longer matches the chosen type,
        # so discard them to prevent mixing (e.g. Vanillin-SMILES + protein logic).
        if st.session_state.get("proc_smiles") or st.session_state.get("proc_db_entry"):
            st.session_state.proc_smiles = None
            st.session_state.proc_db_entry = None
            st.session_state.proc_autofill_msg = (
                "Manueller Override aktiv: erkannte Molekül-Repräsentation wurde "
                "verworfen, da der Molekültyp manuell geändert wurde."
            )
        # Auch Industrie-Defaults für die anderen Felder anwenden.
        _apply_type_defaults()

    def _autofill_from_entry(entry, name):
        """Apply a recognised DB entry's type/subtype/SMILES/representation to
        the form fields. Shared by the dropdown selection and the custom
        free-text path."""
        try:
            mt = entry.get("molecule_type") or "small_molecule"
            # Internal type code is stored directly (the selectbox uses
            # internal codes as option values).
            st.session_state.proc_molecule_type_v2 = mt
            subtype = entry.get("molecule_subtype")
            st.session_state.proc_molecule_subtype_v2 = subtype
            st.session_state.proc_smiles = entry.get("smiles")
            # Carry the full entry through so downstream stages (PDF, report)
            # can show the proper representation (sequence / external_id) when
            # SMILES is not applicable (proteins, complex peptides).
            st.session_state.proc_db_entry = entry
            st.session_state.proc_last_autofill_name = name
            rep = get_representation_label(entry)
            st.session_state.proc_autofill_msg = (
                f"Molekül erkannt: {mt} / {subtype or 'None'} — {rep}"
            )
        except Exception:
            pass

    def _on_molecule_select():
        """Dropdown selection of a known DB molecule → autofill its metadata.
        Streamlit runs on_change callbacks at the start of the rerun (before
        any widgets are instantiated), so writing the type/subtype session
        keys here is safe."""
        sel = st.session_state.get("proc_molecule_select")
        if not sel or sel == t("custom_molecule_option"):
            return
        try:
            entry = get_entry_by_name(sel)
        except Exception:
            entry = None
        if entry:
            # Picking from the dropdown is an explicit choice → unlock so the
            # new molecule's metadata replaces any prior manual override.
            st.session_state.proc_autofill_locked = False
            _autofill_from_entry(entry, sel)

    def _try_autofill():
        """Free-text custom molecule → best-effort match against the DB."""
        name = str(st.session_state.get("proc_molecule_name") or "").strip()
        if not name:
            return
        # If name changed, allow autofill again
        if name != st.session_state.proc_last_autofill_name:
            st.session_state.proc_autofill_locked = False
        if st.session_state.proc_autofill_locked:
            return
        entry = None
        try:
            entry = get_entry_by_name(name)
        except Exception:
            entry = None
        # partial match fallback
        if not entry:
            try:
                low = name.lower()
                for e in MOLECULE_DATABASE:
                    if low in str(e.get("name", "")).lower():
                        entry = e
                        break
            except Exception:
                entry = None
        if entry:
            _autofill_from_entry(entry, name)

    # Target molecule picker: dropdown of all DB molecules + a free-text
    # "custom" option. The dropdown makes the 79 calibrated molecules one
    # click away (like the old route-comparison picker); the custom option
    # preserves free entry of molecules not yet in the DB (which the P4
    # cluster-fallback path is built to handle).
    _db_names = sorted([m["name"] for m in MOLECULE_DATABASE])
    _CUSTOM_OPT = t("custom_molecule_option")
    _mol_options = [_CUSTOM_OPT] + _db_names
    _sel_mol = st.selectbox(
        t("target_molecule"),
        options=_mol_options,
        index=(_mol_options.index("Vanillin") if "Vanillin" in _mol_options else 0),
        key="proc_molecule_select",
        on_change=_on_molecule_select,
        help=t("target_molecule_help"),
    )
    if _sel_mol == _CUSTOM_OPT:
        molecule_name = st.text_input(
            t("target_molecule_custom_label"),
            value="",
            key="proc_molecule_name",
            on_change=_try_autofill,
        )
    else:
        molecule_name = _sel_mol
        # on_change does not fire on the initial render, so autofill the
        # currently-selected molecule's metadata here too — unless the user
        # has manually overridden type/subtype (autofill locked).
        if (st.session_state.get("proc_last_autofill_name") != _sel_mol
                and not st.session_state.get("proc_autofill_locked")):
            try:
                _e0 = get_entry_by_name(_sel_mol)
            except Exception:
                _e0 = None
            if _e0:
                _autofill_from_entry(_e0, _sel_mol)

    # Molecule type / subtype directly under target molecule.
    # Auto-filled when an entry is recognised; manual change locks the autofill
    # and clears the previously detected SMILES (see _lock_autofill).
    # Internal type key is in English; the user-facing label is translated.
    type_keys = ["small_molecule", "peptide", "protein", "natural_product"]
    type_label_keys = {
        "small_molecule": "type_small_molecule",
        "peptide": "type_peptide",
        "protein": "type_protein",
        "natural_product": "type_natural_product",
    }
    molecule_type = st.selectbox(
        t("molecule_type"),
        options=type_keys,
        format_func=lambda k: t(type_label_keys[k]),
        index=0,
        key="proc_molecule_type_v2",
        on_change=_lock_autofill,
        help=t("molecule_type_help"),
    )

    subtype_options_map = {
        "small_molecule": [(None, "subtype_none"), ("volatile", "subtype_volatile"), ("non_volatile", "subtype_non_volatile")],
        "peptide": [(None, "subtype_none"), ("linear", "subtype_linear"), ("cyclic", "subtype_cyclic")],
        "protein": [(None, "subtype_none"), ("antibody", "subtype_antibody"), ("enzyme", "subtype_enzyme")],
        "natural_product": [(None, "subtype_none"), ("terpene", "subtype_terpene"), ("alkaloid", "subtype_alkaloid")],
    }
    opts = subtype_options_map.get(molecule_type, [(None, "subtype_none")])
    opt_keys = [o[0] for o in opts]
    label_keys = {o[0]: o[1] for o in opts}
    molecule_subtype = st.selectbox(
        t("molecule_subtype"),
        options=opt_keys,
        format_func=lambda v: t(label_keys.get(v, "subtype_none")),
        index=0,
        key="proc_molecule_subtype_v2",
        on_change=_lock_autofill,
        help=t("molecule_subtype_help"),
    )

    # ----- Mode switch: existing process vs. greenfield -----
    st.markdown("---")
    no_existing_process = st.checkbox(
        t("no_existing_process_label"),
        value=False,
        key="proc_no_existing",
        help=t("no_existing_process_help"),
    )
    has_existing_process = not no_existing_process

    # ----- SECTION 1: Aktueller Prozess (only if user has one) -----
    method = None
    number_of_steps = None
    if has_existing_process:
        st.subheader(t("section_1_title"))
        method_keys = ["chemical", "biotechnological", "extraction"]
        method_label_keys = {
            "chemical": "method_chemical",
            "biotechnological": "method_biotechnological",
            "extraction": "method_extraction",
        }
        method = st.selectbox(
            t("current_method_label"),
            options=method_keys,
            format_func=lambda k: t(method_label_keys[k]),
            index=0,
            key="proc_method_v2",
            help=t("current_method_help"),
        )

        number_of_steps = st.number_input(
            t("number_of_steps_label"),
            min_value=1,
            max_value=10,
            value=3,
            step=1,
            key="proc_number_of_steps",
            help=t("number_of_steps_help"),
        )
    else:
        st.info(t("greenfield_info"))

    # ----- SECTION 2: Anforderungen an den Prozess -----
    st.markdown("---")
    st.subheader(t("section_2_title"))
    desired_purity_percent = st.number_input(
        t("purity_label"),
        min_value=80.0,
        max_value=99.999,
        value=95.0,
        step=0.1,
        format="%.3f",
        key="proc_purity_percent",
        help=t("purity_help"),
    )

    scale_kg_per_year = st.number_input(
        t("scale_label"),
        min_value=0.001,
        max_value=1_000_000.0,
        value=1.0,
        step=0.1,
        format="%.3f",
        key="proc_scale_kg_per_year",
        help=t("scale_help"),
    )

    # ----- SECTION 3: Marktbedingungen (Ist-Zustand) -----
    # Feedback-Runde-3: im Greenfield-Modus kann der User die Marktbedingungen
    # NICHT kennen (Lieferanten + Edukt-Preis hängen von der nicht-existenten
    # Route ab). Dann blenden wir die Section komplett aus und nutzen
    # Industrie-Benchmark-Defaults basierend auf dem Molekül-Typ.
    if has_existing_process:
        st.markdown("---")
        st.subheader(t("section_3_title"))
        st.caption(t("section_3_caption"))
        num_qualified_suppliers = st.number_input(
            t("num_suppliers_label"),
            min_value=1,
            max_value=20,
            value=3,
            step=1,
            key="proc_num_suppliers",
            help=t("num_suppliers_help"),
        )
        lead_time_weeks = st.number_input(
            t("lead_time_label"),
            min_value=0.5,
            max_value=52.0,
            value=4.0,
            step=0.5,
            format="%.1f",
            key="proc_lead_time_weeks",
            help=t("lead_time_help"),
        )
        single_region_concentration = st.checkbox(
            t("single_region_label"),
            value=False,
            key="proc_single_region",
            help=t("single_region_help"),
        )

        raw_material_cost_eur_per_kg = st.number_input(
            t("raw_cost_label"),
            min_value=0.1,
            max_value=100_000.0,
            value=50.0,
            step=1.0,
            format="%.2f",
            key="proc_raw_cost_eur_per_kg",
            help=t("raw_cost_help"),
        )
        strict_waste_constraints = st.checkbox(
            t("strict_waste_label"),
            value=False,
            key="proc_strict_waste",
            help=t("strict_waste_help"),
        )
    else:
        # Greenfield: realistische Industrie-Defaults pro Type/Subtype.
        # Anchors basierend auf publizierten Pharma-Benchmarks (vgl. cogs_estimator.py).
        st.markdown("---")
        st.info(t("greenfield_market_info"))
        _market_defaults = {
            ("small_molecule", "non_volatile"): (3, 4.0, False, 20.0, True),   # API-Bulk
            ("small_molecule", "volatile"):     (5, 2.0, False, 1.0,  False),  # Solvent-Commodity
            ("natural_product", "alkaloid"):    (2, 6.0, True,  100.0, True),
            ("natural_product", "terpene"):     (3, 4.0, False, 80.0, False),
            ("peptide", "linear"):              (2, 6.0, False, 1000.0, True),
            ("peptide", "cyclic"):              (2, 8.0, False, 2000.0, True),
            ("protein", "antibody"):            (2, 8.0, False, 5000.0, True),
            ("protein", "enzyme"):              (3, 4.0, False, 15.0, False),
        }
        _key = (molecule_type, molecule_subtype)
        _defaults = _market_defaults.get(_key, (3, 4.0, False, 50.0, True))
        (num_qualified_suppliers,
         lead_time_weeks,
         single_region_concentration,
         raw_material_cost_eur_per_kg,
         strict_waste_constraints) = _defaults

    # ----- SECTION 4: Verfügbare Methodiken & Equipment im Labor -----
    st.markdown("---")
    st.subheader(t("section_4_title"))
    st.caption(t("section_4_caption"))
    has_bioreactor = st.checkbox(
        t("bioreactor_label"),
        value=False,
        key="proc_has_bioreactor",
        help=t("bioreactor_help"),
    )
    st.markdown(t("purification_methods_header"))
    purification_methods = {
        "has_flash_chromatography": st.checkbox(
            t("method_flash_chrom_label"),
            value=False,
            key="proc_has_flash_chromatography",
            help=t("method_flash_chrom_help"),
        ),
        "has_prep_hplc": st.checkbox(
            t("method_prep_hplc_label"),
            value=False,
            key="proc_has_prep_hplc",
            help=t("method_prep_hplc_help"),
        ),
        "has_rotary_evaporator": st.checkbox(
            t("method_rotavap_label"),
            value=False,
            key="proc_has_rotary_evaporator",
            help=t("method_rotavap_help"),
        ),
        "has_lyophilizer": st.checkbox(
            t("method_lyo_label"),
            value=False,
            key="proc_has_lyophilizer",
            help=t("method_lyo_help"),
        ),
        "has_crystallization": st.checkbox(
            t("method_cryst_label"),
            value=False,
            key="proc_has_crystallization",
            help=t("method_cryst_help"),
        ),
        "has_membrane_filtration": st.checkbox(
            t("method_membrane_label"),
            value=False,
            key="proc_has_membrane_filtration",
            help=t("method_membrane_help"),
        ),
        "has_fplc": st.checkbox(
            t("method_fplc_label"),
            value=False,
            key="proc_has_fplc",
            help=t("method_fplc_help"),
        ),
        "has_extraction": st.checkbox(
            t("method_extraction_label"),
            value=False,
            key="proc_has_extraction",
            help=t("method_extraction_help"),
        ),
    }
    # Backwards-compatible flag for the engine: at least one "advanced" method ticked.
    has_advanced_purification = any(
        purification_methods[k]
        for k in ("has_prep_hplc", "has_lyophilizer", "has_fplc", "has_membrane_filtration")
    )

    # Biomolecule-only biophysical inputs (show only for peptide/protein)
    aggregation_risk = None
    folding_complexity = None
    biophysical_stability = None
    # Numeric biophysical inputs — replace the previous low/medium/high dropdowns.
    # Engine vocabulary (low/medium/high) is still used downstream; we map the
    # numeric inputs to those buckets via helper functions in the analyze handler.
    tagg_celsius = None
    tm_celsius = None
    num_disulfides = None
    num_domains = None
    has_ptm = None
    if molecule_type in ("peptide", "protein"):
        st.markdown("---")
        st.subheader(t("section_5_title"))
        tagg_celsius = st.number_input(
            t("tagg_label"),
            min_value=20.0,
            max_value=100.0,
            value=60.0,
            step=0.5,
            format="%.1f",
            key="proc_tagg_celsius",
            help=t("tagg_help"),
        )
        col_d1, col_d2 = st.columns(2)
        with col_d1:
            num_disulfides = st.number_input(
                t("disulfides_label"),
                min_value=0,
                max_value=30,
                value=1,
                step=1,
                key="proc_num_disulfides",
                help=t("disulfides_help"),
            )
        with col_d2:
            num_domains = st.number_input(
                t("domains_label"),
                min_value=1,
                max_value=20,
                value=1,
                step=1,
                key="proc_num_domains",
                help=t("domains_help"),
            )
        has_ptm = st.checkbox(
            t("ptm_label"),
            value=False,
            key="proc_has_ptm",
            help=t("ptm_help"),
        )
        tm_celsius = st.number_input(
            t("tm_label"),
            min_value=20.0,
            max_value=110.0,
            value=65.0,
            step=0.5,
            format="%.1f",
            key="proc_tm_celsius",
            help=t("tm_help"),
        )

    st.markdown("---")
    # Action button: analyze process
    if st.button(t("analyze_button")):
        # Validation
        if not molecule_name or str(molecule_name).strip() == "":
            st.warning(t("warn_target_required"))
            return
        # Step count is only required when the user actually has an existing process.
        if has_existing_process and (not number_of_steps or int(number_of_steps) < 1):
            st.warning(t("warn_steps_required"))
            return

        # Validate biologically/chemically impossible combinations (only for existing processes).
        if has_existing_process and molecule_type == "protein" and method == "chemical":
            st.error(t("error_protein_chemical"))
            return

        # Map numeric inputs to the qualitative buckets the existing engine expects.
        # The numeric values are also passed through (separate keys) so downstream
        # modules (concrete_recommendations, recommended_actions, plausibility,
        # PDF) can show them and reason on actual quantities.
        def _purity_label_from_percent(pct: float) -> str:
            # Engine vocabulary: "standard" | "high" | ">99%"
            if pct >= 99.0:
                return ">99%"
            if pct >= 97.0:
                return "high"
            return "standard"

        def _scale_label_from_kg_per_year(kg_per_year: float) -> str:
            # Engine vocabulary: "lab" | "pilot" | "industrial"
            if kg_per_year >= 1000.0:
                return "industrial"
            if kg_per_year >= 1.0:
                return "pilot"
            return "lab"

        def _cost_label_from_eur_per_kg(eur_per_kg: float) -> str:
            # Engine vocabulary: "low" | "medium" | "high"
            if eur_per_kg > 500.0:
                return "high"
            if eur_per_kg >= 10.0:
                return "medium"
            return "low"

        # --- Biophysical numeric → engine bucket mapping ---
        def _aggregation_label_from_tagg(tagg_c: float) -> str:
            # Engine vocabulary: low / medium / high (low = stable)
            if tagg_c >= 65.0:
                return "low"
            if tagg_c >= 50.0:
                return "medium"
            return "high"

        def _folding_label_from_structure(disulfides: int, domains: int, ptm: bool) -> str:
            # Composite: domains + disulfides + glycosylation
            if domains >= 3 or disulfides >= 3 or (domains >= 2 and disulfides >= 2):
                return "high"
            if domains >= 1 and (disulfides >= 1 or ptm):
                return "medium"
            return "low"

        def _stability_label_from_tm(tm_c: float) -> str:
            # Tm bands matching the legacy help-text
            if tm_c >= 70.0:
                return "high"
            if tm_c >= 50.0:
                return "medium"
            return "low"

        def _availability_label_from_sourcing(num_suppliers: int, lead_weeks: float, single_region: bool) -> str:
            # Composite supply-risk to engine label.
            # "low availability" = high supply risk; "high availability" = robust sourcing.
            if num_suppliers <= 1 or lead_weeks >= 8.0:
                return "low"
            if num_suppliers >= 4 and lead_weeks <= 2.0 and not single_region:
                return "high"
            # Single-region with few suppliers downgrades the score one notch.
            if single_region and num_suppliers <= 2:
                return "low"
            return "medium"

        desired_purity = _purity_label_from_percent(float(desired_purity_percent))
        scale = _scale_label_from_kg_per_year(float(scale_kg_per_year))
        raw_material_cost = _cost_label_from_eur_per_kg(float(raw_material_cost_eur_per_kg))
        raw_material_availability = _availability_label_from_sourcing(
            int(num_qualified_suppliers),
            float(lead_time_weeks),
            bool(single_region_concentration),
        )

        processInput = {
            "molecule_name": str(molecule_name).strip(),
            "has_existing_process": bool(has_existing_process),
            "method": method,
            "number_of_steps": int(number_of_steps) if number_of_steps else None,
            "desired_purity": desired_purity,
            "desired_purity_percent": float(desired_purity_percent),
            "molecule_type": molecule_type,
            "molecule_subtype": molecule_subtype,
            "smiles": st.session_state.get("proc_smiles") or None,
            # Full DB entry (when matched) so downstream consumers can show
            # the proper representation — sequence, external_id, etc. — instead
            # of forcing every molecule into a SMILES field.
            "db_entry": st.session_state.get("proc_db_entry") or None,
            "scale": scale,
            "scale_kg_per_year": float(scale_kg_per_year),
            "raw_material_availability": raw_material_availability,
            "num_qualified_suppliers": int(num_qualified_suppliers),
            "lead_time_weeks": float(lead_time_weeks),
            "single_region_concentration": bool(single_region_concentration),
            "raw_material_cost": raw_material_cost,
            "raw_material_cost_eur_per_kg": float(raw_material_cost_eur_per_kg),
            "strict_waste_constraints": bool(strict_waste_constraints),
            "has_advanced_purification": bool(has_advanced_purification),
            "has_bioreactor": bool(has_bioreactor),
            # Detailed lab capabilities (B5): each method as its own boolean flag.
            "available_methods": {k: bool(v) for k, v in purification_methods.items()},
        }
        # Include biomolecule-only biophysical inputs only when relevant.
        # We pass BOTH the qualitative engine-bucket (low/medium/high) AND the
        # numeric values so downstream consumers (recommended_actions,
        # plausibility, PDF report) can reason on real numbers.
        if molecule_type in ("peptide", "protein"):
            processInput["aggregation_risk"] = _aggregation_label_from_tagg(float(tagg_celsius)) if tagg_celsius is not None else "medium"
            processInput["folding_complexity"] = _folding_label_from_structure(
                int(num_disulfides or 0),
                int(num_domains or 1),
                bool(has_ptm),
            )
            processInput["biophysical_stability"] = _stability_label_from_tm(float(tm_celsius)) if tm_celsius is not None else "medium"
            # Numeric values for traceability + downstream rules
            processInput["tagg_celsius"] = float(tagg_celsius) if tagg_celsius is not None else None
            processInput["tm_celsius"] = float(tm_celsius) if tm_celsius is not None else None
            processInput["num_disulfides"] = int(num_disulfides or 0)
            processInput["num_domains"] = int(num_domains or 1)
            processInput["has_ptm"] = bool(has_ptm)

        # Store into session_state for downstream pages
        st.session_state.input_context = processInput

        # Plausibility check (advisory only — does not block the analysis).
        # Surfaces warnings about unusual step counts, impossible scale/equipment
        # combinations, etc., based on a curated knowledge base.
        plaus_result = check_plausibility_with_sources(processInput, lang=st.session_state.ui_lang)
        plausibility_warnings = plaus_result.get("warnings", [])
        plausibility_sources = plaus_result.get("literature_sources", [])
        if plausibility_warnings:
            st.warning(
                t("plausibility_header") + "\n\n"
                + "\n\n".join(f"• {w}" for w in plausibility_warnings)
                + "\n\n" + t("plausibility_continue")
            )
        st.session_state["plausibility_warnings"] = plausibility_warnings
        st.session_state["plausibility_sources"] = plausibility_sources

        with st.spinner(t("spinner_analyzing")):
            try:
                # Use a fresh engine if DB present to avoid stale cache
                if os.path.exists(DB_CSV):
                    fresh_db = MicrobialDB.from_csv(DB_CSV)
                    local_engine = DecisionEngine(db=fresh_db)
                else:
                    local_engine = engine
            except Exception:
                local_engine = engine

            # Call the decision engine analyze wrapper directly
            try:
                recs = local_engine.analyze_process(st.session_state.input_context)
            except Exception:
                recs = local_engine.recommend(None, target_constraints=st.session_state.input_context, top_n=1)

        if not recs:
            st.warning(t("warn_no_results"))
            st.session_state.recommendations = []
            return

        selected = recs[0]
        st.session_state.recommendations = recs
        st.session_state.selected_strategy = selected

        # Generate blueprint + report using existing generators
        bp = generate_production_blueprint(selected, db=db, alternatives_count=2, safety_margin=0.15)
        st.session_state.blueprint = bp
        st.session_state["pdf_bytes"] = None
        try:
            st.session_state.german_report = generate_report(bp)
        except Exception:
            st.session_state.german_report = explain_blueprint(bp)

        # Build the concrete recommendations + recommended actions blocks (C9 + C10).
        # These are stored separately so the report page can render them inline.
        try:
            st.session_state["concrete_recs"] = build_concrete_recommendations(processInput, lang=st.session_state.ui_lang)
        except Exception:
            st.session_state["concrete_recs"] = None
        try:
            st.session_state["recommended_actions"] = build_recommended_actions(processInput, lang=st.session_state.ui_lang)
        except Exception:
            st.session_state["recommended_actions"] = []

        # Navigate to overview (Step 2). User reviews inputs there, then proceeds to report (Step 3).
        st.session_state.page = "overview"


def show_recommendations_page():
    """Show the recommendations list and allow selecting one strategy to generate a blueprint."""
    # NOTE: Candidate selection step has been disabled. This function now
    # short-circuits to the final report if a blueprint already exists, or
    # informs the user and returns to the start page. We keep the original
    # UI code below commented/disabled to preserve it for now.
    if st.session_state.get("blueprint"):
        # If a blueprint was already generated (e.g. auto-selection), show it.
        show_report_page()
        return
    else:
        st.info(t("candidates_disabled"))
        if st.button(t("back_to_home")):
            st.session_state.page = "start"
        return
    st.title("Empfohlene Kandidaten")
    recs = st.session_state.recommendations or []
    # If no recommendations are present, attempt a one-time reload from CSV and re-run the engine
    if not recs:
        try:
            if os.path.exists(DB_CSV):
                fresh_db = MicrobialDB.from_csv(DB_CSV)
                local_engine = DecisionEngine(db=fresh_db)
                recs_fresh = local_engine.recommend(None, target_constraints=st.session_state.input_context, top_n=10)
                if recs_fresh:
                    st.session_state.recommendations = recs_fresh
                    recs = recs_fresh
        except Exception:
            # ignore and fall through to the informational message
            pass

    if not recs:
        st.info("Keine Empfehlungen vorhanden. Bitte starten Sie eine neue Anfrage.")
        # Helpful debug context for the user (non-technical): show input summary and DB size
        try:
            db_rows = db.all().shape[0]
        except Exception:
            db_rows = "unbekannt"
        st.caption(f"Aktuelle Anfrage: {st.session_state.get('input_context', {})} — DB-Einträge: {db_rows}")
        if st.button(t("back_to_home")):
            st.session_state.page = "start"
        return

    # Present a compact recommendation card per top candidate
    st.markdown("---")
    st.subheader("Top-Empfehlungen")

    # Build a user-friendly table showing method/process and scores
    table_rows = []
    for r in recs:
        table_rows.append({
            "Methode": r.get("microorganism"),
            "Prozess-Typ": r.get("strain"),
            "Zielstoff": r.get("compound_name"),
            "Erwarteter Yield": r.get("expected_yield"),
            "Entscheidungs-Score": r.get("decision_score"),
        })

    df_display = pd.DataFrame(table_rows)
    st.dataframe(df_display)

    st.markdown("---")
    st.subheader("Detaillierte Auswahl")
    options = [f"#{r['rank']} — Methode: {r.get('microorganism')} — Score {r.get('decision_score'):.3f}" for r in recs]
    selected_option = st.selectbox("Wählen Sie eine Strategie zur Blueprint-Erstellung", options=options, index=0, key="ui_selected_option")

    sel_idx = options.index(selected_option)
    sel_strategy = recs[sel_idx]

    st.markdown("---")
    st.subheader("Ausgewählte Empfehlung")
    st.write(f"Empfohlene Methode: **{sel_strategy.get('microorganism')}**")
    st.write(f"Prozess-Typ: **{sel_strategy.get('strain')}**")
    st.write(f"Decision Score: **{sel_strategy.get('decision_score'):.2f}**")

    col_a, col_b = st.columns([1, 1])
    with col_a:
        if st.button("Blueprint erzeugen"):
            st.session_state.selected_strategy = sel_strategy
            bp = generate_production_blueprint(sel_strategy, db=db, alternatives_count=2, safety_margin=0.15)
            st.session_state.blueprint = bp
            st.session_state["pdf_bytes"] = None
            try:
                st.session_state.german_report = generate_report(bp)
            except Exception:
                st.session_state.german_report = explain_blueprint(bp)
            st.session_state.page = "report"
    with col_b:
        if st.button(t("back_to_home")):
            st.session_state.page = "start"


def show_report_page():
    """Show the saved blueprint, the German report preview and provide PDF export and feedback UI."""
    # Step indicator: Step 3 (Results) is active here
    render_step_indicator(active_step=3)

    bp = st.session_state.blueprint
    if bp is None:
        st.warning("Kein Blueprint vorhanden. Bitte erzeugen Sie zuerst einen Blueprint.")
        if st.button("Zurück zu Empfehlungen"):
            # Return to the recommendation form instead of the candidate list
            st.session_state.page = "start"
        return

    # Do not render structured data (JSON/dicts) in the UI. Provide a concise, non-technical summary card.
    rec = bp.get("recommended", {}) if isinstance(bp, dict) else {}
    org = rec.get("microorganism") or "nicht angegeben"
    strain = rec.get("strain") or "nicht angegeben"
    stoffklasse = rec.get("compound_class") or bp.get("operating_parameters", {}).get("compound_class") or "nicht angegeben"
    zielstoff = rec.get("compound_name") or "nicht angegeben"

    # Present a concise status line instead of the detailed table
    # Use the molecule name from the recommended blueprint (if present)
    compound_name = None
    if isinstance(bp, dict):
        compound_name = bp.get("recommended", {}).get("compound_name")
    if not compound_name:
        compound_name = "nicht angegeben"
    # Capitalize for nicer presentation
    compound_name = str(compound_name).capitalize()
    st.markdown(f"### " + t("report_success_title").format(compound=compound_name))

    # ----- Eingabe-Übersicht (zeigt die numerischen User-Werte) -----
    ip = st.session_state.get("input_context") or {}
    if ip:
        st.markdown("---")
        st.subheader(f"📥 {t('your_inputs_header')}")
        # Anforderungen + Preise
        col_a, col_b, col_c = st.columns(3)
        pct = ip.get("desired_purity_percent")
        kg = ip.get("scale_kg_per_year")
        eur = ip.get("raw_material_cost_eur_per_kg")
        if pct is not None:
            col_a.metric(t("metric_purity"), f"{pct:.2f} %")
        if kg is not None:
            if kg >= 1000:
                col_b.metric(t("metric_scale"), f"{kg:,.0f} {t('kg_per_year')}")
            else:
                col_b.metric(t("metric_scale"), f"{kg:.3f} {t('kg_per_year')}")
        if eur is not None:
            col_c.metric(t("metric_cost"), f"{eur:,.2f} €/kg")
        # Sourcing-Profil
        n_supp = ip.get("num_qualified_suppliers")
        lead = ip.get("lead_time_weeks")
        single_region = ip.get("single_region_concentration")
        if any(v is not None for v in (n_supp, lead, single_region)):
            col_d, col_e, col_f = st.columns(3)
            if n_supp is not None:
                col_d.metric(t("metric_suppliers"), f"{int(n_supp)}")
            if lead is not None:
                col_e.metric(t("metric_lead_time"), f"{lead:.1f} {t('weeks')}")
            if single_region is not None:
                col_f.metric(t("metric_geo_concentration"), t("yes") if single_region else t("no"))
        # Biophysical profile (only for peptides/proteins)
        tagg = ip.get("tagg_celsius")
        tm = ip.get("tm_celsius")
        n_disulf = ip.get("num_disulfides")
        n_dom = ip.get("num_domains")
        ptm_flag = ip.get("has_ptm")
        if any(v is not None for v in (tagg, tm, n_disulf, n_dom, ptm_flag)):
            st.markdown(f"**{t('biophysical_inputs_header')}:**")
            col_g, col_h, col_i = st.columns(3)
            if tagg is not None:
                col_g.metric(t("metric_tagg"), f"{tagg:.1f} °C")
            if tm is not None:
                col_h.metric(t("metric_tm"), f"{tm:.1f} °C")
            if n_disulf is not None:
                col_i.metric(t("metric_disulfides"), f"{int(n_disulf)}")
            col_j, col_k, _ = st.columns(3)
            if n_dom is not None:
                col_j.metric(t("metric_domains"), f"{int(n_dom)}")
            if ptm_flag is not None:
                col_k.metric(t("metric_ptm"), t("yes") if ptm_flag else t("no"))

    # -------------------------------------------------------------------
    # Feature D — What-if Sensitivity Analysis
    # -------------------------------------------------------------------
    # Process Scientists want to test "what if I push purity up?" or
    # "what if I scale 10x?" without re-running the whole analysis.
    # An expandable panel shows three sliders (scale, purity, raw-material
    # price) and a live-updated COGS-range plot + tile.
    if ip:
        st.markdown("---")
        with st.expander(t("sensitivity_section_title"), expanded=False):
            st.caption(t("sensitivity_section_caption"))

            try:
                from cogs_estimator import (
                    estimate_cogs as _est_sens,
                    estimate_cogs_whatif as _est_whatif,
                    format_cogs_range as _fmt_sens,
                )
            except Exception:
                _est_sens = None
                _est_whatif = None
                _fmt_sens = None

            if _est_sens is None:
                st.warning(t("sensitivity_not_available"))
            else:
                # Defaults from the user's original inputs
                orig_kg = float(ip.get("scale_kg_per_year") or 100.0)
                orig_pct = float(ip.get("desired_purity_percent") or 99.0)
                orig_eur = float(ip.get("raw_material_cost_eur_per_kg") or 50.0)

                col_s1, col_s2, col_s3 = st.columns(3)
                with col_s1:
                    # Log-scale slider for scale — sliders work on a
                    # log10-axis under the hood so users can sweep across
                    # 5+ orders of magnitude in a single drag.
                    log_kg = st.slider(
                        t("sensitivity_scale_label"),
                        min_value=-3.0, max_value=6.0,
                        value=max(-3.0, min(6.0, __import__("math").log10(max(orig_kg, 1e-3)))),
                        step=0.1,
                        format="%.1f",
                        key="sens_log_kg",
                        help=t("sensitivity_scale_help"),
                    )
                    new_kg = 10.0 ** log_kg
                    st.caption(f"≈ {new_kg:,.1f} {t('kg_per_year')}".replace(",", "."))
                with col_s2:
                    new_pct = st.slider(
                        t("sensitivity_purity_label"),
                        min_value=80.0, max_value=99.99,
                        value=float(orig_pct), step=0.1,
                        format="%.2f",
                        key="sens_pct",
                        help=t("sensitivity_purity_help"),
                    )
                with col_s3:
                    new_eur = st.slider(
                        t("sensitivity_rm_label"),
                        min_value=1.0,
                        max_value=max(10000.0, orig_eur * 3.0),
                        value=float(orig_eur),
                        step=max(1.0, orig_eur * 0.05),
                        format="%.0f",
                        key="sens_eur",
                        help=t("sensitivity_rm_help"),
                    )

                # Build a counterfactual process_input (kept for downstream use)
                pi_alt = dict(ip)
                pi_alt["scale_kg_per_year"] = new_kg
                pi_alt["desired_purity_percent"] = new_pct
                pi_alt["raw_material_cost_eur_per_kg"] = new_eur

                # Simulated COGS via the relative what-if projection: it starts
                # from the PDF value at the original inputs (factors == 1.0) and
                # applies smooth, uncapped elasticities for the slider deltas —
                # so RM / purity / scale moves are actually visible (the capped
                # estimate_cogs plateaus instantly for the 79 override molecules).
                try:
                    cogs_alt = _est_whatif(ip, new_kg, new_pct, new_eur) or {}
                except Exception:
                    cogs_alt = {}
                lo_alt = cogs_alt.get("low_eur_per_kg")
                hi_alt = cogs_alt.get("high_eur_per_kg")

                # Comparison metrics next to the original COGS
                try:
                    cogs_orig = _est_sens(ip) or {}
                except Exception:
                    cogs_orig = {}
                lo_o = cogs_orig.get("low_eur_per_kg") or 0
                hi_o = cogs_orig.get("high_eur_per_kg") or 0

                col_m1, col_m2 = st.columns(2)
                with col_m1:
                    st.markdown(f"**{t('sensitivity_original_label')}**")
                    st.metric(
                        "COGS",
                        _fmt_sens(cogs_orig, lang=st.session_state.ui_lang) if _fmt_sens else f"{lo_o:.1f}–{hi_o:.1f} €/kg",
                    )
                with col_m2:
                    st.markdown(f"**{t('sensitivity_simulated_label')}**")
                    delta_hi = (hi_alt - hi_o) if (hi_alt is not None and hi_o) else None
                    st.metric(
                        "COGS",
                        _fmt_sens(cogs_alt, lang=st.session_state.ui_lang) if _fmt_sens else (
                            f"{lo_alt:.1f}–{hi_alt:.1f} €/kg" if lo_alt is not None else "—"
                        ),
                        delta=(f"{delta_hi:+.0f} €/kg (upper)" if delta_hi is not None else None),
                        delta_color="inverse",
                    )

                # Live COGS-vs-Scale plot — gives the user a feel for the
                # scale-economy curve in this specific process.
                try:
                    import math as _math
                    sample_scales = [10 ** (i * 0.5 - 2) for i in range(0, 17)]  # 0.01 → 10^6
                    rows = []
                    for s in sample_scales:
                        try:
                            # What-if projection: vary scale, hold the slider's
                            # purity/RM — shows the real economy-of-scale curve
                            # instead of the capped (flat) point estimate.
                            cs = _est_whatif(ip, s, new_pct, new_eur) or {}
                            rows.append({
                                "Skala (kg/Jahr)": s,
                                "COGS low": cs.get("low_eur_per_kg"),
                                "COGS high": cs.get("high_eur_per_kg"),
                            })
                        except Exception:
                            pass
                    if rows:
                        import pandas as _pd
                        df_plot = _pd.DataFrame(rows).set_index("Skala (kg/Jahr)")
                        st.markdown(f"**{t('sensitivity_chart_title')}**")
                        st.line_chart(df_plot, height=240)
                        st.caption(t("sensitivity_chart_caption"))
                except Exception:
                    pass

                if st.button(t("sensitivity_reset_button"), key="sens_reset"):
                    # Force-reset to original values by clearing the keys
                    # and letting Streamlit re-render.
                    for _k in ("sens_log_kg", "sens_pct", "sens_eur"):
                        st.session_state.pop(_k, None)
                    st.rerun()

    # ----- Structure analysis (RDKit, when SMILES available) -----
    smiles_for_render = ip.get("smiles")
    rdkit_props = None
    try:
        from rdkit_properties import compute_properties as _rdkit_compute, render_structure_png as _rdkit_render, RDKIT_AVAILABLE as _RDKIT_AVAILABLE
        if _RDKIT_AVAILABLE and smiles_for_render:
            rdkit_props = _rdkit_compute(smiles_for_render)
    except Exception:
        rdkit_props = None
    if rdkit_props:
        st.markdown("---")
        st.subheader(f"🧪 {t('structure_analysis_title')}")
        st.caption(t("structure_caption"))
        col_struct, col_metrics = st.columns([1, 1])
        with col_struct:
            try:
                png = _rdkit_render(smiles_for_render, width=350, height=250)
                if png:
                    st.image(png, caption=t("structure_drawing_caption"), use_container_width=True)
            except Exception:
                pass
        with col_metrics:
            st.markdown(f"**{t('structure_label_mw')}:** {rdkit_props['molecular_weight']:.2f} g/mol")
            st.markdown(f"**{t('structure_label_logp')}:** {rdkit_props['logp']:.2f}")
            st.markdown(f"**{t('structure_label_tpsa')}:** {rdkit_props['tpsa']:.2f} Å²")
            st.markdown(f"**{t('structure_label_arom_rings')}:** {rdkit_props['num_aromatic_rings']}")
            st.markdown(f"**{t('structure_label_hbd')}:** {rdkit_props['num_h_donors']}")
            st.markdown(f"**{t('structure_label_hba')}:** {rdkit_props['num_h_acceptors']}")
            st.markdown(f"**{t('structure_label_rotbonds')}:** {rdkit_props['num_rotatable_bonds']}")
            lip_v = rdkit_props['lipinski_violations']
            if lip_v == 0:
                st.markdown(f"✓ {t('structure_label_lipinski_pass')}")
            else:
                st.markdown(f"⚠ {t('structure_label_lipinski_fail')} ({lip_v} {t('structure_label_lipinski')})")

    # ----- Concrete recommendations (C9) -----
    concrete = st.session_state.get("concrete_recs") or {}
    if concrete:
        st.markdown("---")
        st.subheader(f"📊 {t('concrete_recs_header')}")
        if concrete.get("production_route"):
            st.markdown(f"**{t('concrete_recs_route')}:** {concrete['production_route']}")
        if concrete.get("expected_yield"):
            st.markdown(f"**{t('concrete_recs_yield')}:** {concrete['expected_yield']}")
        if concrete.get("processing_time"):
            st.markdown(f"**{t('concrete_recs_time')}:** {concrete['processing_time']}")
        if concrete.get("expected_final_purity"):
            st.markdown(f"**{t('concrete_recs_purity')}:** {concrete['expected_final_purity']}")
        steps = concrete.get("downstream_steps") or []
        if steps:
            st.markdown(f"**{t('concrete_recs_steps')}:**")
            for i, s in enumerate(steps, 1):
                st.markdown(f"{i}. {s}")
        notes = concrete.get("notes") or []
        for n in notes:
            st.caption(n)
        # Inline literature reference for the production hint, when known
        for src in (concrete.get("literature_sources") or []):
            st.caption(f"📚 {t('source_label_inline')}: {src}")

    # ----- Recommended actions (C10) -----
    actions = st.session_state.get("recommended_actions") or []
    if actions:
        st.markdown("---")
        st.subheader(f"⚙️ {t('actions_header')}")
        st.caption(t("actions_caption"))
        # Map raw effort codes (low/medium/high) to localized labels.
        effort_label_keys = {"low": "effort_low", "medium": "effort_medium", "high": "effort_high"}
        for a in actions:
            with st.expander(a.get("title", "")):
                if a.get("rationale"):
                    st.markdown(f"**{t('actions_rationale')}:** {a['rationale']}")
                # Current vs. optimized state — rendered as a side-by-side block
                cur = a.get("current_state")
                opt = a.get("optimized_state")
                if cur or opt:
                    col_cur, col_opt = st.columns(2)
                    if cur:
                        with col_cur:
                            st.markdown(
                                f"<div style='background-color:rgba(255,140,90,0.08); "
                                f"border-left:3px solid #FF8C5A; padding:10px 12px; border-radius:4px; "
                                f"font-size:0.9rem; line-height:1.4;'>"
                                f"<strong style='color:#E8EEF4;'>{t('actions_current_state')}</strong><br>"
                                f"<span style='color:#C8D2DE;'>{cur}</span>"
                                f"</div>",
                                unsafe_allow_html=True,
                            )
                    if opt:
                        with col_opt:
                            st.markdown(
                                f"<div style='background-color:rgba(61,184,197,0.10); "
                                f"border-left:3px solid #3DB8C5; padding:10px 12px; border-radius:4px; "
                                f"font-size:0.9rem; line-height:1.4;'>"
                                f"<strong style='color:#E8EEF4;'>{t('actions_optimized_state')}</strong><br>"
                                f"<span style='color:#C8D2DE;'>{opt}</span>"
                                f"</div>",
                                unsafe_allow_html=True,
                            )
                    st.markdown("<div style='height:8px;'></div>", unsafe_allow_html=True)
                if a.get("expected_impact"):
                    st.markdown(f"**{t('actions_impact')}:** {a['expected_impact']}")
                if a.get("prerequisites"):
                    st.markdown(f"**{t('actions_prerequisites')}:** {a['prerequisites']}")
                if a.get("effort"):
                    eff_key = effort_label_keys.get(str(a["effort"]).lower(), None)
                    eff_label = t(eff_key) if eff_key else a["effort"]
                    st.markdown(f"**{t('actions_effort')}:** {eff_label}")
                if a.get("literature_source"):
                    st.caption(f"📚 {t('source_label_inline')}: {a['literature_source']}")

    # ----- Aggregated references / literature section -----
    aggregated_sources = []
    seen_sources = set()
    def _add_src(s):
        if s and s not in seen_sources:
            seen_sources.add(s)
            aggregated_sources.append(s)
    cr_for_refs = st.session_state.get("concrete_recs") or {}
    for s in (cr_for_refs.get("literature_sources") or []):
        _add_src(s)
    for action_for_ref in (st.session_state.get("recommended_actions") or []):
        _add_src(action_for_ref.get("literature_source"))
    for s in (st.session_state.get("plausibility_sources") or []):
        _add_src(s)
    if aggregated_sources:
        st.markdown("---")
        st.subheader(f"📚 {t('references_title')}")
        st.caption(t("references_intro"))
        for idx, src in enumerate(aggregated_sources, start=1):
            st.markdown(f"**[{idx}]** {src}")

    st.markdown("---")
    # Entferne inline Anzeige des Fließtext-Berichts. Biete nur Vorschau und Export als PDF an.
    # PDF controls: Vorschau anzeigen & Als PDF exportieren
    report_text = st.session_state.german_report or ""

    # Optional: append the multi-route comparison section to the PDF.
    # Toggling this invalidates any cached PDF so the next preview/export
    # reflects the new choice.
    st.checkbox(
        t("append_comparison_label"),
        key="append_comparison",
        help=t("append_comparison_help"),
        on_change=lambda: st.session_state.pop("pdf_bytes", None),
    )

    col_pdf_1, col_pdf_2 = st.columns([1, 1])
    preview_generated = False

    with col_pdf_1:
        if st.button(t("preview_button")):
            try:
                from pdf_exporter import export_report_pdf
            except Exception:
                st.error("PDF-Export nicht verfügbar: reportlab ist nicht installiert. Bitte `pip install reportlab` ausführen.")
            else:
                try:
                    # Only generate if not already cached in session_state
                    if not st.session_state.get("pdf_bytes"):
                        with st.spinner("PDF-Vorschau wird erstellt …"):
                            pdf_extras = {
                                "concrete_recs": st.session_state.get("concrete_recs"),
                                "recommended_actions": st.session_state.get("recommended_actions"),
                                "plausibility_warnings": st.session_state.get("plausibility_warnings"),
                                "plausibility_sources": st.session_state.get("plausibility_sources"),
                            }
                            if st.session_state.get("append_comparison"):
                                pdf_extras["include_comparison"] = True
                                pdf_extras["comparison_input"] = st.session_state.get("input_context")
                            pdf_bytes = export_report_pdf(
                                report_text, bp, None,
                                title=t("pdf_title_default"),
                                extras=pdf_extras,
                                lang=st.session_state.ui_lang,
                            )
                            if not pdf_bytes:
                                raise ValueError("Keine PDF-Bytes erhalten.")
                            st.session_state["pdf_bytes"] = pdf_bytes
                    # Indicate success via the persistent preview area below
                    preview_generated = True
                except Exception as e:
                    st.error(f"PDF-Vorschau fehlgeschlagen: {e}")

    with col_pdf_2:
        if st.button(t("export_pdf_button")):
            # If preview already generated we can reuse bytes cached in session_state
            try:
                pdf_bytes = st.session_state.get("pdf_bytes")
            except Exception:
                pdf_bytes = None

            if pdf_bytes:
                st.download_button(label=t("download_pdf_button"), data=pdf_bytes, file_name="helixar_produktionsbericht.pdf", mime="application/pdf")
            else:
                try:
                    from pdf_exporter import export_report_pdf
                except Exception:
                    st.error("PDF-Export nicht verfügbar: reportlab ist nicht installiert. Bitte `pip install reportlab` ausführen.")
                else:
                    try:
                        # Generate in-memory bytes and offer download
                        pdf_extras = {
                            "concrete_recs": st.session_state.get("concrete_recs"),
                            "recommended_actions": st.session_state.get("recommended_actions"),
                            "plausibility_warnings": st.session_state.get("plausibility_warnings"),
                            "plausibility_sources": st.session_state.get("plausibility_sources"),
                        }
                        if st.session_state.get("append_comparison"):
                            pdf_extras["include_comparison"] = True
                            pdf_extras["comparison_input"] = st.session_state.get("input_context")
                        pdf_bytes = export_report_pdf(
                            report_text, bp, None,
                            title=t("pdf_title_default"),
                            extras=pdf_extras,
                            lang=st.session_state.ui_lang,
                        )
                        st.session_state["pdf_bytes"] = pdf_bytes
                        st.success("PDF wurde erzeugt.")
                        st.download_button(label=t("download_pdf_button"), data=pdf_bytes, file_name="helixar_produktionsbericht.pdf", mime="application/pdf")
                    except Exception as e:
                        st.error(f"PDF-Export fehlgeschlagen: {e}")

    # Persistent preview area: show only a link to open the cached PDF in a new tab.
    # We intentionally avoid embedding or rendering an inline iframe to prevent
    # the large white preview block; users can open the PDF in a new tab instead.
    cached = st.session_state.get("pdf_bytes")
    if cached:
        try:
                        # Render a tiny HTML component that will create a Blob from the
                        # Base64 PDF and open it in a new tab when the user clicks the
                        # link. Running the JS inside a small component iframe is more
                        # reliable for window.open and avoids about:blank.
                        b64 = base64.b64encode(cached).decode("utf-8")
                        html = f"""
<html>
    <body>
        <a id='openPdf' href='#'>Vorschau in neuem Tab öffnen</a>
        <script>
            (function(){{
                const b64 = "{b64}";
                document.getElementById('openPdf').addEventListener('click', function(e){{
                    e.preventDefault();
                    try {{
                        const byteCharacters = atob(b64);
                        const byteNumbers = new Array(byteCharacters.length);
                        for (let i = 0; i < byteCharacters.length; i++) {{
                            byteNumbers[i] = byteCharacters.charCodeAt(i);
                        }}
                        const byteArray = new Uint8Array(byteNumbers);
                        const blob = new Blob([byteArray], {{type: 'application/pdf'}});
                        const url = URL.createObjectURL(blob);
                        window.open(url, '_blank');
                    }} catch(err) {{
                        alert('Vorschau konnte nicht geöffnet werden: ' + err);
                    }}
                }});
            }})();
        </script>
    </body>
</html>
"""
                        # Small height; only the link will be visible.
                        st.components.v1.html(html, height=60)
        except Exception as e:
            st.error(f"PDF-Vorschau konnte nicht erstellt werden: {e}")

    # Note: Inline feedback removed from the report page per UX requirements.

    if st.button("Zurück zu Empfehlungen"):
        # Navigate back to the recommendation *form* (Empfehlung generieren / Startseite)
        # The candidate-selection page is currently disabled, so returning to
        # the form avoids immediately re-entering the report again.
        st.session_state.page = "start"


def show_landing_page():
    """Hero homescreen: title + DNA helix on the right, then 3 info cards and a start CTA."""
    # ---- HERO ROW: text on the left, DNA helix on the right ----
    col_text, col_dna = st.columns([3, 2])
    with col_text:
        st.markdown(
            f"<div style='margin-top:30px;'>"
            f"<div class='hx-hero-title'>{t('hero_title')}</div>"
            f"<div class='hx-hero-subtitle'>{t('hero_subtitle')}</div>"
            f"</div>",
            unsafe_allow_html=True,
        )
    with col_dna:
        st.markdown(render_dna_helix_svg(height_px=320), unsafe_allow_html=True)

    st.markdown("<div style='height:36px;'></div>", unsafe_allow_html=True)

    # ---- THREE INFO CARDS ----
    c1, c2, c3 = st.columns(3, gap="medium")
    with c1:
        st.markdown(
            f"<div class='hx-card'>"
            f"<div class='hx-card-icon'>🧪</div>"
            f"<div class='hx-card-title'>{t('card_smart_analysis_title')}</div>"
            f"<div class='hx-card-body'>{t('card_smart_analysis_body')}</div>"
            f"</div>",
            unsafe_allow_html=True,
        )
    with c2:
        st.markdown(
            f"<div class='hx-card'>"
            f"<div class='hx-card-icon'>📊</div>"
            f"<div class='hx-card-title'>{t('card_optimized_results_title')}</div>"
            f"<div class='hx-card-body'>{t('card_optimized_results_body')}</div>"
            f"</div>",
            unsafe_allow_html=True,
        )
    with c3:
        st.markdown(
            f"<div class='hx-card'>"
            f"<div class='hx-card-icon'>📤</div>"
            f"<div class='hx-card-title'>{t('card_export_share_title')}</div>"
            f"<div class='hx-card-body'>{t('card_export_share_body')}</div>"
            f"</div>",
            unsafe_allow_html=True,
        )

    st.markdown("<div style='height:36px;'></div>", unsafe_allow_html=True)

    # ---- START NEW ANALYSIS CARD ----
    st.markdown(
        f"<div class='hx-card' style='padding:28px 26px;'>"
        f"<div class='hx-card-title' style='font-size:1.25rem;'>{t('start_card_title')}</div>"
        f"<div class='hx-card-body' style='margin-bottom:18px;'>{t('start_card_body')}</div>"
        f"</div>",
        unsafe_allow_html=True,
    )
    # The button itself sits below the descriptive card so styling is consistent
    if st.button(t("start_card_button"), type="primary", use_container_width=False, key="landing_start_btn"):
        st.session_state.page = "start"
        st.rerun()


def render_step_indicator(active_step: int):
    """Render the 1-2-3 step indicator. active_step is 1-indexed (1, 2, or 3).
    Steps before `active_step` show as 'done' (filled, primary color).
    """
    labels = (t("step_inputs"), t("step_analysis"), t("step_results"))
    parts = ["<div class='hx-steps'>"]
    for i, label in enumerate(labels, start=1):
        if i < active_step:
            cls = "hx-step done"
        elif i == active_step:
            cls = "hx-step active"
        else:
            cls = "hx-step"
        parts.append(
            f"<div class='{cls}'><div class='hx-step-circle'>{i}</div>"
            f"<div class='hx-step-label'>{label}</div></div>"
        )
        if i < len(labels):
            line_cls = "hx-step-line done" if i < active_step else "hx-step-line"
            parts.append(f"<div class='{line_cls}'></div>")
    parts.append("</div>")
    st.markdown("".join(parts), unsafe_allow_html=True)


def render_dna_helix_svg(height_px: int = 420) -> str:
    """Generate an inline SVG of a stylised DNA double helix with subtle particles.

    Two strands cross sinusoidally with cross-bars (base pairs) between them.
    Decorative particles around the helix add the 'tech' feel from the mockup.
    Brand colors: marine blue (#1E3A5F) + turquoise (#3DB8C5).
    Transparent background, no external assets.
    """
    import math
    import random
    rnd = random.Random(42)  # fixed seed for deterministic decoration
    width = 200
    height = height_px
    cx = width / 2
    amp = 60     # horizontal swing
    turns = 3.2  # number of helix turns over full height
    n = 80       # samples per strand

    def strand_path(phase: float) -> str:
        pts = []
        for j in range(n + 1):
            t_ = j / n
            y = t_ * height
            x = cx + amp * math.sin(2 * math.pi * turns * t_ + phase)
            pts.append(f"{x:.1f},{y:.1f}")
        return "M " + " L ".join(pts)

    s1 = strand_path(0.0)
    s2 = strand_path(math.pi)

    # Base-pair crossbars between strands at sampled points
    bars = []
    for j in range(0, n + 1, 5):
        t_ = j / n
        y = t_ * height
        x1 = cx + amp * math.sin(2 * math.pi * turns * t_)
        x2 = cx + amp * math.sin(2 * math.pi * turns * t_ + math.pi)
        # Only render when the two strands are visually separated enough
        if abs(x1 - x2) > 6:
            bars.append(
                f"<line x1='{x1:.1f}' y1='{y:.1f}' x2='{x2:.1f}' y2='{y:.1f}' "
                f"stroke='#3DB8C5' stroke-width='1.2' opacity='0.45' />"
            )

    # Decorative particles around the helix
    particles = []
    for _ in range(28):
        px = rnd.uniform(-30, width + 30)
        py = rnd.uniform(0, height)
        pr = rnd.uniform(1.0, 2.6)
        op = rnd.uniform(0.25, 0.85)
        color = "#3DB8C5" if rnd.random() < 0.7 else "#5BCAD7"
        particles.append(
            f"<circle cx='{px:.1f}' cy='{py:.1f}' r='{pr:.1f}' "
            f"fill='{color}' opacity='{op:.2f}' />"
        )

    return f"""
    <svg viewBox='0 0 {width} {height}' xmlns='http://www.w3.org/2000/svg'
         style='display:block; margin:0 auto;' width='{width}' height='{height}'>
      <defs>
        <linearGradient id='hxStrand1' x1='0%' y1='0%' x2='0%' y2='100%'>
          <stop offset='0%' stop-color='#3DB8C5' />
          <stop offset='100%' stop-color='#1E3A5F' />
        </linearGradient>
        <linearGradient id='hxStrand2' x1='0%' y1='0%' x2='0%' y2='100%'>
          <stop offset='0%' stop-color='#1E3A5F' />
          <stop offset='100%' stop-color='#3DB8C5' />
        </linearGradient>
      </defs>
      {''.join(particles)}
      {''.join(bars)}
      <path d='{s1}' fill='none' stroke='url(#hxStrand1)' stroke-width='3.5' stroke-linecap='round' />
      <path d='{s2}' fill='none' stroke='url(#hxStrand2)' stroke-width='3.5' stroke-linecap='round' />
    </svg>
    """


def show_settings_page():
    """Settings page — scaffolded for future expansion (theme, defaults, export, ...)."""
    st.title(f"⚙️  {t('settings_page_title')}")
    st.caption(t("settings_page_intro"))
    st.markdown("---")

    # ---- Section: Language ----
    st.subheader(t("settings_section_language"))
    _opts = list(SUPPORTED_LANGUAGES)
    _label_map = {"de": "Deutsch / German", "en": "English / Englisch"}

    # Current session language
    current_index = _opts.index(st.session_state.get("ui_lang", DEFAULT_LANGUAGE))
    chosen = st.selectbox(
        t("settings_current_language"),
        options=_opts,
        index=current_index,
        format_func=lambda code: _label_map.get(code, code.upper()),
        key="settings_lang_current",
        help=t("settings_current_language_help"),
    )
    if chosen != st.session_state.get("ui_lang"):
        st.session_state.ui_lang = chosen
        st.rerun()

    # Default language for next session (stored in session_state so it survives navigation)
    default_index = _opts.index(st.session_state.get("ui_lang_default", DEFAULT_LANGUAGE))
    default_chosen = st.selectbox(
        t("settings_default_language"),
        options=_opts,
        index=default_index,
        format_func=lambda code: _label_map.get(code, code.upper()),
        key="settings_lang_default",
        help=t("settings_default_language_help"),
    )
    if default_chosen != st.session_state.get("ui_lang_default"):
        st.session_state.ui_lang_default = default_chosen
        st.success(t("settings_saved"))

    # ---- Future sections ----
    # TODO: Add additional sections here (theme preference, default purity / scale,
    # PDF preferences, data-export options, etc.). Pattern:
    #   st.markdown("---")
    #   st.subheader(t("settings_section_<name>"))
    #   ... widgets ...


def show_overview_page():
    """Step 2 — review inputs before triggering the report.
    Shows a compact, sectioned summary of every value the user entered."""
    render_step_indicator(active_step=2)
    st.title(f"📋  {t('overview_title')}")
    st.caption(t("overview_intro"))
    st.markdown("---")

    ip = st.session_state.get("input_context") or {}
    if not ip:
        st.warning(t("warn_no_results"))
        if st.button(t("overview_back")):
            st.session_state.page = "start"
            st.rerun()
        return

    type_label_map = {
        "small_molecule": t("type_small_molecule"),
        "peptide": t("type_peptide"),
        "protein": t("type_protein"),
        "natural_product": t("type_natural_product"),
    }
    method_label_map = {
        "chemical": t("method_chemical"),
        "biotechnological": t("method_biotechnological"),
        "extraction": t("method_extraction"),
    }

    # ---- Molecule ----
    st.subheader(t("overview_section_molecule"))
    col1, col2 = st.columns(2)
    col1.markdown(f"**{t('target_molecule')}:** {ip.get('molecule_name', '—')}")
    col2.markdown(f"**{t('molecule_type')}:** {type_label_map.get(ip.get('molecule_type', ''), '—')}")
    if ip.get("molecule_subtype"):
        col1.markdown(f"**{t('molecule_subtype')}:** {ip['molecule_subtype']}")

    # ---- Process ----
    if ip.get("has_existing_process"):
        st.markdown("---")
        st.subheader(t("overview_section_process"))
        col1, col2 = st.columns(2)
        col1.markdown(f"**{t('current_method_label')}:** {method_label_map.get(ip.get('method', ''), ip.get('method', '—'))}")
        if ip.get("number_of_steps") is not None:
            col2.markdown(f"**{t('number_of_steps_label')}:** {ip['number_of_steps']}")

    # ---- Requirements ----
    st.markdown("---")
    st.subheader(t("overview_section_requirements"))
    col1, col2 = st.columns(2)
    pct = ip.get("desired_purity_percent")
    kg = ip.get("scale_kg_per_year")
    if pct is not None:
        col1.metric(t("metric_purity"), f"{pct:.2f} %")
    if kg is not None:
        if kg >= 1000:
            col2.metric(t("metric_scale"), f"{kg:,.0f} {t('kg_per_year')}")
        else:
            col2.metric(t("metric_scale"), f"{kg:.3f} {t('kg_per_year')}")

    # ---- Market conditions ----
    st.markdown("---")
    st.subheader(t("overview_section_market"))
    col1, col2, col3 = st.columns(3)
    n_supp = ip.get("num_qualified_suppliers")
    lead = ip.get("lead_time_weeks")
    eur = ip.get("raw_material_cost_eur_per_kg")
    if n_supp is not None:
        col1.metric(t("metric_suppliers"), f"{int(n_supp)}")
    if lead is not None:
        col2.metric(t("metric_lead_time"), f"{lead:.1f} {t('weeks')}")
    if eur is not None:
        col3.metric(t("metric_cost"), f"{eur:,.2f} €/kg")
    single_region = ip.get("single_region_concentration")
    if single_region is not None:
        st.markdown(f"**{t('metric_geo_concentration')}:** {t('yes') if single_region else t('no')}")
    if ip.get("strict_waste_constraints"):
        st.markdown(f"⚠️  {t('strict_waste_label')}")

    # ---- Available methods ----
    methods = ip.get("available_methods") or {}
    if methods:
        st.markdown("---")
        st.subheader(t("overview_section_methods"))
        method_keys_to_label = {
            "has_flash_chromatography": t("method_flash_chrom_label"),
            "has_prep_hplc": t("method_prep_hplc_label"),
            "has_rotary_evaporator": t("method_rotavap_label"),
            "has_lyophilizer": t("method_lyo_label"),
            "has_crystallization": t("method_cryst_label"),
            "has_membrane_filtration": t("method_membrane_label"),
            "has_fplc": t("method_fplc_label"),
            "has_extraction": t("method_extraction_label"),
        }
        active = [method_keys_to_label[k] for k, v in methods.items() if v and k in method_keys_to_label]
        if ip.get("has_bioreactor"):
            active.insert(0, t("bioreactor_label"))
        if active:
            for m in active:
                st.markdown(f"- ✓ {m}")
        else:
            st.caption("—")

    # ---- Biophysical (only if peptide/protein) ----
    if ip.get("molecule_type") in ("peptide", "protein"):
        if any(ip.get(k) is not None for k in ("tagg_celsius", "tm_celsius", "num_disulfides", "num_domains", "has_ptm")):
            st.markdown("---")
            st.subheader(t("overview_section_biophys"))
            col1, col2, col3 = st.columns(3)
            if ip.get("tagg_celsius") is not None:
                col1.metric(t("metric_tagg"), f"{ip['tagg_celsius']:.1f} °C")
            if ip.get("tm_celsius") is not None:
                col2.metric(t("metric_tm"), f"{ip['tm_celsius']:.1f} °C")
            if ip.get("num_disulfides") is not None:
                col3.metric(t("metric_disulfides"), f"{int(ip['num_disulfides'])}")
            col4, col5, _ = st.columns(3)
            if ip.get("num_domains") is not None:
                col4.metric(t("metric_domains"), f"{int(ip['num_domains'])}")
            if ip.get("has_ptm") is not None:
                col5.metric(t("metric_ptm"), t("yes") if ip["has_ptm"] else t("no"))

    # ---- Navigation ----
    st.markdown("---")
    col_back, col_continue = st.columns([1, 1])
    with col_back:
        if st.button(t("overview_back"), use_container_width=True):
            st.session_state.page = "start"
            st.rerun()
    with col_continue:
        if st.button(t("overview_continue"), use_container_width=True, type="primary"):
            st.session_state.page = "report"
            st.rerun()


def show_comparison_page(_skip_title: bool = False):
    """Process Comparison — select a molecule and produce a dedicated
    comparison PDF that scores all viable production routes
    side-by-side and ends with a quantified recommendation.

    When called from the recommendation page's mode toggle, set
    `_skip_title=True` to avoid rendering the duplicate page title
    (the parent page already shows the app title and step indicator).
    """
    if not _skip_title:
        st.title(t("comparison_page_title"))
    st.caption(t("comparison_page_intro"))
    st.markdown("<div style='height:14px;'></div>", unsafe_allow_html=True)

    # Molecule picker — uses the same database as the main flow.
    try:
        from molecules_db import MOLECULE_DATABASE
        mol_names = sorted([m["name"] for m in MOLECULE_DATABASE])
    except Exception:
        mol_names = []
    if not mol_names:
        st.warning(t("comparison_no_molecules"))
        return

    sel_name = st.selectbox(
        t("comparison_select_molecule"),
        options=mol_names,
        index=mol_names.index("Vanillin") if "Vanillin" in mol_names else 0,
        help=t("comparison_select_molecule_help"),
    )

    db_entry = next((m for m in MOLECULE_DATABASE if m["name"] == sel_name), None)
    if not db_entry:
        st.error(t("comparison_molecule_not_found"))
        return

    st.markdown("---")
    st.subheader(t("comparison_inputs_header"))

    col_a, col_b = st.columns(2)
    with col_a:
        purity_pct = st.number_input(
            t("purity_label"),
            min_value=80.0, max_value=99.999, value=99.0, step=0.1,
            format="%.2f", key="cmp_purity",
            help=t("purity_help"),
        )
    with col_b:
        scale_kg = st.number_input(
            t("scale_label"),
            min_value=0.001, max_value=1_000_000.0, value=100.0, step=1.0,
            format="%.3f", key="cmp_scale",
            help=t("scale_help"),
        )

    col_c, col_d = st.columns(2)
    with col_c:
        rm_cost = st.number_input(
            t("metric_cost"),
            min_value=0.5, max_value=10_000.0, value=50.0, step=1.0,
            format="%.2f", key="cmp_rm_cost",
            help=t("comparison_rm_help"),
        )
    with col_d:
        n_steps = st.number_input(
            t("comparison_steps_label"),
            min_value=1, max_value=100, value=5, step=1,
            key="cmp_n_steps",
            help=t("comparison_steps_help"),
        )

    st.markdown("---")
    st.markdown(
        f"**{t('comparison_preview_label')}:** "
        f"{db_entry.get('molecule_type', '?')}/{db_entry.get('molecule_subtype', '?')} · "
        f"{db_entry.get('smiles') or db_entry.get('external_id') or '—'}"
    )

    # Build the process_input dict the same way the comparison engine
    # expects it.
    process_input = {
        "molecule_name": sel_name,
        "molecule_type": db_entry.get("molecule_type"),
        "molecule_subtype": db_entry.get("molecule_subtype"),
        "smiles": db_entry.get("smiles"),
        "external_id": db_entry.get("external_id"),
        "method": "chemical",  # placeholder — comparison runs all viable methods
        "number_of_steps": int(n_steps),
        "desired_purity_percent": float(purity_pct),
        "scale_kg_per_year": float(scale_kg),
        "raw_material_cost_eur_per_kg": float(rm_cost),
        "num_qualified_suppliers": 3,
        "lead_time_weeks": 4.0,
        "single_region_concentration": False,
        "strict_waste_constraints": False,
        "has_bioreactor": True,
        "has_advanced_purification": True,
        "available_methods": {
            "has_flash_chromatography": True, "has_prep_hplc": True,
            "has_rotary_evaporator": True, "has_lyophilizer": True,
            "has_crystallization": True, "has_membrane_filtration": True,
            "has_fplc": True, "has_extraction": True,
        },
    }

    st.markdown("<div style='height:10px;'></div>", unsafe_allow_html=True)

    if st.button(t("comparison_generate_button"), type="primary"):
        try:
            from comparison_report import export_comparison_pdf
            with st.spinner(t("comparison_generating_spinner")):
                pdf_bytes = export_comparison_pdf(
                    process_input, lang=st.session_state.get("ui_lang", "de")
                )
            st.session_state["cmp_pdf_bytes"] = pdf_bytes
            st.session_state["cmp_pdf_molecule"] = sel_name
            st.success(t("comparison_success"))
        except Exception as e:
            st.error(f"{t('comparison_error_prefix')}: {e}")

    # Show download button + new-tab preview link if a comparison PDF
    # has been generated. Mirrors the pattern used on the recommendation
    # report page: a tiny HTML component opens the PDF as a Blob in a
    # new browser tab so users get an inline preview without an
    # embedded iframe.
    pdf_bytes_dl = st.session_state.get("cmp_pdf_bytes")
    pdf_mol_dl = st.session_state.get("cmp_pdf_molecule")
    if pdf_bytes_dl and pdf_mol_dl:
        st.markdown("<div style='height:8px;'></div>", unsafe_allow_html=True)
        col_dl, col_pv = st.columns(2)
        with col_dl:
            st.download_button(
                label=t("comparison_download_button").format(molecule=pdf_mol_dl),
                data=pdf_bytes_dl,
                file_name=f"Helixar_Prozessvergleich_{pdf_mol_dl}.pdf",
                mime="application/pdf",
                type="secondary",
                use_container_width=True,
            )
        with col_pv:
            try:
                b64 = base64.b64encode(pdf_bytes_dl).decode("utf-8")
                _preview_label = t("comparison_preview_button").format(molecule=pdf_mol_dl)
                html = (
                    "<html><body style='margin:0; font-family: -apple-system, "
                    "BlinkMacSystemFont, sans-serif;'>"
                    "<a id='openCmpPdf' href='#' style='display:inline-block; "
                    "padding:8px 14px; background:#F0F2F6; color:#1E3A5F; "
                    "border-radius:6px; text-decoration:none; font-weight:600; "
                    "font-size:14px;'>"
                    + _preview_label +
                    "</a>"
                    "<script>(function(){"
                    f"const b64=\"{b64}\";"
                    "document.getElementById('openCmpPdf').addEventListener('click',function(e){"
                    "e.preventDefault();try{"
                    "const bc=atob(b64);const bn=new Array(bc.length);"
                    "for(let i=0;i<bc.length;i++){bn[i]=bc.charCodeAt(i);}"
                    "const ba=new Uint8Array(bn);"
                    "const blob=new Blob([ba],{type:'application/pdf'});"
                    "const url=URL.createObjectURL(blob);window.open(url,'_blank');"
                    "}catch(err){alert('Vorschau konnte nicht geöffnet werden: '+err);}"
                    "});})();</script>"
                    "</body></html>"
                )
                st.components.v1.html(html, height=46)
            except Exception as e:
                st.error(f"{t('comparison_preview_error_prefix')}: {e}")

    # Inline preview of the comparison table (visible without download).
    st.markdown("---")
    st.subheader(t("comparison_inline_preview_header"))
    try:
        from route_comparison import compare_routes, METHODS_DE, METHODS_EN
        cmp = compare_routes(process_input)
        if not cmp.rows:
            st.info(t("comparison_no_routes"))
        else:
            import pandas as _pd
            lang = st.session_state.get("ui_lang", "de")
            methods_table = METHODS_EN if lang == "en" else METHODS_DE
            tier_de = {
                "efficiency": {"low": "Kritisch", "medium": "Optimierbar", "high": "Robust", "very high": "Best-in-class"},
                "cost":       {"low": "Commodity-Niveau", "medium": "Mid-Tier-API", "high": "Specialty-Tier", "very high": "Biologika-Niveau"},
                "risk":       {"low": "Etabliert", "medium": "Akzeptabel", "high": "Überwachungsbedürftig", "very high": "Kritisch"},
            }
            tier_en = {
                "efficiency": {"low": "Critical", "medium": "Optimisable", "high": "Robust", "very high": "Best-in-class"},
                "cost":       {"low": "Commodity-tier", "medium": "Mid-tier API", "high": "Specialty-tier", "very high": "Biologic-tier"},
                "risk":       {"low": "Established", "medium": "Acceptable", "high": "Monitoring-required", "very high": "Critical"},
            }
            tier_table = tier_en if lang == "en" else tier_de
            def _tier(label, metric):
                key = str(label or "medium").lower()
                key = {"niedrig": "low", "mittel": "medium", "hoch": "high", "sehr hoch": "very high"}.get(key, key)
                return tier_table.get(metric, {}).get(key, str(label).title())

            rows = []
            for r in cmp.rows:
                cogs_str = (f"{r.cogs_low_eur_per_kg:.0f}–{r.cogs_high_eur_per_kg:.0f}"
                            if r.cogs_low_eur_per_kg is not None else "—")
                rows.append({
                    "★" if lang != "en" else "Rank":     "★" if r.decision_rank == 1 else f"{r.decision_rank}.",
                    t("comparison_col_method"):           methods_table.get(r.method, r.method),
                    t("comparison_col_cost"):             _tier(r.cost_label, "cost"),
                    t("comparison_col_risk"):             _tier(r.risk_label, "risk"),
                    t("comparison_col_efficiency"):       _tier(r.efficiency_label, "efficiency"),
                    t("comparison_col_cogs"):             cogs_str,
                    t("comparison_col_yield"):            r.expected_yield or "—",
                    t("comparison_col_steps"):            f"{r.n_production_steps}" if r.has_concrete_steps
                                                          else ("generic" if lang == "en" else "generisch"),
                })
            st.dataframe(_pd.DataFrame(rows), hide_index=True, use_container_width=True)

            # Engine reasoning below the table
            top_method = methods_table.get(cmp.recommended_method, cmp.recommended_method)
            reason = (cmp.recommendation_reason_en if lang == "en"
                      else cmp.recommendation_reason_de)
            st.markdown(
                f"**{t('comparison_engine_recommendation')}:** "
                f":blue[**{top_method}**]"
            )
            if reason:
                st.caption(reason)
    except Exception as e:
        st.warning(f"{t('comparison_inline_error_prefix')}: {e}")


def show_demos_page():
    """Pre-built demo scenarios — visual card grid (mockup style)."""
    st.title(t("demo_page_title"))
    st.caption(t("demo_page_intro"))
    st.markdown("<div style='height:18px;'></div>", unsafe_allow_html=True)

    # Try to import RDKit for 2D structure rendering. Fallback gracefully.
    try:
        from rdkit_properties import render_structure_png as _rdk_render, RDKIT_AVAILABLE as _RDK_AVAIL
    except Exception:
        _rdk_render = None
        _RDK_AVAIL = False

    # Generic antibody silhouette (inline SVG) for protein demos when RDKit cannot render
    def _antibody_svg() -> str:
        return """
        <svg viewBox='0 0 200 160' xmlns='http://www.w3.org/2000/svg'
             style='display:block; margin:0 auto;'>
          <g stroke='#3DB8C5' stroke-width='4' stroke-linecap='round' fill='none'>
            <!-- Two heavy chains forming Y-shape -->
            <line x1='65' y1='20' x2='100' y2='80' />
            <line x1='135' y1='20' x2='100' y2='80' />
            <line x1='100' y1='80' x2='100' y2='140' />
            <!-- Light chains (shorter, parallel to upper part of heavy chains) -->
            <line x1='80' y1='35' x2='95' y2='65' stroke='#1E3A5F' />
            <line x1='120' y1='35' x2='105' y2='65' stroke='#1E3A5F' />
            <!-- Disulfide / hinge dot -->
            <circle cx='100' cy='80' r='3' fill='#3DB8C5' stroke='none' />
          </g>
        </svg>
        """

    # Map each demo to a category label (for the small tag chip)
    type_to_chip = {
        "small_molecule": "Small molecule",
        "peptide": "Peptide",
        "protein": "Antibody" if True else "Protein",
        "natural_product": "Natural product",
    }

    cols = st.columns(len(DEMOS), gap="medium")
    for col, demo in zip(cols, DEMOS):
        inp = demo.get("inputs", {})
        smiles = None
        # Try to find a SMILES via molecules_db lookup
        try:
            from molecules_db import get_entry_by_name
            entry = get_entry_by_name(inp.get("proc_molecule_name", ""))
            if entry:
                smiles = entry.get("smiles")
        except Exception:
            pass

        mtype = inp.get("proc_molecule_type_v2") or "small_molecule"
        chip_label = type_to_chip.get(mtype, "Molecule")

        with col:
            # Card content as a single HTML block (chip + image placeholder + title + subtitle)
            st.markdown(
                f"<div class='hx-demo-card'>"
                f"<div class='hx-demo-tag'>{chip_label}</div>",
                unsafe_allow_html=True,
            )
            # Render structure
            if smiles and _RDK_AVAIL and _rdk_render:
                try:
                    png = _rdk_render(smiles, width=320, height=180)
                    if png:
                        st.image(png, use_container_width=True)
                    else:
                        st.markdown(_antibody_svg(), unsafe_allow_html=True)
                except Exception:
                    st.markdown(_antibody_svg(), unsafe_allow_html=True)
            else:
                # No SMILES (large protein / antibody) — show generic antibody silhouette
                st.markdown(_antibody_svg(), unsafe_allow_html=True)

            # Title + subtitle (story preview)
            st.markdown(
                f"<div style='font-weight:600; font-size:1rem; color:#E8EEF4; margin-top:10px;'>{t(demo['title_key'])}</div>"
                f"<div style='font-size:0.85rem; color:#8AA0B8; line-height:1.4; margin-top:4px;'>{t(demo['story_key'])}</div>"
                f"</div>",
                unsafe_allow_html=True,
            )
            st.markdown("<div style='height:8px;'></div>", unsafe_allow_html=True)
            if st.button(t("demo_load_button"), key=f"load_demo_{demo['id']}", use_container_width=True):
                # Copy demo inputs into session_state for the form widgets
                for k, v in demo["inputs"].items():
                    st.session_state[k] = v
                # Reset analysis artefacts
                st.session_state["recommendations"] = []
                st.session_state["selected_strategy"] = None
                st.session_state["blueprint"] = None
                st.session_state["german_report"] = None
                st.session_state["concrete_recs"] = None
                st.session_state["recommended_actions"] = []
                st.session_state["plausibility_warnings"] = []
                st.session_state["pdf_bytes"] = None
                st.session_state.page = "start"
                st.rerun()


def _looks_like_email(s: str) -> bool:
    """Very lightweight email-shape check (not full RFC 5322).

    Empty string is considered valid because the email field is optional.
    """
    if not s or not str(s).strip():
        return True
    s = str(s).strip()
    if " " in s:
        return False
    if s.count("@") != 1:
        return False
    local, _, domain = s.partition("@")
    if not local or not domain or "." not in domain:
        return False
    return True


def show_feedback_page():
    """Feedback page with structured form that emails Lorenz via Formsubmit.co."""
    from feedback_submission import submit_feedback as _submit_feedback

    st.title(t("feedback_page_title"))
    st.markdown(t("feedback_page_intro"))
    st.caption(t("feedback_optional_hint"))

    # Category options — labels are localised, the value sent to email is the
    # short English key so Lorenz sees a stable category name in his inbox.
    category_options = [
        ("bug", t("feedback_category_bug")),
        ("feature", t("feedback_category_feature")),
        ("content", t("feedback_category_content")),
        ("ux", t("feedback_category_ux")),
        ("other", t("feedback_category_other")),
    ]
    category_keys = [k for k, _ in category_options]
    category_label_for = {k: lbl for k, lbl in category_options}

    with st.form("feedback_form", clear_on_submit=False):
        col1, col2 = st.columns(2)
        with col1:
            fb_name = st.text_input(
                t("feedback_label_name"),
                value="",
                placeholder=t("feedback_placeholder_name"),
            )
        with col2:
            fb_company = st.text_input(
                t("feedback_label_company"),
                value="",
                placeholder=t("feedback_placeholder_company"),
            )

        fb_email = st.text_input(
            t("feedback_label_email"),
            value="",
            placeholder=t("feedback_placeholder_email"),
        )

        fb_category_key = st.selectbox(
            t("feedback_label_category"),
            options=category_keys,
            format_func=lambda k: category_label_for[k],
            index=0,
        )

        fb_text = st.text_area(
            t("feedback_label_text"),
            value="",
            placeholder=t("feedback_placeholder_text"),
            height=180,
        )

        submitted = st.form_submit_button(t("feedback_submit_button"))

    if submitted:
        # 1) Required-text check
        if not fb_text or not str(fb_text).strip():
            st.warning(t("feedback_warn_empty_text"))
            return

        # 2) Lightweight email-shape check (only if user filled it)
        if fb_email and not _looks_like_email(fb_email):
            st.warning(t("feedback_email_invalid"))
            return

        # 3) Send via Formsubmit. We keep the category key (stable) plus the
        #    localised label in parentheses so Lorenz sees both.
        category_display = f"{fb_category_key} ({category_label_for[fb_category_key]})"
        active_lang = st.session_state.get("language", DEFAULT_LANGUAGE)

        with st.spinner(t("feedback_sending")):
            result = _submit_feedback(
                feedback_text=fb_text,
                category=category_display,
                name=fb_name,
                company=fb_company,
                user_email=fb_email,
                app_language=active_lang,
            )

        if result.get("ok"):
            st.success(t("feedback_submit_success"))
        else:
            err_code = result.get("message", "")
            if err_code == "needs_activation":
                # Friendly path: the mail service is set up but waiting for the
                # one-time activation click in Lorenz's inbox.
                st.warning(t("feedback_needs_activation"))
            elif err_code in ("timeout", "connection_error"):
                st.error(t("feedback_submit_error_network"))
            else:
                st.error(t("feedback_submit_error_generic"))
            # Show a tiny diagnostic line so we can debug if needed (not alarming)
            st.caption(
                f"(Status: {result.get('http_status') or '—'} · {err_code})"
            )


# Sidebar header: logo + wordmark, then icon-driven navigation.
with st.sidebar:
    if _logo_for_tab:
        # Compact logo in the sidebar header
        col_logo, col_brand = st.columns([1, 2])
        with col_logo:
            st.image(_logo_for_tab, width=44)
        with col_brand:
            st.markdown(
                "<div style='padding-top:6px; font-size:18px; font-weight:700; "
                "color:#E8EEF4; letter-spacing:0.5px;'>HELIXAR <span style='color:#3DB8C5;'>AI</span></div>",
                unsafe_allow_html=True,
            )
    else:
        st.markdown(
            "<div style='font-size:18px; font-weight:700; color:#E8EEF4; padding:8px 0;'>"
            "HELIXAR <span style='color:#3DB8C5;'>AI</span></div>",
            unsafe_allow_html=True,
        )
    st.markdown("<div style='height:8px;'></div>", unsafe_allow_html=True)

# Icon-prefixed navigation labels
_nav_keys = ["home", "generate", "demos", "feedback", "settings"]
_nav_label_keys = {
    "home": "nav_home",
    "generate": "nav_generate",
    "demos": "nav_demos",
    "feedback": "nav_feedback",
    "settings": "nav_settings",
}
_nav_icons = {
    "home": "🏠",
    "generate": "✨",
    "demos": "🧬",
    "feedback": "💬",
    "settings": "⚙️",
}

page = st.sidebar.radio(
    t("nav_label"),
    options=_nav_keys,
    format_func=lambda k: f"{_nav_icons[k]}  {t(_nav_label_keys[k])}",
    label_visibility="collapsed",
)

# Render the global language switch in the top area BEFORE each page renders.
_render_lang_switch()

# --- Page dispatcher: sidebar wins when changed; otherwise honor explicit nav.
# Why two sources: buttons (e.g. "Demo laden", "Bericht anzeigen") set
# st.session_state.page directly to navigate into / between flow steps. But the
# user must also be able to leave that flow via the sidebar. We detect a
# sidebar click by comparing the current sidebar selection to the last-known
# one — if it changed, the user explicitly clicked the sidebar and we discard
# the button-driven explicit_page.
_last_sidebar = st.session_state.get("_last_sidebar_page")
if _last_sidebar is not None and page != _last_sidebar:
    # Sidebar selection changed → user clicked sidebar → reset button-driven page
    st.session_state["page"] = None
st.session_state["_last_sidebar_page"] = page

explicit_page = st.session_state.get("page")
if explicit_page in ("recommendations", "report", "start", "overview"):
    if explicit_page == "recommendations":
        show_recommendations_page()
    elif explicit_page == "report":
        show_report_page()
    elif explicit_page == "start":
        show_start_page()
    elif explicit_page == "overview":
        show_overview_page()
else:
    # fallback to sidebar-driven routing using internal page codes
    if page == "home":
        show_landing_page()
    elif page == "generate":
        show_start_page()
    elif page == "demos":
        show_demos_page()
    elif page == "feedback":
        show_feedback_page()
    elif page == "settings":
        show_settings_page()
    else:
        show_landing_page()

# Footer / Hinweise
st.markdown("---")
st.markdown(t("footer_disclaimer"))
with st.expander(t("notices_title")):
    st.markdown(f"**{t('notice_compliance_title')}** — {t('notice_compliance_body')}")
    st.markdown(f"**{t('notice_data_title')}** — {t('notice_data_body')}")
    st.markdown(f"**{t('notice_ip_title')}** — {t('notice_ip_body')}")
    st.markdown(f"**{t('notice_liability_title')}** — {t('notice_liability_body')}")
# ...existing code...