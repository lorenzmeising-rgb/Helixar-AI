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
from molecules_db import get_entry_by_name, MOLECULE_DATABASE

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

st.set_page_config(page_title="Helixar · Produktions-Blueprint", layout="centered")


# --- Session state initialization --------------------------------------------
if "page" not in st.session_state:
    # Landing page on app load; separate from the recommendation form/page
    st.session_state.page = "landing"
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
    st.title("Helixar · Produktions-Blueprint-Generator")
    st.caption("Rule-based decision support for process optimization (non-operational).")

    # New: process-optimization input form
    st.subheader("Describe Your Process")
    st.markdown("Describe your existing production process so we can analyze it and suggest improvements.")
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
    if "proc_autofill_msg" not in st.session_state:
        st.session_state.proc_autofill_msg = ""

    def _lock_autofill():
        st.session_state.proc_autofill_locked = True

    def _try_autofill():
        name = str(st.session_state.get("proc_molecule_name") or "").strip()
        if not name:
            return
        # If name changed, allow autofill again
        if name != st.session_state.proc_last_autofill_name:
            st.session_state.proc_autofill_locked = False
        if st.session_state.proc_autofill_locked:
            return
        # exact lookup
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
            try:
                mt = entry.get("molecule_type") or "small_molecule"
                label_map = {
                    "small_molecule": "Small molecule",
                    "peptide": "Peptide",
                    "protein": "Protein",
                    "natural_product": "Natural product",
                }
                st.session_state.proc_molecule_type_label = label_map.get(mt, "Small molecule")
                subtype = entry.get("molecule_subtype")
                subtype_map_all = {
                    "small_molecule": {"volatile": "Volatile", "non_volatile": "Non-volatile"},
                    "peptide": {"linear": "Linear", "cyclic": "Cyclic"},
                    "protein": {"antibody": "Antibody", "enzyme": "Enzyme"},
                    "natural_product": {"terpene": "Terpene", "alkaloid": "Alkaloid"},
                }
                label = subtype_map_all.get(mt, {}).get(subtype)
                if label:
                    st.session_state.proc_molecule_subtype_label = label
                else:
                    st.session_state.proc_molecule_subtype_label = "None"
                st.session_state.proc_smiles = entry.get("smiles")
                st.session_state.proc_last_autofill_name = name
                st.session_state.proc_autofill_msg = f"Molecule recognised: {mt} / {subtype or 'None'}"
            except Exception:
                pass

    molecule_name = st.text_input("Target Molecule", value="", key="proc_molecule_name", on_change=_try_autofill)

    # 2. Current production method
    method_options = {"chemical": "chemical", "biotechnological": "biotechnological", "extraction": "extraction"}
    selected_method_label = st.selectbox("Current Production Method", list(method_options.keys()), index=0, key="proc_method_label")
    method = method_options[selected_method_label]

    # 3. Process complexity (number_of_steps)
    number_of_steps = st.number_input("Number of Synthesis Steps", min_value=1, max_value=10, value=3, step=1, key="proc_number_of_steps")

    # 4. Purity requirement
    purity_map = {"standard": "standard", "high": "high", "very high": ">99%"}
    selected_purity_label = st.selectbox("Purity Requirement", list(purity_map.keys()), index=0, key="proc_purity_label")
    desired_purity = purity_map[selected_purity_label]

    # 5. Production scale
    scale_map = {"lab": "lab", "pilot": "pilot", "industrial": "industrial"}
    selected_scale_label = st.selectbox("Production Scale", list(scale_map.keys()), index=0, key="proc_scale_label")
    scale = scale_map[selected_scale_label]

    # 6/9. Critical infrastructure (booleans)
    has_advanced_purification = st.checkbox("Advanced purification equipment available", value=False, key="proc_has_advanced_purification")
    has_bioreactor = st.checkbox("Bioreactor available", value=False, key="proc_has_bioreactor")

    # 8. Raw materials
    raw_avail_map = {"low": "low", "medium": "medium", "high": "high"}
    selected_raw_avail = st.selectbox("Raw material availability", list(raw_avail_map.keys()), index=1, key="proc_raw_avail")
    raw_material_availability = raw_avail_map[selected_raw_avail]
    raw_cost_map = {"low": "low", "medium": "medium", "high": "high"}
    selected_raw_cost = st.selectbox("Raw material cost", list(raw_cost_map.keys()), index=1, key="proc_raw_cost")
    raw_material_cost = raw_cost_map[selected_raw_cost]

    # 8. Waste / toxicity constraints
    strict_waste_constraints = st.checkbox("Strict waste handling constraints", value=False, key="proc_strict_waste")

    # NEW: Molecule type (small_molecule, peptide, protein, natural_product)
    molec_type_map = {
        "Small molecule": "small_molecule",
        "Peptide": "peptide",
        "Protein": "protein",
        "Natural product": "natural_product",
    }
    selected_molec_label = st.selectbox("Molecule Type", list(molec_type_map.keys()), index=0, key="proc_molecule_type_label", on_change=_lock_autofill)
    molecule_type = molec_type_map[selected_molec_label]

    # Optional: Molecule subtype (depends on molecule_type)
    subtype_options_map = {
        "small_molecule": [("None", None), ("Volatile", "volatile"), ("Non-volatile", "non_volatile")],
        "peptide": [("None", None), ("Linear", "linear"), ("Cyclic", "cyclic")],
        "protein": [("None", None), ("Antibody", "antibody"), ("Enzyme", "enzyme")],
        "natural_product": [("None", None), ("Terpene", "terpene"), ("Alkaloid", "alkaloid")],
    }
    opts = subtype_options_map.get(molecule_type, [("None", None)])
    opt_labels = [o[0] for o in opts]
    sel_sub_label = st.selectbox("Molecule Subtype (optional)", opt_labels, index=0, key="proc_molecule_subtype_label", on_change=_lock_autofill)
    # map back to internal value
    label_to_val = {o[0]: o[1] for o in opts}
    molecule_subtype = label_to_val.get(sel_sub_label)

    # Biomolecule-only biophysical inputs (show only for peptide/protein)
    aggregation_risk = None
    folding_complexity = None
    biophysical_stability = None
    if molecule_type in ("peptide", "protein"):
        st.markdown("---")
        st.subheader("Biophysical properties (only for peptides & proteins)")
        ar_map = {"low": "low", "medium": "medium", "high": "high"}
        sel_ar = st.selectbox("Aggregation risk", list(ar_map.keys()), index=1, key="proc_aggregation_risk")
        aggregation_risk = ar_map[sel_ar]
        fc_map = {"low": "low", "medium": "medium", "high": "high"}
        sel_fc = st.selectbox("Folding complexity", list(fc_map.keys()), index=1, key="proc_folding_complexity")
        folding_complexity = fc_map[sel_fc]
        bs_map = {"low": "low", "medium": "medium", "high": "high"}
        sel_bs = st.selectbox("Biophysical stability", list(bs_map.keys()), index=1, key="proc_biophysical_stability")
        biophysical_stability = bs_map[sel_bs]

    st.markdown("---")
    # Action button: analyze process
    if st.button("Analyze Process"):
        # Validation
        if not molecule_name or str(molecule_name).strip() == "":
            st.warning("Please provide the target molecule (required).")
            return
        if not number_of_steps or int(number_of_steps) < 1:
            st.warning("Please provide the number of synthesis steps (required).")
            return

        processInput = {
            "molecule_name": str(molecule_name).strip(),
            "method": method,
            "number_of_steps": int(number_of_steps),
            "desired_purity": desired_purity,
            "molecule_type": molecule_type,
            "molecule_subtype": molecule_subtype,
            "smiles": st.session_state.get("proc_smiles") or None,
            "scale": scale,
            "raw_material_availability": raw_material_availability,
            "raw_material_cost": raw_material_cost,
            "strict_waste_constraints": bool(strict_waste_constraints),
            "has_advanced_purification": bool(has_advanced_purification),
            "has_bioreactor": bool(has_bioreactor),
        }
        # Include biomolecule-only biophysical inputs only when relevant
        if molecule_type in ("peptide", "protein"):
            processInput["aggregation_risk"] = aggregation_risk
            processInput["folding_complexity"] = folding_complexity
            processInput["biophysical_stability"] = biophysical_stability

        # Store into session_state for downstream pages
        st.session_state.input_context = processInput

        with st.spinner("Prozess wird analysiert …"):
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
            st.warning("Analyse konnte keine Ergebnisse erzeugen. Bitte prüfen Sie die Eingaben.")
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

        # Navigate to report
        st.session_state.page = "report"


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
        st.info("Die Kandidatenauswahl wurde deaktiviert. Bitte nutzen Sie die Startseite, um eine neue Anfrage zu stellen.")
        if st.button("Zurück zur Startseite"):
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
        if st.button("Zurück zur Startseite"):
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
        if st.button("Zurück zur Startseite"):
            st.session_state.page = "start"


def show_report_page():
    """Show the saved blueprint, the German report preview and provide PDF export and feedback UI."""
    st.title("Produktions-Blueprint & Bericht")
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
    st.markdown(f"### Erfolgreiche Produktionsempfehlung für {compound_name} generiert")

    # Entferne inline Anzeige des Fließtext-Berichts. Biete nur Vorschau und Export als PDF an.
    # PDF controls: Vorschau anzeigen & Als PDF exportieren
    report_text = st.session_state.german_report or ""

    col_pdf_1, col_pdf_2 = st.columns([1, 1])
    preview_generated = False

    with col_pdf_1:
        if st.button("Vorschau anzeigen"):
            try:
                from pdf_exporter import export_report_pdf
            except Exception:
                st.error("PDF-Export nicht verfügbar: reportlab ist nicht installiert. Bitte `pip install reportlab` ausführen.")
            else:
                try:
                    # Only generate if not already cached in session_state
                    if not st.session_state.get("pdf_bytes"):
                        with st.spinner("PDF-Vorschau wird erstellt …"):
                            pdf_bytes = export_report_pdf(report_text, bp, None, title="Helixar Produktionsbericht")
                            if not pdf_bytes:
                                raise ValueError("Keine PDF-Bytes erhalten.")
                            st.session_state["pdf_bytes"] = pdf_bytes
                    # Indicate success via the persistent preview area below
                    preview_generated = True
                except Exception as e:
                    st.error(f"PDF-Vorschau fehlgeschlagen: {e}")

    with col_pdf_2:
        if st.button("Als PDF exportieren"):
            # If preview already generated we can reuse bytes cached in session_state
            try:
                pdf_bytes = st.session_state.get("pdf_bytes")
            except Exception:
                pdf_bytes = None

            if pdf_bytes:
                st.download_button(label="PDF herunterladen", data=pdf_bytes, file_name="helixar_produktionsbericht.pdf", mime="application/pdf")
            else:
                try:
                    from pdf_exporter import export_report_pdf
                except Exception:
                    st.error("PDF-Export nicht verfügbar: reportlab ist nicht installiert. Bitte `pip install reportlab` ausführen.")
                else:
                    try:
                        # Generate in-memory bytes and offer download
                        pdf_bytes = export_report_pdf(report_text, bp, None, title="Helixar Produktionsbericht")
                        st.session_state["pdf_bytes"] = pdf_bytes
                        st.success("PDF wurde erzeugt.")
                        st.download_button(label="PDF herunterladen", data=pdf_bytes, file_name="helixar_produktionsbericht.pdf", mime="application/pdf")
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
    """Simple landing / start page with logo and slogan placeholders."""
    # Centered landing page: show logo (if available) centered, then slogan beneath
    # Prefer a project-root image file (place the provided logo as 'helixar_logo.png' or 'logo.png')
    logo_candidates = [
        os.path.join(BASE_DIR, "helixar_logo.png"),
        os.path.join(BASE_DIR, "logo.png"),
        os.path.join(BASE_DIR, "helixar_logo.jpg"),
        os.path.join(BASE_DIR, "logo.jpg"),
    ]
    logo_path = None
    for p in logo_candidates:
        if os.path.exists(p):
            logo_path = p
            break

    st.markdown("<div style='text-align:center; margin-top: 40px;'>", unsafe_allow_html=True)
    if logo_path:
        # Use a three-column layout and place the image in the center column to ensure centering
        col1, col2, col3 = st.columns([1, 2, 1])
        with col2:
            st.image(logo_path, width=280)
    else:
        st.markdown("<div style='font-size:48px; font-weight:600;'>LOGO</div>", unsafe_allow_html=True)

    # Slogan (German) centered under the logo
    st.markdown("<div style='margin-top:12px; font-size:18px; color:#555;'>Von Molekül zu Produktionsstrategie</div>", unsafe_allow_html=True)
    st.markdown("</div>", unsafe_allow_html=True)
    st.markdown("---")
    st.markdown("Wählen Sie im Menü 'Empfehlung generieren', um zu starten.")


def show_feedback_page():
    """A lightweight feedback page, separate from reports/blueprints."""
    st.title("Feedback geben")
    user_feedback = st.text_area("Ihr Feedback", value="", placeholder="Ihre Nachricht an das Team")
    if st.button("Feedback absenden"):
        # For now, do not connect to blueprints or reports; simply acknowledge receipt.
        if not user_feedback or str(user_feedback).strip() == "":
            st.warning("Bitte geben Sie zuerst Feedback ein.")
        else:
            st.success("Danke — Ihr Feedback wurde empfangen.")
            # Optionally, store or forward feedback in a future iteration.


# Sidebar navigation (appears on all pages) - simple radio navigation
st.sidebar.title("Helixar")
st.sidebar.markdown("---")
# Direct radio-driven navigation. Keep labels and order as requested.
page = st.sidebar.radio("Navigation", ["Startseite", "Empfehlung generieren", "Feedback geben"]) 

# --- Page dispatcher: prefer explicit session_state navigation from buttons ---
# If a button set st.session_state.page (e.g. "recommendations" or "report"),
# honor that and render the corresponding page. Otherwise fall back to the
# sidebar radio selection (`page`) so the user can still navigate there.
explicit_page = st.session_state.get("page")
# Only honor explicit_page when a button actually requested a different view.
# Treat the default 'landing' value as non-decisive so the sidebar selection works.
if explicit_page in ("recommendations", "report", "start"):
    if explicit_page == "recommendations":
        show_recommendations_page()
    elif explicit_page == "report":
        show_report_page()
    elif explicit_page == "start":
        show_start_page()
else:
    # fallback to sidebar-driven routing
    if page == "Startseite":
        show_landing_page()
    elif page == "Empfehlung generieren":
        show_start_page()
    elif page == "Feedback geben":
        show_feedback_page()
    else:
        show_landing_page()

# Footer / Hinweise
st.markdown("---")
st.markdown("Hinweis: Dieses Tool bietet nicht-operatives Entscheidungs‑Support. Es ersetzt keine Laborvalidierung oder regulatorische Beratung.")
# ...existing code...