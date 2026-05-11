"""Pre-built demo scenarios for sales / pitch use.

Each scenario is a dict with:
  - id (str): stable identifier
  - title_key, story_key, outcome_key: i18n keys for localized text
  - inputs (dict): a complete set of values that maps 1:1 to the form fields
                   in show_start_page(). Loading a demo means copying these
                   values into st.session_state under the matching widget keys.

Demo selection criteria:
  - Vanillin biotech: classic small-molecule example, well-known to industry
  - Trastuzumab biotech: shows the antibody / protein representation logic
  - Aspirin chemical: pharma small-molecule with clean chemical synthesis path

Adding a new demo: append a dict to DEMOS, add the i18n keys, and the new
scenario shows up in the Demo page automatically.
"""

from typing import Dict, Any, List


DEMOS: List[Dict[str, Any]] = [
    {
        "id": "vanillin_biotech",
        "title_key": "demo_vanillin_title",
        "story_key": "demo_vanillin_story",
        "outcome_key": "demo_vanillin_outcome",
        "inputs": {
            # Molecule block
            "proc_molecule_name": "Vanillin",
            "proc_molecule_type_v2": "small_molecule",
            "proc_molecule_subtype_v2": "non_volatile",
            # Mode + Section 1
            "proc_no_existing": False,
            "proc_method_v2": "biotechnological",
            "proc_number_of_steps": 3,
            # Section 2: Anforderungen
            "proc_purity_percent": 98.5,
            "proc_scale_kg_per_year": 100.0,
            # Section 3: Marktbedingungen
            "proc_num_suppliers": 3,
            "proc_lead_time_weeks": 4.0,
            "proc_single_region": False,
            "proc_raw_cost_eur_per_kg": 50.0,
            "proc_strict_waste": False,
            # Section 4: Equipment
            "proc_has_bioreactor": True,
            "proc_has_flash_chromatography": False,
            "proc_has_prep_hplc": False,
            "proc_has_rotary_evaporator": True,
            "proc_has_lyophilizer": False,
            "proc_has_crystallization": True,
            "proc_has_membrane_filtration": False,
            "proc_has_fplc": False,
            "proc_has_extraction": True,
        },
    },
    {
        "id": "trastuzumab_biotech",
        "title_key": "demo_trastuzumab_title",
        "story_key": "demo_trastuzumab_story",
        "outcome_key": "demo_trastuzumab_outcome",
        "inputs": {
            "proc_molecule_name": "Trastuzumab",
            "proc_molecule_type_v2": "protein",
            "proc_molecule_subtype_v2": "antibody",
            "proc_no_existing": False,
            "proc_method_v2": "biotechnological",
            "proc_number_of_steps": 5,
            "proc_purity_percent": 99.5,
            "proc_scale_kg_per_year": 200.0,
            "proc_num_suppliers": 2,
            "proc_lead_time_weeks": 8.0,
            "proc_single_region": False,
            "proc_raw_cost_eur_per_kg": 5000.0,
            "proc_strict_waste": True,
            "proc_has_bioreactor": True,
            "proc_has_flash_chromatography": False,
            "proc_has_prep_hplc": False,
            "proc_has_rotary_evaporator": False,
            "proc_has_lyophilizer": True,
            "proc_has_crystallization": False,
            "proc_has_membrane_filtration": True,
            "proc_has_fplc": True,
            "proc_has_extraction": False,
            # Biophysical (peptide/protein only) — numeric IgG defaults
            "proc_tagg_celsius": 68.0,        # typical antibody Tagg
            "proc_tm_celsius": 75.0,           # robust IgG melting point
            "proc_num_disulfides": 16,         # IgG standard: 12 intrachain + 4 interchain
            "proc_num_domains": 12,            # IgG standard: 4 heavy × 2 + 2 light × 2
            "proc_has_ptm": True,              # glycosylation on Fc
        },
    },
    {
        "id": "aspirin_chemical",
        "title_key": "demo_aspirin_title",
        "story_key": "demo_aspirin_story",
        "outcome_key": "demo_aspirin_outcome",
        "inputs": {
            "proc_molecule_name": "Aspirin",
            "proc_molecule_type_v2": "small_molecule",
            "proc_molecule_subtype_v2": "non_volatile",
            "proc_no_existing": False,
            "proc_method_v2": "chemical",
            "proc_number_of_steps": 1,
            "proc_purity_percent": 99.0,
            "proc_scale_kg_per_year": 5000.0,
            "proc_num_suppliers": 5,
            "proc_lead_time_weeks": 2.0,
            "proc_single_region": False,
            "proc_raw_cost_eur_per_kg": 8.0,
            "proc_strict_waste": True,
            "proc_has_bioreactor": False,
            "proc_has_flash_chromatography": True,
            "proc_has_prep_hplc": False,
            "proc_has_rotary_evaporator": True,
            "proc_has_lyophilizer": False,
            "proc_has_crystallization": True,
            "proc_has_membrane_filtration": False,
            "proc_has_fplc": False,
            "proc_has_extraction": True,
        },
    },
]


def get_demo_by_id(demo_id: str):
    """Return the demo dict matching the given id, or None."""
    for d in DEMOS:
        if d.get("id") == demo_id:
            return d
    return None
