# Structured molecule database for analysis and reporting
# Each entry contains: name, molecule_type, molecule_subtype, smiles

MOLECULE_DATABASE = [
    # SMALL MOLECULES - volatile
    {"name": "Ethanol", "molecule_type": "small_molecule", "molecule_subtype": "volatile", "smiles": "CCO"},
    {"name": "Acetone", "molecule_type": "small_molecule", "molecule_subtype": "volatile", "smiles": "CC(=O)C"},
    {"name": "Ethyl acetate", "molecule_type": "small_molecule", "molecule_subtype": "volatile", "smiles": "CCOC(=O)C"},
    {"name": "Methanol", "molecule_type": "small_molecule", "molecule_subtype": "volatile", "smiles": "CO"},
    {"name": "Isopropanol", "molecule_type": "small_molecule", "molecule_subtype": "volatile", "smiles": "CC(C)O"},

    # SMALL MOLECULES - non_volatile
    {"name": "Vanillin", "molecule_type": "small_molecule", "molecule_subtype": "non_volatile", "smiles": "COc1cc(C=O)ccc1O"},
    {"name": "Ibuprofen", "molecule_type": "small_molecule", "molecule_subtype": "non_volatile", "smiles": "CC(C)Cc1ccc(cc1)C@@HC(=O)O"},
    {"name": "Glucose", "molecule_type": "small_molecule", "molecule_subtype": "non_volatile", "smiles": "C(C1C(C(C(C(O1)O)O)O)O)O"},
    {"name": "Aspirin", "molecule_type": "small_molecule", "molecule_subtype": "non_volatile", "smiles": "CC(=O)Oc1ccccc1C(=O)O"},


    # PEPTIDES (SIMPLIFIED SMILES) - linear
    {"name": "Glutathione", "molecule_type": "peptide", "molecule_subtype": "linear", "smiles": "NCC(=O)NCC(=O)NCC(=O)O"},
    {"name": "Leu-enkephalin", "molecule_type": "peptide", "molecule_subtype": "linear", "smiles": "NCC(=O)NCC(=O)NCC(=O)NCC(=O)O"},
    {"name": "Met-enkephalin", "molecule_type": "peptide", "molecule_subtype": "linear", "smiles": "NCC(=O)NCC(=O)NCC(=O)NCC(=O)O"},
    {"name": "Bradykinin", "molecule_type": "peptide", "molecule_subtype": "linear", "smiles": "NCC(=O)NCC(=O)NCC(=O)NCC(=O)NCC(=O)O"},
    {"name": "Angiotensin II", "molecule_type": "peptide", "molecule_subtype": "linear", "smiles": "NCC(=O)NCC(=O)NCC(=O)NCC(=O)NCC(=O)O"},

    # PEPTIDES - cyclic (simplified/truncated SMILES)
    {"name": "Cyclosporine", "molecule_type": "peptide", "molecule_subtype": "cyclic", "smiles": "NCC(=O)NCC(=O)NCC(=O)1"},
    {"name": "Gramicidin S", "molecule_type": "peptide", "molecule_subtype": "cyclic", "smiles": "NCC(=O)NCC(=O)1"},
    {"name": "Bacitracin", "molecule_type": "peptide", "molecule_subtype": "cyclic", "smiles": "NCC(=O)NCC(=O)1"},
    {"name": "Vancomycin", "molecule_type": "peptide", "molecule_subtype": "cyclic", "smiles": "NCC(=O)NCC(=O)1"},
    {"name": "Daptomycin", "molecule_type": "peptide", "molecule_subtype": "cyclic", "smiles": "NCC(=O)NCC(=O)1"},

    # PROTEINS - antibodies (simplified core signatures)
    {"name": "Adalimumab", "molecule_type": "protein", "molecule_subtype": "antibody", "smiles": "NCC(=O)NCC(=O)NCC(=O)NCC(=O)"},
    {"name": "Trastuzumab", "molecule_type": "protein", "molecule_subtype": "antibody", "smiles": "NCC(=O)NCC(=O)NCC(=O)NCC(=O)"},
    {"name": "Rituximab", "molecule_type": "protein", "molecule_subtype": "antibody", "smiles": "NCC(=O)NCC(=O)NCC(=O)NCC(=O)"},
    {"name": "Bevacizumab", "molecule_type": "protein", "molecule_subtype": "antibody", "smiles": "NCC(=O)NCC(=O)NCC(=O)NCC(=O)"},
    {"name": "Infliximab", "molecule_type": "protein", "molecule_subtype": "antibody", "smiles": "NCC(=O)NCC(=O)NCC(=O)NCC(=O)"},

    # PROTEINS - enzymes (simplified)
    {"name": "Amylase", "molecule_type": "protein", "molecule_subtype": "enzyme", "smiles": "NCC(=O)NCC(=O)NCC(=O)"},
    {"name": "Lipase", "molecule_type": "protein", "molecule_subtype": "enzyme", "smiles": "NCC(=O)NCC(=O)NCC(=O)"},
    {"name": "Protease", "molecule_type": "protein", "molecule_subtype": "enzyme", "smiles": "NCC(=O)NCC(=O)NCC(=O)"},
    {"name": "Cellulase", "molecule_type": "protein", "molecule_subtype": "enzyme", "smiles": "NCC(=O)NCC(=O)NCC(=O)"},
    {"name": "Lactase", "molecule_type": "protein", "molecule_subtype": "enzyme", "smiles": "NCC(=O)NCC(=O)NCC(=O)"},

    # NATURAL PRODUCTS - terpenes
    {"name": "Linalool", "molecule_type": "natural_product", "molecule_subtype": "terpene", "smiles": "CC(C)=CCCC@HCO"},
    {"name": "Citral", "molecule_type": "natural_product", "molecule_subtype": "terpene", "smiles": "CC(C)=CCC/C(C)=C/C=O"},
    {"name": "Limonene", "molecule_type": "natural_product", "molecule_subtype": "terpene", "smiles": "CC1=CCC(CC1)C(=C)C"},
    {"name": "Menthol", "molecule_type": "natural_product", "molecule_subtype": "terpene", "smiles": "CC(C)C1CCC(C(C)C)CC1O"},
    {"name": "Geraniol", "molecule_type": "natural_product", "molecule_subtype": "terpene", "smiles": "CC(C)=CCC/C(C)=C/CO"},
]


# Helper lookups
def _lower(x: str) -> str:
    try:
        return str(x).strip().lower()
    except Exception:
        return ""


def get_entry_by_name(name: str):
    key = _lower(name)
    for e in MOLECULE_DATABASE:
        if _lower(e.get("name", "")) == key:
            return e
    return None


def get_smiles_for(name: str) -> str:
    e = get_entry_by_name(name)
    if e and e.get("smiles"):
        return e.get("smiles")
    return "Not available"


def list_all():
    return list(MOLECULE_DATABASE)
