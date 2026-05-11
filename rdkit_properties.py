"""Cheminformatics wrapper around RDKit.

Design goals:
  1. Real structure analysis (MW, logP, TPSA, aromatic rings, H-bond
     donors/acceptors) via RDKit when SMILES is available.
  2. Graceful fallback: if RDKit is not installed (e.g. on a minimal local
     environment) or if the SMILES is invalid, the wrapper returns None for
     all properties — callers should handle that case (typically by falling
     back to the legacy string-heuristic in decision_engine).
  3. No silent crashes: every entry point is wrapped in try/except.

Why wrap RDKit instead of using it directly:
  - Single import location → easy to mock in tests
  - Controlled error surface (no rdkit.Chem.rdchem.AtomKekulizeException
    bubbling up to the engine)
  - Clear "is RDKit available" check via RDKIT_AVAILABLE constant
"""

from typing import Optional, Dict, Any
import io
import logging


# --- Try-import RDKit. Failure is non-fatal: the module still loads, callers
# --- can check RDKIT_AVAILABLE before calling structure-aware features.
try:
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Crippen, AllChem, Draw
    from rdkit import RDLogger
    # Silence RDKit's noisy warnings about parsing edge cases.
    RDLogger.DisableLog("rdApp.*")
    RDKIT_AVAILABLE = True
except Exception as _rdkit_import_error:  # pragma: no cover
    RDKIT_AVAILABLE = False
    Chem = None  # type: ignore
    Descriptors = None  # type: ignore
    Crippen = None  # type: ignore
    AllChem = None  # type: ignore
    Draw = None  # type: ignore


# --- Property extraction ---

def compute_properties(smiles: Optional[str]) -> Optional[Dict[str, Any]]:
    """Compute physicochemical properties from a SMILES string.

    Returns a dict with the following numeric / boolean fields, or None when:
      - RDKit is not available
      - SMILES is empty / None
      - SMILES cannot be parsed by RDKit (invalid syntax)

    Output dict keys:
      - molecular_weight: float, g/mol
      - logp: float, octanol/water log partition (Crippen method)
      - tpsa: float, topological polar surface area, Å²
      - num_aromatic_rings: int
      - num_aromatic_atoms: int
      - num_h_donors: int   (Lipinski H-bond donors)
      - num_h_acceptors: int (Lipinski H-bond acceptors)
      - num_rotatable_bonds: int
      - num_heavy_atoms: int
      - num_stereocenters: int
      - lipinski_violations: int (number of Rule-of-5 violations)
      - is_aromatic: bool

    All values are deterministic for a given SMILES (no random sampling).
    """
    if not RDKIT_AVAILABLE or not smiles or not isinstance(smiles, str):
        return None
    s = smiles.strip()
    if not s:
        return None
    try:
        mol = Chem.MolFromSmiles(s)
    except Exception:
        return None
    if mol is None:
        return None

    try:
        mw = float(Descriptors.MolWt(mol))
        logp = float(Crippen.MolLogP(mol))
        tpsa = float(Descriptors.TPSA(mol))
        n_arom_rings = int(Descriptors.NumAromaticRings(mol))
        n_arom_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
        n_hbd = int(Descriptors.NumHDonors(mol))
        n_hba = int(Descriptors.NumHAcceptors(mol))
        n_rot = int(Descriptors.NumRotatableBonds(mol))
        n_heavy = int(mol.GetNumHeavyAtoms())
        try:
            stereo = Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False)
            n_stereo = len(stereo)
        except Exception:
            n_stereo = 0

        # Lipinski Rule of 5
        violations = sum([
            mw > 500,
            logp > 5,
            n_hbd > 5,
            n_hba > 10,
        ])

        return {
            "molecular_weight": round(mw, 2),
            "logp": round(logp, 2),
            "tpsa": round(tpsa, 2),
            "num_aromatic_rings": n_arom_rings,
            "num_aromatic_atoms": n_arom_atoms,
            "num_h_donors": n_hbd,
            "num_h_acceptors": n_hba,
            "num_rotatable_bonds": n_rot,
            "num_heavy_atoms": n_heavy,
            "num_stereocenters": n_stereo,
            "lipinski_violations": int(violations),
            "is_aromatic": n_arom_rings > 0,
        }
    except Exception:
        logging.exception("rdkit_properties.compute_properties failed for SMILES: %r", s)
        return None


# --- Qualitative buckets derived from numeric properties ---
# These map the precise RDKit numbers back to the qualitative
# vocabulary the existing engine uses (low / medium / high), so we can
# slot RDKit in without rewriting the entire engine.

def derive_engine_properties(rdkit_props: Optional[Dict[str, Any]]) -> Optional[Dict[str, str]]:
    """Translate RDKit numerics into the qualitative buckets the engine uses.

    Returns a dict shaped like the legacy `_extract_smiles_properties` output:
      - purification_difficulty: low | medium | high
      - stability: low | medium | high
      - toxicity: low | medium | high
      - complexity: low | medium | high
      - raw_material_complexity: low | medium | high
      - scalability: low | medium | high
      - base_toxicity: low | medium | high
      - aromatic_rings: int
      - long_aliphatic_chain: bool
      - polar_group_count: int
      - ester_present, amide_present, crystallization_potential: bool

    Plus the raw numeric fields under their original names — the engine can
    use either qualitative or quantitative as it likes.
    """
    if not rdkit_props:
        return None
    p = rdkit_props

    mw = p["molecular_weight"]
    logp = p["logp"]
    tpsa = p["tpsa"]
    n_rot = p["num_rotatable_bonds"]
    n_arom = p["num_aromatic_rings"]
    n_heavy = p["num_heavy_atoms"]
    n_hbd = p["num_h_donors"]
    n_hba = p["num_h_acceptors"]
    n_stereo = p["num_stereocenters"]
    lip_viol = p["lipinski_violations"]

    # Complexity heuristic — heavy atoms + rotatable bonds + stereocenters
    if n_heavy > 50 or n_rot > 12 or n_stereo > 4:
        complexity = "high"
    elif n_heavy > 25 or n_rot > 6 or n_stereo > 1:
        complexity = "medium"
    else:
        complexity = "low"

    # Purification difficulty — driven by logP, polar groups, MW
    # Highly polar (low logP, high TPSA) → easier to handle aqueously
    # Highly lipophilic (high logP) → may require organic solvents
    # Many H-bond donors+acceptors → harder chromatography
    polar_score = n_hbd + n_hba
    if logp > 5 or (logp < 0 and polar_score > 8) or n_arom > 4:
        purification_difficulty = "high"
    elif (1 < logp < 4) and polar_score < 6 and n_arom <= 2:
        purification_difficulty = "low"
    else:
        purification_difficulty = "medium"

    # Stability — aromatic rings stabilize, many esters/rotatable bonds destabilize
    # Rough proxy: aromatic count + low rotatable bonds = stable
    if n_arom >= 1 and n_rot < 8:
        stability = "high"
    elif n_rot > 12 or n_heavy > 40:
        stability = "low"
    else:
        stability = "medium"

    # Toxicity — heuristic: halogens raise toxicity; we infer from atom symbols
    # We don't have direct atom access here; use logP as a rough proxy
    # (very lipophilic compounds tend to bioaccumulate).
    if logp > 5:
        toxicity = "high"
    elif logp < 0:
        toxicity = "low"
    else:
        toxicity = "medium"

    raw_material_complexity = complexity  # same heuristic
    # Scalability — Lipinski violations + complexity
    if lip_viol >= 2 or complexity == "high":
        scalability = "low"
    elif lip_viol == 0 and complexity == "low":
        scalability = "high"
    else:
        scalability = "medium"

    # Crystallization tends to work for rigid aromatic structures with
    # moderate rotatable-bond count
    crystallization_potential = (n_arom >= 1) and (n_rot <= 6)

    out = {
        "purification_difficulty": purification_difficulty,
        "stability": stability,
        "toxicity": toxicity,
        "complexity": complexity,
        "raw_material_complexity": raw_material_complexity,
        "scalability": scalability,
        "base_toxicity": toxicity,
        "aromatic_rings": n_arom,
        "long_aliphatic_chain": (n_rot > 8 and n_arom == 0),
        "polar_group_count": polar_score,
        "ester_present": False,   # could detect via SMARTS but keep simple
        "amide_present": False,
        "crystallization_potential": crystallization_potential,
        # Numeric values for richer reporting:
        "rdkit_mw": mw,
        "rdkit_logp": logp,
        "rdkit_tpsa": tpsa,
        "rdkit_h_donors": n_hbd,
        "rdkit_h_acceptors": n_hba,
        "rdkit_rotatable_bonds": n_rot,
        "rdkit_lipinski_violations": lip_viol,
    }
    return out


# --- Structure rendering ---

def render_structure_png(smiles: Optional[str], width: int = 350, height: int = 250) -> Optional[bytes]:
    """Render a SMILES string to a PNG byte-string.

    Returns the PNG bytes, or None if RDKit is unavailable / SMILES is
    invalid. Callers can pass the bytes to st.image() or to reportlab's
    Image flowable directly.
    """
    if not RDKIT_AVAILABLE or not smiles or not isinstance(smiles, str):
        return None
    s = smiles.strip()
    if not s:
        return None
    try:
        mol = Chem.MolFromSmiles(s)
        if mol is None:
            return None
        # MolToImage returns a PIL image — convert to PNG bytes.
        img = Draw.MolToImage(mol, size=(width, height))
        buf = io.BytesIO()
        img.save(buf, format="PNG")
        buf.seek(0)
        return buf.read()
    except Exception:
        logging.exception("rdkit_properties.render_structure_png failed for SMILES: %r", s)
        return None
