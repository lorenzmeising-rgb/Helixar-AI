"""Curated molecule database for Helixar AI.

Schema (every entry MAY include any subset of these keys):
  - name (str, required)
  - molecule_type (str, required): small_molecule | peptide | protein | natural_product
  - molecule_subtype (str, optional)
  - smiles (str | None): canonical SMILES, or None when not meaningfully representable
                        (proteins, large/complex peptides, antibodies)
  - sequence (str, optional): amino-acid sequence in single-letter code, or descriptive
                              chain notation (e.g. "γ-Glu-Cys-Gly" for Glutathione)
  - external_id (str, optional): authoritative reference, e.g. "UniProt:P04626",
                                 "DrugBank:DB00072", "ChEBI:16856"
  - representation_note (str, optional): explanation when SMILES is None or
                                         when the representation is unusual
  - molecular_weight_kda (float, optional): for proteins/large biomolecules

Design intent:
  Earlier versions of this file used identical fake SMILES strings for all
  antibodies, all enzymes, and most peptides. That is dishonest and immediately
  visible to any chemist. The current version is honest about what can and
  cannot be represented as SMILES:
    - Small molecules: real SMILES.
    - Small peptides (≤3 AA) where SMILES is well-defined: real SMILES.
    - Larger peptides: SMILES omitted, sequence provided.
    - Proteins / antibodies / enzymes: SMILES omitted (≥10⁴ atoms is impractical
      to represent), external IDs provided instead (UniProt for the protein
      itself, DrugBank for therapeutic antibodies — also lists the target).
"""

from typing import Optional, Dict, Any, List


MOLECULE_DATABASE: List[Dict[str, Any]] = [
    # ---------------- SMALL MOLECULES — volatile ----------------
    {"name": "Ethanol", "molecule_type": "small_molecule", "molecule_subtype": "volatile",
     "smiles": "CCO", "external_id": "ChEBI:16236"},
    {"name": "Acetone", "molecule_type": "small_molecule", "molecule_subtype": "volatile",
     "smiles": "CC(=O)C", "external_id": "ChEBI:15347"},
    {"name": "Ethyl acetate", "molecule_type": "small_molecule", "molecule_subtype": "volatile",
     "smiles": "CCOC(=O)C", "external_id": "ChEBI:27750"},
    {"name": "Methanol", "molecule_type": "small_molecule", "molecule_subtype": "volatile",
     "smiles": "CO", "external_id": "ChEBI:17790"},
    {"name": "Isopropanol", "molecule_type": "small_molecule", "molecule_subtype": "volatile",
     "smiles": "CC(C)O", "external_id": "ChEBI:17824"},

    # ---------------- SMALL MOLECULES — non-volatile ----------------
    {"name": "Vanillin", "molecule_type": "small_molecule", "molecule_subtype": "non_volatile",
     "smiles": "COc1cc(C=O)ccc1O", "external_id": "ChEBI:18346"},
    # Ibuprofen: previous SMILES had a malformed "C@@H" without brackets. Fixed.
    {"name": "Ibuprofen", "molecule_type": "small_molecule", "molecule_subtype": "non_volatile",
     "smiles": "CC(C)Cc1ccc(cc1)[C@@H](C)C(=O)O", "external_id": "ChEBI:5855"},
    # Glucose: stereospecified α-D-glucose (most common reference form for fermentation).
    {"name": "Glucose", "molecule_type": "small_molecule", "molecule_subtype": "non_volatile",
     "smiles": "OC[C@H]1O[C@H](O)[C@H](O)[C@@H](O)[C@@H]1O", "external_id": "ChEBI:17234"},
    {"name": "Aspirin", "molecule_type": "small_molecule", "molecule_subtype": "non_volatile",
     "smiles": "CC(=O)Oc1ccccc1C(=O)O", "external_id": "ChEBI:15365"},

    # ---------------- PEPTIDES — linear ----------------
    # Glutathione: small enough (3 AA) to have a well-defined, hand-verifiable SMILES.
    {"name": "Glutathione", "molecule_type": "peptide", "molecule_subtype": "linear",
     "smiles": "N[C@@H](CCC(=O)N[C@@H](CS)C(=O)NCC(=O)O)C(=O)O",
     "sequence": "γ-Glu-Cys-Gly",
     "external_id": "ChEBI:16856"},
    # Enkephalins: 5 AA. Hand-SMILES feasible but error-prone; sequence is the
    # standard representation in pharmacology.
    {"name": "Leu-enkephalin", "molecule_type": "peptide", "molecule_subtype": "linear",
     "smiles": None,
     "sequence": "YGGFL",
     "external_id": "ChEBI:64628",
     "representation_note": "5-AA opioid pentapeptide; represented by amino-acid sequence."},
    {"name": "Met-enkephalin", "molecule_type": "peptide", "molecule_subtype": "linear",
     "smiles": None,
     "sequence": "YGGFM",
     "external_id": "ChEBI:64624",
     "representation_note": "5-AA opioid pentapeptide; represented by amino-acid sequence."},
    {"name": "Bradykinin", "molecule_type": "peptide", "molecule_subtype": "linear",
     "smiles": None,
     "sequence": "RPPGFSPFR",
     "external_id": "UniProt:P01042",
     "representation_note": "9-AA vasoactive peptide; represented by amino-acid sequence."},
    {"name": "Angiotensin II", "molecule_type": "peptide", "molecule_subtype": "linear",
     "smiles": None,
     "sequence": "DRVYIHPF",
     "external_id": "ChEBI:2719",
     "representation_note": "8-AA peptide hormone; represented by amino-acid sequence."},

    # ---------------- PEPTIDES — cyclic ----------------
    # Cyclic peptides with non-natural residues, glycosylation or branching are
    # too complex for hand-written SMILES. We intentionally omit SMILES and
    # reference the authoritative external entry instead.
    {"name": "Cyclosporine", "molecule_type": "peptide", "molecule_subtype": "cyclic",
     "smiles": None,
     "external_id": "ChEBI:4031",
     "representation_note": "11-residue cyclic peptide with N-methylated and non-proteinogenic AAs; SMILES non-trivial — see ChEBI."},
    {"name": "Gramicidin S", "molecule_type": "peptide", "molecule_subtype": "cyclic",
     "smiles": None,
     "sequence": "cyclo(VOLfPVOLfP)",  # D-Phe at lowercase positions
     "external_id": "ChEBI:5306",
     "representation_note": "Cyclic decapeptide antibiotic with two D-Phe residues."},
    {"name": "Bacitracin", "molecule_type": "peptide", "molecule_subtype": "cyclic",
     "smiles": None,
     "external_id": "ChEBI:3066",
     "representation_note": "Mixture of cyclic peptides (Bacitracin A predominant); structure includes thiazoline ring and D-AAs."},
    {"name": "Vancomycin", "molecule_type": "peptide", "molecule_subtype": "cyclic",
     "smiles": None,
     "external_id": "ChEBI:28001",
     "representation_note": "Tricyclic glycopeptide antibiotic with chlorinated aromatic residues — structure too complex for inline SMILES."},
    {"name": "Daptomycin", "molecule_type": "peptide", "molecule_subtype": "cyclic",
     "smiles": None,
     "external_id": "DrugBank:DB00080",
     "representation_note": "Cyclic lipopeptide antibiotic with decanoyl tail — see DrugBank."},

    # ---------------- PROTEINS — antibodies ----------------
    # Antibodies are ~150 kDa proteins (~1300 AAs) — SMILES is not a sensible
    # representation. We use DrugBank IDs (the standard pharmaceutical reference)
    # and identify the molecular target via UniProt for downstream context.
    {"name": "Adalimumab", "molecule_type": "protein", "molecule_subtype": "antibody",
     "smiles": None,
     "external_id": "DrugBank:DB00051",
     "representation_note": "Fully human IgG1 antibody targeting TNF-α (UniProt:P01375).",
     "molecular_weight_kda": 148.0},
    {"name": "Trastuzumab", "molecule_type": "protein", "molecule_subtype": "antibody",
     "smiles": None,
     "external_id": "DrugBank:DB00072",
     "representation_note": "Humanized IgG1 antibody targeting HER2 / ERBB2 (UniProt:P04626).",
     "molecular_weight_kda": 145.5},
    {"name": "Rituximab", "molecule_type": "protein", "molecule_subtype": "antibody",
     "smiles": None,
     "external_id": "DrugBank:DB00073",
     "representation_note": "Chimeric IgG1 antibody targeting CD20 (UniProt:P11836).",
     "molecular_weight_kda": 144.0},
    {"name": "Bevacizumab", "molecule_type": "protein", "molecule_subtype": "antibody",
     "smiles": None,
     "external_id": "DrugBank:DB00112",
     "representation_note": "Humanized IgG1 antibody targeting VEGF-A (UniProt:P15692).",
     "molecular_weight_kda": 149.0},
    {"name": "Infliximab", "molecule_type": "protein", "molecule_subtype": "antibody",
     "smiles": None,
     "external_id": "DrugBank:DB00065",
     "representation_note": "Chimeric IgG1 antibody targeting TNF-α (UniProt:P01375).",
     "molecular_weight_kda": 144.0},

    # ---------------- PROTEINS — enzymes ----------------
    # Generic enzyme names (Amylase, Lipase, ...) cover many isoforms across
    # organisms. We anchor each entry to one well-characterized, commercially
    # relevant variant and reference its UniProt.
    {"name": "Amylase", "molecule_type": "protein", "molecule_subtype": "enzyme",
     "smiles": None,
     "external_id": "UniProt:P04745",
     "representation_note": "Reference: human salivary α-amylase (AMY1A). Industrial production typically uses Bacillus or fungal homologs.",
     "molecular_weight_kda": 57.7},
    {"name": "Lipase", "molecule_type": "protein", "molecule_subtype": "enzyme",
     "smiles": None,
     "external_id": "UniProt:P16233",
     "representation_note": "Reference: human pancreatic lipase (PNLIP). Industrial workhorse: Candida antarctica lipase B (CALB).",
     "molecular_weight_kda": 51.2},
    {"name": "Protease", "molecule_type": "protein", "molecule_subtype": "enzyme",
     "smiles": None,
     "external_id": "UniProt:P00761",
     "representation_note": "Reference: bovine trypsin (well-characterized). 'Protease' is a broad family — specify isoform for production planning.",
     "molecular_weight_kda": 23.3},
    {"name": "Cellulase", "molecule_type": "protein", "molecule_subtype": "enzyme",
     "smiles": None,
     "external_id": "UniProt:P62694",
     "representation_note": "Reference: Trichoderma reesei endoglucanase Cel7B. Cellulase preparations are typically multi-enzyme cocktails.",
     "molecular_weight_kda": 48.2},
    {"name": "Lactase", "molecule_type": "protein", "molecule_subtype": "enzyme",
     "smiles": None,
     "external_id": "UniProt:P09848",
     "representation_note": "Reference: human lactase-phlorizin hydrolase (LCT). Industrial lactase from Aspergillus or Kluyveromyces.",
     "molecular_weight_kda": 218.0},

    # ---------------- NATURAL PRODUCTS — terpenes ----------------
    # Terpenes are small molecules; real SMILES are tractable. Earlier entries
    # had malformed stereo descriptors (e.g. "C@HCO" instead of "[C@H](CO)").
    {"name": "Linalool", "molecule_type": "natural_product", "molecule_subtype": "terpene",
     "smiles": "CC(=CCC[C@](C)(O)C=C)C", "external_id": "ChEBI:17580"},
    {"name": "Citral", "molecule_type": "natural_product", "molecule_subtype": "terpene",
     "smiles": "CC(=CCC/C(=C/C=O)/C)C", "external_id": "ChEBI:16980",
     "representation_note": "Citral is a mixture of geranial (E) and neral (Z); SMILES shows geranial."},
    {"name": "Limonene", "molecule_type": "natural_product", "molecule_subtype": "terpene",
     "smiles": "CC(=C)[C@@H]1CCC(=CC1)C", "external_id": "ChEBI:15384",
     "representation_note": "(R)-(+)-Limonene shown; commercial limonene is often (R)-enantiomer (citrus origin)."},
    {"name": "Menthol", "molecule_type": "natural_product", "molecule_subtype": "terpene",
     "smiles": "CC(C)[C@@H]1CC[C@@H](C)C[C@H]1O", "external_id": "ChEBI:15409",
     "representation_note": "(–)-Menthol (1R,2S,5R) — pharmacopoeia standard form."},
    {"name": "Geraniol", "molecule_type": "natural_product", "molecule_subtype": "terpene",
     "smiles": "CC(=CCC/C(=C/CO)/C)C", "external_id": "ChEBI:17447"},

    # ============================================================
    # EXTENSION SET (added 2026-05): biotech-/industry-relevant
    # molecules across all categories. Roughly doubles the database.
    # ============================================================

    # ---------------- SMALL MOLECULES — volatile (extension) ----------------
    {"name": "Butanol", "molecule_type": "small_molecule", "molecule_subtype": "volatile",
     "smiles": "CCCCO", "external_id": "ChEBI:28885",
     "representation_note": "n-Butanol; biotech via Clostridium ABE-Fermentation."},
    {"name": "Isobutanol", "molecule_type": "small_molecule", "molecule_subtype": "volatile",
     "smiles": "CC(C)CO", "external_id": "ChEBI:46645",
     "representation_note": "Advanced biofuel; engineered E. coli / S. cerevisiae."},
    {"name": "2,3-Butanediol", "molecule_type": "small_molecule", "molecule_subtype": "volatile",
     "smiles": "CC(O)C(C)O", "external_id": "ChEBI:16812",
     "representation_note": "Plattform-Chemikalie; biotech via Klebsiella, Bacillus, Paenibacillus."},
    {"name": "Ethyl lactate", "molecule_type": "small_molecule", "molecule_subtype": "volatile",
     "smiles": "CCOC(=O)C(C)O", "external_id": "ChEBI:39458",
     "representation_note": "Grünes Lösungsmittel aus biotech-Milchsäure + Ethanol."},
    {"name": "Diacetyl", "molecule_type": "small_molecule", "molecule_subtype": "volatile",
     "smiles": "CC(=O)C(=O)C", "external_id": "ChEBI:16583",
     "representation_note": "Buttriger Aroma-Stoff; biotech via Lactobacillus / Lactococcus."},

    # ---------------- SMALL MOLECULES — non-volatile (extension) ----------------
    {"name": "Citric acid", "molecule_type": "small_molecule", "molecule_subtype": "non_volatile",
     "smiles": "OC(=O)CC(O)(CC(=O)O)C(=O)O", "external_id": "ChEBI:30769",
     "representation_note": "Massenstoff, ~2 Mio t/Jahr biotech via A. niger."},
    {"name": "Lactic acid", "molecule_type": "small_molecule", "molecule_subtype": "non_volatile",
     "smiles": "C[C@H](O)C(=O)O", "external_id": "ChEBI:422",
     "representation_note": "Bioplastik-Vorstufe (PLA); biotech via Lactobacillus, B. coagulans."},
    {"name": "Succinic acid", "molecule_type": "small_molecule", "molecule_subtype": "non_volatile",
     "smiles": "OC(=O)CCC(=O)O", "external_id": "ChEBI:15741",
     "representation_note": "Plattform-Chemikalie (DOE Top 12); biotech via A. succinogenes, B. succiniciproducens."},
    {"name": "Itaconic acid", "molecule_type": "small_molecule", "molecule_subtype": "non_volatile",
     "smiles": "C(=C)(C(=O)O)CC(=O)O", "external_id": "ChEBI:30838",
     "representation_note": "Grüne Monomer-Vorstufe; biotech via A. terreus."},
    {"name": "Paracetamol", "molecule_type": "small_molecule", "molecule_subtype": "non_volatile",
     "smiles": "CC(=O)Nc1ccc(O)cc1", "external_id": "ChEBI:46195",
     "representation_note": "Acetaminophen / N-Acetyl-p-aminophenol; klassische chemische Synthese im Massenmaßstab."},

    # ---------------- PEPTIDES — linear (extension) ----------------
    # All ≥29 AA — SMILES omitted; sequence + DrugBank reference.
    {"name": "Glucagon", "molecule_type": "peptide", "molecule_subtype": "linear",
     "smiles": None,
     "sequence": "HSQGTFTSDYSKYLDSRRAQDFVQWLMNT",
     "external_id": "DrugBank:DB00040",
     "representation_note": "29-AA Pankreashormon; rekombinant in E. coli / Hefe.",
     "molecular_weight_kda": 3.48},
    {"name": "Calcitonin", "molecule_type": "peptide", "molecule_subtype": "linear",
     "smiles": None,
     "sequence": "CSNLSTCVLGKLSQELHKLQTYPRTNTGSGTP-NH2 (Lachs-Calcitonin)",
     "external_id": "DrugBank:DB00017",
     "representation_note": "32-AA Hormon (intern. Disulfid 1-7); rekombinant oder SPPS.",
     "molecular_weight_kda": 3.43},
    {"name": "Liraglutide", "molecule_type": "peptide", "molecule_subtype": "linear",
     "smiles": None,
     "sequence": "HAEGTFTSDVSSYLEGQAAK(γGlu-C16)EFIAWLVRGRG",
     "external_id": "DrugBank:DB06655",
     "representation_note": "GLP-1-Analog mit C16-Fettsäure-Modifikation; semisynthetisch (rekombinante Backbone + chemische Acylierung).",
     "molecular_weight_kda": 3.75},
    {"name": "Teriparatide", "molecule_type": "peptide", "molecule_subtype": "linear",
     "smiles": None,
     "sequence": "SVSEIQLMHNLGKHLNSMERVEWLRKKLQDVHNF (PTH 1-34)",
     "external_id": "DrugBank:DB06285",
     "representation_note": "N-terminales Fragment des Parathormons (1-34); rekombinant in E. coli.",
     "molecular_weight_kda": 4.12},
    {"name": "Exenatide", "molecule_type": "peptide", "molecule_subtype": "linear",
     "smiles": None,
     "sequence": "HGEGTFTSDLSKQMEEEAVRLFIEWLKNGGPSSGAPPPS",
     "external_id": "DrugBank:DB01276",
     "representation_note": "GLP-1-Mimetikum (Exendin-4); synthetisch via SPPS.",
     "molecular_weight_kda": 4.19},

    # ---------------- PEPTIDES — cyclic (extension) ----------------
    # Complex cyclic peptides — SMILES omitted, external_id authoritative.
    {"name": "Octreotide", "molecule_type": "peptide", "molecule_subtype": "cyclic",
     "smiles": None,
     "sequence": "fCFwKTCT-ol (Disulfid 2-7)",
     "external_id": "DrugBank:DB00104",
     "representation_note": "8-AA cyclisches Somatostatin-Analog mit D-Phe + D-Trp; SPPS.",
     "molecular_weight_kda": 1.02},
    {"name": "Polymyxin B", "molecule_type": "peptide", "molecule_subtype": "cyclic",
     "smiles": None,
     "external_id": "DrugBank:DB00781",
     "representation_note": "Cyclisches Lipopeptid-Antibiotikum (10 AA + Fettsäure); fermentativ via Bacillus polymyxa."},
    {"name": "Colistin", "molecule_type": "peptide", "molecule_subtype": "cyclic",
     "smiles": None,
     "external_id": "DrugBank:DB00803",
     "representation_note": "Polymyxin E; cyclisches Lipopeptid-Reserveantibiotikum gegen multiresistente Gram-negative."},
    {"name": "Caspofungin", "molecule_type": "peptide", "molecule_subtype": "cyclic",
     "smiles": None,
     "external_id": "DrugBank:DB00520",
     "representation_note": "Echinocandin-Antimykotikum; semisynthetisch aus Pneumocandin B0 (Glarea lozoyensis)."},
    {"name": "Anidulafungin", "molecule_type": "peptide", "molecule_subtype": "cyclic",
     "smiles": None,
     "external_id": "DrugBank:DB06637",
     "representation_note": "Echinocandin-Antimykotikum; semisynthetisch aus Echinocandin B (A. nidulans)."},

    # ---------------- PROTEINS — antibodies (extension) ----------------
    {"name": "Pembrolizumab", "molecule_type": "protein", "molecule_subtype": "antibody",
     "smiles": None,
     "external_id": "DrugBank:DB09037",
     "representation_note": "Humanisierter IgG4 anti-PD-1 (UniProt:Q15116); CHO-Zellen. Top-1 Onkologie-Verkauf weltweit (Keytruda).",
     "molecular_weight_kda": 149.0},
    {"name": "Nivolumab", "molecule_type": "protein", "molecule_subtype": "antibody",
     "smiles": None,
     "external_id": "DrugBank:DB09035",
     "representation_note": "Vollhumaner IgG4 anti-PD-1 (UniProt:Q15116); CHO-Zellen (Opdivo).",
     "molecular_weight_kda": 146.0},
    {"name": "Cetuximab", "molecule_type": "protein", "molecule_subtype": "antibody",
     "smiles": None,
     "external_id": "DrugBank:DB00002",
     "representation_note": "Chimärer IgG1 anti-EGFR (UniProt:P00533); SP2/0-Mausmyelom-Zellen (Erbitux).",
     "molecular_weight_kda": 152.0},
    {"name": "Pertuzumab", "molecule_type": "protein", "molecule_subtype": "antibody",
     "smiles": None,
     "external_id": "DrugBank:DB06366",
     "representation_note": "Humanisierter IgG1 anti-HER2 (UniProt:P04626); CHO-Zellen (Perjeta) — komplementär zu Trastuzumab.",
     "molecular_weight_kda": 148.0},
    {"name": "Daratumumab", "molecule_type": "protein", "molecule_subtype": "antibody",
     "smiles": None,
     "external_id": "DrugBank:DB09331",
     "representation_note": "Vollhumaner IgG1 anti-CD38 (UniProt:P28907); CHO-Zellen (Darzalex), Multiples Myelom.",
     "molecular_weight_kda": 148.0},

    # ---------------- PROTEINS — enzymes (extension) ----------------
    {"name": "Glucoamylase", "molecule_type": "protein", "molecule_subtype": "enzyme",
     "smiles": None,
     "external_id": "UniProt:P69328",
     "representation_note": "A. niger Glucoamylase (GA1); industrieller Workhorse zur Stärkehydrolyse → Glucose-Sirup.",
     "molecular_weight_kda": 70.0},
    {"name": "Phytase", "molecule_type": "protein", "molecule_subtype": "enzyme",
     "smiles": None,
     "external_id": "UniProt:P34754",
     "representation_note": "A. niger Phytase (PhyA); Tierfutter-Additiv (verbessert Phosphor-Verfügbarkeit).",
     "molecular_weight_kda": 52.0},
    {"name": "Glucose isomerase", "molecule_type": "protein", "molecule_subtype": "enzyme",
     "smiles": None,
     "external_id": "UniProt:P24300",
     "representation_note": "Streptomyces rubiginosus; immobilisiert für High-Fructose-Corn-Syrup (HFCS) — eines der großvolumigsten industriellen Enzyme.",
     "molecular_weight_kda": 43.0},
    {"name": "Pectinase", "molecule_type": "protein", "molecule_subtype": "enzyme",
     "smiles": None,
     "external_id": "UniProt:P26213",
     "representation_note": "A. niger Endo-Polygalacturonase; Saftherstellung, Klärung von Wein/Saft.",
     "molecular_weight_kda": 36.0},
    {"name": "Asparaginase", "molecule_type": "protein", "molecule_subtype": "enzyme",
     "smiles": None,
     "external_id": "UniProt:P00805",
     "representation_note": "E. coli L-Asparaginase II (Tetramer ~140 kDa); Krebstherapie bei akuter lymphatischer Leukämie (ALL).",
     "molecular_weight_kda": 34.0},

    # ---------------- NATURAL PRODUCTS — terpenes (extension) ----------------
    {"name": "α-Pinene", "molecule_type": "natural_product", "molecule_subtype": "terpene",
     "smiles": "CC1=CC[C@H]2CC1[C@@]2(C)C", "external_id": "ChEBI:28660",
     "representation_note": "(1S,5S)-α-Pinen; Hauptbestandteil von Kiefernöl, biotech via engineered S. cerevisiae möglich."},
    {"name": "β-Carotene", "molecule_type": "natural_product", "molecule_subtype": "terpene",
     "smiles": "CC(=CCC/C(C)=C/C=C/C(C)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/C1=C(C)CCCC1(C)C)C",
     "external_id": "ChEBI:17579",
     "representation_note": "C40 Tetraterpen, Provitamin A; biotech via Blakeslea trispora oder engineered Yarrowia lipolytica."},
    {"name": "Astaxanthin", "molecule_type": "natural_product", "molecule_subtype": "terpene",
     "smiles": None,
     "external_id": "ChEBI:40968",
     "representation_note": "Carotinoid, C40H52O4; biotech via Haematococcus pluvialis oder Phaffia rhodozyma — Aquakultur-Pigment, MW ~597 g/mol; volle SMILES bei ChEBI."},
    {"name": "Squalene", "molecule_type": "natural_product", "molecule_subtype": "terpene",
     "smiles": "CC(=CCC/C(C)=C/CC/C(C)=C/CC/C=C(C)/CC/C=C(C)/CC/C=C(C)/C)C",
     "external_id": "ChEBI:15440",
     "representation_note": "Lineares Triterpen; biotech via engineered S. cerevisiae als Adjuvans/Kosmetik-Vorstufe."},
    {"name": "Artemisinin", "molecule_type": "natural_product", "molecule_subtype": "terpene",
     "smiles": None,
     "external_id": "ChEBI:223316",
     "representation_note": "Sesquiterpen-Lacton mit Endoperoxidbrücke; semisynthetisch aus engineered Hefe-Vorstufe Artemisininsäure (Sanofi/Amyris-Verfahren). Malaria-Wirkstoff. Volle SMILES bei ChEBI."},

    # ---------------- NATURAL PRODUCTS — alkaloids (NEW category) ----------------
    {"name": "Caffeine", "molecule_type": "natural_product", "molecule_subtype": "alkaloid",
     "smiles": "CN1C=NC2=C1C(=O)N(C)C(=O)N2C", "external_id": "ChEBI:27732",
     "representation_note": "1,3,7-Trimethylxanthin; primär aus Tee/Kaffeebohnen extrahiert (entkoffeinierungs-Nebenprodukt) oder synthetisch."},
    {"name": "Morphine", "molecule_type": "natural_product", "molecule_subtype": "alkaloid",
     "smiles": "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5",
     "external_id": "ChEBI:17303",
     "representation_note": "Phenanthren-Alkaloid aus Papaver somniferum; Extraktion aus Schlafmohn-Schlafmilch (Latex) oder aus dem Strohextrakt."},
    {"name": "Codeine", "molecule_type": "natural_product", "molecule_subtype": "alkaloid",
     "smiles": "CN1CC[C@]23c4c5ccc(OC)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5",
     "external_id": "ChEBI:16714",
     "representation_note": "3-O-Methylmorphin; meist halbsynthetisch durch Methylierung von Morphin (höhere Wirtschaftlichkeit als Direktextraktion)."},
    {"name": "Atropine", "molecule_type": "natural_product", "molecule_subtype": "alkaloid",
     "smiles": "CN1[C@H]2CC[C@@H]1C[C@H](OC(=O)C(CO)c1ccccc1)C2",
     "external_id": "ChEBI:16708",
     "representation_note": "Tropan-Alkaloid (rac. Hyoscyamin); Extraktion aus Atropa belladonna oder Datura. Anticholinergikum."},
    {"name": "Quinine", "molecule_type": "natural_product", "molecule_subtype": "alkaloid",
     "smiles": "COc1ccc2nccc([C@H](O)[C@@H]3C[C@H]4CC[N@@]3C[C@@H]4C=C)c2c1",
     "external_id": "ChEBI:15854",
     "representation_note": "Chinolin-Alkaloid aus Cinchona-Rinde; klassische Malaria-Therapie, immer noch hauptsächlich extrahiert (Totalsynthese unwirtschaftlich)."},
]


# ---------- Helper lookups ----------

def _norm(x: str) -> str:
    try:
        return str(x).strip().lower()
    except Exception:
        return ""


def get_entry_by_name(name: str) -> Optional[Dict[str, Any]]:
    """Return the full entry dict matching `name` (case-insensitive), or None."""
    key = _norm(name)
    for e in MOLECULE_DATABASE:
        if _norm(e.get("name", "")) == key:
            return e
    return None


def get_smiles_for(name: str) -> str:
    """Return the SMILES for the named molecule, or 'Not available' when SMILES
    is not a meaningful representation (proteins, complex peptides).

    Kept as a string return for backwards compatibility with the rest of the
    codebase; callers that need richer info should use get_entry_by_name().
    """
    e = get_entry_by_name(name)
    if not e:
        return "Not available"
    s = e.get("smiles")
    if s:
        return s
    return "Not available"


def get_representation_label(entry: Dict[str, Any]) -> str:
    """Build a short, human-readable description of how a molecule is represented.

    Useful for UI hints and PDF caption lines. Examples:
      - "SMILES: COc1cc(C=O)ccc1O"
      - "Sequence: RPPGFSPFR (UniProt:P01042)"
      - "DrugBank:DB00072 — humanized IgG1 antibody targeting HER2 (~145 kDa)"
    """
    if not entry:
        return "Nicht verfügbar"
    smi = entry.get("smiles")
    seq = entry.get("sequence")
    ext = entry.get("external_id")
    note = entry.get("representation_note")
    mw = entry.get("molecular_weight_kda")

    if smi:
        ext_part = f" ({ext})" if ext else ""
        return f"SMILES: {smi}{ext_part}"
    if seq:
        ext_part = f" ({ext})" if ext else ""
        return f"Sequenz: {seq}{ext_part}"
    # Protein / complex case — lead with the external ID and the note.
    parts: List[str] = []
    if ext:
        parts.append(ext)
    if note:
        parts.append(note)
    if mw:
        parts.append(f"~{mw:.1f} kDa")
    return " — ".join(parts) if parts else "Repräsentation nicht verfügbar"


def list_all() -> List[Dict[str, Any]]:
    return list(MOLECULE_DATABASE)
