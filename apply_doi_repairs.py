"""
apply_doi_repairs.py
====================
Wendet die Ergebnisse von repair_dois.py auf die Code-Dateien an.

Logik:
  - HIGH-Vorschläge: alte DOI → neue DOI ersetzen
  - MEDIUM-Vorschläge: nur die explizit akzeptierten (MANUAL_ACCEPTS) ersetzen,
    der Rest wird wie LOW/NONE behandelt
  - LOW + NONE + abgelehnte MEDIUM: DOI aus dem literature_source-String entfernen
    (Klartext-Zitierung bleibt erhalten)

Backups: alle modifizierten Dateien werden als <name>.bak.<timestamp> gesichert.

Idempotent: wenn DOI schon korrekt ersetzt, passiert nichts.
"""

from __future__ import annotations

import csv
import re
import shutil
import sys
import time
from pathlib import Path
from typing import Dict, List, Tuple

REPO = Path(__file__).resolve().parent

FILES_TO_MODIFY = [
    REPO / "concrete_recommendations.py",
    REPO / "molecule_hints_extension.py",
]

CSV_IN = REPO / "docs" / "doi_repair_proposals.csv"

# ----------------------------------------------------------------------
# Meine manuell-reviewten MEDIUM-Akzeptanzen (Lorenz-Review)
# ----------------------------------------------------------------------
# Schlüssel: (molecule, method, original_doi)
# Wert: True = akzeptiert (neue DOI anwenden), False/missing = ablehnen (löschen)

MANUAL_ACCEPTS: Dict[Tuple[str, str, str], bool] = {
    ("vanillin",       "extraction",        "10.1080/10408390701764235"): True,   # Sinha 2008 Vanilla extraction
    ("angiotensin ii", "chemical",          "10.1084/jem.103.3.295"):    True,   # Skeggs 1956 Hypertensin II
    ("gramicidin s",   "biotechnological",  "10.1016/S1074-5521(97)90019-1"): True,  # Marahiel 1997 verwandtes Paper
    ("daptomycin",     "biotechnological",  "10.1086/648676"):           True,   # Eisenstein 2010
    ("menthol",        "extraction",        "10.1007/10_2014_283"):      True,   # Lange 2015 p-Menthane
    ("atropine",       "extraction",        "10.1002/pca.728"):          True,   # Berkov 2003 Datura
}


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------

def _escape_for_regex(s: str) -> str:
    """Escape DOI string for use in regex."""
    return re.escape(s)


def remove_doi_from_text(text: str, doi: str) -> Tuple[str, int]:
    """Entfernt eine DOI-Erwähnung aus einem Text.

    WICHTIG: arbeitet auf BELIEBIGEN Eingaben (auch ganzen File-Inhalten);
    macht NUR DOI-bezogene Ersetzungen, nicht globale Whitespace-Bereinigung.
    Das war ein katastrophaler Bug in der Vorversion — globales re.sub auf
    "\\s+" hat alle Newlines aus der Python-Datei gefressen.

    Erkennt mehrere Schreibweisen, NUR INNERHALB derselben Zeile (mit "[^\\n]"
    statt "\\s"), damit der Dateistruktur niemals geschadet wird.

    Returns (neuer_text, anzahl_ersetzungen).
    """
    d = _escape_for_regex(doi)
    # ' *' statt '\s*' — keine Newlines fressen!
    patterns = [
        rf" *\( *doi *: *{d} *\)",                  # " (doi:DOI)"
        rf" *\( *DOI *: *{d} *\)",                  # " (DOI: DOI)"
        rf" *\( *{d} *\)",                          # " (DOI)" naked
        rf" *doi *: *{d}",                          # " doi:DOI" ohne Klammern
        rf" *DOI *: *{d}",                          # " DOI:DOI"
        rf" *{d}",                                  # raw DOI
    ]
    total = 0
    for pat in patterns:
        new_text, n = re.subn(pat, "", text, flags=re.IGNORECASE)
        if n > 0:
            text = new_text
            total += n
    # KEIN globales Whitespace-Cleanup mehr.
    return text, total


def replace_doi_in_text(text: str, old_doi: str, new_doi: str) -> Tuple[str, int]:
    """Ersetzt eine alte DOI durch eine neue. Case-insensitive."""
    d = _escape_for_regex(old_doi)
    new_text, n = re.subn(d, new_doi, text, flags=re.IGNORECASE)
    return new_text, n


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------

def main() -> int:
    if not CSV_IN.exists():
        print(f"FEHLT: {CSV_IN}")
        return 2

    # Lese Vorschläge
    with CSV_IN.open() as f:
        proposals = list(csv.DictReader(f))

    # Sortiere in Aktionen
    to_replace: List[Dict] = []
    to_delete: List[Dict] = []
    skipped: List[Dict] = []

    for p in proposals:
        conf = p["confidence"]
        old = p["original_doi"]
        new = p["proposed_doi"]
        key = (p["molecule"], p["method"], old)

        if conf == "HIGH":
            if new and new.lower() != old.lower():
                to_replace.append(p)
            else:
                skipped.append({**p, "reason": "Crossref-Search lieferte selbe DOI — vermutlich originale war doch ok, kein Change"})
        elif conf == "MEDIUM":
            if MANUAL_ACCEPTS.get(key) and new:
                to_replace.append(p)
            else:
                to_delete.append(p)
        elif conf in ("LOW", "NONE"):
            to_delete.append(p)

    print("Helixar AI — Apply DOI Repairs")
    print()
    print(f"Plan:")
    print(f"  Replace:  {len(to_replace)}")
    print(f"  Delete:   {len(to_delete)}")
    print(f"  Skip:     {len(skipped)}")
    print()

    # Backups
    timestamp = time.strftime("%Y%m%d_%H%M%S")
    print("Backups:")
    for f in FILES_TO_MODIFY:
        if f.exists():
            backup = f.with_suffix(f".bak.{timestamp}")
            shutil.copy2(f, backup)
            print(f"  {f.name} → {backup.name}")
    print()

    # Lese Code-Dateien
    contents: Dict[Path, str] = {}
    for f in FILES_TO_MODIFY:
        if f.exists():
            contents[f] = f.read_text(encoding="utf-8")

    # Wende Replacements an
    total_replace_count = 0
    print("Replace-Operations:")
    for p in to_replace:
        old, new = p["original_doi"], p["proposed_doi"]
        if not old or not new:
            continue
        applied = False
        for f, text in contents.items():
            new_text, n = replace_doi_in_text(text, old, new)
            if n > 0:
                contents[f] = new_text
                total_replace_count += n
                applied = True
        marker = "✅" if applied else "⚠️ nicht gefunden"
        print(f"  {marker} {p['molecule']:<25}/{p['method']:<18} {old} → {new}")

    print()

    # Wende Deletions an
    total_delete_count = 0
    print("Delete-Operations:")
    for p in to_delete:
        old = p["original_doi"]
        if not old:
            continue
        applied = False
        for f, text in contents.items():
            new_text, n = remove_doi_from_text(text, old)
            if n > 0:
                contents[f] = new_text
                total_delete_count += n
                applied = True
        marker = "✅" if applied else "⚠️ nicht gefunden"
        print(f"  {marker} {p['molecule']:<25}/{p['method']:<18} {old}")

    print()

    # Schreibe Dateien zurück
    for f, text in contents.items():
        f.write_text(text, encoding="utf-8")

    print("=" * 60)
    print(f"FERTIG.")
    print(f"  DOIs ersetzt:  {total_replace_count}")
    print(f"  DOIs gelöscht: {total_delete_count}")
    print(f"  Übersprungen:  {len(skipped)}")
    print()
    print(f"Backups gesichert mit Suffix .bak.{timestamp}")
    print(f"Nächster Schritt: python3 validate_sources.py erneut laufen lassen")

    return 0


if __name__ == "__main__":
    sys.exit(main())
