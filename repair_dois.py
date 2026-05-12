"""
repair_dois.py
==============
Findet für jeden FAILED-DOI aus dem Validator-Report die wahrscheinlich
RICHTIGE DOI via Crossref-Search-API.

Workflow:
  1. Lese docs/sources_validation.csv (vom Validator erzeugt)
  2. Für jeden FAILED-Eintrag: extrahiere Autor + Titel + Jahr aus
     dem zitierten Text
  3. Crossref-Search-API → bestes Match abrufen
  4. Score: wie gut passt das Match zu unserer Zitierung?
  5. 3 Buckets:
     HIGH      — sehr wahrscheinlich richtig (Title-Match > 80 %, Autor stimmt)
     MEDIUM    — plausibel, aber manuelle Review nötig
     LOW       — kein gutes Match, DOI wird entfernt

Output:
  docs/doi_repair_proposals.csv  — Vorschläge mit Confidence
  docs/doi_repair_summary.md     — menschlich lesbar
"""

from __future__ import annotations

import csv
import json
import os
import re
import sys
import time
import urllib.parse
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

REPO = Path(__file__).resolve().parent
os.chdir(REPO)

DOCS_DIR = REPO / "docs"
CACHE_DIR = REPO / ".cache"

USER_AGENT = "Helixar-AI-DOI-Repair/1.0 (mailto:lorenzmeising@icloud.com)"
HEADERS = {"User-Agent": USER_AGENT}

CROSSREF_SEARCH = "https://api.crossref.org/works"


# ---------------------------------------------------------------------------
# Cache (gemeinsam mit Validator nutzbar)
# ---------------------------------------------------------------------------

class Cache:
    def __init__(self, path: Path):
        self.path = path
        self.data: Dict[str, Any] = {}
        if path.exists():
            try:
                with path.open() as f:
                    self.data = json.load(f)
            except Exception:
                self.data = {}

    def get(self, key): return self.data.get(key)

    def set(self, key, value):
        self.data[key] = value
        try:
            with self.path.open("w") as f:
                json.dump(self.data, f, ensure_ascii=False, indent=2)
        except Exception:
            pass


CACHE = Cache(CACHE_DIR / "api_responses.json")


def _http_get_json(url: str, params: Dict[str, Any]) -> Optional[Dict]:
    import requests
    cache_key = f"GET::{url}::{json.dumps(params, sort_keys=True)}"
    cached = CACHE.get(cache_key)
    if cached is not None:
        return cached
    for attempt in range(3):
        try:
            r = requests.get(url, params=params, headers=HEADERS, timeout=20)
            if r.status_code == 200:
                data = r.json()
                CACHE.set(cache_key, data)
                return data
            time.sleep(1.0 + attempt)
        except Exception:
            time.sleep(1.0 + attempt)
    return None


# ---------------------------------------------------------------------------
# Zitierungs-Text parsen
# ---------------------------------------------------------------------------

@dataclass
class CitationContext:
    """Was wir aus dem cited_text um den FAILED-DOI extrahieren."""
    failed_doi: str
    author_surname: str = ""
    year: Optional[int] = None
    journal_hint: str = ""
    title_keywords: str = ""
    raw_segment: str = ""        # Der Text-Block der zur DOI gehört


YEAR_RE = re.compile(r"\b(19|20)\d{2}\b")
# "Hansen et al.," / "Smith & Jones," / "Levine,"
AUTHOR_RE = re.compile(r"^([A-Z][a-zA-Zäöüß\-']+)(?:\s+et\s+al\.?|,)", re.IGNORECASE)


def _extract_segment_around_doi(cited_text: str, doi: str) -> str:
    """Findet den Text-Abschnitt der die DOI zitiert (zwischen ';' Trennern)."""
    # Häufige Trenner zwischen Quellen: "; " oder ".)"
    # Splitte konservativ
    segments = re.split(r";\s+(?=[A-Z])", cited_text)
    for seg in segments:
        if doi in seg:
            return seg.strip()
    return cited_text


def parse_citation_context(failed_doi: str, cited_text: str) -> CitationContext:
    """Aus dem cited_text das Segment um den FAILED-DOI extrahieren + parsen."""
    ctx = CitationContext(failed_doi=failed_doi)
    seg = _extract_segment_around_doi(cited_text, failed_doi)
    ctx.raw_segment = seg

    # Author: erstes Wort
    m = AUTHOR_RE.match(seg)
    if m:
        ctx.author_surname = m.group(1)

    # Year
    ym = YEAR_RE.search(seg)
    if ym:
        try:
            ctx.year = int(ym.group(0))
        except Exception:
            pass

    # Journal: zwischen Author und Year liegt meist das Journal
    if ctx.author_surname and ctx.year:
        between = seg
        # Schneide vor Year
        idx_year = seg.find(str(ctx.year))
        if idx_year > 0:
            between = seg[:idx_year]
        # Schneide nach Author
        idx_auth_end = between.find(",")
        if idx_auth_end > 0:
            between = between[idx_auth_end + 1:]
        ctx.journal_hint = between.strip().rstrip(",").strip()

    # Title: nach Year, vor DOI oder vor "(doi:" / "—"
    after_year = seg
    if ctx.year:
        i = seg.find(str(ctx.year))
        if i >= 0:
            after_year = seg[i + 4:]
    # Schneide vor DOI oder typischen Markers
    title_part = after_year
    for marker in ["(doi:", " doi:", "(http", "—"]:
        idx = title_part.find(marker)
        if idx >= 0:
            title_part = title_part[:idx]
            break
    # Reinige: führende "—" / ":" / "—" entfernen
    title_part = title_part.strip(" —-:,;").strip()
    ctx.title_keywords = title_part[:120]   # max 120 chars
    return ctx


# ---------------------------------------------------------------------------
# Crossref-Search
# ---------------------------------------------------------------------------

def search_crossref(ctx: CitationContext) -> List[Dict]:
    """Crossref-Search-API mit Autor + Titel-Keywords + Jahr."""
    params: Dict[str, Any] = {"rows": 5}
    if ctx.author_surname:
        params["query.author"] = ctx.author_surname
    if ctx.title_keywords:
        params["query.title"] = ctx.title_keywords
    elif ctx.journal_hint:
        params["query.container-title"] = ctx.journal_hint
    # Year-Filter (Crossref unterstützt from-pub-date / until-pub-date)
    if ctx.year:
        params["filter"] = f"from-pub-date:{ctx.year-1},until-pub-date:{ctx.year+1}"
    data = _http_get_json(CROSSREF_SEARCH, params=params)
    if not data:
        return []
    return (data.get("message") or {}).get("items") or []


# ---------------------------------------------------------------------------
# Match-Scoring: wie gut passt ein Crossref-Result zu unserer Zitierung?
# ---------------------------------------------------------------------------

def _normalize_for_compare(s: str) -> str:
    return re.sub(r"[^a-z0-9]+", " ", s.lower()).strip()


def _word_overlap(a: str, b: str) -> float:
    """Anteil gemeinsamer Wörter (>=4 Buchstaben), normalisiert."""
    a_words = set(w for w in _normalize_for_compare(a).split() if len(w) >= 4)
    b_words = set(w for w in _normalize_for_compare(b).split() if len(w) >= 4)
    if not a_words:
        return 0.0
    return len(a_words & b_words) / len(a_words)


@dataclass
class RepairProposal:
    molecule: str
    method: str
    original_doi: str
    proposed_doi: str
    confidence: str             # HIGH | MEDIUM | LOW | NONE
    score: int                  # 0-100
    proposed_title: str
    proposed_authors: str
    proposed_year: Optional[int]
    proposed_journal: str
    notes: str
    raw_citation_segment: str


def score_candidate(cand: Dict, ctx: CitationContext) -> Tuple[int, str]:
    """Vergibt einen Score 0-100 + Begründung wie gut ein Kandidat passt."""
    score = 0
    notes = []

    # Autor-Match (+30)
    cand_authors = cand.get("author") or []
    if cand_authors:
        first_family = (cand_authors[0].get("family") or "").lower()
        if first_family and ctx.author_surname.lower() in first_family:
            score += 30
            notes.append("Autor-Match")
        elif first_family and ctx.author_surname.lower() == first_family.split()[0]:
            score += 30
            notes.append("Autor-Match")
        else:
            notes.append(f"Autor-Mismatch (gefunden: {first_family})")

    # Jahr-Match (+20)
    cand_year = None
    issued = cand.get("issued") or {}
    dp = (issued.get("date-parts") or [[None]])[0]
    if dp and dp[0]:
        try:
            cand_year = int(dp[0])
        except Exception:
            pass
    if cand_year and ctx.year:
        if cand_year == ctx.year:
            score += 20
            notes.append("Jahr-Match")
        elif abs(cand_year - ctx.year) == 1:
            score += 10
            notes.append(f"Jahr off-by-one ({cand_year} vs {ctx.year})")
        else:
            notes.append(f"Jahr-Mismatch ({cand_year} vs {ctx.year})")

    # Titel-Match (+40, höchstes Gewicht)
    cand_titles = cand.get("title") or []
    cand_title = cand_titles[0] if cand_titles else ""
    overlap = _word_overlap(ctx.title_keywords, cand_title)
    title_score = int(40 * overlap)
    score += title_score
    if overlap >= 0.6:
        notes.append(f"Titel-Match {int(overlap*100)} %")
    elif overlap >= 0.3:
        notes.append(f"Titel-Teilmatch {int(overlap*100)} %")
    else:
        notes.append(f"Titel-Mismatch {int(overlap*100)} %")

    # Journal-Match (+10)
    cand_container = (cand.get("container-title") or [""])[0]
    if ctx.journal_hint and cand_container:
        if _word_overlap(ctx.journal_hint, cand_container) >= 0.5:
            score += 10
            notes.append("Journal-Match")

    return score, " | ".join(notes)


def classify_confidence(score: int) -> str:
    if score >= 70:
        return "HIGH"
    if score >= 45:
        return "MEDIUM"
    if score >= 25:
        return "LOW"
    return "NONE"


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    csv_in = DOCS_DIR / "sources_validation.csv"
    if not csv_in.exists():
        print(f"FEHLT: {csv_in}. Erst validate_sources.py laufen lassen.")
        return 2

    print("Helixar AI — DOI Repair via Crossref-Search")
    print()

    # Lese FAILED + auch PARTIAL (deren Crossref-Match < 2/3)
    targets: List[Tuple[str, str, str, str]] = []
    with csv_in.open() as f:
        for row in csv.DictReader(f):
            status = row.get("overall_status", "")
            match_score = int(row.get("crossref_match_score", "0") or "0")
            # FAILED → immer; PARTIAL nur wenn match_score < 2
            if status == "FAILED" or (status == "PARTIAL" and match_score < 2):
                targets.append((
                    row["molecule"], row["method"], row["doi"], row["cited_text"]
                ))

    print(f"Reparatur-Kandidaten: {len(targets)} Einträge")
    print()

    proposals: List[RepairProposal] = []
    for i, (mol, meth, doi, cited) in enumerate(targets, 1):
        ctx = parse_citation_context(doi, cited)
        cands = search_crossref(ctx)
        best: Optional[Dict] = None
        best_score = 0
        best_notes = ""
        for c in cands:
            s, n = score_candidate(c, ctx)
            if s > best_score:
                best, best_score, best_notes = c, s, n

        if best:
            p_doi = best.get("DOI") or ""
            p_titles = best.get("title") or []
            p_authors = best.get("author") or []
            author_names = " · ".join(
                f"{a.get('family','')}{', '+a.get('given','') if a.get('given') else ''}"
                for a in p_authors[:5]
            )
            p_year = None
            issued = best.get("issued") or {}
            dp = (issued.get("date-parts") or [[None]])[0]
            if dp and dp[0]:
                try:
                    p_year = int(dp[0])
                except Exception:
                    pass
            p_journal = (best.get("container-title") or [""])[0]
        else:
            p_doi = ""
            p_titles = [""]
            author_names = ""
            p_year = None
            p_journal = ""
            best_notes = "Crossref-Search lieferte keine Kandidaten"

        confidence = classify_confidence(best_score)
        proposals.append(RepairProposal(
            molecule=mol,
            method=meth,
            original_doi=doi,
            proposed_doi=p_doi,
            confidence=confidence,
            score=best_score,
            proposed_title=p_titles[0] if p_titles else "",
            proposed_authors=author_names,
            proposed_year=p_year,
            proposed_journal=p_journal,
            notes=best_notes,
            raw_citation_segment=ctx.raw_segment[:200],
        ))

        emoji = {"HIGH": "🟢", "MEDIUM": "🟡", "LOW": "🟠", "NONE": "🔴"}[confidence]
        print(f"  [{i:3d}/{len(targets)}] {emoji} {mol[:22]:<22}/{meth[:18]:<18} "
              f"{doi[:35]:<35} → score={best_score} ({confidence})")
        time.sleep(0.1)

    # Output
    csv_out = DOCS_DIR / "doi_repair_proposals.csv"
    md_out = DOCS_DIR / "doi_repair_summary.md"

    if proposals:
        with csv_out.open("w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(asdict(proposals[0]).keys()))
            w.writeheader()
            for p in proposals:
                w.writerow(asdict(p))

    # Buckets
    by_conf = {"HIGH": [], "MEDIUM": [], "LOW": [], "NONE": []}
    for p in proposals:
        by_conf[p.confidence].append(p)

    lines = [
        "# DOI-Reparatur-Vorschläge",
        "",
        f"**Gesamt:** {len(proposals)} Kandidaten geprüft",
        "",
        f"- 🟢 **HIGH** (auto-anwenden empfohlen): {len(by_conf['HIGH'])}",
        f"- 🟡 **MEDIUM** (manuelle Review): {len(by_conf['MEDIUM'])}",
        f"- 🟠 **LOW** (zweifelhaft, DOI eher löschen): {len(by_conf['LOW'])}",
        f"- 🔴 **NONE** (kein Crossref-Match — DOI löschen): {len(by_conf['NONE'])}",
        "",
        "---",
        "",
        "## 🟢 HIGH-Confidence — auto-anwenden",
        "",
    ]
    for p in by_conf["HIGH"]:
        lines.append(f"### {p.molecule}/{p.method} (score {p.score})")
        lines.append(f"- Original DOI: `{p.original_doi}`")
        lines.append(f"- Vorgeschlagen: `{p.proposed_doi}`")
        lines.append(f"- Title: {p.proposed_title}")
        lines.append(f"- Autoren: {p.proposed_authors}")
        lines.append(f"- Jahr/Journal: {p.proposed_year} · {p.proposed_journal}")
        lines.append(f"- Notes: {p.notes}")
        lines.append("")

    lines.extend([
        "---",
        "",
        "## 🟡 MEDIUM-Confidence — kurz prüfen",
        "",
    ])
    for p in by_conf["MEDIUM"]:
        lines.append(f"### {p.molecule}/{p.method} (score {p.score})")
        lines.append(f"- Original DOI: `{p.original_doi}`")
        lines.append(f"- Vorgeschlagen: `{p.proposed_doi}`")
        lines.append(f"- Title: {p.proposed_title}")
        lines.append(f"- Autoren: {p.proposed_authors}")
        lines.append(f"- Jahr/Journal: {p.proposed_year} · {p.proposed_journal}")
        lines.append(f"- Notes: {p.notes}")
        lines.append(f"- Original-Zitat: {p.raw_citation_segment}")
        lines.append("")

    lines.extend([
        "---",
        "",
        "## 🟠 LOW + 🔴 NONE — DOI löschen, nur Klartext-Zitierung behalten",
        "",
    ])
    for p in by_conf["LOW"] + by_conf["NONE"]:
        lines.append(f"- **{p.molecule}/{p.method}**: `{p.original_doi}` — {p.notes}")

    md_out.write_text("\n".join(lines), encoding="utf-8")

    print()
    print("=" * 60)
    print(f"FERTIG. {len(proposals)} Vorschläge:")
    print(f"  🟢 HIGH:   {len(by_conf['HIGH']):3d}")
    print(f"  🟡 MEDIUM: {len(by_conf['MEDIUM']):3d}")
    print(f"  🟠 LOW:    {len(by_conf['LOW']):3d}")
    print(f"  🔴 NONE:   {len(by_conf['NONE']):3d}")
    print()
    print(f"Reports:")
    print(f"  {csv_out}")
    print(f"  {md_out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
