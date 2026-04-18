# ...existing code...
import pandas as pd
from typing import Optional, Tuple, Dict, Any, List


COLUMNS = [
    "microorganism",
    "strain",
    "compound_class",
    "compound_name",
    "temperature_range",  # stored as tuple (min, max) or None
    "ph_range",           # stored as tuple (min, max) or None
    "oxygen_level",
    "expected_yield",     # numeric (e.g. g/L) or None
    "confidence_score",   # float 0..1
    "source_reference",
]


def _normalize_range(value: Optional[Any]) -> Optional[Tuple[Optional[float], Optional[float]]]:
    if value is None:
        return None
    if isinstance(value, (list, tuple)) and len(value) == 2:
        try:
            return (None if value[0] is None else float(value[0]),
                    None if value[1] is None else float(value[1]))
        except ValueError:
            return None
    # attempt parsing simple string "min-max"
    if isinstance(value, str):
        parts = value.split("-")
        if len(parts) == 2:
            try:
                return (None if parts[0].strip() == "" else float(parts[0]),
                        None if parts[1].strip() == "" else float(parts[1]))
            except ValueError:
                return None
    return None


def _normalize_confidence(value: Any) -> float:
    try:
        v = float(value)
    except Exception:
        v = 0.0
    # clamp 0..1
    return max(0.0, min(1.0, v))


class MicrobialDB:
    """In-memory pandas-backed table for microbial production knowledge."""

    def __init__(self, df: Optional[pd.DataFrame] = None):
        if df is not None:
            # ensure required columns exist (add missing columns with None/defaults)
            missing = [c for c in COLUMNS if c not in df.columns]
            if missing:
                for c in missing:
                    df[c] = None
            self._df = df.copy()
        else:
            self._df = pd.DataFrame(columns=COLUMNS)
            # ensure proper dtypes for sortable columns
            self._df["confidence_score"] = self._df["confidence_score"].astype(float)

    def add_entry(self, entry: Dict[str, Any]) -> None:
        """Add a new entry dict. Keys should match COLUMNS. Missing keys allowed (set to None)."""
        row = {c: entry.get(c, None) for c in COLUMNS}
        row["temperature_range"] = _normalize_range(row["temperature_range"])
        row["ph_range"] = _normalize_range(row["ph_range"])
        row["confidence_score"] = _normalize_confidence(row["confidence_score"])
        # try cast expected_yield to float if provided
        if row["expected_yield"] is not None:
            try:
                row["expected_yield"] = float(row["expected_yield"])
            except Exception:
                row["expected_yield"] = None
        self._df = pd.concat([self._df, pd.DataFrame([row])], ignore_index=True)
        # ensure numeric type for confidence_score column
        self._df["confidence_score"] = pd.to_numeric(self._df["confidence_score"], errors="coerce").fillna(0.0)

    def query_by_compound_class(self, compound_class: str) -> pd.DataFrame:
        """Return rows matching compound_class (case-insensitive)."""
        mask = self._df["compound_class"].astype(str).str.lower() == str(compound_class).lower()
        return self._df[mask].copy()

    def query_by_microorganism(self, microorganism: str) -> pd.DataFrame:
        """Return rows matching microorganism (case-insensitive)."""
        mask = self._df["microorganism"].astype(str).str.lower() == str(microorganism).lower()
        return self._df[mask].copy()

    def top_candidates(self, n: int = 10) -> pd.DataFrame:
        """Return top n rows ranked by confidence_score (descending)."""
        return self._df.sort_values(by="confidence_score", ascending=False).head(n).copy()

    def all(self) -> pd.DataFrame:
        """Return full table copy."""
        return self._df.copy()

    def to_csv(self, path: str) -> None:
        """Save table to CSV. Ranges are serialized as strings."""
        export = self._df.copy()
        export["temperature_range"] = export["temperature_range"].apply(lambda r: "" if r is None else f"{r[0]}-{r[1]}")
        export["ph_range"] = export["ph_range"].apply(lambda r: "" if r is None else f"{r[0]}-{r[1]}")
        export.to_csv(path, index=False)

    @classmethod
    def from_csv(cls, path: str) -> "MicrobialDB":
        df = pd.read_csv(path)
        # parse ranges back into tuples where possible
        def parse_range_cell(x):
            if pd.isna(x) or x == "":
                return None
            s = str(x).strip()
            parts = s.split("-")
            if len(parts) == 2:
                try:
                    return (None if parts[0].strip() == "" else float(parts[0]),
                            None if parts[1].strip() == "" else float(parts[1]))
                except Exception:
                    return s  # behalte den beschreibenden String, wenn er nicht als Zahlenpaar parst
            return s  # beschreibende Werte (z.B. "moderat", "breit") unverändert lassen

        # Ensure expected legacy columns exist; if missing, add with None so downstream
        # consumers don't fail when CSV is sparse/minimal.
        for c in COLUMNS:
            if c not in df.columns:
                df[c] = None

        if "temperature_range" in df.columns:
            df["temperature_range"] = df["temperature_range"].apply(parse_range_cell)
        if "ph_range" in df.columns:
            df["ph_range"] = df["ph_range"].apply(parse_range_cell)
        if "confidence_score" in df.columns:
            df["confidence_score"] = pd.to_numeric(df["confidence_score"], errors="coerce").fillna(0.0)
        return cls(df)
# ...existing code...