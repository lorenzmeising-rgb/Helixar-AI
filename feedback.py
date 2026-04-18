# ...existing code...
import os
from typing import Dict, Any, Optional, List
from datetime import datetime

import pandas as pd

FEEDBACK_FILENAME = "feedback_history.csv"


def _history_path(provided_path: Optional[str]) -> str:
    if provided_path:
        return provided_path
    base = os.path.dirname(__file__)
    return os.path.join(base, FEEDBACK_FILENAME)


def _clamp01(x: float) -> float:
    return max(0.0, min(1.0, float(x)))


def load_feedback_history(path: Optional[str] = None) -> pd.DataFrame:
    """Load feedback history CSV if it exists, else return empty DataFrame with standard columns."""
    p = _history_path(path)
    cols = [
        "timestamp",
        "microorganism",
        "strain",
        "compound_class",
        "compound_name",
        "outcome",
        "notes",
        "previous_confidence",
        "updated_confidence",
    ]
    if os.path.exists(p):
        try:
            return pd.read_csv(p)
        except Exception:
            return pd.DataFrame(columns=cols)
    return pd.DataFrame(columns=cols)


def _append_history_row(row: Dict[str, Any], path: Optional[str] = None) -> None:
    p = _history_path(path)
    df = load_feedback_history(p)
    df = pd.concat([df, pd.DataFrame([row])], ignore_index=True)
    df.to_csv(p, index=False)


def submit_feedback(
    db,
    selected_strategy: Dict[str, Any],
    outcome: str,
    notes: Optional[str] = None,
    history_path: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Record user feedback and update confidence_score(s) in the provided MicrobialDB instance.

    Parameters:
      - db: MicrobialDB instance
      - selected_strategy: dict identifying the strategy (expects keys like microorganism, strain, compound_name, compound_class)
      - outcome: one of 'success', 'partial', 'failure'
      - notes: optional free-text
      - history_path: optional path to store feedback history CSV

    Behavior:
      - Maps outcomes to simple, transparent confidence adjustments:
          'success' -> +0.10
          'partial' -> +0.03
          'failure' -> -0.20
      - Applies the adjustment to all matching rows in the DB (microorganism+strain+compound where provided).
      - Clamps confidence_score to [0.0, 1.0].
      - Appends one history row per updated DB row with timestamp, previous and updated confidence.
    Returns a concise summary dict.
    """
    outcome_norm = str(outcome).strip().lower()
    if outcome_norm not in ("success", "partial", "failure"):
        return {"ok": False, "message": "invalid outcome; use 'success', 'partial', or 'failure'."}

    deltas = {"success": 0.10, "partial": 0.03, "failure": -0.20}
    delta = float(deltas[outcome_norm])

    # Build a tolerant match mask: require fields present in selected_strategy to match exactly.
    df = db._df  # access internal DataFrame to update confidence; MicrobialDB is pandas-backed
    mask = pd.Series([True] * len(df), index=df.index)
    for key in ("microorganism", "strain", "compound_name", "compound_class"):
        val = selected_strategy.get(key)
        if val is not None:
            mask = mask & (df[key].astype(str) == str(val))

    matched_idx = df[mask].index.tolist()
    if not matched_idx:
        return {"ok": False, "message": "no matching DB entries found", "matched": 0}

    # Ensure numeric confidence column
    df["confidence_score"] = pd.to_numeric(df["confidence_score"], errors="coerce").fillna(0.0)

    updated_rows: List[int] = []
    for idx in matched_idx:
        prev = float(df.at[idx, "confidence_score"])
        new = _clamp01(prev + delta)
        df.at[idx, "confidence_score"] = new
        updated_rows.append(idx)

        history_row = {
            "timestamp": datetime.utcnow().isoformat() + "Z",
            "microorganism": df.at[idx, "microorganism"],
            "strain": df.at[idx, "strain"],
            "compound_class": df.at[idx, "compound_class"],
            "compound_name": df.at[idx, "compound_name"],
            "outcome": outcome_norm,
            "notes": notes or "",
            "previous_confidence": prev,
            "updated_confidence": new,
        }
        _append_history_row(history_row, history_path)

    # Persist updated DataFrame back into the MicrobialDB instance
    db._df = df

    return {
        "ok": True,
        "message": "feedback applied",
        "matched": len(matched_idx),
        "updated_indices": updated_rows,
        "delta": delta,
        "history_file": _history_path(history_path),
    }
# ...existing code...