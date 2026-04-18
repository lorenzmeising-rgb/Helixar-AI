import math
from typing import Dict, Any, List, Optional, Union, Tuple

import pandas as pd
import logging

from database import MicrobialDB
from molecules_db import get_smiles_for


class DecisionEngine:
    """Rule-based decision engine for recommending microbial production strategies.

    Usage:
      engine = DecisionEngine(db=MicrobialDB.from_csv("db.csv"))
      recommendations = engine.recommend("antibiotic", {"time": 12, "budget": 50000, "risk_tolerance": "medium"})
    """

    def __init__(self, db: Optional[MicrobialDB] = None, db_path: Optional[str] = None):
        if db is not None:
            self.db = db
        elif db_path is not None:
            self.db = MicrobialDB.from_csv(db_path)
        else:
            raise ValueError("Provide either a MicrobialDB instance or a db_path to a CSV.")

    @staticmethod
    def _parse_risk_tolerance(rt: Optional[Union[str, float, int]]) -> float:
        """Normalize risk tolerance to float in [0,1].
        0.0 = very risk averse (favor confidence), 1.0 = very risk seeking (favor yield).
        Accepts 'low'|'medium'|'high' or numeric 0..1.
        """
        if rt is None:
            return 0.5
        if isinstance(rt, (float, int)):
            return float(max(0.0, min(1.0, float(rt))))
        s = str(rt).strip().lower()
        if s in ("low", "risk-averse", "risk_averse"):
            return 0.0
        if s in ("high", "risk-seeking", "risk_seeking"):
            return 1.0
        return 0.5

    @staticmethod
    def _normalize_series_for_scoring(s: pd.Series) -> pd.Series:
        """Scale numeric series to 0..1 for scoring. Missing or constant series map to 0.0..0.0 => 0.5."""
        numeric = pd.to_numeric(s, errors="coerce").fillna(float("nan"))
        if numeric.dropna().empty:
            return pd.Series([0.0] * len(s), index=s.index)
        mn = numeric.min()
        mx = numeric.max()
        if math.isclose(mn, mx):
            return pd.Series([0.5] * len(s), index=s.index)
        return ((numeric - mn) / (mx - mn)).fillna(0.0)

    def _apply_filters(self, df: pd.DataFrame, target_constraints: Dict[str, Any]) -> pd.DataFrame:
        """Apply optional pre-scoring filters based on target_constraints.

        - If scale == 'industrial' and 'scalability_score' exists, keep rows with scalability_score > 0.5
        - If sustainability == 'high' and 'type' exists, prefer rows where type contains 'yeast' or 'biolog'

        If filters remove all rows, returns the original df (fallback).
        """
        if df.empty:
            return df

        filtered = df

        # industrial scale filter
        try:
            scale_req = str(target_constraints.get("scale", "")).strip().lower()
        except Exception:
            scale_req = ""
        if scale_req == "industrial" and "scalability_score" in df.columns:
            try:
                keep = pd.to_numeric(df["scalability_score"], errors="coerce").fillna(0.0) > 0.5
                filtered = filtered.loc[keep]
            except Exception:
                pass

        # sustainability preference (prefer biological / yeast processes)
        try:
            sustain_req = str(target_constraints.get("sustainability", "")).strip().lower()
        except Exception:
            sustain_req = ""
        if sustain_req == "high":
            # attempt to match biological/process_type fields
            process_cols = [c for c in ("process_type", "type", "method") if c in df.columns]
            matched = pd.Series([False] * len(filtered), index=filtered.index)
            try:
                for c in process_cols:
                    matched = matched | filtered[c].astype(str).str.contains(r"biolog|yeast|cellul|ferment", case=False, na=False)
                preferred = filtered.loc[matched]
                if not preferred.empty:
                    filtered = preferred
            except Exception:
                pass

        # preferred_method filter (optional)
        try:
            pref = str(target_constraints.get("preferred_method", "")).strip().lower()
        except Exception:
            pref = ""
        if pref and pref != "none":
            method_cols = [c for c in ("method", "process_method", "process_type") if c in df.columns]
            if method_cols:
                try:
                    mask = pd.Series([False] * len(filtered), index=filtered.index)
                    for c in method_cols:
                        mask = mask | filtered[c].astype(str).str.contains(pref, case=False, na=False)
                    filtered_pref = filtered.loc[mask]
                    if not filtered_pref.empty:
                        filtered = filtered_pref
                except Exception:
                    pass

        if filtered.empty:
            return df
        return filtered

    @staticmethod
    def _extract_smiles_properties(smiles: Optional[str]) -> Dict[str, str]:
        """Light-weight SMILES property extractor (no external libs).

        Infers simple, explainable properties relevant to process optimization:
          - purification_difficulty (low|medium|high)
          - stability (low|medium|high)
          - toxicity (low|medium|high)
          - complexity (low|medium|high)
          - raw_material_complexity (low|medium|high)
          - scalability (low|medium|high)

        Uses straightforward string-pattern heuristics; intentionally simple
        and deterministic for traceability.
        """
        # defaults
        props: Dict[str, Any] = {
            "purification_difficulty": "medium",
            "stability": "medium",
            "toxicity": "medium",
            "complexity": "medium",
            "raw_material_complexity": "medium",
            "scalability": "medium",
            "base_toxicity": "medium",
            # enhanced structural signals
            "aromatic_rings": 0,
            "long_aliphatic_chain": False,
            "polar_group_count": 0,
            "ester_present": False,
            "amide_present": False,
            "crystallization_potential": False,
        }

        if not smiles or not isinstance(smiles, str) or smiles.strip() == "":
            return props

        s = smiles.strip()
        try:
            # aromaticity proxy: lowercase 'c1' or 'c2' patterns indicate aromatic ring presence
            arom = s.count("c1") + s.count("c2") + s.count("c3")
            if arom > 0 or ("ar" in s.lower()):
                props["aromatic_rings"] = arom if arom > 0 else 1
                # aromatic systems often confer higher chemical stability
                props["stability"] = "high"
                # aromatic compounds often crystallize well
                props["crystallization_potential"] = True

            # long aliphatic chains -> hydrophobicity and purification difficulty
            if "CCCC" in s or s.count("C") >= 12:
                props["long_aliphatic_chain"] = True
                props["complexity"] = "high"
                props["purification_difficulty"] = "high"

            # polar groups: O, N, OH, NH, COOH increase solubility but may require different separations
            polar = 0
            polar += s.count("O")
            polar += s.count("N")
            # count explicit OH/NH patterns
            if "OH" in s or "NH" in s:
                polar += 1
            props["polar_group_count"] = polar
            if polar >= 2:
                # more polar groups usually improve aqueous solubility
                props["purification_difficulty"] = props.get("purification_difficulty") or "medium"

            # ester / amide presence increases hydrolysis risk
            if "C(=O)O" in s or "COC(=O)" in s or "OC(=O)" in s:
                props["ester_present"] = True
            if "C(=O)N" in s or "NC(=O)" in s or "C(=O)N(" in s:
                props["amide_present"] = True

            # simple toxicity heuristics
            if "Cl" in s or "Br" in s or "I" in s:
                props["toxicity"] = "high"
            elif "C=O" in s and props.get("polar_group_count", 0) < 2:
                props["toxicity"] = "medium"

            # raw material complexity linked to structural complexity
            if props.get("complexity") == "high" or props.get("polar_group_count", 0) > 4:
                props["raw_material_complexity"] = "high"

            # scalability heuristics: low if instability or high complexity
            if props.get("stability") == "low" or props.get("complexity") == "high":
                props["scalability"] = "low"

            # final base toxicity
            props["base_toxicity"] = props.get("toxicity", "medium")
        except Exception:
            pass

        return props

    def _compute_score_breakdown(
        self,
        p: Dict[str, Any],
        properties: Dict[str, Any],
        smiles_candidate: Optional[str],
        molecule_type: str,
    ) -> Tuple[Dict[str, float], Dict[str, float], str]:
        """Compute a reproducible score breakdown (cost and risk contributions) and a short weighting explanation.

        This is a deterministic, explainable summary used for human-readable explanations and does not replace
        the main numeric scoring which must be preserved.
        """
        try:
            # weights
            try:
                w_smiles = self.__class__._get_smiles_weight_static(molecule_type)
            except Exception:
                # fallback
                w_smiles = 0.5
            w_process = 1.0
            w_economic = 1.2
            w_downstream = 1.1
            w_biophysical = 1.3 if molecule_type in ("protein", "peptide") else 0.0

            # inputs
            steps = int(p.get("number_of_steps") or 0)
            method = str(p.get("method") or "").lower()
            desired_purity = str(p.get("desired_purity") or "standard").lower()
            scale = str(p.get("scale") or "lab").lower()
            raw_cost = str(p.get("raw_material_cost") or "medium").lower()
            raw_avail = str(p.get("raw_material_availability") or "medium").lower()
            strict_waste = bool(p.get("strict_waste_constraints"))
            stype = (p.get("molecule_subtype") or "").lower()

            complexity_map = {"low": 1, "medium": 2, "high": 3, "very high": 4}
            stability_map = {"high": 1, "medium": 2, "low": 3}
            comp_num = complexity_map.get(properties.get("complexity", "medium"), 2)
            instability_num = stability_map.get(properties.get("stability", "medium"), 2)

            scoreBreakdown = {"smiles": 0.0, "process": 0.0, "economic": 0.0, "downstream": 0.0, "biophysical": 0.0}
            riskBreakdown = {"smiles": 0.0, "process": 0.0, "economic": 0.0, "downstream": 0.0, "biophysical": 0.0}

            # smiles
            delta_smiles_cost = float(comp_num * w_smiles)
            delta_smiles_risk = float(instability_num * w_smiles)
            if properties.get("base_toxicity") == "high" or properties.get("toxicity") == "high":
                delta_smiles_cost += 2 * w_smiles
                delta_smiles_risk += 2 * w_smiles
            if properties.get("complexity") == "high":
                delta_smiles_cost += 2 * w_smiles
            if smiles_candidate and isinstance(smiles_candidate, str) and len(smiles_candidate) > 100:
                delta_smiles_cost += 2 * w_smiles
                delta_smiles_risk += 2 * w_smiles
            scoreBreakdown["smiles"] = round(delta_smiles_cost, 3)
            riskBreakdown["smiles"] = round(delta_smiles_risk, 3)

            # process
            delta_process_cost = float(steps * 1.2 * w_process)
            delta_process_risk = 0.0
            if scale == "industrial":
                delta_process_cost += 3 * w_process
                delta_process_risk += 3 * w_process
            if method in ("chemical", "chem", "chemical synthesis"):
                delta_process_cost += 2 * w_process
            if method in ("biotechnological", "biotech", "biotechnological synthesis"):
                delta_process_risk += 2 * w_process
            if method in ("extraction", "extract", "extraction-based"):
                delta_process_risk += 1 * w_process
            scoreBreakdown["process"] = round(delta_process_cost, 3)
            riskBreakdown["process"] = round(delta_process_risk, 3)

            # economic
            delta_econ_cost = 0.0
            delta_econ_risk = 0.0
            if raw_cost == "high":
                delta_econ_cost += 4 * w_economic
            if raw_avail == "low":
                delta_econ_risk += 3 * w_economic
            if strict_waste:
                delta_econ_cost += 2 * w_economic
                delta_econ_risk += 1 * w_economic
            scoreBreakdown["economic"] = round(delta_econ_cost, 3)
            riskBreakdown["economic"] = round(delta_econ_risk, 3)

            # downstream
            delta_down_cost = 0.0
            delta_down_risk = 0.0
            if properties.get("purification_difficulty") in ("high", "very high"):
                delta_down_cost += 3 * w_downstream
            if desired_purity in (">99%", "very high"):
                delta_down_cost += 3 * w_downstream
                delta_down_risk += 2 * w_downstream
            if stype == "antibody":
                delta_down_cost += 4
                delta_down_risk += 3
            scoreBreakdown["downstream"] = round(delta_down_cost, 3)
            riskBreakdown["downstream"] = round(delta_down_risk, 3)

            # biophysical
            delta_bio_cost = 0.0
            delta_bio_risk = 0.0
            if w_biophysical > 0:
                agg = str(p.get("aggregation_risk") or "").lower()
                fold = str(p.get("folding_complexity") or "").lower()
                bst = str(p.get("biophysical_stability") or "").lower()
                if agg == "high":
                    delta_bio_risk += 4 * w_biophysical
                    delta_bio_cost += 2 * w_biophysical
                if fold == "high":
                    delta_bio_risk += 3 * w_biophysical
                if bst == "low":
                    delta_bio_risk += 4 * w_biophysical
                    delta_bio_cost += 2 * w_biophysical
            scoreBreakdown["biophysical"] = round(delta_bio_cost, 3)
            riskBreakdown["biophysical"] = round(delta_bio_risk, 3)

            weighting_explanation = (
                f"weights: smiles={w_smiles}, process={w_process}, economic={w_economic}, downstream={w_downstream}, biophysical={w_biophysical}"
            )

            return scoreBreakdown, riskBreakdown, weighting_explanation
        except Exception:
            return {"smiles": 0.0, "process": 0.0, "economic": 0.0, "downstream": 0.0, "biophysical": 0.0}, {
                "smiles": 0.0,
                "process": 0.0,
                "economic": 0.0,
                "downstream": 0.0,
                "biophysical": 0.0,
            }, "weights unavailable"

    @staticmethod
    def _get_smiles_weight_static(mtype: str) -> float:
        try:
            key = str(mtype).lower()
        except Exception:
            return 0.5
        if key == "small_molecule":
            return 1.0
        if key == "peptide":
            return 0.6
        if key == "protein":
            return 0.3
        if key == "natural_product":
            return 0.8
        return 0.5

    def classifyProcess(self, inp: Dict[str, Any], properties: Dict[str, Any]) -> str:
        """Lightweight process-class classifier.

        Deterministically assigns a process_class string based on inputs and SMILES-derived properties.
        This is intentionally conservative and only used to apply very small score adjustments.
        """
        try:
            mt = str(inp.get("molecule_type") or properties.get("molecule_type") or "").lower()
        except Exception:
            mt = ""
        try:
            subtype = str(inp.get("molecule_subtype") or "").lower()
        except Exception:
            subtype = ""
        try:
            pur = str(properties.get("purification_difficulty") or "").lower()
        except Exception:
            pur = ""
        try:
            comp = str(properties.get("complexity") or "").lower()
        except Exception:
            comp = ""

        # Small molecule branches
        if mt == "small_molecule":
            if subtype == "volatile":
                return "volatile_small_molecule"
            # commodity: low complexity, easy purification
            if comp in ("low",) and pur in ("low", "medium"):
                return "commodity_small_molecule"
            if pur in ("high", "very high"):
                return "difficult_purification_small_molecule"
            if comp == "high":
                return "complex_small_molecule"
            return "small_molecule_standard"

        # Peptides
        if mt == "peptide":
            if subtype == "cyclic":
                return "cyclic_peptide"
            return "linear_peptide"

        # Proteins / biologics
        if mt == "protein":
            if subtype == "antibody":
                return "antibody"
            if subtype == "enzyme":
                return "enzyme"
            if comp == "high" or pur in ("high", "very high"):
                return "high_complexity_biologic"
            return "protein_standard"

        # Natural products
        if mt == "natural_product":
            if pur in ("high", "very high") or str(properties.get("downstream_complexity") or "").lower() in ("very high",):
                return "natural_product_complex"
            return "natural_product_standard"

        return "unknown"
        
        # If still empty, return empty result
        if candidates is None or (hasattr(candidates, 'empty') and candidates.empty):
            try:
                print("DB CONTENT: none or empty")
                print("--- DEBUG END ---\n")
            except Exception:
                pass
            return []

        # Apply optional pre-scoring filters (non-destructive). This may narrow down candidates
        # according to constraints like scale or sustainability preference.
        filtered = self._apply_filters(candidates, target_constraints)

        # Critical safety check: if filters removed all candidates, restore full DB copy
        try:
            if filtered is None or (hasattr(filtered, 'empty') and filtered.empty):
                print("ERROR: Candidates empty after filtering - restoring original dataset")
                candidates = self.db.all()
                filtered = candidates.copy()
        except Exception:
            pass

        # (target molecule matching already applied above). No-op here to avoid double-filtering.

        # annotate rows with input infrastructure for downstream reasons/compat checks
        infra = target_constraints.get("infrastructure") if isinstance(target_constraints, dict) else None
        try:
            if infra is not None:
                # ensure list-like
                filtered = filtered.copy()
                filtered["_input_infrastructure"] = [infra] * len(filtered)
            else:
                filtered = filtered.copy()
                filtered["_input_infrastructure"] = [None] * len(filtered)
        except Exception:
            filtered = filtered.copy()
            filtered["_input_infrastructure"] = [None] * len(filtered)

        # chemical infrastructure penalty if chemical_plant not available
        try:
            if isinstance(infra, (list, tuple)):
                if "chemical_plant" not in [i.lower() for i in infra]:
                    # mark penalty for chemical-method rows
                    mask = pd.Series([False] * len(filtered), index=filtered.index)
                    for col in ("method", "process_type", "type"):
                        if col in filtered.columns:
                            mask = mask | filtered[col].astype(str).str.contains(r"chemical|chem", case=False, na=False)
                    filtered["_chemical_penalty"] = mask.astype(float) * 0.05
                else:
                    filtered["_chemical_penalty"] = [0.0] * len(filtered)
            else:
                filtered["_chemical_penalty"] = [0.0] * len(filtered)
        except Exception:
            filtered["_chemical_penalty"] = [0.0] * len(filtered)

        # Final debug prints before scoring
        try:
            print("FINAL CANDIDATES:", filtered["compound_name"].astype(str).tolist())
            print("--- DEBUG END ---\n")
        except Exception:
            pass

        scores = self._build_score(filtered, risk, target_constraints)
        filtered = filtered.copy()
        filtered["decision_score"] = scores

        # Constraint-driven small nudges (deterministic, lightweight). These adjustments
        # are intentionally small so they don't break existing behavior but let user
        # constraints influence ranking.
        try:
            scale_flag = str(target_constraints.get("scale", "")).strip().lower()
        except Exception:
            scale_flag = ""
        try:
            sustain_flag = str(target_constraints.get("sustainability", "")).strip().lower()
        except Exception:
            sustain_flag = ""
        try:
            pref_method = str(target_constraints.get("preferred_method", "")).strip().lower()
        except Exception:
            pref_method = ""

        # normalized helper series (safe)
        norm_scal = self._normalize_series_for_scoring(filtered.get("scalability_score", pd.Series([0] * len(filtered))))
        norm_sust = self._normalize_series_for_scoring(filtered.get("sustainability_score", pd.Series([0] * len(filtered))))
        norm_cost = self._normalize_series_for_scoring(filtered.get("cost_score", pd.Series([0] * len(filtered))))

        # small additive nudges
        if sustain_flag == "high":
            # add up to +0.06 proportional to sustainability score
            filtered["decision_score"] = (filtered["decision_score"] + 0.06 * norm_sust).clip(0.0, 1.0)
        if scale_flag == "industrial":
            # add up to +0.06 proportional to scalability
            filtered["decision_score"] = (filtered["decision_score"] + 0.06 * norm_scal).clip(0.0, 1.0)

        # preferred_method scoring bonus moved into _build_score to centralize logic

        ranked = filtered.sort_values(by="decision_score", ascending=False).head(top_n)

        notes = []
        if "time" in target_constraints or "budget" in target_constraints:
            notes.append("Time/budget constraints provided but no explicit time/cost in DB; recommendations favor higher score.")
        notes.append(f"Risk tolerance normalized to {risk:.2f}; scoring weights => yield-focused: {0.4 + 0.4 * risk:.2f}, confidence-focused: {1 - (0.4 + 0.4 * risk):.2f}")

        results: List[Dict[str, Any]] = []
        for rank, (_, row) in enumerate(ranked.iterrows(), start=1):
            reasons = self._build_reasons(row)
            # include a copy of the input constraints/context for downstream consumers
            input_ctx = {k: target_constraints.get(k) for k in ("target_molecule", "application", "scale", "purity", "sustainability", "preferred_method", "infrastructure") if k in target_constraints}

            # Map new molecule-based DB structure to legacy fields for backward compatibility
            row_dict = row.to_dict()
            try:
                method_val = row.get("method") or row.get("process_type")
            except Exception:
                method_val = None
            try:
                process_type_val = row.get("process_type")
            except Exception:
                process_type_val = None

            # Friendly labels (optional, non-destructive)
            method_label_map = {
                "chemical": "Chemical Process",
                "biotech": "Biotechnological Process",
                "extraction": "Natural Extraction",
            }

            mapped_method = None
            if method_val is not None:
                try:
                    mv = str(method_val).strip().lower()
                    mapped_method = method_label_map.get(mv, method_val)
                except Exception:
                    mapped_method = method_val

            # apply fallbacks
            row_dict["microorganism"] = mapped_method if mapped_method else (method_val if method_val is not None else "unknown")
            row_dict["strain"] = process_type_val if process_type_val is not None else "unknown"

            results.append(
                {
                    "rank": rank,
                    "microorganism": row_dict["microorganism"],
                    "strain": row_dict["strain"],
                    "compound_class": row.get("compound_class"),
                    "compound_name": row.get("compound_name"),
                    "expected_yield": row.get("expected_yield"),
                    "confidence_score": row.get("confidence_score"),
                    "temperature_range": row.get("temperature_range"),
                    "ph_range": row.get("ph_range"),
                    "oxygen_level": row.get("oxygen_level"),
                    "source_reference": row.get("source_reference"),
                    "decision_score": float(row.get("decision_score", 0.0)),
                    "reasons": reasons,
                    "input_context": input_ctx,
                    "notes": "; ".join(notes),
                }
            )

        return results

    def analyze_process(self, process_input: Dict[str, Any]) -> List[Dict[str, Any]]:
        """Public wrapper to analyze an existing process.

        This forwards to the existing recommend(...) path which detects
        process-like inputs and returns a single synthetic strategy containing
        an 'analysis' dict. Returns the same list structure as recommend.
        """
        try:
            # Prefer direct analyzeProcess implementation if available
            analysis = self.analyzeProcess(process_input)
            # Wrap into the same list structure as recommend (selected strategy)
            selected = {
                "microorganism": process_input.get("method") or "process",
                "strain": f"steps:{process_input.get('number_of_steps')}",
                "compound_class": None,
                "compound_name": process_input.get("molecule_name") or process_input.get("target_molecule") or "unknown",
                "expected_yield": analysis.get("efficiency_score"),
                "confidence_score": 1.0,
                "decision_score": analysis.get("efficiency_score"),
                "source_reference": "ProcessOptimizationEngine",
                "analysis": analysis,
                "input_context": process_input,
            }
            return [selected]
        except Exception:
            logging.exception("analyze_process wrapper failed")
            return []

    def analyzeProcess(self, inp: Dict[str, Any], properties: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
        """Analyze an existing process and return deterministic analysis results.

        Returns a dict with keys: efficiency, efficiency_score, cost, cost_score,
        risk, toxicity (and final_toxicity for backward compatibility), issues, improvements, properties
        """
        try:
            # Normalize input and apply defaults
            p = {k: v for k, v in (inp or {}).items()}
            molecule = str(p.get("molecule_name") or p.get("target_molecule") or "").lower().strip()
            # default properties
            if properties is None:
                molecule_props = {
                    "vanillin": {"complexity": "medium", "purification_difficulty": "medium", "stability": "high", "base_toxicity": "low"},
                    "linalool": {"complexity": "medium", "purification_difficulty": "medium", "stability": "medium", "base_toxicity": "low"},
                    "citral": {"complexity": "medium", "purification_difficulty": "medium", "stability": "low", "base_toxicity": "medium"},
                }
                properties = molecule_props.get(molecule, {"complexity": "medium", "purification_difficulty": "medium", "stability": "medium", "base_toxicity": "medium"})

            # Derive structural properties from SMILES when available (or from known molecule names).
            # Find SMILES: prefer explicit user-supplied SMILES, else lookup in central molecule DB
            smiles_candidate = None
            try:
                smiles_candidate = p.get("smiles")
            except Exception:
                smiles_candidate = None
            if not smiles_candidate:
                try:
                    smiles_candidate = get_smiles_for(molecule)
                    if not smiles_candidate or smiles_candidate == "Not available":
                        smiles_candidate = None
                except Exception:
                    smiles_candidate = None

            # Extract SMILES-derived properties and merge (SMILES-derived override defaults)
            struct_props = {}
            try:
                struct_props = self._extract_smiles_properties(smiles_candidate) or {}
                # merge without losing previously-provided keys (SMILES-derived wins)
                merged = dict(properties or {})
                merged.update(struct_props or {})
                properties = merged
            except Exception:
                # keep original properties on failure
                struct_props = struct_props or {}

            # --- Molecule type classification layer ---
            # Read molecule type from input (default: small_molecule)
            molecule_type = str(p.get("molecule_type") or "small_molecule").lower()

            def getMoleculeTypeProperties(t: str) -> Dict[str, Any]:
                if not t:
                    t = "small_molecule"
                key = str(t).lower()
                if key == "small_molecule":
                    return {
                        "purification_difficulty": "low",
                        "stability": "high",
                        "cost_driver": "energy",
                        "downstream_complexity": "low",
                    }
                if key == "peptide":
                    return {
                        "purification_difficulty": "medium",
                        "stability": "medium",
                        "cost_driver": "purification",
                        "downstream_complexity": "medium",
                    }
                if key == "protein":
                    return {
                        "purification_difficulty": "high",
                        "stability": "low",
                        "cost_driver": "chromatography",
                        "downstream_complexity": "high",
                    }
                if key == "natural_product":
                    return {
                        "purification_difficulty": "very high",
                        "stability": "medium",
                        "cost_driver": "purification",
                        "downstream_complexity": "very high",
                    }
                return {
                    "purification_difficulty": "medium",
                    "stability": "medium",
                    "cost_driver": "unknown",
                    "downstream_complexity": "medium",
                }

            def worst(a: Optional[str], b: Optional[str]) -> str:
                order = ["low", "medium", "high", "very high"]
                try:
                    ia = order.index(a) if a in order else 1
                except Exception:
                    ia = 1
                try:
                    ib = order.index(b) if b in order else 1
                except Exception:
                    ib = 1
                return order[max(ia, ib)]

            type_props = getMoleculeTypeProperties(molecule_type)
            # Build merged properties where molecule type and SMILES-derived properties are combined
            smiles_props = struct_props or {}
            merged_props = dict(properties or {})
            merged_props["purification_difficulty"] = worst(
                type_props.get("purification_difficulty"), smiles_props.get("purification_difficulty", merged_props.get("purification_difficulty", "medium"))
            )
            merged_props["stability"] = worst(
                type_props.get("stability"), smiles_props.get("stability", merged_props.get("stability", "medium"))
            )
            # toxicity and complexity primarily from SMILES-derived heuristics when available
            if smiles_props.get("toxicity") is not None:
                merged_props["toxicity"] = smiles_props.get("toxicity")
            if smiles_props.get("complexity") is not None:
                merged_props["complexity"] = smiles_props.get("complexity")
            # cost_driver and downstream complexity come from molecule type
            merged_props["cost_driver"] = type_props.get("cost_driver")
            merged_props["downstream_complexity"] = type_props.get("downstream_complexity")
            # keep base_toxicity if present
            if smiles_props.get("base_toxicity") is not None:
                merged_props["base_toxicity"] = smiles_props.get("base_toxicity")

            # expose molecule_type and type_properties for transparency
            merged_props["molecule_type"] = molecule_type
            merged_props["type_properties"] = type_props

            properties = merged_props

            # --- SMILES weight helper: scale the influence of structural SMILES heuristics
            def getSmilesWeight(mtype: str) -> float:
                try:
                    key = str(mtype).lower()
                except Exception:
                    return 0.5
                if key == "small_molecule":
                    return 1.0
                if key == "peptide":
                    return 0.6
                if key == "protein":
                    return 0.3
                if key == "natural_product":
                    return 0.8
                return 0.5

            # compute weight for later SMILES-driven adjustments
            smiles_weight = getSmilesWeight(molecule_type)

            # --- Molecule subtype layer (optional) ---
            molecule_subtype = (p.get("molecule_subtype") or None)

            def getSubtypeProperties(t: str, s: Optional[str]) -> Dict[str, Any]:
                props: Dict[str, Any] = {}
                if not t:
                    return props
                key = str(t).lower()
                sub = str(s).lower() if s else None
                # small_molecule subtypes
                if key == "small_molecule":
                    if sub == "volatile":
                        props["purification_difficulty"] = "low"
                        props["cost_driver"] = "distillation"
                        props["notes"] = "Volatile compounds can be separated via distillation"
                    elif sub == "non_volatile":
                        props["purification_difficulty"] = "medium"
                        props["cost_driver"] = "crystallization"
                # peptide subtypes
                if key == "peptide":
                    if sub == "linear":
                        props["purification_difficulty"] = "medium"
                        props["stability"] = "medium"
                    elif sub == "cyclic":
                        props["purification_difficulty"] = "high"
                        props["stability"] = "high"
                        props["notes"] = "Cyclic peptides are more stable but harder to purify"
                # protein subtypes
                if key == "protein":
                    if sub == "antibody":
                        props["purification_difficulty"] = "very high"
                        props["cost_driver"] = "chromatography"
                        props["notes"] = "Antibody purification is highly complex and expensive"
                    elif sub == "enzyme":
                        props["purification_difficulty"] = "high"
                        props["cost_driver"] = "chromatography"
                # natural product subtypes
                if key == "natural_product":
                    if sub == "terpene":
                        props["purification_difficulty"] = "medium"
                        props["cost_driver"] = "distillation"
                    elif sub == "alkaloid":
                        props["purification_difficulty"] = "high"
                        props["cost_driver"] = "extraction"
                return props

            subtype_props = getSubtypeProperties(molecule_type, molecule_subtype)
            # Merge subtype properties: worst for purification/stability; subtype cost_driver overrides
            try:
                if subtype_props:
                    properties["purification_difficulty"] = worst(properties.get("purification_difficulty"), subtype_props.get("purification_difficulty", properties.get("purification_difficulty")))
                    properties["stability"] = worst(properties.get("stability"), subtype_props.get("stability", properties.get("stability")))
                    if subtype_props.get("cost_driver"):
                        properties["cost_driver"] = subtype_props.get("cost_driver")
                    if subtype_props.get("notes"):
                        subtype_notes.append(subtype_props.get("notes"))
                    # Subtype-driven cost/issue adjustments
                    cd = subtype_props.get("cost_driver")
                    if cd == "chromatography":
                        cost = "very high"
                        if "Chromatography-intensive purification increases cost" not in issues:
                            issues.append("Chromatography-intensive purification increases cost")
                    if cd == "distillation":
                        if "Distillation-based separation likely feasible" not in issues:
                            issues.append("Distillation-based separation likely feasible")
                    if subtype_props.get("purification_difficulty") == "very high":
                        if "Subtype indicates extremely difficult purification" not in issues:
                            issues.append("Subtype indicates extremely difficult purification")
                        cost = "very high"
                    # Subtype-driven improvements
                    stype = molecule_subtype
                    if stype == "volatile":
                        if "Leverage fractional distillation for efficient recovery of volatiles" not in improvements:
                            improvements.append("Leverage fractional distillation for efficient recovery of volatiles; evaluate reflux ratios and entrainers to avoid azeotropes")
                    if stype == "cyclic":
                        if "Develop advanced chromatographic methods for cyclic peptides (optimize resin and gradient)" not in improvements:
                            improvements.append("Develop advanced chromatographic methods for cyclic peptides (optimize resin, gradient and desolvation) to handle conformational heterogeneity")
                    if stype == "antibody":
                        if "Optimize chromatography loading, residence time and column regeneration to reduce resin and buffer consumption" not in improvements:
                            improvements.append("Optimize chromatography loading, residence time and column regeneration to reduce resin and buffer consumption")
                        # subtype antibody is expensive and risky - ensure scoring reflects that later
                        # (we add conservative adjustments in the subtype scoring block below)
            except Exception:
                pass

            # Initialize issue/improvement lists early so structural checks can append
            issues: List[str] = []
            improvements: List[str] = []
            subtype_notes: List[str] = []

            # initialize variables
            efficiency = "high"
            cost = "medium"
            risk = "low"
            toxicity = properties.get("base_toxicity", "medium")

            # helpers to bump level
            levels = ["low", "medium", "high", "very high"]

            def bump(level: str, steps_up: int = 1) -> str:
                try:
                    idx = levels.index(level)
                except ValueError:
                    idx = 1
                return levels[min(idx + steps_up, len(levels) - 1)]

            # read inputs safely
            try:
                steps = int(p.get("number_of_steps") or 0)
            except Exception:
                steps = 0
            method = str(p.get("method") or "").lower()
            desired_purity = str(p.get("desired_purity") or "standard").lower()
            scale = str(p.get("scale") or "lab").lower()
            raw_avail = str(p.get("raw_material_availability") or "medium").lower()
            raw_cost = str(p.get("raw_material_cost") or "medium").lower()
            strict_waste = bool(p.get("strict_waste_constraints"))
            has_bioreactor = bool(p.get("has_bioreactor"))
            has_advanced_pur = bool(p.get("has_advanced_purification"))

            # --- New internal scoring system (decision-driven) ---
            # costScore and riskScore accumulate drivers; efficiencyScore starts high
            costScore = 0
            riskScore = 0
            efficiencyScore = 10

            # --- Method impact: ensure method influences cost/risk/efficiency ---
            try:
                if method in ("chemical", "chem", "chemical synthesis"):
                    # chemical syntheses tend to add energy-driven costs
                    if "Chemical synthesis may involve energy-intensive steps" not in issues:
                        issues.append("Chemical synthesis may involve energy-intensive steps")
                    # reflect moderate cost pressure in costScore
                    try:
                        costScore += 2
                    except Exception:
                        pass
                if method in ("biotechnological", "biotech", "biotechnological synthesis"):
                    if "Biotech processes require careful control of biological systems" not in issues:
                        issues.append("Biotech processes require careful control of biological systems")
                    # biotech often increases operational risk
                    risk = bump(risk, 1)
                if method in ("extraction", "extract", "extraction-based"):
                    if "Extraction processes may have low yield and high variability" not in issues:
                        issues.append("Extraction processes may have low yield and high variability")
                    efficiency = "low"
            except Exception:
                pass

            # Step 4: Efficiency logic
            if steps > 5:
                efficiency = "low"
                issues.append("High number of synthesis steps reduces efficiency")
            elif steps > 3:
                efficiency = "medium"
            else:
                efficiency = "high"

            # STEP 3: Number of steps direct impact (ensure clear effect)
            if steps > 5:
                # enforce cost and efficiency penalties for very long routes
                cost = "high"
                efficiency = "low"
                if "High number of synthesis steps increases cost and complexity" not in issues:
                    issues.append("High number of synthesis steps increases cost and complexity")

            # SMILES-derived issue generation and suggested improvements (early)
            try:
                # Use structural signals to produce concrete, chemistry-aware statements
                arom = int(properties.get("aromatic_rings") or 0)
                long_chain = bool(properties.get("long_aliphatic_chain"))
                polar_ct = int(properties.get("polar_group_count") or 0)
                ester = bool(properties.get("ester_present"))
                amide = bool(properties.get("amide_present"))

                if properties.get("purification_difficulty") == "high":
                    reason = []
                    if long_chain:
                        reason.append("long hydrophobic alkyl chain -> increased solvent use and mixed-mode separations")
                    if arom:
                        reason.append("aromatic system -> strong π–π interactions can complicate chromatographic separations")
                    if polar_ct >= 3:
                        reason.append("multiple polar groups -> solubility-driven workup constraints")
                    msg = "Difficult purification expected due to: " + ", ".join(reason) if reason else "Difficult downstream separations inferred from structure"
                    if msg not in issues:
                        issues.append(msg)
                    # concrete improvement
                    if "Evaluate crystallization or tailored chromatography protocols" not in improvements:
                        improvements.append("Evaluate crystallization or tailored chromatography protocols (optimize solvent, seed, and temperature) to reduce HPLC time and solvent consumption")

                # Toxicity leads to explicit waste handling/neutralization needs
                if (properties.get("toxicity") == "high") or (properties.get("base_toxicity") == "high"):
                    tox_msg = "Potential hazardous fragments (halogens or reactive carbonyls) increase waste treatment and disposal requirements"
                    if tox_msg not in issues:
                        issues.append(tox_msg)
                    if "Design safer reagent choices and include effluent neutralization in process flow" not in improvements:
                        improvements.append("Design safer reagent choices and include effluent neutralization and segregation steps to meet waste regulations")

                # Structural complexity guidance
                if properties.get("complexity") == "high":
                    comp_msg = "High structural complexity (rings/branches) increases synthetic steps and purification burden"
                    if comp_msg not in issues:
                        issues.append(comp_msg)
                    if "Assess convergent synthesis or protecting-group minimization to reduce step count" not in improvements:
                        improvements.append("Assess convergent synthesis or protecting-group minimization to reduce step count and improve overall yield")

                # Aromatic-specific note
                if arom:
                    arom_msg = f"Aromatic rings detected (count={arom}) — increases chemical stability and often enables crystallization-based purification"
                    if arom_msg not in issues:
                        issues.append(arom_msg)
                    if "Explore crystallization protocols leveraging aromaticity for solid-form separation" not in improvements:
                        improvements.append("Explore crystallization protocols leveraging aromaticity for solid-form separation to lower downstream chromatography load")

                # Long-chain (hydrophobic) consequences
                if long_chain:
                    lc_msg = "Long aliphatic chain detected — increased hydrophobicity will raise solvent consumption and complicate aqueous workups"
                    if lc_msg not in issues:
                        issues.append(lc_msg)
                    if "Consider phase-separation or flash chromatography solvent systems targeting hydrophobic species" not in improvements:
                        improvements.append("Consider phase-separation or flash chromatography solvent systems targeting hydrophobic species; evaluate greener solvents to reduce cost")

                # Polar groups and hydrolysis-sensitive motifs
                if ester or amide:
                    hydro_msg = "Ester/amide motifs detected — hydrolysis risk under acidic/basic conditions can affect stability during workup"
                    if hydro_msg not in issues:
                        issues.append(hydro_msg)
                    if "Avoid harsh pH during workup; consider protecting labile groups or ester-to-amide exchange strategies" not in improvements:
                        improvements.append("Avoid harsh pH during workup; consider protecting labile groups or enzymatic/chemoselective transformations to improve stability")
            except Exception:
                pass

            # Step 5: Decision-driven scoring system
            # Initialize scores (decision-driven core)
            costScore = 0
            riskScore = 0
            efficiencyScore = 10  # start high, reduce with problems

            # BASE INPUT EFFECTS
            # Number of steps
            try:
                if steps > 5:
                    costScore += 3
                    riskScore += 2
                    efficiencyScore -= 3
                    if "High number of synthesis steps increases cost and complexity" not in issues:
                        issues.append("High number of synthesis steps increases cost and complexity")
                elif steps > 3:
                    costScore += 1
                    efficiencyScore -= 1
            except Exception:
                pass

            # Purity
            try:
                if desired_purity in (">99%", "very high"):
                    costScore += 3
                    riskScore += 2
                    if "High purity significantly increases downstream complexity" not in issues:
                        issues.append("High purity significantly increases downstream complexity")
                elif desired_purity == "high":
                    costScore += 2
            except Exception:
                pass

            # Scale
            try:
                if scale == "industrial":
                    costScore += 3
                    riskScore += 3
                elif scale == "pilot":
                    riskScore += 1
                elif scale == "lab":
                    # small risk reduction for lab
                    riskScore = max(0, riskScore - 1)
            except Exception:
                pass

            # Raw materials
            try:
                if raw_cost == "high":
                    costScore += 3
                if raw_avail == "low":
                    riskScore += 3
                    if "Limited raw material availability may affect scalability" not in issues:
                        issues.append("Limited raw material availability may affect scalability")
            except Exception:
                pass

            # Waste constraints
            try:
                if strict_waste:
                    costScore += 2
                    riskScore += 2
                    if "Strict waste regulations increase process cost and complexity" not in issues:
                        issues.append("Strict waste regulations increase process cost and complexity")
            except Exception:
                pass

            # METHOD EFFECTS (ensure method influences scores)
            try:
                if method in ("chemical", "chem", "chemical synthesis"):
                    costScore += 1
                if method in ("biotechnological", "biotech", "biotechnological synthesis"):
                    riskScore += 2
                if method in ("extraction", "extract", "extraction-based"):
                    efficiencyScore -= 2
                    riskScore += 1
            except Exception:
                pass

            # SMILES / STRUCTURE EFFECTS
            try:
                # Scale SMILES-derived effects by the molecule-type-aware weight
                w = smiles_weight if 'smiles_weight' in locals() else 0.5
                if properties.get("purification_difficulty") in ("high", "very high"):
                    try:
                        costScore += int(round(3 * w))
                    except Exception:
                        costScore += 3
                if properties.get("stability") == "low":
                    try:
                        riskScore += int(round(3 * w))
                    except Exception:
                        riskScore += 3
                if properties.get("toxicity") == "high" or properties.get("base_toxicity") == "high":
                    try:
                        costScore += int(round(2 * w))
                        riskScore += int(round(2 * w))
                    except Exception:
                        costScore += 2
                        riskScore += 2
                if properties.get("complexity") == "high":
                    try:
                        costScore += int(round(2 * w))
                    except Exception:
                        costScore += 2
                # Structural complexity signal: very long SMILES strings indicate large/complex molecules
                try:
                    if smiles_candidate and isinstance(smiles_candidate, str) and len(smiles_candidate) > 100:
                        riskScore += 2
                        costScore += 2
                        if "Large and complex molecular structure increases process difficulty" not in issues:
                            issues.append("Large and complex molecular structure increases process difficulty")
                except Exception:
                    pass
            except Exception:
                pass

            # MOLECULE TYPE EFFECTS
            try:
                if molecule_type == "protein":
                    costScore += 3
                    riskScore += 3
                elif molecule_type == "natural_product":
                    costScore += 3
                    riskScore += 2
                elif molecule_type == "peptide":
                    costScore += 1
                    riskScore += 1
            except Exception:
                pass

            # SUBTYPE EFFECTS
            try:
                stype = (p.get("molecule_subtype") or "")
                if stype == "antibody":
                    costScore += 4
                    riskScore += 3
                    # ensure additional antibody-specific signal (keeps behavior but guarantees impact)
                    try:
                        costScore += 3
                        riskScore += 3
                    except Exception:
                        pass
                if stype == "cyclic":
                    costScore += 2
                if stype == "non_volatile":
                    costScore += 2
                if stype == "volatile":
                    # volatile small-molecule subtypes often reduce purification cost (distillation)
                    try:
                        costScore = max(0, costScore - 2)
                    except Exception:
                        pass
            except Exception:
                pass

            # BIOMOLECULE-SPECIFIC (only peptides & proteins)
            try:
                isBiomolecule = molecule_type in ("peptide", "protein")
                # for biomolecules, structural SMILES heuristics are less predictive; warn when weight is low
                try:
                    if isBiomolecule and (smiles_weight if 'smiles_weight' in locals() else 0.5) < 0.5:
                        if "For large biomolecules, structural representation has limited predictive power compared to biophysical behavior" not in issues:
                            issues.append("For large biomolecules, structural representation has limited predictive power compared to biophysical behavior")
                except Exception:
                    pass
                if isBiomolecule:
                    agg = str(p.get("aggregation_risk") or "").lower()
                    fold = str(p.get("folding_complexity") or "").lower()
                    bst = str(p.get("biophysical_stability") or "").lower()
                    if agg == "high":
                        riskScore += 4
                        costScore += 2
                        efficiencyScore -= 3
                        if "High aggregation risk may significantly impact downstream processing" not in issues:
                            issues.append("High aggregation risk may significantly impact downstream processing")
                        if "Optimize conditions to reduce aggregation (e.g. buffer or temperature control)" not in improvements:
                            improvements.append("Optimize conditions to reduce aggregation (e.g. buffer or temperature control)")
                    elif agg == "medium":
                        issues.append("Moderate aggregation risk may affect purification efficiency")
                    if fold == "high":
                        riskScore += 3
                        efficiencyScore -= 2
                        if "Complex folding increases process sensitivity and instability" not in issues:
                            issues.append("Complex folding increases process sensitivity and instability")
                        if "Stabilize folding conditions to improve process robustness" not in improvements:
                            improvements.append("Stabilize folding conditions to improve process robustness")
                    elif fold == "medium":
                        issues.append("Folding complexity may impact reproducibility and yield")
                    if bst == "low":
                        riskScore += 4
                        costScore += 2
                        efficiencyScore -= 4
                        if "Low molecular stability may lead to degradation during processing" not in issues:
                            issues.append("Low molecular stability may lead to degradation during processing")
                        if "Adjust process conditions to improve molecular stability" not in improvements:
                            improvements.append("Adjust process conditions to improve molecular stability")
                    elif bst == "medium":
                        issues.append("Moderate stability may require careful process control")
            except Exception:
                pass

            # INTERACTION EFFECTS
            try:
                if scale == "industrial" and steps > 5:
                    costScore += 3
                    riskScore += 3
                if desired_purity in (">99%", "very high") and properties.get("purification_difficulty") in ("high", "very high"):
                    costScore += 4
                    riskScore += 3
                if isBiomolecule and str(p.get("aggregation_risk") or "").lower() == "high" and scale == "industrial":
                    riskScore += 4
            except Exception:
                pass

            # --- Score breakdown & explanation augmentation ---
                # Track contributions per category so outputs are explainable
                try:
                    # weights (recompute if needed)
                    try:
                        w_smiles = getSmilesWeight(molecule_type)
                    except Exception:
                        w_smiles = 0.5
                    w_process = 1.0
                    w_economic = 1.2
                    w_downstream = 1.1
                    w_biophysical = 1.3 if molecule_type in ("protein", "peptide") else 0.0

                    # SMILES-derived numeric proxies
                    complexity_map = {"low": 1, "medium": 2, "high": 3, "very high": 4}
                    stability_map = {"high": 1, "medium": 2, "low": 3}
                    comp_level = properties.get("complexity", "medium")
                    stab_level = properties.get("stability", "medium")
                    comp_num = complexity_map.get(comp_level, 2)
                    instability_num = stability_map.get(stab_level, 2)

                    scoreBreakdown = {"smiles": 0.0, "process": 0.0, "economic": 0.0, "downstream": 0.0, "biophysical": 0.0}
                    riskBreakdown = {"smiles": 0.0, "process": 0.0, "economic": 0.0, "downstream": 0.0, "biophysical": 0.0}

                    # smiles contributions
                    delta_smiles_cost = comp_num * w_smiles
                    delta_smiles_risk = instability_num * w_smiles
                    try:
                        print(f"[DEBUG] comp_num={comp_num}, instability_num={instability_num}, w_smiles={w_smiles}")
                        print(f"[DEBUG] delta_smiles_cost={delta_smiles_cost}, delta_smiles_risk={delta_smiles_risk}")
                    except Exception:
                        pass
                    if properties.get("base_toxicity") == "high" or properties.get("toxicity") == "high":
                        delta_smiles_cost += 2 * w_smiles
                        delta_smiles_risk += 2 * w_smiles
                    if properties.get("complexity") == "high":
                        delta_smiles_cost += 2 * w_smiles
                    if smiles_candidate and isinstance(smiles_candidate, str) and len(smiles_candidate) > 100:
                        delta_smiles_cost += 2 * w_smiles
                        delta_smiles_risk += 2 * w_smiles
                    scoreBreakdown["smiles"] += float(delta_smiles_cost)
                    riskBreakdown["smiles"] += float(delta_smiles_risk)

                    # process contributions (steps, scale, method)
                    delta_process_cost = steps * 1.2 * w_process
                    delta_process_risk = 0.0
                    if scale == "industrial":
                        delta_process_cost += 3 * w_process
                        delta_process_risk += 3 * w_process
                    # method cost/risk nudges
                    if method in ("chemical", "chem", "chemical synthesis"):
                        delta_process_cost += 2 * w_process
                    if method in ("biotechnological", "biotech", "biotechnological synthesis"):
                        delta_process_risk += 2 * w_process
                    if method in ("extraction", "extract", "extraction-based"):
                        delta_process_risk += 1 * w_process
                    scoreBreakdown["process"] += float(delta_process_cost)
                    riskBreakdown["process"] += float(delta_process_risk)

                    # economic contributions
                    delta_econ_cost = 0.0
                    delta_econ_risk = 0.0
                    if raw_cost == "high":
                        delta_econ_cost += 4 * w_economic
                    if raw_avail == "low":
                        delta_econ_risk += 3 * w_economic
                    if strict_waste:
                        delta_econ_cost += 2 * w_economic
                        delta_econ_risk += 1 * w_economic
                    scoreBreakdown["economic"] += float(delta_econ_cost)
                    riskBreakdown["economic"] += float(delta_econ_risk)

                    # downstream contributions
                    delta_down_cost = 0.0
                    delta_down_risk = 0.0
                    if properties.get("purification_difficulty") in ("high", "very high"):
                        delta_down_cost += 3 * w_downstream
                    if desired_purity in (">99%", "very high"):
                        delta_down_cost += 3 * w_downstream
                        delta_down_risk += 2 * w_downstream
                    # subtype-driven downstream costs (antibody etc.)
                    stype = (p.get("molecule_subtype") or "").lower()
                    if stype == "antibody":
                        delta_down_cost += 4
                        delta_down_risk += 3
                    scoreBreakdown["downstream"] += float(delta_down_cost)
                    riskBreakdown["downstream"] += float(delta_down_risk)

                    # biophysical contributions
                    delta_bio_cost = 0.0
                    delta_bio_risk = 0.0
                    if w_biophysical > 0:
                        agg = str(p.get("aggregation_risk") or "").lower()
                        fold = str(p.get("folding_complexity") or "").lower()
                        bst = str(p.get("biophysical_stability") or "").lower()
                        if agg == "high":
                            delta_bio_risk += 4 * w_biophysical
                            delta_bio_cost += 2 * w_biophysical
                        if fold == "high":
                            delta_bio_risk += 3 * w_biophysical
                        if bst == "low":
                            delta_bio_risk += 4 * w_biophysical
                            delta_bio_cost += 2 * w_biophysical
                    scoreBreakdown["biophysical"] += float(delta_bio_cost)
                    riskBreakdown["biophysical"] += float(delta_bio_risk)

                    # expose numeric breakdowns
                    score_breakdown = scoreBreakdown
                    risk_breakdown = riskBreakdown
                except Exception as e:
                    # fallback - record printed exception for debugging in terminal
                    try:
                        print("[DEBUG] score breakdown calculation failed:", repr(e))
                    except Exception:
                        pass
                    score_breakdown = {"smiles": 0.0, "process": 0.0, "economic": 0.0, "downstream": 0.0, "biophysical": 0.0}
                    risk_breakdown = {"smiles": 0.0, "process": 0.0, "economic": 0.0, "downstream": 0.0, "biophysical": 0.0}

                # ensure efficiencyScore within 0..10 (STEP 5)
                try:
                    efficiencyScore = max(0, min(10, int(round(efficiencyScore))))
                except Exception:
                    efficiencyScore = max(0, min(10, efficiencyScore if isinstance(efficiencyScore, (int, float)) else 10))

                # Baseline risk for biomolecules (STEP 7)
                try:
                    if molecule_type == "peptide":
                        riskScore += 2
                        risk_breakdown["biophysical"] = risk_breakdown.get("biophysical", 0) + 2
                    if molecule_type == "protein":
                        riskScore += 3
                        risk_breakdown["biophysical"] = risk_breakdown.get("biophysical", 0) + 3
                except Exception:
                    pass

                # Human-readable cost drivers (STEP 3)
                try:
                    costDrivers_new: List[str] = []
                    if score_breakdown.get("process", 0) > 2:
                        costDrivers_new.append("Multiple synthesis steps increase operational complexity")
                    if score_breakdown.get("economic", 0) > 2:
                        costDrivers_new.append("High raw material cost contributes significantly to total cost")
                    if score_breakdown.get("downstream", 0) > 2:
                        costDrivers_new.append("Challenging purification increases downstream cost")
                    if score_breakdown.get("smiles", 0) > 2:
                        costDrivers_new.append("Molecular structure increases process complexity")
                except Exception:
                    costDrivers_new = []

                # Structure-specific insight (STEP 8)
                try:
                    insights = []
                    try:
                        arom_ct = int(properties.get("aromatic_rings", 0) or 0)
                    except Exception:
                        arom_ct = 0
                    if arom_ct > 0:
                        insights.append(f"Aromatic rings detected (count={arom_ct}) — supports crystallization screening and potential solid-form separation")
                    if properties.get("long_aliphatic_chain"):
                        insights.append("Long hydrophobic chain increases solvent demand and suggests flash chromatography or phase-separation strategies")
                    if properties.get("ester_present"):
                        insights.append("Ester functionality present — hydrolysis risk during aqueous workup; control pH and consider protecting strategies")
                    if properties.get("amide_present"):
                        insights.append("Amide groups present — monitor stereochemical integrity and hydrolytic stability during prolonged processing")
                    try:
                        polar_ct = int(properties.get("polar_group_count", 0) or 0)
                    except Exception:
                        polar_ct = 0
                    if polar_ct >= 3:
                        insights.append("Multiple polar groups suggest reversed-phase HPLC or polarity-based partitioning for purification")
                    structure_insight = "; ".join(insights) if insights else None
                except Exception:
                    structure_insight = None

                # Downstream detail (STEP 9)
                try:
                    downstream_detail = None
                    if stype == "volatile":
                        downstream_detail = "Distillation is likely a viable and cost-efficient separation method"
                    elif stype == "non_volatile":
                        downstream_detail = "Crystallization or solvent-based purification is expected"
                    elif molecule_type == "protein":
                        downstream_detail = "Chromatography-based purification is required due to molecular complexity"
                    else:
                        # fallback to earlier downstream insights if available
                        downstream_detail = "; ".join(downstreamInsights) if downstreamInsights else None
                except Exception:
                    downstream_detail = None

                # Enhance weighting explanation (STEP 10)
                try:
                    weighting_explanation = (
                        "Final assessment is based on weighted contributions from molecular structure (SMILES), process complexity, economic constraints, "
                        "and downstream/biophysical requirements. See score_breakdown for numeric contributions."
                    )
                except Exception:
                    weighting_explanation = "weights unavailable"

                # Convert scores into qualitative levels
            def mapScore(score: int) -> str:
                try:
                    s = float(score)
                except Exception:
                    s = 0.0
                if s < 5:
                    return "low"
                if s < 10:
                    return "medium"
                if s < 15:
                    return "high"
                return "very high"

            # Final qualitative assignments driven by scores
            # Before mapping to qualitative labels, run process-class adjustments (small +/-1..3)
            try:
                process_class = self.classifyProcess(p, properties)
            except Exception:
                process_class = "unknown"

            # Apply only small, incremental adjustments to numeric scores based on process class
            try:
                if process_class == "commodity_small_molecule":
                    costScore = max(0, costScore - 2)
                    riskScore = max(0, riskScore - 2)
                elif process_class == "volatile_small_molecule":
                    costScore = max(0, costScore - 1)
                elif process_class == "complex_small_molecule":
                    costScore += 1
                    riskScore += 1
                elif process_class == "difficult_purification_small_molecule":
                    costScore += 2
                elif process_class == "linear_peptide":
                    riskScore += 1
                    costScore += 1
                elif process_class == "cyclic_peptide":
                    riskScore += 2
                    costScore += 1
                elif process_class == "enzyme":
                    riskScore += 1
                    costScore += 1
                elif process_class == "antibody":
                    riskScore += 2
                    costScore += 2
                elif process_class == "high_complexity_biologic":
                    riskScore += 2
                    costScore += 2
                elif process_class == "natural_product_complex":
                    costScore += 1
                    riskScore += 1
            except Exception:
                pass

            # record process class for transparency
            try:
                process_class_explanation = "Additional adjustments based on recognized process class"
            except Exception:
                process_class_explanation = "Process class adjustment applied"

            cost = mapScore(costScore)
            risk = mapScore(riskScore)
            if efficiencyScore >= 8:
                efficiency = "high"
            elif efficiencyScore >= 5:
                efficiency = "medium"
            else:
                efficiency = "low"

            # numeric proxy for efficiency (0..1)
            try:
                efficiency_score = max(0.0, min(1.0, float(efficiencyScore) / 10.0))
            except Exception:
                efficiency_score = 0.5

            # (cost already mapped via mapScore)

            # Adapt cost logic based on molecule type-driven cost_driver
            try:
                cd = str(properties.get("cost_driver") or "").lower()
                if cd == "chromatography":
                    issues.append("Chromatography likely a major cost driver")
                elif cd == "energy":
                    issues.append("Energy-intensive separation steps likely dominate cost")
            except Exception:
                pass

            # Extremely challenging purification raises cost sharply
            if properties.get("purification_difficulty") == "very high":
                cost = "very high"
                if "Extremely challenging purification expected" not in issues:
                    issues.append("Extremely challenging purification expected")

            # Step 6: Toxicity logic
            tox_order = {"low": 0, "medium": 1, "high": 2, "very high": 3}
            tox_level = tox_order.get(properties.get("base_toxicity", "medium"), 1)
            if method == "chemical" and steps > 5:
                tox_level = min(tox_level + 1, 3)
            if scale == "industrial":
                tox_level = min(tox_level + 1, 3)
            if strict_waste:
                issues.append("Process may conflict with waste regulations")
            tox_rev = {0: "low", 1: "medium", 2: "high", 3: "very high"}
            toxicity = tox_rev.get(tox_level, "medium")

            # Step 7: Risk logic
            if steps > 5:
                risk = "medium"
            if properties.get("stability") == "low":
                risk = "high"
                issues.append("Low molecular stability increases process risk")
            if scale == "industrial" and efficiency == "low":
                risk = "high"

            # STEP 5: Scale impact on risk/cost
            try:
                if scale == "lab":
                    # lab scale tends to be lower operational risk: reduce one level safely
                    try:
                        idx = levels.index(risk) if risk in levels else 1
                    except Exception:
                        idx = 1
                    risk = levels[max(0, idx - 1)]
                elif scale == "pilot":
                    risk = bump(risk, 1)
                elif scale == "industrial":
                    # industrial increases both cost and risk
                    if "Scaling introduces operational and mixing challenges" not in issues:
                        issues.append("Scaling introduces operational and mixing challenges")
                    cost = "high"
                    risk = bump(risk, 2)
            except Exception:
                pass

            # Molecule-type specific risk adjustments (record issues; scores are authoritative)
            try:
                if molecule_type == "protein":
                    if "Protein stability and aggregation risk must be considered" not in issues:
                        issues.append("Protein stability and aggregation risk must be considered")
                if molecule_type == "natural_product":
                    if "Complex mixture of similar compounds increases process risk" not in issues:
                        issues.append("Complex mixture of similar compounds increases process risk")
            except Exception:
                pass

            # --- Biophysical inputs (ONLY for peptides & proteins) ---
            try:
                is_biomolecule = molecule_type in ("peptide", "protein")
                agg_risk = None
                fold_complex = None
                bio_stability = None
                if is_biomolecule:
                    agg_risk = str(p.get("aggregation_risk") or "").lower() or None
                    fold_complex = str(p.get("folding_complexity") or "").lower() or None
                    bio_stability = str(p.get("biophysical_stability") or "").lower() or None
                    # expose in properties for transparency
                    properties["aggregation_risk"] = agg_risk
                    properties["folding_complexity"] = fold_complex
                    properties["biophysical_stability"] = bio_stability

                    # Aggregation effects (record issues/improvements; scoring handled earlier)
                    if agg_risk == "high":
                        if "High aggregation risk may significantly impact downstream processing" not in issues:
                            issues.append("High aggregation risk may significantly impact downstream processing")
                        if "Optimize conditions to reduce aggregation (e.g. buffer or temperature control)" not in improvements:
                            improvements.append("Optimize conditions to reduce aggregation (e.g. buffer or temperature control)")
                    elif agg_risk == "medium":
                        if "Moderate aggregation risk may affect purification efficiency" not in issues:
                            issues.append("Moderate aggregation risk may affect purification efficiency")

                    # Folding complexity (record only)
                    if fold_complex == "high":
                        if "Complex folding increases process sensitivity and instability" not in issues:
                            issues.append("Complex folding increases process sensitivity and instability")
                        if "Stabilize folding conditions to improve process robustness" not in improvements:
                            improvements.append("Stabilize folding conditions to improve process robustness")
                    elif fold_complex == "medium":
                        if "Folding complexity may impact reproducibility and yield" not in issues:
                            issues.append("Folding complexity may impact reproducibility and yield")

                    # Biophysical stability (record only)
                    if bio_stability == "low":
                        if "Low molecular stability may lead to degradation during processing" not in issues:
                            issues.append("Low molecular stability may lead to degradation during processing")
                        if "Adjust process conditions to improve molecular stability" not in improvements:
                            improvements.append("Adjust process conditions to improve molecular stability")
                    elif bio_stability == "medium":
                        if "Moderate stability may require careful process control" not in issues:
                            issues.append("Moderate stability may require careful process control")
            except Exception:
                pass

            # STEP 12: Interaction effects
            try:
                if scale == "industrial" and steps > 5:
                    if "Industrial scale combined with long synthesis route drives costs very high" not in issues:
                        issues.append("Industrial scale combined with long synthesis route drives costs very high")
                if desired_purity in ("very high", ">99%") and properties.get("purification_difficulty") in ("high", "very high"):
                    if "Very high purity targets with difficult purification significantly increase cost and risk" not in issues:
                        issues.append("Very high purity targets with difficult purification significantly increase cost and risk")
                if raw_cost == "high" and steps > 5:
                    if "Expensive raw materials combined with many steps drive unit economics unfavorably" not in issues:
                        issues.append("Expensive raw materials combined with many steps drive unit economics unfavorably")
            except Exception:
                pass

            # Step 8: Infrastructure checks
            if method in ("biotechnological", "biotech") and not has_bioreactor:
                issues.append("Biotech process without bioreactor limits feasibility")
                risk = "high"

            if desired_purity in (">99%", "very high") and not has_advanced_pur:
                issues.append("Lack of advanced purification may limit achievable purity")
                risk = bump(risk, 1) if risk in levels else risk

            # Step 9: Raw material logic
            costDrivers = []
            if raw_avail == "low":
                msg = "Limited raw material availability — may cause frequent batch failures or sourcing delays affecting scale-up"
                if msg not in issues:
                    issues.append(msg)
                cost = bump(cost, 1)
                costDrivers.append("Limited raw material availability (supply risk)")
            if raw_cost == "high":
                msg = "High-cost precursors directly increase per-unit material cost and can dominate COGS"
                if msg not in issues:
                    issues.append(msg)
                costDrivers.append("Expensive raw materials (high COGS)")

            # Add other cost drivers identified earlier with precise wording
            if steps > 5:
                costDrivers.append(f"Multiple synthesis steps ({steps}) increase operational costs (reagents, labour, unit ops)")
            if desired_purity in ("high", ">99%", "very high"):
                costDrivers.append("Very high purity targets increase downstream separation cost (longer HPLC/extra polishing steps)")
            if properties.get("purification_difficulty") in ("high", "very high"):
                # reference structural reasons when possible
                pd_reasons = []
                if properties.get("long_aliphatic_chain"):
                    pd_reasons.append("hydrophobic chain")
                if int(properties.get("aromatic_rings") or 0) > 0:
                    pd_reasons.append("aromatic system")
                if properties.get("polar_group_count", 0) > 0:
                    pd_reasons.append("multiple polar groups")
                reason_text = " (" + ", ".join(pd_reasons) + ")" if pd_reasons else ""
                costDrivers.append(f"Difficult purification increases chromatography/HPLC time and solvent consumption{reason_text}")
            if strict_waste:
                costDrivers.append("Strict waste handling and disposal requirements (compliance cost)")

            # deduplicate costDrivers preserving order
            seen_cd = set()
            dedup_costDrivers = []
            for d in costDrivers:
                if d not in seen_cd:
                    dedup_costDrivers.append(d)
                    seen_cd.add(d)
            costDrivers = dedup_costDrivers

            # If the computed cost level is at least 'medium' (costScore >=5) ensure 1-3 concrete cost drivers exist
            try:
                if (not costDrivers) and isinstance(costScore, (int, float)) and costScore >= 5:
                    inferred = []
                    if steps and steps > 2:
                        inferred.append(f"Number of synthesis steps ({steps}) increases reagent and unit-op costs")
                    if properties.get("purification_difficulty") in ("high", "very high"):
                        inferred.append("Difficult purification (chromatography/HPLC) increases solvent, resin and labour costs")
                    if raw_cost == "high":
                        inferred.append("High-cost starting materials increase COGS")
                    # fallback to structural reason
                    if not inferred:
                        if properties.get("complexity") == "high":
                            inferred.append("High molecular complexity increases synthesis and purification burden")
                        else:
                            inferred.append("Structural features increase downstream purification time and solvent use")
                    # limit to 3
                    costDrivers = inferred[:3]
            except Exception:
                pass

            # Add molecule-type specific suggested improvements and downstream insights
            downstreamInsights: List[str] = []
            try:
                # Small molecules: choose separation mode based on subtype and structural signals
                if molecule_type == "small_molecule":
                    if molecule_subtype == "volatile":
                        # distillation favorable for volatiles
                        if "Leverage fractional distillation for volatile separation" not in improvements:
                            improvements.append("Leverage fractional distillation for volatile separation to reduce solvent usage and enable continuous recovery")
                        downstreamInsights.append("Distillation likely feasible due to volatility; consider azeotrope management and reflux optimization")
                    else:
                        # non-volatile: crystallization or solvent-based
                        if properties.get("crystallization_potential"):
                            if "Implement crystallization protocols targeting solid form selection" not in improvements:
                                improvements.append("Implement crystallization protocols targeting solid form selection (solvent screening, seeding) to reduce chromatography load")
                            downstreamInsights.append("Crystallization or solvent-based purification expected due to aromatic/rigid structure")
                        else:
                            if "Design solvent-based extraction and flash-chromatography workflows" not in improvements:
                                improvements.append("Design solvent-based extraction and flash-chromatography workflows with solvent minimization strategies")
                            downstreamInsights.append("Crystallization or solvent-based purification expected for non-volatile small molecules")

                # Peptides: typically require HPLC; advise preparative methods
                elif molecule_type == "peptide":
                    if "Develop preparative reversed-phase HPLC methods and scale transfer plans" not in improvements:
                        improvements.append("Develop preparative reversed-phase HPLC methods and scale-transfer plans (optimize gradient, column loading, and solvent recovery)")
                    downstreamInsights.append("Peptide downstreams often require preparative HPLC (reverse-phase) and orthogonal polishing steps; plan for solvent recycling")

                # Proteins: chromatography centric, antibodies need affinity steps
                elif molecule_type == "protein":
                    if molecule_subtype == "antibody":
                        if "Include affinity capture (Protein A) followed by polishing chromatography in process train" not in improvements:
                            improvements.append("Include affinity capture (e.g., Protein A) followed by polishing chromatography (ion exchange, HIC) to reach purity targets; expect high resin and buffer costs")
                        downstreamInsights.append("Antibody downstream dominated by affinity capture and multi-step chromatography; plan for resin reuse and CIP cycles")
                    else:
                        if "Design multi-modal chromatography train for enzyme purification (capture + polish)" not in improvements:
                            improvements.append("Design multi-modal chromatography train for enzyme purification (capture + polish) and assess feedstock clarity to minimise column fouling")
                        downstreamInsights.append("Protein downstreams require chromatography-based capture and polishing; ensure robust column regeneration and buffer systems")

                # Natural products: mixtures requiring selective separation
                elif molecule_type == "natural_product":
                    if "Develop selective extraction and fractionation followed by targeted chromatography" not in improvements:
                        improvements.append("Develop selective extraction and fractionation followed by targeted chromatography and consider TR/GC for terpenes")
                    downstreamInsights.append("Natural products often require selective extraction, fractionation and targeted chromatographic separation due to complex mixtures")
            except Exception:
                pass

            # Include any subtype-specific notes collected earlier
            try:
                if subtype_notes:
                    downstreamInsights.extend([n for n in subtype_notes if n])
            except Exception:
                pass

            # Step 10: Generate improvements and attach potential cost impact
            # Estimate savings potential (low/medium/high)
            savingsPotential = "low"
            if cost == "high":
                savingsPotential = "high"
            if steps > 5:
                savingsPotential = "high"
            if desired_purity in (">99%", "very high"):
                # very high purity makes savings harder (but still medium)
                savingsPotential = savingsPotential if savingsPotential == "high" else "medium"

            def impact_label(level: str) -> str:
                return str(level).upper()

            if steps > 5:
                improvements.append(f"Reduce number of synthesis steps to improve efficiency and lower cost → Potential cost reduction: {impact_label('HIGH')}")

            # STEP 14: Improvements tailored to inputs
            try:
                if desired_purity in ("high", ">99%", "very high"):
                    if "Develop targeted purification train (crystallization, HPLC, polishing) to meet purity with lowest cost" not in improvements:
                        improvements.append("Develop targeted purification train (crystallization, HPLC, polishing) to meet purity with lowest cost")
                if scale == "industrial":
                    if "Simplify process design for improved scalability" not in improvements:
                        improvements.append("Simplify process design for improved scalability")
            except Exception:
                pass

            if properties.get("purification_difficulty") == "high" and desired_purity != "standard":
                improvements.append(f"Consider alternative purification methods such as crystallization → Potential cost reduction: {impact_label('MEDIUM')}")

            if cost in ("high", "very high"):
                improvements.append(f"Evaluate alternative suppliers for high-cost precursors and reduce protecting-group use to lower COGS → Potential cost reduction: {impact_label('MEDIUM')}")

            if risk == "high":
                improvements.append(f"Implement process control improvements (robust pH/temperature control, inline monitoring) and simplify critical unit ops to reduce failure modes → Potential cost reduction: {impact_label('MEDIUM')}")

            if strict_waste:
                improvements.append(f"Switch to less toxic reagents or reduce hazardous waste streams → Potential cost reduction: {impact_label('MEDIUM')}")

            # --- Trade-off engine: combine structural properties with process inputs ---
            # Start with adjusted defaults
            adjustedCost = cost
            adjustedRisk = risk
            adjustedPurification = properties.get("purification_difficulty", "medium")

            tradeoffs: List[str] = []

            # Purity × Purification
            if desired_purity not in ("standard", "") and properties.get("purification_difficulty") == "high":
                adjustedCost = "very high"
                adjustedRisk = "high"
                if "High purity requirement combined with difficult purification significantly increases cost" not in issues:
                    issues.append("High purity requirement combined with difficult purification significantly increases cost")
                tradeoffs.append("High purity requirement combined with difficult purification significantly increases cost and risk")

            # Steps × Complexity
            if steps > 5 and properties.get("complexity") == "high":
                adjustedCost = "very high"
                adjustedRisk = "high"
                if "Complex molecule combined with many synthesis steps leads to inefficient process" not in issues:
                    issues.append("Complex molecule combined with many synthesis steps leads to inefficient process")
                tradeoffs.append("Complex molecule plus many steps reduces efficiency and increases cost")

            # Scale × Stability
            if scale == "industrial" and properties.get("stability") == "low":
                adjustedRisk = "high"
                if "Low molecular stability may cause major issues during industrial scale-up" not in issues:
                    issues.append("Low molecular stability may cause major issues during industrial scale-up")
                tradeoffs.append("Low molecular stability increases scale-up risk and may require additional controls")

            # Scale × Steps
            if scale == "industrial" and steps > 5:
                adjustedCost = "very high"
                if "High number of steps is not suitable for industrial scale" not in issues:
                    issues.append("High number of steps is not suitable for industrial scale")
                tradeoffs.append("High step count reduces scalability and increases cost at industrial scale")

            # Toxicity × Waste constraints
            if properties.get("toxicity") == "high" and strict_waste:
                adjustedCost = "very high"
                adjustedRisk = "high"
                if "Toxic structure combined with strict waste constraints increases cost and complexity" not in issues:
                    issues.append("Toxic structure combined with strict waste constraints increases cost and complexity")
                tradeoffs.append("High inherent toxicity combined with strict waste rules raises compliance and disposal cost")

            # Raw material × Complexity
            if raw_cost == "high" and properties.get("complexity") == "high":
                adjustedCost = "very high"
                if "Complex molecule likely requires expensive raw materials" not in issues:
                    issues.append("Complex molecule likely requires expensive raw materials")
                tradeoffs.append("High raw material cost combined with structural complexity drives unit economics worsely")

            # Infrastructure × Requirements
            if desired_purity in ("very high", ">99%") and not has_advanced_pur:
                adjustedRisk = "high"
                if "Purity target cannot be efficiently achieved with current infrastructure" not in issues:
                    issues.append("Purity target cannot be efficiently achieved with current infrastructure")
                tradeoffs.append("Purity targets above standard require advanced purification assets or increase risk/cost")

            # Add some generic tradeoffs based on inputs
            if desired_purity not in ("standard", ""):
                tradeoffs.append("Higher purity targets increase purification cost and process complexity")
            if steps > 5:
                tradeoffs.append("Reducing steps may impact yield or product quality and requires route redesign")
            if properties.get("complexity") == "high":
                tradeoffs.append("Simplifying structure or route may reduce cost but requires alternative synthesis strategies")

            # Apply adjusted cost/risk back to main variables for final reporting
            cost = adjustedCost
            risk = adjustedRisk

            # Numeric proxies
            eff_score_map = {"low": 0.25, "medium": 0.5, "high": 0.85}
            cost_score_map = {"low": 0.25, "medium": 0.5, "high": 0.85, "very high": 0.95}
            # Deduplicate improvements and issues while preserving order
            seen_impr = set()
            dedup_impr: List[str] = []
            for im in improvements:
                if im not in seen_impr:
                    dedup_impr.append(im)
                    seen_impr.add(im)

            seen_iss = set()
            dedup_iss: List[str] = []
            for it in issues:
                if it not in seen_iss:
                    dedup_iss.append(it)
                    seen_iss.add(it)

            # Post-process textual outputs: remove redundancy, force specificity, and map SMILES/process signals to concrete language
            try:
                # topics for coarse semantic deduplication
                _topics = ["crystalliz", "hydrolys", "aggregation", "distill", "affinity", "convergent", "protect", "hplc", "solvent", "resin", "scale-up", "mixing", "stability"]
                def _topic_of(s: str):
                    ls = s.lower()
                    for t in _topics:
                        if t in ls:
                            return t
                    return None

                def _merge_unique(lst: List[str]) -> List[str]:
                    out: List[str] = []
                    seen_topics: List[str] = []
                    for s in lst:
                        t = _topic_of(s) or s.lower()
                        if t in seen_topics:
                            # if a similar-topic message already exists, prefer the more specific one (longer)
                            existing_idx = next((i for i, e in enumerate(out) if (_topic_of(e) or e.lower()) == t), None)
                            if existing_idx is not None:
                                # replace if current string is more specific
                                if len(s) > len(out[existing_idx]):
                                    out[existing_idx] = s
                            continue
                        out.append(s)
                        seen_topics.append(t)
                    return out

                dedup_iss = _merge_unique(dedup_iss)
                dedup_impr = _merge_unique(dedup_impr)

                # Remove hedging language and force more declarative, specific phrasing
                def _sharpen_text(s: str) -> str:
                    t = s
                    replacements = [
                        (" may involve ", " includes "),
                        (" may be ", " is likely to be "),
                        (" may cause ", " can cause "),
                        (" may increase ", " increases "),
                        (" may affect ", " affects "),
                        (" can be ", " is "),
                        (" may ", " is likely to "),
                    ]
                    for a, b in replacements:
                        try:
                            t = t.replace(a, b)
                        except Exception:
                            pass
                    # trim excessive whitespace
                    return " ".join(t.split())

                dedup_iss = [ _sharpen_text(x) for x in dedup_iss ]
                dedup_impr = [ _sharpen_text(x) for x in dedup_impr ]
            except Exception:
                pass

            # Build cause→effect explanations for each high-level metric
            try:
                eff_reasons = []
                if steps > 3:
                    eff_reasons.append(f"{steps} synthesis steps increase unit operations and cycle time")
                if properties.get("complexity") == "high":
                    eff_reasons.append("structural complexity increases handling and purification time")
                if method in ("extraction", "extract"):
                    eff_reasons.append("extraction workflows can reduce yield/throughput")
                if isBiomolecule and str(p.get("aggregation_risk") or "").lower() == "high":
                    eff_reasons.append("high aggregation risk reduces downstream throughput and increases batch failures")
                efficiency_explanation = "; ".join(eff_reasons) if eff_reasons else "Efficiency driven by moderate step count and tractable purification requirements"

                cost_reasons = []
                # costDrivers is already a list of descriptive drivers
                for cd in costDrivers:
                    cost_reasons.append(cd)
                if properties.get("purification_difficulty") == "high":
                    cost_reasons.append("structurally-driven purification (chromatography/HPLC) increases solvent and resin costs")
                if raw_cost == "high":
                    cost_reasons.append("expensive starting materials increase COGS")
                cost_explanation = "; ".join(cost_reasons) if cost_reasons else "Cost driven by standard materials and moderate downstream operations"

                risk_reasons = []
                if properties.get("stability") == "low":
                    risk_reasons.append("low molecular stability increases degradation risk during processing and storage")
                if isBiomolecule and str(p.get("folding_complexity") or "").lower() == "high":
                    risk_reasons.append("complex folding increases sensitivity to process parameters and risk of misfolding")
                if strict_waste:
                    risk_reasons.append("tight waste regulations increase compliance risk and process constraints")
                if scale == "industrial" and steps > 5:
                    risk_reasons.append("industrial scale-up with many steps multiplies operational risk")
                risk_explanation = "; ".join(risk_reasons) if risk_reasons else "Risk is manageable given stable structure and controllable process parameters"
            except Exception:
                efficiency_explanation = "Efficiency explanation not available"
                cost_explanation = "Cost explanation not available"
                risk_explanation = "Risk explanation not available"

            # Ensure a basic weighting_explanation exists for transparency
            if 'weighting_explanation' not in locals():
                try:
                    if 'weights' not in locals():
                        wsm = getSmilesWeight(molecule_type)
                        weights_local = {"smiles": wsm, "process": 1.0, "economic": 1.2, "downstream": 1.1, "biophysical": 1.3 if molecule_type in ("protein","peptide") else 0.0}
                        weighting_explanation = f"weights: smiles={weights_local['smiles']}, process={weights_local['process']}, economic={weights_local['economic']}, downstream={weights_local['downstream']}, biophysical={weights_local['biophysical']}"
                    else:
                        weighting_explanation = f"weights: smiles={weights.get('smiles')}, process={weights.get('process')}, economic={weights.get('economic')}, downstream={weights.get('downstream')}, biophysical={weights.get('biophysical')}"
                except Exception:
                    weighting_explanation = "weights unavailable"

            # Final pass: ensure structure_insight reflects SMILES/features and is molecule-type aware
            try:
                # For biomolecules prefer biophysical insights; for small molecules rely on SMILES-derived signals
                _insights_final = []
                if molecule_type in ('peptide', 'protein'):
                    # Prioritize aggregation/folding/stability signals
                    agg = str(p.get('aggregation_risk') or '').lower()
                    fold = str(p.get('folding_complexity') or '').lower()
                    bst = str(p.get('biophysical_stability') or '').lower()
                    if agg:
                        _insights_final.append(f'Aggregation risk={agg} — aggregation reduces yield and increases polishing/formulation cost')
                    if fold:
                        _insights_final.append(f'Folding complexity={fold} — complex folding increases refolding and process development cost')
                    if bst:
                        _insights_final.append(f'Biophysical stability={bst} — low stability increases degradation risk and lowers batch recovery')
                else:
                    try:
                        arom_ct = int(properties.get('aromatic_rings', 0) or 0)
                    except Exception:
                        arom_ct = 0
                    if arom_ct > 0:
                        _insights_final.append(f'Aromatic rings detected (count={arom_ct}) — supports crystallization screening to lower HPLC demand')
                    if properties.get('long_aliphatic_chain'):
                        _insights_final.append('Long hydrophobic chain increases solvent demand and suggests flash chromatography or phase-separation strategies')
                    if properties.get('ester_present'):
                        _insights_final.append('Ester functionality present — hydrolysis sensitivity requires controlled pH during workup')
                    if properties.get('amide_present'):
                        # avoid highlighting trivial amide motifs for biomolecules
                        if molecule_type not in ('peptide', 'protein'):
                            _insights_final.append('Amide groups present — monitor stereochemical integrity during prolonged processing')
                    try:
                        polar_ct = int(properties.get('polar_group_count', 0) or 0)
                    except Exception:
                        polar_ct = 0
                    if polar_ct >= 3:
                        _insights_final.append('Multiple polar groups suggest reversed-phase HPLC or polarity-based partitioning for purification')
                structure_insight_final = '; '.join(_insights_final) if _insights_final else (structure_insight if 'structure_insight' in locals() else None)
            except Exception:
                structure_insight_final = structure_insight if 'structure_insight' in locals() else None

            # --- Final consistency and policy enforcement (textual outputs only) ---
            try:
                # Replace remaining generic phrases with cause->effect statements
                def _remove_generic_phrases(lst: List[str]) -> List[str]:
                    out: List[str] = []
                    for s in lst:
                        t = s
                        # targeted replacements to remove vague language
                        replacements = [
                            ("requires careful control of biological systems", "has aggregation or contamination risks that reduce yield and complicate purification"),
                            ("requires careful control", "has process-sensitive parameters that, if uncontrolled, increase degradation and batch failure rates"),
                            ("may introduce challenges", "increases variability and failure modes during scale-up"),
                            ("introduces operational and mixing challenges", "increases residence-time variability and mixing-dependent degradation during scale-up"),
                        ]
                        for a, b in replacements:
                            try:
                                if a in t:
                                    t = t.replace(a, b)
                            except Exception:
                                pass
                        out.append(t)
                    return out

                dedup_iss = _remove_generic_phrases(dedup_iss)
                dedup_impr = _remove_generic_phrases(dedup_impr)
            except Exception:
                pass

            # Ensure cost driver rule: if cost level is medium or higher, provide 1-3 concrete, specific drivers
            try:
                costDrivers_final = costDrivers if isinstance(costDrivers, list) else []
                # interpret cost levels: include 'medium' and above
                if str(cost).lower() in ("medium", "high", "very high"):
                    # remove generic/SMILES-driven drivers for biomolecules unless biophysical drivers are absent
                    if molecule_type in ("peptide", "protein"):
                        biophys_drivers: List[str] = []
                        agg = str(p.get("aggregation_risk") or "").lower()
                        fold = str(p.get("folding_complexity") or "").lower()
                        bst = str(p.get("biophysical_stability") or "").lower()
                        if agg == "high":
                            biophys_drivers.append("High aggregation risk increases polishing chromatography and formulation costs")
                        if bst == "low":
                            biophys_drivers.append("Low biophysical stability increases degradation risk and rework cost")
                        if fold == "high":
                            biophys_drivers.append("Complex folding increases process development and refolding operation costs")
                        # prefer biophysical drivers if present
                        if biophys_drivers:
                            costDrivers_final = biophys_drivers[:3]
                    # otherwise construct drivers from process/structure
                    if not costDrivers_final:
                        inferred = []
                        if steps and steps > 2:
                            inferred.append(f"Multiple synthesis steps ({steps}) increase reagent, labour and unit-operation costs")
                        if properties.get("purification_difficulty") in ("high", "very high"):
                            inferred.append("Difficult purification (chromatography/HPLC) increases solvent, resin and labour costs")
                        if raw_cost == "high":
                            inferred.append("High-cost starting materials increase COGS and raw material procurement expenses")
                        # structural fallback (avoid SMILES minutiae for peptides/proteins)
                        if not inferred:
                            if properties.get("complexity") == "high":
                                inferred.append("High molecular complexity increases synthesis steps and purification burden")
                            else:
                                inferred.append("Structural features increase downstream purification time and solvent consumption")
                        costDrivers_final = inferred[:3]
                else:
                    costDrivers_final = costDrivers if isinstance(costDrivers, list) else []
            except Exception:
                costDrivers_final = costDrivers if isinstance(costDrivers, list) else []

            # replace cost drivers in final output
            try:
                costDrivers = costDrivers_final
            except Exception:
                pass

            # For peptides/proteins: downweight SMILES-derived text and prioritize biophysical cause->effect phrasing
            try:
                if molecule_type in ("peptide", "protein"):
                    # filter out trivial SMILES-driven items from dedup_iss and dedup_impr
                    def _is_trivial_smiles_item(s: str) -> bool:
                        low = s.lower()
                        trivial_markers = ["amide", "ester", "aromatic rings detected", "crystallization"]
                        return any(m in low for m in trivial_markers)

                    dedup_iss = [s for s in dedup_iss if not _is_trivial_smiles_item(s)]
                    # ensure biophysical issues are present and specific
                    biophys_add: List[str] = []
                    agg = str(p.get("aggregation_risk") or "").lower()
                    fold = str(p.get("folding_complexity") or "").lower()
                    bst = str(p.get("biophysical_stability") or "").lower()
                    if agg == "high":
                        biophys_add.append("High aggregation risk reduces yield and complicates downstream polishing and filtration")
                    if bst == "low":
                        biophys_add.append("Low biophysical stability increases degradation during hold times and lowers batch recovery")
                    if fold == "high":
                        biophys_add.append("Complex folding increases sensitivity to process conditions and raises refolding/workup costs")
                    # prepend these so they are prioritized
                    dedup_iss = biophys_add + dedup_iss
            except Exception:
                pass

            # Final refinement: remove generic fallbacks, merge related items, enforce cost-driver counts,
            # align wording to molecule type, and make trade-offs and recommendations concrete.
            try:
                # Helper: topic mapping for merging related items
                topic_map = {
                    'steps': ['step', 'synthesis', 'multi-step', 'step count'],
                    'purification': ['purificat', 'hplc', 'chromatograph', 'polish', 'distillation', 'crystalliz'],
                    'raw_materials': ['raw material', 'starting material', 'precursor', 'cog'],
                    'aggregation': ['aggregat', 'aggregate'],
                    'stability': ['stabil', 'degrad'],
                    'folding': ['fold', 'refold'],
                    'solvent': ['solvent', 'phase-separat', 'partition'],
                    'scale': ['scale', 'industrial', 'pilot', 'mixing', 'residence-time'],
                    'ester': ['ester'],
                    'aromatic': ['aromatic', 'crystalliz'],
                }

                def _topic_of_text(s: str) -> Optional[str]:
                    low = s.lower()
                    for t, kws in topic_map.items():
                        for kw in kws:
                            if kw in low:
                                return t
                    return None

                # Merge related messages into one strong statement per topic
                def _merge_by_topic(lst: List[str]) -> List[str]:
                    grouped: Dict[str, List[str]] = {}
                    others: List[str] = []
                    for s in lst:
                        t = _topic_of_text(s)
                        if t:
                            grouped.setdefault(t, []).append(s)
                        else:
                            others.append(s)
                    merged: List[str] = []
                    for t, items in grouped.items():
                        # create a consolidated sentence depending on topic
                        if t == 'steps':
                            merged.append('High synthesis step count increases operational complexity, extends cycle time, and limits scalability')
                        elif t == 'purification':
                            merged.append('Difficult purification (HPLC/chromatography or distillation) increases solvent, buffer and resin costs and extends downstream cycle time')
                        elif t == 'raw_materials':
                            merged.append('Expensive or scarce starting materials increase COGS and procurement risk')
                        elif t == 'aggregation':
                            merged.append('High aggregation risk reduces yield and increases polishing and formulation costs')
                        elif t == 'stability':
                            merged.append('Low biophysical stability increases degradation during processing and reduces batch recovery')
                        elif t == 'folding':
                            merged.append('Complex folding increases sensitivity to process conditions and raises development and refolding costs')
                        elif t == 'solvent':
                            merged.append('Hydrophobicity and solvent demand increase solvent consumption and complicate phase separations')
                        elif t == 'scale':
                            merged.append('Industrial scale-up increases mixing and residence-time variability that can drive degradation and batch failures')
                        elif t == 'ester':
                            merged.append('Ester groups introduce hydrolysis sensitivity under variable pH, increasing degradation risk during workup')
                        elif t == 'aromatic':
                            merged.append('Aromatic rings enable crystallization-based purification opportunities that can reduce HPLC burden')
                        else:
                            # fallback: prefer the longest, most specific item
                            merged.append(max(items, key=len))
                    # append other, non-topic items but avoid duplicates
                    for s in others:
                        if s not in merged:
                            merged.append(s)
                    return merged

                dedup_iss = _merge_by_topic(dedup_iss)
                dedup_impr = _merge_by_topic(dedup_impr)

                # Align messages to molecule type: remove domain-mismatched suggestions
                def _filter_by_molecule_type(lst: List[str], mtype: str) -> List[str]:
                    out: List[str] = []
                    for s in lst:
                        low = s.lower()
                        if mtype in ('peptide', 'protein'):
                            # drop small-molecule specific actions
                            if any(k in low for k in ('distillation', 'reflux', 'entrainer', 'fractional distillation')):
                                continue
                        if mtype == 'small_molecule':
                            # drop biologic-specific suggestions
                            if any(k in low for k in ('protein a', 'affinity capture', 'bioreactor', 'resin reuse', 'cip')):
                                continue
                        out.append(s)
                    return out

                dedup_iss = _filter_by_molecule_type(dedup_iss, molecule_type)
                dedup_impr = _filter_by_molecule_type(dedup_impr, molecule_type)

                # Strengthen SMILES-based insights: ensure structural features are expressed as process implications
                smi_insights: List[str] = []
                try:
                    # For biomolecules, avoid relying on SMILES minutiae; prioritize biophysical inputs
                    if molecule_type not in ('peptide', 'protein'):
                        try:
                            arom_ct = int(properties.get('aromatic_rings', 0) or 0)
                        except Exception:
                            arom_ct = 0
                        if arom_ct > 0:
                            smi_insights.append(f'Aromatic rings detected (count={arom_ct}) — enable crystallization screening to lower HPLC demand')
                        if properties.get('ester_present'):
                            smi_insights.append('Ester functionality present — hydrolysis sensitivity requires controlled pH during workup to avoid degradation')
                        if properties.get('long_aliphatic_chain'):
                            smi_insights.append('Long hydrophobic chain increases solvent consumption and complicates aqueous workups')
                        try:
                            polar_ct = int(properties.get('polar_group_count', 0) or 0)
                        except Exception:
                            polar_ct = 0
                        if polar_ct >= 3:
                            smi_insights.append('Multiple polar groups suggest reversed-phase HPLC or polarity-based partitioning for purification')
                    # Prepend SMILES insights if not redundant
                    for si in reversed(smi_insights):
                        if si not in dedup_iss:
                            dedup_iss.insert(0, si)
                except Exception:
                    pass

                # Enforce cost driver counts: HIGH/VERY HIGH -> 2-3 drivers, MEDIUM -> at least 1
                try:
                    costDrivers_list = costDrivers if isinstance(costDrivers, list) else ([] if not costDrivers else [costDrivers])
                    desired_min = 1 if str(cost).lower() == 'medium' else 2 if str(cost).lower() in ('high', 'very high') else 0
                    if len(costDrivers_list) < desired_min:
                        inf: List[str] = []
                        if steps and steps > 2:
                            inf.append(f'Multiple synthesis steps ({steps}) increase reagent, labour and unit-operation costs')
                        if properties.get('purification_difficulty') in ('high', 'very high'):
                            inf.append('Difficult purification (chromatography/HPLC) increases solvent, resin and labour costs')
                        if raw_cost == 'high':
                            inf.append('High-cost starting materials increase COGS and procurement expenses')
                        if efficiencyScore < 5:
                            inf.append('Low process efficiency increases unit cost through yield loss and rework')
                        # for biomolecules prefer biophysical drivers
                        if molecule_type in ('peptide', 'protein'):
                            agg = str(p.get('aggregation_risk') or '').lower()
                            bst = str(p.get('biophysical_stability') or '').lower()
                            if agg == 'high' and 'High aggregation risk' not in inf:
                                inf.insert(0, 'High aggregation risk increases polishing and formulation costs')
                            if bst == 'low' and 'Low biophysical stability' not in inf:
                                inf.insert(0, 'Low biophysical stability increases degradation risk and rework cost')
                        # fill up to 3 drivers
                        for x in inf:
                            if x not in costDrivers_list:
                                costDrivers_list.append(x)
                            if len(costDrivers_list) >= 3:
                                break
                    # final trim to 3
                    costDrivers = costDrivers_list[:3]
                except Exception:
                    pass

                # Improve tradeoffs: make explicit cause->effect tradeoffs
                try:
                    tradeoffs_refined: List[str] = []
                    for tr in tradeoffs:
                        t = tr.lower()
                        if 'purity' in t and 'difficult' in t:
                            tradeoffs_refined.append('Very high purity targets require extended polishing (HPLC/chromatography), increasing solvent/buffer consumption and cycle time')
                        elif 'steps' in t and 'complex' in t:
                            tradeoffs_refined.append('Reducing synthesis steps lowers operational cost and cycle time but may require novel reagents or convergent routes that increase raw material expenses')
                        elif 'stability' in t and 'industrial' in t:
                            tradeoffs_refined.append('Low molecular stability at industrial scale increases degradation risk, requiring tighter control and chilled holds which raise CAPEX/OPEX')
                        else:
                            # keep original but ensure it is concise
                            tradeoffs_refined.append(tr)
                    tradeoffs = tradeoffs_refined[:5]
                except Exception:
                    pass

                # Make recommendations actionable and concise: ensure each improvement ties to a problem and an effect
                def _actionable(impr: str) -> str:
                    s = impr
                    # common patterns and stronger rewrites
                    if 'crystallization' in s.lower():
                        return 'Develop seeded crystallization and polymorph screening to reduce HPLC load and solvent use'
                    if 'prepare' in s.lower() and 'hplc' in s.lower():
                        return 'Develop preparative RP-HPLC method with column sizing and solvent recycle to enable scale-up'
                    if 'affinity capture' in s.lower() or 'protein a' in s.lower():
                        return 'Implement affinity capture (Protein A) followed by IEX/HIC polishing to achieve purity while controlling resin costs via reuse and CIP'
                    if 'aggregation' in s.lower():
                        return 'Introduce formulation screens, low-shear handling and inline aggregate monitoring to reduce aggregate-related yield loss'
                    if 'refold' in s.lower() or 'fold' in s.lower():
                        return 'Add controlled refolding steps with gradual denaturant removal and defined acceptance criteria to improve yield'
                    # fallback: trim and keep
                    return ' '.join(s.split())

                dedup_impr = [_actionable(x) for x in dedup_impr]

                # Remove improvements that merely restate issues
                try:
                    filtered_impr: List[str] = []
                    for imp in dedup_impr:
                        imp_low = imp.lower()
                        is_redundant = any(imp_low in iss.lower() or iss.lower() in imp_low for iss in dedup_iss)
                        if not is_redundant:
                            filtered_impr.append(imp)
                    dedup_impr = filtered_impr
                except Exception:
                    pass

                # Remove weak language: convert 'may/could/potentially' to stronger phrasing where safe
                def _strengthen_lang(lst: List[str]) -> List[str]:
                    out: List[str] = []
                    for s in lst:
                        t = s.replace(' may ', ' is likely to ').replace(' could ', ' can ').replace(' potentially ', ' ')
                        # remove duplicate whitespace
                        t = ' '.join(t.split())
                        out.append(t)
                    return out

                dedup_iss = _strengthen_lang(dedup_iss)
                dedup_impr = _strengthen_lang(dedup_impr)

                # Keep outputs concise: cap counts
                dedup_iss = dedup_iss[:8]
                dedup_impr = dedup_impr[:8]
                costDrivers = costDrivers[:3]
            except Exception:
                # leave originals if anything fails
                pass

            # Replace remaining specific generic biotech phrase with explicit cause→effect
            try:
                dedup_iss = [
                    ("Biotechnological production increases contamination and aggregation risk, requiring sterile bioreactor operations and validated media which raise operational cost" if (s.lower().strip().startswith('biotech processes require careful control') or s.lower().strip().startswith('biotech process require careful control')) else s)
                    for s in dedup_iss
                ]
            except Exception:
                pass

            # Normalize improvement phrasing (remove arrows and make short)
            try:
                dedup_impr = [s.replace('→', '-').replace('Potential cost reduction', 'Expected cost reduction') for s in dedup_impr]
            except Exception:
                pass

            analysis = {
                "efficiency": efficiency,
                "efficiency_score": eff_score_map.get(efficiency, 0.5),
                "cost": cost,
                "costLevel": cost,
                "costScore": costScore,
                "cost_score": cost_score_map.get(cost, 0.5),
                "costDrivers": costDrivers,
                "savingsPotential": savingsPotential,
                "risk": risk,
                "toxicity": toxicity,
                "final_toxicity": toxicity,
                "issues": dedup_iss,
                "improvements": dedup_impr,
                "efficiency_explanation": efficiency_explanation,
                "cost_explanation": costDrivers_new if 'costDrivers_new' in locals() else cost_explanation,
                "score_breakdown": score_breakdown if 'score_breakdown' in locals() else {"smiles":0.0,"process":0.0,"economic":0.0,"downstream":0.0,"biophysical":0.0},
                "risk_breakdown": risk_breakdown if 'risk_breakdown' in locals() else {"smiles":0.0,"process":0.0,"economic":0.0,"downstream":0.0,"biophysical":0.0},
                "structure_insight": structure_insight_final,
                "downstream_detail": downstream_detail if 'downstream_detail' in locals() else None,
                "risk_explanation": risk_explanation,
                # numeric raw scores and weighting transparency
                "numeric_scores": {
                    "costScore": costScore,
                    "riskScore": riskScore,
                    "efficiencyScore": efficiencyScore,
                },
                "weighting_explanation": weighting_explanation if 'weighting_explanation' in locals() else None,
                "properties": properties,
                "process_class": process_class if 'process_class' in locals() else None,
                "process_class_explanation": process_class_explanation if 'process_class_explanation' in locals() else None,
                "tradeoffs": tradeoffs,
                "downstream_insights": downstreamInsights,
                "molecule_type": molecule_type,
            }

            return analysis
        except Exception:
            logging.exception("analyzeProcess failed")
            return {
                "efficiency": "medium",
                "efficiency_score": 0.5,
                "cost": "medium",
                "costLevel": "medium",
                "costScore": 0,
                "cost_score": 0.5,
                "costDrivers": [],
                "savingsPotential": "low",
                "risk": "medium",
                "toxicity": "medium",
                "final_toxicity": "medium",
                "issues": [],
                "improvements": [],
                "properties": properties if properties is not None else {},
                "downstream_insights": [],
                "molecule_type": str((inp or {}).get('molecule_type') or 'small_molecule').lower(),
            }