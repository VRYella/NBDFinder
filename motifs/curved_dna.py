"""
Category 1: Curved DNA Detection Module
=======================================
This module implements detection algorithms for curved DNA structures, focusing on poly(A)/poly(T) tracts and their periodic phasing, as supported by scientific literature.
- Poly(dA:dT) and poly(T) tracts generate intrinsic curvature via minor-groove narrowing and a hydration spine, with cooperative behavior that consolidates for tracts ≥4–6 bp and strengthens with length.
- In-phase (~10–11 bp) spacing of multiple tracts reinforces macroscopic curvature (see: Travers & Muskhelishvili, 2007, Nat Rev Microbiol; Segal & Widom, 2009, Nat Rev Genet).

Author: Dr. Venkata Rajesh Yella
Updated: 2024
"""

import re
from .shared_utils import wrap, calculate_conservation_score

# --- Utility scoring functions for global array scoring ---

def _phase_score(centers, target=10.5, tol=2.5):
    """
    Compute phasing score over tract centers.
    Peaks at target ± tol (default ~10.5±2.5 bp, canonical DNA helical repeat).
    """
    score = 0.0
    if len(centers) < 2: return 0.0
    for i in range(1, len(centers)):
        spacing = abs(centers[i] - centers[i-1])
        # Gaussian-like peak centered at target
        score += max(0, 1.0 - ((spacing - target) / tol)**2)
    return score

def _length_step_bonus(region_len):
    """
    Bonus for array region length. +1 for ≥50 bp, +2 for ≥100 bp.
    """
    return 1.0*(region_len >= 50) + 1.0*(region_len >= 100)

def _non_overlap_selection(motif_list):
    """
    Select non-overlapping intervals from a motif list by descending score.
    Standard greedy approach for interval selection.
    """
    sorted_motifs = sorted(motif_list, key=lambda m: m['Score'], reverse=True)
    selected = []
    covered = set()
    for motif in sorted_motifs:
        region = set(range(motif['Start'], motif['End']+1))
        if not region & covered:
            selected.append(motif)
            covered |= region
    return selected

def find_polyA_polyT_tracts(seq: str, min_len: int = 7) -> list:
    """
    Find contiguous A or T tracts of minimum length.
    Returns: list of (start, end, tract_seq)
    """
    results = []
    i = 0; n = len(seq)
    while i < n:
        if seq[i] == 'A' or seq[i] == 'T':
            ch = seq[i]; start = i
            while i < n and seq[i] == ch:
                i += 1
            if i - start >= min_len:
                results.append((start, i-1, seq[start:i]))
        else:
            i += 1
    return results

def curvature_score(seq):
    """
    Calculate curvature score for a DNA segment based on AT-richness and tract organization.
    - Heuristic for local tracts: length & AT-content.
    """
    if not seq: return 0.0
    at_frac = (seq.count('A') + seq.count('T')) / len(seq)
    runs = re.findall(r"(A+|T+)", seq)
    run_bonus = sum(len(r)**0.5 for r in runs)
    return len(seq) * (1.0 + at_frac) + 0.5 * run_bonus

def find_global_curved_polyA_polyT(seq: str, min_tract_len: int = 3, min_repeats: int = 3, min_spacing: int = 8, max_spacing: int = 12, min_score: float = 6.0):
    """
    Detect globally curved DNA with periodic, phased poly(A)/poly(T) tracts.
    Scientific basis: in-phase arrays of A/T tracts enforce curvature (phasing ~10–11 bp).
    Strategy:
      - Seed with A/T tracts (min_tract_len) to capture short contributors, but scoring favors longer tracts and better phasing.
      - Group consecutive tracts whose center spacings lie within [min_spacing, max_spacing].
      - Score: tract_load (sum of tract lengths), tract_count, phasing, length bonus.
    Returns: (list of motif dicts, list of region tuples)
    """
    tracts = find_polyA_polyT_tracts(seq, min_tract_len)
    results = []; apr_regions = []
    for i in range(len(tracts) - min_repeats + 1):
        group = [tracts[i]]
        for j in range(1, min_repeats):
            prev_center = (tracts[i + j - 1][0] + tracts[i + j - 1][1]) // 2
            curr_center = (tracts[i + j][0] + tracts[i + j][1]) // 2
            spacing = curr_center - prev_center
            if min_spacing <= spacing <= max_spacing:
                group.append(tracts[i + j])
            else:
                break
        if len(group) >= min_repeats:
            gstart = group[0][0]; gend = group[-1][1]
            motif_seq = seq[gstart:gend+1]
            centers = [((s+e)//2) for s, e, _ in group]
            tract_len_sum = sum((e - s + 1) for s, e, _ in group)
            tract_count = len(group)
            phs = _phase_score(centers, target=10.5, tol=2.5)
            length_bonus = _length_step_bonus(gend - gstart + 1)
            # Deterministic global score (no heavy composition dependency)
            score = (tract_len_sum / 10.0) + (1.0 * tract_count) + (2.0 * phs) + length_bonus
            if score >= min_score:
                conservation_result = calculate_conservation_score(motif_seq, "Curved DNA")
                conservation_score = conservation_result["enrichment_score"]
                motif = {
                    "Sequence Name": "",
                    "Class": "Curved DNA",
                    "Subtype": "Global curvature",
                    "Start": gstart + 1,
                    "End": gend + 1,
                    "Length": gend - gstart + 1,
                    "Sequence": wrap(motif_seq),
                    "Score": float(score),
                    "Conservation_Score": float(conservation_score),
                    "Conservation_P_Value": float(conservation_result["p_value"]),
                    "Conservation_Significance": conservation_result["significance"],
                    # Tract count and array spacing for interpretability
                    "Arms/Repeat Unit/Copies": f"Tracts={tract_count}",
                    "Spacer": f"Spacing={min_spacing}-{max_spacing}"
                }
                results.append(motif)
                apr_regions.append((motif["Start"], motif["End"]))
    return results, apr_regions

def find_local_curved_polyA_polyT(seq: str, apr_regions: list, min_len: int = 7) -> list:
    """
    Detect local curved DNA with isolated poly(A)/poly(T) tracts.
    Strategy:
      - Report tracts of length ≥min_len that are outside global APR regions.
      - Score emphasizes tract length (curvature_score).
    """
    results = []
    tracts = find_polyA_polyT_tracts(seq, min_len)
    for start, end, tract_seq in tracts:
        s, e = start + 1, end + 1
        if not any(r_start <= s <= r_end or r_start <= e <= r_end for r_start, r_end in apr_regions):
            conservation_result = calculate_conservation_score(tract_seq, "Curved DNA")
            conservation_score = conservation_result["enrichment_score"]
            results.append({
                "Sequence Name": "",
                "Class": "Curved DNA",
                "Subtype": "Local Curvature",
                "Start": s,
                "End": e,
                "Length": len(tract_seq),
                "Sequence": wrap(tract_seq),
                "Score": float(curvature_score(tract_seq)),
                "Conservation_Score": float(conservation_score),
                "Conservation_P_Value": float(conservation_result["p_value"]),
                "Conservation_Significance": conservation_result["significance"],
                "Arms/Repeat Unit/Copies": "",
                "Spacer": ""
            })
    return results

def find_curved_DNA(seq: str) -> list:
    """
    Main function to detect curved DNA structures.
    Returns both global (phased arrays) and local (isolated tracts) curvature calls.
    Enhancements:
      1. Phasing-aware scoring for global arrays (target ~10.5 bp).
      2. Deterministic length/tract-count terms for reproducibility.
      3. Final non-overlap selection across global and local calls for high-confidence intervals.
    """
    global_results, apr_regions = find_global_curved_polyA_polyT(seq)
    local_results = find_local_curved_polyA_polyT(seq, apr_regions)
    # Merge and retain only top-scoring, non-overlapping intervals
    return _non_overlap_selection(global_results + local_results)
