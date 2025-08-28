"""
Category 1: Curved DNA Detection Module
=======================================
This module implements detection algorithms for curved DNA structures, focusing on poly(A)/poly(T) tracts and their periodic phasing.

- Poly(dA:dT) and poly(T) tracts generate intrinsic curvature via minor-groove narrowing and hydration spine formation.
- Tracts ≥4–6 bp contribute cooperatively; longer tracts (≥7 bp) act as rigid motifs.
- In-phase (~10–11 bp) arrays of ≥3 tracts reinforce global curvature (phased A-tract curvature).
- Isolated long tracts (≥7 bp) confer local curvature if not overlapping with global arrays.

Author: Dr. Venkata Rajesh Yella
Updated: 2024 (with annotations)
"""

import re
from .shared_utils import wrap, calculate_conservation_score

# --- Utility scoring functions for global curvature arrays ---

def _phase_score(centers, target=10.5, tol=2.5):
    """
    Compute phasing score for tract centers.
    Peaks when tract spacing ≈ canonical DNA helical repeat (10–11 bp).
    """
    score = 0.0
    if len(centers) < 2:
        return 0.0
    for i in range(1, len(centers)):
        spacing = abs(centers[i] - centers[i-1])
        # Gaussian-like penalty away from target spacing
        score += max(0, 1.0 - ((spacing - target) / tol)**2)
    return score

def _length_step_bonus(region_len):
    """
    Length bonus for extended phased array regions:
    +1 for ≥50 bp, +2 for ≥100 bp.
    """
    return 1.0 * (region_len >= 50) + 1.0 * (region_len >= 100)

def _non_overlap_selection(motif_list):
    """
    Select non-overlapping motifs, prioritizing by descending Score.
    Standard greedy interval scheduling.
    """
    sorted_motifs = sorted(motif_list, key=lambda m: m['Score'], reverse=True)
    selected, covered = [], set()
    for motif in sorted_motifs:
        region = set(range(motif['Start'], motif['End']+1))
        if not region & covered:  # keep only if not overlapping
            selected.append(motif)
            covered |= region
    return selected

# --- Core detection functions ---

def find_polyA_polyT_tracts(seq: str, min_len: int = 7) -> list:
    """
    Detect contiguous A or T tracts ≥ min_len.
    Returns: list of (start, end, tract_seq)
    """
    results, i, n = [], 0, len(seq)
    while i < n:
        if seq[i] in ("A", "T"):
            ch, start = seq[i], i
            while i < n and seq[i] == ch:
                i += 1
            if i - start >= min_len:
                results.append((start, i-1, seq[start:i]))
        else:
            i += 1
    return results

def curvature_score(seq):
    """
    Heuristic curvature score for local tracts:
    Combines AT-content, run length, and tract size.
    """
    if not seq:
        return 0.0
    at_frac = (seq.count("A") + seq.count("T")) / len(seq)
    runs = re.findall(r"(A+|T+)", seq)
    run_bonus = sum(len(r)**0.5 for r in runs)
    return len(seq) * (1.0 + at_frac) + 0.5 * run_bonus

def find_global_curved_polyA_polyT(
    seq: str,
    min_tract_len: int = 4,   # allow shorter seeds, arrays require ≥3 tracts
    min_repeats: int = 3,     # 3–6 phased tracts → global curvature
    min_spacing: int = 6,
    max_spacing: int = 15,
    min_score: float = 4.0
):
    """
    Detect phased arrays of poly(A)/poly(T) tracts (global curvature).
    Criteria:
      - ≥3 tracts, each ≥4 bp
      - Phased at ~10–11 bp (± tolerance)
      - Arrays of 3–6 tracts → global curvature signal
    Returns: (motif_list, region_list)
    """
    tracts = find_polyA_polyT_tracts(seq, min_tract_len)
    results, apr_regions = [], []

    for i in range(len(tracts)):
        group = [tracts[i]]
        # Build phased array by chaining tracts with proper spacing
        for j in range(i+1, len(tracts)):
            prev_center = (group[-1][0] + group[-1][1]) // 2
            curr_center = (tracts[j][0] + tracts[j][1]) // 2
            spacing = curr_center - prev_center

            if min_spacing <= spacing <= max_spacing:
                group.append(tracts[j])
            elif len(group) >= min_repeats:
                break

        # Score arrays with ≥3 tracts
        if len(group) >= min_repeats:
            gstart, gend = group[0][0], group[-1][1]
            motif_seq = seq[gstart:gend+1]
            centers = [(s+e)//2 for s, e, _ in group]
            tract_len_sum = sum((e - s + 1) for s, e, _ in group)
            tract_count = len(group)
            phs = _phase_score(centers, target=10.5, tol=3.0)
            length_bonus = _length_step_bonus(gend - gstart + 1)

            score = (tract_len_sum / 8.0) + (1.5 * tract_count) + (3.0 * phs) + length_bonus

            if score >= min_score:
                conservation_result = calculate_conservation_score(motif_seq, "Curved DNA")
                conservation_score = conservation_result["enrichment_score"]
                motif = {
                    "Sequence Name": "",
                    "Class": "Curved DNA",
                    "Subtype": "Global curvature (phased A-tracts)",
                    "Start": gstart+1,
                    "End": gend+1,
                    "Length": gend - gstart + 1,
                    "Sequence": wrap(motif_seq),
                    "Score": float(score),
                    "Conservation_Score": float(conservation_score),
                    "Conservation_P_Value": float(conservation_result["p_value"]),
                    "Conservation_Significance": conservation_result["significance"],
                    "Arms/Repeat Unit/Copies": f"Tracts={tract_count}",
                    "Spacer": f"Spacing={min_spacing}-{max_spacing}"
                }
                results.append(motif)
                apr_regions.append((motif["Start"], motif["End"]))
    return results, apr_regions

def find_local_curved_polyA_polyT(seq: str, apr_regions: list, min_len: int = 7) -> list:
    """
    Detect isolated poly(A)/poly(T) tracts (local curvature).
    Criteria:
      - Tracts ≥7 bp (rigid motif threshold)
      - Exclude regions overlapping phased arrays
    Returns: motif list
    """
    results, tracts = [], find_polyA_polyT_tracts(seq, min_len)
    for start, end, tract_seq in tracts:
        s, e = start+1, end+1
        # Ensure no overlap with global regions
        if not any(r_start <= s <= r_end or r_start <= e <= r_end for r_start, r_end in apr_regions):
            conservation_result = calculate_conservation_score(tract_seq, "Curved DNA")
            conservation_score = conservation_result["enrichment_score"]
            results.append({
                "Sequence Name": "",
                "Class": "Curved DNA",
                "Subtype": "Local curvature (long A-tract)",
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
    Master function: detects both global (phased arrays) and local (isolated long tracts) curvature.
    Uses greedy selection to report non-overlapping high-confidence motifs.
    """
    global_results, apr_regions = find_global_curved_polyA_polyT(seq)
    local_results = find_local_curved_polyA_polyT(seq, apr_regions)
    return _non_overlap_selection(global_results + local_results)
