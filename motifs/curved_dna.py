"""
Category 1: Curved DNA Detection Module
======================================

This module implements detection algorithms for curved DNA structures,
including global and local poly(A)/poly(T) tracts that cause DNA bending.

Scientific Basis:
- A-tracts and T-tracts cause DNA bending through narrowed minor groove
- Curvature score reflects biological propensity for curvature
- Periodic spacing of A/T tracts enhances curvature

Author: Dr. Venkata Rajesh Yella
Updated: 2024
"""

import re
from .shared_utils import wrap, calculate_conservation_score

def find_polyA_polyT_tracts(seq: str, min_len: int = 7) -> list:
    """Find contiguous A or T tracts of minimum length."""
    results = []
    i = 0
    n = len(seq)
    while i < n:
        if seq[i] == 'A' or seq[i] == 'T':
            ch = seq[i]
            start = i
            while i < n and seq[i] == ch:
                i += 1
            if i - start >= min_len:
                results.append((start, i-1, seq[start:i]))
        else:
            i += 1
    return results

def curvature_score(seq):
    """
    Calculate curvature score based on AT-richness and tract organization.
    
    Scientific Basis: A-tracts and T-tracts cause DNA bending through 
    narrowed minor groove. Score reflects biological propensity for curvature.
    
    Parameters:
    seq (str): DNA sequence to score
    
    Returns:
    float: Curvature propensity score
    """
    # Raw score: length scaled by AT-bias and mild periodicity bonus based on A/T tracts
    if not seq:
        return 0.0
    at_frac = (seq.count('A') + seq.count('T')) / len(seq)
    # count segments of mono-base runs (A or T)
    runs = re.findall(r"(A+|T+)", seq)
    run_bonus = sum(len(r)**0.5 for r in runs)  # diminishing returns
    return len(seq) * (1.0 + at_frac) + 0.5 * run_bonus

def find_global_curved_polyA_polyT(seq: str, min_tract_len: int = 3, min_repeats: int = 3, min_spacing: int = 8, max_spacing: int = 12, min_score: int = 6) -> tuple:
    """Find globally curved DNA with periodic poly(A)/poly(T) tracts."""
    tracts = find_polyA_polyT_tracts(seq, min_tract_len)
    results = []
    apr_regions = []
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
            motif_seq = seq[group[0][0]:group[-1][1]+1]
            score = curvature_score(motif_seq)
            if score >= min_score:
                # Calculate conservation score
                conservation_result = calculate_conservation_score(motif_seq, "Curved DNA")
                conservation_score = conservation_result["enrichment_score"]
                
                motif = {
                    "Sequence Name": "",
                    "Class": "Curved DNA",
                    "Subtype": "Global curvature",
                    "Start": group[0][0] + 1,
                    "End": group[-1][1] + 1,
                    "Length": group[-1][1] - group[0][0] + 1,
                    "Sequence": wrap(motif_seq),
                    "Score": float(score),
                    "Conservation_Score": float(conservation_score),
                    "Conservation_P_Value": float(conservation_result["p_value"]),
                    "Conservation_Significance": conservation_result["significance"],
                    "Arms/Repeat Unit/Copies": "",
                    "Spacer": ""
                }
                results.append(motif)
                apr_regions.append((motif["Start"], motif["End"]))
    return results, apr_regions

def find_local_curved_polyA_polyT(seq: str, apr_regions: list, min_len: int = 7) -> list:
    """Find local curved DNA with isolated poly(A)/poly(T) tracts."""
    results = []
    tracts = find_polyA_polyT_tracts(seq, min_len)
    for start, end, tract_seq in tracts:
        s, e = start + 1, end + 1
        if not any(r_start <= s <= r_end or r_start <= e <= r_end for r_start, r_end in apr_regions):
            # Calculate conservation score
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
    
    Returns both global and local curved DNA structures.
    Enhanced scoring system focuses on:
    1. Curvature propensity based on A/T content and tract organization
    2. Conservation analysis using composition-preserving shuffles
    3. Distinction between global periodic and local isolated structures
    """
    global_results, apr_regions = find_global_curved_polyA_polyT(seq)
    local_results = find_local_curved_polyA_polyT(seq, apr_regions)
    return global_results + local_results