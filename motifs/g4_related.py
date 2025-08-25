"""
Category 8: G4-related Detection Module - Modernized with Exact G4Hunter
========================================================================

This module implements detection algorithms for all G-quadruplex related structures
using exact G4Hunter algorithm and literature-standard patterns from:
- QuadBase2
- G4Hunter (Bedrat et al. NAR 2016)  
- Comprehensive G4 classification framework

MODERNIZED G4 FAMILY DETECTION:
-------------------------------
- Uses exact G4Hunter scoring algorithm (not fast approximation)
- Implements 9 distinct G4 subclasses with literature-standard patterns
- Follows priority-based prediction order for overlapping motifs
- Non-overlapping selection based on priority, score, and structural factors

G4 SUBCLASSES IMPLEMENTED:
--------------------------
1. Multimeric/Tandem: two canonical G4s within 20 nt
2. Split/Bipartite: two G3+L1-7 blocks separated by long linker (20-50 nt)  
3. Bulged: G4s with single-base bulges in G-runs (up to 3 nt between Gs)
4. Canonical G3+L1–7: four G-runs (≥3 Gs), loops 1–7 nt
5. Long-Loop G3+L8–12: canonical G4 with at least one loop 8–12 nt
6. Extended G3+L1–12: canonical G4 with loops up to 12 nt
7. Two-Tetrad G2L1–12: four G-runs of 2 Gs, loops 1–12 nt
8. Generalized G2+L1–12: four G-runs (≥2 Gs), loops 1–12 nt
9. G-Triplex: three G-runs (≥3 Gs), loops up to 15 nt

SCIENTIFIC BASIS:
-----------------
    - G-quadruplexes are four-stranded structures formed by guanine-rich sequences
    - Various topologies based on loop length, G-run composition, and structural features
    - Formation potential based on exact G4Hunter scoring methodology
    - Priority-based classification prevents redundant overlapping predictions

REFERENCES:
-----------
    - Bedrat et al. (2016) NAR (G4Hunter algorithm)
    - Hänsel-Hertsch et al. (2017) Nat Rev Mol Cell Biol
    - Chambers et al. (2015) Nat Biotechnol (G4 experimental validation)
    - QuadBase2 database for comprehensive G4 classification

Author: Dr. Venkata Rajesh Yella
Updated: 2024 with modernized G4 family detection framework
"""

import re
from math import log2, sqrt

def mean(values): 
    """Mean utility for lists, guarded for empty lists."""
    return sum(values) / len(values) if values else 0

# Import shared utilities and the modernized G4 family module
from .shared_utils import (
    wrap, calculate_conservation_score, overlapping_finditer,
    calculate_structural_factor
)
from .g4_family_modernized import (
    candidates_from, non_overlapping, g4hunter_score as exact_g4hunter_score,
    PRIORITY as G4_PRIORITY_MAP
)

# =============================================================
# MODERNIZED G4 FAMILY DETECTION FUNCTIONS
# =============================================================

def convert_to_nbdfinder_format(candidates, sequence_name=""):
    """Convert modernized G4 candidates to NBDFinder format"""
    results = []
    for cand in candidates:
        # Map subclass names to NBDFinder conventions
        subclass_mapping = {
            "Multimeric/Tandem": "Multimeric G4",
            "Split/Bipartite": "Bipartite G4", 
            "Bulged": "Bulged G4",
            "Canonical G3+L1–7": "Canonical G4",
            "Long-Loop G3+L8–12": "Relaxed G4",
            "Extended G3+L1–12": "Relaxed G4",
            "Two-Tetrad G2L1–12": "Imperfect G4",
            "Generalized G2+L1–12": "Imperfect G4",
            "G-Triplex": "G-Triplex intermediate"
        }
        
        subtype = subclass_mapping.get(cand['subclass'], cand['subclass'])
        motif_seq = cand['seq']
        g4_score = cand['score']
        
        # Calculate structural factor and conservation
        structural_factor = calculate_structural_factor(motif_seq, subtype.lower().replace(" ", "_"), [])
        conservation_result = calculate_conservation_score(motif_seq, subtype)
        formation_data = get_g4_formation_category(g4_score)
        
        results.append({
            "Sequence Name": sequence_name,
            "Class": "G-Quadruplex Family", 
            "Subtype": subtype,
            "Start": cand['start'] + 1,  # 1-based indexing
            "End": cand['end'],
            "Length": len(motif_seq),
            "Sequence": wrap(motif_seq),
            "ScoreMethod": "G4Hunter_Exact",
            "Score": float(g4_score * len(motif_seq) * structural_factor),
            "G4Hunter_Score": float(g4_score),
            "Structural_Factor": round(structural_factor, 3),
            "Formation_Category": formation_data["category"],
            "Conservation_Score": float(conservation_result["enrichment_score"]),
            "Conservation_P_Value": float(conservation_result["p_value"]),
            "Priority": cand['priority'],
            "Arms/Repeat Unit/Copies": "",
            "Spacer": ""
        })
    return results

def find_all_g4_motifs(seq, use_non_overlapping=True, sequence_name=""):
    """
    Find all G4 motifs using modernized G4 family detection framework.
    
    Args:
        seq: DNA sequence
        use_non_overlapping: If True, apply priority-based non-overlapping selection
        sequence_name: Name for the sequence
        
    Returns:
        List of G4 motifs in NBDFinder format
    """
    # Get all G4 candidates using modernized patterns and exact G4Hunter scoring
    candidates = candidates_from(seq, strand='+')
    
    # Apply non-overlapping selection if requested
    if use_non_overlapping:
        candidates = non_overlapping(candidates)
    
    # Convert to NBDFinder format
    return convert_to_nbdfinder_format(candidates, sequence_name)

# =============================================================
# LEGACY COMPATIBILITY FUNCTIONS
# =============================================================

def g4hunter_score(seq):
    """
    Exact G4Hunter scoring algorithm.
    Replaces the previous unified_hunter_score with exact G4Hunter implementation.
    """
    return exact_g4hunter_score(seq)
# Legacy compatibility function - using modernized G4 detection
def find_g4_unified(seq, motif_subtype="Canonical G4", pattern=r"G{3,}\w{1,7}G{3,}\w{1,7}G{3,}\w{1,7}G{3,}", g4_threshold=1.0, min_length=12):
    """
    Unified G4 detection framework - updated to use modernized G4 family detection
    """
    # Use modernized G4 detection and filter by subtype
    all_results = find_all_g4_motifs(seq, use_non_overlapping=False)
    
    # Filter by requested subtype
    filtered_results = [r for r in all_results if r['Subtype'] == motif_subtype]
    
    # Apply additional filtering based on legacy parameters
    final_results = []
    for result in filtered_results:
        if (result['G4Hunter_Score'] >= g4_threshold and 
            result['Length'] >= min_length):
            final_results.append(result)
    
    return final_results

# Categorize G4 formation potential based on experimental thresholds.
def get_g4_formation_category(g4hunter_score):
    """Categorize G4 formation potential based on experimental thresholds."""
    if g4hunter_score >= 1.5:
        return {"category": "High Formation Potential", "threshold": "≥ 1.5", 
                "experimental_evidence": "Strong", "formation_probability": "85-95%",
                "stability": "High", "color": "#d32f2f"}
    elif g4hunter_score >= 1.0:
        return {"category": "Moderate Formation Potential", "threshold": "1.0 - 1.5",
                "experimental_evidence": "Moderate", "formation_probability": "60-85%", 
                "stability": "Moderate", "color": "#f57c00"}
    elif g4hunter_score < 1.0:
        return {"category": "Low Formation Potential", "threshold": "< 1.0",
                "experimental_evidence": "Weak/Variable", "formation_probability": "10-60%",
                "stability": "Low", "color": "#388e3c"}
    else:
        return {"category": "Unknown", "threshold": "N/A", "experimental_evidence": "Unknown",
                "formation_probability": "Unknown", "stability": "Unknown", "color": "#666666"}

# G4Hunter scoring algorithm - using exact implementation
def g4hunter_score(seq):
    """
    Exact G4Hunter scoring algorithm.
    Replaces the previous unified_hunter_score with exact G4Hunter implementation.
    """
    return exact_g4hunter_score(seq)

# G4 structural factor calculation using unified framework.
def g4_structural_factor(motif_seq, motif_type="canonical"):
    """
    Calculate structural factor using unified framework for consistency.
    - motif_seq: sequence of detected motif
    - motif_type: string describing the G4 subtype (default: canonical)
    """
    # Extract loop lengths between G-runs for G4 structures
    g_run_spans = [match.span() for match in re.finditer(r"G{3,}", motif_seq)]
    if len(g_run_spans) >= 2:
        loops = [g_run_spans[i+1][0] - g_run_spans[i][1] for i in range(len(g_run_spans)-1)]
    else:
        loops = []
    return calculate_structural_factor(motif_seq, motif_type, loops)

# Updated G4 detection functions using modernized G4 family framework

def find_gquadruplex(seq):
    """
    Find canonical G-quadruplex structures using modernized framework.
    """
    all_g4s = find_all_g4_motifs(seq, use_non_overlapping=False)
    return [g4 for g4 in all_g4s if g4['Subtype'] == 'Canonical G4']

def find_relaxed_gquadruplex(seq):
    """
    Find relaxed G-quadruplex structures using modernized framework.
    """
    all_g4s = find_all_g4_motifs(seq, use_non_overlapping=False)
    return [g4 for g4 in all_g4s if g4['Subtype'] == 'Relaxed G4']

def find_bulged_gquadruplex(seq):
    """
    Find bulged G-quadruplex structures using modernized framework.
    """
    all_g4s = find_all_g4_motifs(seq, use_non_overlapping=False)
    return [g4 for g4 in all_g4s if g4['Subtype'] == 'Bulged G4']

def find_imperfect_gquadruplex(seq):
    """
    Find imperfect G-quadruplex structures using modernized framework.
    """
    all_g4s = find_all_g4_motifs(seq, use_non_overlapping=False)
    return [g4 for g4 in all_g4s if g4['Subtype'] == 'Imperfect G4']

def find_multimeric_gquadruplex(seq):
    """
    Find multimeric G-quadruplex structures using modernized framework.
    """
    all_g4s = find_all_g4_motifs(seq, use_non_overlapping=False)
    return [g4 for g4 in all_g4s if g4['Subtype'] == 'Multimeric G4']

def find_bipartite_gquadruplex(seq):
    """
    Find bipartite G-quadruplex structures using modernized framework.
    """
    all_g4s = find_all_g4_motifs(seq, use_non_overlapping=False)
    return [g4 for g4 in all_g4s if g4['Subtype'] == 'Bipartite G4']

def find_gtriplex(seq):
    """
    Find G-triplex structures using modernized framework.
    """
    all_g4s = find_all_g4_motifs(seq, use_non_overlapping=False)
    return [g4 for g4 in all_g4s if g4['Subtype'] == 'G-Triplex intermediate']