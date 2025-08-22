"""
Category 8: G4-related Detection Module
=======================================

This module implements detection algorithms for all G-quadruplex related structures
using the unified NBDFinder computational framework while preserving G4Hunter
biological accuracy.

UNIFIED FRAMEWORK INTEGRATION:
------------------------------
All G4 detection functions now use the unified framework through find_g4_unified(),
which maintains the exact G4Hunter scoring methodology while providing:
    - Standardized computational pipeline
    - Consistent scoring array generation  
    - Enhanced region detection and filtering
    - Preserved formation potential classification
    - Maintained structural factor calculations

SCIENTIFIC BASIS:
-----------------
    - G-quadruplexes are four-stranded structures formed by guanine-rich sequences
    - Various topologies: canonical, relaxed, bulged, bipartite, multimeric, imperfect
    - G-triplex structures are three-stranded intermediates
    - Formation potential based on experimental G4Hunter thresholds (Bedrat et al., 2016)

BIOLOGICAL ACCURACY PRESERVATION:
---------------------------------
    - Preserves all original G4Hunter scoring: G/C run identification (capped at 4), 
      average score calculation, formation category thresholds (≥1.5 High, ≥1.0 Moderate, <1.0 Low),
      structural factor calculations for stability assessment

REFERENCES:
-----------
    - Bedrat et al. (2016) NAR (G4Hunter algorithm)
    - Hänsel-Hertsch et al. (2017) Nat Rev Mol Cell Biol
    - Chambers et al. (2015) Nat Biotechnol (G4 experimental validation)

Author: Dr. Venkata Rajesh Yella
Updated: 2024 with unified framework integration
"""

import re
from math import log2, sqrt

def mean(values): 
    """Mean utility for lists, guarded for empty lists."""
    return sum(values) / len(values) if values else 0

# Import shared utilities from the unified NBDFinder framework
from .shared_utils import (
    wrap, calculate_conservation_score, overlapping_finditer,
    unified_hunter_score, calculate_structural_factor
)

def find_g4_unified(seq, motif_subtype="Canonical G4", pattern=r"G{3,}\w{1,7}G{3,}\w{1,7}G{3,}\w{1,7}G{3,}", g4_threshold=1.0, min_length=12):
    """
    Unified G4 detection framework
    """
    results = []
    for m in overlapping_finditer(pattern, seq):
        motif_seq = m.group(0)
        if len(motif_seq) < min_length:
            continue
        g4_score = g4hunter_score(motif_seq)
        if g4_score >= g4_threshold:
            structural_factor = g4_structural_factor(motif_seq, motif_subtype.lower().replace(" ", "_"))
            conservation_result = calculate_conservation_score(motif_seq, motif_subtype)
            formation_data = get_g4_formation_category(g4_score)
            results.append({
                "Sequence Name": "", "Class": "G-Quadruplex Family", "Subtype": motif_subtype,
                "Start": m.start() + 1, "End": m.end(), "Length": len(motif_seq),
                "Sequence": wrap(motif_seq), "ScoreMethod": "G4_UnifiedHunter_raw",
                "Score": float(g4_score * len(motif_seq) * structural_factor),
                "G4Hunter_Score": float(g4_score), "Structural_Factor": round(structural_factor, 3),
                "Formation_Category": formation_data["category"],
                "Conservation_Score": float(conservation_result["enrichment_score"]),
                "Conservation_P_Value": float(conservation_result["p_value"]),
                "Arms/Repeat Unit/Copies": "", "Spacer": ""
            })
    return results

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

# G4Hunter scoring algorithm using unified computational framework.
def g4hunter_score(seq):
    """
    G4Hunter scoring using unified computational framework.
    Maintains original G4Hunter algorithm accuracy while using consistent logic structure.
    """
    return unified_hunter_score(seq, target_base='G', complementary_base='C')

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

# Canonical quadruplex detection
def find_gquadruplex(seq):
    """
    Find canonical G-quadruplex structures using unified framework.
    - Uses G4Hunter scoring and structural assessment.
    """
    results = []
    pattern = r"G{3,}\w{1,7}G{3,}\w{1,7}G{3,}\w{1,7}G{3,}"
    for m in overlapping_finditer(pattern, seq):
        motif_seq = m.group(0)
        g4_score = g4hunter_score(motif_seq)
        if g4_score >= 1.0:  # G4Hunter canonical threshold
            structural_factor = g4_structural_factor(motif_seq, "canonical")
            conservation_result = calculate_conservation_score(motif_seq, "Canonical G4")
            formation_data = get_g4_formation_category(g4_score)
            results.append({
                "Sequence Name": "", "Class": "G-Quadruplex Family", "Subtype": "Canonical G4",
                "Start": m.start() + 1, "End": m.end(), "Length": len(motif_seq),
                "Sequence": wrap(motif_seq), "ScoreMethod": "G4_UnifiedHunter_raw",
                "Score": float(g4_score * len(motif_seq) * structural_factor),
                "G4Hunter_Score": float(g4_score), "Structural_Factor": round(structural_factor, 3),
                "Formation_Category": formation_data["category"],
                "Conservation_Score": float(conservation_result["enrichment_score"]),
                "Conservation_P_Value": float(conservation_result["p_value"]),
                "Arms/Repeat Unit/Copies": "", "Spacer": ""
            })
    return results

# Relaxed quadruplex detection (longer loops)
def find_relaxed_gquadruplex(seq):
    """
    Find relaxed G-quadruplex structures with longer loops using unified framework.
    """
    return find_g4_unified(
        seq,
        motif_subtype="Relaxed G4", 
        pattern=r"G{3,}\w{1,12}G{3,}\w{1,12}G{3,}\w{1,12}G{3,}",
        g4_threshold=0.8,
        min_length=12
    )

# Bulged quadruplex detection
def find_bulged_gquadruplex(seq):
    """
    Find bulged G-quadruplex structures using unified framework.
    """
    return find_g4_unified(
        seq,
        motif_subtype="Bulged G4",
        pattern=r"G{2,}\w{1,3}G{1,}\w{1,7}G{3,}\w{1,7}G{3,}\w{1,7}G{3,}",
        g4_threshold=0.9,
        min_length=12
    )

# Imperfect quadruplex detection (G-to-A substitutions)
def find_imperfect_gquadruplex(seq):
    """
    Find imperfect G-quadruplex structures with G-to-A substitutions using unified framework.
    """
    return find_g4_unified(
        seq,
        motif_subtype="Imperfect G4",
        pattern=r"[GA]{3,}\w{1,7}[GA]{3,}\w{1,7}[GA]{3,}\w{1,7}[GA]{3,}",
        g4_threshold=0.7,
        min_length=12
    )

# Multimeric quadruplex detection
def find_multimeric_gquadruplex(seq):
    """
    Find multimeric G-quadruplex structures using unified framework.
    """
    return find_g4_unified(
        seq,
        motif_subtype="Multimeric G4",
        pattern=r"(G{3,}\w{1,12}){4,}",
        g4_threshold=1.2,
        min_length=15
    )

# Bipartite quadruplex detection
def find_bipartite_gquadruplex(seq):
    """
    Find bipartite G-quadruplex structures using unified framework.
    """
    return find_g4_unified(
        seq,
        motif_subtype="Bipartite G4", 
        pattern=r"G{3,}\w{1,7}G{3,}\w{20,50}G{3,}\w{1,7}G{3,}",
        g4_threshold=0.8,
        min_length=15
    )

# G-triplex detection: three G runs, experimental intermediate
def find_gtriplex(seq):
    """
    Find G-triplex structures using specialized scoring (not covered by G4Hunter).
    """
    results = []
    pattern = r"G{3,}\w{1,15}G{3,}\w{1,15}G{3,}"
    for m in overlapping_finditer(pattern, seq):
        motif_seq = m.group(0)
        g_runs = re.findall(r"G{3,}", motif_seq)
        if len(g_runs) >= 3:
            g_content = motif_seq.count('G') / len(motif_seq)
            score = len(motif_seq) * g_content * len(g_runs)
            conservation_result = calculate_conservation_score(motif_seq, "G-Triplex")
            results.append({
                "Sequence Name": "",
                "Class": "G-Quadruplex Family",
                "Subtype": "G-Triplex intermediate",
                "Start": m.start() + 1,
                "End": m.end(),
                "Length": len(motif_seq),
                "Sequence": wrap(motif_seq),
                "Score": float(score),
                "G_Run_Count": len(g_runs),
                "G_Content": round(g_content, 3),
                "Conservation_Score": float(conservation_result["enrichment_score"]),
                "ScoreMethod": "G_Content_Run_Count",
                "Arms/Repeat Unit/Copies": "",
                "Spacer": ""
            })
    return results