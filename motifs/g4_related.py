"""
Category 8: G4-related Detection Module
======================================

This module implements detection algorithms for all G-quadruplex related structures
including canonical G4, variants, and G-triplex structures.

Scientific Basis:
- G-quadruplexes are four-stranded structures formed by guanine-rich sequences
- Various topologies: canonical, relaxed, bulged, bipartite, multimeric, imperfect
- G-triplex structures are three-stranded intermediates

References:
- Bedrat et al. (2016) NAR (G4Hunter algorithm)
- Hänsel-Hertsch et al. (2017) Nat Rev Mol Cell Biol

Author: Dr. Venkata Rajesh Yella
Updated: 2024
"""

import re

def mean(values):
    """Simple mean calculation."""
    return sum(values) / len(values) if values else 0
from .shared_utils import wrap, calculate_conservation_score, overlapping_finditer

def get_g4_formation_category(g4hunter_score):
    """
    Categorize G4 formation potential based on experimental data thresholds.
    
    Scientific Basis: Experimental validation shows different G4 formation
    propensities based on G4Hunter scores (Bedrat et al., 2016).
    
    Parameters:
    g4hunter_score (float): G4Hunter mean score
    
    Returns:
    dict: Category information with experimental formation data
    """
    if g4hunter_score >= 1.5:
        return {
            "category": "High Formation Potential",
            "threshold": "≥ 1.5",
            "experimental_evidence": "Strong",
            "formation_probability": "85-95%",
            "stability": "High",
            "color": "#d32f2f"  # Red
        }
    elif 1.0 <= g4hunter_score < 1.5:
        return {
            "category": "Moderate Formation Potential", 
            "threshold": "1.0 - 1.5",
            "experimental_evidence": "Moderate",
            "formation_probability": "60-85%",
            "stability": "Moderate",
            "color": "#f57c00"  # Orange
        }
    elif g4hunter_score < 1.0:
        return {
            "category": "Low Formation Potential",
            "threshold": "< 1.0",
            "experimental_evidence": "Weak/Variable",
            "formation_probability": "10-60%",
            "stability": "Low",
            "color": "#388e3c"  # Green
        }
    else:
        return {
            "category": "Unknown",
            "threshold": "N/A",
            "experimental_evidence": "Unknown",
            "formation_probability": "Unknown",
            "stability": "Unknown", 
            "color": "#666666"  # Gray
        }

def g4hunter_score(seq):
    """
    G4Hunter scoring function: assign scores based on G/C runs and average.
    
    Scientific Basis: G4Hunter computes G vs C bias to predict G-quadruplex
    formation potential. Positive scores indicate G-quadruplex propensity.
    Reference: Bedrat et al. (2016) Nucleic Acids Research
    
    Parameters:
    seq (str): DNA sequence to score
    
    Returns:
    float: G4Hunter mean score (run-based G vs C bias)
    """
    seq = seq.upper()
    scores = []
    i = 0
    while i < len(seq):
        # G-run: assign +score for each G in the run
        if seq[i] == 'G':
            run = 1
            while i + run < len(seq) and seq[i + run] == 'G':
                run += 1
            s = min(run, 4)
            scores += [s] * run
            i += run
        # C-run: assign -score for each C in the run
        elif seq[i] == 'C':
            run = 1
            while i + run < len(seq) and seq[i + run] == 'C':
                run += 1
            s = -min(run, 4)
            scores += [s] * run
            i += run
        # A/T: assign 0
        else:
            scores.append(0)
            i += 1
    return sum(scores) / len(scores) if scores else 0.0

def g4_structural_factor(motif_seq, motif_type="canonical"):
    """
    Calculate structural factor based on G-run architecture, loop lengths, and motif type.
    Accounts for loop lengths, number of G-runs, bulges, and architecture stability.
    """
    g_runs = re.findall(r"G{3,}", motif_seq)
    num_g_runs = len(g_runs)
    
    # Base structural factor
    structural_factor = 1.0
    
    # G-run count bonus (more G-runs = more stable)
    if num_g_runs >= 4:
        structural_factor += 0.1 * (num_g_runs - 4)  # Bonus for extra G-runs
    
    # Loop length analysis
    if num_g_runs >= 2:
        # Extract loops between G-runs
        pattern = r"G{3,}([ATGC]*?)G{3,}"
        loops = re.findall(pattern, motif_seq)
        if loops:
            avg_loop_len = np.mean([len(loop) for loop in loops])
            # Optimal loop lengths (1-7 nt) get bonus, longer loops get penalty
            if avg_loop_len <= 7:
                structural_factor += 0.1  # Compact loops bonus
            else:
                structural_factor -= 0.05 * (avg_loop_len - 7)  # Long loop penalty
    
    # Motif-type specific adjustments
    if motif_type == "bipartite":
        structural_factor += 0.2  # Bipartite architecture bonus
    elif motif_type == "multimeric":
        structural_factor += 0.15  # Multimeric architecture bonus
    elif motif_type == "bulged":
        structural_factor -= 0.1   # Bulges reduce stability
    elif motif_type == "relaxed":
        structural_factor -= 0.05  # Relaxed criteria slight penalty
    
    return max(0.1, structural_factor)  # Ensure positive factor

# G4 Detection Functions (simplified for space)
def find_gquadruplex(seq):
    """Find canonical G-quadruplex structures."""
    results = []
    pattern = r"G{3,}\w{1,7}G{3,}\w{1,7}G{3,}\w{1,7}G{3,}"
    for m in overlapping_finditer(pattern, seq):
        motif_seq = m.group(0)
        g4_score = g4hunter_score(motif_seq)
        if g4_score >= 1.0:  # G4Hunter threshold
            conservation_result = calculate_conservation_score(motif_seq, "Canonical G4")
            formation_data = get_g4_formation_category(g4_score)
            
            results.append({
                "Sequence Name": "",
                "Class": "Canonical G4",
                "Subtype": "G4_Canonical",
                "Start": m.start() + 1,
                "End": m.end(),
                "Length": len(motif_seq),
                "Sequence": wrap(motif_seq),
                "Score": float(g4_score * len(motif_seq)),
                "G4Hunter_Score": float(g4_score),
                "Formation_Category": formation_data["category"],
                "Conservation_Score": float(conservation_result["enrichment_score"]),
                "Conservation_P_Value": float(conservation_result["p_value"]),
                "Arms/Repeat Unit/Copies": "",
                "Spacer": ""
            })
    return results

def find_relaxed_gquadruplex(seq):
    """Find relaxed G-quadruplex structures with longer loops."""
    results = []
    pattern = r"G{3,}\w{1,12}G{3,}\w{1,12}G{3,}\w{1,12}G{3,}"
    for m in overlapping_finditer(pattern, seq):
        motif_seq = m.group(0)
        g4_score = g4hunter_score(motif_seq)
        if g4_score >= 0.8:  # Lower threshold for relaxed
            conservation_result = calculate_conservation_score(motif_seq, "Relaxed G4")
            
            results.append({
                "Sequence Name": "",
                "Class": "Relaxed G4",
                "Subtype": "G4_Relaxed",
                "Start": m.start() + 1,
                "End": m.end(),
                "Length": len(motif_seq),
                "Sequence": wrap(motif_seq),
                "Score": float(g4_score * len(motif_seq)),
                "G4Hunter_Score": float(g4_score),
                "Conservation_Score": float(conservation_result["enrichment_score"]),
                "Arms/Repeat Unit/Copies": "",
                "Spacer": ""
            })
    return results

def find_bulged_gquadruplex(seq):
    """Find bulged G-quadruplex structures."""
    results = []
    pattern = r"G{2,}\w{1,3}G{1,}\w{1,7}G{3,}\w{1,7}G{3,}\w{1,7}G{3,}"
    for m in overlapping_finditer(pattern, seq):
        motif_seq = m.group(0)
        g4_score = g4hunter_score(motif_seq)
        if g4_score >= 0.9:
            conservation_result = calculate_conservation_score(motif_seq, "Bulged G4")
            
            results.append({
                "Sequence Name": "",
                "Class": "Bulged G4",
                "Subtype": "G4_Bulged",
                "Start": m.start() + 1,
                "End": m.end(),
                "Length": len(motif_seq),
                "Sequence": wrap(motif_seq),
                "Score": float(g4_score * len(motif_seq)),
                "G4Hunter_Score": float(g4_score),
                "Conservation_Score": float(conservation_result["enrichment_score"]),
                "Arms/Repeat Unit/Copies": "",
                "Spacer": ""
            })
    return results

def find_imperfect_gquadruplex(seq):
    """Find imperfect G-quadruplex structures."""
    results = []
    pattern = r"G{2,}\w{1,12}G{2,}\w{1,12}G{2,}\w{1,12}G{2,}"
    for m in overlapping_finditer(pattern, seq):
        motif_seq = m.group(0)
        g4_score = g4hunter_score(motif_seq)
        if g4_score >= 0.7:
            conservation_result = calculate_conservation_score(motif_seq, "Imperfect G4")
            
            results.append({
                "Sequence Name": "",
                "Class": "Imperfect G4",
                "Subtype": "G4_Imperfect",
                "Start": m.start() + 1,
                "End": m.end(),
                "Length": len(motif_seq),
                "Sequence": wrap(motif_seq),
                "Score": float(g4_score * len(motif_seq)),
                "G4Hunter_Score": float(g4_score),
                "Conservation_Score": float(conservation_result["enrichment_score"]),
                "Arms/Repeat Unit/Copies": "",
                "Spacer": ""
            })
    return results

def find_multimeric_gquadruplex(seq):
    """Find multimeric G-quadruplex structures."""
    results = []
    pattern = r"(G{3,}\w{1,12}){4,}"
    for m in overlapping_finditer(pattern, seq):
        motif_seq = m.group(0)
        g4_score = g4hunter_score(motif_seq)
        if g4_score >= 1.2:
            conservation_result = calculate_conservation_score(motif_seq, "Multimeric G4")
            
            results.append({
                "Sequence Name": "",
                "Class": "Multimeric G4",
                "Subtype": "G4_Multimeric",
                "Start": m.start() + 1,
                "End": m.end(),
                "Length": len(motif_seq),
                "Sequence": wrap(motif_seq),
                "Score": float(g4_score * len(motif_seq)),
                "G4Hunter_Score": float(g4_score),
                "Conservation_Score": float(conservation_result["enrichment_score"]),
                "Arms/Repeat Unit/Copies": "",
                "Spacer": ""
            })
    return results

def find_bipartite_gquadruplex(seq):
    """Find bipartite G-quadruplex structures."""
    results = []
    pattern = r"G{3,}\w{1,7}G{3,}\w{20,50}G{3,}\w{1,7}G{3,}"
    for m in overlapping_finditer(pattern, seq):
        motif_seq = m.group(0)
        g4_score = g4hunter_score(motif_seq)
        if g4_score >= 0.8:
            conservation_result = calculate_conservation_score(motif_seq, "Bipartite G4")
            
            results.append({
                "Sequence Name": "",
                "Class": "Bipartite G4",
                "Subtype": "G4_Bipartite",
                "Start": m.start() + 1,
                "End": m.end(),
                "Length": len(motif_seq),
                "Sequence": wrap(motif_seq),
                "Score": float(g4_score * len(motif_seq)),
                "G4Hunter_Score": float(g4_score),
                "Conservation_Score": float(conservation_result["enrichment_score"]),
                "Arms/Repeat Unit/Copies": "",
                "Spacer": ""
            })
    return results

def find_gtriplex(seq):
    """Find G-triplex structures."""
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
                "Class": "G-Triplex",
                "Subtype": "G_Triplex",
                "Start": m.start() + 1,
                "End": m.end(),
                "Length": len(motif_seq),
                "Sequence": wrap(motif_seq),
                "Score": float(score),
                "G_Run_Count": len(g_runs),
                "G_Content": round(g_content, 3),
                "Conservation_Score": float(conservation_result["enrichment_score"]),
                "Arms/Repeat Unit/Copies": "",
                "Spacer": ""
            })
    return results