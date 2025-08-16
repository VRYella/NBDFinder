"""
Category 5: R-loop Detection Module
==================================

This module implements detection algorithms for R-loop structures
formed by RNA-DNA hybrids that displace the non-template DNA strand.

Scientific Basis:
- R-loops form when RNA-DNA hybrids displace the non-template DNA strand
- Require G-rich initiating zone (RIZ) and G-rich extending zone (REZ)
- Associated with transcription, DNA damage, and genetic instability

References:
- Aguilera & García-Muse (2012) Mol Cell
- Ginno et al. (2012) Mol Cell

Author: Dr. Venkata Rajesh Yella
Updated: 2024
"""

import re
from .shared_utils import wrap, calculate_conservation_score, gc_content

# RLFS models for R-loop forming sequences
RLFS_MODELS = {
    "m1": r"G{3,}[ATGC]{1,10}?G{3,}(?:[ATGC]{1,10}?G{3,}){1,}",
    "m2": r"G{4,}(?:[ATGC]{1,10}?G{4,}){1,}",
}

def advanced_rloop_score(riz_seq, rez_seq, w1=50.0, w2=10.0, alpha=0.25):
    """
    Advanced R-loop stability scoring: (GC_fraction × W1 + G_run_count × W2) × length^α
    
    Enhanced Scoring System:
    1. GC fraction component for R-loop stability
    2. G-run count component for initiation potential
    3. Length scaling factor to balance length dominance
    
    Parameters:
    riz_seq (str): R-loop initiating zone sequence
    rez_seq (str): R-loop extending zone sequence
    w1 (float): Weight for GC fraction component
    w2 (float): Weight for G-run count component
    alpha (float): Length scaling exponent
    
    Returns:
    float: Combined R-loop stability score
    """
    if not riz_seq and not rez_seq:
        return 0.0
    
    # Combine RIZ and REZ sequences
    combined_seq = riz_seq + rez_seq
    total_length = len(combined_seq)
    
    if total_length == 0:
        return 0.0
    
    # GC fraction component
    gc_count = combined_seq.count('G') + combined_seq.count('C')
    gc_fraction = gc_count / total_length
    
    # G-run count component (runs of 3+ consecutive G's)
    g_runs = len(re.findall(r"G{3,}", combined_seq))
    
    # Length scaling factor to temper length dominance while retaining scale
    length_factor = (total_length ** alpha)
    
    # Combined score formula
    score = (gc_fraction * w1 + g_runs * w2) * length_factor
    
    return score

def find_rez_advanced(seq, start_pos, max_search_len=2000, min_window=100, step=50, min_gc=40):
    """
    Advanced REZ detection: find the best GC-rich downstream region with optimized scoring.
    
    Enhanced Algorithm:
    1. Sliding window approach for optimal REZ detection
    2. GC content and length-based scoring
    3. Flexible window sizing for best REZ identification
    
    Parameters:
    seq (str): DNA sequence
    start_pos (int): Starting position for REZ search
    max_search_len (int): Maximum search length
    min_window (int): Minimum window size
    step (int): Step size for sliding window
    min_gc (float): Minimum GC content threshold
    
    Returns:
    dict: Best REZ region or None if not found
    """
    if start_pos >= len(seq):
        return None
    
    best_rez = None
    best_score = 0
    search_end = min(len(seq), start_pos + max_search_len)
    
    # Use sliding window approach to find optimal REZ
    for window_start in range(start_pos, search_end - min_window + 1, step):
        for window_size in range(min_window, min(max_search_len, search_end - window_start) + 1, step):
            window_end = window_start + window_size
            if window_end > len(seq):
                break
                
            window_seq = seq[window_start:window_end]
            gc_content_val = gc_content(window_seq)
            
            if gc_content_val >= min_gc:
                # Score this REZ candidate
                window_score = gc_content_val * len(window_seq) * 0.1  # Weight by GC content and length
                
                if window_score > best_score:
                    best_score = window_score
                    best_rez = {
                        'seq': window_seq,
                        'start': window_start - start_pos,  # Relative to start_pos
                        'end': window_end - start_pos,
                        'length': len(window_seq),
                        'gc_content': gc_content_val
                    }
    
    return best_rez

def find_rlfs(seq, models=("m1", "m2"), min_total_length=100):
    """
    Advanced RLFS detection using QmRLFS models with enhanced scoring.
    
    Enhanced Scoring System:
    1. Multiple model patterns for R-loop forming sequences
    2. Combined RIZ and REZ scoring for stability assessment
    3. Conservation analysis for evolutionary significance
    4. Minimum length thresholds for biological relevance
    
    Parameters:
    seq (str): DNA sequence
    models (tuple): RLFS model patterns to use
    min_total_length (int): Minimum total R-loop length
    
    Returns:
    list: List of R-loop motif dictionaries
    """
    if len(seq) < min_total_length:
        return []
    
    results = []
    
    for model_name in models:
        pattern = RLFS_MODELS.get(model_name)
        if not pattern:
            continue
            
        for m in re.finditer(pattern, seq, re.IGNORECASE):
            riz_seq = m.group(0)
            riz_end_pos = m.end()
            
            # Find downstream REZ
            rez = find_rez_advanced(seq, riz_end_pos)
            
            if rez:
                total_length = len(riz_seq) + rez['length']
                if total_length >= min_total_length:
                    # Calculate R-loop stability score
                    score = advanced_rloop_score(riz_seq, rez['seq'])
                    
                    # Calculate conservation score
                    full_rloop_seq = riz_seq + rez['seq']
                    conservation_result = calculate_conservation_score(full_rloop_seq, "R-Loop")
                    conservation_score = conservation_result["enrichment_score"]
                    
                    results.append({
                        "Sequence Name": "",
                        "Class": "R-Loop",
                        "Subtype": f"RLFS_{model_name}",
                        "Start": m.start() + 1,
                        "End": riz_end_pos + rez['end'],
                        "Length": total_length,
                        "Sequence": wrap(full_rloop_seq),
                        "ScoreMethod": "RLFS_RIZ_REZ_raw",
                        "Score": float(score),
                        "RIZ_Length": len(riz_seq),
                        "REZ_Length": rez['length'],
                        "REZ_GC_Content": rez['gc_content'],
                        "Conservation_Score": float(conservation_score),
                        "Conservation_P_Value": float(conservation_result["p_value"]),
                        "Conservation_Significance": conservation_result["significance"],
                        "Arms/Repeat Unit/Copies": f"RIZ={len(riz_seq)};REZ={rez['length']}",
                        "Spacer": ""
                    })
    
    return results