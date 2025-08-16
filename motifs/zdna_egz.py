"""
Category 6: Z-DNA and eGZ Detection Module
=========================================

This module implements detection algorithms for Z-DNA and extruded-G (eGZ) structures.

Scientific Basis:
- Z-DNA is a left-handed double helix favored by alternating purine-pyrimidine sequences
- eGZ structures form from CGG repeat expansions with extruded guanines
- Associated with genetic instability and disease mechanisms

References:
- Rich & Zhang (2003) Nat Rev Genet
- Wang et al. (1979) Nature
- Usdin & Woodford (2007) DNA Repair

Author: Dr. Venkata Rajesh Yella
Updated: 2024
"""

import re
import numpy as np
from .shared_utils import wrap, calculate_conservation_score

def zdna_dinucleotide_weights(seq, 
                             gc_weight=7.0, 
                             gt_ca_weight=1.25, 
                             at_base_weight=0.5,
                             at_consecutive_penalty=(-5.0, -100.0),
                             mismatch_penalty_type="linear",
                             mismatch_penalty_base=3,
                             mismatch_penalty_delta=3):
    """
    Advanced Z-DNA dinucleotide weight array with optimized penalties.
    
    Scoring System:
    - GC/CG dinucleotides: High weight (canonical Z-DNA)
    - GT/TG and AC/CA: Moderate weight (Z-DNA-favorable)
    - AT/TA: Base weight with consecutive penalties
    - Unknown bases: Mismatch penalties
    
    Parameters:
    seq (str): DNA sequence
    gc_weight (float): Weight for GC/CG dinucleotides
    gt_ca_weight (float): Weight for GT/TG/AC/CA dinucleotides
    at_base_weight (float): Base weight for AT/TA dinucleotides
    at_consecutive_penalty (tuple): Penalties for consecutive AT dinucleotides
    
    Returns:
    np.array: Dinucleotide weights for Z-DNA propensity
    """
    if len(seq) < 2:
        return np.array([])
    
    weights = np.empty(len(seq) - 1, dtype=float)
    consecutive_at_count = 0
    consecutive_mismatch_count = 0
    
    for i in range(len(seq) - 1):
        dinuc = seq[i:i+2].upper()
        
        if dinuc in ("GC", "CG"):
            # High weight for canonical Z-DNA dinucleotides
            weights[i] = gc_weight
            consecutive_at_count = 0
            consecutive_mismatch_count = 0
            
        elif dinuc in ("GT", "TG", "AC", "CA"):
            # Moderate weight for Z-DNA-favorable dinucleotides
            weights[i] = gt_ca_weight
            consecutive_at_count = 0
            consecutive_mismatch_count = 0
            
        elif dinuc in ("AT", "TA"):
            # AT base weight with consecutive penalty
            base_weight = at_base_weight
            if consecutive_at_count >= 4:  # Start penalties after 4 consecutive AT
                if consecutive_at_count < 6:
                    base_weight += at_consecutive_penalty[0]
                else:
                    base_weight += at_consecutive_penalty[1]
            weights[i] = base_weight
            consecutive_at_count += 1
            consecutive_mismatch_count = 0
            
        else:
            # Mismatch penalty (non-standard bases or N's)
            consecutive_mismatch_count += 1
            consecutive_at_count = 0
            
            if mismatch_penalty_type == "exponential":
                penalty = -(mismatch_penalty_base ** min(consecutive_mismatch_count, 10))
            else:  # linear
                penalty = -mismatch_penalty_base - mismatch_penalty_delta * (consecutive_mismatch_count - 1)
            
            weights[i] = penalty
    
    return weights

def kadane_maximum_subarray(weights, min_length=12):
    """
    Kadane's algorithm for maximum subarray sum, optimized for Z-DNA detection.
    
    Enhanced Algorithm:
    1. Finds all potential Z-DNA segments using modified Kadane's algorithm
    2. Filters segments by minimum length and score criteria
    3. Removes overlapping segments, keeping highest scoring ones
    
    Parameters:
    weights (np.array): Dinucleotide weights array
    min_length (int): Minimum segment length
    
    Returns:
    list: List of (start, end, score) tuples for significant segments
    """
    if len(weights) < min_length:
        return []
    
    segments = []
    n = len(weights)
    
    # Find all potential Z-DNA segments using modified Kadane's
    for start in range(n - min_length + 1):
        max_ending_here = 0
        max_so_far = float('-inf')
        current_start = start
        best_start = start
        best_end = start
        
        for end in range(start, n):
            max_ending_here += weights[end]
            
            if max_ending_here > max_so_far:
                max_so_far = max_ending_here
                best_start = current_start
                best_end = end
            
            if max_ending_here < 0:
                max_ending_here = 0
                current_start = end + 1
        
        # Only keep segments that meet minimum criteria
        segment_length = best_end - best_start + 1
        if max_so_far > 0 and segment_length >= min_length:
            segments.append((best_start, best_end, max_so_far))
    
    # Remove overlapping segments, keeping the highest scoring ones
    segments.sort(key=lambda x: x[2], reverse=True)  # Sort by score descending
    non_overlapping = []
    used_positions = set()
    
    for start, end, score in segments:
        segment_positions = set(range(start, end + 1))
        if not segment_positions.intersection(used_positions):
            non_overlapping.append((start, end, score))
            used_positions.update(segment_positions)
    
    return non_overlapping

def find_zdna(seq, threshold=50, min_length=12, **kwargs):
    """
    Find Z-DNA using advanced dinucleotide weights and Kadane's maximum subarray.
    
    Enhanced Scoring System:
    1. Dinucleotide weight calculation with sequence-specific penalties
    2. Kadane's maximum subarray algorithm for optimal segment detection
    3. Conservation analysis using composition-preserving shuffles
    4. GC content and dinucleotide composition metrics
    
    Parameters:
    seq (str): DNA sequence
    threshold (float): Minimum score threshold for detection
    min_length (int): Minimum motif length
    **kwargs: Additional parameters for dinucleotide weights
    
    Returns:
    list: List of Z-DNA motif dictionaries
    """
    seq = seq.upper()
    if len(seq) < min_length:
        return []
    
    # Calculate dinucleotide weights
    weights = zdna_dinucleotide_weights(seq, **kwargs)
    
    # Find maximum subarrays using Kadane's algorithm
    segments = kadane_maximum_subarray(weights, min_length)
    
    motifs = []
    for start_idx, end_idx, raw_score in segments:
        if raw_score >= threshold:
            # Adjust indices for dinucleotide to sequence mapping
            seq_start = start_idx
            seq_end = end_idx + 1  # +1 because dinucleotide represents position i to i+1
            
            motif_seq = seq[seq_start:seq_end + 1]
            
            # Calculate additional scoring metrics
            gc_content = (motif_seq.count('G') + motif_seq.count('C')) / len(motif_seq)
            gc_cg_dinucs = sum(1 for i in range(len(motif_seq)-1) 
                              if motif_seq[i:i+2] in ['GC', 'CG'])
            
            # Calculate conservation score
            conservation_result = calculate_conservation_score(motif_seq, "Z-DNA")
            conservation_score = conservation_result["enrichment_score"]
            
            motifs.append({
                "Sequence Name": "",
                "Class": "Z-DNA",
                "Subtype": "Z-Seeker_Kadane",
                "Start": seq_start + 1,  # 1-based indexing
                "End": seq_end + 1,
                "Length": len(motif_seq),
                "Sequence": wrap(motif_seq),
                "ScoreMethod": "Kadane_MaxSubarray_raw",
                "Score": float(raw_score),
                "GC_Content": round(gc_content, 3),
                "GC_CG_Dinucleotides": gc_cg_dinucs,
                "Conservation_Score": float(conservation_score),
                "Conservation_P_Value": float(conservation_result["p_value"]),
                "Conservation_Significance": conservation_result["significance"],
                "Arms/Repeat Unit/Copies": "",
                "Spacer": ""
            })
    
    return motifs

def find_egz_motif(seq):
    """
    Enhanced eGZ (Extruded-G) Motif Detection
    
    Scientific Basis: CGG repeat expansions form extruded-G structures with left-handed
    conformations (Usdin & Woodford, DNA Repair 2007). Associated with fragile X syndrome
    and other repeat expansion disorders.
    
    Enhanced Scoring System:
    1. Copy number-based scoring with G-content bias
    2. Conservation analysis for repeat stability
    3. Clinical relevance thresholds (3+ repeats for sensitivity)
    
    Parameters:
    seq (str): DNA sequence
    
    Returns:
    list: List of eGZ motif dictionaries
    """
    pattern = re.compile(r'(CGG){3,}', re.IGNORECASE)  # Lowered from 4 to 3 for sensitivity
    results = []
    for m in pattern.finditer(seq):
        motif_seq = m.group(0)
        n_repeats = len(motif_seq) // 3
        
        # Enhanced scoring: copies × unit_length × G-bias factor
        # CGG repeats tie to extruded-G/left-handed conformations with strong G-content emphasis
        g_frac = motif_seq.count('G') / len(motif_seq)
        score = n_repeats * 3 * (1.0 + 2.0*g_frac)
        
        # Calculate conservation score
        conservation_result = calculate_conservation_score(motif_seq, "eGZ")
        conservation_score = conservation_result["enrichment_score"]
        
        results.append({
            "Sequence Name": "",
            "Family": "Double-stranded",
            "Class": "eGZ (Extruded-G)",
            "Subtype": "CGG_Expansion",
            "Start": m.start() + 1,
            "End": m.end(),
            "Length": len(motif_seq),
            "Sequence": wrap(motif_seq),
            "ScoreMethod": "eGZ_CGG_Expansion_raw",
            "Score": float(score),
            "CGG_Repeats": n_repeats,
            "G_Fraction": round(g_frac, 3),
            "Conservation_Score": float(conservation_score),
            "Conservation_P_Value": float(conservation_result["p_value"]),
            "Conservation_Significance": conservation_result["significance"],
            "Arms/Repeat Unit/Copies": f"Unit=CGG;Copies={n_repeats}",
            "Spacer": ""
        })
    return results