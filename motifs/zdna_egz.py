"""
Category 6: Z-DNA and eGZ Detection Module
==========================================

This module implements detection algorithms for Z-DNA and extruded-G (eGZ) structures
using the unified NBDFinder computational framework while preserving the advanced
Kadane's maximum subarray algorithm for optimal region identification.

UNIFIED FRAMEWORK INTEGRATION:
------------------------------
- Original dinucleotide weight calculation methodology
- Preserved Kadane's maximum subarray algorithm for region optimization
- Maintained threshold-based filtering and confidence scoring
- Enhanced scoring with GC content and alternating pattern analysis

SCIENTIFIC BASIS:
-----------------
- Z-DNA: favored by alternating purine-pyrimidine sequences (GC/CG/GT/TG/AC/CA)
- eGZ: formed from CGG repeat expansions with extruded guanines

BIOLOGICAL ACCURACY PRESERVATION:
---------------------------------
- Dinucleotide weight assignments: GC/CG (+7), GT/TG/AC/CA (+1.25), AT/TA (0)
- Advanced Kadane's with overlap resolution
- Confidence scoring, length and composition-based filtering

REFERENCES:
-----------
- Rich & Zhang (2003) Nat Rev Genet
- Wang et al. (1979) Nature (Z-DNA structure discovery)
- Usdin & Woodford (2007) DNA Repair (eGZ and repeat expansions)

Author: Dr. Venkata Rajesh Yella
Updated: 2024 with unified framework integration
"""

import re
from math import log2, sqrt
from .shared_utils import (wrap, calculate_conservation_score,
                          unified_hunter_score, calculate_structural_factor)

# Dinucleotide weights for Z-DNA
def zdna_dinucleotide_weights(seq, gc_weight=7.0, gt_ca_weight=1.25, at_base_weight=0.5,
                             at_consecutive_penalty=(-5.0, -100.0), mismatch_penalty_base=3):
    """Advanced Z-DNA dinucleotide weight array with optimized penalties."""
    if len(seq) < 2: return []
    weights, consecutive_at_count = [0.0] * (len(seq) - 1), 0
    for i in range(len(seq) - 1):
        dinuc = seq[i:i+2].upper()
        if dinuc in ("GC", "CG"):
            weights[i] = gc_weight; consecutive_at_count = 0
        elif dinuc in ("GT", "TG", "AC", "CA"):
            weights[i] = gt_ca_weight; consecutive_at_count = 0
        elif dinuc in ("AT", "TA"):
            base_weight = at_base_weight
            if consecutive_at_count >= 4:
                base_weight += at_consecutive_penalty[0] if consecutive_at_count < 6 else at_consecutive_penalty[1]
            weights[i] = base_weight; consecutive_at_count += 1
        else:
            weights[i] = -mismatch_penalty_base; consecutive_at_count = 0
    return weights

# Optimized Kadane's algorithm for maximum subarray
def kadane_maximum_subarray(weights, min_length=12):
    """Kadane's algorithm for maximum subarray sum, optimized for Z-DNA detection."""
    if len(weights) < min_length: return []
    segments, n = [], len(weights)
    for start in range(n - min_length + 1):
        max_ending_here = 0; max_so_far = float('-inf')
        current_start = best_start = best_end = start
        for end in range(start, n):
            max_ending_here += weights[end]
            if max_ending_here > max_so_far:
                max_so_far = max_ending_here; best_start = current_start; best_end = end
            if max_ending_here < 0:
                max_ending_here = 0; current_start = end + 1
        region_length = best_end - best_start + 1
        if max_so_far > 0 and region_length >= min_length:
            segments.append((best_start, best_end, max_so_far))
    segments.sort(key=lambda x: x[2], reverse=True)
    non_overlapping, used_positions = [], set()
    for start, end, score in segments:
        segment_positions = set(range(start, end + 1))
        if not segment_positions.intersection(used_positions):
            non_overlapping.append((start, end, score))
            used_positions.update(segment_positions)
    return non_overlapping

# Z-DNA detection using unified framework and Kadane's
def find_zdna(seq, threshold=50, min_length=12, **kwargs):
    """
    Find Z-DNA using unified computational framework with Kadane's maximum subarray.
    Maintains Kadane's algorithm while incorporating unified scoring for consistency.
    """
    seq = seq.upper()
    if len(seq) < min_length: return []
    weights = zdna_dinucleotide_weights(seq, **kwargs)
    segments = kadane_maximum_subarray(weights, min_length)
    motifs = []
    for start_idx, end_idx, raw_score in segments:
        if raw_score >= threshold:
            seq_start, seq_end = start_idx, end_idx + 1
            motif_seq = seq[seq_start:seq_end + 1]
            unified_score = unified_hunter_score(motif_seq, target_base='G', complementary_base='A')
            structural_factor = calculate_structural_factor(motif_seq, "Z-DNA")
            combined_score = (raw_score + abs(unified_score) * len(motif_seq)) / 2 * structural_factor
            gc_content = (motif_seq.count('G') + motif_seq.count('C')) / len(motif_seq)
            gc_cg_dinucs = sum(1 for i in range(len(motif_seq)-1) if motif_seq[i:i+2] in ['GC', 'CG'])
            conservation_result = calculate_conservation_score(motif_seq, "Z-DNA")
            motifs.append({
                "Sequence Name": "", "Class": "Z-DNA", "Subtype": "Z-Seeker_Kadane",
                "Start": seq_start + 1, "End": seq_end + 1, "Length": len(motif_seq),
                "Sequence": wrap(motif_seq), "ScoreMethod": "Kadane_UnifiedFramework_raw", 
                "Score": float(combined_score), "Kadane_Score": float(raw_score),
                "Unified_Score": float(unified_score), "Structural_Factor": round(structural_factor, 3),
                "GC_Content": round(gc_content, 3), "GC_CG_Dinucleotides": gc_cg_dinucs,
                "Conservation_Score": float(conservation_result["enrichment_score"]),
                "Conservation_P_Value": float(conservation_result["p_value"]),
                "Conservation_Significance": conservation_result["significance"],
                "Arms/Repeat Unit/Copies": "", "Spacer": ""
            })
    return motifs

# eGZ motif detection for CGG repeat expansions
def find_egz_motif(seq):
    """Enhanced eGZ (Extruded-G) motif detection for CGG repeat expansions."""
    pattern = re.compile(r'(CGG){3,}', re.IGNORECASE)  # 3+ repeats for sensitivity
    results = []
    for m in pattern.finditer(seq):
        motif_seq = m.group(0)
        n_repeats = len(motif_seq) // 3
        g_frac = motif_seq.count('G') / len(motif_seq)
        score = n_repeats * 3 * (1.0 + 2.0 * g_frac)
        conservation_result = calculate_conservation_score(motif_seq, "eGZ")
        results.append({
            "Sequence Name": "", "Family": "Double-stranded", "Class": "Z-DNA", "Subtype": "eGZ (Extruded-G) DNA",
            "Start": m.start() + 1, "End": m.end(), "Length": len(motif_seq), "Sequence": wrap(motif_seq),
            "ScoreMethod": "eGZ_CGG_Expansion_raw", "Score": float(score), "CGG_Repeats": n_repeats,
            "G_Fraction": round(g_frac, 3), "Conservation_Score": float(conservation_result["enrichment_score"]),
            "Conservation_P_Value": float(conservation_result["p_value"]),
            "Conservation_Significance": conservation_result["significance"],
            "Arms/Repeat Unit/Copies": f"Unit=CGG;Copies={n_repeats}", "Spacer": ""
        })
    return results