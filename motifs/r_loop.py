"""
Category 5: R-loop Detection Module
==================================

This module implements detection algorithms for R-loop structures using the
unified NBDFinder computational framework while preserving RLFS+REZ biological
accuracy for RNA-DNA hybrid formation prediction.

UNIFIED FRAMEWORK INTEGRATION:
-------------------------------
R-loop detection now integrates with the unified framework while maintaining:
    - Original RLFS pattern matching for RIZ (R-loop Initiating Zone) identification
    - Advanced REZ (R-loop Extending Zone) detection with GC optimization
    - Preserved stability scoring formula: (GC_fraction × W1 + G_runs × W2) × length^α
    - Bipartite R-loop structure analysis (RIZ + REZ)

SCIENTIFIC BASIS:
-----------------
    - R-loops form when RNA-DNA hybrids displace the non-template DNA strand
    - Require G-rich initiating zone (RIZ) and G-rich extending zone (REZ)
    - Associated with transcription, DNA damage, and genetic instability
    - Formation depends on GC skew and G-run density

BIOLOGICAL ACCURACY PRESERVATION:
---------------------------------
    - QmRLFS model patterns (m1, m2) for RIZ identification
    - Advanced REZ sliding window optimization with GC content thresholds
    - Original stability scoring with empirically derived weights (W1=50, W2=10, α=0.25)
    - Length scaling to prevent bias toward extremely long sequences

REFERENCES:
-----------
    - Aguilera & García-Muse (2012) Mol Cell
    - Ginno et al. (2012) Mol Cell  
    - Crossley et al. (2019) Mol Cell (RLFS methodology)

Author: Dr. Venkata Rajesh Yella
Updated: 2024 with unified framework integration
"""

import re
from .shared_utils import (wrap, calculate_conservation_score, gc_content)

# RLFS models for R-loop forming sequences
RLFS_MODELS = {
    "m1": r"G{3,}[ATGC]{1,10}?G{3,}(?:[ATGC]{1,10}?G{3,}){1,}",
    "m2": r"G{4,}(?:[ATGC]{1,10}?G{4,}){1,}",
}

# Advanced R-loop stability scoring using unified computational framework.
def advanced_rloop_score(riz_seq, rez_seq, w1=50.0, w2=10.0, alpha=0.25):
    """
    Combines Hunter-style scoring with R-loop specific factors.
    - Consistent G-bias scoring (Hunter-style)
    - Combines GC fraction and G-run count components
    - Length scaling to balance dominance
    """
    if not riz_seq and not rez_seq: return 0.0
    combined_seq = riz_seq + rez_seq; total_length = len(combined_seq)
    if total_length == 0: return 0.0
    from .shared_utils import unified_hunter_score, calculate_structural_factor
    hunter_score = unified_hunter_score(combined_seq, target_base='G', complementary_base='C')
    structural_factor = calculate_structural_factor(combined_seq, "R-loop")
    gc_count = combined_seq.count('G') + combined_seq.count('C')
    gc_fraction = gc_count / total_length
    g_runs = len(re.findall(r"G{3,}", combined_seq))
    length_factor = (total_length ** alpha)
    traditional_score = (gc_fraction * w1 + g_runs * w2) * length_factor
    unified_score = hunter_score * total_length * structural_factor
    combined_score = (traditional_score + unified_score) / 2
    return combined_score

# Advanced REZ detection: find the best GC-rich downstream region with optimized scoring.
def find_rez_advanced(seq, start_pos, max_search_len=2000, min_window=50, step=25, min_gc=50):  # Literature standard ≥50% GC (PMID: 22243696)
    """
    Sliding window approach for optimal REZ detection (GC-content and length-based scoring).
    """
    if start_pos >= len(seq): return None
    best_rez, best_score = None, 0
    search_end = min(len(seq), start_pos + max_search_len)
    for window_start in range(start_pos, search_end - min_window + 1, step):
        for window_size in range(min_window, min(max_search_len, search_end - window_start) + 1, step):
            window_end = window_start + window_size
            if window_end > len(seq): break
            window_seq = seq[window_start:window_end]
            gc_content_val = gc_content(window_seq)
            if gc_content_val >= min_gc:
                window_score = gc_content_val * len(window_seq) * 0.1
                if window_score > best_score:
                    best_score = window_score
                    best_rez = {
                        'seq': window_seq,
                        'start': window_start - start_pos,
                        'end': window_end - start_pos,
                        'length': len(window_seq),
                        'gc_content': gc_content_val
                    }
    return best_rez

# Advanced RLFS detection using unified framework with preserved biological accuracy.
def find_rlfs(seq, models=("m1", "m2"), min_total_length=100):  # Literature standard minimum 100bp (PMID: 30318411)
    """
    Maintains RLFS+REZ detection for biological accuracy, integrates unified scoring/reporting.
    """
    if len(seq) < min_total_length: return []
    results = []
    for model_name in models:
        pattern = RLFS_MODELS.get(model_name)
        if not pattern: continue
        for m in re.finditer(pattern, seq, re.IGNORECASE):
            riz_seq = m.group(0); riz_end_pos = m.end()
            rez = find_rez_advanced(seq, riz_end_pos)
            if rez:
                total_length = len(riz_seq) + rez['length']
                if total_length >= min_total_length:
                    score = advanced_rloop_score(riz_seq, rez['seq'])
                    full_rloop_seq = riz_seq + rez['seq']
                    conservation_result = calculate_conservation_score(full_rloop_seq, "R-loop")
                    conservation_score = conservation_result["enrichment_score"]
                    results.append({
                        "Sequence Name": "",
                        "Class": "R-loop",
                        "Subtype": "R-loop",
                        "Start": m.start() + 1,
                        "End": riz_end_pos + rez['end'],
                        "Length": total_length,
                        "Sequence": wrap(full_rloop_seq),
                        "ScoreMethod": "RLFS_UnifiedFramework_raw",
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