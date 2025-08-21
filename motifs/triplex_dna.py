"""
Category 4: Triplex DNA Detection Module
=======================================

This module implements detection algorithms for triplex DNA (H-DNA) and sticky DNA structures.

Scientific Basis:
- Triplex DNA forms at homopurine-homopyrimidine mirror repeats
- Sticky DNA structures form from GAA/TTC repeats in Friedreich's ataxia
- Triple helix formation through Hoogsteen hydrogen bonding

References:
- Frank-Kamenetskii & Mirkin (1995) Annu Rev Biochem
- Usdin (2008) Chromosome Res
- Campuzano et al. (1996) Science

Author: Dr. Venkata Rajesh Yella
Updated: 2024
"""

import re
from .shared_utils import wrap, calculate_conservation_score

def purine_fraction(seq):
    """Calculate fraction of purines (A, G) in sequence."""
    return (seq.count('A') + seq.count('G')) / max(1, len(seq))

def pyrimidine_fraction(seq):
    """Calculate fraction of pyrimidines (C, T) in sequence."""
    return (seq.count('C') + seq.count('T')) / max(1, len(seq))

def find_hdna(seq):
    """
    Find H-DNA (triplex DNA) structures at mirror repeats.
    
    Enhanced Scoring System:
    1. Mirror length scaling with homopurine/pyrimidine enrichment
    2. Spacer length penalty for structural flexibility
    3. Conservation analysis for evolutionary significance
    
    Parameters:
    seq (str): DNA sequence
    
    Returns:
    list: List of H-DNA motif dictionaries
    """
    results = []
    n = len(seq)
    min_score_threshold = 25.0  # Minimum score to avoid excessive low-quality matches
    used_positions = set()  # Track used positions to prevent excessive overlap
    
    for rep_len in range(10, min(101, n//2)):
        for spacer in range(0, 9):
            pattern = re.compile(rf"(([ATGC]{{{rep_len}}})[ATGC]{{{spacer}}}\2)", re.IGNORECASE)
            for m in pattern.finditer(seq):
                repeat = m.group(2)
                mirror_start = m.start()
                mirror_end = m.end()
                
                # Skip if this region significantly overlaps with an already found motif
                current_positions = set(range(mirror_start, mirror_end))
                overlap_ratio = len(current_positions.intersection(used_positions)) / len(current_positions)
                if overlap_ratio > 0.5:  # Skip if >50% overlap
                    continue
                
                full_seq = seq[mirror_start:mirror_end]
                pur_frac = purine_fraction(full_seq)
                pyr_frac = pyrimidine_fraction(full_seq)
                is_triplex = (pur_frac >= 0.9 or pyr_frac >= 0.9)
                # Enhanced scoring: mirror_length × homopurine/pyrimidine_enrichment - spacer_penalty
                homogeneity = max(pur_frac, pyr_frac)
                score = len(full_seq) * (1.0 + 1.5*homogeneity) - spacer * 1.0
                
                # Only keep high-scoring matches to avoid excessive low-quality results
                if score >= min_score_threshold:
                    # Calculate conservation score
                    conservation_result = calculate_conservation_score(full_seq, "Triplex")
                    conservation_score = conservation_result["enrichment_score"]
                    
                    results.append({
                        "Sequence Name": "",
                        "Class": "Triplex",
                        "Subtype": "Triplex" if is_triplex else "Triplex",
                        "Start": mirror_start + 1,
                        "End": mirror_end,
                        "Length": len(full_seq),
                        "Spacer": spacer,
                        "Sequence": wrap(full_seq),
                        "ScoreMethod": "Triplex_Homogeneity_raw",
                        "Score": float(score),
                        "PurineFrac": round(pur_frac, 3),
                        "PyrimidineFrac": round(pyr_frac, 3),
                        "Homogeneity": round(homogeneity, 3),
                        "Conservation_Score": float(conservation_score),
                        "Conservation_P_Value": float(conservation_result["p_value"]),
                        "Conservation_Significance": conservation_result["significance"],
                        "Arms/Repeat Unit/Copies": f"Arms={rep_len}",
                        "Spacer": str(spacer)
                    })
                    # Mark these positions as used
                    used_positions.update(current_positions)
    return results

def find_sticky_dna(seq):
    """
    Enhanced Sticky DNA Detection with Improved Sensitivity
    
    Scientific Basis: GAA/TTC triplet repeats form sticky DNA structures through 
    inter-strand purine-purine interactions. Pathogenic thresholds: >59 repeats 
    for Friedreich's ataxia (Campuzano et al., Science 1996).
    
    Enhanced Scoring System:
    1. Copy number × unit length × AT-richness bonus
    2. Pathogenic threshold marking (≥59 repeats)
    3. Conservation analysis for repeat stability
    
    Parameters:
    seq (str): DNA sequence
    
    Returns:
    list: List of sticky DNA motif dictionaries
    """
    motifs = []
    seq = seq.replace('\n','').replace(' ','').upper()
    
    # Enhanced pattern with lower threshold to capture more expansions
    pattern = r"(?:GAA){6,}|(?:TTC){6,}"  # Lowered threshold for better sensitivity
    
    for m in re.finditer(pattern, seq):
        repeat_len = len(m.group())
        repeat_count = repeat_len // 3
        
        # Enhanced scoring: copies × unit_length × AT_richness_bonus
        # Very high thresholds (≥59 repeats) mark pathogenic ranges
        at_frac = (m.group().count('A') + m.group().count('T')) / repeat_len
        score = repeat_count * 3 * (1.0 + 0.5*at_frac)
        
        # Mark pathogenic ranges
        pathogenic = repeat_count >= 59
        subtype = "GAA_TTC_Pathogenic" if pathogenic else "GAA_TTC_Repeat"
        
        # Calculate conservation score
        conservation_result = calculate_conservation_score(m.group(), "Sticky DNA")
        conservation_score = conservation_result["enrichment_score"]
        
        motifs.append({
            "Sequence Name": "",
            "Class": "Triplex",
            "Subtype": "sticky DNA",
            "Start": m.start() + 1,
            "End": m.end(),
            "Length": repeat_len,
            "RepeatCount": repeat_count,
            "Sequence": wrap(m.group()),
            "ScoreMethod": "GAA_TTC_Pathogenic_raw",
            "Score": float(score),
            "AT_Fraction": round(at_frac, 3),
            "Pathogenic": pathogenic,
            "Conservation_Score": float(conservation_score),
            "Conservation_P_Value": float(conservation_result["p_value"]),
            "Conservation_Significance": conservation_result["significance"],
            "Arms/Repeat Unit/Copies": f"Unit={'GAA' if 'GAA' in m.group() else 'TTC'};Copies={repeat_count}",
            "Spacer": ""
        })
    return motifs