"""
Category 7: i-motif and AC-motif Detection Module
================================================

This module implements detection algorithms for i-motif and AC-motif structures.

Scientific Basis:
- i-motifs are C-rich structures complementary to G-quadruplexes
- AC-motifs show alternating A/C patterns forming non-canonical structures
- Both are pH-dependent and involved in gene regulation

References:
- Zeraati et al. (2018) Nat Chem
- Kocsis et al. (2021) Int J Mol Sci
- Varizhuk et al. (2019) Biochimie

Author: Dr. Venkata Rajesh Yella
Updated: 2024
"""

import re
from .shared_utils import wrap, calculate_conservation_score, overlapping_finditer

# i-motif scoring - unified Hunter-style framework
def imotif_score(seq):
    """
    Unified Hunter-style scoring for i-Motifs using consistent computational logic.
    Uses the same core framework as G4Hunter but optimized for C-rich structures.
    """
    from .shared_utils import unified_hunter_score
    return unified_hunter_score(seq, target_base='C', complementary_base='G')

# i-motif detection - unified framework with proper constraints
def find_imotif(seq):
    """
    Find i-motif structures using unified computational logic and proper scientific constraints.
    Implements canonical and relaxed i-motif detection based on current literature.
    
    Scientific Basis:
    - Canonical i-motif: C3+N1-7 loops (Zeraati et al. Nat Chem 2018)
    - Relaxed i-motif: C3+N1-12 loops (extended for biological relevance)
    - Uses unified Hunter-style scoring for consistency
    """
    from .shared_utils import calculate_structural_factor
    results = []
    
    # Canonical i-motif pattern: C3+ with loops 1-7 nucleotides
    canonical_pattern = r"(?=(C{3,}[ATGC]{1,7}C{3,}[ATGC]{1,7}C{3,}[ATGC]{1,7}C{3,}))"
    for m in overlapping_finditer(canonical_pattern, seq):
        motif_seq = m.group(1)
        im_score = imotif_score(motif_seq)
        if im_score >= 0.4:  # Threshold for i-motif propensity
            # Analyze loop lengths
            c_run_spans = [match.span() for match in re.finditer(r"C{3,}", motif_seq)]
            loops = [c_run_spans[i+1][0] - c_run_spans[i][1] for i in range(len(c_run_spans)-1)]
            
            # Calculate structural factor
            structural_factor = calculate_structural_factor(motif_seq, "Canonical i-motif", loops)
            
            conservation_result = calculate_conservation_score(motif_seq, "i-Motif")
            c_runs = [len(r) for r in re.findall(r"C{3,}", motif_seq)]
            c_fraction = motif_seq.count('C') / len(motif_seq)
            
            results.append({
                "Sequence Name": "", "Class": "i-motif family", "Subtype": "Canonical i-motif",
                "Start": m.start() + 1, "End": m.start() + len(motif_seq), "Length": len(motif_seq),
                "Sequence": wrap(motif_seq), "ScoreMethod": "iM_UnifiedHunter_raw", 
                "Score": float(im_score * len(motif_seq) * structural_factor), 
                "iMotif_Mean": float(im_score), "Structural_Factor": round(structural_factor, 3),
                "C_Run_Count": len(c_runs), "C_Run_Sum": sum(c_runs) if c_runs else 0,
                "C_Fraction": round(c_fraction, 3), "Conservation_Score": float(conservation_result["enrichment_score"]),
                "Conservation_P_Value": float(conservation_result["p_value"]),
                "Conservation_Significance": conservation_result["significance"],
                "Loop_Lengths": loops, "Arms/Repeat Unit/Copies": "", "Spacer": ""
            })
    
    # Relaxed i-motif pattern: C3+ with loops 8-12 nucleotides
    relaxed_pattern = r"(?=(C{3,}[ATGC]{8,12}C{3,}[ATGC]{8,12}C{3,}[ATGC]{8,12}C{3,}))"
    for m in overlapping_finditer(relaxed_pattern, seq):
        motif_seq = m.group(1)
        im_score = imotif_score(motif_seq)
        if im_score >= 0.3:  # Lower threshold for relaxed structures
            # Analyze loop lengths
            c_run_spans = [match.span() for match in re.finditer(r"C{3,}", motif_seq)]
            loops = [c_run_spans[i+1][0] - c_run_spans[i][1] for i in range(len(c_run_spans)-1)]
            
            # Calculate structural factor
            structural_factor = calculate_structural_factor(motif_seq, "Relaxed i-motif", loops)
            
            conservation_result = calculate_conservation_score(motif_seq, "i-Motif")
            c_runs = [len(r) for r in re.findall(r"C{3,}", motif_seq)]
            c_fraction = motif_seq.count('C') / len(motif_seq)
            
            results.append({
                "Sequence Name": "", "Class": "i-motif family", "Subtype": "Relaxed i-motif",
                "Start": m.start() + 1, "End": m.start() + len(motif_seq), "Length": len(motif_seq),
                "Sequence": wrap(motif_seq), "ScoreMethod": "iM_UnifiedHunter_raw", 
                "Score": float(im_score * len(motif_seq) * structural_factor), 
                "iMotif_Mean": float(im_score), "Structural_Factor": round(structural_factor, 3),
                "C_Run_Count": len(c_runs), "C_Run_Sum": sum(c_runs) if c_runs else 0,
                "C_Fraction": round(c_fraction, 3), "Conservation_Score": float(conservation_result["enrichment_score"]),
                "Conservation_P_Value": float(conservation_result["p_value"]),
                "Conservation_Significance": conservation_result["significance"],
                "Loop_Lengths": loops, "Arms/Repeat Unit/Copies": "", "Spacer": ""
            })
    
    return results

# AC-motif scoring - optimized
def ac_motif_score(seq):
    """Enhanced AC-motif scoring using alternating A/C bias analysis."""
    if len(seq) == 0: return 0.0
    
    # Count A and C runs of length 3+
    a3_runs = len(re.findall(r"A{3,}", seq))
    c3_runs = len(re.findall(r"C{3,}", seq))
    ac_fraction = (seq.count('A') + seq.count('C')) / len(seq)
    
    # Boundary and alternation bonuses
    boundary_score = (2 if seq.startswith('AAA') else 0) + (2 if seq.endswith(('CCC', 'AAA')) else 0)
    alternation_bonus = min(a3_runs, c3_runs) * 1.5 if a3_runs >= 1 and c3_runs >= 3 else 0
    
    return (len(seq) * ac_fraction * 0.5) + boundary_score + alternation_bonus + (a3_runs + c3_runs) * 2

# AC-motif detection - based on Hur et al. NAR 2021 
def find_ac_motifs(seq):
    """
    Find AC-motif structures based on Hur et al. NAR 2021 patterns.
    
    Scientific Implementation:
    - Based on reviewed regular expression from Hur et al. NAR 2021
    - Applies correct pattern and filtering logic for AC-motifs
    - Uses unified computational framework for consistency
    
    References:
    - Hur et al. NAR 2021 (AC motif patterns and formation requirements)
    """
    results, used_positions = [], set()
    
    # Enhanced patterns for AC-motifs based on Hur et al. NAR 2021
    # Pattern 1: Classical AC-motif with A-rich and C-rich segments
    # More scientifically accurate pattern based on literature
    patterns = [
        # Primary pattern: A3+ tract followed by C3+ tract patterns
        r"A{3,}[ATGC]{1,15}C{3,}[ATGC]{1,15}C{3,}[ATGC]{1,15}C{3,}",
        # Secondary pattern: AC-rich regions with specific constraints
        r"[AC]{4}[ATGC]{0,3}[AC]{4}[ATGC]{0,3}[AC]{4}[ATGC]{0,3}[AC]{4}",
        # Tertiary pattern: Alternating AC motifs
        r"(A{2,4}C{2,4}){3,}",
    ]
    
    for pattern_idx, pattern in enumerate(patterns):
        for m in re.finditer(pattern, seq, re.IGNORECASE):
            motif_seq = m.group(0)
            
            # Check for overlap with existing motifs
            current_positions = set(range(m.start(), m.end()))
            if len(current_positions.intersection(used_positions)) / len(current_positions) > 0.3:
                continue
            
            # Enhanced AC-motif scoring using unified framework
            score = ac_motif_score(motif_seq)
            
            # Apply more stringent filtering based on Hur et al. criteria
            a_content = motif_seq.count('A') / len(motif_seq)
            c_content = motif_seq.count('C') / len(motif_seq)
            ac_content = a_content + c_content
            
            # Minimum requirements based on literature
            if score < 10.0 or ac_content < 0.6 or len(motif_seq) < 15:
                continue
                
            # Calculate detailed metrics
            a3_runs = len(re.findall(r"A{3,}", motif_seq))
            c3_runs = len(re.findall(r"C{3,}", motif_seq))
            conservation_result = calculate_conservation_score(motif_seq, "AC-Motif")
            
            # Calculate structural factor for AC-motifs
            from .shared_utils import calculate_structural_factor
            structural_factor = calculate_structural_factor(motif_seq, "AC-motif")
            
            subtype = f"AC-motif_Pattern{pattern_idx + 1}"
            
            results.append({
                "Sequence Name": "", "Class": "i-motif family", "Subtype": "AC-motif",
                "Start": m.start() + 1, "End": m.end(), "Length": len(motif_seq),
                "Sequence": wrap(motif_seq), "ScoreMethod": "AC_UnifiedFramework_raw", 
                "Score": float(score * structural_factor),
                "AC_Fraction": round(ac_content, 3), "A_Fraction": round(a_content, 3),
                "C_Fraction": round(c_content, 3), "A3_Runs": a3_runs, "C3_Runs": c3_runs,
                "Structural_Factor": round(structural_factor, 3),
                "Pattern_Type": subtype,
                "Conservation_Score": float(conservation_result["enrichment_score"]),
                "Conservation_P_Value": float(conservation_result["p_value"]),
                "Conservation_Significance": conservation_result["significance"],
                "Arms/Repeat Unit/Copies": f"A3={a3_runs};C3={c3_runs}", "Spacer": ""
            })
            used_positions.update(current_positions)
    
    return results