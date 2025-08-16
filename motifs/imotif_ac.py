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

def imotif_score(seq):
    """
    G4Hunter-style scoring for i-Motifs: C vs G bias calculation with reverse logic.
    
    Scientific Basis: i-Motifs are C-rich structures complementary to G-quadruplexes.
    This scoring system uses G4Hunter run-based logic but with reversed polarity.
    Reference: Zeraati et al. (2018) Nat Chem; G4Hunter methodology adapted for i-motifs
    
    Enhanced Scoring System:
    1. C-run scoring with positive values (reverse of G4Hunter)
    2. G-run scoring with negative values for competition
    3. Run length scaling up to 4 bases maximum
    
    Parameters:
    seq (str): DNA sequence to score
    
    Returns:
    float: i-Motif Hunter mean score (run-based C vs G bias)
    """
    seq = seq.upper()
    scores = []
    i = 0
    while i < len(seq):
        # C-run: assign +score for each C in the run (reverse of G4Hunter)
        if seq[i] == 'C':
            run = 1
            while i + run < len(seq) and seq[i + run] == 'C':
                run += 1
            s = min(run, 4)
            scores += [s] * run
            i += run
        # G-run: assign -score for each G in the run (reverse of G4Hunter)
        elif seq[i] == 'G':
            run = 1
            while i + run < len(seq) and seq[i + run] == 'G':
                run += 1
            s = -min(run, 4)
            scores += [s] * run
            i += run
        # A/T: assign 0
        else:
            scores.append(0)
            i += 1
    return sum(scores) / len(scores) if scores else 0.0

def find_imotif(seq):
    """
    Find i-motif structures using overlapping pattern detection.
    
    Enhanced Scoring System:
    1. G4Hunter-style i-motif scoring for C vs G bias
    2. Loop length analysis for structural classification
    3. Composite scoring combining bias and length
    4. Conservation analysis for evolutionary significance
    
    Parameters:
    seq (str): DNA sequence
    
    Returns:
    list: List of i-motif motif dictionaries
    """
    results = []
    pattern = r"(?=(C{3,}\w{1,12}C{3,}\w{1,12}C{3,}\w{1,12}C{3,}))"
    for m in overlapping_finditer(pattern, seq):
        motif_seq = m.group(1)
        im_score = imotif_score(motif_seq)
        if im_score >= 0.4:  # Threshold for i-motif propensity (positive C bias)
            # Analyze loop lengths for subtype classification
            c_run_spans = [match.span() for match in re.finditer(r"C{3,}", motif_seq)]
            loops = []
            for i in range(len(c_run_spans)-1):
                loop_start = c_run_spans[i][1]
                loop_end = c_run_spans[i+1][0]
                loops.append(loop_end - loop_start)
            
            # Classify based on loop characteristics
            if loops and all(1 <= l <= 7 for l in loops):
                subtype = "Canonical_iMotif"
            elif loops and any(8 <= l <= 12 for l in loops):
                subtype = "LongLoop_iMotif"
            else:
                subtype = "Other_iMotif"
            
            # Calculate conservation score
            conservation_result = calculate_conservation_score(motif_seq, "i-Motif")
            conservation_score = conservation_result["enrichment_score"]
            
            # Calculate composite score: iM_score × length (similar to G4Hunter approach)
            composite_score = im_score * len(motif_seq)
            
            # Additional scoring details
            c_runs = [len(r) for r in re.findall(r"C{3,}", motif_seq)]
            c_fraction = motif_seq.count('C') / len(motif_seq)
            
            results.append({
                "Sequence Name": "",
                "Class": "i-Motif",
                "Subtype": subtype,
                "Start": m.start() + 1,
                "End": m.start() + len(motif_seq),
                "Length": len(motif_seq),
                "Sequence": wrap(motif_seq),
                "ScoreMethod": "iM_G4HunterStyle_raw",
                "Score": float(composite_score),
                "iMotif_Mean": float(im_score),
                "C_Run_Count": len(c_runs),
                "C_Run_Sum": sum(c_runs) if c_runs else 0,
                "C_Fraction": round(c_fraction, 3),
                "Conservation_Score": float(conservation_score),
                "Conservation_P_Value": float(conservation_result["p_value"]),
                "Conservation_Significance": conservation_result["significance"],
                "Loop_Lengths": loops if loops else [],
                "Arms/Repeat Unit/Copies": "",
                "Spacer": ""
            })
    return results

def ac_motif_score(seq):
    """
    Enhanced AC-motif scoring using alternating A/C bias analysis.
    
    Scientific Basis: AC-motifs show alternating purine-pyrimidine patterns that can
    form non-canonical structures. Score based on A/C alternation and structural potential.
    Reference: Kocsis et al. (2021) Int J Mol Sci; AC-motif structural analysis
    
    Enhanced Scoring System:
    1. A/C content scaling with sequence length
    2. Boundary recognition for proper AC-motif structure
    3. Alternation pattern bonuses for canonical patterns
    4. Run count bonuses for structural elements
    
    Parameters:
    seq (str): DNA sequence to score
    
    Returns:
    float: AC-motif propensity score
    """
    if len(seq) == 0:
        return 0.0
    
    # Count A and C runs of length 3+
    a3_runs = len(re.findall(r"A{3,}", seq))
    c3_runs = len(re.findall(r"C{3,}", seq))
    
    # Calculate A/C fraction
    ac_fraction = (seq.count('A') + seq.count('C')) / len(seq)
    
    # Boundary bonus for proper AC-motif structure
    boundary_score = 0
    if seq.startswith('AAA'):
        boundary_score += 2
    if seq.endswith('CCC') or seq.endswith('AAA'):
        boundary_score += 2
    
    # Alternation pattern bonus - check for spacing between A3 and C3 runs
    alternation_bonus = 0
    if a3_runs >= 1 and c3_runs >= 3:  # Canonical AC-motif pattern
        alternation_bonus = min(a3_runs, c3_runs) * 1.5
    
    # Final score: length-normalized AC content + structural bonuses
    score = (len(seq) * ac_fraction * 0.5) + boundary_score + alternation_bonus + (a3_runs + c3_runs) * 2
    
    return score

def find_ac_motifs(seq):
    """
    Find AC-motif structures using pattern matching with enhanced scoring.
    
    Enhanced Scoring System:
    1. Multiple pattern recognition for AC-motif variants
    2. A3/C3 run analysis for structural classification
    3. Overlap prevention for quality control
    4. Conservation analysis for evolutionary significance
    
    Parameters:
    seq (str): DNA sequence
    
    Returns:
    list: List of AC-motif motif dictionaries
    """
    results = []
    used_positions = set()
    
    # Enhanced patterns for AC-motifs with varying structures
    patterns = [
        r"AAA[ATGC]{1,20}CCC[ATGC]{1,20}CCC[ATGC]{1,20}CCC",  # Strict A3-C3 consensus
        r"[AC]{15,40}"  # Relaxed AC-rich regions
    ]
    
    for pattern_idx, pattern in enumerate(patterns):
        for m in re.finditer(pattern, seq, re.IGNORECASE):
            motif_seq = m.group(0)
            
            # Check for overlap with existing motifs
            current_positions = set(range(m.start(), m.end()))
            overlap_ratio = len(current_positions.intersection(used_positions)) / len(current_positions)
            if overlap_ratio > 0.3:  # Skip if >30% overlap
                continue
                
            score = ac_motif_score(motif_seq)
            
            # Only keep high-scoring matches
            if score < 15.0:  # Minimum score threshold
                continue
            
            # Calculate detailed metrics
            a3_runs = len(re.findall(r"A{3,}", motif_seq))
            c3_runs = len(re.findall(r"C{3,}", motif_seq))
            ac_fraction = (motif_seq.count('A') + motif_seq.count('C')) / len(motif_seq)
            
            # Calculate conservation score
            conservation_result = calculate_conservation_score(motif_seq, "AC-Motif")
            conservation_score = conservation_result["enrichment_score"]
            
            # Determine subtype based on structure
            if motif_seq.startswith('AAA') and c3_runs >= 3:
                subtype = "A3-C3_Consensus"
            elif motif_seq.startswith('CCC') and a3_runs >= 1:
                subtype = "C3-A3_Consensus"
            elif pattern_idx == 0:
                subtype = "Strict_AC_Motif"
            else:
                subtype = "Relaxed_AC_Motif"
            
            results.append({
                "Sequence Name": "",
                "Class": "AC-Motif",
                "Subtype": subtype,
                "Start": m.start() + 1,
                "End": m.end(),
                "Length": len(motif_seq),
                "Sequence": wrap(motif_seq),
                "ScoreMethod": "AC_StructuralAnalysis_raw",
                "Score": float(score),
                "AC_Fraction": round(ac_fraction, 3),
                "A3_Runs": a3_runs,
                "C3_Runs": c3_runs,
                "Conservation_Score": float(conservation_score),
                "Conservation_P_Value": float(conservation_result["p_value"]),
                "Conservation_Significance": conservation_result["significance"],
                "Arms/Repeat Unit/Copies": f"A3={a3_runs};C3={c3_runs}",
                "Spacer": ""
            })
            
            # Mark positions as used
            used_positions.update(current_positions)
    
    return results