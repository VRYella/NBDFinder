"""
Category 2: Hairpin and Cruciform Detection Module
=================================================

This module implements detection algorithms for cruciform structures
formed by inverted repeats that can extrude into four-way junctions.

Scientific Basis:
- Cruciform structures form at palindromic sequences (inverted repeats)
- AT-rich sequences easier to extrude due to lower melting temperature
- Important for gene regulation and recombination

Author: Dr. Venkata Rajesh Yella
Updated: 2024
"""

from .shared_utils import wrap, calculate_conservation_score, reverse_complement

def find_cruciform(seq):
    """
    Find cruciform structures formed by inverted repeats.
    
    Enhanced Scoring System:
    1. Arm length scaling with AT-richness bonus
    2. Spacer length penalty for flexibility
    3. Conservation analysis for evolutionary significance
    
    Parameters:
    seq (str): DNA sequence
    
    Returns:
    list: List of cruciform motif dictionaries
    """
    results = []
    n = len(seq)
    for i in range(n - 2*10):
        for arm_len in range(10, min(101, (n-i)//2)):
            for spacer_len in range(0, 4):
                arm = seq[i:i+arm_len]
                rev_arm = reverse_complement(arm)
                mid = i + arm_len + spacer_len
                if mid + arm_len > n:
                    continue
                candidate = seq[mid:mid+arm_len]
                if candidate == rev_arm:
                    full = seq[i:mid+arm_len]
                    # Enhanced scoring: arm_len × AT_richness_bonus - spacer_penalty
                    # AT-rich sequences easier to extrude
                    at_frac = (arm.count('A') + arm.count('T')) / arm_len
                    score = arm_len * (1.0 + 0.5*at_frac) - 2.0*spacer_len
                    
                    # Calculate conservation score
                    conservation_result = calculate_conservation_score(full, "Cruciform")
                    conservation_score = conservation_result["enrichment_score"]
                    
                    results.append({
                        "Sequence Name": "",
                        "Class": "Cruciform",
                        "Subtype": f"Inverted_Repeat_spacer{spacer_len}",
                        "Start": i+1,
                        "End": mid+arm_len,
                        "Length": len(full),
                        "Sequence": wrap(full),
                        "ScoreMethod": "IR_AT_Enhanced_raw",
                        "Score": float(score),
                        "Arm_Length": arm_len,
                        "AT_Fraction": round(at_frac, 3),
                        "Conservation_Score": float(conservation_score),
                        "Conservation_P_Value": float(conservation_result["p_value"]),
                        "Conservation_Significance": conservation_result["significance"],
                        "Arms/Repeat Unit/Copies": f"Arms={arm_len}",
                        "Spacer": str(spacer_len)
                    })
    return results