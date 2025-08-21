"""
Category 3: Slipped DNA Detection Module
======================================

This module implements detection algorithms for slipped DNA structures,
including direct repeats and short tandem repeats (STRs).

Scientific Basis:
- DNA slippage occurs during replication at direct repeats and STRs
- Forms looped-out structures that can lead to genetic instability
- AT-rich repeats are more prone to slippage due to lower melting temperature

References:
- Kunkel & Bebenek (2000) Annu Rev Biochem
- Wells (2007) Trends Biochem Sci

Author: Dr. Venkata Rajesh Yella
Updated: 2024
"""

from .shared_utils import wrap, calculate_conservation_score

# Slipped DNA detection - optimized for performance and reduced overlaps
def find_slipped_dna(seq):
    """Find slipped DNA structures including direct repeats and STRs with enhanced scoring."""
    results = []
    min_len_dr, max_len_dr = 10, 300
    min_score_threshold = 25.0  # Minimum score to avoid excessive low-quality matches
    used_positions = set()  # Track used positions to prevent excessive overlap
    
    # Direct repeats - optimized to reduce overlapping matches
    for i in range(len(seq) - min_len_dr * 2 + 1):
        if i in used_positions: continue  # Skip if position already covered
            
        for l in range(min_len_dr, min(max_len_dr+1, (len(seq)-i)//2+1)):
            repeat = seq[i:i+l]
            if seq[i+l:i+2*l] == repeat:
                # Enhanced scoring: unit_len × composition weight (AT-rich repeats more flexible)
                at_frac = (repeat.count('A') + repeat.count('T')) / max(1, len(repeat))
                score = 2*l * (1.0 + 0.5*at_frac)
                
                # Only keep high-scoring matches and prevent excessive overlap
                if score >= min_score_threshold:
                    current_positions = set(range(i, i+2*l))
                    overlap_ratio = len(current_positions.intersection(used_positions)) / len(current_positions)
                    
                    if overlap_ratio < 0.3:  # Allow some overlap but not excessive
                        conservation_result = calculate_conservation_score(repeat+repeat, "Slipped DNA")
                        
                        results.append({
                            "Sequence Name": "", "Class": "Slipped DNA", "Subtype": "Slipped DNA [Direct Repeat]",
                            "Start": i+1, "End": i+2*l, "Length": 2*l, "Sequence": wrap(repeat+repeat),
                            "ScoreMethod": "DR_Composition_raw", "Score": float(score), "AT_Fraction": round(at_frac, 3),
                            "Conservation_Score": float(conservation_result["enrichment_score"]),
                            "Conservation_P_Value": float(conservation_result["p_value"]),
                            "Conservation_Significance": conservation_result["significance"],
                            "Arms/Repeat Unit/Copies": f"UnitLen={l};Copies=2", "Spacer": ""
                        })
                        used_positions.update(range(i, i+l))  # Mark core positions as used
                        break  # Found a good match at this position, move to next
    
    # Short Tandem Repeats (STRs) - Enhanced algorithm
    min_unit_str = 1
    max_unit_str = 6
    min_reps_str = 5
    min_len_str = 15
    i = 0
    n = len(seq)
    
    while i < n - min_unit_str * min_reps_str + 1:
        found = False
        for unit in range(min_unit_str, max_unit_str+1):
            if i + unit * min_reps_str > n:
                continue
            repeat_unit = seq[i:i+unit]
            if 'n' in repeat_unit.lower():
                continue
            reps = 1
            while (i + reps*unit + unit <= n and seq[i + reps*unit:i + (reps+1)*unit] == repeat_unit):
                reps += 1
            if reps >= min_reps_str and reps*unit >= min_len_str:
                remainder = 0
                rs = i + reps*unit
                re_idx = rs
                while (re_idx < n and seq[re_idx] == repeat_unit[re_idx % unit]):
                    remainder += 1
                    re_idx += 1
                full_len = reps*unit + remainder
                gc_frac = (repeat_unit.count('G') + repeat_unit.count('C')) / max(1, len(repeat_unit))
                score = full_len * (1.0 + 0.3*gc_frac) * (reps ** 0.5)
                
                # Calculate conservation score
                str_seq = seq[i:i + full_len]
                conservation_result = calculate_conservation_score(str_seq, "Slipped DNA")
                conservation_score = conservation_result["enrichment_score"]
                
                results.append({
                    "Sequence Name": "",
                    "Class": "Slipped DNA",
                    "Subtype": "Slipped DNA [STR]",
                    "Start": i+1,
                    "End": i + full_len,
                    "Length": full_len,
                    "Unit": repeat_unit,
                    "Copies": reps,
                    "Sequence": wrap(str_seq),
                    "ScoreMethod": "STR_Enhanced_raw",
                    "Score": float(score),
                    "GC_Fraction": round(gc_frac, 3),
                    "Conservation_Score": float(conservation_score),
                    "Conservation_P_Value": float(conservation_result["p_value"]),
                    "Conservation_Significance": conservation_result["significance"],
                    "Arms/Repeat Unit/Copies": f"Unit={repeat_unit};Copies={reps}",
                    "Spacer": ""
                })
                i = i + full_len - 1
                found = True
                break
        if not found:
            i += 1
    
    return results