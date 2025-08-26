"""
Category 3:  Slipped DNA Detection Module
=================================================
This module implements detection algorithms for slipped DNA structures, including direct repeats and STRs

Scientific Basis:

DNA slippage occurs during replication at direct repeats and STRs, creating looped-out slip-outs that can expand or contract repeat length.
[Replication slippage overview: Viguera et al., EMBO J 2001; Slipped-strand mispairing concept]

Slippage frequency scales with repeat copy number and is higher for shorter repeat units; perfect-identity detection with deterministic thresholds 
yields reproducible catalogs. [STR mutation/copy-size dependence: Gymrek et al., PNAS 2022]

Composition (AT/GC) is retained as annotation only; scoring emphasizes unit length, copy number, and total length to avoid composition bias 
while remaining scientifically grounded. [Catalog-style non-B motif annotation: nBMST/Non-B DB]

References:
Kunkel & Bebenek (2000) Annu Rev Biochem — replication fidelity and slippage context.
Wells (2007) Trends Biochem Sci — non-B DNA and genome instability overview.
Viguera et al. (2001) EMBO J — polymerase pausing/dissociation in slippage between direct repeats.
Non-B DNA Motif Search Tool/Non-B DB (Nucleic Acids Res, 2011/2010) — deterministic non-B motif cataloging.
Gymrek et al. (2022) PNAS — STR mutation/selection processes (copy number, unit size effects).

Internal deterministic scoring helpers
Rationale (annotation):
- Unit-length weight: shorter units are more slippage-prone; approximate with l^(-alpha) discretized
for reproducibility (keeps scores small and comparable). [STR instability vs. unit size]
- Copy-number bonus: monotonic, capped; higher copy number increases slippage likelihood. [STR copy effect]
- Length thresholds (30/50/100 nt): simple, transparent steps capturing length-dependent instability
observed across non-B motifs and repeats; avoids overfitting. [Catalog practice/length sensitivity]
"""
from .shared_utils import wrap, calculate_conservation_score

def find_slipped_dna_advanced(seq):
    """
    Detects slipped DNA motifs (DR, STR) in deterministic, mismatch-free fashion.
    Each call is scored and annotated with composition and conservation metrics.
    Returns: List of non-overlapping, catalog-compliant motif dictionaries.
    """
    results, preliminary = [], []; n = len(seq)
    # --- DR parameters ---
    min_len_dr, max_len_dr, max_spacer = 10, 300, 5
    # --- STR parameters ---
    min_unit_str, max_unit_str, min_len_str = 1, 6, 15
    min_reps_by_unit = {1:6, 2:5, 3:4, 4:3, 5:3, 6:3}
    min_score_threshold = 10.0

    # --- Direct Repeats (DR) block: strict catalog-mismatch-free search ---
    for i in range(0, max(0, n - min_len_dr * 2 + 1)):
        if seq[i].lower() == 'n': continue  # skip ambiguous base
        made_call = False
        max_l_here = min(max_len_dr, (n - i) // 2)
        for l in range(max_l_here, min_len_dr - 1, -1):
            repeat = seq[i:i+l]
            if 'n' in repeat.lower(): continue
            max_sp = min(max_spacer, n - (i + 2*l))
            for sp in range(0, max_sp + 1):
                j = i + l + sp
                if j + l > n: continue
                if seq[j:j+l] != repeat: continue
                # Found DR; catalog extension logic
                copies, remainder, end = 2, 0, j + l
                if sp == 0:
                    k = j + l
                    while k + l <= n and seq[k:k+l] == repeat: copies += 1; k += l
                    rs, re = i, k
                    while re < n and seq[rs] == seq[re]: remainder += 1; rs += 1; re += 1
                    end = k if remainder == 0 else re
                else:
                    rs, re = i, j + l
                    while re < n and seq[rs] == seq[re]: remainder += 1; rs += 1; re += 1
                    end = j + l if remainder == 0 else re
                total_len = end - i
                score = _score_DR(l, copies, total_len, sp)
                if score >= min_score_threshold:
                    dr_seq = seq[i:end]
                    conservation_result = calculate_conservation_score(dr_seq, "Slipped DNA")
                    preliminary.append({
                        "Sequence Name": "",
                        "Class": "Slipped DNA",
                        "Subtype": "Slipped DNA [Direct Repeat]",
                        "Start": i+1,
                        "End": end,
                        "Length": total_len,
                        "Sequence": wrap(dr_seq),
                        "ScoreMethod": "DR_PerfectBlock_v1",
                        "Score": float(score),
                        "AT_Fraction": round((repeat.count('A')+repeat.count('T'))/max(1,len(repeat)), 3),
                        "Conservation_Score": float(conservation_result.get("enrichment_score", 0.0)),
                        "Conservation_P_Value": float(conservation_result.get("p_value", 1.0)),
                        "Conservation_Significance": conservation_result.get("significance", ""),
                        "Arms/Repeat Unit/Copies": f"UnitLen={l};Copies={copies}",
                        "Spacer": str(sp)
                    }); made_call = True; break
            if made_call: break

    # --- STR block: strict mismatch-free search with catalog copy/length threshold ---
    i = 0
    while i <= n - min_len_str:
        if seq[i].lower() == 'n': i += 1; continue
        found = False
        for unit in range(min_unit_str, max_unit_str+1):
            if i + unit > n: break
            repeat_unit = seq[i:i+unit]
            if 'n' in repeat_unit.lower(): continue
            
            # Exclude GAA/TTC patterns and their variants - these should be handled by sticky DNA detection
            if (repeat_unit.upper() in ['GAA', 'TTC', 'AAG', 'AGG', 'GGA', 'TCT', 'CTC', 'CTT'] and 
                any(x in seq[i:i+50].upper() for x in ['GAAGAA', 'TTCTTC'])):
                continue
            reps, j = 1, i+unit
            while j + unit <= n and seq[j:j+unit] == repeat_unit: reps += 1; j += unit
            remainder, rs, re_idx = 0, i, j
            while re_idx < n and seq[rs] == seq[re_idx]: remainder += 1; rs += 1; re_idx += 1
            full_len = reps*unit + remainder
            if reps >= min_reps_by_unit.get(unit, 999) and full_len >= min_len_str:
                str_seq = seq[i:i+full_len]
                score = _score_STR(unit, reps, full_len)
                if score >= min_score_threshold:
                    conservation_result = calculate_conservation_score(str_seq, "Slipped DNA")
                    results.append({
                        "Sequence Name": "",
                        "Class": "Slipped DNA",
                        "Subtype": "Slipped DNA [STR]",
                        "Start": i+1,
                        "End": i+full_len,
                        "Length": full_len,
                        "Unit": repeat_unit,
                        "Copies": reps,
                        "Sequence": wrap(str_seq),
                        "ScoreMethod": "STR_PerfectBlock_v1",
                        "Score": float(score),
                        "GC_Fraction": round((repeat_unit.count('G')+repeat_unit.count('C'))/max(1,len(repeat_unit)), 3),
                        "Conservation_Score": float(conservation_result.get("enrichment_score", 0.0)),
                        "Conservation_P_Value": float(conservation_result.get("p_value", 1.0)),
                        "Conservation_Significance": conservation_result.get("significance", ""),
                        "Arms/Repeat Unit/Copies": f"Unit={repeat_unit};Copies={reps}",
                        "Spacer": ""
                    })
                i = i + max(full_len, unit); found = True; break
        if not found: i += 1

    # --- Global non-overlap: ensures one motif per region (catalog standard) ---
    if preliminary:
        all_hits = preliminary + results
        results = _non_overlap_selection(all_hits)
    return results

def _score_DR(l, copies, total_len, spacer):
    """
    Scoring for Direct Repeats per catalog rules.
    Emphasizes unit length, copy number, and total length.
    """
    # Unit-length weight: shorter units are more slippage-prone
    unit_weight = max(1.0, 5.0 / l) if l > 0 else 1.0
    
    # Copy-number bonus: monotonic, capped
    copy_bonus = min(copies * 2.0, 10.0)
    
    # Length thresholds: simple, transparent steps
    length_bonus = 0.0
    if total_len >= 100:
        length_bonus = 3.0
    elif total_len >= 50:
        length_bonus = 2.0
    elif total_len >= 30:
        length_bonus = 1.0
    
    # Spacer penalty
    spacer_penalty = spacer * 0.5
    
    return unit_weight + copy_bonus + length_bonus - spacer_penalty

def _score_STR(unit, reps, full_len):
    """
    Scoring for Short Tandem Repeats per catalog rules.
    Emphasizes unit length, repeat count, and total length.
    """
    # Unit-length weight: shorter units are more slippage-prone (discretized)
    unit_weight = {1: 4.0, 2: 3.0, 3: 2.5, 4: 2.0, 5: 1.5, 6: 1.0}.get(unit, 1.0)
    
    # Copy-number bonus: monotonic, capped
    copy_bonus = min(reps * 1.5, 12.0)
    
    # Length thresholds: simple, transparent steps
    length_bonus = 0.0
    if full_len >= 100:
        length_bonus = 3.0
    elif full_len >= 50:
        length_bonus = 2.0
    elif full_len >= 30:
        length_bonus = 1.0
    
    return unit_weight + copy_bonus + length_bonus

def _non_overlap_selection(motif_list):
    """
    Select non-overlapping intervals from a motif list by descending score.
    Standard greedy approach for interval selection.
    """
    sorted_motifs = sorted(motif_list, key=lambda m: m['Score'], reverse=True)
    selected = []
    covered = set()
    for motif in sorted_motifs:
        region = set(range(motif['Start'], motif['End']+1))
        if not region & covered:
            selected.append(motif)
            covered |= region
    return selected

def find_slipped_dna(seq):
    """
    Wrapper function for backwards compatibility.
    Calls find_slipped_dna_advanced with the sequence.
    """
    return find_slipped_dna_advanced(seq)
