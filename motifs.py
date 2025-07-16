import re
import numpy as np
from typing import List, Dict, Tuple
from collections import defaultdict, Counter
import random

# ========== UTILS FUNCTIONS ==========
def parse_fasta(fasta_str: str) -> str:
    return "".join([line.strip() for line in fasta_str.split('\n') if not line.startswith(">")]).upper().replace(" ", "").replace("U", "T")

def wrap(seq: str, width: int = 60) -> str:
    return "\n".join([seq[i:i+width] for i in range(0, len(seq), width)])

def gc_content(seq: str) -> float:
    return 100 * (seq.count('G') + seq.count('C')) / max(1, len(seq))

def reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ATGC", "TACG"))[::-1]

def is_palindrome(seq: str) -> bool:
    return seq == reverse_complement(seq)

def g4hunter_score(seq: str) -> float:
    scores = []
    for c in seq.upper():
        if c == 'G':
            scores.append(1)
        elif c == 'C':
            scores.append(-1)
        else:
            scores.append(0)
    return np.mean(scores) if scores else 0

def percentileofscore(a, score, kind='rank'):
    a = np.asarray(a)
    if len(a) == 0: return 0.0
    if kind == 'rank': return (sum(a <= score) / len(a) * 100)
    elif kind == 'strict': return (sum(a < score) / len(a)) * 100
    elif kind == 'weak': return (sum(a <= score) / len(a)) * 100
    elif kind == 'mean': return (sum(a < score) + sum(a <= score)) / (2 * len(a)) * 100
    else: raise ValueError("kind must be 'rank', 'strict', 'weak' or 'mean'")

# ========== MOTIF DETECTION ==========
def overlapping_finditer(pattern, seq):
    for m in re.compile(pattern, re.IGNORECASE).finditer(seq):
        yield m

def all_motifs(seq):
    """
    Identifies all non-B DNA motifs in the input sequence, including:
    - Curved DNA motifs (Global: phased A/T tracts, Local: A7/T7, no TA step)
    - Z-DNA, slipped DNA, R-loops, cruciforms, triplexes, G-quadruplexes, i-motifs, hybrids, sticky DNA, etc.

    References:
    - Crothers DM et al. (1992) DNA bending by A-tracts. Science.
    - Brukner I et al. (1995) Curvature of DNA: phasing of A-tracts. J Biomol Struct Dyn.
    - Trifonov EN (1980) Sequence-dependent deformational anisotropy of chromatin DNA.
    """
    if not seq or not re.match("^[ATGC]+$", seq, re.IGNORECASE):
        return []
    seq = seq.upper()
    results = (
        find_curved_DNA(seq) +
        find_zdna(seq) +
        find_slipped_dna(seq) +
        find_rlfs(seq) +
        find_cruciform(seq) +
        find_hdna(seq) +
        find_gtriplex(seq) +
        find_gquadruplex(seq) +
        find_relaxed_gquadruplex(seq) +
        find_bulged_gquadruplex(seq) +
        find_bipartite_gquadruplex(seq) +
        find_multimeric_gquadruplex(seq) +
        find_imotif(seq) +
        find_hybrids(seq) +
        find_sticky_dna(seq)
    )
    return [m for m in results if validate_motif(m, len(seq))]
#####################################################################################################################
#####################################################################################################################
import re

def wrap(seq, width=60):
    """Wraps sequence for display purposes (optional)."""
    return '\n'.join([seq[i:i+width] for i in range(0, len(seq), width)])

def find_AT_tracts_no_TA(seq: str, min_len: int) -> list:
    """
    Finds runs of A or T bases (≥min_len) with NO TA steps.
    Returns list of (start, end, sequence).
    """
    tracts = []
    i = 0
    while i < len(seq):
        if seq[i] in "AT":
            start = i
            tract_seq = seq[i]
            i += 1
            ta_found = False
            while i < len(seq) and seq[i] in "AT":
                if seq[i-1:i+1] == 'TA':
                    ta_found = True
                tract_seq += seq[i]
                i += 1
            if len(tract_seq) >= min_len and not ta_found:
                tracts.append((start, i-1, tract_seq))
        else:
            i += 1
    return tracts

def find_AT_tracts_relaxed(seq: str, min_len: int, max_TA: int = 0) -> list:
    """
    Finds runs of A or T bases (≥min_len), allowing at most max_TA TA steps per tract.
    Returns list of (start, end, sequence).
    NOTE: For relaxed version, max_TA=0 (NO TA allowed).
    """
    tracts = []
    i = 0
    while i < len(seq):
        if seq[i] in "AT":
            start = i
            tract_seq = seq[i]
            i += 1
            ta_count = 0
            while i < len(seq) and seq[i] in "AT":
                if seq[i-1:i+1] == 'TA':
                    ta_count += 1
                tract_seq += seq[i]
                i += 1
            if len(tract_seq) >= min_len and ta_count <= max_TA:
                tracts.append((start, i-1, tract_seq))
        else:
            i += 1
    return tracts

def curvature_score(seq):
    """
    Scores a curved DNA motif:
    Each AA, TT, or AT step = 1 point.
    Example: AAA = 2 (AA, AA), ATTA = 2 (AT, TT).
    """
    score = 0
    for i in range(len(seq) - 1):
        pair = seq[i:i+2]
        if pair in ["AA", "TT", "AT"]:
            score += 1
    return score

def find_global_curved(seq: str, min_tract_len: int = 3, min_repeats: int = 3, min_spacing: int = 8, max_spacing: int = 12, min_score: int = 6) -> tuple:
    """
    Detects global curved DNA motifs:
      - ≥3 A/T tracts (≥3 bases, NO TA step)
      - Phased at 8–12 bp center-to-center
      - Motif score (AA/TT/AT steps) must be ≥ min_score (default 6)
    """
    tracts = find_AT_tracts_no_TA(seq, min_tract_len)
    results = []
    apr_regions = []
    for i in range(len(tracts) - min_repeats + 1):
        group = [tracts[i]]
        for j in range(1, min_repeats):
            prev_center = (tracts[i + j - 1][0] + tracts[i + j - 1][1]) // 2
            curr_center = (tracts[i + j][0] + tracts[i + j][1]) // 2
            spacing = curr_center - prev_center
            if min_spacing <= spacing <= max_spacing:
                group.append(tracts[i + j])
            else:
                break
        if len(group) >= min_repeats:
            motif_seq = seq[group[0][0]:group[-1][1]+1]
            score = curvature_score(motif_seq)
            if score >= min_score:
                motif = {
                    "Class": "Curved_DNA",
                    "Subtype": "Global_Curved_Strict",
                    "Start": group[0][0] + 1,
                    "End": group[-1][1] + 1,
                    "Length": group[-1][1] - group[0][0] + 1,
                    "Sequence": wrap(motif_seq),
                    "ScoreMethod": "Strict: curvature step score",
                    "Score": score
                }
                results.append(motif)
                apr_regions.append((motif["Start"], motif["End"]))
    return results, apr_regions

def find_local_curved(seq: str, apr_regions: list, min_len: int = 7) -> list:
    """
    Detects relaxed local curved DNA motifs:
      - A/T tract of ≥7 bases (NO TA step), not overlapping global curved regions
    """
    results = []
    tracts = find_AT_tracts_relaxed(seq, min_len, max_TA=0)  # NO TA allowed even for relaxed
    for start, end, tract_seq in tracts:
        s, e = start + 1, end + 1
        if not any(r_start <= s <= r_end or r_start <= e <= r_end for r_start, r_end in apr_regions):
            results.append({
                "Class": "Curved_DNA",
                "Subtype": "Local_Curved_Relaxed",
                "Start": s,
                "End": e,
                "Length": len(tract_seq),
                "Sequence": wrap(tract_seq),
                "ScoreMethod": "Relaxed: tract length",
                "Score": len(tract_seq)
            })
    return results

def find_curved_DNA(seq: str) -> list:
    """
    Main function for curved DNA motif detection.
    - Global: ≥3 A/T tracts (≥3 bases, NO TA step), spaced 8–12 bp (center-to-center),
      curvature step score ≥ 6
    - Local: A/T tracts (≥7 bases, NO TA step), not overlapping global motifs
    """
    global_results, apr_regions = find_global_curved(seq)
    local_results = find_local_curved(seq, apr_regions)
    return global_results + local_results

############################################ 1. Curved DNA Motif: Code End  ########################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#


############################################ 2. Z-DNA  ############################################################

import numpy as np

def zdna_seeker_scoring_array(
    seq,
    GC_weight=7.0, AT_weight=0.5, GT_weight=1.25, AC_weight=1.25,
    consecutive_AT_scoring=(0.5, 0.5, 0.5, 0.5, 0.0, 0.0, -5.0, -100.0),
    mismatch_penalty_type="linear",
    mismatch_penalty_starting_value=3,
    mismatch_penalty_linear_delta=3,
    cadence_reward=0.0
):
    """
    Computes a Z-DNA propensity score for every dinucleotide in the input sequence.

    Scoring is based on the type of dinucleotide:
      - 'GC'/'CG': Strong Z-DNA former, highest score.
      - 'GT'/'TG' and 'AC'/'CA': Intermediate scores.
      - 'AT'/'TA': Weak Z-DNA former, penalized more for consecutive repeats.
      - All other dinucleotides: Penalized with mismatch penalty (linear or exponential).

    Parameters:
        seq (str): DNA sequence (A/C/G/T).
        GC_weight, AT_weight, GT_weight, AC_weight (float): Weights for dinucleotide types.
        consecutive_AT_scoring (tuple): Additional penalties for consecutive AT/TA.
        mismatch_penalty_type (str): "linear" or "exponential" penalty for mismatches.
        mismatch_penalty_starting_value (float): Starting penalty for mismatches.
        mismatch_penalty_linear_delta (float): Increment for linear penalty.
        cadence_reward (float): Additional score for canonical dinucleotides.

    Returns:
        scoring_array (np.ndarray): Array of scores for each position (dinucleotide).
    """
    scoring_array = np.empty(len(seq) - 1, dtype=float)  # Score for each dinucleotide
    mismatches_counter = 0       # Tracks consecutive mismatches
    consecutive_AT_counter = 0   # Tracks consecutive AT/TA dinucleotides

    for i in range(len(seq) - 1):
        t = seq[i:i+2].upper()   # Current dinucleotide
        if t in ("GC", "CG"):
            # Strongest Z-DNA-forming dinucleotides
            scoring_array[i] = GC_weight
            mismatches_counter = 0
            consecutive_AT_counter = 0
        elif t in ("GT", "TG"):
            scoring_array[i] = GT_weight
            mismatches_counter = 0
            consecutive_AT_counter = 0
        elif t in ("AC", "CA"):
            scoring_array[i] = AC_weight
            mismatches_counter = 0
            consecutive_AT_counter = 0
        elif t in ("AT", "TA"):
            # Apply additional penalty for consecutive AT/TA
            adjusted_weight = AT_weight
            if consecutive_AT_counter < len(consecutive_AT_scoring):
                adjusted_weight += consecutive_AT_scoring[consecutive_AT_counter]
            else:
                adjusted_weight += consecutive_AT_scoring[-1]  # Use last value if too many
            scoring_array[i] = adjusted_weight
            consecutive_AT_counter += 1
            mismatches_counter = 0
        else:
            # Penalty for non-canonical dinucleotides (mismatch)
            mismatches_counter += 1
            consecutive_AT_counter = 0
            if mismatch_penalty_type == "exponential":
                # Exponential penalty for long stretches
                scoring_array[i] = -mismatch_penalty_starting_value ** mismatches_counter \
                                   if mismatches_counter < 15 else -32000
            elif mismatch_penalty_type == "linear":
                # Linear penalty increases with length of mismatch stretch
                scoring_array[i] = -mismatch_penalty_starting_value \
                                   - mismatch_penalty_linear_delta * (mismatches_counter - 1)
            else:
                scoring_array[i] = -10  # Default penalty

        # Optional reward for canonical dinucleotide "cadence"
        if t in ("GC", "CG", "GT", "TG", "AC", "CA", "AT", "TA"):
            scoring_array[i] += cadence_reward

    return scoring_array

def wrap(seq, width=50):
    """
    Formats a DNA sequence into lines of specified width.

    Used for displaying long motifs in a readable way.

    Parameters:
        seq (str): The DNA sequence to format.
        width (int): The length of each line.

    Returns:
        str: Formatted sequence.
    """
    return '\n'.join(seq[i:i+width] for i in range(0, len(seq), width))

def find_zdna(
    seq,
    threshold=50,
    drop_threshold=50, #A region must have a total Z-DNA score of at least 50 to be reported as a motif.It’s a biological cutoff.
    GC_weight=7.0, AT_weight=0.5, GT_weight=1.25, AC_weight=1.25,
    consecutive_AT_scoring=(0.5, 0.5, 0.5, 0.5, 0.0, 0.0, -5.0, -100.0),
    mismatch_penalty_type="linear",
    mismatch_penalty_starting_value=3,
    mismatch_penalty_linear_delta=3,
    cadence_reward=0.0
):
    """
    Identifies candidate Z-DNA-forming regions ("motifs") in the input sequence.

    Uses a modified Kadane's algorithm to find subarrays (contiguous regions)
    where the cumulative Z-DNA score exceeds the specified threshold.
    Motif is reported when the score drops by at least drop_threshold, or at the end.

    Parameters:
        seq (str): DNA sequence.
        threshold (float): Minimum score for motif detection.
        drop_threshold (float): Minimum score drop to segment motifs.
        All other parameters: See zdna_seeker_scoring_array.

    Returns:
        motifs (list of dict): Each dict contains motif annotation:
            - Class: always "Z-DNA"
            - Subtype: always "Z-Seeker"
            - Start: 1-based motif start coordinate
            - End: 1-based motif end coordinate (inclusive)
            - Length: motif length
            - Sequence: wrapped motif sequence
            - ScoreMethod: "Z-Seeker Weighted"
            - Score: motif score (float, 2 decimals)
    """
    seq = seq.upper()  # Ensure uppercase
    scoring = zdna_seeker_scoring_array(
        seq,
        GC_weight=GC_weight, AT_weight=AT_weight,
        GT_weight=GT_weight, AC_weight=AC_weight,
        consecutive_AT_scoring=consecutive_AT_scoring,
        mismatch_penalty_type=mismatch_penalty_type,
        mismatch_penalty_starting_value=mismatch_penalty_starting_value,
        mismatch_penalty_linear_delta=mismatch_penalty_linear_delta,
        cadence_reward=cadence_reward
    )
    motifs = []

    # Find motifs using modified Kadane's algorithm
    start_idx = 0
    max_ending_here = scoring[0]
    current_max = 0
    candidate = None
    end_idx = 1

    for i in range(1, len(scoring)):
        num = scoring[i]
        # Start new region if current score is higher than continuing previous
        if num >= max_ending_here + num:
            start_idx = i
            end_idx = i + 1
            max_ending_here = num
        else:
            max_ending_here += num
            end_idx = i + 1

        # If score above threshold, mark as candidate
        if max_ending_here >= threshold and (candidate is None or current_max < max_ending_here):
            candidate = (start_idx, end_idx, max_ending_here)
            current_max = max_ending_here

        # If score drops enough, record motif and reset
        if candidate and (max_ending_here < 0 or current_max - max_ending_here >= drop_threshold):
            s, e, score = candidate
            motifs.append({
                "Class": "Z-DNA",
                "Subtype": "Z-Seeker",
                "Start": s + 1,           # 1-based indexing for output
                "End": e + 1,             # inclusive end (covers all input bases)
                "Length": e - s + 1,
                "Sequence": wrap(seq[s:e+1]),
                "ScoreMethod": "Z-Seeker Weighted",
                "Score": f"{score:.2f}",
            })
            candidate = None
            max_ending_here = current_max = 0

    # If there's a leftover candidate at the end, record it
    if candidate:
        s, e, score = candidate
        motifs.append({
            "Class": "Z-DNA",
            "Subtype": "Z-Seeker",
            "Start": s + 1,
            "End": e + 1,
            "Length": e - s + 1,
            "Sequence": wrap(seq[s:e+1]),
            "ScoreMethod": "Z-Seeker Weighted",
            "Score": f"{score:.2f}",
        })
    return motifs
###################################################################################
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

########################################### 3. Slipped-DNA  ############################################################

def wrap(seq, width=60):
    """Wraps sequence for display (optional, adjust/remove as needed)."""
    return '\n'.join(seq[i:i+width] for i in range(0, len(seq), width))

def find_slipped_dna(seq):
    """
    Finds slipped-strand DNA motifs: direct repeats and short tandem repeats (STRs/microsatellites).
    Returns a list of dictionaries with details for each motif.
    """
    results = []
    # Direct repeats (>=10bp unit, 2 copies, up to 300bp unit)
    min_len_dr = 10
    max_len_dr = 300
    for i in range(len(seq) - min_len_dr * 2 + 1):
        for l in range(min_len_dr, min(max_len_dr+1, (len(seq)-i)//2+1)):
            repeat = seq[i:i+l]
            if seq[i+l:i+2*l] == repeat:
                results.append({
                    "Class": "Slipped_DNA",
                    "Subtype": "Direct_Repeat",
                    "Start": i+1,
                    "End": i+2*l,
                    "Length": 2*l,
                    "Sequence": wrap(repeat+repeat),
                    "ScoreMethod": "nBST_DR",
                    "Score": f"{min(1.0, l/300):.2f}"
                })

    # Short Tandem Repeats (microsatellites: 1-6bp unit, >=5 copies, at least 15bp array)
    min_unit_str = 1
    max_unit_str = 6
    min_reps_str = 5
    min_len_str = 15
    for unit in range(min_unit_str, max_unit_str+1):
        for i in range(len(seq) - unit * min_reps_str + 1):
            repeat_unit = seq[i:i+unit]
            if 'n' in repeat_unit.lower():
                continue  # skip ambiguous bases
            reps = 1
            while (i + reps*unit + unit <= len(seq) and
                   seq[i + reps*unit:i + (reps+1)*unit] == repeat_unit):
                reps += 1
            if reps >= min_reps_str and reps*unit >= min_len_str:
                # Optionally check for partial repeat at end
                remainder = 0
                rs = i + reps*unit
                re = rs
                while (re < len(seq) and seq[re] == repeat_unit[re % unit]):
                    remainder += 1
                    re += 1
                results.append({
                    "Class": "Slipped_DNA",
                    "Subtype": "STR",
                    "Start": i+1,
                    "End": i + reps*unit + remainder,
                    "Length": reps*unit + remainder,
                    "Unit": repeat_unit,
                    "Copies": reps,
                    "Sequence": wrap(seq[i:i + reps*unit + remainder]),
                    "ScoreMethod": "nBST_STR",
                    "Score": f"{min(1.0, reps/20):.2f}"
                })
                # Skip overlapping STRs
                i = i + reps*unit + remainder - 1

    return results

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++#

############################################ 4.R-loop  ############################################################

import re

def wrap(seq, width=60):
    return '\n'.join(seq[i:i+width] for i in range(0, len(seq), width))

def gc_content(seq):
    gc = seq.count('G') + seq.count('C')
    return (gc / len(seq)) * 100 if seq else 0

# RLFS regex models based on QmRLFS-finder
RLFS_MODELS = {
    "m1": r"G{3,}[ATGC]{1,10}?G{3,}(?:[ATGC]{1,10}?G{3,}){1,}",
    "m2": r"G{4,}(?:[ATGC]{1,10}?G{4,}){1,}",
}

def find_rlfs(seq, models=("m1", "m2")):
    if len(seq) < 100:
        return []
    results = []
    for model_name in models:
        pattern = RLFS_MODELS[model_name]
        for m in re.finditer(pattern, seq, re.IGNORECASE):
            riz_seq = m.group(0)
            if gc_content(riz_seq) < 50:
                continue
            rez = find_rez_max(seq, m.end())
            if rez:
                rez_seq = rez['seq']
                stability = min(1.0, 0.6 * (gc_content(riz_seq + rez_seq) / 100) +
                                0.4 * (len(re.findall(r"G{3,}", riz_seq + rez_seq)) / 5))
                results.append({
                    "Class": "R-Loop",
                    "Subtype": f"RLFS_{model_name}",
                    "Start": m.start() + 1,
                    "End": m.start() + len(riz_seq) + rez['end'],
                    "Length": len(riz_seq) + rez['end'],
                    "Sequence": wrap(riz_seq + rez_seq),
                    "ScoreMethod": "QmRLFS_Thermo",
                    "Score": f"{stability:.2f}"
                })
    return results

def find_rez_max(seq, start_pos, max_len=2000, step=100, min_gc=40):
    # Find the maximal REZ region downstream of RIZ with >=min_gc%
    max_window = ""
    for win_start in range(start_pos, min(len(seq), start_pos + max_len), step):
        win_end = min(win_start + step, len(seq))
        window = seq[win_start:win_end]
        if gc_content(window) >= min_gc and len(window) > len(max_window):
            max_window = window
    if max_window:
        return {'seq': max_window, 'end': len(max_window)}
    return None
####################################++++++++++++++++++++++++++++++++##########################################

#################################5. Cruciform ###################################################################
def find_cruciform(seq):
    """Detect cruciform-forming palindromes with a spacer (0–3 bases) and arm (10–100 bases)."""
    results = []
    for i in range(len(seq) - 2*10):  # minimum required for two arms
        for arm_len in range(10, 101):  # arm length from 10 to 100
            for spacer_len in range(0, 4):  # spacer length from 0 to 3
                arm = seq[i:i+arm_len]
                rev_arm = reverse_complement(arm)
                mid = i + arm_len + spacer_len
                if mid + arm_len > len(seq): continue
                candidate = seq[mid:mid+arm_len]
                if candidate == rev_arm:
                    full = seq[i:mid+arm_len]
                    score = min(1.0, (arm_len / 100) + ((arm.count('A') + arm.count('T')) / arm_len * 0.3))
                    results.append({
                        "Class": "Cruciform",
                        "Subtype": f"Inverted_Repeat_spacer{spacer_len}",
                        "Start": i+1,
                        "End": mid+arm_len,
                        "Length": len(full),
                        "Sequence": wrap(full),
                        "ScoreMethod": "nBST_IR",
                        "Score": f"{score:.2f}"
                    })
    return results

####++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##################################
################################################6. Triplex DNA #######################################################
import re

def purine_fraction(seq):
    """Return fraction of purines (A/G) in the sequence."""
    return (seq.count('A') + seq.count('G')) / max(1, len(seq))

def pyrimidine_fraction(seq):
    """Return fraction of pyrimidines (C/T) in the sequence."""
    return (seq.count('C') + seq.count('T')) / max(1, len(seq))

def wrap(seq, width=60):
    """Wrap sequence output to fixed width (for display)."""
    return '\n'.join(seq[i:i+width] for i in range(0, len(seq), width))

def find_hdna(seq):
    """Detects mirror repeats suitable for triplex (H-DNA) formation with spacer 0–8 nt."""
    results = []
    n = len(seq)
    # Repeat length: 10–100nt, Spacer: 0–8nt
    for rep_len in range(10, min(101, n//2)):
        for spacer in range(0, 9):  # Spacer 0–8 nt inclusive
            # Pattern: (repeat)[spacer](repeat)
            pattern = re.compile(rf"(?=(([ATGC]{{{rep_len}}})[ATGC]{{{spacer}}}\2))", re.IGNORECASE)
            for m in pattern.finditer(seq):
                repeat = m.group(2)
                mirror_start = m.start()
                mirror_end = mirror_start + 2*rep_len + spacer
                if mirror_end > n:
                    continue
                full_seq = seq[mirror_start:mirror_end]
                pur_frac = purine_fraction(full_seq)
                pyr_frac = pyrimidine_fraction(full_seq)
                is_triplex = (pur_frac >= 0.9 or pyr_frac >= 0.9)
                results.append({
                    "Class": "Triplex_DNA" if is_triplex else "Mirror_Repeat",
                    "Subtype": "Triplex_Motif" if is_triplex else "Mirror_Repeat",
                    "Start": mirror_start + 1,
                    "End": mirror_end,
                    "Length": len(full_seq),
                    "Spacer": spacer,
                    "Sequence": wrap(full_seq),
                    "PurineFrac": round(pur_frac, 2),
                    "PyrimidineFrac": round(pyr_frac, 2)
                })
    return results


####++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##################################
################################################7. Sticky DNA #######################################################
import re

def find_sticky_dna(seq):
    """
    Detects uninterrupted (GAA)n or (TTC)n tracts where n = 59–270 (potential sticky DNA partners).
    Sticky DNA forms only when two such tracts are present (on separate molecules).
    References:
      - Sakamoto N et al., Molecular Cell, 1999.
      - Potaman VN et al., Nucleic Acids Res, 2004.
      - Grabczyk E et al., Biochemistry, 2000.
    """
    motifs = []
    # Simple regex: matches uninterrupted GAA or TTC repeats, n = 59–270
    pattern = r"(?:GAA){59,270}|(?:TTC){59,270}"
    for m in re.finditer(pattern, seq):
        repeat_count = len(m.group()) // 3
        motifs.append({
            "Class": "Sticky_DNA",
            "Subtype": "GAA_TTC_Repeat",
            "Start": m.start() + 1,
            "End": m.end(),
            "Length": len(m.group()),
            "RepeatCount": repeat_count,
            "Sequence": m.group(),  # or wrap(m.group()) if you want wrapped output
            "ScoreMethod": "Sakamoto1999",
            "Score": f"{min(1.0, repeat_count/270):.2f}",
            "References": (
                "Sakamoto N et al., Molecular Cell, 1999; "
                "Potaman VN et al., Nucleic Acids Res, 2004; "
                "Grabczyk E et al., Biochemistry, 2000"
            )
        })
    return motifs


#######################################++++++++++##########################################################
####++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++##################################
################################################8. G-Triplex DNA #######################################################

def find_gtriplex(seq):
    results = []
    for m in overlapping_finditer(r"(?=(G{3,}\w{1,7}G{3,}\w{1,7}G{3,}))", seq):
        seq_frag = m.group(1)
        g_runs = [len(r) for r in re.findall(r"G{3,}", seq_frag)]
        if len(g_runs) < 3: continue
        loops = [len(l) for l in re.findall(r"G{3,}(\w{1,7})G{3,}", seq_frag)]
        score = min(1.0, sum(g_runs)/15 + sum(1/l if l > 0 else 0.5 for l in loops)/3)
        results.append({
            "Class": "G-Triplex", "Subtype": "Three_G-Runs", "Start": m.start()+1,
            "End": m.start()+len(seq_frag), "Length": len(seq_frag), "Sequence": wrap(seq_frag),
            "ScoreMethod": "G3_Stability", "Score": f"{score:.2f}"})
    return results

def find_gquadruplex(seq):
    return [{
        "Class": "G4", "Subtype": "Canonical_G4", "Start": m.start()+1, "End": m.end(),
        "Length": len(m.group(1)), "Sequence": wrap(m.group(1)), "ScoreMethod": "G4Hunter_v2",
        "Score": f"{g4hunter_score(m.group(1)):.2f}"}
        for m in overlapping_finditer(r"(G{3,}\w{1,7}G{3,}\w{1,7}G{3,}\w{1,7}G{3,})", seq)
        if g4hunter_score(m.group(1)) >= 1.2]

def find_relaxed_gquadruplex(seq):
    return [{
        "Class": "G4", "Subtype": "Relaxed_G4", "Start": m.start()+1, "End": m.end(),
        "Length": len(m.group(1)), "Sequence": wrap(m.group(1)), "ScoreMethod": "G4Hunter_LongLoop",
        "Score": f"{g4hunter_score(m.group(1))*0.8:.2f}"}
        for m in overlapping_finditer(r"(G{3,}\w{1,12}G{3,}\w{1,12}G{3,}\w{1,12}G{3,})", seq)
        if g4hunter_score(m.group(1)) >= 0.8]

def find_bulged_gquadruplex(seq):
    return [{
        "Class": "G4", "Subtype": "Bulged_G4", "Start": m.start()+1, "End": m.end(),
        "Length": len(m.group(1)), "Sequence": wrap(m.group(1)), "ScoreMethod": "G4Hunter_Bulge",
        "Score": f"{g4hunter_score(m.group(1))*0.7:.2f}"}
        for m in overlapping_finditer(r"(G{3,}\w{0,3}G{3,}\w{0,3}G{3,}\w{0,3}G{3,})", seq)
        if len(re.findall(r"G{3,}", m.group(1))) >= 4]

def find_bipartite_gquadruplex(seq):
    results = []
    for m in overlapping_finditer(r"(G{3,}\w{1,30}G{3,}\w{1,30}G{3,}\w{1,30}G{3,}\w{1,30}G{3,}\w{1,30}G{3,})", seq):
        seq_frag = m.group(1)
        if len(re.findall(r"G{3,}", seq_frag)) < 6: continue
        unit1, unit2 = seq_frag[:len(seq_frag)//2], seq_frag[len(seq_frag)//2:]
        score = max(g4hunter_score(unit1), g4hunter_score(unit2)) * 0.9
        if score >= 0.9:
            results.append({
                "Class": "G4", "Subtype": "Bipartite_G4", "Start": m.start()+1, "End": m.end(),
                "Length": len(seq_frag), "Sequence": wrap(seq_frag), "ScoreMethod": "Bipartite_Score",
                "Score": f"{score:.2f}"})
    return results

def find_multimeric_gquadruplex(seq):
    return [{
        "Class": "G4", "Subtype": "Multimeric_G4", "Start": m.start()+1, "End": m.end(),
        "Length": len(m.group(1)), "Sequence": wrap(m.group(1)), "ScoreMethod": "G4Hunter_Multimer",
        "Score": f"{g4hunter_score(m.group(1))*1.2:.2f}"}
        for m in overlapping_finditer(r"(G{3,}\w{1,12}){4,}", seq)
        if g4hunter_score(m.group(1)) >= 1.0]

def find_imotif(seq):
    return [{
        "Class": "i-Motif", "Subtype": "C_Quadruplex", "Start": m.start()+1,
        "End": m.start()+len(m.group(1)), "Length": len(m.group(1)), "Sequence": wrap(m.group(1)),
        "ScoreMethod": "iM_Stability", "Score": f"{imotif_score(m.group(1)):.2f}"}
        for m in overlapping_finditer(r"(?=(C{3,}\w{0,7}C{3,}\w{0,7}C{3,}\w{0,7}C{3,}))", seq)
        if imotif_score(m.group(1)) >= 0.7]

def imotif_score(seq):
    c_runs = [len(r) for r in re.findall(r"C{3,}", seq)]
    return min(1.0, sum(c_runs)/16 + (seq.count('C')/len(seq))*0.5) if len(c_runs) >= 4 else 0

def find_hybrids(seq):
    g4_regions = [(m.start()+1, m.end(), m.group(1)) for m in overlapping_finditer(r"(G{3,}\w{1,12}G{3,}\w{1,12}G{3,}\w{1,12}G{3,})", seq)]
    im_regions = [(m.start()+1, m.start()+len(m.group(1)), m.group(1)) for m in overlapping_finditer(r"(C{3,}\w{0,7}C{3,}\w{0,7}C{3,}\w{0,7}C{3,})", seq)]
    return [{
        "Class": "Hybrid", "Subtype": "G4_iM_Overlap", "Start": min(g_start, c_start),
        "End": max(g_end, c_end), "Length": max(g_end, c_end)-min(g_start, c_start)+1,
        "Sequence": wrap(seq[min(g_start,c_start)-1:max(g_end,c_end)]),
        "ScoreMethod": "Hybrid_Score",
        "Score": f"{min(1.0, (g4hunter_score(g_seq)+imotif_score(c_seq))/2*1.1):.2f}"}
        for g_start, g_end, g_seq in g4_regions
        for c_start, c_end, c_seq in im_regions
        if (g_start <= c_end) and (c_start <= g_end)]

def find_hotspots(motif_hits, seq_len, window=100, min_count=3):
    hotspots = []
    positions = [(hit['Start'], hit['End']) for hit in motif_hits]
    for i in range(0, seq_len - window + 1):
        region_start, region_end = i+1, i+window
        count = sum(s <= region_end and e >= region_start for s,e in positions)
        if count >= min_count:
            motifs_in_region = [m for m in motif_hits if m['Start'] <= region_end and m['End'] >= region_start]
            type_div = len({m['Subtype'] for m in motifs_in_region})
            hotspots.append({
                "RegionStart": region_start, "RegionEnd": region_end,
                "MotifCount": count, "TypeDiversity": type_div,
                "Score": f"{min(1.0, count/10 + type_div/5):.2f}"})
    return merge_hotspots(hotspots)

def merge_hotspots(hotspots):
    if not hotspots: return []
    merged = [hotspots[0]]
    for current in hotspots[1:]:
        last = merged[-1]
        if current['RegionStart'] <= last['RegionEnd']:
            last['RegionEnd'] = max(last['RegionEnd'], current['RegionEnd'])
            last['MotifCount'] += current['MotifCount']
            last['TypeDiversity'] = max(last['TypeDiversity'], current['TypeDiversity'])
            last['Score'] = f"{min(1.0, float(last['Score']) + float(current['Score'])):.2f}"
        else:
            merged.append(current)
    return merged

# ========== NEW: NON-OVERLAPPING MOTIF SELECTION ==========

def select_best_nonoverlapping_motifs(motifs: List[Dict], motif_priority: List[str] = None) -> List[Dict]:
    """
    Given a list of motif dicts, return a non-overlapping list, keeping the best motif in overlapping regions.
    Priority: motif type > best score > longest motif.
    motif_priority: list from highest to lowest priority, using motif 'Subtype' field.
    """
    if motif_priority is None:
        motif_priority = [
            'Multimeric_G4', 'Bipartite_G4', 'Dimeric_G4', 'Canonical_G4', 
            'Relaxed_G4', 'Non_canonical_G4', 'Bulged_G4', 'Three_G-Runs'
        ]
    # Map subtype to rank (lower is higher priority)
    subtype_rank = {subtype: i for i, subtype in enumerate(motif_priority)}
    def motif_key(m):
        # Lower rank is better; if subtype not in priority, assign lowest priority
        rank = subtype_rank.get(m.get('Subtype'), len(subtype_rank))
        try:
            score = float(m.get('Score', 0))
        except ValueError:
            score = 0.0
        length = m.get('Length', 0)
        return (rank, -score, -length)  # Lower rank, higher score, longer

    # Sort all motifs by best-to-worst (for tie-breaking)
    sorted_motifs = sorted(motifs, key=motif_key)
    selected = []
    occupied = set()
    for m in sorted_motifs:
        region = set(range(m['Start'], m['End']+1))
        if occupied.isdisjoint(region):
            selected.append(m)
            occupied.update(region)
    return selected

# ========== END NEW SECTION ==========

def validate_motif(motif, seq_length):
    """
    Checks that a motif dict has required keys and valid coordinates.
    Returns True if the motif is valid, False otherwise.
    """
    required_keys = ["Class", "Subtype", "Start", "End", "Length", "Sequence"]
    # Check all required fields exist
    if not all(key in motif for key in required_keys):
        return False
    # Check valid coordinates
    if not (1 <= motif["Start"] <= motif["End"] <= seq_length):
        return False
    # Check non-empty sequence
    if len(motif["Sequence"].replace('\n', '')) == 0:
        return False
    return True
