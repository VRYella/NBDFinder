import re
import numpy as np
from typing import List, Dict, Tuple
from collections import defaultdict, Counter
import random

# ========= UTILITY FUNCTIONS =========

def parse_fasta(fasta_str: str) -> str:
    return "".join([line.strip() for line in fasta_str.split('\n') if not line.startswith(">")]).upper().replace(" ", "").replace("U", "T")

def wrap(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))

def gc_content(seq: str) -> float:
    gc = seq.count('G') + seq.count('C')
    return (gc / max(1, len(seq))) * 100 if seq else 0

def reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ATGC", "TACG"))[::-1]

def is_palindrome(seq: str) -> bool:
    return seq == reverse_complement(seq)

def percentileofscore(a, score, kind='rank'):
    a = np.asarray(a)
    if len(a) == 0: return 0.0
    if kind == 'rank': return (sum(a <= score) / len(a) * 100)
    elif kind == 'strict': return (sum(a < score) / len(a)) * 100
    elif kind == 'weak': return (sum(a <= score) / len(a)) * 100
    elif kind == 'mean': return (sum(a < score) + sum(a <= score)) / (2 * len(a)) * 100
    else: raise ValueError("kind must be 'rank', 'strict', 'weak' or 'mean'")

def overlapping_finditer(pattern, seq):
    regex = re.compile(pattern, re.IGNORECASE)
    pos = 0
    while pos < len(seq):
        m = regex.search(seq, pos)
        if not m:
            break
        yield m
        pos = m.start() + 1

# ========= K-MER CONSERVATION SCORING =========

def shuffle_sequence(seq: str, num_shuffles: int = 100) -> List[str]:
    """Generate shuffled versions of the sequence for conservation analysis."""
    shuffled_seqs = []
    seq_list = list(seq)
    for _ in range(num_shuffles):
        shuffled = seq_list.copy()
        random.shuffle(shuffled)
        shuffled_seqs.append(''.join(shuffled))
    return shuffled_seqs

def kmer_conservation(seq: str, k: int = 3, num_shuffles: int = 100) -> float:
    """Calculate k-mer conservation score by comparing to shuffled sequences."""
    if len(seq) < k:
        return 0.0
    
    # Count k-mers in original sequence
    original_kmers = defaultdict(int)
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        original_kmers[kmer] += 1
    
    # Calculate diversity in original (Shannon entropy-like measure)
    total_kmers = sum(original_kmers.values())
    if total_kmers == 0:
        return 0.0
    
    original_diversity = -sum((count/total_kmers) * np.log2(count/total_kmers) 
                             for count in original_kmers.values() if count > 0)
    
    # Calculate average diversity in shuffled sequences
    shuffled_diversities = []
    shuffled_seqs = shuffle_sequence(seq, num_shuffles)
    
    for shuffled_seq in shuffled_seqs:
        shuffled_kmers = defaultdict(int)
        for i in range(len(shuffled_seq) - k + 1):
            kmer = shuffled_seq[i:i+k]
            shuffled_kmers[kmer] += 1
        
        shuffled_total = sum(shuffled_kmers.values())
        if shuffled_total > 0:
            shuffled_diversity = -sum((count/shuffled_total) * np.log2(count/shuffled_total) 
                                    for count in shuffled_kmers.values() if count > 0)
            shuffled_diversities.append(shuffled_diversity)
    
    if not shuffled_diversities:
        return 0.0
    
    avg_shuffled_diversity = np.mean(shuffled_diversities)
    
    # Conservation score: lower diversity relative to shuffled = higher conservation
    if avg_shuffled_diversity == 0:
        return 1.0 if original_diversity == 0 else 0.0
    
    conservation_score = max(0, (avg_shuffled_diversity - original_diversity) / avg_shuffled_diversity)
    return min(1.0, conservation_score)

def classify_conservation(score: float) -> str:
    """Classify conservation score into categories."""
    if score >= 0.8:
        return "High"
    elif score >= 0.6:
        return "Medium"
    elif score >= 0.4:
        return "Low"
    else:
        return "None"

def add_conservation_metrics(motif: Dict, seq: str) -> Dict:
    """Add k-mer conservation metrics to a motif dictionary."""
    motif_seq = seq[motif['Start']-1:motif['End']]
    
    # Calculate conservation for different k-mer sizes
    k3_score = kmer_conservation(motif_seq, k=3)
    k4_score = kmer_conservation(motif_seq, k=4)
    k5_score = kmer_conservation(motif_seq, k=5)
    
    # Average conservation score
    avg_conservation = (k3_score + k4_score + k5_score) / 3
    
    motif['Conservation_Score'] = f"{avg_conservation:.3f}"
    motif['Conservation_Class'] = classify_conservation(avg_conservation)
    motif['K3_Conservation'] = f"{k3_score:.3f}"
    motif['K4_Conservation'] = f"{k4_score:.3f}"
    motif['K5_Conservation'] = f"{k5_score:.3f}"
    
    return motif

def validate_motif_counts(motifs: List[Dict], seq: str) -> bool:
    """Validate that motif positions are within sequence bounds and non-overlapping where expected."""
    seq_len = len(seq)
    
    for motif in motifs:
        # Check bounds
        if motif['Start'] < 1 or motif['End'] > seq_len:
            return False
        
        # Check that start <= end
        if motif['Start'] > motif['End']:
            return False
        
        # Check that length matches
        expected_length = motif['End'] - motif['Start'] + 1
        if motif['Length'] != expected_length:
            return False
    
    return True

# ========== 1. CURVED DNA (Strict PolyA/PolyT Only) ==========

def find_polyA_polyT_tracts(seq: str, min_len: int = 7) -> list:
    results = []
    i = 0
    n = len(seq)
    while i < n:
        if seq[i] == 'A':
            start = i
            while i < n and seq[i] == 'A':
                i += 1
            if i - start >= min_len:
                results.append((start, i-1, seq[start:i]))
        elif seq[i] == 'T':
            start = i
            while i < n and seq[i] == 'T':
                i += 1
            if i - start >= min_len:
                results.append((start, i-1, seq[start:i]))
        else:
            i += 1
    return results

def curvature_score(seq):
    return len(seq)

def find_global_curved_polyA_polyT(seq: str, min_tract_len: int = 3, min_repeats: int = 3, min_spacing: int = 8, max_spacing: int = 12, min_score: int = 6) -> Tuple[list, list]:
    tracts = find_polyA_polyT_tracts(seq, min_tract_len)
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
                    "Class": "Curved DNA",
                    "Subtype": "Global Curved PolyA/PolyT",
                    "Start": group[0][0] + 1,
                    "End": group[-1][1] + 1,
                    "Length": group[-1][1] - group[0][0] + 1,
                    "Sequence": wrap(motif_seq),
                    "ScoreMethod": "PolyA/PolyT curvature tract score",
                    "Score": score
                }
                results.append(motif)
                apr_regions.append((motif["Start"], motif["End"]))
    return results, apr_regions

def find_local_curved_polyA_polyT(seq: str, apr_regions: list, min_len: int = 7) -> list:
    results = []
    tracts = find_polyA_polyT_tracts(seq, min_len)
    for start, end, tract_seq in tracts:
        s, e = start + 1, end + 1
        if not any(r_start <= s <= r_end or r_start <= e <= r_end for r_start, r_end in apr_regions):
            results.append({
                "Class": "Curved DNA",
                "Subtype": "Local Curved PolyA/PolyT",
                "Start": s,
                "End": e,
                "Length": len(tract_seq),
                "Sequence": wrap(tract_seq),
                "ScoreMethod": "PolyA/PolyT tract length",
                "Score": len(tract_seq)
            })
    return results

def find_curved_DNA(seq: str) -> list:
    global_results, apr_regions = find_global_curved_polyA_polyT(seq)
    local_results = find_local_curved_polyA_polyT(seq, apr_regions)
    return global_results + local_results

# ========== 2. Z-DNA ==========

def zdna_seeker_scoring_array(
    seq,
    GC_weight=7.0, AT_weight=0.5, GT_weight=1.25, AC_weight=1.25,
    consecutive_AT_scoring=(0.5, 0.5, 0.5, 0.5, 0.0, 0.0, -5.0, -100.0),
    mismatch_penalty_type="linear",
    mismatch_penalty_starting_value=3,
    mismatch_penalty_linear_delta=3,
    cadence_reward=0.0
):
    scoring_array = np.empty(len(seq) - 1, dtype=float)
    mismatches_counter = 0
    consecutive_AT_counter = 0
    for i in range(len(seq) - 1):
        t = seq[i:i+2].upper()
        if t in ("GC", "CG"):
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
            adjusted_weight = AT_weight
            if consecutive_AT_counter < len(consecutive_AT_scoring):
                adjusted_weight += consecutive_AT_scoring[consecutive_AT_counter]
            else:
                adjusted_weight += consecutive_AT_scoring[-1]
            scoring_array[i] = adjusted_weight
            consecutive_AT_counter += 1
            mismatches_counter = 0
        else:
            mismatches_counter += 1
            consecutive_AT_counter = 0
            if mismatch_penalty_type == "exponential":
                scoring_array[i] = -mismatch_penalty_starting_value ** mismatches_counter if mismatches_counter < 15 else -32000
            elif mismatch_penalty_type == "linear":
                scoring_array[i] = -mismatch_penalty_starting_value - mismatch_penalty_linear_delta * (mismatches_counter - 1)
            else:
                scoring_array[i] = -10
        if t in ("GC", "CG", "GT", "TG", "AC", "CA", "AT", "TA"):
            scoring_array[i] += cadence_reward
    return scoring_array

def find_zdna(
    seq,
    threshold=50,
    drop_threshold=50,
    GC_weight=7.0, AT_weight=0.5, GT_weight=1.25, AC_weight=1.25,
    consecutive_AT_scoring=(0.5, 0.5, 0.5, 0.5, 0.0, 0.0, -5.0, -100.0),
    mismatch_penalty_type="linear",
    mismatch_penalty_starting_value=3,
    mismatch_penalty_linear_delta=3,
    cadence_reward=0.0
):
    seq = seq.upper()
    if len(seq) < 12:
        return []
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
    start_idx = 0
    max_ending_here = scoring[0]
    current_max = 0
    candidate = None
    end_idx = 1
    for i in range(1, len(scoring)):
        num = scoring[i]
        if num >= max_ending_here + num:
            start_idx = i
            end_idx = i + 1
            max_ending_here = num
        else:
            max_ending_here += num
            end_idx = i + 1
        if max_ending_here >= threshold and (candidate is None or current_max < max_ending_here):
            candidate = (start_idx, end_idx, max_ending_here)
            current_max = max_ending_here
        if candidate and (max_ending_here < 0 or current_max - max_ending_here >= drop_threshold):
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
            candidate = None
            max_ending_here = current_max = 0
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

# ========== 2a. eGZ-MOTIF ==========

def find_egz_motif(seq):
    pattern = re.compile(r'(CGG){4,}', re.IGNORECASE)
    results = []
    for m in pattern.finditer(seq):
        motif_seq = m.group(0)
        n_repeats = len(motif_seq) // 3
        score = min(1.0, n_repeats/12)
        results.append({
            "Class": "eGZ (Extruded-G)",
            "Subtype": "Repeat_normalized",
            "Start": m.start() + 1,
            "End": m.end(),
            "Length": len(motif_seq),
            "Sequence": wrap(motif_seq),
            "ScoreMethod": "Repeat normalized",
            "Score": f"{score:.2f}",
            "CGG Repeats": n_repeats
        })
    return results

# ========== 3. SLIPPED DNA (DR, STR) ==========

def find_slipped_dna(seq):
    results = []
    min_len_dr = 10
    max_len_dr = 300
    for i in range(len(seq) - min_len_dr * 2 + 1):
        for l in range(min_len_dr, min(max_len_dr+1, (len(seq)-i)//2+1)):
            repeat = seq[i:i+l]
            if seq[i+l:i+2*l] == repeat:
                results.append({
                    "Class": "Slipped DNA",
                    "Subtype": "Direct Repeat",
                    "Start": i+1,
                    "End": i+2*l,
                    "Length": 2*l,
                    "Sequence": wrap(repeat+repeat),
                    "ScoreMethod": "nBST DR",
                    "Score": f"{min(1.0, l/300):.2f}"
                })
    min_unit_str = 1
    max_unit_str = 6
    min_reps_str = 5
    min_len_str = 15
    i = 0
    while i < len(seq) - min_unit_str * min_reps_str + 1:
        found = False
        for unit in range(min_unit_str, max_unit_str+1):
            if i + unit * min_reps_str > len(seq):
                continue
            repeat_unit = seq[i:i+unit]
            if 'n' in repeat_unit.lower():
                continue
            reps = 1
            while (i + reps*unit + unit <= len(seq) and
                   seq[i + reps*unit:i + (reps+1)*unit] == repeat_unit):
                reps += 1
            if reps >= min_reps_str and reps*unit >= min_len_str:
                remainder = 0
                rs = i + reps*unit
                re_idx = rs
                while (re_idx < len(seq) and seq[re_idx] == repeat_unit[re_idx % unit]):
                    remainder += 1
                    re_idx += 1
                results.append({
                    "Class": "Slipped DNA",
                    "Subtype": "Short Tandem Repeat",
                    "Start": i+1,
                    "End": i + reps*unit + remainder,
                    "Length": reps*unit + remainder,
                    "Sequence": wrap(seq[i:i + reps*unit + remainder]),
                    "ScoreMethod": "nBST STR",
                    "Score": f"{min(1.0, reps/20):.2f}"
                })
                i = i + reps*unit + remainder - 1
                found = True
                break
        if not found:
            i += 1
    return results

# ========== 4. R-LOOP ==========

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
                    "Subtype": f"RLFS {model_name}",
                    "Start": m.start() + 1,
                    "End": m.start() + len(riz_seq) + rez['end'],
                    "Length": len(riz_seq) + rez['end'],
                    "Sequence": wrap(riz_seq + rez_seq),
                    "ScoreMethod": "QmRLFS Thermo",
                    "Score": f"{stability:.2f}"
                })
    return results

def find_rez_max(seq, start_pos, max_len=2000, step=100, min_gc=40):
    max_window = ""
    for win_start in range(start_pos, min(len(seq), start_pos + max_len), step):
        win_end = min(win_start + step, len(seq))
        window = seq[win_start:win_end]
        if gc_content(window) >= min_gc and len(window) > len(max_window):
            max_window = window
    if max_window:
        return {'seq': max_window, 'end': len(max_window)}
    return None

# ========== 5. CRUCIFORM (INVERTED REPEAT) ==========

def find_cruciform(seq):
    results = []
    for i in range(len(seq) - 2*10):
        for arm_len in range(10, min(101, (len(seq)-i)//2)):
            for spacer_len in range(0, 4):
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
                        "Subtype": f"Inverted Repeat Spacer {spacer_len}",
                        "Start": i+1,
                        "End": mid+arm_len,
                        "Length": len(full),
                        "Sequence": wrap(full),
                        "ScoreMethod": "nBST IR",
                        "Score": f"{score:.2f}"
                    })
    return results

# ========== 6. TRIPLEX DNA (H-DNA, MIRROR REPEAT) ==========

def purine_fraction(seq):
    return (seq.count('A') + seq.count('G')) / max(1, len(seq))

def pyrimidine_fraction(seq):
    return (seq.count('C') + seq.count('T')) / max(1, len(seq))

def find_hdna(seq):
    results = []
    n = len(seq)
    for rep_len in range(10, min(101, n//2)):
        for spacer in range(0, 9):
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
                    "Class": "Triplex DNA" if is_triplex else "Mirror Repeat",
                    "Subtype": "Triplex Motif" if is_triplex else "Mirror Repeat",
                    "Start": mirror_start + 1,
                    "End": mirror_end,
                    "Length": len(full_seq),
                    "Sequence": wrap(full_seq),
                    "ScoreMethod": "Purine/Pyrimidine Fraction",
                    "Score": f"{max(pur_frac, pyr_frac):.2f}"
                })
    return results

# ========== 7. STICKY DNA ==========

def find_sticky_dna(seq):
    motifs = []
    seq = seq.replace('\n','').replace(' ','').upper()
    pattern = r"(?:GAA){59,}|(?:TTC){59,}"
    for m in re.finditer(pattern, seq):
        repeat_count = len(m.group()) // 3
        motifs.append({
            "Class": "Sticky DNA",
            "Subtype": "GAA/TTC Repeat",
            "Start": m.start() + 1,
            "End": m.end(),
            "Length": len(m.group()),
            "Repeat Count": repeat_count,
            "Sequence": m.group(),
            "ScoreMethod": "Sakamoto 1999",
            "Score": f"{min(1.0, repeat_count/270):.2f}",
        })
    return motifs

# ========== 8. G-TRIPLEX DNA & G4 VARIANTS ==========

def g4hunter_score(seq):
    scores = []
    for c in seq.upper():
        if c == 'G':
            scores.append(1)
        elif c == 'C':
            scores.append(-1)
        else:
            scores.append(0)
    return np.mean(scores) if scores else 0

def find_multimeric_gquadruplex(seq):
    results = []
    pattern = r"(G{3,}\w{1,12}){4,}"
    for m in overlapping_finditer(pattern, seq):
        motif_seq = m.group(0)
        if g4hunter_score(motif_seq) >= 1.0:
            results.append({
                "Class": "Multimeric G4",
                "Subtype": "Multimeric G4",
                "Start": m.start()+1,
                "End": m.end(),
                "Length": len(motif_seq),
                "Sequence": wrap(motif_seq),
                "ScoreMethod": "G4Hunter Multimer",
                "Score": f"{g4hunter_score(motif_seq)*1.2:.2f}"})
    return results

def find_bipartite_gquadruplex(seq):
    results = []
    pattern = r"(G{3,}\w{1,30}G{3,}\w{1,30}G{3,}\w{1,30}G{3,}\w{10,30}G{3,}\w{1,30}G{3,}\w{1,30}G{3,}\w{1,30}G{3,})"
    for m in overlapping_finditer(pattern, seq):
        motif_seq = m.group(1)
        if len(re.findall(r"G{3,}", motif_seq)) < 8:
            continue
        half = len(motif_seq)//2
        unit1, unit2 = motif_seq[:half], motif_seq[half:]
        score = max(g4hunter_score(unit1), g4hunter_score(unit2)) * 0.9
        if score >= 0.9:
            results.append({
                "Class": "Bipartite G4",
                "Subtype": "Bipartite G4",
                "Start": m.start()+1,
                "End": m.end(),
                "Length": len(motif_seq),
                "Sequence": wrap(motif_seq),
                "ScoreMethod": "Bipartite Score",
                "Score": f"{score:.2f}"})
    return results

def find_gquadruplex(seq):
    pattern = r"(G{3,}\w{1,7}G{3,}\w{1,7}G{3,}\w{1,7}G{3,})"
    results = []
    for m in overlapping_finditer(pattern, seq):
        motif_seq = m.group(1)
        score = g4hunter_score(motif_seq)
        if score >= 1.2:
            results.append({
                "Class": "G4",
                "Subtype": "Canonical G4",
                "Start": m.start()+1,
                "End": m.end(),
                "Length": len(motif_seq),
                "Sequence": wrap(motif_seq),
                "ScoreMethod": "G4Hunter v2",
                "Score": f"{score:.2f}"})
    return results

def find_relaxed_gquadruplex(seq):
    pattern = r"(G{3,}\w{8,12}G{3,}\w{8,12}G{3,}\w{8,12}G{3,})"
    results = []
    for m in overlapping_finditer(pattern, seq):
        motif_seq = m.group(1)
        score = g4hunter_score(motif_seq)
        if score >= 0.8:
            results.append({
                "Class": "Relaxed G4",
                "Subtype": "Relaxed G4",
                "Start": m.start()+1,
                "End": m.end(),
                "Length": len(motif_seq),
                "Sequence": wrap(motif_seq),
                "ScoreMethod": "G4Hunter LongLoop",
                "Score": f"{score*0.8:.2f}"})
    return results

def find_bulged_gquadruplex(seq):
    pattern = r"(G{3,}\w{0,3}G{3,}\w{0,3}G{3,}\w{0,3}G{3,})"
    results = []
    for m in overlapping_finditer(pattern, seq):
        motif_seq = m.group(1)
        if len(re.findall(r"G{3,}", motif_seq)) >= 4:
            score = g4hunter_score(motif_seq)
            results.append({
                "Class": "Bulged G4",
                "Subtype": "Bulged G4",
                "Start": m.start()+1,
                "End": m.end(),
                "Length": len(motif_seq),
                "Sequence": wrap(motif_seq),
                "ScoreMethod": "G4Hunter Bulge",
                "Score": f"{score*0.7:.2f}"})
    return results

def find_imperfect_gquadruplex(seq):
    pattern = r"(G{2,3}\w{1,7}G{3,}\w{1,7}G{3,}\w{1,7}G{3,})"\
              r"|(G{3,}\w{1,7}G{2,3}\w{1,7}G{3,}\w{1,7}G{3,})"\
              r"|(G{3,}\w{1,7}G{3,}\w{1,7}G{2,3}\w{1,7}G{3,})"\
              r"|(G{3,}\w{1,7}G{3,}\w{1,7}G{3,}\w{1,7}G{2,3})"
    results = []
    for m in overlapping_finditer(pattern, seq):
        motif_seq = m.group(0)
        score = g4hunter_score(motif_seq)
        if score >= 1.2:
            results.append({
                "Class": "G4",
                "Subtype": "Imperfect G4",
                "Start": m.start()+1,
                "End": m.end(),
                "Length": len(motif_seq),
                "Sequence": wrap(motif_seq),
                "ScoreMethod": "G4Hunter Imperfect",
                "Score": f"{score:.2f}"})
    return results

def find_gtriplex(seq):
    pattern = r"(G{3,}\w{1,7}G{3,}\w{1,7}G{3,})"
    results = []
    for m in overlapping_finditer(pattern, seq):
        motif_seq = m.group(1)
        g_runs = [len(r) for r in re.findall(r"G{3,}", motif_seq)]
        if len(g_runs) < 3: continue
        loops = [len(l) for l in re.findall(r"G{3,}(\w{1,7})G{3,}", motif_seq)]
        score = min(1.0, sum(g_runs)/15 + sum(1/l if l > 0 else 0.5 for l in loops)/3)
        results.append({
            "Class": "G-Triplex",
            "Subtype": "Three G-Runs",
            "Start": m.start()+1,
            "End": m.end(),
            "Length": len(motif_seq),
            "Sequence": wrap(motif_seq),
            "ScoreMethod": "G3 Stability",
            "Score": f"{score:.2f}"})
    return results

def find_nonoverlapping_g4_variants(seq):
    mask = [False]*len(seq)
    results = []
    finders = [
        find_multimeric_gquadruplex,
        find_bipartite_gquadruplex,
        find_gquadruplex,
        find_relaxed_gquadruplex,
        find_bulged_gquadruplex,
        find_imperfect_gquadruplex,
        find_gtriplex
    ]
    for finder in finders:
        motifs = finder(seq)
        for motif in motifs:
            region = range(motif["Start"]-1, motif["End"])
            if not any(mask[i] for i in region):
                results.append(motif)
                for i in region:
                    mask[i] = True
    return sorted(results, key=lambda x: x["Start"])

# ========== 9. I-motif &  VARIANTS ==========

def find_imotif(seq):
    results = []
    pattern = r"(?=(C{3,}\w{1,12}C{3,}\w{1,12}C{3,}\w{1,12}C{3,}))"
    for m in overlapping_finditer(pattern, seq):
        motif_seq = m.group(1)
        score = imotif_g4hunter_score(motif_seq)
        if score >= 0.7:
            c_run_spans = [match.span() for match in re.finditer(r"C{3,}", motif_seq)]
            loops = []
            for i in range(len(c_run_spans)-1):
                loop_start = c_run_spans[i][1]
                loop_end = c_run_spans[i+1][0]
                loops.append(loop_end - loop_start)
            if loops and all(1 <= l <= 7 for l in loops):
                subtype = "Canonical i-Motif"
            elif loops and any(8 <= l <= 12 for l in loops):
                subtype = "Long Loop i-Motif"
            else:
                subtype = "Other i-Motif"
            results.append({
                "Class": "i-Motif",
                "Subtype": subtype,
                "Start": m.start() + 1,
                "End": m.start() + len(motif_seq),
                "Length": len(motif_seq),
                "Sequence": wrap(motif_seq),
                "ScoreMethod": "iM G4HunterStyle",
                "Score": f"{score:.2f}"
            })
    return results

def imotif_g4hunter_score(seq):
    """G4Hunter-style scoring for i-motifs (C-rich sequences)."""
    scores = []
    for c in seq.upper():
        if c == 'C':
            scores.append(1)
        elif c == 'G':
            scores.append(-1)
        else:
            scores.append(0)
    return abs(np.mean(scores)) if scores else 0

# ========== 10. AC-MOTIF ==========

def find_ac_motifs(seq):
    pattern = re.compile(
        r"(?=(?:A{3}[ACGT]{4,6}C{3}[ACGT]{4,6}C{3}[ACGT]{4,6}C{3}|"
        r"C{3}[ACGT]{4,6}C{3}[ACGT]{4,6}C{3}[ACGT]{4,6}A{3}))",
        re.IGNORECASE
    )
    results = []
    for m in pattern.finditer(seq):
        motif_seq = m.group(0).upper()
        score = ac_motif_g4hunter_score(motif_seq)
        if score >= 0.5:  # Apply threshold for G4Hunter-style scoring
            results.append({
                "Class": "AC-Motif",
                "Subtype": "Consensus",
                "Start": m.start() + 1,
                "End": m.start() + len(motif_seq),
                "Length": len(motif_seq),
                "Sequence": wrap(motif_seq),
                "ScoreMethod": "AC G4HunterStyle",
                "Score": f"{score:.2f}"
            })
    return results

def ac_motif_g4hunter_score(seq):
    """G4Hunter-style scoring for AC-motifs."""
    scores = []
    for c in seq.upper():
        if c == 'A':
            scores.append(1)
        elif c == 'C':
            scores.append(1)
        elif c == 'T':
            scores.append(-0.5)
        elif c == 'G':
            scores.append(-0.5)
        else:
            scores.append(0)
    return abs(np.mean(scores)) if scores else 0

# ========== HYBRID MOTIFS ==========

def find_hybrids(motifs, seq):
    events = []
    for idx, m in enumerate(motifs):
        events.append((m['Start'], 'start', idx))
        events.append((m['End'] + 1, 'end', idx))
    events.sort()
    active = set()
    region_start = None
    results = []
    for pos, typ, idx in events:
        if typ == 'start':
            active.add(idx)
            if len(active) == 2:
                region_start = pos
        elif typ == 'end':
            if len(active) == 2:
                region_end = pos - 1
                involved_idxs = list(active)
                involved_classes = {motifs[i]['Class'] for i in involved_idxs}
                if len(involved_classes) >= 2:
                    region_motifs = [motifs[i] for i in involved_idxs]
                    results.append({
                        "Class": "Hybrid",
                        "Subtype": "Overlap of " + " & ".join(sorted(involved_classes)),
                        "Start": region_start,
                        "End": region_end,
                        "Length": region_end - region_start + 1,
                        "Motif Classes": sorted(involved_classes),
                        "Contributing Motifs": region_motifs,
                        "ScoreMethod": "Hybrid Overlap",
                        "Score": f"{min(1.0, len(involved_classes)/5 + len(region_motifs)/10):.2f}",
                        "Sequence": seq[region_start-1:region_end]
                    })
            active.discard(idx)
    return results

# ========== HOTSPOTS / NON-B DNA CLUSTERS ==========

def find_hotspots(motif_hits, seq_len, window=100, min_count=3):
    hotspots = []
    positions = [(hit['Start'], hit['End']) for hit in motif_hits]
    for i in range(0, seq_len - window + 1):
        region_start, region_end = i+1, i+window
        count = sum(s <= region_end and e >= region_start for s,e in positions)
        if count >= min_count:
            motifs_in_region = [m for m in motif_hits if m['Start'] <= region_end and m['End'] >= region_start]
            type_div = len({m['Subtype'] for m in motifs_in_region})
            seq_region = motif_hits[0]['Sequence'] if motif_hits else ""
            hotspots.append({
                "Class": "Non-B DNA Clusters",
                "Subtype": "Hotspot",
                "Start": region_start,
                "End": region_end,
                "Length": region_end - region_start + 1,
                "Sequence": "",  
                "ScoreMethod": "Hotspot",
                "Score": f"{min(1.0, count/10 + type_div/5):.2f}",
                "Motif Count": count,
                "Type Diversity": type_div
            })
    return merge_hotspots(hotspots)

def merge_hotspots(hotspots):
    if not hotspots: return []
    merged = [hotspots[0]]
    for current in hotspots[1:]:
        last = merged[-1]
        if current['Start'] <= last['End']:
            last['End'] = max(last['End'], current['End'])
            last['Length'] = last['End'] - last['Start'] + 1
            last['Motif Count'] += current['Motif Count']
            last['Type Diversity'] = max(last['Type Diversity'], current['Type Diversity'])
            last['Score'] = f"{min(1.0, float(last['Score']) + float(current['Score'])):.2f}"
        else:
            merged.append(current)
    return merged

# ========== NON-OVERLAPPING MOTIF SELECTION ==========

def select_best_nonoverlapping_motifs(motifs: list, motif_priority: list = None) -> list:
    if motif_priority is None:
        motif_priority = [
            'Multimeric G4', 'Bipartite G4', 'Dimeric G4', 'Canonical G4',
            'Relaxed G4', 'Non-canonical G4', 'Bulged G4', 'Three G-Runs'
        ]
    subtype_rank = {subtype: i for i, subtype in enumerate(motif_priority)}
    def motif_key(m):
        rank = subtype_rank.get(m.get('Subtype'), len(subtype_rank))
        try:
            score = float(m.get('Score', 0))
        except ValueError:
            score = 0.0
        length = m.get('Length', 0)
        return (m.get('Class', ''), rank, -score, -length)
    sorted_motifs = sorted(motifs, key=motif_key)
    selected = []
    occupied_per_class = dict()
    for m in sorted_motifs:
        motif_class = m.get('Class', 'Other')
        region = set(range(m['Start'], m['End']+1))
        if motif_class not in occupied_per_class:
            occupied_per_class[motif_class] = set()
        if occupied_per_class[motif_class].isdisjoint(region):
            selected.append(m)
            occupied_per_class[motif_class].update(region)
    return selected

def validate_motif(motif, seq_length):
    required_keys = ["Class", "Subtype", "Start", "End", "Length", "Sequence"]
    if not all(key in motif for key in required_keys):
        return False
    if not (1 <= motif["Start"] <= motif["End"] <= seq_length):
        return False
    if len(motif["Sequence"].replace('\n', '')) == 0:
        return False
    return True

# ========== MASTER MOTIF DISCOVERY ==========

def all_motifs(seq, selected_classes=None, nonoverlap=False, report_hotspots=False):
    """Find all motifs in sequence. Returns both overlapping and non-overlapping results.
    
    Args:
        seq: DNA sequence string
        selected_classes: List of motif classes to analyze (None = all)
        nonoverlap: If True, return only non-overlapping results
        report_hotspots: Whether to include hotspot analysis
        
    Returns:
        If nonoverlap=False: Dictionary with 'overlapping' and 'non_overlapping' keys
        If nonoverlap=True: List of non-overlapping motifs (backward compatibility)
    """
    if not seq or not re.match("^[ATGC]+$", seq, re.IGNORECASE):
        return [] if nonoverlap else {"overlapping": [], "non_overlapping": []}
    
    seq = seq.upper()
    
    # List of all motif finders with display names
    motif_finders = [
        ("Sticky DNA", find_sticky_dna),
        ("Curved DNA", find_curved_DNA),
        ("Z-DNA", find_zdna),
        ("eGZ (Extruded-G)", find_egz_motif),
        ("Slipped DNA", find_slipped_dna),
        ("R-Loop", find_rlfs),
        ("Cruciform", find_cruciform),
        ("Triplex DNA", find_hdna),
        ("G-Triplex", find_gtriplex),
        ("G4", find_gquadruplex),
        ("Relaxed G4", find_relaxed_gquadruplex),
        ("Bulged G4", find_bulged_gquadruplex),
        ("Bipartite G4", find_bipartite_gquadruplex),
        ("Multimeric G4", find_multimeric_gquadruplex),
        ("i-Motif", find_imotif),
        ("AC-Motif", find_ac_motifs),
    ]
    
    # Handle selection logic
    run_all = False
    if selected_classes:
        # Hybrid or Hotspot selection means run all
        if "Hybrid" in selected_classes or "Non-B DNA Clusters" in selected_classes:
            run_all = True
    
    # Run motif finders according to selection
    overlapping_motifs = []
    for name, func in motif_finders:
        if run_all or not selected_classes or name in selected_classes:
            motifs = func(seq)
            # Add conservation metrics to each motif
            for motif in motifs:
                motif = add_conservation_metrics(motif, seq)
            overlapping_motifs += motifs
    
    # Validate motifs
    overlapping_motifs = [m for m in overlapping_motifs if validate_motif(m, len(seq))]
    
    # Create non-overlapping selection for hybrid and cluster analysis
    non_overlapping_motifs = select_best_nonoverlapping_motifs(overlapping_motifs)
    
    # Hybrid and Hotspot cluster analysis uses non-overlapping motifs
    if run_all or (selected_classes and "Hybrid" in selected_classes):
        hybrid_motifs = find_hybrids(non_overlapping_motifs, seq)
        overlapping_motifs += hybrid_motifs
        non_overlapping_motifs += hybrid_motifs
    
    if run_all or (selected_classes and "Non-B DNA Clusters" in selected_classes):
        cluster_motifs = find_hotspots(non_overlapping_motifs, len(seq))
        overlapping_motifs += cluster_motifs
        non_overlapping_motifs += cluster_motifs
    
    # Return format depends on nonoverlap parameter for backward compatibility
    if nonoverlap:
        return non_overlapping_motifs
    else:
        return {
            "overlapping": overlapping_motifs,
            "non_overlapping": non_overlapping_motifs
        }

# Convenience function for getting both overlapping and non-overlapping results
def discover_all_motifs(seq, selected_classes=None):
    """Discover all motifs and return both overlapping and non-overlapping results.
    
    Returns:
        Dictionary with 'overlapping' and 'non_overlapping' keys
    """
    return all_motifs(seq, selected_classes, nonoverlap=False)

# ========== END OF FILE ==========
