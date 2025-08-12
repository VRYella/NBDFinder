import re
import numpy as np

# =========================
# Basic sequence utilities
# =========================

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

def overlapping_finditer(pattern, seq):
    regex = re.compile(pattern, re.IGNORECASE)
    pos = 0
    while pos < len(seq):
        m = regex.search(seq, pos)
        if not m:
            break
        yield m
        pos = m.start() + 1

# =========================
# Curved DNA (PolyA/PolyT) with improved raw scoring
# =========================

def find_polyA_polyT_tracts(seq: str, min_len: int = 7) -> list:
    results = []
    i = 0
    n = len(seq)
    while i < n:
        if seq[i] == 'A' or seq[i] == 'T':
            ch = seq[i]
            start = i
            while i < n and seq[i] == ch:
                i += 1
            if i - start >= min_len:
                results.append((start, i-1, seq[start:i]))
        else:
            i += 1
    return results

def curvature_score(seq):
    # Raw score: length scaled by AT-bias and mild periodicity bonus based on A/T tracts
    if not seq:
        return 0.0
    at_frac = (seq.count('A') + seq.count('T')) / len(seq)
    # count segments of mono-base runs (A or T)
    runs = re.findall(r"(A+|T+)", seq)
    run_bonus = sum(len(r)**0.5 for r in runs)  # diminishing returns
    return len(seq) * (1.0 + at_frac) + 0.5 * run_bonus

def find_global_curved_polyA_polyT(seq: str, min_tract_len: int = 3, min_repeats: int = 3, min_spacing: int = 8, max_spacing: int = 12, min_score: int = 6) -> tuple:
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
                    "Sequence Name": "",
                    "Class": "Curved_DNA",
                    "Subtype": "Global_Curved_Strict_PolyA_or_PolyT",
                    "Start": group[0][0] + 1,
                    "End": group[-1][1] + 1,
                    "Length": group[-1][1] - group[0][0] + 1,
                    "Sequence": wrap(motif_seq),
                    "Score": float(score),
                    "Arms/Repeat Unit/Copies": "",
                    "Spacer": ""
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
                "Sequence Name": "",
                "Class": "Curved_DNA",
                "Subtype": "Local_Curved_Strict_PolyA_or_PolyT",
                "Start": s,
                "End": e,
                "Length": len(tract_seq),
                "Sequence": wrap(tract_seq),
                "Score": float(curvature_score(tract_seq)),
                "Arms/Repeat Unit/Copies": "",
                "Spacer": ""
            })
    return results

def find_curved_DNA(seq: str) -> list:
    global_results, apr_regions = find_global_curved_polyA_polyT(seq)
    local_results = find_local_curved_polyA_polyT(seq, apr_regions)
    return global_results + local_results

# =========================
# Z-DNA seeker (raw scoring retained, no normalization)
# =========================

def zdna_seeker_scoring_array(seq, GC_weight=7.0, AT_weight=0.5, GT_weight=1.25, AC_weight=1.25,
        consecutive_AT_scoring=(0.5, 0.5, 0.5, 0.5, 0.0, 0.0, -5.0, -100.0),
        mismatch_penalty_type="linear",
        mismatch_penalty_starting_value=3,
        mismatch_penalty_linear_delta=3,
        cadence_reward=0.0):
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
                scoring_array[i] = - (mismatch_penalty_starting_value ** mismatches_counter if mismatches_counter < 15 else 32000.0)
            elif mismatch_penalty_type == "linear":
                scoring_array[i] = -mismatch_penalty_starting_value - mismatch_penalty_linear_delta * (mismatches_counter - 1)
            else:
                scoring_array[i] = -10.0
        if t in ("GC", "CG", "GT", "TG", "AC", "CA", "AT", "TA"):
            scoring_array[i] += cadence_reward
    return scoring_array

def find_zdna(seq, threshold=50, drop_threshold=50, GC_weight=7.0, AT_weight=0.5, GT_weight=1.25, AC_weight=1.25,
        consecutive_AT_scoring=(0.5, 0.5, 0.5, 0.5, 0.0, 0.0, -5.0, -100.0),
        mismatch_penalty_type="linear",
        mismatch_penalty_starting_value=3,
        mismatch_penalty_linear_delta=3,
        cadence_reward=0.0):
    seq = seq.upper()
    if len(seq) < 12:
        return []
    scoring = zdna_seeker_scoring_array(seq, GC_weight=GC_weight, AT_weight=AT_weight,
        GT_weight=GT_weight, AC_weight=AC_weight,
        consecutive_AT_scoring=consecutive_AT_scoring,
        mismatch_penalty_type=mismatch_penalty_type,
        mismatch_penalty_starting_value=mismatch_penalty_starting_value,
        mismatch_penalty_linear_delta=mismatch_penalty_linear_delta,
        cadence_reward=cadence_reward)
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
                "Sequence Name": "",
                "Class": "Z-DNA",
                "Subtype": "Z-Seeker",
                "Start": s + 1,
                "End": e + 1,
                "Length": e - s + 1,
                "Sequence": wrap(seq[s:e+1]),
                "Score": float(score),
                "Arms/Repeat Unit/Copies": "",
                "Spacer": ""
            })
            candidate = None
            max_ending_here = current_max = 0
    if candidate:
        s, e, score = candidate
        motifs.append({
            "Sequence Name": "",
            "Class": "Z-DNA",
            "Subtype": "Z-Seeker",
            "Start": s + 1,
            "End": e + 1,
            "Length": e - s + 1,
            "Sequence": wrap(seq[s:e+1]),
            "Score": float(score),
            "Arms/Repeat Unit/Copies": "",
            "Spacer": ""
        })
    return motifs

# =========================
# eGZ (extruded-G) CGG repeats
# =========================

def find_egz_motif(seq):
    pattern = re.compile(r'(CGG){4,}', re.IGNORECASE)
    results = []
    for m in pattern.finditer(seq):
        motif_seq = m.group(0)
        n_repeats = len(motif_seq) // 3
        # Raw score: repeats * unit_len * G-bias
        g_frac = motif_seq.count('G') / len(motif_seq)
        score = n_repeats * 3 * (1.0 + 2.0*g_frac)
        results.append({
            "Sequence Name": "",
            "Family": "Double-stranded",
            "Class": "Z-DNA",
            "Subclass": "eGZ (extruded-G)",
            "Start": m.start() + 1,
            "End": m.end(),
            "Length": len(motif_seq),
            "Sequence": wrap(motif_seq),
            "ScoreMethod": "Repeat_raw",
            "Score": float(score),
            "CGG_Repeats": n_repeats,
            "Arms/Repeat Unit/Copies": f"Unit=CGG;Copies={n_repeats}",
            "Spacer": ""
        })
    return results

# =========================
# Slipped DNA (Direct repeats and STR)
# =========================

def find_slipped_dna(seq):
    results = []
    min_len_dr = 10
    max_len_dr = 300
    # Direct repeats
    for i in range(len(seq) - min_len_dr * 2 + 1):
        for l in range(min_len_dr, min(max_len_dr+1, (len(seq)-i)//2+1)):
            repeat = seq[i:i+l]
            if seq[i+l:i+2*l] == repeat:
                # Raw score: length with composition weight (AT-rich direct repeats more flexible)
                at_frac = (repeat.count('A') + repeat.count('T')) / max(1, len(repeat))
                score = 2*l * (1.0 + 0.5*at_frac)
                results.append({
                    "Sequence Name": "",
                    "Class": "Slipped_DNA",
                    "Subtype": "Direct_Repeat",
                    "Start": i+1,
                    "End": i+2*l,
                    "Length": 2*l,
                    "Sequence": wrap(repeat+repeat),
                    "ScoreMethod": "DR_raw",
                    "Score": float(score),
                    "Arms/Repeat Unit/Copies": f"UnitLen={l};Copies=2",
                    "Spacer": ""
                })
    # STRs
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
                results.append({
                    "Sequence Name": "",
                    "Class": "Slipped_DNA",
                    "Subtype": "STR",
                    "Start": i+1,
                    "End": i + full_len,
                    "Length": full_len,
                    "Unit": repeat_unit,
                    "Copies": reps,
                    "Sequence": wrap(seq[i:i + full_len]),
                    "ScoreMethod": "STR_raw",
                    "Score": float(score),
                    "Arms/Repeat Unit/Copies": f"Unit={repeat_unit};Copies={reps}",
                    "Spacer": ""
                })
                i = i + full_len - 1
                found = True
                break
        if not found:
            i += 1
    return results

# =========================
# R-Loop prediction (RLFS models) with raw stability
# =========================

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
                concat = riz_seq + rez_seq
                g_runs = len(re.findall(r"G{3,}", concat))
                # Raw stability: GC fraction weight + G-run density scaled by length
                gc_frac = gc_content(concat) / 100.0
                score = (gc_frac * 50.0 + g_runs * 10.0) * (len(concat) ** 0.25)
                results.append({
                    "Sequence Name": "",
                    "Class": "R-Loop",
                    "Subtype": f"RLFS_{model_name}",
                    "Start": m.start() + 1,
                    "End": m.start() + len(riz_seq) + rez['end'],
                    "Length": len(riz_seq) + rez['end'],
                    "Sequence": wrap(concat),
                    "ScoreMethod": "QmRLFS_raw",
                    "Score": float(score),
                    "Arms/Repeat Unit/Copies": "",
                    "Spacer": ""
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

# =========================
# Cruciform (Inverted repeats)
# =========================

def find_cruciform(seq):
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
                    # Raw score: arm length with AT-rich bonus minus spacer penalty
                    at_frac = (arm.count('A') + arm.count('T')) / arm_len
                    score = arm_len * (1.0 + 0.5*at_frac) - spacer_len * 2.0
                    results.append({
                        "Sequence Name": "",
                        "Class": "Cruciform",
                        "Subtype": f"Inverted_Repeat_spacer{spacer_len}",
                        "Start": i+1,
                        "End": mid+arm_len,
                        "Length": len(full),
                        "Sequence": wrap(full),
                        "ScoreMethod": "IR_raw",
                        "Score": float(score),
                        "Arms/Repeat Unit/Copies": f"Arms={arm_len}",
                        "Spacer": str(spacer_len)
                    })
    return results

# =========================
# Triplex / Mirror repeats (H-DNA)
# =========================

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
                # Raw score: mirror length with homopurine/pyrimidine enrichment and spacer penalty
                homogeneity = max(pur_frac, pyr_frac)
                score = len(full_seq) * (1.0 + 1.5*homogeneity) - spacer * 1.0
                results.append({
                    "Sequence Name": "",
                    "Class": "Triplex_DNA" if is_triplex else "Mirror_Repeat",
                    "Subtype": "Triplex_Motif" if is_triplex else "Mirror_Repeat",
                    "Start": mirror_start + 1,
                    "End": mirror_end,
                    "Length": len(full_seq),
                    "Spacer": spacer,
                    "Sequence": wrap(full_seq),
                    "PurineFrac": round(pur_frac, 2),
                    "PyrimidineFrac": round(pyr_frac, 2),
                    "Score": float(score),
                    "Arms/Repeat Unit/Copies": f"Arms={rep_len}",
                    "Spacer": str(spacer)
                })
    return results

# =========================
# Sticky DNA (GAA/TTC long repeats)
# =========================

def find_sticky_dna(seq):
    motifs = []
    seq = seq.replace('\n','').replace(' ','').upper()
    pattern = r"(?:GAA){59,}|(?:TTC){59,}"
    for m in re.finditer(pattern, seq):
        repeat_len = len(m.group())
        repeat_count = repeat_len // 3
        # Raw score: repeat_count * unit_length with A/T bias
        at_frac = (m.group().count('A') + m.group().count('T')) / repeat_len
        score = repeat_count * 3 * (1.0 + 0.5*at_frac)
        motifs.append({
            "Sequence Name": "",
            "Class": "Sticky_DNA",
            "Subtype": "GAA_TTC_Repeat",
            "Start": m.start() + 1,
            "End": m.end(),
            "Length": repeat_len,
            "RepeatCount": repeat_count,
            "Sequence": wrap(m.group()),
            "ScoreMethod": "Sakamoto1999_raw",
            "Score": float(score),
            "Arms/Repeat Unit/Copies": f"Unit={'GAA' if 'GAA' in m.group() else 'TTC'};Copies={repeat_count}",
            "Spacer": ""
        })
    return motifs

# =========================
# G4Hunter and G-quadruplex variants (raw scaling)
# =========================

def g4hunter_score(seq):
    scores = []
    for c in seq.upper():
        if c == 'G':
            scores.append(1)
        elif c == 'C':
            scores.append(-1)
        else:
            scores.append(0)
    return np.mean(scores) if scores else 0.0

def find_multimeric_gquadruplex(seq):
    results = []
    pattern = r"(G{3,}\w{1,12}){4,}"
    for m in overlapping_finditer(pattern, seq):
        motif_seq = m.group(0)
        g4h = g4hunter_score(motif_seq)
        if g4h >= 0.5:
            score = (g4h * len(motif_seq)) * 1.2
            results.append({
                "Sequence Name": "",
                "Class": "G4",
                "Subtype": "Multimeric_G4",
                "Start": m.start()+1,
                "End": m.end(),
                "Length": len(motif_seq),
                "Sequence": wrap(motif_seq),
                "ScoreMethod": "G4Hunter_Multimer_raw",
                "Score": float(score),
                "Arms/Repeat Unit/Copies": "",
                "Spacer": ""
            })
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
        score = max(g4hunter_score(unit1), g4hunter_score(unit2)) * len(motif_seq) * 0.9
        if score > 0:
            results.append({
                "Sequence Name": "",
                "Class": "G4",
                "Subtype": "Bipartite_G4",
                "Start": m.start()+1,
                "End": m.end(),
                "Length": len(motif_seq),
                "Sequence": wrap(motif_seq),
                "ScoreMethod": "Bipartite_raw",
                "Score": float(score),
                "Arms/Repeat Unit/Copies": "",
                "Spacer": ""
            })
    return results

def find_gquadruplex(seq):
    pattern = r"(G{3,}\w{1,7}G{3,}\w{1,7}G{3,}\w{1,7}G{3,})"
    results = []
    for m in overlapping_finditer(pattern, seq):
        motif_seq = m.group(1)
        g4h = g4hunter_score(motif_seq)
        score = g4h * len(motif_seq)  # raw
        if g4h >= 0.8:
            results.append({
                "Sequence Name": "",
                "Class": "G4",
                "Subtype": "Canonical_G4",
                "Start": m.start()+1,
                "End": m.end(),
                "Length": len(motif_seq),
                "Sequence": wrap(motif_seq),
                "ScoreMethod": "G4Hunter_v2_raw",
                "Score": float(score),
                "Arms/Repeat Unit/Copies": "",
                "Spacer": ""
            })
    return results

def find_relaxed_gquadruplex(seq):
    pattern = r"(G{3,}\w{8,12}G{3,}\w{8,12}G{3,}\w{8,12}G{3,})"
    results = []
    for m in overlapping_finditer(pattern, seq):
        motif_seq = m.group(1)
        g4h = g4hunter_score(motif_seq)
        score = g4h * len(motif_seq) * 0.8
        if g4h >= 0.5:
            results.append({
                "Sequence Name": "",
                "Class": "G4",
                "Subtype": "Relaxed_G4",
                "Start": m.start()+1,
                "End": m.end(),
                "Length": len(motif_seq),
                "Sequence": wrap(motif_seq),
                "ScoreMethod": "G4Hunter_LongLoop_raw",
                "Score": float(score),
                "Arms/Repeat Unit/Copies": "",
                "Spacer": ""
            })
    return results

def find_bulged_gquadruplex(seq):
    pattern = r"(G{3,}\w{0,3}G{3,}\w{0,3}G{3,}\w{0,3}G{3,})"
    results = []
    for m in overlapping_finditer(pattern, seq):
        motif_seq = m.group(1)
        if len(re.findall(r"G{3,}", motif_seq)) >= 4:
            score = g4hunter_score(motif_seq) * len(motif_seq) * 0.7
            if score > 0:
                results.append({
                    "Sequence Name": "",
                    "Class": "G4",
                    "Subtype": "Bulged_G4",
                    "Start": m.start()+1,
                    "End": m.end(),
                    "Length": len(motif_seq),
                    "Sequence": wrap(motif_seq),
                    "ScoreMethod": "G4Hunter_Bulge_raw",
                    "Score": float(score),
                    "Arms/Repeat Unit/Copies": "",
                    "Spacer": ""
                })
    return results

def find_imperfect_gquadruplex(seq):
    pattern = r"(G{2,3}\w{1,7}G{3,}\w{1,7}G{3,}\w{1,7}G{3,})"\
              r"|(G{3,}\w{1,7}G{2,3}\w{1,7}G{3,}\w{1,7}G{3,})"\
              r"|(G{3,}\w{1,7}G{3,}\w{1,7}G{2,3}\w{1,7}G{3,})"\
              r"|(G{3,}\w{1,7}G{3,}\w{1,7}G{3,}\w{1,7}G{2,3})"
    results = []
    for m in overlapping_finditer(pattern, seq):
        motif_seq = m.group(0)
        g4h = g4hunter_score(motif_seq)
        score = g4h * len(motif_seq)  # raw
        if g4h >= 0.7:
            results.append({
                "Sequence Name": "",
                "Class": "G4",
                "Subtype": "Imperfect_G4",
                "Start": m.start()+1,
                "End": m.end(),
                "Length": len(motif_seq),
                "Sequence": wrap(motif_seq),
                "ScoreMethod": "G4Hunter_Imperfect_raw",
                "Score": float(score),
                "Arms/Repeat Unit/Copies": "",
                "Spacer": ""
            })
    return results

# =========================
# G-triplex
# =========================

def find_gtriplex(seq):
    pattern = r"(G{3,}\w{1,7}G{3,}\w{1,7}G{3,})"
    results = []
    for m in overlapping_finditer(pattern, seq):
        motif_seq = m.group(1)
        g_runs = [len(r) for r in re.findall(r"G{3,}", motif_seq)]
        if len(g_runs) < 3:
            continue
        loops = [len(l) for l in re.findall(r"G{3,}(\w{1,7})G{3,}", motif_seq)]
        loop_term = sum(1/l if l > 0 else 0.5 for l in loops)
        score = (sum(g_runs) * 2.0) + (loop_term * 5.0)
        results.append({
            "Sequence Name": "",
            "Class": "G-Triplex",
            "Subtype": "Three_G-Runs",
            "Start": m.start()+1,
            "End": m.end(),
            "Length": len(motif_seq),
            "Sequence": wrap(motif_seq),
            "ScoreMethod": "G3_raw",
            "Score": float(score),
            "Arms/Repeat Unit/Copies": "",
            "Spacer": ""
        })
    return results

# =========================
# i-Motif
# =========================

def imotif_score(seq):
    c_runs = [len(r) for r in re.findall(r"C{3,}", seq)]
    if len(c_runs) < 4 or len(seq) == 0:
        return 0.0
    c_fraction = seq.count('C') / len(seq)
    c_run_spans = [match.span() for match in re.finditer(r"C{3,}", seq)]
    loops = []
    for i in range(len(c_run_spans)-1):
        loop_start = c_run_spans[i][1]
        loop_end = c_run_spans[i+1][0]
        loops.append(loop_end - loop_start)
    # Raw score: sum of C-run sizes plus compact loop bonuses and C-fraction
    loop_bonus = sum(1.0/(l+1) for l in loops) if loops else 0.5
    return (sum(c_runs) * 1.0) + (c_fraction * len(seq) * 0.5) + (loop_bonus * 3.0)

def find_imotif(seq):
    results = []
    pattern = r"(?=(C{3,}\w{1,12}C{3,}\w{1,12}C{3,}\w{1,12}C{3,}))"
    for m in overlapping_finditer(pattern, seq):
        motif_seq = m.group(1)
        score = imotif_score(motif_seq)
        if score > 0:
            c_run_spans = [match.span() for match in re.finditer(r"C{3,}", motif_seq)]
            loops = []
            for i in range(len(c_run_spans)-1):
                loop_start = c_run_spans[i][1]
                loop_end = c_run_spans[i+1][0]
                loops.append(loop_end - loop_start)
            if loops and all(1 <= l <= 7 for l in loops):
                subtype = "Canonical_iMotif"
            elif loops and any(8 <= l <= 12 for l in loops):
                subtype = "LongLoop_iMotif"
            else:
                subtype = "Other_iMotif"
            results.append({
                "Sequence Name": "",
                "Class": "i-Motif",
                "Subtype": subtype,
                "Start": m.start() + 1,
                "End": m.start() + len(motif_seq),
                "Length": len(motif_seq),
                "Sequence": wrap(motif_seq),
                "ScoreMethod": "iM_raw",
                "Score": float(score),
                "Arms/Repeat Unit/Copies": "",
                "Spacer": ""
            })
    return results

# =========================
# AC-motifs (consensus)
# =========================

def find_ac_motifs(seq):
    pattern = re.compile(
        r"(?=(?:A{3}[ACGT]{4,6}C{3}[ACGT]{4,6}C{3}[ACGT]{4,6}C{3}|"
        r"C{3}[ACGT]{4,6}C{3}[ACGT]{4,6}C{3}[ACGT]{4,6}A{3}))",
        re.IGNORECASE
    )
    results = []
    for m in pattern.finditer(seq):
        motif_seq = m.group(0).upper()
        # Raw score: length times boundary run emphasis (A3 and C3 presence)
        boundary_bonus = (3 if motif_seq.startswith('AAA') else 0) + (3 if motif_seq.endswith('AAA') else 0)
        c_runs = len(re.findall(r"C{3}", motif_seq))
        score = len(motif_seq) + boundary_bonus + c_runs * 2.0
        results.append({
            "Sequence Name": "",
            "Class": "AC-Motif",
            "Subtype": "Consensus",
            "Start": m.start() + 1,
            "End": m.start() + len(motif_seq),
            "Length": len(motif_seq),
            "Sequence": wrap(motif_seq),
            "ScoreMethod": "PatternMatch_raw",
            "Score": float(score),
            "Arms/Repeat Unit/Copies": "",
            "Spacer": ""
        })
    return results

# =========================
# Hybrid overlaps and hotspots (augmented fields)
# =========================

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
                    score = sum(float(m.get("Score", 0.0)) for m in region_motifs) * 0.1
                    results.append({
                        "Sequence Name": motifs[0].get("Sequence Name", ""),
                        "Class": "Hybrid",
                        "Subtype": "_".join(sorted(involved_classes)) + "_Overlap",
                        "Start": region_start,
                        "End": region_end,
                        "Length": region_end - region_start + 1,
                        "MotifClasses": sorted(involved_classes),
                        "ContributingMotifs": region_motifs,
                        "ScoreMethod": "HybridOverlap_raw",
                        "Score": float(score),
                        "Sequence": wrap(seq[region_start-1:region_end]),
                        "Arms/Repeat Unit/Copies": "",
                        "Spacer": ""
                    })
            active.discard(idx)
    return results

def find_hotspots(motif_hits, seq_len, window=100, min_count=3):
    hotspots = []
    positions = [(hit['Start'], hit['End']) for hit in motif_hits]
    for i in range(0, seq_len - window + 1):
        region_start, region_end = i+1, i+window
        count = sum(s <= region_end and e >= region_start for s,e in positions)
        if count >= min_count:
            motifs_in_region = [m for m in motif_hits if m['Start'] <= region_end and m['End'] >= region_start]
            type_div = len({m['Subtype'] for m in motifs_in_region})
            total_score = sum(float(m.get("Score", 0.0)) for m in motifs_in_region)
            hotspots.append({
                "Sequence Name": motif_hits[0].get("Sequence Name", "") if motif_hits else "",
                "Class": "Non-B DNA Clusters",
                "Subtype": "Hotspot",
                "Start": region_start,
                "End": region_end,
                "Length": region_end - region_start + 1,
                "Sequence": "",
                "ScoreMethod": "Hotspot_raw",
                "Score": float(total_score),
                "MotifCount": count,
                "TypeDiversity": type_div,
                "Arms/Repeat Unit/Copies": "",
                "Spacer": ""
            })
    return merge_hotspots(hotspots)

def merge_hotspots(hotspots):
    if not hotspots:
        return []
    merged = [hotspots[0]]
    for current in hotspots[1:]:
        last = merged[-1]
        if current['Start'] <= last['End']:
            last['End'] = max(last['End'], current['End'])
            last['Length'] = last['End'] - last['Start'] + 1
            last['MotifCount'] += current['MotifCount']
            last['TypeDiversity'] = max(last['TypeDiversity'], current['TypeDiversity'])
            last['Score'] = float(last['Score']) + float(current['Score'])
        else:
            merged.append(current)
    return merged

# =========================
# Selection, validation, stats
# =========================

def select_best_nonoverlapping_motifs(motifs: list, motif_priority: list = None) -> list:
    if motif_priority is None:
        motif_priority = [
            'Multimeric_G4', 'Bipartite_G4', 'Dimeric_G4', 'Canonical_G4',
            'Relaxed_G4', 'Non_canonical_G4', 'Bulged_G4', 'Three_G-Runs'
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

def get_basic_stats(seq, motifs=None):
    seq = seq.upper()
    length = len(seq)
    gc = gc_content(seq)
    at = (seq.count('A') + seq.count('T')) / length * 100 if length else 0
    stats = {
        "Length": length,
        "GC%": round(gc, 2),
        "AT%": round(at, 2),
        "A": seq.count('A'),
        "T": seq.count('T'),
        "G": seq.count('G'),
        "C": seq.count('C'),
    }
    if motifs is not None:
        covered = set()
        for m in motifs:
            covered.update(range(m['Start'], m['End']))
        coverage_pct = (len(covered) / length * 100) if length else 0
        stats["Motif Coverage %"] = round(coverage_pct, 2)
    return stats

# =========================
# Aggregator
# =========================

def all_motifs(seq, nonoverlap=False, report_hotspots=False, sequence_name="Sequence"):
    if not seq or not re.match("^[ATGC]+$", seq, re.IGNORECASE):
        return []
    seq = seq.upper()
    motif_list = (
        find_sticky_dna(seq) +
        find_curved_DNA(seq) +
        find_zdna(seq) +
        find_egz_motif(seq) +
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
        find_ac_motifs(seq)
    )
    # Validate and standardize fields
    motif_list = [m for m in motif_list if validate_motif(m, len(seq))]
    # Add hybrids
    motif_list += find_hybrids(motif_list, seq)
    # De-overlap per class if asked
    if nonoverlap:
        motif_list = select_best_nonoverlapping_motifs(motif_list)
    # Hotspots appended if asked
    if report_hotspots:
        motif_list += find_hotspots(motif_list, len(seq))
    # Add Sequence Name and ensure ordered keys exist
    for m in motif_list:
        m["Sequence Name"] = sequence_name
        # Ensure mandatory ordered fields exist and not missing
        if "Arms/Repeat Unit/Copies" not in m:
            m["Arms/Repeat Unit/Copies"] = ""
        if "Spacer" not in m:
            m["Spacer"] = ""
        if "Score" in m:
            try:
                m["Score"] = float(m["Score"])
            except Exception:
                pass
    return motif_list

# =========================
# Utility: formatted output rows in the exact requested order
# =========================

def format_motif_rows(motifs):
    ordered = []
    for m in motifs:
        row = {
            "Sequence Name": m.get("Sequence Name", ""),
            "Class": m.get("Class", ""),
            "Subtype": m.get("Subtype", m.get("Subclass", "")),
            "Start": m.get("Start", ""),
            "End": m.get("End", ""),
            "Length": m.get("Length", ""),
            "Sequence": m.get("Sequence", ""),
            "Score": m.get("Score", ""),
            "Arms/Repeat Unit/Copies": m.get("Arms/Repeat Unit/Copies", ""),
            "Spacer": m.get("Spacer", "")
        }
        ordered.append(row)
    return ordered

# =========================
