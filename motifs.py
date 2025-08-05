import re
import numpy as np

def wrap(seq, width=60):
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))

def gc_content(seq):
    gc = seq.count('G') + seq.count('C')
    return (gc / max(1, len(seq))) * 100 if seq else 0

def reverse_complement(seq):
    return seq.translate(str.maketrans("ATGC", "TACG"))[::-1]

def overlapping_finditer(pattern, seq):
    regex = re.compile(pattern, re.IGNORECASE)
    pos = 0
    while pos < len(seq):
        m = regex.search(seq, pos)
        if not m:
            break
        yield m
        pos = m.start() + 1

# --------- SCORING FUNCTIONS ---------
def score_curved(motif):
    # Best practice: tract length/40, max 1.0
    return round(min(1.0, motif['Length']/40), 2)

def score_slipped(motif):
    # Direct repeat: length/300; STR: copies/20
    if motif['Subtype'] == "Direct_Repeat":
        return round(min(1.0, motif['Length']/300), 2)
    else:
        return round(min(1.0, motif.get('Copies', 1)/20), 2)

def score_cruciform(motif):
    # Arm length/100 + AT%*0.3, max 1.0
    seq = motif['Sequence'].replace('\n', '')
    arm_len = motif['Length']//2
    at_frac = (seq.count('A') + seq.count('T')) / len(seq)
    return round(min(1.0, arm_len/100 + at_frac*0.3), 2)

def score_mirror_triplex(motif):
    # max(purineFrac, pyrimidineFrac) + length/50 - spacer/10
    score = max(motif.get('PurineFrac', 0), motif.get('PyrimidineFrac', 0))
    score += min(1.0, motif['Length']/50)
    score -= motif.get('Spacer', 0)/10
    return round(min(2.0, score), 2)

def score_gtriplex(motif):
    g_runs = motif.get('G_runs', 3)
    avg_loop = motif.get('AverageLoopLength', 0)
    score = g_runs/3 - avg_loop/7 + min(1.0, motif['Length']/30)
    return round(min(2.0, score), 2)

def score_sticky(motif):
    return round(min(1.0, motif.get('RepeatCount', 0)/200), 2)

def score_imotif(motif):
    seq = motif['Sequence'].replace('\n', '')
    c_runs = len(re.findall(r"C{3,}", seq))
    c_fraction = seq.count('C') / len(seq) if seq else 0
    c_run_spans = [match.span() for match in re.finditer(r"C{3,}", seq)]
    loops = []
    for i in range(len(c_run_spans)-1):
        loop_start = c_run_spans[i][1]
        loop_end = c_run_spans[i+1][0]
        loops.append(loop_end - loop_start)
    loop_score = sum(1/(l+1) for l in loops) / max(1, len(loops)) if loops else 0.5
    score = min(1.0, c_runs/6 + c_fraction*0.5 + loop_score*0.3)
    return round(score, 2)

def score_acmotif(motif):
    return 1.0

def score_hybrid(motif):
    diversity = len(set(motif.get('MotifClasses', [])))
    motif_count = len(motif.get('ContributingMotifs', []))
    return round(min(1.0, diversity/5 + motif_count/10), 2)

def score_hotspot(motif):
    return round(min(1.0, motif.get('MotifCount', 0)/10 + motif.get('TypeDiversity', 0)/5), 2)

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

# --------- MOTIF FINDERS ---------
def find_curved_DNA(seq):
    def find_polyA_polyT_tracts(seq, min_len=7):
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
    def find_global(seq):
        tracts = find_polyA_polyT_tracts(seq, 3)
        results = []
        apr_regions = []
        min_repeats, min_spacing, max_spacing = 3, 8, 12
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
                motif = {
                    "Class": "Curved_DNA",
                    "Subtype": "Global_Curved_Strict_PolyA_or_PolyT",
                    "Start": group[0][0] + 1,
                    "End": group[-1][1] + 1,
                    "Length": group[-1][1] - group[0][0] + 1,
                    "Sequence": wrap(motif_seq),
                    "ScoreMethod": "Strict: PolyA/PolyT curvature tract score"
                }
                motif["Score"] = score_curved(motif)
                results.append(motif)
                apr_regions.append((motif["Start"], motif["End"]))
        return results, apr_regions
    def find_local(seq, apr_regions):
        results = []
        tracts = find_polyA_polyT_tracts(seq, 7)
        for start, end, tract_seq in tracts:
            s, e = start + 1, end + 1
            if not any(r_start <= s <= r_end or r_start <= e <= r_end for r_start, r_end in apr_regions):
                motif = {
                    "Class": "Curved_DNA",
                    "Subtype": "Local_Curved_Strict_PolyA_or_PolyT",
                    "Start": s,
                    "End": e,
                    "Length": len(tract_seq),
                    "Sequence": wrap(tract_seq),
                    "ScoreMethod": "Strict: PolyA/PolyT tract length"
                }
                motif["Score"] = score_curved(motif)
                results.append(motif)
        return results
    global_results, apr_regions = find_global(seq)
    local_results = find_local(seq, apr_regions)
    return global_results + local_results

def find_slipped_dna(seq):
    results = []
    min_len_dr = 10
    max_len_dr = 300
    for i in range(len(seq) - min_len_dr * 2 + 1):
        for l in range(min_len_dr, min(max_len_dr+1, (len(seq)-i)//2+1)):
            repeat = seq[i:i+l]
            if seq[i+l:i+2*l] == repeat:
                motif = {
                    "Class": "Slipped_DNA",
                    "Subtype": "Direct_Repeat",
                    "Start": i+1,
                    "End": i+2*l,
                    "Length": 2*l,
                    "Sequence": wrap(repeat+repeat),
                    "ScoreMethod": "nBST_DR"
                }
                motif["Score"] = score_slipped(motif)
                results.append(motif)
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
                motif = {
                    "Class": "Slipped_DNA",
                    "Subtype": "STR",
                    "Start": i+1,
                    "End": i + reps*unit + remainder,
                    "Length": reps*unit + remainder,
                    "Unit": repeat_unit,
                    "Copies": reps,
                    "Sequence": wrap(seq[i:i + reps*unit + remainder]),
                    "ScoreMethod": "nBST_STR"
                }
                motif["Score"] = score_slipped(motif)
                results.append(motif)
                i = i + reps*unit + remainder - 1
                found = True
                break
        if not found:
            i += 1
    return results

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
                    motif = {
                        "Class": "Cruciform",
                        "Subtype": f"Inverted_Repeat_spacer{spacer_len}",
                        "Start": i+1,
                        "End": mid+arm_len,
                        "Length": len(full),
                        "Sequence": wrap(full),
                        "ScoreMethod": "nBST_IR"
                    }
                    motif["Score"] = score_cruciform(motif)
                    results.append(motif)
    return results

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
                pur_frac = (full_seq.count('A') + full_seq.count('G')) / max(1, len(full_seq))
                pyr_frac = (full_seq.count('C') + full_seq.count('T')) / max(1, len(full_seq))
                is_triplex = (pur_frac >= 0.9 or pyr_frac >= 0.9)
                motif = {
                    "Class": "Triplex_DNA" if is_triplex else "Mirror_Repeat",
                    "Subtype": "Triplex_Motif" if is_triplex else "Mirror_Repeat",
                    "Start": mirror_start + 1,
                    "End": mirror_end,
                    "Length": len(full_seq),
                    "Spacer": spacer,
                    "Sequence": wrap(full_seq),
                    "PurineFrac": round(pur_frac, 2),
                    "PyrimidineFrac": round(pyr_frac, 2)
                }
                motif["Score"] = score_mirror_triplex(motif)
                results.append(motif)
    return results

def find_gtriplex(seq):
    pattern = r"(G{3,}\w{1,7}G{3,}\w{1,7}G{3,})"
    results = []
    for m in overlapping_finditer(pattern, seq):
        motif_seq = m.group(1)
        g_runs = [len(r) for r in re.findall(r"G{3,}", motif_seq)]
        if len(g_runs) < 3: continue
        loops = [len(l) for l in re.findall(r"G{3,}(\w{1,7})G{3,}", motif_seq)]
        motif = {
            "Class": "G-Triplex",
            "Subtype": "Three_G-Runs",
            "Start": m.start()+1,
            "End": m.end(),
            "Length": len(motif_seq),
            "Sequence": wrap(motif_seq),
            "G_runs": len(g_runs),
            "AverageLoopLength": np.mean(loops) if loops else 0
        }
        motif["Score"] = score_gtriplex(motif)
        results.append(motif)
    return results

def find_sticky_dna(seq):
    motifs = []
    seq = seq.replace('\n','').replace(' ','').upper()
    pattern = r"(?:GAA){59,}|(?:TTC){59,}"
    for m in re.finditer(pattern, seq):
        repeat_count = len(m.group()) // 3
        motif = {
            "Class": "Sticky_DNA",
            "Subtype": "GAA_TTC_Repeat",
            "Start": m.start() + 1,
            "End": m.end(),
            "Length": len(m.group()),
            "RepeatCount": repeat_count,
            "Sequence": m.group(),
            "ScoreMethod": "Sakamoto1999"
        }
        motif["Score"] = score_sticky(motif)
        motifs.append(motif)
    return motifs

def find_imotif(seq):
    results = []
    pattern = r"(?=(C{3,}\w{1,12}C{3,}\w{1,12}C{3,}\w{1,12}C{3,}))"
    for m in overlapping_finditer(pattern, seq):
        motif_seq = m.group(1)
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
        motif = {
            "Class": "i-Motif",
            "Subtype": subtype,
            "Start": m.start() + 1,
            "End": m.start() + len(motif_seq),
            "Length": len(motif_seq),
            "Sequence": wrap(motif_seq),
            "ScoreMethod": "iM_G4HunterStyle"
        }
        motif["Score"] = score_imotif(motif)
        results.append(motif)
    return results

def find_ac_motifs(seq):
    pattern = re.compile(
        r"(?=(?:A{3}[ACGT]{4,6}C{3}[ACGT]{4,6}C{3}[ACGT]{4,6}C{3}|"
        r"C{3}[ACGT]{4,6}C{3}[ACGT]{4,6}C{3}[ACGT]{4,6}A{3}))",
        re.IGNORECASE
    )
    results = []
    for m in pattern.finditer(seq):
        motif_seq = m.group(0).upper()
        motif = {
            "Class": "AC-Motif",
            "Subtype": "Consensus",
            "Start": m.start() + 1,
            "End": m.start() + len(motif_seq),
            "Length": len(motif_seq),
            "Sequence": wrap(motif_seq),
            "ScoreMethod": "PatternMatch"
        }
        motif["Score"] = score_acmotif(motif)
        results.append(motif)
    return results

# --- G4, Z-DNA, RLFS, EGZ (original logic, original scoring) ---
def find_gquadruplex(seq):
    pattern = r"(G{3,}\w{1,7}G{3,}\w{1,7}G{3,}\w{1,7}G{3,})"
    results = []
    for m in overlapping_finditer(pattern, seq):
        motif_seq = m.group(1)
        score = g4hunter_score(motif_seq)
        if score >= 1.2:
            results.append({
                "Class": "G4",
                "Subtype": "Canonical_G4",
                "Start": m.start()+1,
                "End": m.end(),
                "Length": len(motif_seq),
                "Sequence": wrap(motif_seq),
                "ScoreMethod": "G4Hunter_v2",
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
                "Class": "G4",
                "Subtype": "Relaxed_G4",
                "Start": m.start()+1,
                "End": m.end(),
                "Length": len(motif_seq),
                "Sequence": wrap(motif_seq),
                "ScoreMethod": "G4Hunter_LongLoop",
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
                "Class": "G4",
                "Subtype": "Bulged_G4",
                "Start": m.start()+1,
                "End": m.end(),
                "Length": len(motif_seq),
                "Sequence": wrap(motif_seq),
                "ScoreMethod": "G4Hunter_Bulge",
                "Score": f"{score*0.7:.2f}"})
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
                "Class": "G4",
                "Subtype": "Bipartite_G4",
                "Start": m.start()+1,
                "End": m.end(),
                "Length": len(motif_seq),
                "Sequence": wrap(motif_seq),
                "ScoreMethod": "Bipartite_Score",
                "Score": f"{score:.2f}"})
    return results

def find_multimeric_gquadruplex(seq):
    results = []
    pattern = r"(G{3,}\w{1,12}){4,}"
    for m in overlapping_finditer(pattern, seq):
        motif_seq = m.group(0)
        if g4hunter_score(motif_seq) >= 1.0:
            results.append({
                "Class": "G4",
                "Subtype": "Multimeric_G4",
                "Start": m.start()+1,
                "End": m.end(),
                "Length": len(motif_seq),
                "Sequence": wrap(motif_seq),
                "ScoreMethod": "G4Hunter_Multimer",
                "Score": f"{g4hunter_score(motif_seq)*1.2:.2f}"})
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
                "Subtype": "Imperfect_G4",
                "Start": m.start()+1,
                "End": m.end(),
                "Length": len(motif_seq),
                "Sequence": wrap(motif_seq),
                "ScoreMethod": "G4Hunter_Imperfect",
                "Score": f"{score:.2f}"})
    return results

def find_zdna(seq, threshold=50, drop_threshold=50, GC_weight=7.0, AT_weight=0.5, GT_weight=1.25, AC_weight=1.25,
        consecutive_AT_scoring=(0.5, 0.5, 0.5, 0.5, 0.0, 0.0, -5.0, -100.0),
        mismatch_penalty_type="linear",
        mismatch_penalty_starting_value=3,
        mismatch_penalty_linear_delta=3,
        cadence_reward=0.0):
    seq = seq.upper()
    if len(seq) < 12:
        return []
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
                    scoring_array[i] = -mismatch_penalty_starting_value ** mismatches_counter if mismatches_counter < 15 else -32000
                elif mismatch_penalty_type == "linear":
                    scoring_array[i] = -mismatch_penalty_starting_value - mismatch_penalty_linear_delta * (mismatches_counter - 1)
                else:
                    scoring_array[i] = -10
            if t in ("GC", "CG", "GT", "TG", "AC", "CA", "AT", "TA"):
                scoring_array[i] += cadence_reward
        return scoring_array
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

def find_egz_motif(seq):
    pattern = re.compile(r'(CGG){4,}', re.IGNORECASE)
    results = []
    for m in pattern.finditer(seq):
        motif_seq = m.group(0)
        n_repeats = len(motif_seq) // 3
        score = min(1.0, n_repeats/12)
        results.append({
            "Family": "Double-stranded",
            "Class": "Z-DNA",
            "Subclass": "eGZ (extruded-G)",
            "Start": m.start() + 1,
            "End": m.end(),
            "Length": len(motif_seq),
            "Sequence": wrap(motif_seq),
            "ScoreMethod": "Repeat_normalized",
            "Score": f"{score:.2f}",
            "CGG_Repeats": n_repeats
        })
    return results

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

# --- Hybrid/Hotspot ---
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
                    motif = {
                        "Class": "Hybrid",
                        "Subtype": "_".join(sorted(involved_classes)) + "_Overlap",
                        "Start": region_start,
                        "End": region_end,
                        "Length": region_end - region_start + 1,
                        "MotifClasses": sorted(involved_classes),
                        "ContributingMotifs": region_motifs,
                        "ScoreMethod": "HybridOverlap",
                        "Sequence": seq[region_start-1:region_end]
                    }
                    motif["Score"] = score_hybrid(motif)
                    results.append(motif)
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
            seq_region = motif_hits[0]['Sequence'] if motif_hits else ""
            motif = {
                "Class": "Non-B DNA Clusters",
                "Subtype": "Hotspot",
                "Start": region_start,
                "End": region_end,
                "Length": region_end - region_start + 1,
                "Sequence": "",
                "ScoreMethod": "Hotspot",
                "MotifCount": count,
                "TypeDiversity": type_div
            }
            motif["Score"] = score_hotspot(motif)
            hotspots.append(motif)
    # merge hotspots
    if not hotspots: return []
    merged = [hotspots[0]]
    for current in hotspots[1:]:
        last = merged[-1]
        if current['Start'] <= last['End']:
            last['End'] = max(last['End'], current['End'])
            last['Length'] = last['End'] - last['Start'] + 1
            last['MotifCount'] += current['MotifCount']
            last['TypeDiversity'] = max(last['TypeDiversity'], current['TypeDiversity'])
            last['Score'] = score_hotspot(last)
        else:
            merged.append(current)
    return merged

def validate_motif(motif, seq_length):
    required_keys = ["Class", "Subtype", "Start", "End", "Length", "Sequence"]
    if not all(key in motif for key in required_keys):
        return False
    if not (1 <= motif["Start"] <= motif["End"] <= seq_length):
        return False
    if len(motif["Sequence"].replace('\n', '')) == 0:
        return False
    return True

def all_motifs(seq, nonoverlap=False, report_hotspots=False):
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
        find_imperfect_gquadruplex(seq) +
        find_imotif(seq) +
        find_ac_motifs(seq)
    )
    motif_list = [m for m in motif_list if validate_motif(m, len(seq))]
    motif_list += find_hybrids(motif_list, seq)
    if nonoverlap:
        motif_list = select_best_nonoverlapping_motifs(motif_list)
    if report_hotspots:
        motif_list += find_hotspots(motif_list, len(seq))
    return motif_list

# (You may want to include select_best_nonoverlapping_motifs, get_basic_stats, etc. from your repo for completeness.)

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

def all_motifs(seq, nonoverlap=False, report_hotspots=False):
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
    motif_list = [m for m in motif_list if validate_motif(m, len(seq))]
    motif_list += find_hybrids(motif_list, seq)
    if nonoverlap:
        motif_list = select_best_nonoverlapping_motifs(motif_list)
    if report_hotspots:
        motif_list += find_hotspots(motif_list, len(seq))
    return motif_list
