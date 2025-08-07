import re
import numpy as np
import random

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

def scientific_conservation_score(seq, motif_seq, n_shuffles=100):
    k = min(8, len(motif_seq))
    if k < 4: k = 4
    motif_kmer = motif_seq[:k]
    obs = seq.count(motif_kmer)
    if obs == 0:
        return 0.0, 1.0, "not significant"
    shuf_counts = []
    seq_list = list(seq)
    for _ in range(n_shuffles):
        random.shuffle(seq_list)
        shuf_seq = ''.join(seq_list)
        shuf_counts.append(shuf_seq.count(motif_kmer))
    exp = np.mean(shuf_counts)
    score = np.log2(obs / (exp+1e-6))
    pvalue = sum(count >= obs for count in shuf_counts) / n_shuffles
    if pvalue > 0.05:
        signif = "not significant"
    elif score < 0.5:
        signif = "low"
    elif score < 1.5:
        signif = "medium"
    elif score < 2.5:
        signif = "high"
    else:
        signif = "very high"
    return round(score,3), round(pvalue,4), signif

def validate_motif(motif, seq_length):
    required_keys = ["Class", "Subtype", "Start", "End", "Length", "Sequence"]
    if not all(key in motif for key in required_keys):
        return False
    if not (1 <= motif["Start"] <= motif["End"] <= seq_length):
        return False
    if len(motif["Sequence"].replace('\n', '')) == 0:
        return False
    return True

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
            cons_score, pval, cons_sig = scientific_conservation_score(seq, motif_seq)
            motif = {
                "Class": "Curved_DNA",
                "Subtype": "Global_Curved_Strict_PolyA_or_PolyT",
                "Start": group[0][0] + 1,
                "End": group[-1][1] + 1,
                "Length": group[-1][1] - group[0][0] + 1,
                "Sequence": wrap(motif_seq),
                "Stability Score": score,
                "Significance": "high" if score > 10 else "medium",
                "Conservation score": cons_score,
                "Conservation p-value": pval,
                "Conservation significance": cons_sig,
                "Repeat/arm/Contributing region": motif_seq,
                "Spacer": None
            }
            results.append(motif)
            apr_regions.append((motif["Start"], motif["End"]))
    return results, apr_regions

def find_local_curved_polyA_polyT(seq: str, apr_regions: list, min_len: int = 7) -> list:
    results = []
    tracts = find_polyA_polyT_tracts(seq, min_len)
    for start, end, tract_seq in tracts:
        s, e = start + 1, end + 1
        cons_score, pval, cons_sig = scientific_conservation_score(seq, tract_seq)
        if not any(r_start <= s <= r_end or r_start <= e <= r_end for r_start, r_end in apr_regions):
            results.append({
                "Class": "Curved_DNA",
                "Subtype": "Local_Curved_Strict_PolyA_or_PolyT",
                "Start": s,
                "End": e,
                "Length": len(tract_seq),
                "Sequence": wrap(tract_seq),
                "Stability Score": len(tract_seq),
                "Significance": "medium" if len(tract_seq) > 10 else "low",
                "Conservation score": cons_score,
                "Conservation p-value": pval,
                "Conservation significance": cons_sig,
                "Repeat/arm/Contributing region": tract_seq,
                "Spacer": None
            })
    return results

def find_curved_DNA(seq: str) -> list:
    global_results, apr_regions = find_global_curved_polyA_polyT(seq)
    local_results = find_local_curved_polyA_polyT(seq, apr_regions)
    return global_results + local_results

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
                    region_seq = seq[region_start-1:region_end]
                    cons_score, pval, cons_sig = scientific_conservation_score(seq, region_seq)
                    results.append({
                        "Class": "Hybrid",
                        "Subtype": "_".join(sorted(involved_classes)) + "_Overlap",
                        "Start": region_start,
                        "End": region_end,
                        "Length": region_end - region_start + 1,
                        "Sequence": wrap(region_seq),
                        "Stability Score": len(region_seq),
                        "Significance": "high" if len(involved_classes) > 2 else "medium",
                        "Conservation score": cons_score,
                        "Conservation p-value": pval,
                        "Conservation significance": cons_sig,
                        "Repeat/arm/Contributing region": None,
                        "Spacer": None,
                        "MotifClasses": sorted(involved_classes),
                        "ContributingMotifs": region_motifs
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
            seq_region = ''.join([m['Sequence'].replace('\n', '') for m in motifs_in_region])
            cons_score, pval, cons_sig = scientific_conservation_score(seq_region, seq_region)
            hotspots.append({
                "Class": "Non-B DNA Clusters",
                "Subtype": "Hotspot",
                "Start": region_start,
                "End": region_end,
                "Length": region_end - region_start + 1,
                "Sequence": wrap(seq_region),
                "Stability Score": count,
                "Significance": "high" if type_div > 3 else "medium",
                "Conservation score": cons_score,
                "Conservation p-value": pval,
                "Conservation significance": cons_sig,
                "Repeat/arm/Contributing region": None,
                "Spacer": None,
                "MotifCount": count,
                "TypeDiversity": type_div
            })
    return merge_hotspots(hotspots)

def find_zdna(seq, sequence_name=""):
    seq = seq.upper()
    if len(seq) < 12:
        return []
    # Use your existing zdna_seeker_scoring_array function for window scoring
    scoring = zdna_seeker_scoring_array(seq)
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
        if max_ending_here >= 50 and (candidate is None or current_max < max_ending_here):
            candidate = (start_idx, end_idx, max_ending_here)
            current_max = max_ending_here
        if candidate and (max_ending_here < 0 or current_max - max_ending_here >= 50):
            s, e, score = candidate
            motif_seq = seq[s:e+1]
            cons_score, pval, cons_sig = scientific_conservation_score(seq, motif_seq)
            motifs.append({
                "SequenceName": sequence_name,
                "Class": "Z-DNA",
                "Subtype": "Z-Seeker",
                "Start": s + 1,
                "End": e + 1,
                "Length": e - s + 1,
                "sequence": wrap(motif_seq),
                "Stability Score": score,
                "Significance": "high" if score > 80 else "medium",
                "Conservation score": cons_score,
                "Conservation p-value": pval,
                "Conservation significance": cons_sig,
                "Repeat/arm/Contributing region": motif_seq,
                "Spacer": None
            })
            candidate = None
            max_ending_here = current_max = 0
    if candidate:
        s, e, score = candidate
        motif_seq = seq[s:e+1]
        cons_score, pval, cons_sig = scientific_conservation_score(seq, motif_seq)
        motifs.append({
            "SequenceName": sequence_name,
            "Class": "Z-DNA",
            "Subtype": "Z-Seeker",
            "Start": s + 1,
            "End": e + 1,
            "Length": e - s + 1,
            "sequence": wrap(motif_seq),
            "Stability Score": score,
            "Significance": "high" if score > 80 else "medium",
            "Conservation score": cons_score,
            "Conservation p-value": pval,
            "Conservation significance": cons_sig,
            "Repeat/arm/Contributing region": motif_seq,
            "Spacer": None
        })
    return motifs
def find_mirror_repeat(seq, sequence_name=""):
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
                if not is_triplex:  # Only non-triplex mirror repeats
                    cons_score, pval, cons_sig = scientific_conservation_score(seq, full_seq)
                    results.append({
                        "SequenceName": sequence_name,
                        "Class": "Mirror_Repeat",
                        "Subtype": "Non-triplex_Mirror_Repeat",
                        "Start": mirror_start + 1,
                        "End": mirror_end,
                        "Length": len(full_seq),
                        "sequence": wrap(full_seq),
                        "Stability Score": len(full_seq),
                        "Significance": "medium" if len(full_seq) > 20 else "low",
                        "Conservation score": cons_score,
                        "Conservation p-value": pval,
                        "Conservation significance": cons_sig,
                        "Repeat/arm/Contributing region": repeat,
                        "Spacer": spacer
                    })
    return results

def find_str(seq, sequence_name=""):
    min_unit_str = 1
    max_unit_str = 6
    min_reps_str = 5
    min_len_str = 15
    results = []
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
            while (i + reps*unit + unit <= len(seq) and seq[i + reps*unit:i + (reps+1)*unit] == repeat_unit):
                reps += 1
            if reps >= min_reps_str and reps*unit >= min_len_str:
                remainder = 0
                rs = i + reps*unit
                re_idx = rs
                while (re_idx < len(seq) and seq[re_idx] == repeat_unit[re_idx % unit]):
                    remainder += 1
                    re_idx += 1
                motif_seq = seq[i:i + reps*unit + remainder]
                cons_score, pval, cons_sig = scientific_conservation_score(seq, motif_seq)
                results.append({
                    "SequenceName": sequence_name,
                    "Class": "STR",
                    "Subtype": "Short_Tandem_Repeat",
                    "Start": i+1,
                    "End": i + reps*unit + remainder,
                    "Length": reps*unit + remainder,
                    "sequence": wrap(motif_seq),
                    "Stability Score": reps,
                    "Significance": "high" if reps > 8 else "medium",
                    "Conservation score": cons_score,
                    "Conservation p-value": pval,
                    "Conservation significance": cons_sig,
                    "Repeat/arm/Contributing region": repeat_unit,
                    "Spacer": None
                })
                i = i + reps*unit + remainder - 1
                found = True
                break
        if not found:
            i += 1
    return results

def find_multimeric_gquadruplex(seq, sequence_name=""):
    pattern = r"(G{3,}\w{1,12}){4,}"
    results = []
    for m in overlapping_finditer(pattern, seq):
        motif_seq = m.group(0)
        score = g4hunter_score(motif_seq)
        if score >= 1.0:
            cons_score, pval, cons_sig = scientific_conservation_score(seq, motif_seq)
            results.append({
                "SequenceName": sequence_name,
                "Class": "G4",
                "Subtype": "Multimeric_G4",
                "Start": m.start()+1,
                "End": m.end(),
                "Length": len(motif_seq),
                "sequence": wrap(motif_seq),
                "Stability Score": round(score*1.2,2),
                "Significance": "high" if score > 1.2 else "medium",
                "Conservation score": cons_score,
                "Conservation p-value": pval,
                "Conservation significance": cons_sig,
                "Repeat/arm/Contributing region": None,
                "Spacer": None
            })
    return results






def merge_hotspots(hotspots):
    if not hotspots: return []
    merged = [hotspots[0]]
    for current in hotspots[1:]:
        last = merged[-1]
        if current['Start'] <= last['End']:
            last['End'] = max(last['End'], current['End'])
            last['Length'] = last['End'] - last['Start'] + 1
            last['MotifCount'] += current['MotifCount']
            last['TypeDiversity'] = max(last['TypeDiversity'], current['TypeDiversity'])
            last['Stability Score'] += current['Stability Score']
            last['Conservation score'] = max(last['Conservation score'], current['Conservation score'])
        else:
            merged.append(current)
    return merged

# ... (Add all other motif finders here, updating their outputs to contain conservation scores as above) ...

def all_motifs(seq, nonoverlap=False, report_hotspots=False):
    if not seq or not re.match("^[ATGC]+$", seq, re.IGNORECASE):
        return []
    seq = seq.upper()
    motif_list = (
        find_curved_DNA(seq) +
        # ... add all other motif finders here ...
    )
    motif_list = [m for m in motif_list if validate_motif(m, len(seq))]
    motif_list += find_hybrids(motif_list, seq)
    if report_hotspots:
        motif_list += find_hotspots(motif_list, len(seq))
    return motif_list
