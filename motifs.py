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
    if not seq or not re.match("^[ATGC]+$", seq, re.IGNORECASE): return []
    seq = seq.upper()
    results = (find_apr(seq) + find_zdna(seq) + find_slipped_dna(seq) + find_rlfs(seq) + 
               find_cruciform(seq) + find_hdna(seq) + find_gtriplex(seq) + find_gquadruplex(seq) + 
               find_relaxed_gquadruplex(seq) + find_bulged_gquadruplex(seq) + 
               find_bipartite_gquadruplex(seq) + find_multimeric_gquadruplex(seq) + 
               find_imotif(seq) + find_hybrids(seq) + find_sticky_dna(seq))
    return [m for m in results if validate_motif(m, len(seq))]

def validate_motif(motif, seq_len):
    return (0 < motif['Start'] <= motif['End'] <= seq_len and 
            motif['Length'] == (motif['End'] - motif['Start'] + 1) and 
            re.match("^[ATGC]+$", motif['Sequence']))

def find_apr(seq):
    """Detect A-Phased Repeats with curvature scoring"""
    results = []
    apr_regions = []
    for m in overlapping_finditer(r"(?=((A{3,6}[ATGC]{2,5}){3,}))", seq):
        seq_frag = m.group(1)
        a_count = seq_frag.count('A') + seq_frag.count('T')
        score = min(1.0, 0.2 * len(re.findall(r"A{3,6}", seq_frag)) + (a_count / len(seq_frag)))
        motif = {
            "Class": "Curved_DNA", "Subtype": "A-Phased_Repeat", "Start": m.start()+1,
            "End": m.start()+len(seq_frag), "Length": len(seq_frag), "Sequence": wrap(seq_frag),
            "ScoreMethod": "Brukner_Curvature", "Score": f"{score:.2f}"
        }
        results.append(motif)
        apr_regions.append((motif['Start'], motif['End']))
    
    # Append non-overlapping local curved regions
    results += find_local_curved(seq, apr_regions)
    return results

def find_local_curved(seq, apr_regions):
    """Detect local curved DNA (A7/T7) not overlapping with APR"""
    results = []
    for m in overlapping_finditer(r"(A{7,}|T{7,})", seq):
        start, end = m.start()+1, m.end()
        if not any(s <= start <= e or s <= end <= e for s, e in apr_regions):
            results.append({
                "Class": "Curved_DNA", "Subtype": "Local_Curved", "Start": start,
                "End": end, "Length": len(m.group()), "Sequence": wrap(m.group()),
                "ScoreMethod": "A/T_Stretch", "Score": f"{min(1.0, len(m.group())/10):.2f}"})
    return results


def find_zdna(seq):
    results = []
    for m in overlapping_finditer(r"(?=(([AT][CG]){6,}))", seq):
        seq_frag = m.group(1)
        cg = seq_frag.count('CG') + seq_frag.count('GC')
        score = min(1.0, (cg/5) + (len(seq_frag)/15))
        results.append({
            "Class": "Z-DNA", "Subtype": "Z-Forming_Repeat", "Start": m.start()+1,
            "End": m.start()+len(seq_frag), "Length": len(seq_frag), "Sequence": wrap(seq_frag),
            "ScoreMethod": "Z-Seeker", "Score": f"{score:.2f}"})
    return results

def find_slipped_dna(seq):
    pattern = r"(?=(([ATGC]{10,300})\1))"
    return [{
        "Class": "Slipped_DNA", "Subtype": "Direct_Repeat", "Start": m.start()+1,
        "End": m.start()+len(m.group(1))*2, "Length": len(m.group(1))*2,
        "Sequence": wrap(m.group(1)+m.group(1)), "ScoreMethod": "nBST_DR",
        "Score": f"{min(1.0, len(m.group(1))/300):.2f}"} for m in overlapping_finditer(pattern, seq)]

def find_rlfs(seq):
    if len(seq) < 100: return []
    results = []
    for m in re.finditer(r"(G{3,}[ATGC]{1,10}G{3,}[ATGC]{1,10}G{4,})", seq, re.IGNORECASE):
        riz_seq = m.group(1)
        if gc_content(riz_seq) < 50: continue
        for rez in find_rez(seq, m.end()):
            stability = min(1.0, 0.6 * (gc_content(riz_seq + rez['seq']) / 100) + 
                            0.4 * (len(re.findall(r"G{3,}", riz_seq + rez['seq'])) / 5))
            results.append({
                "Class": "R-Loop", "Subtype": "RLFS", "Start": m.start()+1,
                "End": m.start()+len(riz_seq)+rez['end'], "Length": len(riz_seq)+rez['end'],
                "Sequence": wrap(riz_seq+rez['seq']), "ScoreMethod": "QmRLFS_Thermo",
                "Score": f"{stability:.2f}"})
    return results

def find_rez(seq, start_pos, max_len=500):
    window = seq[start_pos:start_pos+max_len]
    return [{'seq': window, 'end': len(window)}] if gc_content(window) >= 40 else []

def find_cruciform(seq):
    """Detect cruciform-forming inverted repeats (non-palindromic with spacer)"""
    results = []
    pattern = r"(?=([ATGC]{6,})([ATGC]{0,10})\2[::-1])"  # Pseudocode pattern

    # Better: Use a loop to find inverted repeats manually
    for i in range(len(seq) - 12):
        for arm_len in range(6, 15):  # arms of 6 to 14
            arm = seq[i:i+arm_len]
            spacer = seq[i+arm_len:i+arm_len+5]
            rev = reverse_complement(arm)
            end = i + arm_len + 5 + arm_len
            if end <= len(seq) and seq[i+arm_len+5:end] == rev:
                full_seq = seq[i:end]
                score = min(1.0, (arm_len / 15) + ((arm.count('A') + arm.count('T')) / arm_len * 0.3))
                results.append({
                    "Class": "Cruciform", "Subtype": "Inverted_Repeat", "Start": i+1,
                    "End": end, "Length": end - i, "Sequence": wrap(full_seq),
                    "ScoreMethod": "nBST_IR", "Score": f"{score:.2f}"
                })
    return results

def find_hdna(seq):
    return [{
        "Class": "Triplex_DNA", "Subtype": "H-DNA_Pyrimidine", "Start": m.start()+1,
        "End": m.start()+len(m.group(1)), "Length": len(m.group(1)), "Sequence": wrap(m.group(1)),
        "ScoreMethod": "Triplex_Propensity", "Score": f"{min(1.0, len(m.group(2))/10 + 0.3):.2f}"}
        for m in overlapping_finditer(r"(?=(([CT]{5,})\s*?\2))", seq) if len(m.group(2)) >= 5]

def find_sticky_dna(seq):
    return [{
        "Class": "Sticky_DNA", "Subtype": "GAA_TTC_Repeat", "Start": m.start()+1,
        "End": m.end(), "Length": len(m.group()), "Sequence": wrap(m.group()),
        "ScoreMethod": "Potaman_Score", "Score": f"{min(1.0, len(m.group())/30):.2f}"}
        for m in overlapping_finditer(r"(GAA){3,}|(TTC){3,}", seq)]

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
