import re
import numpy as np
from collections import defaultdict, Counter
import random

# ========== UTILS FUNCTIONS ==========
def parse_fasta(fasta_str): return "".join([line.strip() for line in fasta_str.split('\n') if not line.startswith(">")]).upper().replace(" ", "").replace("U", "T")
def wrap(seq, width=60): return "\n".join([seq[i:i+width] for i in range(0, len(seq), width)])
def gc_content(seq): return 100 * (seq.count('G') + seq.count('C')) / max(1, len(seq))
def reverse_complement(seq): return seq.translate(str.maketrans("ATGC", "TACG"))[::-1]
def is_palindrome(seq): return seq == reverse_complement(seq)
def calculate_tm(seq): return 2*(seq.count('A')+seq.count('T')) + 4*(seq.count('G')+seq.count('C'))) if len(seq) < 14 else 64.9 + 41*(seq.count('G')+seq.count('C')-16.4)/len(seq)
def shuffle_sequence(seq): return ''.join(random.sample(seq, len(seq)))
def g4hunter_score(seq):
    scores = [3 if b=='G' else -1 if b=='C' else 0 for b in seq.upper()]
    return sum(scores)/len(scores) if scores else 0

# ========== MOTIF DETECTION ==========
def overlapping_finditer(pattern, seq):
    for m in re.compile(pattern, re.IGNORECASE).finditer(seq): yield m

def all_motifs(seq):
    if not seq or not re.match("^[ATGC]+$", seq, re.IGNORECASE): return []
    seq = seq.upper(); motifs = []
    for finder in [find_curved_dna, find_zdna, find_slipped_dna, find_rlfs, find_cruciform, 
                  find_triplex, find_gtriplex, find_gquadruplex, find_imotif, find_hybrids]:
        motifs.extend(finder(seq))
    return [m for m in motifs if validate_motif(m, len(seq))]

def validate_motif(motif, seq_len):
    return (0 < motif['Start'] <= motif['End'] <= seq_len and 
            motif['Length'] == (motif['End'] - motif['Start'] + 1) and 
            re.match("^[ATGC]+$", motif['Sequence']))

# 1. Curved DNA (APRs and local curvature)
def find_curved_dna(seq):
    results = []
    # APR detection (nBSTAPR logic)
    for m in overlapping_finditer(r"(?=((A{3,11}[ATGC]{2,5}){3,}))", seq):
        apr_seq = m.group(1); a_count = apr_seq.count('A') + apr_seq.count('T')
        score = min(1.0, 0.2 * len(re.findall(r"A{3,11}", apr_seq)) + (a_count/len(apr_seq)))
        results.append({
            "Class": "Curved_DNA", "Subtype": "A-Phased_Repeat", "Start": m.start()+1,
            "End": m.start()+len(apr_seq), "Length": len(apr_seq), "Sequence": wrap(apr_seq),
            "ScoreMethod": "Brukner_Curvature", "Score": f"{score:.2f}"
        })
    # Local curvature (non-overlapping)
    for m in overlapping_finditer(r"(?=(A{7}|T{7}))", seq):
        if not any(m.start() < r['End'] and m.end() > r['Start'] for r in results):
            results.append({
                "Class": "Curved_DNA", "Subtype": "Local_Curvature", "Start": m.start()+1,
                "End": m.end(), "Length": len(m.group(1)), "Sequence": wrap(m.group(1)),
                "ScoreMethod": "A/T_Run", "Score": "1.00"
            })
    return results

# 2. Z-DNA (Z-Seeker scoring)
def find_zdna(seq):
    results = []
    for m in overlapping_finditer(r"(?=(([GC][GC]){6,}))", seq):
        seq_frag = m.group(1); cg_count = seq_frag.count('CG') + seq_frag.count('GC')
        score = min(1.0, (len(seq_frag)/20 + cg_count/10))
        results.append({
            "Class": "Z-DNA", "Subtype": "Z-Forming_Repeat", "Start": m.start()+1,
            "End": m.start()+len(seq_frag), "Length": len(seq_frag), "Sequence": wrap(seq_frag),
            "ScoreMethod": "Z-Seeker", "Score": f"{score:.2f}"
        })
    return results

# 3. Slipped DNA (nBST logic)
def find_slipped_dna(seq):
    return [{
        "Class": "Slipped_DNA", "Subtype": "Direct_Repeat", "Start": m.start()+1,
        "End": m.start()+len(m.group(1))), "Length": len(m.group(1))), 
        "Sequence": wrap(m.group(1)), "ScoreMethod": "nBST",
        "Score": f"{min(1.0, len(m.group(1))/50):.2f}"
    } for m in overlapping_finditer(r"(?=(([ATGC]{10,300})(.{0,100}?)\2))", seq)]

# 4. R-Loops (Improved QmRLFS)
def find_rlfs(seq):
    if len(seq) < 100: return []
    results = []
    for m in re.finditer(r"(G{3,}[ATGC]{1,10}G{3,}[ATGC]{1,10}G{4,})", seq, re.IGNORECASE):
        riz_seq = m.group(1); window = seq[m.end():m.end()+300]
        if gc_content(riz_seq + window) > 50:
            score = min(1.0, 0.6*gc_content(riz_seq+window)/100 + 0.4*(len(re.findall(r"G{3,}", riz_seq+window))/5))
            results.append({
                "Class": "R-Loop", "Subtype": "RLFS", "Start": m.start()+1,
                "End": m.end()+len(window), "Length": len(riz_seq)+len(window),
                "Sequence": wrap(riz_seq+window), "ScoreMethod": "QmRLFS_Thermo",
                "Score": f"{score:.2f}"
            })
    return results

# 5. Cruciform/Hairpin
def find_cruciform(seq):
    return [{
        "Class": "Cruciform", "Subtype": "Inverted_Repeat", "Start": m.start()+1,
        "End": m.start()+(2*len(m.group(2))), "Length": 2*len(m.group(2))),
        "Sequence": wrap(m.group(1)), "ScoreMethod": "nBSTAPR",
        "Score": f"{min(1.0, len(m.group(2))/15 + (m.group(2).count('A')+m.group(2).count('T'))/len(m.group(2))*0.3):.2f}"
    } for m in overlapping_finditer(r"(?=(([ATGC]{6,})\s*?\2))", seq) if len(m.group(2)) >= 6]

# 6. Triplex (H-DNA) and Sticky DNA
def find_triplex(seq):
    results = [{
        "Class": "Triplex_DNA", "Subtype": "H-DNA", "Start": m.start()+1,
        "End": m.start()+len(m.group(1))), "Length": len(m.group(1))),
        "Sequence": wrap(m.group(1)), "ScoreMethod": "nBST",
        "Score": f"{min(1.0, len(m.group(2))/15):.2f}"
    } for m in overlapping_finditer(r"(?=(([CT]{10,})(.{0,100}?)\2))", seq)]
    # Sticky DNA (GAA/TTC repeats)
    results.extend([{
        "Class": "Triplex_DNA", "Subtype": "Sticky_DNA", "Start": m.start()+1,
        "End": m.end(), "Length": len(m.group(1))), "Sequence": wrap(m.group(1))),
        "ScoreMethod": "GAA_Repeat", "Score": f"{min(1.0, len(m.group(1))/30):.2f}"
    } for m in overlapping_finditer(r"(?=((GAA|TTC){3,}))", seq)])
    return results

# 7. G-Triplex
def find_gtriplex(seq):
    return [{
        "Class": "G-Triplex", "Subtype": "Three_G-Runs", "Start": m.start()+1,
        "End": m.start()+len(m.group(1))), "Length": len(m.group(1))),
        "Sequence": wrap(m.group(1)), "ScoreMethod": "G3_Stability",
        "Score": f"{min(1.0, len(re.findall(r"G{3,}", m.group(1)))/3):.2f}"
    } for m in overlapping_finditer(r"(?=(G{3,}\w{1,7}G{3,}\w{1,7}G{3,}))", seq)]

# 8. G-Quadruplex variants
def find_gquadruplex(seq):
    results = []
    # Multimeric > Bipartite > Canonical > Relaxed > Bulged
    for pattern, subtype in [
        (r"(G{3,}\w{1,12}){4,}", "Multimeric_G4"),
        (r"(G{3,}\w{1,30}){6,}", "Bipartite_G4"),
        (r"G{3,}\w{1,7}G{3,}\w{1,7}G{3,}\w{1,7}G{3,}", "Canonical_G4"),
        (r"G{3,}\w{1,12}G{3,}\w{1,12}G{3,}\w{1,12}G{3,}", "Relaxed_G4"),
        (r"G{3,}(\w{0,3}G{3,}){3,}", "Bulged_G4")
    ]:
        for m in overlapping_finditer(pattern, seq):
            score = g4hunter_score(m.group(1)); threshold = 0.8 if "Relaxed" in subtype else 1.0 if "Canonical" in subtype else 0.9
            if score >= threshold and not any(m.start() < r['End'] and m.end() > r['Start'] for r in results):
                results.append({
                    "Class": "G4", "Subtype": subtype, "Start": m.start()+1,
                    "End": m.end(), "Length": len(m.group(1))), "Sequence": wrap(m.group(1))),
                    "ScoreMethod": "G4Hunter", "Score": f"{score:.2f}"
                })
    return results

# 9. i-Motif
def find_imotif(seq):
    return [{
        "Class": "i-Motif", "Subtype": "C_Quadruplex", "Start": m.start()+1,
        "End": m.start()+len(m.group(1))), "Length": len(m.group(1))),
        "Sequence": wrap(m.group(1)), "ScoreMethod": "iM_Stability",
        "Score": f"{min(1.0, len(re.findall(r"C{3,}", m.group(1)))/4):.2f}"
    } for m in overlapping_finditer(r"(?=(C{3,}\w{0,7}C{3,}\w{0,7}C{3,}\w{0,7}C{3,}))", seq)]

# 10. Hybrid Motifs
def find_hybrids(seq):
    hybrids = []
    g4s = [(m.start()+1, m.end(), m.group(1)) for m in overlapping_finditer(r"(G{3,}\w{1,12}G{3,}\w{1,12}G{3,}\w{1,12}G{3,})", seq)]
    ims = [(m.start()+1, m.start()+len(m.group(1)), m.group(1)) for m in overlapping_finditer(r"(C{3,}\w{0,7}C{3,}\w{0,7}C{3,}\w{0,7}C{3,})", seq)]
    for g_start, g_end, g_seq in g4s:
        for c_start, c_end, c_seq in ims:
            if (g_start <= c_end) and (c_start <= g_end):
                overlap_seq = seq[min(g_start,c_start)-1:max(g_end,c_end)]
                score = min(1.0, (g4hunter_score(g_seq) + len(re.findall(r"C{3,}", c_seq))/4)/2 *1.1)
                hybrids.append({
                    "Class": "Hybrid", "Subtype": "G4_iM", "Start": min(g_start,c_start),
                    "End": max(g_end,c_end), "Length": max(g_end,c_end)-min(g_start,c_start)+1,
                    "Sequence": wrap(overlap_seq), "ScoreMethod": "Hybrid_Score",
                    "Score": f"{score:.2f}"
                })
    return hybrids

# 11. Hotspot Detection
def find_hotspots(motifs, seq_len, window=100, min_count=3):
    hotspots = []; pos = [(m['Start'], m['End']) for m in motifs]
    for i in range(0, seq_len - window + 1):
        region = (i+1, i+window)
        count = sum(s <= region[1] and e >= region[0] for s,e in pos)
        if count >= min_count:
            motifs_in_region = [m for m in motifs if m['Start'] <= region[1] and m['End'] >= region[0]]
            diversity = len({m['Subtype'] for m in motifs_in_region})
            score = min(1.0, count/10 + diversity/5)
            hotspots.append({
                "RegionStart": region[0], "RegionEnd": region[1],
                "MotifCount": count, "TypeDiversity": diversity,
                "Score": f"{score:.2f}", "Coverage": f"{100*sum(m['End']-m['Start']+1 for m in motifs_in_region)/window:.1f}%"
            })
    # Merge overlapping hotspots
    merged = []
    for h in sorted(hotspots, key=lambda x: x['RegionStart']):
        if merged and h['RegionStart'] <= merged[-1]['RegionEnd']:
            last = merged[-1]
            last['RegionEnd'] = max(last['RegionEnd'], h['RegionEnd'])
            last['MotifCount'] += h['MotifCount']
            last['TypeDiversity'] = max(last['TypeDiversity'], h['TypeDiversity'])
            last['Score'] = f"{min(1.0, float(last['Score']) + float(h['Score'])):.2f}"
            last['Coverage'] = f"{100*(int(last['Coverage'][:-1])*window/100 + sum(m['End']-m['Start']+1 for m in motifs if m['Start'] >= last['RegionStart'] and m['End'] <= h['RegionEnd']))/(h['RegionEnd']-last['RegionStart']+1):.1f}%"
        else: merged.append(h)
    return merged
