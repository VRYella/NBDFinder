import re
import numpy as np
from collections import defaultdict, Counter
import random

# ========== UTILITY FUNCTIONS ==========
def parse_fasta(fasta_str):
    return "".join([line.strip() for line in fasta_str.split('\n') 
                   if not line.startswith(">")]).upper().replace(" ", "").replace("U", "T")

def wrap(seq, width=60):
    return "\n".join([seq[i:i+width] for i in range(0, len(seq), width)])

def gc_content(seq):
    return 100 * (seq.count('G') + seq.count('C')) / max(1, len(seq))

def reverse_complement(seq):
    return seq.translate(str.maketrans("ATGC", "TACG"))[::-1]

def is_palindrome(seq):
    return seq == reverse_complement(seq)

def calculate_tm(seq):
    if len(seq) < 14:
        return 2*(seq.count('A')+seq.count('T')) + 4*(seq.count('G')+seq.count('C'))
    else:
        return 64.9 + 41*(seq.count('G')+seq.count('C')-16.4)/len(seq)

def shuffle_sequence(seq):
    return ''.join(random.sample(seq, len(seq)))

def g4hunter_score(seq):
    scores = [3 if b=='G' else -1 if b=='C' else 0 for b in seq.upper()]
    return sum(scores)/len(scores) if scores else 0

# ========== MOTIF DETECTION FUNCTIONS ==========
def overlapping_finditer(pattern, seq):
    for m in re.compile(pattern, re.IGNORECASE).finditer(seq): 
        yield m

def validate_motif(motif, seq_len):
    return (0 < motif['Start'] <= motif['End'] <= seq_len and 
            motif['Length'] == (motif['End'] - motif['Start'] + 1) and 
            re.match("^[ATGC]+$", motif['Sequence']))

def all_motifs(seq):
    """Main function to detect all non-B DNA motifs in a sequence"""
    if not seq or not re.match("^[ATGC]+$", seq, re.IGNORECASE):
        return []
    
    seq = seq.upper()
    motifs = []
    detection_functions = [
        find_curved_dna, find_zdna, find_slipped_dna, find_rlfs,
        find_cruciform, find_triplex, find_gtriplex,
        find_gquadruplex, find_imotif, find_hybrids
    ]
    
    for finder in detection_functions:
        try:
            motifs.extend(finder(seq))
        except Exception as e:
            print(f"Warning: {finder.__name__} failed: {str(e)}")
            continue
            
    return [m for m in motifs if validate_motif(m, len(seq))]

# ========== INDIVIDUAL MOTIF DETECTORS ==========
def find_curved_dna(seq):
    results = []
    # APR detection (nBSTAPR logic)
    for m in overlapping_finditer(r"(?=((A{3,11}[ATGC]{2,5}){3,}))", seq):
        apr_seq = m.group(1)
        a_count = apr_seq.count('A') + apr_seq.count('T')
        score = min(1.0, 0.2 * len(re.findall(r"A{3,11}", apr_seq)) + (a_count/len(apr_seq)))
        results.append({
            "Class": "Curved_DNA",
            "Subtype": "A-Phased_Repeat",
            "Start": m.start()+1,
            "End": m.start()+len(apr_seq),
            "Length": len(apr_seq),
            "Sequence": wrap(apr_seq),
            "ScoreMethod": "Brukner_Curvature",
            "Score": f"{score:.2f}"
        })
    
    # Local curvature (non-overlapping)
    for m in overlapping_finditer(r"(?=(A{7}|T{7}))", seq):
        if not any(m.start() < r['End'] and m.end() > r['Start'] for r in results):
            results.append({
                "Class": "Curved_DNA",
                "Subtype": "Local_Curvature",
                "Start": m.start()+1,
                "End": m.end(),
                "Length": len(m.group(1)),
                "Sequence": wrap(m.group(1)),
                "ScoreMethod": "A/T_Run",
                "Score": "1.00"
            })
    return results

# [Rest of your individual motif detection functions...]
# (Keep the implementations of find_zdna, find_slipped_dna, etc. as they were,
#  just ensure consistent formatting and fix any syntax errors like unmatched parentheses)

def find_hotspots(motifs, seq_len, window=100, min_count=3):
    """Identify genomic regions with high motif density"""
    if not motifs or seq_len <= 0:
        return []
        
    hotspots = []
    positions = [(m['Start'], m['End']) for m in motifs]
    
    for i in range(0, seq_len - window + 1):
        region_start = i + 1
        region_end = i + window
        count = sum(s <= region_end and e >= region_start for s, e in positions)
        
        if count >= min_count:
            motifs_in_region = [m for m in motifs 
                              if m['Start'] <= region_end and m['End'] >= region_start]
            diversity = len({m['Subtype'] for m in motifs_in_region})
            score = min(1.0, count/10 + diversity/5)
            
            hotspots.append({
                "RegionStart": region_start,
                "RegionEnd": region_end,
                "MotifCount": count,
                "TypeDiversity": diversity,
                "Score": f"{score:.2f}",
                "Coverage": f"{100*sum(m['End']-m['Start']+1 for m in motifs_in_region)/window:.1f}%"
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
            cov = (int(last['Coverage'][:-1])*window/100 + 
                   sum(m['End']-m['Start']+1 for m in motifs 
                       if m['Start'] >= last['RegionStart'] and m['End'] <= h['RegionEnd']))
            last['Coverage'] = f"{100*cov/(h['RegionEnd']-last['RegionStart']+1):.1f}%"
        else:
            merged.append(h)
            
    return merged
