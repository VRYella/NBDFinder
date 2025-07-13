import re
import numpy as np
from typing import List, Dict, Iterator
from utils import (wrap, gc_content, reverse_complement, 
                  is_palindrome, calculate_tm, kmer_conservation,
                  motif_conservation)

# Constants
MIN_CRUCIFORM_SCORE = 0.7
MIN_G4_SCORE = 1.0
MIN_IMOTIF_SCORE = 0.7
MIN_ZDNA_SCORE = 0.8

def overlapping_finditer(pattern: str, seq: str) -> Iterator[re.Match]:
    """Find overlapping matches"""
    regex = re.compile(pattern, re.IGNORECASE)
    for m in regex.finditer(seq):
        yield m

def all_motifs(seq: str) -> List[Dict]:
    """Main detection function with non-overlapping priority"""
    if not seq or not re.match("^[ATGC]+$", seq, re.IGNORECASE):
        return []
    
    seq = seq.upper()
    conservation_scores = kmer_conservation(seq)
    
    # Detection priority order
    detectors = [
        find_gquadruplex,
        find_cruciform,
        find_imotif,
        find_zdna,
        find_hdna,
        find_rlfs,
        find_apr,
        find_bent_dna,
        find_slipped_dna
    ]
    
    # Track covered positions
    covered = np.zeros(len(seq), dtype=bool)
    results = []
    
    for detector in detectors:
        for motif in detector(seq):
            start, end = motif['Start']-1, motif['End']-1
            
            # Skip overlaps
            if np.any(covered[start:end+1]):
                continue
                
            # Add conservation score
            motif['Conservation'] = motif_conservation(motif['Sequence'], conservation_scores)
            results.append(motif)
            covered[start:end+1] = True
    
    return [m for m in results if validate_motif(m, len(seq))]

def validate_motif(motif: Dict, seq_len: int) -> bool:
    """Validate motif structure"""
    return (0 < motif['Start'] <= motif['End'] <= seq_len and
            motif['Length'] == (motif['End'] - motif['Start'] + 1) and
            re.match("^[ATGC]+$", motif['Sequence']))

# --- Motif Detection Functions ---

def find_cruciform(seq: str) -> List[Dict]:
    """Detect cruciform structures with strict requirements"""
    pattern = r"(?=(([ATGC]{6,20})([ATGC]{0,4}?)([ATGC]{6,20})))"
    results = []
    
    for m in overlapping_finditer(pattern, seq):
        stem1, loop, stem2 = m.group(2), m.group(3), m.group(4)
        
        if reverse_complement(stem1) != stem2:
            continue
            
        stem_len = len(stem1)
        total_len = len(m.group(1))
        at_content = (stem1.count('A') + stem1.count('T') + stem2.count('A') + stem2.count('T')) / (2*stem_len)
        loop_penalty = min(1.0, len(loop)/4)
        score = min(1.0, (stem_len/20) + (at_content*0.4) - (loop_penalty*0.2))
        
        if score >= MIN_CRUCIFORM_SCORE:
            results.append({
                "Class": "Cruciform",
                "Subtype": "Inverted_Repeat",
                "Start": m.start() + 1,
                "End": m.start() + total_len,
                "Length": total_len,
                "Sequence": wrap(m.group(1)),
                "ScoreMethod": "Cruciform_Stability_v2",
                "Score": f"{score:.2f}"
            })
    return results

def find_gquadruplex(seq: str) -> List[Dict]:
    """Detect G-quadruplexes with strict requirements"""
    pattern = r"(G{3,5})([ATGC]{1,7})(G{3,5})([ATGC]{1,7})(G{3,5})([ATGC]{1,7})(G{3,5})"
    results = []
    
    for m in re.finditer(pattern, seq):
        g_runs = [m.group(i) for i in range(1, 8, 2)]
        loops = [m.group(i) for i in range(2, 7, 2)]
        
        g_score = sum(min(4, len(r)) for r in g_runs)/12
        loop_penalty = sum(1/(len(l)+1) for l in loops)/3
        score = min(1.0, 0.7*g_score + 0.3*(1-loop_penalty))
        
        if score >= MIN_G4_SCORE:
            full_seq = ''.join(g_runs + loops[:3])
            results.append({
                "Class": "G4",
                "Subtype": "Canonical_G4",
                "Start": m.start() + 1,
                "End": m.start() + len(full_seq),
                "Length": len(full_seq),
                "Sequence": wrap(full_seq),
                "ScoreMethod": "G4Hunter_Strict",
                "Score": f"{score:.2f}"
            })
    return results

def find_imotif(seq: str) -> List[Dict]:
    """Detect i-Motifs with C-rich requirements"""
    pattern = r"(?=(C{3,7})([ATGC]{0,7})(C{3,7})([ATGC]{0,7})(C{3,7})([ATGC]{0,7})(C{3,7}))"
    results = []
    
    for m in overlapping_finditer(pattern, seq):
        c_runs = [m.group(i) for i in range(1, 8, 2)]
        loops = [m.group(i) for i in range(2, 7, 2)]
        
        c_score = sum(min(4, len(r)) for r in c_runs)/12
        loop_penalty = sum(1/(len(l)+1) for l in loops)/3
        ph = 5.5  # Assumed pH for i-Motif stability
        score = min(1.0, (0.6*c_score + 0.2*(1-loop_penalty) + 0.2*(7.0 - ph)/2))
        
        if score >= MIN_IMOTIF_SCORE:
            full_seq = ''.join(c_runs + loops[:3])
            results.append({
                "Class": "i-Motif",
                "Subtype": "C_Quadruplex",
                "Start": m.start() + 1,
                "End": m.start() + len(full_seq),
                "Length": len(full_seq),
                "Sequence": wrap(full_seq),
                "ScoreMethod": "iM_Stability_pH",
                "Score": f"{score:.2f}"
            })
    return results

# [Additional motif detection functions with similar improvements...]

def find_hotspots(seq: str, motifs: List[Dict], window: int = 100, min_count: int = 3) -> List[Dict]:
    """Identify enriched motif regions"""
    positions = [(m['Start'], m['End']) for m in motifs]
    hotspots = []
    
    for i in range(len(seq) - window + 1):
        start, end = i+1, i+window
        count = sum(s <= end and e >= start for s,e in positions)
        
        if count >= min_count:
            region_motifs = [m for m in motifs if m['Start'] <= end and m['End'] >= start]
            diversity = len({m['Subtype'] for m in region_motifs})
            cons_score = np.mean([m.get('Conservation', 0) for m in region_motifs])
            
            hotspots.append({
                "RegionStart": start,
                "RegionEnd": end,
                "MotifCount": count,
                "TypeDiversity": diversity,
                "AvgConservation": f"{cons_score:.2f}",
                "Score": f"{min(1.0, count/10 + diversity/5 + cons_score/3):.2f}"
            })
    
    return hotspots
