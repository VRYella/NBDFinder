"""
Category 9: Hybrid Structures Detection Module
==============================================

This module implements detection algorithms for hybrid structures
formed by overlapping non-B DNA motifs.

Scientific Basis:
- Hybrid structures form when multiple non-B DNA motifs overlap
- Can create complex regulatory elements and hotspots
- Important for understanding combinatorial effects

Author: Dr. Venkata Rajesh Yella
Updated: 2024
"""

from .shared_utils import wrap

def find_hybrids(motifs, seq):
    """
    Find hybrid structures formed by overlapping motifs.
    
    Enhanced Algorithm:
    1. Event-based overlap detection using sorted start/end events
    2. Multi-motif overlap identification with class diversity requirements
    3. Combined scoring from contributing motifs
    
    Parameters:
    motifs (list): List of detected motifs
    seq (str): Original DNA sequence
    
    Returns:
    list: List of hybrid structure motif dictionaries
    """
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