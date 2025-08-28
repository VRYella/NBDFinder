"""
Category 9: Hybrid Structures Detection Module
==============================================

This module implements detection algorithms for hybrid structures
formed by overlapping non-B DNA motifs.

Scientific Basis:
-----------------
- Hybrid structures form when multiple non-B DNA motifs overlap
- Overlapping motifs can generate higher-order DNA conformations
- These hybrids often act as complex regulatory elements and hotspots
- They are implicated in chromatin remodeling, nucleosome exclusion,
  transcription regulation, replication stalling, and genome instability

Biological Insight:
-------------------
- Individual motifs (e.g., G-quadruplex, Z-DNA, A-tracts) are important,
  but when they overlap, the combined structure may have emergent properties
- For example, a G-quadruplex overlapping with Z-DNA can stabilize unusual
  supercoiling domains, or an A-tract plus Mirror-repeat may alter local stiffness
- Detecting hybrids helps identify *regulatory hotspots* and *fragile regions*
  where non-B DNA contributes to disease, mutation, or gene expression regulation

Author: Dr. Venkata Rajesh Yella
Updated: 2024 (annotated for biological interpretation)
"""

from .shared_utils import wrap

def find_hybrids(motifs, seq):
    """
    Find hybrid structures with detailed subclass breakdown analysis.
    
    Algorithmic Approach:
    ---------------------
    1. Event-based overlap detection (start/end tracking of motifs)
    2. Identify regions where ≥2 different motif *classes* overlap
    3. Build breakdown by Class (e.g., G-quadruplex, Z-DNA) and Subtype
       (e.g., Canonical PQS, Left-handed Z-DNA, Inverted Mirror)
    4. Score hybrid based on contributing motifs, overlap strength,
       and potential interaction bonus
    5. Report each hybrid with sequence fragment and contributing motifs
    
    Biological Interpretation:
    --------------------------
    - Hybrids detected here represent **combinatorial non-B DNA zones**
    - Overlap degree: how tightly structures coexist in the same DNA window
    - Interaction strength: likelihood that overlapping motifs *cooperate*
      (e.g., one stabilizing the other, or creating local rigidity)
    - High-scoring hybrids = *regulatory hotspots* where nucleosome
      organization, transcription factor binding, or polymerase stalling
      may be strongly influenced.
    
    Parameters:
    -----------
    motifs (list): List of detected motifs from upstream analysis
    seq (str): Original DNA sequence
    
    Returns:
    --------
    list: List of hybrid structure dictionaries with rich annotation
    """
    
    # Build a timeline of motif "start" and "end" events
    events = []
    for idx, m in enumerate(motifs):
        events.append((m['Start'], 'start', idx))
        events.append((m['End'] + 1, 'end', idx))
    events.sort()
    
    active = set()         # currently overlapping motifs
    region_start = None
    results = []
    
    # Walk through events in order
    for pos, typ, idx in events:
        if typ == 'start':
            active.add(idx)
            
            # A hybrid region begins when we first hit ≥2 overlaps
            if len(active) == 2:
                region_start = pos
        
        elif typ == 'end':
            if len(active) >= 2:
                region_end = pos - 1
                involved_idxs = list(active)
                
                # Collect motifs in this overlapping region
                involved_motifs = [motifs[i] for i in involved_idxs]
                subclass_breakdown = {}
                class_breakdown = {}
                
                # Count by class and subclass
                for motif in involved_motifs:
                    motif_class = motif.get('Class', 'Unknown')
                    motif_subtype = motif.get('Subtype', 'Unknown')
                    
                    class_breakdown[motif_class] = class_breakdown.get(motif_class, 0) + 1
                    subclass_key = f"{motif_class}:{motif_subtype}"
                    subclass_breakdown[subclass_key] = subclass_breakdown.get(subclass_key, 0) + 1
                
                # Only proceed if ≥2 different motif classes overlap
                if len(class_breakdown) >= 2:
                    total_length = region_end - region_start + 1
                    
                    # Overlap degree = hybrid size relative to largest motif
                    overlap_degree = total_length / max(m['Length'] for m in involved_motifs)
                    
                    # Base score = sum of contributing motif scores
                    base_score = sum(float(m.get("Score", 0.0)) for m in involved_motifs)
                    
                    # Interaction bonus = reward for multiple overlaps
                    # Biological analogy: more motifs = more structural stress/rigidity
                    interaction_bonus = len(involved_motifs) * 0.2
                    final_score = base_score * (1 + interaction_bonus) * overlap_degree
                    
                    # Sort for clean labeling
                    sorted_classes = sorted(class_breakdown.keys())
                    sorted_subtypes = sorted(subclass_breakdown.keys())
                    
                    subtype_detail = "_".join([f"{k}({v})" for k, v in sorted(class_breakdown.items())])
                    subclass_detail = "; ".join([f"{k.split(':')[1]}({v})" for k, v in sorted(subclass_breakdown.items())])
                    
                    # Biological hotspot record
                    results.append({
                        "Sequence Name": motifs[0].get("Sequence Name", ""),
                        "Class": "Hybrid",
                        "Subtype": f"Hybrid_{subtype_detail}",
                        "Start": region_start,
                        "End": region_end,
                        "Length": total_length,
                        "MotifClasses": sorted_classes,
                        "ClassBreakdown": class_breakdown,
                        "SubclassBreakdown": subclass_breakdown,
                        "SubclassDetail": subclass_detail,
                        "ContributingMotifs": involved_motifs,
                        "OverlapDegree": round(overlap_degree, 3),   # Co-localization tightness
                        "InteractionStrength": round(interaction_bonus, 3),  # Motif synergy
                        "ScoreMethod": "HybridOverlap_SubclassAnalysis_raw",
                        "Score": float(final_score),
                        "Sequence": wrap(seq[region_start-1:region_end]),
                        "Arms/Repeat Unit/Copies": f"Classes={len(class_breakdown)};Subtypes={len(subclass_breakdown)}",
                        "Spacer": ""
                    })
            # Remove motif as its end is reached
            active.discard(idx)
    return results


def analyze_longest_hybrid(hybrid_motifs, all_motifs):
    """
    Identify the *longest* hybrid structure and evaluate its significance.
    
    Biological significance:
    ------------------------
    - Longer hybrid regions with more overlapping motifs = greater chance
      of disrupting canonical B-DNA and forming stable alternative structures
    - Such regions can act as *genomic instability hotspots* (linked to
      cancer translocations, fragile sites, repeat expansion disorders)
    - They also overlap with functional elements like promoters, replication
      origins, or enhancers where DNA flexibility is critical
    
    Parameters:
    -----------
    hybrid_motifs (list): List of detected hybrids
    all_motifs (list): All detected motifs (to count overlaps)
    
    Returns:
    --------
    dict: Longest hybrid summary with overlaps and class diversity
    """
    if not hybrid_motifs:
        return None
    
    best_hybrid = None
    max_length = 0
    max_overlap_count = 0
    
    for hybrid in hybrid_motifs:
        length = hybrid.get('Length', 0)
        overlap_count = 0
        h_start, h_end = hybrid.get('Start', 0), hybrid.get('End', 0)
        
        # Count how many motifs overlap this hybrid (excluding hybrids themselves)
        for motif in all_motifs:
            if motif.get('Class') == 'Hybrid': 
                continue
            m_start, m_end = motif.get('Start', 0), motif.get('End', 0)
            if max(h_start, m_start) <= min(h_end, m_end):
                overlap_count += 1
        
        # Pick the longest; break ties with overlap count
        if length > max_length or (length == max_length and overlap_count > max_overlap_count):
            best_hybrid = hybrid
            max_length = length
            max_overlap_count = overlap_count
    
    if not best_hybrid:
        return None
    
    return {
        'hybrid': best_hybrid,
        'length': max_length,
        'overlap_count': max_overlap_count,
        'sequence_name': best_hybrid.get('Sequence Name', ''),
        'start': best_hybrid.get('Start', 0),
        'end': best_hybrid.get('End', 0),
        'classes': best_hybrid.get('MotifClasses', []),
        'summary': f"Length={max_length}bp, Overlaps={max_overlap_count}, Classes={'+'.join(best_hybrid.get('MotifClasses', []))}"
    }
