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
    Find hybrid structures with detailed subclass breakdown analysis.
    
    Enhanced Algorithm:
    1. Event-based overlap detection using sorted start/end events
    2. Multi-motif overlap identification with detailed subclass analysis
    3. Subclass-specific scoring and reporting instead of aggregation
    4. Comprehensive hybrid structure characterization
    
    Parameters:
    motifs (list): List of detected motifs
    seq (str): Original DNA sequence
    
    Returns:
    list: List of hybrid structure motif dictionaries with subclass breakdown
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
            if len(active) >= 2:
                region_end = pos - 1
                involved_idxs = list(active)
                
                # Detailed subclass analysis
                involved_motifs = [motifs[i] for i in involved_idxs]
                subclass_breakdown = {}
                class_breakdown = {}
                
                for motif in involved_motifs:
                    motif_class = motif.get('Class', 'Unknown')
                    motif_subtype = motif.get('Subtype', 'Unknown')
                    
                    # Count by class
                    class_breakdown[motif_class] = class_breakdown.get(motif_class, 0) + 1
                    
                    # Count by subclass
                    subclass_key = f"{motif_class}:{motif_subtype}"
                    subclass_breakdown[subclass_key] = subclass_breakdown.get(subclass_key, 0) + 1
                
                # Only proceed if we have multiple different classes
                if len(class_breakdown) >= 2:
                    # Calculate overlap degree and interaction scores
                    total_length = region_end - region_start + 1
                    overlap_degree = total_length / max(m['Length'] for m in involved_motifs)
                    
                    # Combined scoring from contributing motifs with structural interaction factor
                    base_score = sum(float(m.get("Score", 0.0)) for m in involved_motifs)
                    interaction_bonus = len(involved_motifs) * 0.2  # Bonus for multiple interactions
                    final_score = base_score * (1 + interaction_bonus) * overlap_degree
                    
                    # Create detailed subtype identifier
                    sorted_classes = sorted(class_breakdown.keys())
                    sorted_subtypes = sorted(subclass_breakdown.keys())
                    
                    subtype_detail = "_".join([f"{k}({v})" for k, v in sorted(class_breakdown.items())])
                    subclass_detail = "; ".join([f"{k.split(':')[1]}({v})" for k, v in sorted(subclass_breakdown.items())])
                    
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
                        "OverlapDegree": round(overlap_degree, 3),
                        "InteractionStrength": round(interaction_bonus, 3),
                        "ScoreMethod": "HybridOverlap_SubclassAnalysis_raw",
                        "Score": float(final_score),
                        "Sequence": wrap(seq[region_start-1:region_end]),
                        "Arms/Repeat Unit/Copies": f"Classes={len(class_breakdown)};Subtypes={len(subclass_breakdown)}",
                        "Spacer": ""
                    })
            active.discard(idx)
    return results

def analyze_longest_hybrid(hybrid_motifs, all_motifs):
    """
    Analyze hybrid motifs to find the longest region with highest overlap count.
    
    Biological significance: Longer hybrid regions with more overlaps indicate complex
    regulatory hotspots with potential for enhanced genomic instability and function.
    
    Parameters:
    hybrid_motifs (list): List of detected hybrid motifs
    all_motifs (list): All detected motifs for overlap calculation
    
    Returns:
    dict: Analysis results with longest hybrid details or None if no hybrids
    """
    if not hybrid_motifs: return None
    
    best_hybrid = None; max_length = 0; max_overlap_count = 0
    
    for hybrid in hybrid_motifs:
        length = hybrid.get('Length', 0); overlap_count = 0
        h_start, h_end = hybrid.get('Start', 0), hybrid.get('End', 0)
        
        # Count overlapping motifs (excluding hybrid class to avoid self-counting)
        for motif in all_motifs:
            if motif.get('Class') == 'Hybrid': continue
            m_start, m_end = motif.get('Start', 0), motif.get('End', 0)
            if max(h_start, m_start) <= min(h_end, m_end): overlap_count += 1
        
        # Select based on length first, then overlap count as tiebreaker
        if length > max_length or (length == max_length and overlap_count > max_overlap_count):
            best_hybrid = hybrid; max_length = length; max_overlap_count = overlap_count
    
    if not best_hybrid: return None
    
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