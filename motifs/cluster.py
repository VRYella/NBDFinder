"""
Category 10: Non-B DNA Clusters Detection Module
===============================================

This module implements detection algorithms for non-B DNA clusters
and hotspots where multiple motifs co-localize.

Scientific Basis:
- Clustering of non-B DNA motifs creates regulatory hotspots
- High motif density regions are associated with genetic instability
- Important for understanding genome organization and disease mechanisms

Author: Dr. Venkata Rajesh Yella
Updated: 2024
"""

def find_hotspots(motif_hits, seq_len, window=100, min_count=3):
    """
    Find hotspots with high density of non-B DNA motifs.
    
    Enhanced Algorithm:
    1. Sliding window approach for motif density calculation
    2. Motif count and type diversity scoring
    3. Hotspot merging to avoid redundant overlapping regions
    
    Parameters:
    motif_hits (list): List of detected motifs
    seq_len (int): Length of the sequence
    window (int): Window size for hotspot detection
    min_count (int): Minimum motif count for hotspot
    
    Returns:
    list: List of hotspot motif dictionaries
    """
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
    """
    Merge overlapping hotspots to avoid redundancy.
    
    Parameters:
    hotspots (list): List of hotspot regions
    
    Returns:
    list: List of merged hotspot regions
    """
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

def select_best_nonoverlapping_motifs(motifs: list, motif_priority: list = None) -> list:
    """
    Select best non-overlapping motifs using priority-based selection.
    
    Enhanced Algorithm:
    1. Sort motifs by priority class and score
    2. Greedy selection avoiding significant overlaps
    3. Maintains highest quality motifs per region
    
    Parameters:
    motifs (list): List of all detected motifs
    motif_priority (list): Priority order for motif classes
    
    Returns:
    list: List of non-overlapping high-quality motifs
    """
    if motif_priority is None:
        motif_priority = [
            'Multimeric_G4', 'Bipartite_G4', 'Imperfect_G4', 'Canonical_G4',
            'Relaxed_G4', 'Bulged_G4', 'G-Triplex', 'i-Motif', 'AC-Motif',
            'Z-DNA', 'eGZ (Extruded-G)', 'Curved_DNA', 'Cruciform', 'Slipped_DNA',
            'Triplex_DNA', 'R-Loop', 'Sticky_DNA', 'Hybrid', 'Non-B DNA Clusters'
        ]
    
    # Create priority mapping
    priority_map = {motif: i for i, motif in enumerate(motif_priority)}
    
    # Sort by priority (lower index = higher priority), then by score descending
    def sort_key(motif):
        class_name = motif.get('Class', '')
        subtype = motif.get('Subtype', '')
        # Try class first, then subtype
        priority = priority_map.get(class_name, priority_map.get(subtype, len(motif_priority)))
        score = float(motif.get('Score', 0))
        return (priority, -score)
    
    sorted_motifs = sorted(motifs, key=sort_key)
    
    selected = []
    used_positions = set()
    
    for motif in sorted_motifs:
        start = motif.get('Start', 0)
        end = motif.get('End', 0)
        motif_positions = set(range(start, end + 1))
        
        # Check for significant overlap (>50% of motif)
        overlap = len(motif_positions.intersection(used_positions))
        overlap_ratio = overlap / len(motif_positions) if motif_positions else 0
        
        if overlap_ratio <= 0.5:  # Allow up to 50% overlap
            selected.append(motif)
            used_positions.update(motif_positions)
    
    return selected