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
    Find hotspots with detailed subclass breakdown analysis.
    
    Enhanced Algorithm:
    1. Sliding window approach for motif density calculation
    2. Detailed subclass-specific counting and analysis
    3. Subclass diversity scoring and comprehensive reporting
    4. Hotspot merging with subclass preservation
    
    Parameters:
    motif_hits (list): List of detected motifs
    seq_len (int): Length of the sequence
    window (int): Window size for hotspot detection
    min_count (int): Minimum motif count for hotspot
    
    Returns:
    list: List of hotspot motif dictionaries with subclass breakdown
    """
    hotspots = []
    positions = [(hit['Start'], hit['End']) for hit in motif_hits]
    
    for i in range(0, seq_len - window + 1, window // 4):  # Sliding window with overlap
        region_start, region_end = i+1, i+window
        
        # Find motifs in current window
        motifs_in_region = [m for m in motif_hits if m['Start'] <= region_end and m['End'] >= region_start]
        count = len(motifs_in_region)
        
        if count >= min_count:
            # Detailed subclass analysis
            class_breakdown = {}
            subclass_breakdown = {}
            subtype_breakdown = {}
            
            for motif in motifs_in_region:
                motif_class = motif.get('Class', 'Unknown')
                motif_subtype = motif.get('Subtype', 'Unknown')
                
                # Count by class
                class_breakdown[motif_class] = class_breakdown.get(motif_class, 0) + 1
                
                # Count by subclass (Class:Subtype combination)
                subclass_key = f"{motif_class}:{motif_subtype}"
                subclass_breakdown[subclass_key] = subclass_breakdown.get(subclass_key, 0) + 1
                
                # Count by subtype only
                subtype_breakdown[motif_subtype] = subtype_breakdown.get(motif_subtype, 0) + 1
            
            # Calculate diversity indices
            class_diversity = len(class_breakdown)
            subclass_diversity = len(subclass_breakdown)
            subtype_diversity = len(subtype_breakdown)
            
            # Enhanced scoring considering subclass complexity
            total_score = sum(float(m.get("Score", 0.0)) for m in motifs_in_region)
            diversity_bonus = (class_diversity * 10 + subclass_diversity * 5 + subtype_diversity * 2)
            final_score = total_score + diversity_bonus
            
            # Create detailed subclass description
            class_detail = "; ".join([f"{k}({v})" for k, v in sorted(class_breakdown.items())])
            subclass_detail = "; ".join([f"{k.split(':')[1]}({v})" for k, v in sorted(subclass_breakdown.items())])
            
            # Determine cluster complexity level
            if class_diversity >= 5:
                complexity_level = "High"
            elif class_diversity >= 3:
                complexity_level = "Moderate" 
            else:
                complexity_level = "Low"
                
            hotspots.append({
                "Sequence Name": motif_hits[0].get("Sequence Name", "") if motif_hits else "",
                "Class": "Non-B DNA Clusters",
                "Subtype": f"Cluster_{complexity_level}_Complexity",
                "Start": region_start,
                "End": region_end,
                "Length": region_end - region_start + 1,
                "Sequence": "",
                "ScoreMethod": "Cluster_SubclassAnalysis_raw",
                "Score": float(final_score),
                "MotifCount": count,
                "ClassDiversity": class_diversity,
                "SubclassDiversity": subclass_diversity,
                "SubtypeDiversity": subtype_diversity,
                "ClassBreakdown": class_breakdown,
                "SubclassBreakdown": subclass_breakdown,
                "SubtypeBreakdown": subtype_breakdown,
                "ClassDetail": class_detail,
                "SubclassDetail": subclass_detail,
                "ComplexityLevel": complexity_level,
                "DiversityScore": diversity_bonus,
                "Arms/Repeat Unit/Copies": f"Classes={class_diversity};Subclasses={subclass_diversity}",
                "Spacer": ""
            })
    
    return merge_hotspots_with_subclass_preservation(hotspots)

def merge_hotspots_with_subclass_preservation(hotspots):
    """
    Merge overlapping hotspots while preserving detailed subclass information.
    
    Parameters:
    hotspots (list): List of hotspot regions with subclass breakdown
    
    Returns:
    list: List of merged hotspot regions with preserved subclass analysis
    """
    if not hotspots:
        return []
    
    # Sort by start position
    hotspots_sorted = sorted(hotspots, key=lambda x: x['Start'])
    merged = [hotspots_sorted[0]]
    
    for current in hotspots_sorted[1:]:
        last = merged[-1]
        
        # Check for overlap
        if current['Start'] <= last['End']:
            # Merge hotspots with subclass preservation
            last['End'] = max(last['End'], current['End'])
            last['Length'] = last['End'] - last['Start'] + 1
            last['MotifCount'] += current['MotifCount']
            
            # Merge class breakdowns
            for class_name, count in current['ClassBreakdown'].items():
                last['ClassBreakdown'][class_name] = last['ClassBreakdown'].get(class_name, 0) + count
            
            # Merge subclass breakdowns
            for subclass_name, count in current['SubclassBreakdown'].items():
                last['SubclassBreakdown'][subclass_name] = last['SubclassBreakdown'].get(subclass_name, 0) + count
            
            # Merge subtype breakdowns
            for subtype_name, count in current['SubtypeBreakdown'].items():
                last['SubtypeBreakdown'][subtype_name] = last['SubtypeBreakdown'].get(subtype_name, 0) + count
            
            # Recalculate diversity metrics
            last['ClassDiversity'] = len(last['ClassBreakdown'])
            last['SubclassDiversity'] = len(last['SubclassBreakdown'])
            last['SubtypeDiversity'] = len(last['SubtypeBreakdown'])
            
            # Update complexity level
            if last['ClassDiversity'] >= 5:
                last['ComplexityLevel'] = "High"
            elif last['ClassDiversity'] >= 3:
                last['ComplexityLevel'] = "Moderate"
            else:
                last['ComplexityLevel'] = "Low"
            
            last['Subtype'] = f"Cluster_{last['ComplexityLevel']}_Complexity"
            
            # Recalculate diversity score and total score
            last['DiversityScore'] = (last['ClassDiversity'] * 10 + 
                                    last['SubclassDiversity'] * 5 + 
                                    last['SubtypeDiversity'] * 2)
            last['Score'] = float(last['Score']) + float(current['Score'])
            
            # Update detail strings
            last['ClassDetail'] = "; ".join([f"{k}({v})" for k, v in sorted(last['ClassBreakdown'].items())])
            last['SubclassDetail'] = "; ".join([f"{k.split(':')[1]}({v})" for k, v in sorted(last['SubclassBreakdown'].items())])
            last['Arms/Repeat Unit/Copies'] = f"Classes={last['ClassDiversity']};Subclasses={last['SubclassDiversity']}"
            
        else:
            merged.append(current)
    
    return merged

def merge_hotspots(hotspots):
    """
    Legacy merge function for backward compatibility.
    """
    return merge_hotspots_with_subclass_preservation(hotspots)

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
            'Multimeric G4', 'Bipartite G4', 'Imperfect G4', 'Canonical G4',
            'Relaxed G4', 'Bulged G4', 'G-Triplex intermediate', 'Canonical i-motif', 'AC-motif',
            'Z-DNA', 'eGZ (Extruded-G) DNA', 'Curved DNA', 'Cruciform DNA', 'Slipped DNA',
            'Triplex', 'R-loop', 'sticky DNA', 'Hybrid', 'Non-B DNA cluster regions'
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