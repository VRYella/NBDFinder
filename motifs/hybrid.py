"""
Category 9: Hybrid Structures Detection Module
==============================================

This module implements detection algorithms for hybrid structures
formed by overlapping non-B DNA motifs, collapsing nested overlaps
so that only the largest hybrid region per locus is retained.

Author: Dr. Venkata Rajesh Yella
Updated: 2025 (nested collapse + biological annotations)
"""

from .shared_utils import wrap

def find_hybrids(motifs, seq):
    """
    Find hybrid structures with detailed subclass breakdown analysis.
    Collapse nested overlaps: keep only the largest hybrid per locus.
    
    Biological Interpretation:
    --------------------------
    - Hybrids represent *combinatorial non-B DNA zones*
    - Larger hybrids = stronger candidates for *regulatory hotspots* or *fragile regions*
    - Collapsing ensures redundant nested hybrids do not inflate results
    """

    # Build timeline of motif events
    events = []
    for idx, m in enumerate(motifs):
        events.append((m['Start'], 'start', idx))
        events.append((m['End'] + 1, 'end', idx))
    events.sort()

    active = set()
    region_start = None
    raw_results = []

    for pos, typ, idx in events:
        if typ == 'start':
            active.add(idx)
            if len(active) == 2:
                region_start = pos
        elif typ == 'end':
            if len(active) >= 2:
                region_end = pos - 1
                involved_idxs = list(active)
                involved_motifs = [motifs[i] for i in involved_idxs]

                class_breakdown, subclass_breakdown = {}, {}
                for motif in involved_motifs:
                    motif_class = motif.get('Class', 'Unknown')
                    motif_subtype = motif.get('Subtype', 'Unknown')
                    class_breakdown[motif_class] = class_breakdown.get(motif_class, 0) + 1
                    subclass_key = f"{motif_class}:{motif_subtype}"
                    subclass_breakdown[subclass_key] = subclass_breakdown.get(subclass_key, 0) + 1

                # Require ≥2 classes for a hybrid
                if len(class_breakdown) >= 2:
                    total_length = region_end - region_start + 1
                    overlap_degree = total_length / max(m['Length'] for m in involved_motifs)
                    base_score = sum(float(m.get("Score", 0.0)) for m in involved_motifs)
                    interaction_bonus = len(involved_motifs) * 0.2
                    final_score = base_score * (1 + interaction_bonus) * overlap_degree

                    sorted_classes = sorted(class_breakdown.keys())
                    subtype_detail = "_".join([f"{k}({v})" for k, v in sorted(class_breakdown.items())])
                    subclass_detail = "; ".join([f"{k.split(':')[1]}({v})" for k, v in sorted(subclass_breakdown.items())])

                    raw_results.append({
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
                        "Spacer": "",
                        "Annotation": (
                            "This hybrid region spans multiple non-B DNA classes, "
                            "indicating a potential hotspot for chromatin remodeling, "
                            "replication stalling, or mutagenesis. Longer spans with "
                            "diverse motif composition have higher biological impact."
                        )
                    })
            active.discard(idx)

    # Collapse nested hybrids: retain only the largest per locus
    if not raw_results:
        return []

    raw_results.sort(key=lambda x: (x['Start'], -(x['Length'])))
    collapsed = []
    current = None

    for r in raw_results:
        if current is None:
            current = r
        else:
            if r['Start'] >= current['Start'] and r['End'] <= current['End']:
                # nested inside -> skip
                continue
            else:
                collapsed.append(current)
                current = r
    if current:
        collapsed.append(current)

    return collapsed

def analyze_longest_hybrid(hybrid_motifs, all_motifs):
    """
    Identify the *longest* hybrid region (after collapsing nested overlaps).
    Provides annotation of biological significance.
    """
    if not hybrid_motifs:
        return None

    best_hybrid = max(hybrid_motifs, key=lambda h: (h.get('Length', 0), len(h.get('MotifClasses', []))))

    h_start, h_end = best_hybrid['Start'], best_hybrid['End']
    overlap_count = 0
    for motif in all_motifs:
        if motif.get('Class') == 'Hybrid':
            continue
        m_start, m_end = motif.get('Start', 0), motif.get('End', 0)
        if max(h_start, m_start) <= min(h_end, m_end):
            overlap_count += 1

    return {
        'hybrid': best_hybrid,
        'length': best_hybrid['Length'],
        'overlap_count': overlap_count,
        'sequence_name': best_hybrid.get('Sequence Name', ''),
        'start': h_start,
        'end': h_end,
        'classes': best_hybrid.get('MotifClasses', []),
        'summary': (
            f"Longest hybrid: {best_hybrid['Length']}bp, {overlap_count} overlapping motifs, "
            f"Classes={'+'.join(best_hybrid.get('MotifClasses', []))}. "
            "Such long hybrids with multiple motif classes often correspond to fragile genomic regions, "
            "sites of transcriptional regulation, or DNA replication stress zones."
        )
    }
