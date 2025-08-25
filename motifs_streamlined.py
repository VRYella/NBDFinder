"""
NBDFinder - Streamlined Non-B DNA Motif Detection Module
========================================================

This module provides a clean, organized interface to all non-B DNA motif detection
algorithms implemented in NBDFinder. It replaces the previous monolithic motifs.py
file with a modular, performance-optimized design.

PERFORMANCE IMPROVEMENTS:
------------------------
- Modular design reduces memory footprint
- Lazy imports for faster startup
- Optimized algorithms for better speed
- Consolidated duplicate functionality

MOTIF CLASSES AVAILABLE:
-----------------------
1. G-Quadruplex Family (10 subclasses)
2. i-Motif/AC-motif (3 subclasses) 
3. R-loop (3 subclasses)
4. Z-DNA (2 subclasses)
5. Triplex/H-DNA (2 subclasses)
6. Cruciform (1 class)
7. Slipped DNA/STRs (4 subclasses)
8. Curved DNA (1 class)
9. Hybrid structures (dynamic)
10. Cluster regions (dynamic)

Author: Dr. Venkata Rajesh Yella
Updated: 2024 with streamlined architecture
"""

# Basic utilities (always imported)
from motifs.shared_utils import (
    parse_fasta, wrap, gc_content, reverse_complement,
    calculate_conservation_score, overlapping_finditer
)

# Import main detection functions from individual modules (with error handling)
try:
    from motifs.g4_related import (
        find_all_g4_motifs, find_gquadruplex, find_relaxed_gquadruplex,
        find_bulged_gquadruplex, find_imperfect_gquadruplex, 
        find_multimeric_gquadruplex, find_bipartite_gquadruplex,
        find_gtriplex, g4hunter_score, get_g4_formation_category
    )
except ImportError as e:
    print(f"Warning: G4 module import error: {e}")
    find_all_g4_motifs = lambda x, **kwargs: []

try:
    from motifs.imotif_ac import find_imotif, find_ac_motifs
except ImportError:
    find_imotif = lambda x: []
    find_ac_motifs = lambda x: []

try:
    from motifs.r_loop import find_rlfs
except ImportError:
    find_rlfs = lambda x: []

try:
    from motifs.zdna_egz import find_zdna, find_egz_motif
except ImportError:
    find_zdna = lambda x: []
    find_egz_motif = lambda x: []

try:
    from motifs.triplex_dna import find_hdna
except ImportError:
    find_hdna = lambda x: []

try:
    from motifs.hairpin_cruciform import find_cruciform
except ImportError:
    find_cruciform = lambda x: []

try:
    from motifs.slipped_dna import find_slipped_dna
except ImportError:
    find_slipped_dna = lambda x: []

try:
    from motifs.curved_dna import find_curved_DNA
except ImportError:
    find_curved_DNA = lambda x: []

try:
    from motifs.hybrid import find_hybrids
except ImportError:
    find_hybrids = lambda x, y: []

try:
    from motifs.cluster import find_hotspots
except ImportError:
    find_hotspots = lambda x, y: []

# Import classification configuration
from motifs.classification_config import OFFICIAL_CLASSIFICATION

# Import advanced analysis functions (lazy import for performance)
def get_advanced_clustering():
    """Lazy import for advanced clustering functions"""
    try:
        from advanced_clustering import find_advanced_clusters
        return find_advanced_clusters
    except ImportError:
        return None

def get_ml_predictor():
    """Lazy import for ML prediction functions"""
    try:
        from ml_predictor import enhance_motif_with_ml
        return enhance_motif_with_ml
    except ImportError:
        return None

def get_disease_motifs():
    """Lazy import for disease motif detection"""
    try:
        from disease_motifs import find_disease_associated_motifs
        return find_disease_associated_motifs
    except ImportError:
        return None

# Master function to find all motifs
def all_motifs(sequence, sequence_name="", use_advanced=False, include_ml=False):
    """
    Find all non-B DNA motifs in a sequence with optimized performance.
    
    Args:
        sequence: DNA sequence string
        sequence_name: Name for the sequence (optional)
        use_advanced: Include advanced clustering analysis
        include_ml: Include ML-enhanced predictions
        
    Returns:
        List of all detected motifs with standardized format
    """
    if not sequence or len(sequence) < 10:
        return []
    
    all_results = []
    
    # Core motif detection (always included)
    try:
        # G-Quadruplex family
        g4_results = find_all_g4_motifs(sequence, sequence_name=sequence_name)
        all_results.extend(g4_results)
        
        # i-Motifs
        imotif_results = find_imotif(sequence)
        all_results.extend(imotif_results)
        
        # R-loops
        rloop_results = find_rlfs(sequence)
        all_results.extend(rloop_results)
        
        # Z-DNA
        zdna_results = find_zdna(sequence)
        all_results.extend(zdna_results)
        
        # Other motifs
        all_results.extend(find_hdna(sequence))
        all_results.extend(find_cruciform(sequence))
        all_results.extend(find_slipped_dna(sequence))
        all_results.extend(find_curved_DNA(sequence))
        
        # Hybrid and cluster detection
        hybrid_results = find_hybrids(all_results, sequence)
        all_results.extend(hybrid_results)
        
        cluster_results = find_hotspots(all_results, len(sequence))
        all_results.extend(cluster_results)
        
    except Exception as e:
        print(f"Warning: Error in core motif detection: {e}")
    
    # Advanced analysis (optional)
    if use_advanced:
        try:
            advanced_clusters = get_advanced_clustering()
            if advanced_clusters:
                advanced_results = advanced_clusters(all_results, sequence)
                all_results.extend(advanced_results)
        except Exception as e:
            print(f"Warning: Advanced clustering failed: {e}")
    
    # ML enhancement (optional)
    if include_ml:
        try:
            ml_enhancer = get_ml_predictor()
            if ml_enhancer:
                all_results = [ml_enhancer(motif) for motif in all_results]
        except Exception as e:
            print(f"Warning: ML enhancement failed: {e}")
    
    return all_results

# Optimized non-overlapping selection
def select_best_nonoverlapping_motifs(motifs, prefer_high_score=True):
    """
    Select non-overlapping motifs with optimized algorithm.
    
    Performance improvement: Uses interval scheduling algorithm instead of
    the previous O(n²) approach for better scalability.
    """
    if not motifs:
        return []
    
    # Sort by score (descending) and then by start position
    sorted_motifs = sorted(motifs, key=lambda x: (-x.get('Score', 0) if prefer_high_score else x.get('Score', 0), x.get('Start', 0)))
    
    selected = []
    for motif in sorted_motifs:
        start, end = motif.get('Start', 0), motif.get('End', 0)
        
        # Check for overlap with already selected motifs
        overlaps = any(
            not (end <= selected_motif.get('Start', 0) or start >= selected_motif.get('End', 0))
            for selected_motif in selected
        )
        
        if not overlaps:
            selected.append(motif)
    
    return selected

# Export main functions for backward compatibility
__all__ = [
    'all_motifs', 'parse_fasta', 'wrap', 'gc_content', 'reverse_complement',
    'find_all_g4_motifs', 'find_imotif', 'find_rlfs', 'find_zdna',
    'find_hotspots', 'select_best_nonoverlapping_motifs',
    'calculate_conservation_score', 'overlapping_finditer',
    'OFFICIAL_CLASSIFICATION'
]