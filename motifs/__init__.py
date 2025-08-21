"""
NBDFinder Motifs Package - Reorganized by Scientific Classification
================================================================

This package organizes non-B DNA motif detection into 10 scientific categories
for improved maintainability and specialized scoring systems.

Categories:
1. Curved DNA
2. Hairpin and Cruciform  
3. Slipped DNA
4. Triplex DNA
5. R-loop
6. Z-DNA and eGZ
7. i-motif and AC-motif
8. G4 related including G triplex
9. Hybrid
10. Non-B DNA Clusters

Author: Dr. Venkata Rajesh Yella
Updated: 2024
"""

# Import all functions to maintain backward compatibility
from .shared_utils import *
from .curved_dna import *
from .hairpin_cruciform import *
from .slipped_dna import *
from .triplex_dna import *
from .r_loop import *
from .zdna_egz import *
from .imotif_ac import *
from .g4_related import *
from .hybrid import *
from .cluster import *

# Re-export the main all_motifs function
from .shared_utils import all_motifs, format_motif_rows

__version__ = "2.0.0"
__all__ = [
    # Shared utilities
    'parse_fasta', 'wrap', 'gc_content', 'reverse_complement', 'is_palindrome',
    'overlapping_finditer', 'calculate_conservation_score', 'validate_motif',
    'get_basic_stats', 'all_motifs', 'format_motif_rows', 'clean_subtype_name',
    
    # Category 1: Curved DNA
    'find_curved_DNA', 'curvature_score', 'find_polyA_polyT_tracts',
    'find_global_curved_polyA_polyT', 'find_local_curved_polyA_polyT',
    
    # Category 2: Hairpin and Cruciform
    'find_cruciform',
    
    # Category 3: Slipped DNA
    'find_slipped_dna',
    
    # Category 4: Triplex DNA
    'find_hdna', 'find_sticky_dna', 'purine_fraction', 'pyrimidine_fraction',
    
    # Category 5: R-loop
    'find_rlfs', 'find_rez_advanced', 'advanced_rloop_score',
    
    # Category 6: Z-DNA and eGZ
    'find_zdna', 'find_egz_motif', 'zdna_dinucleotide_weights', 'kadane_maximum_subarray',
    
    # Category 7: i-motif and AC-motif
    'find_imotif', 'find_ac_motifs', 'imotif_score', 'ac_motif_score',
    
    # Category 8: G4 related
    'find_gquadruplex', 'find_relaxed_gquadruplex', 'find_bulged_gquadruplex',
    'find_imperfect_gquadruplex', 'find_multimeric_gquadruplex', 'find_bipartite_gquadruplex',
    'find_gtriplex', 'g4hunter_score', 'g4_structural_factor', 'get_g4_formation_category',
    
    # Category 9: Hybrid
    'find_hybrids',
    
    # Category 10: Clusters
    'find_hotspots', 'merge_hotspots', 'select_best_nonoverlapping_motifs'
]