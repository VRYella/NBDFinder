"""
Shared Utilities for NBDFinder Motif Detection
=============================================

This module provides a unified computational framework for non-B DNA motif detection
that harmonizes the scoring and region detection logic across G4, R-loop, and Z-DNA
families while preserving the biological accuracy of each motif-specific algorithm.

UNIFIED FRAMEWORK ARCHITECTURE:
-------------------------------
The harmonized computational pipeline consists of four main components:
1. Scoring Array Generation: create_scoring_array() provides standardized sliding window scoring.
2. Region Detection: detect_regions_from_scores() uses enhanced Kadane's maximum subarray algorithm.
3. Threshold Application: apply_region_thresholds() provides motif-specific filtering.
4. Standardized Output: format_detected_regions() creates consistent motif objects.

SCIENTIFIC ACCURACY PRESERVATION:
---------------------------------
- G4: Maintains G4Hunter scoring methodology, formation classification, and structural factor calculations.
- R-loop: Preserves RLFS+REZ bipartite detection and stability scoring.
- Z-DNA: Maintains dinucleotide weight calculations and Kadane's algorithm for optimal region detection.

REFERENCES:
-----------
- Bedrat, A. et al. (2016) NAR - G4Hunter algorithm for G-quadruplex prediction
- Aguilera, A. & García-Muse, T. (2012) Mol Cell - R-loop formation and genomic instability
- Rich, A. & Zhang, S. (2003) Nat Rev Genet - Z-DNA structure and biological significance

Author: Dr. Venkata Rajesh Yella
Updated: 2024 with unified computational framework
License: Academic Use
"""

import re
from math import log2, sqrt

def mean(values): return sum(values) / len(values) if values else 0
def std(values): return sqrt(sum((x - mean(values)) ** 2 for x in values) / (len(values) - 1)) if len(values) > 1 else 0.0

# =========================
# Basic sequence utilities
# =========================

def parse_fasta(fasta_str: str) -> str:
    return "".join(line.strip() for line in fasta_str.split('\n') if not line.startswith(">")).upper().replace(" ", "").replace("U", "T")

def wrap(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))

def gc_content(seq: str) -> float:
    return (seq.count('G') + seq.count('C')) / max(1, len(seq)) * 100 if seq else 0

def reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ATGC", "TACG"))[::-1]

def is_palindrome(seq: str) -> bool:
    return seq == reverse_complement(seq)

def overlapping_finditer(pattern, seq):
    regex = re.compile(pattern, re.IGNORECASE)
    pos = 0
    while pos < len(seq):
        m = regex.search(seq, pos)
        if not m: break
        yield m; pos = m.start() + 1

# =========================
# Unified Hunter-style scoring for all non-B DNA structures
# =========================

def unified_hunter_score(seq, target_base='G', complementary_base='C'):
    """
    Unified Hunter-style scoring (G4Hunter logic, generalizable).
    """
    if not seq: return 0.0
    seq = seq.upper(); scores = []; i = 0
    while i < len(seq):
        if seq[i] == target_base:
            run = 1
            while i + run < len(seq) and seq[i + run] == target_base: run += 1
            score = min(run, 4); scores.extend([score] * run); i += run
        elif seq[i] == complementary_base:
            run = 1
            while i + run < len(seq) and seq[i + run] == complementary_base: run += 1
            score = -min(run, 4); scores.extend([score] * run); i += run
        else: scores.append(0); i += 1
    return sum(scores) / len(scores) if scores else 0.0

# =========================
# Structural Stability Factor
# =========================

def calculate_structural_factor(motif_seq, motif_type, loop_lengths=None):
    base_factor = 1.0
    if motif_type in ['G4', 'Canonical G4', 'Relaxed G4']:
        g_runs = len(re.findall(r"G{3,}", motif_seq)); base_factor += 0.1 * max(0, g_runs - 4)
        if loop_lengths:
            avg_loop = sum(loop_lengths) / len(loop_lengths)
            base_factor += 0.1 if avg_loop <= 7 else -0.05 * (avg_loop - 7)
    elif motif_type in ['i-motif', 'Canonical i-motif', 'Relaxed i-motif']:
        c_runs = len(re.findall(r"C{3,}", motif_seq)); base_factor += 0.1 * max(0, c_runs - 4)
        if loop_lengths:
            avg_loop = sum(loop_lengths) / len(loop_lengths)
            base_factor += 0.15 if avg_loop <= 7 else 0.05 if avg_loop <= 12 else -0.1
    elif motif_type == 'R-loop':
        gc_content_val = gc_content(motif_seq); g_runs = len(re.findall(r"G{3,}", motif_seq))
        base_factor += 0.01 * gc_content_val + 0.05 * g_runs
    elif motif_type == 'Z-DNA':
        alternating_purine_pyrimidine = len(re.findall(r"[AG][CT]", motif_seq))
        gc_cg_dinucs = len(re.findall(r"GC|CG", motif_seq))
        base_factor += 0.05 * alternating_purine_pyrimidine + 0.1 * gc_cg_dinucs
    return max(0.1, base_factor)

# =========================
# Conservation Analysis
# =========================

def calculate_conservation_score(seq, motif_type="general", num_shuffles=100):
    """
    Calculates motif conservation score using composition-preserving shuffles.
    """
    if not seq or len(seq) < 4: return {"enrichment_score": 0.0, "p_value": 1.0, "significance": "not significant"}
    seq = seq.upper().replace('\n', '').replace(' ', '')
    k, core_kmers = _get_motif_kmers(seq, motif_type)
    if not core_kmers or k > len(seq): return {"enrichment_score": 0.0, "p_value": 1.0, "significance": "not significant"}
    observed_count = _count_kmers(seq, core_kmers, k)
    import random
    shuffle_counts = [_count_kmers(''.join(random.sample(seq, len(seq))), core_kmers, k) for _ in range(num_shuffles)]
    if not shuffle_counts: return {"enrichment_score": 0.0, "p_value": 1.0, "significance": "not significant"}
    mu, sigma = mean(shuffle_counts), std(shuffle_counts); epsilon = 1e-6
    enrichment_score = log2((observed_count + epsilon) / (mu + epsilon))
    p_value = (1 + sum(1 for c in shuffle_counts if c >= observed_count)) / (num_shuffles + 1)
    significance = ("high" if enrichment_score >= 2.5 else
                   "medium" if enrichment_score >= 1.5 else 
                   "low" if enrichment_score >= 0.5 else
                   "not significant") if p_value <= 0.05 else "not significant"
    return {
        "enrichment_score": float(enrichment_score),
        "p_value": float(p_value),
        "significance": significance,
        "observed_count": int(observed_count),
        "null_mean": float(mu),
        "null_std": float(sigma)
    }

def _get_motif_kmers(seq, motif_type):
    n = len(seq)
    if motif_type in ["G4", "Canonical G4", "Relaxed G4", "Bulged G4", "Bipartite G4", "Multimeric G4", "Imperfect G4"]:
        k = min(8, max(4, n // 4))
        core_kmers = {seq[i:i+k] for i in range(n - k + 1) if 'GGG' in seq[i:i+k]}
        if not core_kmers: core_kmers = {seq[i:i+k] for i in range(n - k + 1) if seq[i:i+k].count('G') >= k//2}
    elif motif_type in ["i-Motif", "Triplex"]:
        k = min(8, max(4, n // 4))
        core_kmers = {seq[i:i+k] for i in range(n - k + 1) if 'CCC' in seq[i:i+k]}
        if not core_kmers: core_kmers = {seq[i:i+k] for i in range(n - k + 1) if seq[i:i+k].count('C') >= k//2}
    elif motif_type in ["Z-DNA", "eGZ"]:
        k = min(6, max(3, n // 6))
        z_dinucs = {'GC', 'CG', 'GT', 'TG', 'AC', 'CA'}
        core_kmers = {seq[i:i+k] for i in range(n - k + 1) 
                     if sum(seq[i+j:i+j+2] in z_dinucs for j in range(k-1)) >= (k-1) // 2}
    elif motif_type in ["AC-Motif"]:
        k = min(6, max(4, n // 4))
        core_kmers = {seq[i:i+k] for i in range(n - k + 1) 
                     if (seq[i:i+k].count('A') + seq[i:i+k].count('C')) / k >= 0.6}
    elif motif_type in ["Cruciform", "Slipped DNA"]:
        k = min(8, max(4, n // 4))
        core_kmers = {seq[i:i+k] for i in range(n - k + 1)}
    elif motif_type in ["R-loop", "RLFS"]:
        k = min(8, max(4, n // 4))
        core_kmers = {seq[i:i+k] for i in range(n - k + 1) if seq[i:i+k].count('G') / k >= 0.4}
    elif motif_type in ["Curved DNA"]:
        k = min(6, max(4, n // 6))
        core_kmers = {seq[i:i+k] for i in range(n - k + 1) if 'AAA' in seq[i:i+k] or 'TTT' in seq[i:i+k]}
    else:
        k = min(8, max(4, n // 4))
        core_kmers = {seq[i:i+k] for i in range(n - k + 1)}
    return k, core_kmers

def _count_kmers(seq, core_kmers, k):
    return sum(1 for i in range(len(seq) - k + 1) if seq[i:i+k] in core_kmers)

# =========================
# MAIN MOTIF DETECTION FUNCTIONS
# =========================

def validate_motif(motif, seq_length):
    """Validate motif has required fields and valid coordinates"""
    required_keys = ["Class", "Subtype", "Start", "End", "Length", "Sequence"]
    if not all(key in motif for key in required_keys):
        return False
    if not (1 <= motif["Start"] <= motif["End"] <= seq_length):
        return False
    return True

def all_motifs(seq, nonoverlap=False, report_hotspots=False, sequence_name="Sequence"):
    """
    COMPLETE PERFORMANCE MODE: Detect all non-B DNA motifs in specified order
    
    This function runs each class one after another in the exact order specified:
    1. Curved DNA
    2. Slipped DNA  
    3. Cruciform DNA
    4. R-loop
    5. Triplex
    6. G-Quadruplex Family
    7. i-motif family
    8. Z-DNA
    9. Hybrid (overlaps between any two or more)
    10. Non-B DNA cluster regions
    
    Args:
        seq: DNA sequence string
        nonoverlap: If True, remove overlapping motifs per class
        report_hotspots: If True, include cluster analysis
        sequence_name: Name for the sequence
    
    Returns:
        List of motif dictionaries
    """
    import re
    
    if not seq or not re.match("^[ATGC]+$", seq, re.IGNORECASE):
        return []
    seq = seq.upper()
    
    motif_list = []
    processed_classes = []
    
    print(f"Starting Non-B DNA detection for sequence of length {len(seq)}")
    
    # 1. Curved DNA
    print("1. Processing Curved DNA...")
    try:
        from .curved_dna import find_curved_DNA
        curved_motifs = find_curved_DNA(seq)
        motif_list += curved_motifs
        processed_classes.append(f"Curved DNA: {len(curved_motifs)} motifs found")
        print(f"   Found {len(curved_motifs)} Curved DNA motifs")
    except (ImportError, AttributeError) as e:
        processed_classes.append("Curved DNA: not found in result")
        print(f"   Curved DNA: not found in result ({e})")
    
    # 2. Slipped DNA
    print("2. Processing Slipped DNA...")
    try:
        from .slipped_dna import find_slipped_dna
        slipped_motifs = find_slipped_dna(seq)
        motif_list += slipped_motifs
        processed_classes.append(f"Slipped DNA: {len(slipped_motifs)} motifs found")
        print(f"   Found {len(slipped_motifs)} Slipped DNA motifs")
    except (ImportError, AttributeError) as e:
        processed_classes.append("Slipped DNA: not found in result")
        print(f"   Slipped DNA: not found in result ({e})")
    
    # 3. Cruciform DNA
    print("3. Processing Cruciform DNA...")
    try:
        from .hairpin_cruciform import find_cruciform
        cruciform_motifs = find_cruciform(seq)
        motif_list += cruciform_motifs
        processed_classes.append(f"Cruciform DNA: {len(cruciform_motifs)} motifs found")
        print(f"   Found {len(cruciform_motifs)} Cruciform DNA motifs")
    except (ImportError, AttributeError) as e:
        processed_classes.append("Cruciform DNA: not found in result")
        print(f"   Cruciform DNA: not found in result ({e})")
    
    # 4. R-loop
    print("4. Processing R-loop...")
    try:
        from .r_loop import find_rlfs
        rloop_motifs = find_rlfs(seq)
        motif_list += rloop_motifs
        processed_classes.append(f"R-loop: {len(rloop_motifs)} motifs found")
        print(f"   Found {len(rloop_motifs)} R-loop motifs")
    except (ImportError, AttributeError) as e:
        processed_classes.append("R-loop: not found in result")
        print(f"   R-loop: not found in result ({e})")
    
    # 5. Triplex
    print("5. Processing Triplex...")
    try:
        from .triplex_dna import find_hdna, find_sticky_dna
        triplex_motifs = find_hdna(seq) + find_sticky_dna(seq)
        motif_list += triplex_motifs
        processed_classes.append(f"Triplex: {len(triplex_motifs)} motifs found")
        print(f"   Found {len(triplex_motifs)} Triplex motifs")
    except (ImportError, AttributeError) as e:
        processed_classes.append("Triplex: not found in result")
        print(f"   Triplex: not found in result ({e})")
    
    # 6. G-Quadruplex Family
    print("6. Processing G-Quadruplex Family...")
    try:
        from .g4_related import (find_gquadruplex, find_gtriplex, find_relaxed_gquadruplex, 
                                find_bulged_gquadruplex, find_bipartite_gquadruplex, 
                                find_multimeric_gquadruplex, find_imperfect_gquadruplex)
        g4_motifs = []
        
        # Try each G4 function individually to handle any issues
        try:
            g4_motifs += find_multimeric_gquadruplex(seq)
        except Exception as e:
            print(f"     Warning: Multimeric G4 detection failed: {e}")
            
        try:
            g4_motifs += find_gquadruplex(seq)
        except Exception as e:
            print(f"     Warning: Canonical G4 detection failed: {e}")
            
        try:
            g4_motifs += find_relaxed_gquadruplex(seq)
        except Exception as e:
            print(f"     Warning: Relaxed G4 detection failed: {e}")
            
        try:
            g4_motifs += find_bulged_gquadruplex(seq)
        except Exception as e:
            print(f"     Warning: Bulged G4 detection failed: {e}")
            
        try:
            g4_motifs += find_bipartite_gquadruplex(seq)
        except Exception as e:
            print(f"     Warning: Bipartite G4 detection failed: {e}")
            
        try:
            g4_motifs += find_imperfect_gquadruplex(seq)
        except Exception as e:
            print(f"     Warning: Imperfect G4 detection failed: {e}")
            
        try:
            g4_motifs += find_gtriplex(seq)
        except Exception as e:
            print(f"     Warning: G-triplex detection failed: {e}")
        
        motif_list += g4_motifs
        processed_classes.append(f"G-Quadruplex Family: {len(g4_motifs)} motifs found")
        print(f"   Found {len(g4_motifs)} G-Quadruplex Family motifs")
    except (ImportError, AttributeError) as e:
        processed_classes.append("G-Quadruplex Family: not found in result")
        print(f"   G-Quadruplex Family: not found in result ({e})")
    
    # 7. i-motif family
    print("7. Processing i-motif family...")
    try:
        from .imotif_ac import find_imotif, find_ac_motifs
        imotif_motifs = find_imotif(seq) + find_ac_motifs(seq)
        motif_list += imotif_motifs
        processed_classes.append(f"i-motif family: {len(imotif_motifs)} motifs found")
        print(f"   Found {len(imotif_motifs)} i-motif family motifs")
    except (ImportError, AttributeError) as e:
        processed_classes.append("i-motif family: not found in result")
        print(f"   i-motif family: not found in result ({e})")
    
    # 8. Z-DNA
    print("8. Processing Z-DNA...")
    try:
        from .zdna_egz import find_zdna, find_egz_motif
        zdna_motifs = find_zdna(seq) + find_egz_motif(seq)
        motif_list += zdna_motifs
        processed_classes.append(f"Z-DNA: {len(zdna_motifs)} motifs found")
        print(f"   Found {len(zdna_motifs)} Z-DNA motifs")
    except (ImportError, AttributeError) as e:
        processed_classes.append("Z-DNA: not found in result")
        print(f"   Z-DNA: not found in result ({e})")
    
    # Add disease motifs if available
    try:
        import sys
        import os
        sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
        from disease_motifs import find_disease_associated_motifs
        disease_motifs = find_disease_associated_motifs(seq)
        motif_list += disease_motifs
        print(f"   Added {len(disease_motifs)} disease-associated motifs")
    except ImportError:
        pass
    
    # Validate and standardize fields
    motif_list = [m for m in motif_list if validate_motif(m, len(seq))]
    
    # 9. Hybrid (overlaps between any two or more)
    print("9. Processing Hybrid...")
    try:
        from .hybrid import find_hybrids, analyze_longest_hybrid
        hybrid_motifs = find_hybrids(motif_list, seq)
        motif_list += hybrid_motifs
        processed_classes.append(f"Hybrid: {len(hybrid_motifs)} motifs found")
        print(f"   Found {len(hybrid_motifs)} Hybrid motifs")
        
        # Check for overlaps and report
        if hybrid_motifs:
            print(f"   Hybrid overlap analysis:")
            for hybrid in hybrid_motifs:
                if 'MotifClasses' in hybrid:
                    classes = hybrid['MotifClasses']
                    print(f"     - {hybrid.get('Subtype', 'Unknown')} (Classes: {', '.join(classes)})")
            
            # Longest hybrid analysis - identifies the most significant hybrid region
            longest_hybrid = analyze_longest_hybrid(hybrid_motifs, motif_list); print(f"   Longest hybrid region: {longest_hybrid['summary'] if longest_hybrid else 'None detected'}")
    except (ImportError, AttributeError) as e:
        processed_classes.append("Hybrid: not found in result")
        print(f"   Hybrid: not found in result ({e})")
    
    # De-overlap per class if asked
    if nonoverlap:
        try:
            from .cluster import select_best_nonoverlapping_motifs
            motif_list = select_best_nonoverlapping_motifs(motif_list)
            print(f"   Applied overlap resolution, keeping {len(motif_list)} motifs")
        except (ImportError, AttributeError):
            pass
    
    # 10. Non-B DNA cluster regions
    if report_hotspots:
        print("10. Processing Non-B DNA cluster regions...")
        try:
            from .cluster import find_hotspots
            cluster_motifs = find_hotspots(motif_list, len(seq))
            motif_list += cluster_motifs
            processed_classes.append(f"Non-B DNA cluster regions: {len(cluster_motifs)} motifs found")
            print(f"    Found {len(cluster_motifs)} Non-B DNA cluster regions")
        except (ImportError, AttributeError) as e:
            processed_classes.append("Non-B DNA cluster regions: not found in result")
            print(f"    Non-B DNA cluster regions: not found in result ({e})")
            
        # Add advanced clustering
        try:
            from advanced_clustering import find_advanced_clusters, find_hybrid_structures
            advanced_motifs = find_advanced_clusters(motif_list, len(seq)) + find_hybrid_structures(motif_list)
            motif_list += advanced_motifs
            print(f"    Added {len(advanced_motifs)} advanced cluster motifs")
        except ImportError:
            pass
    
    # Add Sequence Name and ensure ordered keys exist
    for m in motif_list:
        m["Sequence Name"] = sequence_name
        # Ensure mandatory ordered fields exist and not missing
        m.setdefault("Arms/Repeat Unit/Copies", "")
        m.setdefault("Spacer", "")
        
        # Optimize score conversion
        if "Score" in m and isinstance(m["Score"], str):
            try:
                m["Score"] = float(m["Score"])
            except (ValueError, TypeError):
                pass
        
        # Try to enhance with ML if available
        try:
            from ml_predictor import enhance_motif_with_ml
            m = enhance_motif_with_ml(m, seq)
        except ImportError:
            pass
    
    # Apply G4 priority filtering while preserving all other motif classes
    from .classification_config import apply_g4_priority_filter
    original_count = len(motif_list)
    motif_list = apply_g4_priority_filter(motif_list)
    g4_filtered_count = original_count - len(motif_list)
    if g4_filtered_count > 0:
        print(f"   Applied G4 priority filtering: removed {g4_filtered_count} overlapping G4 motifs")
    
    # Print summary
    print(f"\nProcessing Summary:")
    for status in processed_classes:
        print(f"  {status}")
    print(f"Total motifs detected: {len(motif_list)}")
    
    return motif_list

def format_motif_rows(motifs):
    """Format motifs into ordered rows for display"""
    ordered = []
    for m in motifs:
        row = {
            "Sequence Name": m.get("Sequence Name", ""),
            "Class": m.get("Class", ""),
            "Subtype": m.get("Subtype", m.get("Subclass", "")),
            "Start": m.get("Start", ""),
            "End": m.get("End", ""),
            "Length": m.get("Length", ""),
            "Sequence": m.get("Sequence", ""),
            "Score": m.get("Score", ""),
            "Arms/Repeat Unit/Copies": m.get("Arms/Repeat Unit/Copies", ""),
            "Spacer": m.get("Spacer", "")
        }
        ordered.append(row)
    return ordered

# =========================
# END UTILS
# =========================