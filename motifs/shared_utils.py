"""
Shared Utilities for NBDFinder Motif Detection
=============================================

Common utilities used across all motif detection modules.
Includes sequence manipulation, conservation analysis, and validation functions.

Author: Dr. Venkata Rajesh Yella
Updated: 2024
"""

import re

def mean(values):
    """Simple mean calculation without numpy."""
    return sum(values) / len(values) if values else 0

def std(values):
    """Simple standard deviation calculation without numpy."""
    if len(values) <= 1:
        return 0.0
    m = mean(values)
    variance = sum((x - m) ** 2 for x in values) / (len(values) - 1)
    return variance ** 0.5

def log2(x):
    """Simple log2 calculation."""
    import math
    return math.log2(x)

# =========================
# Basic sequence utilities
# =========================

def parse_fasta(fasta_str: str) -> str:
    return "".join([line.strip() for line in fasta_str.split('\n') if not line.startswith(">")]).upper().replace(" ", "").replace("U", "T")

def wrap(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))

def gc_content(seq: str) -> float:
    gc = seq.count('G') + seq.count('C')
    return (gc / max(1, len(seq))) * 100 if seq else 0

def reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ATGC", "TACG"))[::-1]

def is_palindrome(seq: str) -> bool:
    return seq == reverse_complement(seq)

def overlapping_finditer(pattern, seq):
    regex = re.compile(pattern, re.IGNORECASE)
    pos = 0
    while pos < len(seq):
        m = regex.search(seq, pos)
        if not m:
            break
        yield m
        pos = m.start() + 1

# =========================
# Conservation Analysis
# =========================

def calculate_conservation_score(seq, motif_type="general", num_shuffles=100):
    """
    Calculate motif-agnostic conservation score using composition-preserving shuffles.
    
    Scientific Basis: This implements a precise, math-first recipe for conservation scoring
    that can be attached to any non-B DNA motif. Uses k-mer analysis with empirical
    null distribution from composition-preserving shuffles.
    
    Algorithm:
    1. Choose k-mer length based on motif type and sequence
    2. Identify motif-defining k-mers (core structural elements)
    3. Count k-mer occurrences in original sequence
    4. Generate composition-preserving shuffles for null distribution
    5. Calculate log2 enrichment score and empirical p-value
    
    Parameters:
    seq (str): DNA sequence to score
    motif_type (str): Motif type for k-mer selection strategy
    num_shuffles (int): Number of composition-preserving shuffles (default: 100)
    
    Returns:
    dict: Conservation metrics with enrichment score, p-value, and significance
    """
    if not seq or len(seq) < 4:
        return {"enrichment_score": 0.0, "p_value": 1.0, "significance": "not significant"}
    
    seq = seq.upper().replace('\n', '').replace(' ', '')
    n = len(seq)
    
    # Step 1: Choose k-mer length and define motif-specific k-mers
    k, core_kmers = _get_motif_kmers(seq, motif_type)
    
    if not core_kmers or k > n:
        return {"enrichment_score": 0.0, "p_value": 1.0, "significance": "not significant"}
    
    # Step 2: Count observed k-mer occurrences
    observed_count = _count_kmers(seq, core_kmers, k)
    
    # Step 3: Generate composition-preserving shuffles and count k-mers
    import random
    shuffle_counts = []
    
    for _ in range(num_shuffles):
        # Create composition-preserving shuffle
        shuffled_seq = ''.join(random.sample(seq, len(seq)))
        shuffle_count = _count_kmers(shuffled_seq, core_kmers, k)
        shuffle_counts.append(shuffle_count)
    
    # Step 4: Calculate empirical statistics
    if not shuffle_counts:
        return {"enrichment_score": 0.0, "p_value": 1.0, "significance": "not significant"}
    
    mu = mean(shuffle_counts)  # Null mean
    sigma = std(shuffle_counts) if len(shuffle_counts) > 1 else 0.0  # Null std
    
    # Step 5: Calculate log2 enrichment score
    epsilon = 1e-6  # Pseudocount to avoid division by zero
    enrichment_score = log2((observed_count + epsilon) / (mu + epsilon))
    
    # Step 6: Calculate empirical p-value (right-tailed test for enrichment)
    p_value = (1 + sum(1 for c in shuffle_counts if c >= observed_count)) / (num_shuffles + 1)
    
    # Step 7: Determine significance with effect size categories
    if p_value > 0.05:
        significance = "not significant"
    elif enrichment_score >= 2.5:
        significance = "high"
    elif enrichment_score >= 1.5:
        significance = "medium"
    elif enrichment_score >= 0.5:
        significance = "low"
    else:
        significance = "not significant"
    
    return {
        "enrichment_score": float(enrichment_score),
        "p_value": float(p_value),
        "significance": significance,
        "observed_count": int(observed_count),
        "null_mean": float(mu),
        "null_std": float(sigma)
    }


def _get_motif_kmers(seq, motif_type):
    """
    Define motif-specific k-mer selection strategy.
    
    Returns:
    tuple: (k_length, set_of_core_kmers)
    """
    import re
    
    n = len(seq)
    
    # Choose k based on sequence length and motif type
    if motif_type in ["G4", "Canonical G4", "Relaxed G4", "Bulged G4", "Bipartite G4", "Multimeric G4", "Imperfect G4"]:
        # G-quadruplex: Focus on GGG k-mers and G-runs
        k = min(8, max(4, n // 4))
        core_kmers = set()
        # Add all GGG-containing k-mers
        for i in range(n - k + 1):
            kmer = seq[i:i+k]
            if 'GGG' in kmer:
                core_kmers.add(kmer)
        # If no GGG k-mers, use all G-rich k-mers
        if not core_kmers:
            for i in range(n - k + 1):
                kmer = seq[i:i+k]
                if kmer.count('G') >= k//2:
                    core_kmers.add(kmer)
    
    elif motif_type in ["i-Motif", "Triplex"]:
        # i-Motif: Focus on CCC k-mers and C-runs
        k = min(8, max(4, n // 4))
        core_kmers = set()
        for i in range(n - k + 1):
            kmer = seq[i:i+k]
            if 'CCC' in kmer:
                core_kmers.add(kmer)
        if not core_kmers:
            for i in range(n - k + 1):
                kmer = seq[i:i+k]
                if kmer.count('C') >= k//2:
                    core_kmers.add(kmer)
    
    elif motif_type in ["Z-DNA", "eGZ"]:
        # Z-DNA: Focus on alternating dinucleotides
        k = min(6, max(3, n // 6))
        core_kmers = set()
        z_dinucs = ['GC', 'CG', 'GT', 'TG', 'AC', 'CA']
        for i in range(n - k + 1):
            kmer = seq[i:i+k]
            # Count Z-favorable dinucleotides in k-mer
            z_count = sum(kmer[j:j+2] in z_dinucs for j in range(len(kmer)-1))
            if z_count >= (k-1) // 2:  # At least half are Z-favorable
                core_kmers.add(kmer)
    
    elif motif_type in ["AC-Motif"]:
        # AC-motif: Focus on alternating A/C patterns
        k = min(6, max(4, n // 4))
        core_kmers = set()
        for i in range(n - k + 1):
            kmer = seq[i:i+k]
            ac_content = (kmer.count('A') + kmer.count('C')) / k
            if ac_content >= 0.6:  # 60% A/C content
                core_kmers.add(kmer)
    
    elif motif_type in ["Cruciform", "Slipped DNA"]:
        # Cruciform/Slipped: Focus on palindromic/repeat patterns
        k = min(8, max(4, n // 4))
        core_kmers = set()
        for i in range(n - k + 1):
            kmer = seq[i:i+k]
            # Add all k-mers for repeat analysis
            core_kmers.add(kmer)
    
    elif motif_type in ["R-Loop", "RLFS"]:
        # R-loop: Focus on G-rich regions (purine-rich)
        k = min(8, max(4, n // 4))
        core_kmers = set()
        for i in range(n - k + 1):
            kmer = seq[i:i+k]
            g_content = kmer.count('G') / k
            if g_content >= 0.4:  # 40% G content
                core_kmers.add(kmer)
    
    elif motif_type in ["Curved DNA"]:
        # Curved DNA: Focus on A/T tracts
        k = min(6, max(4, n // 6))
        core_kmers = set()
        for i in range(n - k + 1):
            kmer = seq[i:i+k]
            if 'AAA' in kmer or 'TTT' in kmer:
                core_kmers.add(kmer)
    
    else:
        # General case: use all k-mers
        k = min(8, max(4, n // 4))
        core_kmers = set()
        for i in range(n - k + 1):
            core_kmers.add(seq[i:i+k])
    
    return k, core_kmers


def _count_kmers(seq, core_kmers, k):
    """
    Count overlapping occurrences of core k-mers in sequence.
    
    Returns:
    int: Total count of core k-mer occurrences
    """
    count = 0
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i+k]
        if kmer in core_kmers:
            count += 1
    return count

# =========================
# Validation and Utilities
# =========================

def validate_motif(motif, seq_length):
    required_keys = ["Class", "Subtype", "Start", "End", "Length", "Sequence"]
    if not all(key in motif for key in required_keys):
        return False
    if not (1 <= motif["Start"] <= motif["End"] <= seq_length):
        return False
    if len(motif["Sequence"].replace('\n', '')) == 0:
        return False
    return True

def get_basic_stats(seq, motifs=None):
    seq = seq.upper()
    length = len(seq)
    gc = gc_content(seq)
    at = (seq.count('A') + seq.count('T')) / length * 100 if length else 0
    stats = {
        "Length": length,
        "GC%": round(gc, 2),
        "AT%": round(at, 2),
        "A": seq.count('A'),
        "T": seq.count('T'),
        "G": seq.count('G'),
        "C": seq.count('C'),
    }
    if motifs is not None:
        covered = set()
        for m in motifs:
            covered.update(range(m['Start'], m['End']))
        coverage_pct = (len(covered) / length * 100) if length else 0
        stats["Motif Coverage %"] = round(coverage_pct, 2)
    return stats

# =========================
# Main aggregator function
# =========================

def all_motifs(seq, nonoverlap=False, report_hotspots=False, sequence_name="Sequence"):
    """
    Main function to detect all non-B DNA motifs in a sequence.
    Now imports from specialized modules.
    """
    if not seq or not re.match("^[ATGC]+$", seq, re.IGNORECASE):
        return []
    seq = seq.upper()
    
    # Import detection functions from specialized modules
    from .curved_dna import find_curved_DNA
    from .hairpin_cruciform import find_cruciform
    from .slipped_dna import find_slipped_dna
    from .triplex_dna import find_hdna, find_sticky_dna
    from .r_loop import find_rlfs
    from .zdna_egz import find_zdna, find_egz_motif
    from .imotif_ac import find_imotif, find_ac_motifs
    from .g4_related import (find_gquadruplex, find_relaxed_gquadruplex, 
                             find_bulged_gquadruplex, find_bipartite_gquadruplex,
                             find_multimeric_gquadruplex, find_imperfect_gquadruplex,
                             find_gtriplex)
    from .hybrid import find_hybrids
    from .cluster import find_hotspots, select_best_nonoverlapping_motifs
    
    motif_list = (
        find_sticky_dna(seq) +
        find_curved_DNA(seq) +
        find_zdna(seq) +
        find_egz_motif(seq) +
        find_slipped_dna(seq) +
        find_rlfs(seq) +
        find_cruciform(seq) +
        find_hdna(seq) +
        find_gtriplex(seq) +
        find_gquadruplex(seq) +
        find_relaxed_gquadruplex(seq) +
        find_bulged_gquadruplex(seq) +
        find_bipartite_gquadruplex(seq) +
        find_multimeric_gquadruplex(seq) +
        find_imperfect_gquadruplex(seq) +
        find_imotif(seq) +
        find_ac_motifs(seq)
    )
    # Validate and standardize fields
    motif_list = [m for m in motif_list if validate_motif(m, len(seq))]
    # Add hybrids
    motif_list += find_hybrids(motif_list, seq)
    # De-overlap per class if asked
    if nonoverlap:
        motif_list = select_best_nonoverlapping_motifs(motif_list)
    # Hotspots appended if asked
    if report_hotspots:
        motif_list += find_hotspots(motif_list, len(seq))
    # Add Sequence Name and ensure ordered keys exist
    for m in motif_list:
        m["Sequence Name"] = sequence_name
        # Ensure mandatory ordered fields exist and not missing
        if "Arms/Repeat Unit/Copies" not in m:
            m["Arms/Repeat Unit/Copies"] = ""
        if "Spacer" not in m:
            m["Spacer"] = ""
        if "Score" in m:
            try:
                m["Score"] = float(m["Score"])
            except Exception:
                pass
    return motif_list

# =========================
# Utility: formatted output rows in the exact requested order
# =========================

def format_motif_rows(motifs):
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