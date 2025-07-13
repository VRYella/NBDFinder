import re
import numpy as np
from typing import List, Dict, Tuple
from collections import defaultdict, Counter
import random

try:
    from scipy.stats import percentileofscore
except ImportError:
    # Fallback implementation if scipy is not available
    def percentileofscore(a, score, kind='rank'):
        """Simplified percentile calculation without scipy"""
        a = np.asarray(a)
        if len(a) == 0:
            return 0.0
        if kind == 'rank':
            return (sum(a <= score) / len(a) * 100
        elif kind == 'strict':
            return (sum(a < score) / len(a)) * 100
        elif kind == 'weak':
            return (sum(a <= score) / len(a)) * 100
        elif kind == 'mean':
            return (sum(a < score) + sum(a <= score)) / (2 * len(a)) * 100
        else:
            raise ValueError("kind must be 'rank', 'strict', 'weak' or 'mean'")

def parse_fasta(fasta_str: str) -> str:
    """Parse FASTA string to DNA sequence"""
    lines = [line.strip() for line in fasta_str.split('\n') if not line.startswith(">")]
    return "".join(lines).upper().replace(" ", "").replace("U", "T")

def wrap(seq: str, width: int = 60) -> str:
    """Format sequence with line breaks"""
    return "\n".join([seq[i:i+width] for i in range(0, len(seq), width)])

def gc_content(seq: str) -> float:
    """Calculate GC content percentage"""
    seq = seq.upper()
    return 100 * (seq.count('G') + seq.count('C')) / max(1, len(seq))

def reverse_complement(seq: str) -> str:
    """Generate reverse complement"""
    comp = str.maketrans("ATGC", "TACG")
    return seq.translate(comp)[::-1]

def is_palindrome(seq: str) -> bool:
    """Check for perfect palindromes"""
    return seq == reverse_complement(seq)

def calculate_tm(seq: str) -> float:
    """Calculate DNA melting temperature"""
    if len(seq) < 14:
        return 2*(seq.count('A') + seq.count('T')) + 4*(seq.count('G') + seq.count('C'))
    return 64.9 + 41*(seq.count('G') + seq.count('C') - 16.4)/len(seq)

def shuffle_sequence(seq: str) -> str:
    """Create randomized sequence preserving composition"""
    return ''.join(random.sample(seq, len(seq)))

def kmer_conservation(seq: str, k: int = 6, n_shuffles: int = 1000) -> Dict[str, Tuple[float, float]]:
    """
    Calculate k-mer conservation scores
    Returns: {kmer: (log2_ratio, p_value)}
    """
    # Count observed kmers
    kmer_counts = Counter(seq[i:i+k] for i in range(len(seq)-k+1))
    total_kmers = len(seq) - k + 1
    
    # Generate null distribution
    null_counts = defaultdict(list)
    for _ in range(n_shuffles):
        shuffled = shuffle_sequence(seq)
        for kmer, count in Counter(shuffled[i:i+k] for i in range(len(shuffled)-k+1)).items():
            null_counts[kmer].append(count)
    
    # Calculate conservation metrics
    results = {}
    for kmer, observed in kmer_counts.items():
        expected = (1/4)**k * total_kmers
        log2_ratio = np.log2((observed + 1e-6)/(expected + 1e-6))  # Pseudocounts
        
        # Calculate p-value from null distribution
        if kmer in null_counts:
            p_value = 1 - percentileofscore(null_counts[kmer], observed)/100
        else:
            p_value = 1.0
            
        results[kmer] = (log2_ratio, p_value)
    
    return results

def motif_conservation(motif_seq: str, conservation_scores: Dict) -> float:
    """Calculate average conservation for a motif"""
    scores = []
    for i in range(len(motif_seq)-5):
        kmer = motif_seq[i:i+6]
        if kmer in conservation_scores:
            scores.append(conservation_scores[kmer][0])
    return np.mean(scores) if scores else 0.0
