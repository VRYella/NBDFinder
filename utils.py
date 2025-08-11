# ======================== #
#  Nucleotide sequence utilities
# ======================== #
import re
import numpy as np

# Parse FASTA format and clean sequence
def parse_fasta(fasta_str: str) -> str:
    return "".join([line.strip() for line in fasta_str.split('\n') if not line.startswith(">")]).upper().replace(" ", "").replace("U", "T")

# Wrap sequence for display
def wrap(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))

# GC content (%)
def gc_content(seq: str) -> float:
    gc = seq.count('G') + seq.count('C')
    return (gc / max(1, len(seq))) * 100 if seq else 0

# Reverse complement
def reverse_complement(seq: str) -> str:
    return seq.translate(str.maketrans("ATGC", "TACG"))[::-1]

# Palindrome check
def is_palindrome(seq: str) -> bool:
    return seq == reverse_complement(seq)

# Overlapping regex search
def overlapping_finditer(pattern, seq):
    regex = re.compile(pattern, re.IGNORECASE)
    pos = 0
    while pos < len(seq):
        m = regex.search(seq, pos)
        if not m: break
        yield m
        pos = m.start() + 1

# Find polyA/polyT tracts
def find_polyA_polyT_tracts(seq: str, min_len: int = 7) -> list:
    results = []
    i = 0; n = len(seq)
    while i < n:
        if seq[i] == 'A' or seq[i] == 'T':
            start = i
            while i < n and seq[i] == seq[start]: i += 1
            if i - start >= min_len:
                results.append((start, i-1, seq[start:i]))
        else: i += 1
    return results

# Purine/pyrimidine fraction
def purine_fraction(seq): return (seq.count('A') + seq.count('G')) / max(1, len(seq))
def pyrimidine_fraction(seq): return (seq.count('C') + seq.count('T')) / max(1, len(seq))

# G4Hunter score (G/C balance)
def g4hunter_score(seq):
    scores = [1 if c == 'G' else -1 if c == 'C' else 0 for c in seq.upper()]
    return np.mean(scores) if scores else 0

# Validate motif data structure
def validate_motif(motif, seq_length):
    required_keys = ["Class", "Subtype", "Start", "End", "Length", "Sequence"]
    if not all(key in motif for key in required_keys): return False
    if not (1 <= motif["Start"] <= motif["End"] <= seq_length): return False
    if len(motif["Sequence"].replace('\n', '')) == 0: return False
    return True
