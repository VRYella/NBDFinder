import re
import numpy as np

# --- FASTA Parsing ---
def parse_fasta(fasta_str: str) -> str:
    lines = fasta_str.strip().splitlines()
    seq = [line.strip() for line in lines if not line.startswith(">")]
    return "".join(seq).upper().replace(" ", "").replace("U", "T")

# --- Sequence Wrapping ---
def wrap(seq: str, width=60) -> str:
    return "\n".join([seq[i:i+width] for i in range(0, len(seq), width)])

# --- GC Content Calculation ---
def gc_content(seq: str) -> float:
    seq = seq.upper()
    g = seq.count('G')
    c = seq.count('C')
    return 100.0 * (g + c) / max(1, len(seq))

# --- Reverse Complement ---
def reverse_complement(seq: str) -> str:
    comp = str.maketrans("ATGC", "TACG")
    return seq.translate(comp)[::-1]

# --- G4Hunter Score ---
def g4hunter_score(seq: str) -> float:
    vals = []
    i = 0
    while i < len(seq):
        s = seq[i]
        if s == 'G':
            run_len = 1
            while i + run_len < len(seq) and seq[i + run_len] == 'G':
                run_len += 1
            score = min(run_len, 4)
            vals.extend([score] * run_len)
            i += run_len
        elif s == 'C':
            run_len = 1
            while i + run_len < len(seq) and seq[i + run_len] == 'C':
                run_len += 1
            score = -min(run_len, 4)
            vals.extend([score] * run_len)
            i += run_len
        else:
            vals.append(0)
            i += 1
    return np.mean(vals) if vals else 0.0

# --- Z-DNA Seeker Score ---
def zseeker_score(seq: str) -> float:
    dinucs = re.findall(r"(GC|CG|GT|TG|AC|CA)", seq.upper())
    return len(dinucs) / max(1, len(seq) / 2) if len(seq) >= 2 else 0.0

# --- Dinucleotide Entropy (optional advanced use) ---
def dinucleotide_entropy(seq: str) -> float:
    from collections import Counter
    seq = seq.upper()
    dinucs = [seq[i:i+2] for i in range(len(seq) - 1)]
    total = len(dinucs)
    if total == 0:
        return 0.0
    freqs = Counter(dinucs)
    probs = [count / total for count in freqs.values()]
    return -sum(p * np.log2(p) for p in probs)
