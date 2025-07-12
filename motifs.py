import re
from utils import wrap, gc_content, reverse_complement, g4hunter_score, zseeker_score

# Utility function for overlapping pattern matching
def overlapping_finditer(pattern, seq):
    """Find all overlapping matches of a pattern in sequence"""
    regex = re.compile(pattern)
    for m in regex.finditer(seq):
        yield m

# Main function that collects all motifs
def all_motifs(seq):
    """
    Run all motif detection functions and return combined results
    Args:
        seq: DNA sequence to analyze
    Returns:
        List of dictionaries with motif information
    """
    if not seq or not re.match("^[ATGC]+$", seq.upper()):
        return []
    
    return (
        find_apr(seq) + find_bent_dna(seq) + find_zdna(seq) +
        find_slipped_dna(seq) + find_rlfs(seq) + find_cruciform(seq) +
        find_hdna(seq) + find_gtriplex(seq) + find_gquadruplex(seq) +
        find_relaxed_gquadruplex(seq) + find_bulged_gquadruplex(seq) +
        find_bipartite_gquadruplex(seq) + find_multimeric_gquadruplex(seq) +
        find_imotif(seq) + find_hybrids(seq)
    )

# 1. CURVED / FLEXIBLE / LOCAL BENT DNA
def find_apr(seq):
    """Find A-Phased Repeats (APRs) associated with DNA curvature"""
    pattern = r"(?=(?:AAATT){2,})"
    return [
        {
            "Class": "Curved_DNA",
            "Subtype": "A-Phased_Repeat", 
            "Start": m.start()+1,
            "End": m.end(),
            "Length": len(m.group(0)),
            "Sequence": wrap(m.group(0)),
            "ScoreMethod": "NBST",
            "Score": "0"
        }
        for m in overlapping_finditer(pattern, seq)
    ]

def find_bent_dna(seq):
    """Find poly-A/T tracts associated with DNA bending"""
    pattern = r"(?=(A{6,7}|T{6,7}))"
    return [
        {
            "Class": "Bent_DNA",
            "Subtype": "Poly-A/T", 
            "Start": m.start()+1,
            "End": m.end(),
            "Length": len(m.group(1)),
            "Sequence": wrap(m.group(1)),
            "ScoreMethod": "NBST",
            "Score": "0"
        }
        for m in overlapping_finditer(pattern, seq)
    ]

# 2. Z-DNA
def find_zdna(seq):
    """Find Z-DNA forming alternating purine-pyrimidine repeats"""
    pattern = r"(?=((?:CG|GC|GT|TG|AC|CA){6,}))"
    return [
        {
            "Class": "Z-DNA",
            "Subtype": "Purine-Pyrimidine Repeats", 
            "Start": m.start()+1,
            "End": m.end(),
            "Length": len(m.group(1)),
            "Sequence": wrap(m.group(1)),
            "ScoreMethod": "Z-Seeker",
            "Score": f"{zseeker_score(m.group(1)):.2f}"
        }
        for m in overlapping_finditer(pattern, seq)
    ]

# 3. SLIPPED DNA
def find_slipped_dna(seq):
    """Find slipped DNA structures (AT repeats)"""
    pattern = r"(?=((?:AT){6,}))"
    return [
        {
            "Class": "Slipped_DNA",
            "Subtype": "AT Repeats", 
            "Start": m.start()+1,
            "End": m.end(),
            "Length": len(m.group(1)),
            "Sequence": wrap(m.group(1)),
            "ScoreMethod": "NBST",
            "Score": "0"
        }
        for m in overlapping_finditer(pattern, seq)
    ]

# 4. R-LOOPS (QmRLFS-finder logic)
def find_rlfs(seq):
    """Find R-loop forming sequences"""
    results = []
    seq = seq.upper()
    n = len(seq)
    
    # R-loop initiation zones (RIZ)
    model_m1 = r"(?=(G{3,}[ATGC]{1,10}G{3,}[ATGC]{1,10}G{4,}))"
    model_m2 = r"(?=(G{3,}[ATGC]{1,10}G{4,}))"
    
    riz_hits = []
    for pattern, model in [(model_m1, "m1"), (model_m2, "m2")]:
        for match in re.finditer(pattern, seq):
            riz_seq = match.group(1)
            start = match.start(1)
            end = match.end(1)
            if gc_content(riz_seq) >= 50.0:
                riz_hits.append((start, end, model, riz_seq))
    
    # Find R-loop extension zones (REZ)
    for riz_start, riz_end, model, riz_seq in riz_hits:
        for linker_len in range(0, 51):
            rez_start = riz_end + linker_len
            for rez_len in range(100, 2001):
                rez_end = rez_start + rez_len
                if rez_end > n:
                    break
                rez_seq = seq[rez_start:rez_end]
                if gc_content(rez_seq) >= 40.0:
                    results.append({
                        "Class": "R-Loop",
                        "Subtype": f"RLFS_{model}",
                        "Start": riz_start+1,
                        "End": rez_end,
                        "Length": rez_end - riz_start,
                        "Sequence": wrap(seq[riz_start:rez_end]),
                        "ScoreMethod": "QmRLFS",
                        "Score": "1"
                    })
                    break
            if any(r['Start'] == riz_start+1 for r in results):
                break
    return results

# 5. CRUCIFORM DNA / HAIRPIN
def find_cruciform(seq):
    """Find cruciform DNA structures (A-T palindromes)"""
    pattern = r"(?=(A{4,}TTTT))"
    return [
        {
            "Class": "Cruciform_DNA",
            "Subtype": "A-T Palindrome", 
            "Start": m.start()+1,
            "End": m.end(),
            "Length": len(m.group(1)),
            "Sequence": wrap(m.group(1)),
            "ScoreMethod": "NBST",
            "Score": "0"
        }
        for m in overlapping_finditer(pattern, seq)
    ]

# 6. TRIPLEX DNA / H-DNA
def find_hdna(seq):
    """Find H-DNA (triplex DNA) structures"""
    pattern = r"(?=(T{3,}[ATGC]{1,7}A{3,}))"
    return [
        {
            "Class": "Triplex_DNA",
            "Subtype": "H-DNA T-A", 
            "Start": m.start()+1,
            "End": m.end(),
            "Length": len(m.group(1)),
            "Sequence": wrap(m.group(1)),
            "ScoreMethod": "NBST",
            "Score": "0"
        }
        for m in overlapping_finditer(pattern, seq)
    ]

# 7. G-TRIPLEX
def find_gtriplex(seq):
    """Find G-triplex structures"""
    pattern = r"(?=(G{3,}(?:[ATGC]{0,7}G{3,}){2}))"
    return [
        {
            "Class": "G-Triplex",
            "Subtype": "Three G-runs", 
            "Start": m.start()+1,
            "End": m.end(),
            "Length": len(m.group(0)),
            "Sequence": wrap(m.group(0)),
            "ScoreMethod": "G4Hunter",
            "Score": f"{g4hunter_score(m.group(0)):.2f}"
        }
        for m in overlapping_finditer(pattern, seq)
    ]

# 8. G-QUADRUPLEX VARIANTS
def find_gquadruplex(seq):
    """Find canonical G-quadruplex structures"""
    pattern = r"(G{3,}(?:[ATGC]{1,7}G{3,}){3})"
    return [
        {
            "Class": "Quadruplex",
            "Subtype": "Canonical_G-Quadruplex", 
            "Start": m.start()+1,
            "End": m.end(),
            "Length": len(m.group(0)),
            "Sequence": wrap(m.group(0)),
            "ScoreMethod": "G4Hunter",
            "Score": f"{g4hunter_score(m.group(0)):.2f}"
        }
        for m in overlapping_finditer(pattern, seq)
    ]

def find_relaxed_gquadruplex(seq):
    """Find relaxed G-quadruplex structures (longer loops)"""
    pattern = r"(G{3,}(?:[ATGC]{1,12}G{3,}){3})"
    return [
        {
            "Class": "Quadruplex",
            "Subtype": "Relaxed_G-Quadruplex", 
            "Start": m.start()+1,
            "End": m.end(),
            "Length": len(m.group(0)),
            "Sequence": wrap(m.group(0)),
            "ScoreMethod": "G4Hunter",
            "Score": f"{g4hunter_score(m.group(0)):.2f}"
        }
        for m in overlapping_finditer(pattern, seq)
    ]

def find_bulged_gquadruplex(seq):
    """Find bulged G-quadruplex structures"""
    pattern = r"(G{3,}[ATGC]{0,3}G{3,}[ATGC]{0,3}G{3,}[ATGC]{0,3}G{3,})"
    return [
        {
            "Class": "Quadruplex",
            "Subtype": "Bulged_G-Quadruplex", 
            "Start": m.start()+1,
            "End": m.end(),
            "Length": len(m.group(0)),
            "Sequence": wrap(m.group(0)),
            "ScoreMethod": "G4Hunter (bulge)",
            "Score": f"{g4hunter_score(m.group(0)):.2f}"
        }
        for m in overlapping_finditer(pattern, seq)
    ]

def find_bipartite_gquadruplex(seq):
    """Find bipartite G-quadruplex structures"""
    pattern = r"(?=(G{3,}(?:[ATGC]{0,30}G{3,}){3}))"
    return [
        {
            "Class": "Quadruplex",
            "Subtype": "Bipartite_G-Quadruplex", 
            "Start": m.start()+1,
            "End": m.end(),
            "Length": len(m.group(0)),
            "Sequence": wrap(m.group(0)),
            "ScoreMethod": "G4Hunter",
            "Score": f"{g4hunter_score(m.group(0)):.2f}"
        }
        for m in overlapping_finditer(pattern, seq)
    ]

def find_multimeric_gquadruplex(seq):
    """Find multimeric G-quadruplex structures"""
    pattern = r"(?=((G{3,}(?:[ATGC]{0,12}G{3,}){4,})))"
    return [
        {
            "Class": "Quadruplex",
            "Subtype": "Multimeric_G-Quadruplex", 
            "Start": m.start()+1,
            "End": m.end(),
            "Length": len(m.group(1)),
            "Sequence": wrap(m.group(1)),
            "ScoreMethod": "G4Hunter",
            "Score": f"{g4hunter_score(m.group(1)):.2f}"
        }
        for m in overlapping_finditer(pattern, seq)
    ]

# 9. i-MOTIF
def find_imotif(seq):
    """Find i-motif (C-rich) structures"""
    pattern = r"(?=(C{3,}(?:[ATGC]{0,7}C{3,}){3}))"
    return [
        {
            "Class": "i-Motif",
            "Subtype": "C-rich Quadruplex", 
            "Start": m.start()+1,
            "End": m.end(),
            "Length": len(m.group(0)),
            "Sequence": wrap(m.group(0)),
            "ScoreMethod": "G4Hunter",
            "Score": f"{-g4hunter_score(m.group(0).replace('C','G')):.2f}"
        }
        for m in overlapping_finditer(pattern, seq)
    ]

# 10. HYBRID MOTIFS
def find_hybrids(seq):
    """Find hybrid motifs (overlapping G4 and i-motif)"""
    hits = []
    g4 = [(m.start(), m.end()) for m in overlapping_finditer(r"(?=(G{3,}(?:[ATGC]{0,7}G{3,}){3}))", seq)]
    im = [(m.start(), m.end()) for m in overlapping_finditer(r"(?=(C{3,}(?:[ATGC]{0,7}C{3,}){3}))", seq)]
    
    for a_start, a_end in g4:
        for b_start, b_end in im:
            if a_start <= b_end and b_start <= a_end:
                start, end = min(a_start, b_start), max(a_end, b_end)
                hits.append({
                    "Class": "Hybrid",
                    "Subtype": "G4-iMotif", 
                    "Start": start+1,
                    "End": end,
                    "Length": end-start,
                    "Sequence": wrap(seq[start:end]),
                    "ScoreMethod": "Overlap",
                    "Score": "0"
                })
    return hits

# 11. HOTSPOTS / CLUSTERED REGIONS
def find_hotspots(seq, motif_hits, window=100, min_count=3):
    """Find regions with high motif density (hotspots)"""
    n = len(seq)
    hotspots = []
    positions = [(hit['Start'], hit['End']) for hit in motif_hits]
    
    for i in range(1, n - window + 2):
        start, end = i, i + window - 1
        count = sum(mstart <= end and mend >= start for mstart, mend in positions)
        if count >= min_count:
            hotspots.append({
                "RegionStart": start,
                "RegionEnd": end,
                "MotifCount": count
            })
    
    # Remove duplicate regions
    return list({(h['RegionStart'], h['RegionEnd']): h for h in hotspots}.values())
