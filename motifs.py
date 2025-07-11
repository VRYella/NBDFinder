import re
from utils import wrap, gc_content, reverse_complement, g4hunter_score, zseeker_score

def overlapping_finditer(pattern, seq):
    regex = re.compile(pattern)
    for m in regex.finditer(seq):
        yield m

def safe_motif_output(seq, pattern, motif_class, subtype, score_func=None, score_method="None"):
    results = []
    for m in overlapping_finditer(pattern, seq):
        match_seq = m.group(0)
        if not match_seq:
            continue
        start = m.start() + 1
        end = m.end()
        if end <= start:
            continue
        score = score_func(match_seq) if score_func else 0
        results.append(dict(
            Class=motif_class,
            Subtype=subtype,
            Start=start,
            End=end,
            Length=end - start,
            Sequence=wrap(match_seq),
            ScoreMethod=score_method,
            Score=f"{score:.2f}" if isinstance(score, float) else str(score)
        ))
    return results

def find_gquadruplex(seq):
    return safe_motif_output(seq, r"(G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,})",
                             "Quadruplex", "Canonical_G-Quadruplex", g4hunter_score, "G4Hunter")

def find_relaxed_gquadruplex(seq):
    return safe_motif_output(seq, r"(G{3,}[ATGC]{1,12}G{3,}[ATGC]{1,12}G{3,}[ATGC]{1,12}G{3,})",
                             "Quadruplex", "Relaxed_G-Quadruplex", g4hunter_score, "G4Hunter")

def find_bulged_gquadruplex(seq):
    return safe_motif_output(seq, r"(G{3,}[ATGC]{0,3}G{3,}[ATGC]{0,3}G{3,}[ATGC]{0,3}G{3,})",
                             "Quadruplex", "Bulged_G-Quadruplex", g4hunter_score, "G4Hunter (bulge)")

def find_imotif(seq):
    def score_imotif(s): return -g4hunter_score(s.replace('C', 'G'))
    return safe_motif_output(seq, r"(C{3,}[ATGC]{1,7}C{3,}[ATGC]{1,7}C{3,}[ATGC]{1,7}C{3,})",
                             "Quadruplex", "i-Motif", score_imotif, "G4Hunter")

def find_gtriplex(seq):
    return safe_motif_output(seq, r"(G{3,}[ATGC]{1,7}G{3,}[ATGC]{1,7}G{3,})",
                             "Triplex", "G-Triplex", g4hunter_score, "G4Hunter")

def find_bipartite_gquadruplex(seq):
    return safe_motif_output(seq, r"(G{3,}[ATGC]{1,30}G{3,}[ATGC]{1,30}G{3,}[ATGC]{1,30}G{3,})",
                             "Quadruplex", "Bipartite_G-Quadruplex", g4hunter_score, "G4Hunter")

def find_multimeric_gquadruplex(seq):
    return safe_motif_output(seq, r"(G{3,}(?:[ATGC]{1,12}G{3,}){4,})",
                             "Quadruplex", "Multimeric_G-Quadruplex", g4hunter_score, "G4Hunter")

def find_zdna(seq):
    return safe_motif_output(seq, r"((?:CG){6,})",
                             "Z-DNA", "CG_Repeat", zseeker_score, "ZSeeker")

def find_hdna(seq):
    return safe_motif_output(seq, r"(T{3,}[ATGC]{1,7}A{3,})", "H-DNA", "T-A")

def find_sticky_dna(seq):
    return safe_motif_output(seq, r"(CTGCTGCTGCTG)", "Sticky_DNA", "CTG")

def find_slipped_dna(seq):
    return safe_motif_output(seq, r"((?:AT){6,})", "Slipped_DNA", "AT_Slippage")

def find_cruciform(seq):
    return safe_motif_output(seq, r"(A{4,}TTTT)", "Cruciform", "A-T")

def find_bent_dna(seq):
    return safe_motif_output(seq, r"(A{6,7}|T{6,7})", "Bent_DNA", "Poly-A/T")

def find_apr(seq):
    return safe_motif_output(seq, r"((?:AAATT){2,})", "A-Phased_Repeat", "APR")

def find_mirror_repeat(seq):
    return safe_motif_output(seq, r"(ATCGCGAT)", "Mirror_Repeat", "ATCGCGAT")

def find_polyG(seq):
    return safe_motif_output(seq, r"(G{6,})", "Direct_Repeat", "Poly-G")

def find_local_bent(seq):
    return safe_motif_output(seq, r"(A{6,7}|T{6,7})", "Bent_DNA", "Poly-A/T")

def all_motifs(seq):
    motif_funcs = [
        find_gquadruplex, find_relaxed_gquadruplex, find_bulged_gquadruplex,
        find_imotif, find_gtriplex, find_bipartite_gquadruplex, find_multimeric_gquadruplex,
        find_zdna, find_hdna, find_sticky_dna, find_slipped_dna, find_cruciform,
        find_bent_dna, find_apr, find_mirror_repeat, find_polyG, find_local_bent
    ]
    all_hits = []
    for func in motif_funcs:
        all_hits.extend(func(seq))
    return all_hits
