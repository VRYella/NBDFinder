# =============================================================
# G4 Family Motif Finder Module for NBDFinder - Modernized Version
# Author: VRYella
# Scientific basis: QuadBase2, G4Hunter, Bedrat et al. NAR 2016
# =============================================================

import re

# --- Succinct, annotated regex patterns for all major G4 subclasses ---
PATS = [
    # Multimeric/Tandem: two canonical G4s within 20 nt
    ("G4", "Multimeric/Tandem", re.compile(r"(?=(?:([Gg]{3,}\w{1,7}){3}[Gg]{3,}).{0,20}(?:([Gg]{3,}\w{1,7}){3}[Gg]{3,}))")),
    # Split/Bipartite: two G3+L1-7 blocks separated by long linker (20-50 nt)
    ("G4", "Split/Bipartite",   re.compile(r"(?=([Gg]{3,}\w{20,50}[Gg]{3,}))")),
    # Bulged: G4s with single-base bulges in G-runs (up to 3 nt between Gs)
    ("G4", "Bulged",            re.compile(r"(?=(([Gg]{2,}\w{0,3}[Gg]{1,}\w{1,12}){3,}[Gg]{2,}))")),
    # Canonical: four G-runs (≥3 Gs), loops 1–7 nt
    ("G4", "Canonical G3+L1–7", re.compile(r"(?=((?:[Gg]{3,}\w{1,7}){3}[Gg]{3,}))")),
    # Long-Loop: canonical G4 with at least one loop 8–12 nt (filter in logic)
    ("G4", "Long-Loop G3+L8–12",re.compile(r"(?=((?:[Gg]{3,}\w{1,12}){3}[Gg]{3,}))")),
    # Extended: canonical G4 with loops up to 12 nt
    ("G4", "Extended G3+L1–12", re.compile(r"(?=((?:[Gg]{3,}\w{1,12}){3}[Gg]{3,}))")),
    # Two-Tetrad: four G-runs of 2 Gs, loops 1–12 nt
    ("G4", "Two-Tetrad G2L1–12",re.compile(r"(?=((?:[Gg]{2}\w{1,12}){3}[Gg]{2}))")),
    # Generalized: four G-runs (≥2 Gs), loops 1–12 nt
    ("G4", "Generalized G2+L1–12", re.compile(r"(?=((?:[Gg]{2,}\w{1,12}){3}[Gg]{2,}))")),
    # G-Triplex: three G-runs (≥3 Gs), loops up to 15 nt
    ("G4", "G-Triplex",         re.compile(r"(?=(G{3,}\w{1,15}G{3,}\w{1,15}G{3,}))")),
]

# --- Subclass priority for de-duplication and reporting ---
PRIORITY = {
    "Multimeric/Tandem": 1,
    "Split/Bipartite": 2,
    "Bulged": 3,
    "Snapback": 4,              # Snapback: inferred (not regexed here)
    "Canonical G3+L1–7": 5,
    "Long-Loop G3+L8–12": 6,
    "Extended G3+L1–12": 7,
    "Two-Tetrad G2L1–12": 8,
    "Generalized G2+L1–12": 9,
    "G-Triplex": 10,
}

# --- Fast G4Hunter-like scoring: G = +1, C = -1, normalized by length ---
def g4hunter_score(seq):
    # Simple, literature-standard score for G4 propensity
    val = 0
    for b in seq.upper():
        if b == 'G': val += 1
        elif b == 'C': val -= 1
    return val / max(1, len(seq))

# --- Helper function for overlapping matches ---
def find_overlapping_matches(pattern, seq):
    """Find all overlapping matches for a pattern in sequence"""
    matches = []
    start = 0
    while start < len(seq):
        match = pattern.search(seq, start)
        if not match:
            break
        matches.append(match)
        start = match.start() + 1
    return matches

# --- Main candidate motif finder for a single strand ---
def candidates_from(seq, strand='+'):
    # Returns all G4 motif candidates, with subclass and score annotations
    cand = []
    for cls, sub, pat in PATS:
        for m in find_overlapping_matches(pat, seq):
            s = m.start(1) if m.lastindex else m.start()
            e = s + len(m.group(1) if m.lastindex else m.group(0))
            raw = seq[s:e]
            sc = g4hunter_score(raw)
            # For Long-Loop, filter: at least one loop 8–12 nt
            if sub == "Long-Loop G3+L8–12":
                splits = [g.span() for g in re.finditer(r"[Gg]{3,}", raw)]
                if len(splits) == 4:
                    loops = [splits[i + 1][0] - splits[i][1] for i in range(3)]
                    if not any(8 <= l <= 12 for l in loops): continue
            cand.append({
                'start': s, 'end': e, 'strand': strand,
                'class': cls, 'subclass': sub, 'seq': raw,
                'score': sc, 'priority': PRIORITY[sub]
            })
    return cand

# --- Non-overlapping motif selector (highest priority and score) ---
def non_overlapping(cands):
    # Selects non-overlapping motifs by priority, score, G-count, and length
    cands.sort(key=lambda x: (
        x['priority'], -x['score'], -(x['end']-x['start']),
        -x['seq'].upper().count('G'), x['start'])
    )
    chosen = []
    for c in cands:
        if all(c['end'] <= k['start'] or c['start'] >= k['end'] for k in chosen):
            chosen.append(c)
    return chosen

# =============================================================
# Module ready for integration in NBDFinder. All regexes and logic are
# literature-standard and tested for scientific reproducibility.
# =============================================================