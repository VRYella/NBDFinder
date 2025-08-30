"""
Category 8: G4-related Detection Module - All Known Motif Types (Single Window)
===============================================================================
Implements detection of all G-quadruplex related structures using a strictly literature-aligned, reproducible framework.
Includes: canonical, relaxed, bulged, bipartite, multimeric, imperfect (G2), triplex.
References:
- Bedrat et al. (2016) Nucleic Acids Res 44:e70 (G4Hunter)
- Chambers et al. (2015) Nat Biotechnol
- QuadBase2, G4RNA screener, pqsfinder, and computational reviews[1][5].
Author: Dr. Venkata Rajesh Yella, 2025. All scientific upgrade annotations inline.
"""

import re
from math import isfinite
from .shared_utils import wrap, calculate_conservation_score, overlapping_finditer, calculate_structural_factor

def _g4hunter_base_scores(seq):
    """
    Per-base G4Hunter score, see Bedrat et al.
    """
    n = len(seq)
    s = *n
    g_run = 0
    for i, ch in enumerate(seq):
        if ch == 'G':
            g_run += 1
        else:
            g_run = 0
        if ch == 'G':
            s[i] = min(4, g_run)
    c_run = 0
    for i, ch in enumerate(seq):
        if ch == 'C':
            c_run += 1
        else:
            c_run = 0
        if ch == 'C':
            s[i] = -min(4, c_run)
        elif ch in ('A','T') and s[i]==0:
            s[i] = 0
    return s

def g4hunter_score_single_window(seq):
    """
    Return global mean G4Hunter score. No windowing (single region)!
    """
    if not seq:
        return 0.0
    base_scores = _g4hunter_base_scores(seq)
    mean_score = sum(base_scores) / float(len(base_scores)) if base_scores else 0.0
    return float(mean_score)

def _find_g_runs(seq):
    return [m.span() for m in re.finditer(r"G+", seq)]

def _loop_lengths_from_runs(runs):
    return [runs[i+1]-runs[i][6] for i in range(len(runs)-1)] if len(runs)>1 else []

# Canonical G4: Four G-runs ≥3nt, loops 1-7nt
def _is_canonical_G3_L1_7(runs, loops):
    if len(runs)<4 or len(loops)<3: return False
    for k in range(4):
        if (runs[k][6]-runs[k])<3: return False
    return all(1<=L<=7 for L in loops[:3])

# Relaxed (long-loop, extended): Four G-runs ≥3nt, at least one loop 8-12nt, others ≤12nt
def _is_longloop_G3_L8_12(runs, loops):
    if len(runs)<4 or len(loops)<3: return False
    if any((runs[k][6]-runs[k])<3 for k in range(4)): return False
    if not all(1<=L<=12 for L in loops[:3]): return False
    return any(8<=L<=12 for L in loops[:3]) and not _is_canonical_G3_L1_7(runs, loops)

def _is_extended_G3_L1_12(runs, loops):
    if len(runs)<4 or len(loops)<3: return False
    if any((runs[k][6]-runs[k])<3 for k in range(4)): return False
    if not all(1<=L<=12 for L in loops[:3]): return False
    return not _is_canonical_G3_L1_7(runs, loops) and not _is_longloop_G3_L8_12(runs, loops)

# Imperfect: Four G-runs of 2nt, loops ≤12nt
def _is_two_tetrad_G2_L1_12(runs, loops):
    if len(runs)<4 or len(loops)<3: return False
    if any((runs[k][6]-runs[k])!=2 for k in range(4)): return False
    return all(1<=L<=12 for L in loops[:3])

def _is_generalized_G2plus_L1_12(runs, loops):
    if len(runs)<4 or len(loops)<3: return False
    if any((runs[k][6]-runs[k])<2 for k in range(4)): return False
    if not all(1<=L<=12 for L in loops[:3]): return False
    return not _is_two_tetrad_G2_L1_12(runs, loops)

# Bulged: One or more G3+ runs interrupted by a single base (GGGA, GGGT, etc.)
def _is_bulged(seq):
    runs = _find_g_runs(seq)
    loops = _loop_lengths_from_runs(runs)
    if len(runs)<4 or len(loops)<3: return False
    has_bulge = False
    if re.search(r'GGG[ACT]GGG', seq):
        has_bulge = True
    elif re.search(r'GG[ACT]GG', seq) and not re.search(r'GGG[ACT]GGG', seq):
        has_bulge = True
    long_runs = sum(1 for r in runs[:4] if (r[6]-r)>=3)
    return has_bulge and long_runs >= 2

# Triplex: Three G-runs ≥3nt, loops ≤15nt
def _is_triplex(runs, loops):
    if len(runs)<3 or len(loops)<2: return False
    if any((runs[k][6]-runs[k])<3 for k in range(3)): return False
    return all(1<=L<=15 for L in loops[:2])

# Extraction: All motif types!
def _extract_candidates(seq, threshold=1.2):
    candidates = []

    # Canonical/relaxed/extended/bulged/imperfect/generalized (G4 four-run)
    for m in re.finditer(r"(G+)(.{1,15})G+(.{1,15})G+(.{1,15})G+", seq):
        start, end = m.start(), m.end()
        sub = seq[start:end]
        runs = _find_g_runs(sub)
        if len(runs)<4: continue
        loops = _loop_lengths_from_runs(runs)
        g4h = g4hunter_score_single_window(sub)
        # Loop penalty/bonus sequence-only
        penalty = 1.0
        if any(L == 1 for L in loops[:3]): penalty *= 0.75
        elif all(2 <= L <= 4 for L in loops[:3]): penalty *= 1.10
        elif any(L > 12 for L in loops[:3]): penalty *= 0.5
        g4h *= penalty
        if g4h < threshold: continue
        subtype = None
        if _is_bulged(sub): subtype = "Bulged"
        elif _is_canonical_G3_L1_7(runs, loops): subtype = "Canonical G3+L1–7"
        elif _is_longloop_G3_L8_12(runs, loops): subtype = "Long-Loop G3+L8–12"
        elif _is_extended_G3_L1_12(runs, loops): subtype = "Extended G3+L1–12"
        elif _is_two_tetrad_G2_L1_12(runs, loops): subtype = "Two-Tetrad G2L1–12"
        elif _is_generalized_G2plus_L1_12(runs, loops): subtype = "Generalized G2+L1–12"
        if subtype:
            candidates.append({
                "start": start, "end": end, "seq": sub, "subclass": subtype,
                "g4h": g4h, "runs": runs, "loops": loops
            })

    # Bipartite: Two G4 cores separated by long linkers
    for m in re.finditer(r"(G{3,})(.{1,12})G{3,}(.{15,50})G{3,})(.{1,12})G{3,}", seq):
        start, end = m.start(), m.end()
        sub = seq[start:end]
        runs = _find_g_runs(sub)
        loops = _loop_lengths_from_runs(runs)
        g4h = g4hunter_score_single_window(sub)
        # Loop adjustment (bipartite always long linker penalty)
        if g4h < threshold: continue
        candidates.append({
            "start": start, "end": end, "seq": sub, "subclass": "Bipartite/Split",
            "g4h": g4h, "runs": runs, "loops": loops
        })

    # Multimeric/Tandem: Multiple G4 motifs with short linkers
    last_end = None
    multimer_buffer = []
    for m in re.finditer(r"(G{3,})(.{1,7})G{3,}(.{1,7})G{3,}(.{1,7})G{3,}", seq):
        start, end = m.start(), m.end()
        if last_end and (start - last_end) <= 20:
            multimer_buffer[-1]['end'] = end
            multimer_buffer[-1]['seq'] += seq[last_end:end]
        else:
            multimer_buffer.append({'start': start, 'end': end, 'seq': seq[start:end]})
        last_end = end
    for mm in multimer_buffer:
        g4h = g4hunter_score_single_window(mm['seq'])
        if g4h >= threshold:
            candidates.append({
                "start": mm['start'],
                "end": mm['end'],
                "seq": mm['seq'],
                "subclass": "Multimeric/Tandem",
                "g4h": g4h,
                "runs": [],
                "loops": []
            })

    # Triplex: Three G runs (≥3), loops ≤15
    for m in re.finditer(r"(G{3,})(.{1,15})G{3,}(.{1,15})G{3,}", seq):
        start, end = m.start(), m.end()
        sub = seq[start:end]
        runs = _find_g_runs(sub)
        loops = _loop_lengths_from_runs(runs)
        if not _is_triplex(runs, loops): continue
        g4h = g4hunter_score_single_window(sub)
        if g4h < max(1.5, threshold): continue
        candidates.append({
            "start": start, "end": end, "seq": sub, "subclass": "G-Triplex",
            "g4h": g4h, "runs": runs, "loops": loops
        })

    return candidates
# ==============================
# End of fully expanded single-window, sequence-only G4 detection module
# ==============================
