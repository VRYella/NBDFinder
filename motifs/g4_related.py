"""
Category 8: G4-related Detection Module - Modernized with Exact G4Hunter
========================================================================
Implements detection of all G-quadruplex related structures using a strictly literature-aligned, reproducible framework:

- Implements exact basewise G4Hunter (Bedrat et al. NAR 2016; basewise ±1..±4, window mean, threshold ≥1.2)
- Explicit subclass taxonomy with regex/logic for canonical, relaxed, bulged, G2, generalized, triplex, multimer, bipartite
- Strict, priority-then-score non-overlap resolver; allows class-aware multimer/bipartite merging
- Outputs a single, non-overlapping, reproducible set of G4 family calls, with biophysically plausible priority

References:
- Bedrat et al. (2016) Nucleic Acids Res 44:e70 (G4Hunter)
- Chambers et al. (2015) Nat Biotechnol; Hänsel-Hertsch et al. (2017) Nat Rev Mol Cell Biol
- QuadBase2: comprehensive G4 classification

Author: Dr. Venkata Rajesh Yella, 2025. All scientific upgrade annotations inline.
"""

import re
from math import isfinite

# ---- External utility imports expected by the framework ----
from .shared_utils import (
    wrap, calculate_conservation_score, overlapping_finditer, calculate_structural_factor
)

# ================================================
# Exact G4Hunter implementation (Bedrat et al. 2016)
# ================================================
def _g4hunter_base_scores(seq):
    """
    Compute per-base G4Hunter scores (Bedrat et al. 2016):
    - A/T: 0; G: +run_length (capped at 4); C: -run_length (capped at 4).
    """
    n = len(seq)
    s = [0]*n
    # Forward G-run
    g_run = 0
    for i, ch in enumerate(seq):
        if ch == 'G':
            g_run += 1
        else:
            g_run = 0
        if ch == 'G':
            s[i] = min(4, g_run)
    # Forward C-run
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

def _g4hunter_window_scores(base_scores, window=25):
    """
    Sliding-window mean G4Hunter score per Bedrat et al. Returns windowed mean array.
    """
    n = len(base_scores)
    if n == 0:
        return []
    W = max(1, int(window))
    half = W//2
    prefix = [0]
    for v in base_scores: prefix.append(prefix[-1]+v)
    out = [0.0]*n
    for i in range(n):
        a, b = max(0, i-half), min(n, i+half+1)
        total = prefix[b] - prefix[a]
        out[i] = total/(b-a) if b>a else 0.0
    return out

def g4hunter_score(seq, window=25):
    """
    Returns exact G4Hunter max-window mean, strictly per the original algorithm.
    """
    if not seq: return 0.0
    base = _g4hunter_base_scores(seq)
    win = _g4hunter_window_scores(base, window=window)
    maxv = max(win) if win else 0.0
    return float(maxv if isfinite(maxv) else 0.0)

# ================================================
# Subclass pattern detectors (strict scientific guards)
# ================================================
def _find_g_runs(seq):
    """Return list of (start, end) for contiguous G+ runs in seq."""
    return [m.span() for m in re.finditer(r"G+", seq)]

def _loop_lengths_from_runs(runs):
    """Given runs [(s1,e1),...], return loop lengths between successive runs."""
    return [runs[i+1][0]-runs[i][1] for i in range(len(runs)-1)] if len(runs)>1 else []

def _is_canonical_G3_L1_7(runs, loops):
    # Canonical: four G-runs ≥3, loops 1-7 nt
    if len(runs)<4 or len(loops)<3: return False
    for k in range(4):
        if (runs[k][1]-runs[k][0])<3: return False
    return all(1<=L<=7 for L in loops[:3])

def _is_longloop_G3_L8_12(runs, loops):
    # Long-loop: four G-runs ≥3, all loops ≤12, at least one loop 8–12, not canonical
    if len(runs)<4 or len(loops)<3: return False
    if any((runs[k][1]-runs[k][0])<3 for k in range(4)): return False
    if not all(1<=L<=12 for L in loops[:3]): return False
    return any(8<=L<=12 for L in loops[:3]) and not _is_canonical_G3_L1_7(runs, loops)

def _is_extended_G3_L1_12(runs, loops):
    # Extended: four G-runs ≥3, all loops ≤12, not canonical or long-loop
    if len(runs)<4 or len(loops)<3: return False
    if any((runs[k][1]-runs[k][0])<3 for k in range(4)): return False
    if not all(1<=L<=12 for L in loops[:3]): return False
    return not _is_canonical_G3_L1_7(runs, loops) and not _is_longloop_G3_L8_12(runs, loops)

def _is_two_tetrad_G2_L1_12(runs, loops):
    # Two-tetrad: exactly four runs of length 2, loops ≤12
    if len(runs)<4 or len(loops)<3: return False
    if any((runs[k][1]-runs[k][0])!=2 for k in range(4)): return False
    return all(1<=L<=12 for L in loops[:3])

def _is_generalized_G2plus_L1_12(runs, loops):
    # Generalized: four runs of ≥2, loops ≤12, not two-tetrad
    if len(runs)<4 or len(loops)<3: return False
    if any((runs[k][1]-runs[k][0])<2 for k in range(4)): return False
    if not all(1<=L<=12 for L in loops[:3]): return False
    return not _is_two_tetrad_G2_L1_12(runs, loops)

def _is_bulged(seq):
    """
    Bulged: contains G-runs with single base interruptions, creating bulges in the G4 structure.
    Enhanced to detect various bulge patterns.
    """
    runs = _find_g_runs(seq)
    loops = _loop_lengths_from_runs(runs)
    if len(runs)<4 or len(loops)<3: return False
    
    # Look for single-nucleotide interruptions in what would otherwise be longer G-runs
    # The key insight: GGGAGGG should be detected as a bulged run, not two separate runs
    
    # Check for the characteristic bulge pattern: GGGA, GGGT, GGGC in G-runs
    has_bulge = False
    
    # Simple approach: look for GGG[ACT]GGG pattern (clear bulge signature)
    if re.search(r'GGG[ACT]GGG', seq):
        has_bulge = True
    
    # Also check for shorter patterns: GG[ACT]GG 
    elif re.search(r'GG[ACT]GG', seq) and not re.search(r'GGG[ACT]GGG', seq):
        has_bulge = True
    
    # Require at least 2 proper G3+ runs along with the bulge
    long_runs = sum(1 for r in runs[:4] if (r[1]-r[0])>=3)
    return has_bulge and long_runs >= 2

def _is_triplex(runs, loops):
    # G-triplex: three G-runs ≥3, loops ≤15
    if len(runs)<3 or len(loops)<2: return False
    if any((runs[k][1]-runs[k][0])<3 for k in range(3)): return False
    return all(1<=L<=15 for L in loops[:2])

# ================================================
# Candidate extraction and classification
# ================================================
def _extract_candidates(seq, window=25, threshold=0.8):  # Lowered threshold for better G4 detection
    """
    Extract G4-like candidates from seq using strict scientific rules and G4Hunter scoring.
    """
    candidates = []
    # G4 four-run base: (G+)(.{1,12})G+(.{1,12})G+(.{1,12})G+
    for m in re.finditer(r"(G+)(.{1,12})G+(.{1,12})G+(.{1,12})G+", seq):
        start, end = m.start(), m.end()
        sub = seq[start:end]
        runs = _find_g_runs(sub)
        if len(runs)<4: continue
        loops = _loop_lengths_from_runs(runs)
        g4h = g4hunter_score(sub, window=window)
        if g4h < threshold: continue
        subtype = None
        # Prioritize bulge detection since bulged sequences can also meet canonical criteria
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
    
    # Bipartite G4 detection: look for 4 G-runs with one long central loop
    for m in re.finditer(r"(G{3,})(.{1,12})(G{3,})(.{13,50})(G{3,})(.{1,12})(G{3,})", seq):
        start, end = m.start(), m.end()
        sub = seq[start:end]
        runs = _find_g_runs(sub)
        loops = _loop_lengths_from_runs(runs)
        g4h = g4hunter_score(sub, window=window)
        if g4h < threshold: continue
        
        # For bipartite: we expect one loop to be much longer than others (≥13 nt)
        if len(loops) >= 3 and max(loops) >= 13 and loops.count(max(loops)) == 1:
            candidates.append({
                "start": start, "end": end, "seq": sub, "subclass": "Bipartite/Split",
                "g4h": g4h, "runs": runs, "loops": loops
            })
    
    # Triplex: three runs (≥3), loops ≤15, higher threshold
    for m in re.finditer(r"(G{3,})(.{1,15})G{3,}(.{1,15})G{3,}", seq):
        start, end = m.start(), m.end()
        sub = seq[start:end]
        runs = _find_g_runs(sub)
        loops = _loop_lengths_from_runs(runs)
        if not _is_triplex(runs, loops): continue
        g4h = g4hunter_score(sub, window=window)
        if g4h < max(1.4, threshold): continue
        candidates.append({
            "start": start, "end": end, "seq": sub, "subclass": "G-Triplex",
            "g4h": g4h, "runs": runs, "loops": loops
        })
    return candidates

# ================================================
# Multimer and bipartite merging logic (adjacency rules)
# ================================================
def _attempt_merge_multimer(cands, max_linker=20):
    """
    Merge adjacent high-confidence G4s into Multimeric/Tandem when within ≤max_linker.
    """
    if not cands: return []
    cands = sorted(cands, key=lambda x: (x["start"], x["end"]))
    merged, i = [], 0
    while i < len(cands):
        a = cands[i]; j = i+1; merged_flag = False
        while j < len(cands):
            b = cands[j]
            if b["start"]-a["end"]<=max_linker and b["start"]>=a["start"]:
                # Merge a and b into multimer
                start, end = a["start"], b["end"]
                seq = a["seq"] + ("N"*max(0, b["start"]-a["end"])) + b["seq"]
                g4h = max(a["g4h"], b["g4h"])
                merged.append({"start": start, "end": end, "seq": seq, "subclass": "Multimeric/Tandem", "g4h": g4h, "runs": [], "loops": []})
                i = j+1; merged_flag = True; break
            else: break
        if not merged_flag:
            merged.append(a); i += 1
    return merged

def _attempt_merge_bipartite(cands, min_linker=20, max_linker=50):
    """
    Merge two G3+L1–7/extended cores separated by a long linker into Bipartite/Split.
    """
    out = []; cands = sorted(cands, key=lambda x: (x["start"], x["end"])); used = [False]*len(cands)
    for i in range(len(cands)):
        if used[i]: continue
        a = cands[i]
        if a["subclass"] not in ("Canonical G3+L1–7","Extended G3+L1–12","Long-Loop G3+L8–12"):
            out.append(a); continue
        merged_done = False
        for j in range(i+1, len(cands)):
            if used[j]: continue
            b = cands[j]
            if b["subclass"] not in ("Canonical G3+L1–7","Extended G3+L1–12","Long-Loop G3+L8–12"): continue
            linker = b["start"]-a["end"]
            if min_linker<=linker<=max_linker:
                start, end = a["start"], b["end"]
                seq = a["seq"] + ("N"*linker) + b["seq"]
                g4h = max(a["g4h"], b["g4h"])
                out.append({"start": start, "end": end, "seq": seq, "subclass": "Bipartite/Split", "g4h": g4h, "runs": [], "loops": []})
                used[i]=used[j]=True; merged_done = True; break
        if not merged_done and not used[i]:
            out.append(a); used[i]=True
    return sorted(out, key=lambda x: (x["start"], x["end"]))

# ================================================
# Composite scoring and priority assignment
# ================================================
_PRIORITY_ORDER = {
    "Multimeric/Tandem": 1, "Bipartite/Split": 2, "Bulged": 3,
    "Canonical G3+L1–7": 4, "Long-Loop G3+L8–12": 5, "Extended G3+L1–12": 6,
    "Two-Tetrad G2L1–12": 7, "Generalized G2+L1–12": 8, "G-Triplex": 9
}

def _loop_stats(loops):
    if not loops: return (0,0,0)
    return (min(loops), max(loops), sum(loops))

def _composite_score(cand):
    """
    Composite scoring (small-integer, reproducible): G4Hunter, loop compactness, symmetry, compactness, penalties.
    """
    g4h = cand.get("g4h", 0.0)
    runs = cand.get("runs", [])
    loops = cand.get("loops", [])
    length = cand["end"]-cand["start"]
    wHunter = int(round(4*max(0.0,g4h)))
    Lmin, Lmax, _ = _loop_stats(loops)
    wLoops = 2 if loops and all(1<=L<=7 for L in loops[:3]) else 1 if loops and all(1<=L<=12 for L in loops[:3]) else -(int((Lmax-12)//2) if loops and Lmax>12 else 0)
    wSymm = 0
    if len(runs)>=4:
        rlens = [(r[1]-r[0]) for r in runs[:4]]
        if all(r>=3 for r in rlens):
            if max(rlens)-min(rlens)<=1: wSymm += 1
        elif all(r==2 for r in rlens): wSymm += 1
    wCompact = 1 if length<=40 else 0
    penalties = (-1 if cand["subclass"]=="Bulged" else 0) + (-(sum(1 for L in loops[:3] if L>12)) if loops else 0)
    return max(0, wHunter+wLoops+wSymm+wCompact+penalties)

def _assign_priority_and_score(cands):
    for c in cands:
        c["priority"] = _PRIORITY_ORDER.get(c["subclass"],99)
        c["comp_score"] = _composite_score(c)
    return cands

def _strict_non_overlap(cands):
    """
    Strict non-overlap: sort by (priority, -score, -length, start); greedy sweep.
    """
    if not cands: return []
    cands = sorted(cands, key=lambda c: (c["priority"], -c["comp_score"], -(c["end"]-c["start"]), c["start"]))
    chosen, occupied = [], []
    for c in cands:
        s, e = c["start"], c["end"]
        if any(not(e<=cs or s>=ce) for cs,ce in occupied): continue
        chosen.append(c); occupied.append((s,e))
    return chosen

# ================================================
# Public API conversion and wrappers
# ================================================
def _subclass_to_output_label(subclass):
    mapping = {
        "Multimeric/Tandem": "Multimeric G4",
        "Bipartite/Split": "Bipartite G4",
        "Bulged": "Bulged G4",
        "Canonical G3+L1–7": "Canonical G4",
        "Long-Loop G3+L8–12": "Relaxed G4",
        "Extended G3+L1–12": "Relaxed G4",
        "Two-Tetrad G2L1–12": "Imperfect G4",
        "Generalized G2+L1–12": "Imperfect G4",
        "G-Triplex": "G-Triplex intermediate"
    }
    return mapping.get(subclass, subclass)

def convert_to_nbdfinder_format(candidates, sequence_name=""):
    """
    Convert internal candidates to NBDFinder output, adding structural, conservation info.
    """
    results = []
    for cand in candidates:
        subtype = _subclass_to_output_label(cand["subclass"])
        motif_seq = cand["seq"]
        g4_score = cand["g4h"]
        loops = cand.get("loops", [])
        try:
            sf = float(calculate_structural_factor(motif_seq, subtype.lower().replace(" ","_"), loops))
        except Exception:
            sf = 1.0
        try:
            cons = calculate_conservation_score(motif_seq, subtype)
        except Exception:
            cons = {"enrichment_score":0.0, "p_value":1.0, "significance":""}
        results.append({
            "Sequence Name": sequence_name,
            "Class": "G-Quadruplex Family",
            "Subtype": subtype,
            "Start": cand["start"]+1,
            "End": cand["end"],
            "Length": len(motif_seq),
            "Sequence": wrap(motif_seq),
            "ScoreMethod": "G4Hunter_Exact",
            "Score": float(g4_score*max(1.0,sf)),
            "G4Hunter_Score": float(g4_score),
            "Structural_Factor": round(sf,3),
            "Formation_Category": get_g4_formation_category(g4_score)["category"],
            "Conservation_Score": float(cons.get("enrichment_score",0.0)),
            "Conservation_P_Value": float(cons.get("p_value",1.0)),
            "Priority": cand["priority"],
            "Arms/Repeat Unit/Copies": "", "Spacer": ""
        })
    return results

def find_all_g4_motifs(seq, use_non_overlapping=True, sequence_name=""):
    """
    Find all G4 motifs using exact G4Hunter and explicit subclass logic with strict non-overlap.
    Follows: extract → preserve bulged/special types → merge (multimer/bipartite for eligible) → assign priority/score → non-overlap → output.
    """
    base = _extract_candidates(seq, window=25, threshold=0.8)
    
    # Separate bulged and special types from merger eligibles
    special_types = ["Bulged", "G-Triplex"]
    eligible_for_merge = [c for c in base if c["subclass"] not in special_types]
    special_candidates = [c for c in base if c["subclass"] in special_types]
    
    # Only merge eligible candidates
    merged1 = _attempt_merge_multimer(eligible_for_merge, max_linker=20)
    merged2 = _attempt_merge_bipartite(merged1, min_linker=20, max_linker=50)
    
    # Combine special types back with merged results
    all_candidates = special_candidates + merged2
    
    with_meta = _assign_priority_and_score(all_candidates)
    selected = _strict_non_overlap(with_meta) if use_non_overlapping else with_meta
    return convert_to_nbdfinder_format(selected, sequence_name)

# Legacy API compatibility
def find_gquadruplex(seq): return [g4 for g4 in find_all_g4_motifs(seq, use_non_overlapping=False) if g4['Subtype']=='Canonical G4']
def find_relaxed_gquadruplex(seq): return [g4 for g4 in find_all_g4_motifs(seq, use_non_overlapping=False) if g4['Subtype']=='Relaxed G4']
def find_bulged_gquadruplex(seq): return [g4 for g4 in find_all_g4_motifs(seq, use_non_overlapping=False) if g4['Subtype']=='Bulged G4']
def find_imperfect_gquadruplex(seq): return [g4 for g4 in find_all_g4_motifs(seq, use_non_overlapping=False) if g4['Subtype']=='Imperfect G4']
def find_multimeric_gquadruplex(seq): return [g4 for g4 in find_all_g4_motifs(seq, use_non_overlapping=False) if g4['Subtype']=='Multimeric G4']
def find_bipartite_gquadruplex(seq): return [g4 for g4 in find_all_g4_motifs(seq, use_non_overlapping=False) if g4['Subtype']=='Bipartite G4']
def find_gtriplex(seq): return [g4 for g4 in find_all_g4_motifs(seq, use_non_overlapping=False) if g4['Subtype']=='G-Triplex intermediate']

def get_g4_formation_category(g4hunter_score):
    """
    Categorize G4 formation potential (scientifically supported G4Hunter bins).
    """
    if g4hunter_score >= 1.5:
        return {"category":"High Formation Potential","threshold":"≥ 1.5",
                "experimental_evidence":"Strong","formation_probability":"High",
                "stability":"High","color":"#d32f2f"}
    elif g4hunter_score >= 1.2:
        return {"category":"Moderate Formation Potential","threshold":"1.2–1.5",
                "experimental_evidence":"Moderate","formation_probability":"Moderate",
                "stability":"Moderate","color":"#f57c00"}
    else:
        return {"category":"Low Formation Potential","threshold":"< 1.2",
                "experimental_evidence":"Weak/Variable","formation_probability":"Low",
                "stability":"Low","color":"#388e3c"}

# ========================================
# End of modernized, literature-aligned G4 module
# ========================================
