"""
Category 6: Z-DNA and eGZ Detection Module
==========================================

Implements detection algorithms for Z-DNA and eGZ (extruded-G) structures using a
transition-weighted dinucleotide scoring array and a Kadane-style maximum subarray
selection. Parameter defaults and logic follow literature: GC/CG ≈ 7; GT/TG/AC/CA ≈ 1.25;
AT/TA ≈ 0.5 (with heavy penalties for long AT/TA runs); region score threshold ≈ 50.
Non-overlapping calls are enforced, with Z-DNA > eGZ precedence.

Author: Dr. Venkata Rajesh Yella
Updated: 2025
"""

import re

# =============== Transition-weight scoring for Z-DNA ===============
def zdna_dinucleotide_weights(
    seq,
    gc_weight=7.0,  # strong Z-support (GC/CG)
    gt_ca_weight=1.25,  # moderate Z-support (GT/TG/AC/CA)
    at_weight=0.5,  # weak Z-support (AT/TA)
    at_consecutive_penalty=(0.5, 0.5, 0.5, 0.0, 0.0, -5.0, -100.0),  # penalize long AT/TA runs
    mismatch_penalty_start=3.0,
    mismatch_penalty_mode="linear",  # "linear" or "exponential"
    mismatch_penalty_delta=3.0
):
    """
    Build a transition score array of length len(seq)-1:
    GC/CG -> +gc_weight; GT/TG/AC/CA -> +gt_ca_weight; AT/TA -> adjusted by run penalties;
    other dinucleotides -> negative mismatch penalty escalated by streak length.
    Long AT/TA runs or mismatch streaks sharply penalized to prevent Z-calls in hairpin-prone or
    non-alternating tracts.
    """
    s = seq.upper(); n = len(s)
    if n < 2: return []
    w = [0.0] * (n - 1); at_streak = 0; mm_streak = 0
    for i in range(n - 1):
        d = s[i:i+2]
        if d in ("GC", "CG"):
            w[i] = gc_weight; at_streak = 0; mm_streak = 0
        elif d in ("GT", "TG", "AC", "CA"):
            w[i] = gt_ca_weight; at_streak = 0; mm_streak = 0
        elif d in ("AT", "TA"):
            at_streak += 1
            idx = min(at_streak, len(at_consecutive_penalty)) - 1
            w[i] = at_weight + at_consecutive_penalty[idx]
            mm_streak = 0
        else:
            mm_streak += 1
            if mismatch_penalty_mode == "linear":
                penalty = mismatch_penalty_start + mismatch_penalty_delta * (mm_streak - 1)
            else:
                penalty = mismatch_penalty_start * (2 ** (mm_streak - 1))
            w[i] = -abs(penalty); at_streak = 0
    return w

# =============== Kadane-style maximal Z windows ===============
def kadane_max_subarrays(weights, min_transitions=12, drop_threshold=None):
    """
    Identify high-scoring subarrays using a Kadane-like search.
    Returns (start_transition_idx, end_transition_idx, sum_score), non-overlapping, sorted by score.
    """
    n = len(weights)
    if n < min_transitions: return []
    candidates = []
    for start in range(0, n - min_transitions + 1):
        best_sum = float("-inf"); best_end = start; run_sum = 0.0; current_start = start
        for end in range(start, n):
            run_sum += weights[end]
            if run_sum > best_sum: best_sum = run_sum; best_end = end
            if run_sum < 0: run_sum = 0.0; current_start = end + 1
            if drop_threshold is not None and run_sum < -abs(drop_threshold): break
        if best_sum > 0 and (best_end - start + 1) >= min_transitions:
            candidates.append((start, best_end, best_sum))
    candidates.sort(key=lambda x: (-x[2], -(x[1]-x[0]+1), x[0]))
    non_overlapping = []; occupied = []
    for s, e, score in candidates:
        if any(not (e < os or s > oe) for os, oe, _ in occupied): continue
        non_overlapping.append((s, e, score)); occupied.append((s, e, score))
    return non_overlapping

# =============== Z-DNA detector (strict, thresholded) ===============
def _annotate_z_region(seq, t_start, t_end, raw_score):
    """
    Convert transition indexes [t_start..t_end] to base coordinates.
    Returns: base_start, base_end, motif_seq, length, gc_frac, gc_cg_dinucs.
    """
    base_start = t_start; base_end = t_end + 1
    motif_seq = seq[base_start:base_end+1]
    L = len(motif_seq)
    gc_frac = (motif_seq.count("G") + motif_seq.count("C")) / max(1, L)
    gc_cg = sum(1 for i in range(L-1) if motif_seq[i:i+2] in ("GC", "CG"))
    return base_start, base_end, motif_seq, L, gc_frac, gc_cg

def find_zdna(
    seq,
    threshold=35.0,  # lowered threshold for better sensitivity
    min_length_nt=13,  # minimum NT span (>=12 transitions)
    min_transitions=12,
    drop_threshold=35.0,  # early-abandon guard (lowered)
    gc_weight=7.0,
    gt_ca_weight=1.25,
    at_weight=0.5,
    at_consecutive_penalty=(0.5,0.5,0.5,0.0,0.0,-5.0,-100.0),
    mismatch_penalty_start=3.0,
    mismatch_penalty_mode="linear",
    mismatch_penalty_delta=3.0,
    conservation_fn=None,  # optional callable(seq, label) -> dict
    structural_factor_fn=None  # optional callable(seq, subtype_key, loops) -> float
):
    """
    Detect Z-DNA candidate regions using transition-weight scoring and Kadane-style selection.
    - threshold≈50 and weights (GC/CG~7; GT/TG/AC/CA~1.25; AT/TA~0.5 with penalties) separate
      Z-adopting alternating tracts from nonforming sequences.
    - Returns strict non-overlapping set, Score-desc, Length-desc, Start-asc.
    """
    s = seq.upper()
    if len(s) < min_length_nt: return []
    weights = zdna_dinucleotide_weights(
        s, gc_weight, gt_ca_weight, at_weight,
        at_consecutive_penalty, mismatch_penalty_start,
        mismatch_penalty_mode, mismatch_penalty_delta
    )
    segments = kadane_max_subarrays(weights, min_transitions, drop_threshold)
    records = []
    for t_start, t_end, raw in segments:
        b_start, b_end, motif_seq, L, gc_frac, gc_cg = _annotate_z_region(s, t_start, t_end, raw)
        if raw < threshold or L < min_length_nt: continue
        # Structural factor (optional)
        sf = 1.0
        if structural_factor_fn is not None:
            try: sf = float(structural_factor_fn(motif_seq, "z_dna", []))
            except Exception: sf = 1.0
        score = float(raw * max(1.0, sf))
        # Conservation (optional)
        cons = {"enrichment_score": 0.0, "p_value": 1.0, "significance": ""}
        if conservation_fn is not None:
            try: cons = conservation_fn(motif_seq, "Z-DNA")
            except Exception: pass
        records.append({
            "Sequence Name": "",
            "Class": "Z-DNA",
            "Subtype": "Z-DNA",
            "Start": b_start + 1,  # 1-based
            "End": b_end + 1,  # inclusive (matches previous conventions)
            "Length": L,
            "Sequence": motif_seq,
            "ScoreMethod": "Kadane_TransitionArray_v1",
            "Score": score,
            "Kadane_Score": float(raw),
            "Structural_Factor": round(sf, 3),
            "GC_Content": round(gc_frac, 3),
            "GC_CG_Dinucleotides": gc_cg,
            "Conservation_Score": float(cons.get("enrichment_score", 0.0)),
            "Conservation_P_Value": float(cons.get("p_value", 1.0)),
            "Conservation_Significance": cons.get("significance", ""),
            "Arms/Repeat Unit/Copies": "",
            "Spacer": ""
        })
    # Strict non-overlap
    records.sort(key=lambda r: (-float(r["Score"]), -int(r["Length"]), int(r["Start"])))
    chosen, occ = [], []
    for r in records:
        s1, e1 = r["Start"], r["End"]
        if any(not (e1 <= cs or s1 >= ce) for cs, ce in occ): continue
        chosen.append(r); occ.append((s1, e1))
    return chosen

# =============== eGZ (CGG expansions) detector ===============
def find_egz_motif(
    seq,
    min_repeats=3,
    merge_adjacent=True,
    max_gap=3,
    conservation_fn=None,
    structural_factor_fn=None
):
    """
    Detect CGG expansions consistent with extruded-G contexts:
    - Identify (CGG)n with n >= min_repeats.
    - Merge adjacent blocks separated by ≤max_gap (optional).
    - Score scales with repeat count and G-fraction; reported as Z-DNA class subtype eGZ.
    Priority: When resolving overlaps with Z-DNA, Z-DNA regions take precedence.
    """
    s = seq.upper()
    pat = re.compile(r"(?:CGG){%d,}" % max(1, min_repeats))
    raw_hits = [(m.start(), m.end()) for m in pat.finditer(s)]
    # Merge adjacent if needed
    merged = []
    if merge_adjacent and raw_hits:
        cs, ce = raw_hits[0]
        for a, b in raw_hits[1:]:
            if a - ce <= max_gap: ce = b
            else: merged.append((cs, ce)); cs, ce = a, b
        merged.append((cs, ce))
    else:
        merged = raw_hits
    records = []
    for a, b in merged:
        motif_seq = s[a:b]
        L = len(motif_seq)
        copies = L // 3
        g_frac = motif_seq.count("G") / max(1, L)
        base_score = copies * (1.0 + 2.0 * g_frac)
        sf = 1.0
        if structural_factor_fn is not None:
            try: sf = float(structural_factor_fn(motif_seq, "egz", []))
            except Exception: sf = 1.0
        cons = {"enrichment_score": 0.0, "p_value": 1.0, "significance": ""}
        if conservation_fn is not None:
            try: cons = conservation_fn(motif_seq, "eGZ")
            except Exception: pass
        records.append({
            "Sequence Name": "",
            "Family": "Double-stranded",
            "Class": "Z-DNA",
            "Subtype": "eGZ (Extruded-G) DNA",
            "Start": a + 1,
            "End": b,
            "Length": L,
            "Sequence": motif_seq,
            "ScoreMethod": "eGZ_CGG_Expansion_v1",
            "Score": float(base_score * max(1.0, sf)),
            "CGG_Repeats": copies,
            "G_Fraction": round(g_frac, 3),
            "Structural_Factor": round(sf, 3),
            "Conservation_Score": float(cons.get("enrichment_score", 0.0)),
            "Conservation_P_Value": float(cons.get("p_value", 1.0)),
            "Conservation_Significance": cons.get("significance", ""),
            "Arms/Repeat Unit/Copies": f"Unit=CGG;Copies={copies}",
            "Spacer": ""
        })
    # Strict non-overlap within eGZ set
    records.sort(key=lambda r: (-float(r["Score"]), -int(r["Length"]), int(r["Start"])))
    chosen, occ = [], []
    for r in records:
        s1, e1 = r["Start"], r["End"]
        if any(not (e1 <= cs or s1 >= ce) for cs, ce in occ): continue
        chosen.append(r); occ.append((s1, e1))
    return chosen

# =============== Unified resolver (Z-DNA > eGZ) ===============
def find_zdna_and_egz(
    seq,
    zdna_kwargs=None,
    egz_kwargs=None,
    enforce_nonoverlap=True
):
    """
    Run Z-DNA and eGZ and resolve a final, non-overlapping set with priority:
    eGZ (CGG expansions) > Z-DNA (alternating Pu/Py) for CGG-rich sequences.
    Z-DNA > eGZ for other alternating sequences.
    """
    zdna_kwargs = zdna_kwargs or {}
    egz_kwargs = egz_kwargs or {}
    
    # First, identify CGG-rich regions
    cgg_density = seq.count('CGG') / max(1, len(seq) // 3)
    
    e_calls = find_egz_motif(seq, **egz_kwargs)
    z_calls = find_zdna(seq, **zdna_kwargs)
    
    if not enforce_nonoverlap: return z_calls + e_calls
    
    # Merge and enforce priority based on sequence characteristics
    allrecs = []
    
    # For sequences with high CGG density, prioritize eGZ
    if cgg_density > 0.5:  # More than 50% of possible positions are CGG
        for r in e_calls: r["_priority"] = 1; allrecs.append(r)
        for r in z_calls: r["_priority"] = 2; allrecs.append(r)
    else:
        # For other sequences, prioritize Z-DNA
        for r in z_calls: r["_priority"] = 1; allrecs.append(r)
        for r in e_calls: r["_priority"] = 2; allrecs.append(r)
    
    allrecs.sort(key=lambda r: (r["_priority"], -float(r["Score"]), -int(r["Length"]), int(r["Start"])))
    chosen, occ = [], []
    for r in allrecs:
        s1, e1 = r["Start"], r["End"]
        if any(not (e1 <= cs or s1 >= ce) for cs, ce in occ): continue
        r.pop("_priority", None)
        chosen.append(r); occ.append((s1, e1))
    return chosen

# ---- Notes on defaults and tuning ----
# Keep default region threshold near 50 with GC/CG=7, GT/TG/AC/CA=1.25, AT/TA=0.5, and steep AT-run
# penalties; these recapitulate experimental minima for Z tracts. Raise min_transitions or threshold for
# more stringency. For eGZ, raise min_repeats or reduce max_gap. Z-DNA calls always take precedence.
# Conservation and structural-factors are optional hooks.
