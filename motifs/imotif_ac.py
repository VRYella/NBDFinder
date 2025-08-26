"""
Category 7: i-motif and AC-motif Detection Module (2024, rigorous)
==================================================================

Implements exact, strand-aware i-motif and AC-motif detection:
- i-motif: strict C3+ runs, windowed iM-Hunter C-centric scoring, loop class reporting, experimental thresholds, strict non-overlap
- AC-motif: Hur et al. NAR 2021 A/C repeat criteria, composition/length/run requirements, AC-score, strict non-overlap (i-motif priority)

References:
- Zeraati et al. (2018) Nat Chem; Kocsis et al. (2021) Int J Mol Sci; Varizhuk et al. (2019) Biochimie
- Hur et al. (2021) NAR; iM-Seeker (2024) for benchmarking

Author: Dr. Venkata Rajesh Yella (updates by Copilot, 2024)
"""

import re
from .shared_utils import wrap, calculate_conservation_score, overlapping_finditer

# ---- iM-Hunter exact: per-base, C-centric, windowed ----
def _im_base_scores(seq):
    """Per-base iM-Hunter: C: +k, G: –k, A/T: 0; k = run_length capped at 4 [C-rich folding]."""
    n = len(seq); s = [0]*n
    c_run = 0
    for i, ch in enumerate(seq):
        c_run = c_run + 1 if ch == 'C' else 0
        if ch == 'C': s[i] = min(4, c_run)
    g_run = 0
    for i, ch in enumerate(seq):
        g_run = g_run + 1 if ch == 'G' else 0
        if ch == 'G': s[i] = -min(4, g_run)
    return s

def _sliding_mean(arr, window=25):
    """Sliding mean over integer array, window centered at each base."""
    n = len(arr); W = max(1, int(window)); half = W//2
    pref = [0]
    for v in arr: pref.append(pref[-1]+v)
    out = [0.0]*n
    for i in range(n):
        a = max(0,i-half); b = min(n,i+half+1)
        out[i] = (pref[b]-pref[a])/(b-a) if b>a else 0.0
    return out

def imotif_score(seq, window=25):
    """Exact iM-Hunter max windowed score (canonical ≥0.45, relaxed ≥0.3)."""
    if not seq: return 0.0
    base = _im_base_scores(seq)
    win = _sliding_mean(base, window)
    return float(max(win) if win else 0.0)

# ---- Helper: find C-runs, loop lengths ----
def _find_c_runs(seq):
    """Returns list of (start,end) for all C runs in seq."""
    return [m.span() for m in re.finditer(r"C+", seq)]

def _loop_lengths(spans):
    """Given C-run spans, returns list of loop lengths (between runs)."""
    return [spans[i+1][0]-spans[i][1] for i in range(len(spans)-1)]

def _is_canonical_c3_loops_1_7(c_spans, loops):
    """True if ≥4 C3+ runs, first 3 loops all 1–7."""
    if len(c_spans)<4 or len(loops)<3: return False
    if any((e-s)<3 for (s,e) in c_spans[:4]): return False
    return all(1<=L<=7 for L in loops[:3])

def _is_relaxed_c3_loops_1_12(c_spans, loops):
    """True if ≥4 C3+ runs, first 3 loops all 1–12, not canonical."""
    if len(c_spans)<4 or len(loops)<3: return False
    if any((e-s)<3 for (s,e) in c_spans[:4]): return False
    if not all(1<=L<=12 for L in loops[:3]): return False
    return not _is_canonical_c3_loops_1_7(c_spans, loops)

# ---- i-motif detection: canonical & relaxed (non-overlapping) ----
def find_imotif(seq):
    """
    Find strict i-motif structures using experimental C-run/loop/score criteria.
    Returns: list of dicts, canonical/relaxed, non-overlapping, fielded.
    """
    from .shared_utils import calculate_structural_factor
    results, used = [], []
    # Seed: four C-runs, up to 12-nt loops (canonical+relaxed window)
    seed_pat = r"(?=(C+).{1,12}C+.{1,12}C+.{1,12}C+)"
    for m in overlapping_finditer(seed_pat, seq):
        start = m.start()
        window_seq = seq[start:start+len(m.group(0))]
        c_spans = _find_c_runs(window_seq)
        if len(c_spans)<4: continue
        s0,e0 = c_spans[0]
        s3,e3 = c_spans[3]
        motif_seq = window_seq[s0:e3]
        abs_start = start + s0; abs_end = start + e3
        c_runs = _find_c_runs(motif_seq)
        loops = _loop_lengths(c_runs)
        if len(c_runs)<4 or len(loops)<3: continue
        im_score = imotif_score(motif_seq, window=25)
        # Annotate loop profile class for interpretability
        subtype = None
        if _is_canonical_c3_loops_1_7(c_runs, loops) and im_score >= 0.45:
            subtype = "Canonical i-motif"
        elif _is_relaxed_c3_loops_1_12(c_runs, loops) and im_score >= 0.30:
            subtype = "Relaxed i-motif"
        if not subtype: continue
        structural_factor = calculate_structural_factor(motif_seq, subtype, loops)
        conservation_result = calculate_conservation_score(motif_seq, "i-Motif")
        c_run_lengths = [e-s for (s,e) in c_runs]
        c_fraction = motif_seq.count('C')/max(1,len(motif_seq))
        rec = {
            "Sequence Name": "", "Class": "i-motif family", "Subtype": subtype,
            "Start": abs_start+1, "End": abs_end, "Length": abs_end-abs_start,
            "Sequence": wrap(motif_seq), "ScoreMethod": "iM_Hunter_Exact",
            "Score": float(im_score*max(1.0,structural_factor)), "iMotif_Mean": float(im_score),
            "Structural_Factor": round(structural_factor,3), "C_Run_Count": len(c_run_lengths),
            "C_Run_Sum": sum(c_run_lengths), "C_Fraction": round(c_fraction,3),
            "Conservation_Score": float(conservation_result.get("enrichment_score",0.0)),
            "Conservation_P_Value": float(conservation_result.get("p_value",1.0)),
            "Conservation_Significance": conservation_result.get("significance",""),
            "Loop_Lengths": loops, "Arms/Repeat Unit/Copies": "", "Spacer": ""
        }
        results.append(rec)
        used.append((rec["Start"],rec["End"],1 if subtype=="Canonical i-motif" else 2))
    # Strict non-overlap (canonical > relaxed > score > length)
    results = sorted(results, key=lambda r: ((1 if r["Subtype"]=="Canonical i-motif" else 2), -r["Score"], -r["Length"], r["Start"]))
    final_im, occ = [], []
    for r in results:
        s,e = r["Start"], r["End"]
        if any(not(e<=cs or s>=ce) for cs,ce in occ): continue
        final_im.append(r); occ.append((s,e))
    return final_im

# ---- AC-motif detection: Hur et al. NAR 2021 ----
def ac_motif_score(seq):
    """
    AC-motif score: A/C fraction, A3+/C3+ tracts, alternation density. Conservative to avoid outranking i-motif.
    """
    if len(seq)==0: return 0.0
    a3_runs = len(re.findall(r"A{3,}",seq))
    c3_runs = len(re.findall(r"C{3,}",seq))
    ac_fraction = (seq.count('A')+seq.count('C'))/len(seq)
    alt_blocks = len(re.findall(r"(?:AC|CA){3,}",seq))
    return (len(seq)*ac_fraction*0.3)+1.5*(a3_runs+c3_runs)+1.0*alt_blocks

def find_ac_motifs(seq):
    """
    AC-motif (Hur et al. NAR 2021): A/C repeat segments, C3+ anchor, AC-content ≥0.6, length ≥15, score ≥10.
    Returns: non-overlapping set, i-motif priority if overlaps.
    """
    from .shared_utils import calculate_structural_factor
    results, used = [], []
    patterns = [
        r"C{3,}[AC]{6,}C{3,}",           # C3+ anchor - AC-rich - C3+
        r"(?:AC|CA){6,}",                # long alternating blocks
        r"A{3,}[AC]{6,}C{3,}",           # A3+ anchor - AC-rich - C3+
    ]
    for pat in patterns:
        for m in re.finditer(pat, seq):
            motif_seq = m.group(0); L = len(motif_seq)
            if L<15: continue
            a_frac = motif_seq.count('A')/L; c_frac = motif_seq.count('C')/L; ac_frac = a_frac+c_frac
            if ac_frac<0.6: continue
            if not re.search(r"C{3,}", motif_seq): continue
            score = ac_motif_score(motif_seq)
            if score<10.0: continue
            s,e = m.start()+1, m.end()
            if any(not(e<=cs or s>=ce) for cs,ce in used): continue
            conservation_result = calculate_conservation_score(motif_seq, "AC-Motif")
            structural_factor = calculate_structural_factor(motif_seq, "AC-motif", [])
            a3_runs = len(re.findall(r"A{3,}",motif_seq))
            c3_runs = len(re.findall(r"C{3,}",motif_seq))
            results.append({
                "Sequence Name": "", "Class": "i-motif family", "Subtype": "AC-motif",
                "Start": s, "End": e, "Length": L, "Sequence": wrap(motif_seq),
                "ScoreMethod": "AC_UnifiedFramework_raw", "Score": float(score*max(1.0,structural_factor)),
                "AC_Fraction": round(ac_frac,3), "A_Fraction": round(a_frac,3), "C_Fraction": round(c_frac,3),
                "A3_Runs": a3_runs, "C3_Runs": c3_runs, "Structural_Factor": round(structural_factor,3),
                "Pattern_Type": "AC-core", "Conservation_Score": float(conservation_result.get("enrichment_score",0.0)),
                "Conservation_P_Value": float(conservation_result.get("p_value",1.0)),
                "Conservation_Significance": conservation_result.get("significance",""),
                "Arms/Repeat Unit/Copies": f"A3={a3_runs};C3={c3_runs}", "Spacer": ""
            })
            used.append((s,e))
    return results

# ---- Final: strict non-overlap across i-motif/AC-motif ----
def find_imotif_and_ac(seq):
    """
    Run i-motif and AC-motif detectors, yield strict non-overlapping calls by priority:
    1. Canonical i-motif; 2. Relaxed i-motif; 3. AC-motif (Hur et al. NAR 2021)
    """
    im = find_imotif(seq)
    ac = find_ac_motifs(seq)
    def pri(r): return 1 if r["Subtype"]=="Canonical i-motif" else 2 if r["Subtype"]=="Relaxed i-motif" else 3
    allrecs = im + ac
    allrecs = sorted(allrecs, key=lambda r: (pri(r), -float(r["Score"]), -int(r["Length"]), int(r["Start"])))
    chosen, occ = [], []
    for r in allrecs:
        s,e = r["Start"], r["End"]
        if any(not(e<=cs or s>=ce) for cs,ce in occ): continue
        chosen.append(r); occ.append((s,e))
    return chosen
