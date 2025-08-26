"""
Category 4: Triplex DNA Detection Module
=======================================

This module implements detection algorithms for triplex DNA (H-DNA) and sticky DNA structures.

Scientific Basis:
- Intramolecular triplex (H-DNA) forms at homopurine•homopyrimidine mirror repeats via Hoogsteen or reverse Hoogsteen pairing in acidic or crowded conditions; stability rises with longer homopurine/pyrimidine arms and short loops. 
- Sticky DNA is associated with GAA/TTC triplet-repeat expansions (Friedreich’s ataxia), forming inter- or intramolecular complexes with characteristic length dependence.
- Here, detection is sequence-only: enriched purine/pyrimidine mirror architecture with bounded loop size (H-DNA) and long GAA/TTC runs (sticky DNA); scoring is deterministic and mismatch-free for reproducibility.

Author: Dr. Venkata Rajesh Yella
Updated: 2025
"""

import re
from .shared_utils import wrap, calculate_conservation_score

# --------- Utilities: sequence composition and mirrors ---------

def purine_fraction(seq):
    """Fraction of purines (A,G) in sequence (A/G only)."""
    if not seq: return 0.0
    return (seq.count('A') + seq.count('G')) / len(seq)

def pyrimidine_fraction(seq):
    """Fraction of pyrimidines (C,T) in sequence (C/T only)."""
    if not seq: return 0.0
    return (seq.count('C') + seq.count('T')) / len(seq)

def _is_homopurine(seq, thresh=0.9):
    """Return True if sequence is ≥ thresh fraction A/G."""
    return purine_fraction(seq) >= thresh

def _is_homopyrimidine(seq, thresh=0.9):
    """Return True if sequence is ≥ thresh fraction C/T."""
    return pyrimidine_fraction(seq) >= thresh

def _mirror_repeat_slices(seq, arm_len, spacer):
    """
    Yields (start, end, armL, loop, armR) for mirror repeats:
    - Mirror: X (loop) mirror(X), where mirror(X) is the letterwise mirror (not reverse complement).
    - Only mirror (not inverted or palindromic).
    """
    n = len(seq)
    if arm_len <= 0 or spacer < 0: 
        return
    total = 2*arm_len + spacer
    for i in range(0, n - total + 1):
        left = seq[i:i+arm_len]
        loop = seq[i+arm_len:i+arm_len+spacer] if spacer > 0 else ""
        right = seq[i+arm_len+spacer:i+arm_len+spacer+arm_len]
        # Require right to be a perfect mirror of left (mirror, not palindrome)
        if all(left[j] == right[arm_len-1-j] for j in range(arm_len)):
            yield (i, i+total, left, loop, right)

def _triplex_score(arm_len, homogeneity, spacer_len, total_len):
    """
    Deterministic H-DNA score: arm length × homogeneity, penalize loop, small bonus for long mirrors.
    """
    length_step = 1 if total_len >= 50 else 0
    return arm_len * (1.0 + 1.5*homogeneity) - 0.8*spacer_len + length_step

def _non_overlap(records):
    """
    Strict non-overlap: sort by Score desc, Length desc, Start asc.
    Returns the maximal set of non-overlapping motifs.
    """
    recs = sorted(records, key=lambda r: (-float(r["Score"]), -int(r["Length"]), int(r["Start"])))
    chosen, occ = [], []
    for r in recs:
        s, e = r["Start"], r["End"]
        if any(not (e <= cs or s >= ce) for cs, ce in occ):  # if overlap
            continue
        chosen.append(r)
        occ.append((s, e))
    return chosen

# --------- H-DNA (intramolecular triplex) detection ---------

def find_hdna(seq):
    """
    Detect H-DNA (triplex) at homopurine•homopyrimidine mirror repeats.
    - Arms: 10–100 nt, loop: 0–8 nt (≤12 for permissive).
    - Require arm homogeneity: one arm ≥0.9 purine, mirrored arm ≥0.9 pyrimidine, or vice versa.
    - Mirror defined as left == mirror(right).
    - Score: arm length × homogeneity, penalize loop.
    - Strict non-overlap enforced.
    - Output fields compatible with previous API, ScoreMethod updated.
    """
    results = []
    n = len(seq)
    min_arm, max_arm = 10, 100
    max_spacer = 8
    min_total = 20
    min_score_threshold = 25.0

    for arm_len in range(min_arm, min(max_arm, n//2) + 1):
        for spacer in range(0, max_spacer + 1):
            for s, e, left, loop, right in _mirror_repeat_slices(seq, arm_len, spacer):
                total_len = e - s
                if total_len < min_total:
                    continue

                # Homopurine/homopyrimidine arms: either arm is purine-rich and its mirror is pyrimidine-rich
                left_pur, left_pyr = _is_homopurine(left), _is_homopyrimidine(left)
                right_pur, right_pyr = _is_homopurine(right), _is_homopyrimidine(right)
                valid = (left_pur and right_pyr) or (left_pyr and right_pur)
                if not valid:
                    continue

                arms_seq = left + right
                pur_frac = purine_fraction(arms_seq)
                pyr_frac = pyrimidine_fraction(arms_seq)
                homogeneity = max(pur_frac, pyr_frac)
                score = _triplex_score(arm_len, homogeneity, spacer, total_len)
                if score < min_score_threshold:
                    continue

                full_seq = seq[s:e]
                cons = calculate_conservation_score(full_seq, "Triplex")
                results.append({
                    "Sequence Name": "",
                    "Class": "Triplex",
                    "Subtype": "Triplex",
                    "Start": s + 1,
                    "End": e,
                    "Length": total_len,
                    "Spacer": str(spacer),
                    "Sequence": wrap(full_seq),
                    "ScoreMethod": "Triplex_Mirror_v1",
                    "Score": float(score),
                    "PurineFrac": round(pur_frac, 3),
                    "PyrimidineFrac": round(pyr_frac, 3),
                    "Homogeneity": round(homogeneity, 3),
                    "Conservation_Score": float(cons.get("enrichment_score", 0.0)),
                    "Conservation_P_Value": float(cons.get("p_value", 1.0)),
                    "Conservation_Significance": cons.get("significance", ""),
                    "Arms/Repeat Unit/Copies": f"Arms={arm_len}",
                    "Spacer": str(spacer)
                })
    # Strict non-overlap within H-DNA set
    return _non_overlap(results)

# --------- Sticky DNA (GAA/TTC) detection ---------

def find_sticky_dna(seq):
    """
    Detect sticky DNA (GAA/TTC triplet-repeat expansions).
    - Sensitivity: ≥6 copies (18 nt).
    - Pathogenic range annotated at ≥59 copies.
    - Score: proportional to repeat count, mild AT bonus.
    - Strict non-overlap enforced within sticky DNA calls.
    - Output fields compatible with previous API.
    """
    motifs = []
    s = seq.upper()
    pat = re.compile(r"(?:GAA){6,}|(?:TTC){6,}")
    for m in pat.finditer(s):
        motif_seq = m.group(0)
        L = len(motif_seq)
        copies = L // 3
        at_frac = (motif_seq.count('A') + motif_seq.count('T')) / max(1, L)
        score = copies * (1.0 + 0.3*at_frac)
        pathogenic = copies >= 59
        cons = calculate_conservation_score(motif_seq, "Sticky DNA")
        motifs.append({
            "Sequence Name": "",
            "Class": "Triplex",
            "Subtype": "sticky DNA",
            "Start": m.start() + 1,
            "End": m.end(),
            "Length": L,
            "RepeatCount": copies,
            "Sequence": wrap(motif_seq),
            "ScoreMethod": "GAA_TTC_v1",
            "Score": float(score),
            "AT_Fraction": round(at_frac, 3),
            "Pathogenic": pathogenic,
            "Conservation_Score": float(cons.get("enrichment_score", 0.0)),
            "Conservation_P_Value": float(cons.get("p_value", 1.0)),
            "Conservation_Significance": cons.get("significance", ""),
            "Arms/Repeat Unit/Copies": f"Unit={'GAA' if motif_seq.startswith('GAA') else 'TTC'};Copies={copies}",
            "Spacer": ""
        })
    return _non_overlap(motifs)

# --------- Unified resolver: Triplex > sticky DNA, strict non-overlap ---------

def find_triplex_and_sticky(seq):
    """
    Run both detectors and return a strict non-overlapping set:
    - Priority: Triplex (H-DNA mirrors) > sticky DNA (GAA/TTC repeats).
    - No double-counting: if a region qualifies as both, only Triplex is reported.
    """
    triplex = find_hdna(seq)
    sticky = find_sticky_dna(seq)

    allrecs = []
    for r in triplex:
        r["_p"] = 1
        allrecs.append(r)
    for r in sticky:
        r["_p"] = 2
        allrecs.append(r)

    allrecs.sort(key=lambda r: (r["_p"], -float(r["Score"]), -int(r["Length"]), int(r["Start"])))
    chosen, occ = [], []
    for r in allrecs:
        s, e = r["Start"], r["End"]
        if any(not (e <= cs or s >= ce) for cs, ce in occ):  # overlaps with chosen
            continue
        r.pop("_p", None)
        chosen.append(r)
        occ.append((s, e))
    return chosen

# --- End of module ---
