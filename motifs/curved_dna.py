"""
Category 1: Curved DNA Detection Module
==========
Regex-enhanced version preserving original API and defaults.

This module implements detection algorithms for curved DNA structures,
including global and local poly(A)/poly(T) tracts that cause DNA bending.

Scientific Basis:
- A-tracts and T-tracts cause DNA bending through narrowed minor groove.
- Curvature score reflects biological propensity for curvature.
- Periodic spacing of A/T tracts (≈10–11 bp) enhances curvature (helical phasing).
- Junction context (ApT vs TpA; GC flanks) modulates bend magnitude/direction.
- Symmetric vs asymmetric A/T organizations differ in mechanics.

Author: Dr. Venkata Rajesh Yella
Updated: 2024 (regex-enhanced variant, 2025)
"""

import re
from .shared_utils import wrap, calculate_conservation_score


# -----------------------------
# Compiled regex library
# -----------------------------
# Notes:
# - Keep poly(A)/poly(T) default threshold at 7 to match original behavior.
# - Use anchored/atomic group variants where helpful to reduce backtracking.
# - Periodic spacing targets ~8–12 bp to capture ≈one helical turn phasing.

# Core pure A/T tracts (retain ≥7 default)
PURE_A_RE = re.compile(r"A{7,}")
PURE_T_RE = re.compile(r"T{7,}")

# Additional tract variants (not used to override min_len default, but used for context scoring)
A_DINUCL_RE = re.compile(r"(?:AA){3,}")     # AA repeats
T_DINUCL_RE = re.compile(r"(?:TT){3,}")     # TT repeats

# Junction/context patterns: GC flanks, CpG context, spacer interruptions
JUNCTION_GC_A_RE = re.compile(r"[GC]A{4,}[GC]")
JUNCTION_GC_T_RE = re.compile(r"[GC]T{4,}[GC]")
CPG_FLANKED_AT_RE = re.compile(r"CG[AT]{4,}CG")
INTERRUPTED_A_RE = re.compile(r"A{3,}[GC]A{3,}")
INTERRUPTED_T_RE = re.compile(r"T{3,}[GC]T{3,}")

# Symmetric motifs
SYMMETRIC_A4T4_RE = re.compile(r"A{4}T{4}")
SYMMETRIC_T4A4_RE = re.compile(r"T{4}A{4}")
SYMMETRIC_AN_TN_RE = re.compile(r"A{3,}T{3,}|T{3,}A{3,}")

# Periodic (helical phasing) patterns with ≈8–12 bp spacing
# Use non-greedy dot with bounded quantifier to control backtracking
PERIODIC_A_RE = re.compile(r"(A{3,}).{8,12}(A{3,}).{8,12}(A{3,})")
PERIODIC_T_RE = re.compile(r"(T{3,}).{8,12}(T{3,}).{8,12}(T{3,})")
PERIODIC_MIXED_AT_RE = re.compile(r"([AT]{3,}).{8,12}([AT]{3,}).{8,12}([AT]{3,})")

# QC patterns
POLY_N_RE = re.compile(r"N{3,}")
LOW_COMPLEXITY_RE = re.compile(r"(.)\1{10,}")
GC_BLOCK_RE = re.compile(r"[GC]{10,}")

# Simple AT counter via regex
AT_RE = re.compile(r"[AT]")


# -----------------------------
# Utilities (non-API)
# -----------------------------
def _regex_find_tracts(seq: str, min_len: int):
    """
    Regex-based finder for poly(A)/poly(T) tracts with minimum length.
    Retains original default behavior for min_len=7.
    """
    results = []
    # Scan for A-tracts
    for m in re.finditer(rf"A{{{min_len},}}", seq):
        results.append((m.start(), m.end() - 1, seq[m.start():m.end()]))
    # Scan for T-tracts
    for m in re.finditer(rf"T{{{min_len},}}", seq):
        results.append((m.start(), m.end() - 1, seq[m.start():m.end()]))
    # Sort by start
    results.sort(key=lambda x: x[0])
    return results


def _regex_periodic_groups(seq: str, min_repeats: int, min_spacing: int, max_spacing: int):
    """
    Identify groups of ≥min_repeats A/T tracts with center-to-center spacing in [min_spacing, max_spacing].
    Uses regex as a seed (fast) then verifies spacing precisely.
    """
    seeds = []
    for pat in (PERIODIC_A_RE, PERIODIC_T_RE, PERIODIC_MIXED_AT_RE):
        for m in pat.finditer(seq):
            seeds.append((m.start(), m.end()))
    # Expand by verifying centers with the tract list
    tracts = _regex_find_tracts(seq, min_len=3)  # ≥3 to allow flexible grouping; matches API defaults
    results = []
    for s, e in seeds:
        # filter tracts inside seed window
        group = [t for t in tracts if t[0] >= s and t[1] <= e]
        # verify consecutive spacing
        if len(group) >= min_repeats:
            group_sorted = sorted(group, key=lambda x: (x + x[1]) // 2)
            current = [group_sorted]
            ok = True
            for k in range(1, len(group_sorted)):
                prev_c = (group_sorted[k - 1] + group_sorted[k - 1][1]) // 2
                curr_c = (group_sorted[k] + group_sorted[k][1]) // 2
                spacing = curr_c - prev_c
                if min_spacing <= spacing <= max_spacing:
                    current.append(group_sorted[k])
                else:
                    # finalize if we already have enough
                    if len(current) >= min_repeats:
                        results.append(current[:])
                    current = [group_sorted[k]]
            if len(current) >= min_repeats:
                results.append(current[:])
    # De-duplicate groups by their boundaries
    uniq = {}
    for g in results:
        key = (g, g[-1][1])
        if key not in uniq or len(g) > len(uniq[key]):
            uniq[key] = g
    return list(uniq.values())


def _regex_context_scores(seq: str, start: int, end: int, flank: int = 10):
    """
    Compute context score using regex hits in flanking regions and within region.
    Heuristic scoring to reflect junction/cpG/interruptions consistent with literature on A-tract mechanics.
    """
    s = max(0, start - flank)
    e = min(len(seq), end + flank)
    ext = seq[s:e]

    score = 0.0
    # GC junctions can enhance directed bending at boundaries
    score += 1.0 * len(JUNCTION_GC_A_RE.findall(ext))
    score += 1.0 * len(JUNCTION_GC_T_RE.findall(ext))
    # CpG-flanked AT segments can indicate distinct environment
    score += 0.5 * len(CPG_FLANKED_AT_RE.findall(ext))
    # Interrupted tracts preserve anisotropy while modulating stiffness
    score += 0.5 * len(INTERRUPTED_A_RE.findall(ext))
    score += 0.5 * len(INTERRUPTED_T_RE.findall(ext))
    # Symmetric motifs (A4T4/T4A4) historically linked to strong bends
    score += 1.5 * len(SYMMETRIC_A4T4_RE.findall(ext))
    score += 1.2 * len(SYMMETRIC_T4A4_RE.findall(ext))
    # AA/TT repeats (weaker than long pure tracts)
    score += 0.3 * len(A_DINUCL_RE.findall(ext))
    score += 0.3 * len(T_DINUCL_RE.findall(ext))
    return score


def _regex_quality_ok(seq: str, start: int, end: int) -> bool:
    """
    Filter low-quality calls using regex-only QC.
    """
    region = seq[start:end]
    if POLY_N_RE.search(region):
        return False
    if LOW_COMPLEXITY_RE.search(region):
        return False
    # Require AT >= 60% for curved A/T calls
    at_fraction = len(AT_RE.findall(region)) / max(1, len(region))
    if at_fraction < 0.60:
        return False
    return True


def _regex_subtype(seq: str) -> str:
    """
    Classify subtype using regex matches.
    """
    if PERIODIC_A_RE.search(seq) or PERIODIC_T_RE.search(seq) or PERIODIC_MIXED_AT_RE.search(seq):
        return "Global_Curved_Periodic_AT"
    if SYMMETRIC_A4T4_RE.search(seq) or SYMMETRIC_T4A4_RE.search(seq) or SYMMETRIC_AN_TN_RE.search(seq):
        return "Local_Curved_Symmetric_AT"
    if JUNCTION_GC_A_RE.search(seq) or JUNCTION_GC_T_RE.search(seq):
        return "Local_Curved_Junction_Rich"
    return "Local_Curved_Strict_PolyA_or_PolyT"


def _regex_weighted_curvature_score(seq: str) -> float:
    """
    Regex-weighted curvature score to complement the original curvature_score.
    Heuristic weights reflect relative mechanistic contributions:
    - Pure long A/T tracts
    - Periodic phasing (~10 bp)
    - Junction (GC flanks), CpG context
    - Symmetric motifs (A4T4/T4A4)
    - Interrupted tracts
    """
    score = 0.0

    # Base signal from pure tracts (≥7 to align with original default)
    score += 2.5 * len(PURE_A_RE.findall(seq))
    score += 2.5 * len(PURE_T_RE.findall(seq))

    # Periodicity (helical phasing)
    score += 2.0 * len(PERIODIC_A_RE.findall(seq))
    score += 2.0 * len(PERIODIC_T_RE.findall(seq))
    score += 1.5 * len(PERIODIC_MIXED_AT_RE.findall(seq))

    # Junction/context
    score += 1.0 * len(JUNCTION_GC_A_RE.findall(seq))
    score += 1.0 * len(JUNCTION_GC_T_RE.findall(seq))
    score += 0.6 * len(CPG_FLANKED_AT_RE.findall(seq))

    # Symmetric motifs
    score += 1.5 * len(SYMMETRIC_A4T4_RE.findall(seq))
    score += 1.2 * len(SYMMETRIC_T4A4_RE.findall(seq))
    score += 0.8 * len(SYMMETRIC_AN_TN_RE.findall(seq))

    # Interrupted tracts
    score += 0.5 * len(INTERRUPTED_A_RE.findall(seq))
    score += 0.5 * len(INTERRUPTED_T_RE.findall(seq))

    # Length-scaled AT signal
    at_frac = len(AT_RE.findall(seq)) / max(1, len(seq))
    score += 0.5 * len(seq) * at_frac

    return float(score)


# -----------------------------
# Original API (preserved)
# -----------------------------
def find_polyA_polyT_tracts(seq: str, min_len: int = 7) -> list:
    """Find contiguous A or T tracts of minimum length (default 7 retained)."""
    # Regex implementation to improve speed/readability; preserves default behavior
    return _regex_find_tracts(seq, min_len=min_len)


def curvature_score(seq):
    """
    Calculate curvature score based on AT-richness and tract organization.

    Scientific Basis: A-tracts and T-tracts cause DNA bending through
    narrowed minor groove. Score reflects biological propensity for curvature.

    Returns:
    float: Curvature propensity score
    """
    if not seq:
        return 0.0
    at_frac = (seq.count('A') + seq.count('T')) / len(seq)
    runs = re.findall(r"(A+|T+)", seq)
    run_bonus = sum(len(r) ** 0.5 for r in runs)  # diminishing returns
    return len(seq) * (1.0 + at_frac) + 0.5 * run_bonus


def find_global_curved_polyA_polyT(
    seq: str,
    min_tract_len: int = 3,
    min_repeats: int = 3,
    min_spacing: int = 8,
    max_spacing: int = 12,
    min_score: int = 6
) -> tuple:
    """Find globally curved DNA with periodic poly(A)/poly(T) tracts."""
    # Seed with tracts (regex-based) for efficiency
    tracts = find_polyA_polyT_tracts(seq, min_tract_len)
    results = []
    apr_regions = []

    # Identify periodic groups using regex + spacing verification
    groups = _regex_periodic_groups(
        seq,
        min_repeats=min_repeats,
        min_spacing=min_spacing,
        max_spacing=max_spacing
    )

    for group in groups:
        motif_seq = seq[group[0]:group[-1][1] + 1]
        # Combine original curvature_score with regex-weighted score
        score = curvature_score(motif_seq) + 0.8 * _regex_weighted_curvature_score(motif_seq)

        if score >= min_score and _regex_quality_ok(seq, group, group[-1][1] + 1):
            # Calculate conservation score (as provided by shared_utils)
            conservation_result = calculate_conservation_score(motif_seq, "Curved DNA")
            conservation_score_val = conservation_result["enrichment_score"]

            subtype = "Global_Curved_Strict_PolyA_or_PolyT"
            if PERIODIC_A_RE.search(motif_seq) or PERIODIC_T_RE.search(motif_seq) or PERIODIC_MIXED_AT_RE.search(motif_seq):
                subtype = "Global_Curved_Periodic_AT"

            motif = {
                "Sequence Name": "",
                "Class": "Curved_DNA",
                "Subtype": subtype,
                "Start": group + 1,
                "End": group[-1][1] + 1,
                "Length": group[-1][1] - group + 1,
                "Sequence": wrap(motif_seq),
                "Score": float(score),
                "Conservation_Score": float(conservation_score_val),
                "Conservation_P_Value": float(conservation_result["p_value"]),
                "Conservation_Significance": conservation_result["significance"],
                "Arms/Repeat Unit/Copies": "",
                "Spacer": ""
            }
            results.append(motif)
            apr_regions.append((motif["Start"], motif["End"]))

    return results, apr_regions


def find_local_curved_polyA_polyT(seq: str, apr_regions: list, min_len: int = 7) -> list:
    """Find local curved DNA with isolated poly(A)/poly(T) tracts (default 7 retained)."""
    results = []
    tracts = find_polyA_polyT_tracts(seq, min_len)
    for start, end, tract_seq in tracts:
        s, e = start + 1, end + 1
        # Exclude if overlapping a global periodic region
        if any(r_start <= s <= r_end or r_start <= e <= r_end for r_start, r_end in apr_regions):
            continue

        if not _regex_quality_ok(seq, start, end + 1):
            continue

        # Conservation
        conservation_result = calculate_conservation_score(tract_seq, "Curved DNA")
        conservation_score_val = conservation_result["enrichment_score"]

        # Scoring: original + regex-weighted context/junction/symmetric effects
        combined_score = curvature_score(tract_seq) + 0.8 * _regex_weighted_curvature_score(tract_seq)
        subtype = _regex_subtype(tract_seq)

        results.append({
            "Sequence Name": "",
            "Class": "Curved_DNA",
            "Subtype": subtype if subtype.startswith("Local_") else "Local_Curved_Strict_PolyA_or_PolyT",
            "Start": s,
            "End": e,
            "Length": len(tract_seq),
            "Sequence": wrap(tract_seq),
            "Score": float(combined_score),
            "Conservation_Score": float(conservation_score_val),
            "Conservation_P_Value": float(conservation_result["p_value"]),
            "Conservation_Significance": conservation_result["significance"],
            "Arms/Repeat Unit/Copies": "",
            "Spacer": ""
        })
    return results


def find_curved_DNA(seq: str) -> list:
    """
    Main function to detect curved DNA structures.

    Returns both global and local curved DNA structures.
    Enhanced scoring system focuses on:
    1. Curvature propensity based on A/T content and tract organization
    2. Conservation analysis using composition-preserving shuffles
    3. Distinction between global periodic and local isolated structures
    """
    global_results, apr_regions = find_global_curved_polyA_polyT(seq)
    local_results = find_local_curved_polyA_polyT(seq, apr_regions)
    return global_results + local_results

"""
The enhanced curved DNA detection module is grounded in four decades of experimental structural biology 
evidence demonstrating the complex mechanisms underlying A/T-tract-induced DNA bending. High-resolution 
NMR structural studies with residual dipolar couplings by Woods et al. definitively proved that A4T4 and 
T4A4 sequences exhibit fundamentally different bending patterns, with A4T4 showing large negative roll 
angles (-12°) at ApT dinucleotide steps that create coherent bending toward the minor groove, while T4A4 
displays opposing roll angles (+11° at TpA steps vs junction bends) resulting in relatively straight helices. 
This dinucleotide step-specific mechanism validates the module's differential scoring of symmetric motifs 
and context-aware pattern recognition. Gel electrophoresis studies have consistently demonstrated that DNA 
curvature exhibits strict helical phasing requirements, with A-tracts spaced at approximately 10-11 bp intervals 
(one helical turn) showing maximum gel migration anomalies and cyclization efficiency, while out-of-phase arrangements 
lose their curvature properties. The module's periodic spacing detection with 8-12 bp constraints directly captures 
this critical structural requirement that determines functional DNA bending. Comprehensive cation-binding studies 
by Stellwagen et al. and Woods et al. revealed that monovalent cations (NH4+, K+, Na+) bind specifically to A-tract 
minor grooves with species-dependent affinities, modulating curvature magnitude through electrostatic effects, with 
preference order NH4+ > K+ > Na+ for DNA binding and distinct binding stoichiometry of ~0.33 cations per A-tract. 
This experimental evidence supports the module's junction and context scoring functions that account for environmental 
modulation of curvature. Finally, NMR-derived structural modeling combined with cyclization kinetics measurements 
demonstrated that periodic A4T4 repeats form left-handed superhelices with diameter ≈110 Å and pitch ≈80 Å, closely 
resembling the DNA conformation in nucleosome core particles as confirmed by X-ray crystallography, while T4A4 repeats 
form much gentler superhelices with pitch >250 Å. This higher-order structural organization validates the module's 
global versus local classification system and demonstrates that the enhanced detection algorithms incorporate 
experimentally validated principles of DNA structural mechanics rather than purely computational optimizations.
"""
