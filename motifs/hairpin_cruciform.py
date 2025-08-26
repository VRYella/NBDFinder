"""
Category 2: Hairpin and Cruciform Detection Module
==================================================

Scientific Basis:
-----------------
Cruciform structures extrude from palindromic inverted repeats under negative supercoiling in duplex DNA; two opposing intra-strand hairpins form a four-way junction. 
This requires inverted (palindromic) arms, not direct repeats. See reviews and single-molecule assays. (PMC3176155; NAR 2011 “Cruciform structures are a common DNA feature”;
NAR 2011 real-time extrusion)

Short spacers and longer arms favor extrusion; AT-rich arms can lower local melting barriers but are not sufficient alone. AT is retained as annotation; scoring emphasizes
arm length and spacer length for determinism. (PMC3176155; classical plasmid supercoiling experiments)

Computational detection is therefore palindrome-first: X–spacer–X^rc with arm ≥10 bp and short spacer prioritized; results are reported as maximal, non-overlapping 
candidates with transparent scoring.

Notes (scientific/context):
- Cruciforms are extruded four-way junctions from inverted repeats (palindromes), favored by longer arms and short spacers in supercoiled DNA (duplex context). 
(Review and real-time assays)
- This detector enforces perfect inverted-repeat matching: arm == reverse_complement(opposite arm).
- AT_Fraction is retained as annotation only; scoring emphasizes arm length and spacer length for reproducibility.

Suggested citations:
- Cruciform structures review and in vivo evidence: “Cruciform structures are a common DNA feature important for…”, BMC Mol Biol 2011 (PMC3176155).
- Real-time cruciform extrusion assay: “Real-time detection of cruciform extrusion by single-molecule DNA nanomanipulation”, Nucleic Acids Res 2011.
- Mechanistic background linking palindromic IRs to cruciforms and the role of supercoiling: classic plasmid superhelicity studies summarized in reviews above.

Author: Dr. Venkata Rajesh Yella
Updated: 2025
"""

from .shared_utils import wrap, calculate_conservation_score, reverse_complement

# --- Internal deterministic scoring (mismatch-free, palindrome-driven) ---
def _arm_weight(arm_len):
    # Longer arms increase cruciform propensity under supercoiling; use small integer weights
    # Discretize arm influence to keep scores comparable across candidates
    # Example: scale ~ log-like with caps for stability
    if arm_len >= 100: return 8
    if arm_len >= 70: return 7
    if arm_len >= 50: return 6
    if arm_len >= 30: return 5
    if arm_len >= 20: return 4
    if arm_len >= 15: return 3
    if arm_len >= 12: return 2
    return 1 # minimal qualified arms (>=10)

def _spacer_penalty(spacer_len):
    # Short spacers favor cruciform; penalize as spacer increases
    # Keep penalties small to avoid overwhelming arm weight
    if spacer_len == 0: return 0
    if spacer_len <= 3: return 1
    if spacer_len <= 10: return 2
    return 3 # for permissive searches that go beyond very short spacers

def _length_bonus(total_len):
    # Optional small bonus for overall compact structural segment
    if total_len >= 60: return 2
    if total_len >= 40: return 1
    return 0

def _score_cruciform(arm_len, spacer_len, total_len):
    # Deterministic score: arm contribution + length bonus – spacer penalty
    return _arm_weight(arm_len) + _length_bonus(total_len) - _spacer_penalty(spacer_len)

# --- Cruciform detection - palindrome-enforced, maximal, deterministic scoring ---
def find_cruciform(seq):
    """
    Find cruciform structures formed by inverted repeats with enhanced scoring.
    Notes (scientific/context):
    - Cruciforms are extruded four-way junctions from inverted repeats (palindromes), favored by longer arms and short spacers in supercoiled DNA (duplex context).
    (Review and real-time assays)
    - This detector enforces perfect inverted-repeat matching: arm == reverse_complement(opposite arm).
    - AT_Fraction is retained as annotation only; scoring emphasizes arm length and spacer length for reproducibility.
    """
    results = []
    n = len(seq)
    # Stringent defaults: arms >=10 bp, short spacers prioritized (0–3)
    # Loop can be extended to 0–10 for permissive scans, but short spacers are generally more relevant.
    min_arm = 10
    max_arm = 100  # computational bound; can be raised as needed
    min_spacer = 0
    max_spacer = 3  # small spacer emphasis; configurable

    # Practical minimum total length is thus >= 2*10 = 20 bp
    for i in range(0, max(0, n - 2*min_arm)):
        # For each start, try longer arms first to emit maximal palindrome blocks
        max_arm_here = min(max_arm, (n - i)//2)
        arm_called = False

        for arm_len in range(max_arm_here, min_arm - 1, -1):
            arm = seq[i:i+arm_len]
            if 'n' in arm.lower():
                continue  # skip ambiguous arms

            # Spacer loop
            for spacer_len in range(min_spacer, max_spacer + 1):
                mid = i + arm_len + spacer_len
                end_arm = mid + arm_len
                if end_arm > n:
                    continue

                opp = seq[mid:end_arm]
                # PALINDROME ENFORCEMENT: only allow perfect inverted repeats
                if reverse_complement(arm) != opp:
                    continue

                # Candidate cruciform-forming IR detected
                full = seq[i:end_arm]
                total_len = len(full)

                # Deterministic score (no composition); AT retained as annotation
                at_frac = (arm.count('A') + arm.count('T')) / max(1, arm_len)
                score = float(_score_cruciform(arm_len, spacer_len, total_len))

                # Optional conservation annotation (external to sequence-only call)
                conservation_result = calculate_conservation_score(full, "Cruciform")

                results.append({
                    "Sequence Name": "",
                    "Class": "Cruciform DNA",
                    "Subtype": "Cruciform DNA [IR]/HairPin [IR]",
                    "Start": i+1,
                    "End": end_arm,
                    "Length": total_len,
                    "Sequence": wrap(full),

                    # Updated scoring method emphasizes palindrome geometry and spacer
                    "ScoreMethod": "IR_PerfectPalindrome_v1",
                    "Score": score,

                    # Keep AT_Fraction for annotation only (not used in score)
                    "Arm_Length": arm_len,
                    "AT_Fraction": round(at_frac, 3),

                    "Conservation_Score": float(conservation_result.get("enrichment_score", 0.0)),
                    "Conservation_P_Value": float(conservation_result.get("p_value", 1.0)),
                    "Conservation_Significance": conservation_result.get("significance", ""),

                    "Arms/Repeat Unit/Copies": f"Arms={arm_len}",
                    "Spacer": str(spacer_len)
                })

                arm_called = True
                break  # report the longest qualifying arm at this start (maximal block)

            if arm_called:
                break

    return results
