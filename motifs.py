def find_ac_motifs(seq):
    """
    Finds AC-motif consensus sequences in a DNA string.

    Consensus Definition: Hur et al., Nucleic Acids Res., 2021, 49:10150–10165
    - Pattern: A3 N{4-6} C3 N{4-6} C3 N{4-6} C3
      or     : C3 N{4-6} C3 N{4-6} C3 N{4-6} A3
    Returns:
        list of motif dicts.
    """
    pattern = re.compile(
        r"(?=(?:A{3}[ACGT]{4,6}C{3}[ACGT]{4,6}C{3}[ACGT]{4,6}C{3}|"
        r"C{3}[ACGT]{4,6}C{3}[ACGT]{4,6}C{3}[ACGT]{4,6}A{3}))",
        re.IGNORECASE
    )
    results = []
    for m in pattern.finditer(seq):
        motif_seq = m.group(0).upper()
        results.append({
            "Class": "AC-Motif",
            "Subtype": "Consensus",
            "Start": m.start() + 1,
            "End": m.start() + len(motif_seq),
            "Length": len(motif_seq),
            "Sequence": wrap(motif_seq),
            "ScoreMethod": "PatternMatch",
            "Score": "1.0"
        })
    return results
