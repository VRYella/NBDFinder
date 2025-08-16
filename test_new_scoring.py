#!/usr/bin/env python3
"""
Test script to validate the new G4Hunter run-based scoring system and i-motif reverse logic.
"""

import motifs

def test_g4hunter_scoring():
    """Test the new G4Hunter run-based scoring function."""
    print("=" * 70)
    print("Testing G4Hunter Run-Based Scoring")
    print("=" * 70)
    
    test_cases = [
        # (sequence, expected_score, description)
        ("GGG", 3.0, "3 G-run: each G gets score 3"),
        ("GGGG", 4.0, "4 G-run: each G gets score 4"),
        ("GGGGG", 4.0, "5 G-run: each G gets score 4 (capped)"),
        ("CCC", -3.0, "3 C-run: each C gets score -3"),
        ("CCCC", -4.0, "4 C-run: each C gets score -4"),
        ("CCCCC", -4.0, "5 C-run: each C gets score -4 (capped)"),
        ("ATATAT", 0.0, "A/T only: should be 0"),
        ("GGGCCC", 0.0, "Equal G and C runs: should cancel out"),
        ("GGGAT", 1.8, "Mixed: (3+3+3+0+0)/5 = 1.8"),
        ("GGGTTAGGGTTAGGGTTAGGG", None, "Canonical G4 pattern"),
    ]
    
    print("Testing individual scoring patterns:")
    for seq, expected, desc in test_cases:
        score = motifs.g4hunter_score(seq)
        if expected is not None:
            status = "✅" if abs(score - expected) < 0.01 else "❌"
            print(f"  {status} {seq:25} -> {score:.4f} (expected {expected:.4f}) - {desc}")
        else:
            print(f"  ℹ️  {seq:25} -> {score:.4f} - {desc}")
    
    print("\n" + "=" * 70)
    print("Testing i-Motif Reverse Logic Scoring")
    print("=" * 70)
    
    test_cases_imotif = [
        # (sequence, expected_score, description)
        ("CCC", 3.0, "3 C-run: each C gets score 3 (reverse logic)"),
        ("CCCC", 4.0, "4 C-run: each C gets score 4 (reverse logic)"),
        ("CCCCC", 4.0, "5 C-run: each C gets score 4 (capped)"),
        ("GGG", -3.0, "3 G-run: each G gets score -3 (reverse logic)"),
        ("GGGG", -4.0, "4 G-run: each G gets score -4 (reverse logic)"),
        ("GGGGG", -4.0, "5 G-run: each G gets score -4 (capped)"),
        ("ATATAT", 0.0, "A/T only: should be 0"),
        ("CCCGGG", 0.0, "Equal C and G runs: should cancel out"),
        ("CCCTAACCCTAACCCTAACCC", None, "Canonical i-motif pattern"),
    ]
    
    print("Testing individual scoring patterns:")
    for seq, expected, desc in test_cases_imotif:
        score = motifs.imotif_score(seq)
        if expected is not None:
            status = "✅" if abs(score - expected) < 0.01 else "❌"
            print(f"  {status} {seq:25} -> {score:.4f} (expected {expected:.4f}) - {desc}")
        else:
            print(f"  ℹ️  {seq:25} -> {score:.4f} - {desc}")

def test_scoring_comparison():
    """Compare the old vs new scoring behavior."""
    print("\n" + "=" * 70)
    print("Scoring Behavior Comparison")
    print("=" * 70)
    
    test_sequences = [
        "GGGTTAGGGTTAGGGTTAGGG",  # Canonical G4
        "CCCTAACCCTAACCCTAACCC",  # Canonical i-motif
        "GGG",                    # Short G-run
        "CCC",                    # Short C-run
        "GGGGGGGG",              # Long G-run
        "CCCCCCCC",              # Long C-run
    ]
    
    print("Sequence analysis:")
    print(f"{'Sequence':<25} {'G4Hunter':<12} {'i-Motif':<12} {'Notes'}")
    print("-" * 65)
    
    for seq in test_sequences:
        g4_score = motifs.g4hunter_score(seq)
        im_score = motifs.imotif_score(seq)
        
        # Determine motif preference
        if g4_score > 1.0:
            notes = "G4 favorable"
        elif im_score > 1.0:
            notes = "i-motif favorable"
        elif abs(g4_score) < 0.1 and abs(im_score) < 0.1:
            notes = "Neutral"
        else:
            notes = "Mixed"
            
        print(f"{seq:<25} {g4_score:<12.4f} {im_score:<12.4f} {notes}")

if __name__ == "__main__":
    test_g4hunter_scoring()
    test_scoring_comparison()
    print("\n" + "=" * 70)
    print("🎉 All scoring tests completed!")
    print("=" * 70)