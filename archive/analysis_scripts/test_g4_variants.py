#!/usr/bin/env python3
"""
Comprehensive test for all G-quadruplex variants specified in the problem statement.
Tests: Canonical G4, Relaxed G4, Bulged G4, Bipartite G4, Multimeric G4, Imperfect G4
"""

import motifs

def test_g4_variants():
    """Test all G4 variants required by the problem statement."""
    print("="*70)
    print("Testing All G-Quadruplex Variants")
    print("="*70)
    
    # Test cases for each G4 type
    test_cases = [
        ("Canonical G4", "GGGTTAGGGTTAGGGTTAGGG", "find_gquadruplex"),
        ("Relaxed G4", "GGGAAAAAAAAGGGAAAAAAAAGGGAAAAAAAAGGG", "find_relaxed_gquadruplex"),
        ("Bulged G4", "GGGAGGGTGGGTGGG", "find_bulged_gquadruplex"),
        ("Bipartite G4", "GGGAAAGGGAAAGGGAAAGGGAAAAAAAAAAAAGGGAAAGGGAAAGGGAAAGGG", "find_bipartite_gquadruplex"),
        ("Multimeric G4", "GGGAAAGGGAAAGGGAAAGGGAAAGGGAAAGGGAAAGGGAAAGGG", "find_multimeric_gquadruplex"),
        ("Imperfect G4", "GGAAAGGGAAAGGGAAAGGG", "find_imperfect_gquadruplex"),
    ]
    
    for motif_type, sequence, function_name in test_cases:
        print(f"\n{motif_type}:")
        print(f"  Sequence: {sequence}")
        
        # Test individual function
        func = getattr(motifs, function_name)
        individual_results = func(sequence)
        print(f"  Individual function results: {len(individual_results)}")
        
        # Test all_motifs integration
        all_results = motifs.all_motifs(sequence)
        g4_results = [m for m in all_results if motif_type.replace(" ", " ") in m.get('Class', '')]
        print(f"  all_motifs() integration: {len(g4_results)}")
        
        if individual_results:
            for result in individual_results[:1]:  # Show first result
                print(f"    Class: {result.get('Class', 'N/A')}")
                print(f"    Score: {result.get('Score', 'N/A')}")
                print(f"    G4Hunter Mean: {result.get('G4Hunter_Mean', 'N/A')}")
        
        if not individual_results:
            print(f"    ⚠️  No results - may need threshold adjustment")
        else:
            print(f"    ✅ Working correctly")
    
    print(f"\n{'-'*50}")
    print("Testing comprehensive sequence with all G4 types:")
    
    # Comprehensive sequence containing multiple G4 types
    comprehensive_seq = (
        "GGGTTAGGGTTAGGGTTAGGG"  # Canonical
        "AAAAAA"
        "GGGTTTAGGGAAAGGGTTTAGGG"  # Relaxed
        "AAAAAA"
        "GGGAGGGTGGGTGGG"  # Bulged
        "AAAAAA"
        "GGAAAGGGAAAGGGAAAGGG"  # Imperfect
    )
    
    print(f"Sequence length: {len(comprehensive_seq)} bp")
    all_results = motifs.all_motifs(comprehensive_seq)
    
    # Count each G4 type
    g4_counts = {}
    for result in all_results:
        class_name = result.get('Class', 'Other')
        if 'G4' in class_name or 'Canonical' in class_name:
            g4_counts[class_name] = g4_counts.get(class_name, 0) + 1
    
    print("G4 types detected:")
    for g4_type, count in g4_counts.items():
        print(f"  {g4_type}: {count}")
    
    total_g4s = sum(g4_counts.values())
    print(f"Total G4 structures detected: {total_g4s}")
    
    # Test categorization according to problem statement
    print(f"\n{'-'*50}")
    print("Problem Statement Categorization Verification:")
    
    expected_categories = [
        "Canonical G4", "Relaxed G4", "Bulged G4", 
        "Bipartite G4", "Multimeric G4", "Imperfect G4"
    ]
    
    found_categories = set()
    for result in all_results:
        class_name = result.get('Class', '')
        if any(cat in class_name for cat in expected_categories):
            found_categories.add(class_name)
    
    print("Categories found in results:")
    for category in expected_categories:
        status = "✅" if category in found_categories else "❌"
        print(f"  {status} {category}")
    
    missing = set(expected_categories) - found_categories
    if missing:
        print(f"\nMissing categories: {missing}")
    else:
        print(f"\n🎉 All G4 categories successfully implemented!")
    
    print("="*70)

if __name__ == "__main__":
    test_g4_variants()