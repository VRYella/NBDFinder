#!/usr/bin/env python3
"""
Final validation test for all Non-B DNA motif categories as specified in the problem statement.
"""

import motifs

def validate_problem_statement_requirements():
    """Validate all motif categories specified in the problem statement."""
    print("="*80)
    print("PROBLEM STATEMENT VALIDATION TEST")
    print("="*80)
    
    # Create a comprehensive test sequence containing representatives of all categories
    test_sequence = (
        # G-quadruplex-related
        "GGGTTAGGGTTAGGGTTAGGG"  # Canonical G4
        "AAAAAA"
        "GGGAAAAAAAAGGGAAAAAAAAGGGAAAAAAAAGGG"  # Relaxed G4
        "AAAAAA"
        "GGGAGGGTGGGTGGG"  # Bulged G4
        "AAAAAA"
        "GGGAAAGGGAAAGGGAAAGGGAAAAAAAAAAAAGGGAAAGGGAAAGGGAAAGGG"  # Bipartite G4
        "AAAAAA"
        "GGGAAAGGGAAAGGGAAAGGGAAAGGGAAAGGGAAAGGGAAAGGG"  # Multimeric G4
        "AAAAAA"
        "GGAAAGGGAAAGGGAAAGGG"  # Imperfect G4
        "AAAAAA"
        "GGGAAAGGGAAAGGG"  # G-Triplex
        "AAAAAA"
        "CCCTTACCCTTACCCTTACCC"  # i-Motif
        "AAAAAA"
        
        # Helix/curvature
        "CGCGCGCGCGCGCGCGCGCGCGCGCGCG"  # Z-DNA
        "AAAAAA"
        "CGGCGGCGGCGGCGGCGG"  # eGZ (Extruded-G)
        "AAAAAA"
        "AAAAAAATTTTTTTAAAAAAAATTTTTTT"  # Curved DNA
        "AAAAAA"
        "AAACCCAAACCCAAACCCAAA"  # AC-Motif
        "AAAAAA"
        
        # Repeat/junction
        "ATCGATCGATCGATCGATCGATCGATCG"  # Slipped DNA
        "AAAAAA"
        "ATCGATCGATCGAAACGATCGATCGAT"  # Cruciform
        "AAAAAA"
        "GAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAA"  # Sticky DNA
        "AAAAAA"
        "AAAAAAAAGGGGGGGGAAAAAAAA"  # Triplex DNA  
        "AAAAAA"
        
        # Hybrid/cluster (these are detected from overlaps)
        "GGGTTTGGGTTTGGGTTTGGGAAACCCAAACCCAAACCCAAA"  # R-Loop attempt
    )
    
    print(f"Test sequence length: {len(test_sequence)} bp")
    print("\nRunning comprehensive motif detection...")
    
    # Run the analysis (include hotspots to test Non-B DNA Clusters)
    all_results = motifs.all_motifs(test_sequence, report_hotspots=True)
    
    # Organize by categories from problem statement (using actual class names)
    categories = {
        "G-quadruplex-related": ["Canonical G4", "Relaxed G4", "Bulged G4", "Bipartite G4", 
                                "Multimeric G4", "Imperfect G4", "G-Triplex", "i-Motif", "Hybrid"],
        "Helix/curvature": ["Z-DNA", "eGZ (Extruded-G)", "Curved_DNA", "AC-Motif"],
        "Repeat/junction": ["Slipped_DNA", "Cruciform", "Sticky_DNA", "Triplex_DNA"],
        "Hybrid/cluster": ["R-Loop", "Non-B DNA Clusters", "Hybrid"]
    }
    
    print(f"\nTotal motifs detected: {len(all_results)}")
    print("\n" + "="*60)
    print("CATEGORY BREAKDOWN:")
    print("="*60)
    
    for category_name, motif_types in categories.items():
        print(f"\n{category_name.upper()}:")
        print("-" * (len(category_name) + 1))
        
        category_count = 0
        found_types = set()
        
        for result in all_results:
            motif_class = result.get('Class', '')
            if motif_class in motif_types:
                category_count += 1
                found_types.add(motif_class)
        
        for motif_type in motif_types:
            status = "✅" if motif_type in found_types else "❌"
            count = len([r for r in all_results if r.get('Class') == motif_type])
            print(f"  {status} {motif_type}: {count} detected")
        
        print(f"  → Total in category: {category_count}")
    
    # Summary statistics
    print("\n" + "="*60)
    print("SUMMARY STATISTICS:")
    print("="*60)
    
    all_expected_types = []
    for types in categories.values():
        all_expected_types.extend(types)
    
    unique_expected = set(all_expected_types)
    found_types = set(result.get('Class', '') for result in all_results)
    successfully_detected = found_types.intersection(unique_expected)
    
    print(f"Expected motif types: {len(unique_expected)}")
    print(f"Successfully detected: {len(successfully_detected)}")
    print(f"Detection rate: {len(successfully_detected)/len(unique_expected)*100:.1f}%")
    
    missing = unique_expected - found_types
    if missing:
        print(f"Missing types: {missing}")
    
    # Performance check
    import time
    start_time = time.time()
    motifs.all_motifs(test_sequence)
    end_time = time.time()
    
    print(f"Performance: {(end_time - start_time)*1000:.2f} ms for {len(test_sequence)} bp sequence")
    
    print("\n" + "="*80)
    if len(missing) == 0:
        print("🎉 ALL PROBLEM STATEMENT REQUIREMENTS SUCCESSFULLY IMPLEMENTED!")
    else:
        print(f"⚠️  {len(missing)} motif types need additional tuning")
    print("="*80)

if __name__ == "__main__":
    validate_problem_statement_requirements()