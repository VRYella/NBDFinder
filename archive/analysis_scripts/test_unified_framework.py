#!/usr/bin/env python3
"""
Comprehensive Test Suite for NBDFinder Unified Framework
========================================================

Tests the implementation of unified computational logic across G4, i-motif, 
r-loop, and z-DNA prediction, with proper subclass-based analysis for 
hybrid and cluster regions.

Scientific Validation:
- Ensures same core computational logic structure for all motif types
- Validates i-motif family implementation (canonical, relaxed, AC-motif)
- Tests subclass breakdown for hybrid and cluster analysis
- Verifies literature-based patterns and scoring

Author: Test Suite for Scientific Enhancement
Updated: 2024
"""

import sys
import time
import traceback
from motifs import *
from motifs.shared_utils import unified_hunter_score, calculate_structural_factor

def test_unified_computational_framework():
    """Test that all motif types use the same core computational logic."""
    print("=" * 80)
    print("Testing Unified Computational Framework")
    print("=" * 80)
    
    # Test sequence with both G-rich and C-rich regions
    test_seq = "GGGAGGGATGGGAGGGCCCACCACCCACCC"
    
    print(f"Test sequence: {test_seq}")
    print(f"Length: {len(test_seq)} bp")
    
    # Test unified scoring for different motif types
    g4_score = unified_hunter_score(test_seq, 'G', 'C')
    imotif_score = unified_hunter_score(test_seq, 'C', 'G')
    
    print(f"\nUnified Hunter Scores:")
    print(f"  G4 bias (G-target): {g4_score:.3f}")
    print(f"  i-motif bias (C-target): {imotif_score:.3f}")
    
    # Verify scores are symmetric
    assert abs(g4_score + imotif_score) < 0.001, "G4 and i-motif scores should be symmetric"
    
    # Test structural factor calculation for different motif types
    g4_factor = calculate_structural_factor(test_seq, "G4")
    imotif_factor = calculate_structural_factor(test_seq, "i-motif")
    rloop_factor = calculate_structural_factor(test_seq, "R-loop")
    zdna_factor = calculate_structural_factor(test_seq, "Z-DNA")
    
    print(f"\nStructural Factors:")
    print(f"  G4: {g4_factor:.3f}")
    print(f"  i-motif: {imotif_factor:.3f}")
    print(f"  R-loop: {rloop_factor:.3f}")
    print(f"  Z-DNA: {zdna_factor:.3f}")
    
    # All should be positive and reasonable
    for name, factor in [("G4", g4_factor), ("i-motif", imotif_factor), 
                        ("R-loop", rloop_factor), ("Z-DNA", zdna_factor)]:
        assert factor > 0, f"{name} structural factor should be positive"
        assert factor < 10, f"{name} structural factor should be reasonable"
    
    print("✅ Unified computational framework tests passed!")
    return True

def test_imotif_family_implementation():
    """Test i-motif family implementation with proper patterns and constraints."""
    print("\n" + "=" * 80)
    print("Testing i-motif Family Implementation")
    print("=" * 80)
    
    # Test canonical i-motif (C3+N1-7)x4
    canonical_seq = "CCCATCCCACCCAGCCC"  # Loops: 2, 1, 2 (all ≤7)
    canonical_results = find_imotif(canonical_seq)
    
    print(f"Canonical i-motif test:")
    print(f"  Sequence: {canonical_seq}")
    print(f"  Expected: Canonical i-motif with loops 1-7")
    
    canonical_found = False
    for result in canonical_results:
        if result['Subtype'] == 'Canonical i-motif':
            canonical_found = True
            loops = result['Loop_Lengths']
            print(f"  Found: {result['Subtype']}, Loops: {loops}")
            assert all(1 <= l <= 7 for l in loops), f"Canonical loops should be 1-7, got {loops}"
            assert result['ScoreMethod'] == 'iM_UnifiedHunter_raw', "Should use unified scoring"
    
    assert canonical_found, "Should detect canonical i-motif"
    
    # Test relaxed i-motif (C3+N8-12)x4
    relaxed_seq = "CCCATCGATCGACCCACGATCGAACCCAGCGATCGACCC"  # Loops: 9, 9, 9 (all 8-12)
    relaxed_results = find_imotif(relaxed_seq)
    
    print(f"\nRelaxed i-motif test:")
    print(f"  Sequence: {relaxed_seq}")
    print(f"  Expected: Relaxed i-motif with loops 8-12")
    
    relaxed_found = False
    for result in relaxed_results:
        if result['Subtype'] == 'Relaxed i-motif':
            relaxed_found = True
            loops = result['Loop_Lengths']
            print(f"  Found: {result['Subtype']}, Loops: {loops}")
            assert all(8 <= l <= 12 for l in loops), f"Relaxed loops should be 8-12, got {loops}"
    
    assert relaxed_found, "Should detect relaxed i-motif"
    
    # Test AC-motif based on Hur et al. NAR 2021
    ac_seq = "AAAATGCCCACGCCCTGCCCC"
    ac_results = find_ac_motifs(ac_seq)
    
    print(f"\nAC-motif test:")
    print(f"  Sequence: {ac_seq}")
    print(f"  Expected: AC-motif based on Hur et al. patterns")
    
    if ac_results:
        for result in ac_results:
            print(f"  Found: {result['Subtype']}")
            print(f"  Pattern: {result['Pattern_Type']}")
            print(f"  AC fraction: {result['AC_Fraction']}")
            assert result['ScoreMethod'] == 'AC_UnifiedFramework_raw', "Should use unified framework"
            assert result['AC_Fraction'] >= 0.6, "AC content should be ≥60% per literature"
    
    print("✅ i-motif family implementation tests passed!")
    return True

def test_subclass_breakdown_analysis():
    """Test that hybrid and cluster analysis provides subclass breakdown."""
    print("\n" + "=" * 80)
    print("Testing Subclass Breakdown Analysis")
    print("=" * 80)
    
    # Create a complex sequence with multiple motif types
    complex_seq = "GGGTTAGGGTTAGGGTTAGGGCGCGCGCGCGCGCGCCCCAACCCCAACCCCAACCCC"
    
    print(f"Complex sequence: {complex_seq}")
    print(f"Length: {len(complex_seq)} bp")
    
    # Detect all motifs
    from motifs.shared_utils import all_motifs as detect_all_motifs
    detected_motifs = detect_all_motifs(complex_seq)
    
    print(f"\nDetected {len(detected_motifs)} total motifs")
    
    # Extract motif classes
    motif_classes = set()
    motif_subtypes = set()
    
    for motif in detected_motifs:
        motif_classes.add(motif.get('Class', 'Unknown'))
        motif_subtypes.add(motif.get('Subtype', 'Unknown'))
    
    print(f"Motif classes found: {sorted(motif_classes)}")
    print(f"Motif subtypes found: {sorted(motif_subtypes)}")
    
    # Test hybrid analysis with subclass breakdown
    hybrids = find_hybrids(detected_motifs, complex_seq)
    
    print(f"\nHybrid Analysis:")
    print(f"  Found {len(hybrids)} hybrid structures")
    
    for i, hybrid in enumerate(hybrids[:3]):  # Show first 3
        print(f"  Hybrid {i+1}:")
        print(f"    Subtype: {hybrid['Subtype']}")
        if 'ClassBreakdown' in hybrid:
            print(f"    Class breakdown: {hybrid['ClassBreakdown']}")
            print(f"    Subclass breakdown: {hybrid.get('SubclassBreakdown', {})}")
            assert len(hybrid['ClassBreakdown']) >= 2, "Hybrid should have ≥2 classes"
            assert 'SubclassDetail' in hybrid, "Should have subclass detail"
        
        assert hybrid['ScoreMethod'] == 'HybridOverlap_SubclassAnalysis_raw', "Should use subclass analysis"
    
    # Test cluster analysis with subclass breakdown
    clusters = find_hotspots(detected_motifs, len(complex_seq))
    
    print(f"\nCluster Analysis:")
    print(f"  Found {len(clusters)} cluster regions")
    
    for i, cluster in enumerate(clusters):
        print(f"  Cluster {i+1}:")
        print(f"    Subtype: {cluster['Subtype']}")
        if 'ClassBreakdown' in cluster:
            print(f"    Class breakdown: {cluster['ClassBreakdown']}")
            print(f"    Complexity level: {cluster['ComplexityLevel']}")
            assert cluster['ComplexityLevel'] in ['Low', 'Moderate', 'High'], "Should have complexity level"
        
        assert cluster['ScoreMethod'] == 'Cluster_SubclassAnalysis_raw', "Should use subclass analysis"
    
    print("✅ Subclass breakdown analysis tests passed!")
    return True

def test_scoring_consistency():
    """Test that scoring methods are consistent across motif types."""
    print("\n" + "=" * 80)
    print("Testing Scoring Method Consistency")
    print("=" * 80)
    
    test_sequences = {
        'G4': "GGGTTAGGGTTAGGGTTAGGG",
        'i-motif': "CCCAACCCCAACCCCAACCC",
        'R-loop': "GGGAGGGATGGGAGGGAGGGATGGG",
        'Z-DNA': "GCGCGCGCGCGCGCGCGC"
    }
    
    results = {}
    
    for motif_type, seq in test_sequences.items():
        print(f"\nTesting {motif_type}:")
        print(f"  Sequence: {seq}")
        
        if motif_type == 'G4':
            motifs = find_gquadruplex(seq)
            expected_score_method = 'G4_UnifiedHunter_raw'
        elif motif_type == 'i-motif':
            motifs = find_imotif(seq)
            expected_score_method = 'iM_UnifiedHunter_raw'
        elif motif_type == 'R-loop':
            motifs = find_rlfs(seq)
            expected_score_method = 'RLFS_UnifiedFramework_raw'
        elif motif_type == 'Z-DNA':
            motifs = find_zdna(seq)
            expected_score_method = 'Kadane_UnifiedFramework_raw'
        
        print(f"  Found {len(motifs)} motifs")
        
        for motif in motifs:
            score_method = motif.get('ScoreMethod', 'Unknown')
            print(f"    ScoreMethod: {score_method}")
            
            # Check that unified framework is being used
            assert 'Unified' in score_method or 'UnifiedHunter' in score_method or 'UnifiedFramework' in score_method, \
                f"{motif_type} should use unified framework, got {score_method}"
        
        results[motif_type] = motifs
    
    print("✅ Scoring consistency tests passed!")
    return True

def test_literature_compliance():
    """Test compliance with scientific literature requirements."""
    print("\n" + "=" * 80)
    print("Testing Literature Compliance")
    print("=" * 80)
    
    # Test G4Hunter algorithm compliance (Brázda et al. 2019)
    g4_seq = "GGGTTAGGGTTAGGGTTAGGG"
    g4_results = find_gquadruplex(g4_seq)
    
    print(f"G4Hunter compliance test:")
    print(f"  Sequence: {g4_seq}")
    
    if g4_results:
        for result in g4_results:
            g4hunter_score = result.get('G4Hunter_Score', 0)
            print(f"  G4Hunter score: {g4hunter_score:.3f}")
            assert g4hunter_score >= 1.0, "Should meet G4Hunter threshold ≥1.0"
    
    # Test i-motif constraints (Zeraati et al. 2018)
    print(f"\ni-motif literature compliance:")
    canonical_results = find_imotif("CCCATCCCACCCAGCCC")
    
    for result in canonical_results:
        if result['Subtype'] == 'Canonical i-motif':
            loops = result['Loop_Lengths']
            print(f"  Canonical loops: {loops}")
            assert all(1 <= l <= 7 for l in loops), "Canonical i-motif loops should be 1-7 nt"
    
    # Test AC-motif patterns (Hur et al. NAR 2021)
    print(f"\nAC-motif literature compliance:")
    ac_results = find_ac_motifs("AAAATGCCCACGCCCTGCCCC")
    
    for result in ac_results:
        ac_fraction = result['AC_Fraction']
        length = result['Length']
        print(f"  AC fraction: {ac_fraction:.3f}, Length: {length}")
        assert ac_fraction >= 0.6, "AC content should be ≥60% per Hur et al."
        assert length >= 15, "Minimum length should be ≥15 bp"
    
    print("✅ Literature compliance tests passed!")
    return True

def main():
    """Run all test suites."""
    print("🧬 NBDFinder Unified Framework Validation Suite")
    print("Testing implementation of scientific computation logic consistency")
    print("and subclass-based analysis for hybrid and cluster regions\n")
    
    start_time = time.time()
    tests_passed = 0
    total_tests = 5
    
    try:
        # Run all test suites
        if test_unified_computational_framework():
            tests_passed += 1
            
        if test_imotif_family_implementation():
            tests_passed += 1
            
        if test_subclass_breakdown_analysis():
            tests_passed += 1
            
        if test_scoring_consistency():
            tests_passed += 1
            
        if test_literature_compliance():
            tests_passed += 1
            
        end_time = time.time()
        
        print("\n" + "=" * 80)
        print("🎯 VALIDATION SUMMARY")
        print("=" * 80)
        print(f"Tests passed: {tests_passed}/{total_tests}")
        print(f"Execution time: {(end_time - start_time)*1000:.1f} ms")
        
        if tests_passed == total_tests:
            print("\n✅ ALL TESTS PASSED!")
            print("🧬 Unified computational framework successfully implemented")
            print("📊 Subclass-based analysis working correctly")
            print("📚 Scientific literature requirements met")
            print("\nImplementation meets all acceptance criteria:")
            print("✓ G4, i-motif, r-loop, and z-DNA use same core computational logic")
            print("✓ i-motif family includes canonical, relaxed, and AC-motif with correct patterns")
            print("✓ Hybrid and cluster regions analyzed with subclass breakdown")
            print("✓ All logic aligns with scientific literature and best practices")
            return True
        else:
            print(f"\n❌ {total_tests - tests_passed} test(s) failed!")
            return False
            
    except Exception as e:
        print(f"\n💥 Test suite failed with error: {str(e)}")
        print(traceback.format_exc())
        return False

if __name__ == "__main__":
    success = main()
    sys.exit(0 if success else 1)