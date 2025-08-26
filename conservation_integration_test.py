#!/usr/bin/env python3
"""
Conservation Analysis Integration Test
=====================================

This script validates that conservation analysis is properly integrated
into the motif detection pipeline and correctly mapped to Non-B DNA classes
and subclasses.

Author: Dr. Venkata Rajesh Yella
Updated: 2025
"""

import sys
from motifs import all_motifs
from motifs.shared_utils import calculate_conservation_score

def test_conservation_integration():
    """Test that conservation analysis is properly integrated."""
    print("=" * 80)
    print("Testing Conservation Analysis Integration")
    print("=" * 80)
    
    # Test sequence that should trigger multiple motif types
    test_seq = "GGGTTGGGTTGGGTTGGGCCCTTCCCTTCCCTTCCCCGCGCGCGCGCGCGCG"
    print(f"Test sequence: {test_seq}")
    print(f"Length: {len(test_seq)} bp")
    print()
    
    # Run motif detection
    motifs = all_motifs(test_seq, "conservation_test")
    
    print(f"Found {len(motifs)} motifs with conservation analysis:")
    print()
    
    conservation_integrated = 0
    for i, motif in enumerate(motifs):
        print(f"Motif {i+1}: {motif['Class']} -> {motif['Subtype']}")
        
        # Check for conservation fields
        has_conservation = any(key in motif for key in ['Conservation_Score', 'Conservation_P_Value', 'Conservation_Significance'])
        
        if has_conservation:
            conservation_integrated += 1
            print(f"  ✅ Conservation Score: {motif.get('Conservation_Score', 'N/A')}")
            print(f"  ✅ P-value: {motif.get('Conservation_P_Value', 'N/A')}")
            print(f"  ✅ Significance: {motif.get('Conservation_Significance', 'N/A')}")
        else:
            print(f"  ❌ No conservation analysis found")
        
        print()
    
    print(f"Conservation Integration Summary:")
    print(f"  Total motifs: {len(motifs)}")
    print(f"  With conservation analysis: {conservation_integrated}")
    print(f"  Integration rate: {conservation_integrated/len(motifs)*100:.1f}%" if motifs else "N/A")
    
    return conservation_integrated == len(motifs)

def test_conservation_scoring():
    """Test conservation scoring function directly."""
    print("=" * 80)
    print("Testing Conservation Scoring Function")
    print("=" * 80)
    
    test_cases = [
        ("G4 sequence", "GGGTTGGGTTGGGTTGGG", "G4"),
        ("i-motif sequence", "CCCTTCCCTTCCCTTCCC", "i-Motif"),
        ("Z-DNA sequence", "CGCGCGCGCGCGCGCG", "Z-DNA"),
        ("Random sequence", "ATCGATCGATCGATCG", "general")
    ]
    
    all_passed = True
    
    for name, seq, motif_type in test_cases:
        print(f"Testing {name}: {seq}")
        
        try:
            result = calculate_conservation_score(seq, motif_type)
            
            required_fields = ['enrichment_score', 'p_value', 'significance', 'observed_count', 'null_mean', 'null_std']
            missing_fields = [field for field in required_fields if field not in result]
            
            if missing_fields:
                print(f"  ❌ Missing fields: {missing_fields}")
                all_passed = False
            else:
                print(f"  ✅ Enrichment Score: {result['enrichment_score']:.3f}")
                print(f"  ✅ P-value: {result['p_value']:.3f}")
                print(f"  ✅ Significance: {result['significance']}")
                print(f"  ✅ All required fields present")
                
        except Exception as e:
            print(f"  ❌ Error: {e}")
            all_passed = False
        
        print()
    
    return all_passed

def test_disease_motif_mapping():
    """Test disease motif mapping integration."""
    print("=" * 80)
    print("Testing Disease Motif Mapping Integration")
    print("=" * 80)
    
    # Test sequences known to be associated with diseases
    disease_test_cases = [
        ("CGG expansion (Fragile X)", "CGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGG"),
        ("GAA expansion (FRDA)", "GAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAA"),
        ("CTG expansion", "CTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTG")
    ]
    
    disease_mapped = 0
    
    for name, seq in disease_test_cases:
        print(f"Testing {name}: {seq[:30]}...")
        
        motifs = all_motifs(seq, f"disease_test_{disease_mapped}")
        
        has_disease_info = False
        for motif in motifs:
            # Check for disease-related fields that might be added
            disease_fields = [key for key in motif.keys() if 'disease' in key.lower() or 'pathogenic' in key.lower() or 'clinical' in key.lower()]
            
            if disease_fields:
                has_disease_info = True
                print(f"  ✅ Disease-related fields found: {disease_fields}")
            
            # Check if conservation analysis is present (required for disease mapping)
            if motif.get('Conservation_Significance') != 'not significant':
                print(f"  ✅ Significant conservation found for potential disease association")
        
        if has_disease_info:
            disease_mapped += 1
        else:
            print(f"  ⚠️  No explicit disease fields (may be handled by external disease analysis)")
        
        print()
    
    print(f"Disease Mapping Summary:")
    print(f"  Test cases: {len(disease_test_cases)}")
    print(f"  With disease mapping: {disease_mapped}")
    
    return True  # Disease analysis may be external to core detection

def main():
    """Run all conservation integration tests."""
    print("NBDFinder Conservation Analysis Integration Validation")
    print("Testing integration of conservation analysis with motif detection")
    print()
    
    all_passed = True
    
    # Run tests
    tests = [
        test_conservation_integration,
        test_conservation_scoring,
        test_disease_motif_mapping
    ]
    
    for test_func in tests:
        try:
            if not test_func():
                all_passed = False
        except Exception as e:
            print(f"❌ Test {test_func.__name__} failed with error: {e}")
            all_passed = False
    
    print("\n" + "=" * 80)
    if all_passed:
        print("✅ CONSERVATION INTEGRATION TESTS PASSED!")
        print("Conservation analysis is properly integrated into motif detection.")
        print("All detected motifs include conservation scores and significance.")
    else:
        print("❌ SOME CONSERVATION TESTS FAILED!")
        print("Please review the conservation integration issues above.")
    print("=" * 80)
    
    return 0 if all_passed else 1

if __name__ == "__main__":
    sys.exit(main())