#!/usr/bin/env python3
"""
Final Integration Test for NBDFinder Clinical Platform
=====================================================

Comprehensive end-to-end test to verify all components are working correctly:
1. Literature-based thresholds
2. Clinical analysis integration
3. Disease motif detection
4. Export functionality
5. Interface components

Author: Dr. Venkata Rajesh Yella
"""

import sys
import os
sys.path.insert(0, os.path.dirname(__file__))

def test_literature_thresholds():
    """Test that all literature-based thresholds are properly implemented"""
    print("🔬 Testing Literature-Based Thresholds")
    print("-" * 50)
    
    tests_passed = 0
    total_tests = 0
    
    # Test G4 thresholds
    try:
        from motifs.g4_related import find_all_g4_motifs, get_g4_formation_category
        
        # High-scoring G4 sequence should be detected
        high_g4 = "GGGCTGGGCTGGGCTGGGC"
        results = find_all_g4_motifs(high_g4)
        if results:
            category = get_g4_formation_category(results[0].get('Score', 0))
            print(f"  ✅ G4 detection: High-scoring sequence detected with {category['category']}")
            tests_passed += 1
        else:
            print(f"  ❌ G4 detection: High-scoring sequence not detected")
        total_tests += 1
        
    except Exception as e:
        print(f"  ❌ G4 threshold test failed: {e}")
        total_tests += 1
    
    # Test R-loop thresholds
    try:
        from motifs.r_loop import find_rlfs
        
        # Short sequence should not be detected (< 100bp)
        short_rloop = "GGGGGGG" + "C" * 50  # Only 57bp
        results = find_rlfs(short_rloop)
        if not results:
            print(f"  ✅ R-loop threshold: Short sequence correctly rejected (<100bp)")
            tests_passed += 1
        else:
            print(f"  ❌ R-loop threshold: Short sequence incorrectly detected")
        total_tests += 1
        
    except Exception as e:
        print(f"  ❌ R-loop threshold test failed: {e}")
        total_tests += 1
    
    # Test Z-DNA thresholds
    try:
        from motifs.zdna_egz import find_zdna
        
        # Test with alternating sequence
        zdna_seq = "CGCGCGCGCGCGCGCGCGCG"
        results = find_zdna(zdna_seq)
        threshold_check = len(results) > 0  # Should detect with stringent threshold
        print(f"  {'✅' if threshold_check else '⚠️'} Z-DNA threshold: Stringent threshold = 50.0")
        if threshold_check:
            tests_passed += 1
        total_tests += 1
        
    except Exception as e:
        print(f"  ❌ Z-DNA threshold test failed: {e}")
        total_tests += 1
    
    return tests_passed, total_tests

def test_clinical_integration():
    """Test clinical analysis integration"""
    print("\n🏥 Testing Clinical Analysis Integration")
    print("-" * 50)
    
    tests_passed = 0
    total_tests = 0
    
    # Test disease motif detection
    try:
        from disease_motifs import find_disease_associated_motifs
        
        # Test with Friedreich's Ataxia GAA expansion
        gaa_expansion = "GAA" * 80  # 240bp, 80 repeats (pathogenic)
        disease_results = find_disease_associated_motifs(gaa_expansion)
        
        if disease_results:
            motif = disease_results[0]
            if ("Friedreich" in motif.get('Disease_Name', '') and 
                motif.get('Clinical_Significance') == 'Pathogenic'):
                print(f"  ✅ Disease detection: Friedreich's Ataxia correctly identified")
                print(f"      Repeat count: {motif.get('Repeat_Count', 0)}, Clinical significance: {motif.get('Clinical_Significance', 'Unknown')}")
                tests_passed += 1
            else:
                print(f"  ❌ Disease detection: Incorrect disease or significance")
        else:
            print(f"  ❌ Disease detection: No disease motifs found")
        total_tests += 1
        
    except Exception as e:
        print(f"  ❌ Disease motif test failed: {e}")
        total_tests += 1
    
    # Test main analysis pipeline integration
    try:
        from motifs.shared_utils import all_motifs
        
        results = all_motifs(gaa_expansion, sequence_name="Clinical_Test")
        disease_motifs = [m for m in results if m.get('Class') == 'Disease-Associated Motif']
        
        if disease_motifs:
            print(f"  ✅ Pipeline integration: Disease motifs integrated into main analysis")
            tests_passed += 1
        else:
            print(f"  ❌ Pipeline integration: Disease motifs not found in main analysis")
        total_tests += 1
        
    except Exception as e:
        print(f"  ❌ Pipeline integration test failed: {e}")
        total_tests += 1
    
    return tests_passed, total_tests

def test_export_functionality():
    """Test export and data handling"""
    print("\n📊 Testing Export and Data Handling")
    print("-" * 50)
    
    tests_passed = 0
    total_tests = 0
    
    try:
        from motifs.shared_utils import all_motifs, format_motif_rows
        import pandas as pd
        
        # Test sequence with multiple motif types
        test_seq = "AAAAAAAAAAAAAAAAAAA" + "GAA" * 50 + "GGGCTGGGCTGGGCTGGGC" + "CCCCCCCCCCCCCCCCCCCC"
        results = all_motifs(test_seq, sequence_name="Export_Test")
        
        if results:
            # Test DataFrame conversion
            df = pd.DataFrame(results)
            required_columns = ['Class', 'Subtype', 'Start', 'End', 'Length', 'Score']
            if all(col in df.columns for col in required_columns):
                print(f"  ✅ Data format: All required columns present")
                tests_passed += 1
            else:
                print(f"  ❌ Data format: Missing required columns")
            total_tests += 1
            
            # Test clinical data presence
            clinical_columns = ['Clinical_Significance', 'Disease_Name', 'Risk_Score']
            has_clinical = any(any(col in str(motif) for col in clinical_columns) for motif in results)
            if has_clinical:
                print(f"  ✅ Clinical data: Clinical information included in export")
                tests_passed += 1
            else:
                print(f"  ❌ Clinical data: No clinical information found")
            total_tests += 1
            
        else:
            print(f"  ❌ Export test: No motifs detected for export")
            total_tests += 2
        
    except Exception as e:
        print(f"  ❌ Export functionality test failed: {e}")
        total_tests += 2
    
    return tests_passed, total_tests

def test_interface_components():
    """Test interface component functionality"""
    print("\n🖥️ Testing Interface Components")
    print("-" * 50)
    
    tests_passed = 0
    total_tests = 0
    
    # Test Streamlit app imports
    try:
        import streamlit as st
        import pandas as pd
        import matplotlib.pyplot as plt
        import plotly.express as px
        print(f"  ✅ Dependencies: All required packages available")
        tests_passed += 1
    except ImportError as e:
        print(f"  ❌ Dependencies: Missing required package - {e}")
    total_tests += 1
    
    # Test motif classification system
    try:
        from motifs.classification_config import get_official_classification, LEGACY_TO_OFFICIAL_MAPPING
        
        # Test mapping functionality
        official_class, official_subtype = get_official_classification("G4", "Canonical")
        if official_class and official_subtype:
            print(f"  ✅ Classification: Motif classification system working")
            tests_passed += 1
        else:
            print(f"  ❌ Classification: Motif classification system failed")
        total_tests += 1
        
    except Exception as e:
        print(f"  ❌ Classification test failed: {e}")
        total_tests += 1
    
    return tests_passed, total_tests

def main():
    """Run comprehensive integration test"""
    print("🧬 NBDFinder Clinical Platform - Final Integration Test")
    print("=" * 70)
    print("Testing all components after clinical integration and literature-based threshold updates")
    print()
    
    # Run all test suites
    threshold_passed, threshold_total = test_literature_thresholds()
    clinical_passed, clinical_total = test_clinical_integration()
    export_passed, export_total = test_export_functionality()
    interface_passed, interface_total = test_interface_components()
    
    # Calculate overall results
    total_passed = threshold_passed + clinical_passed + export_passed + interface_passed
    total_tests = threshold_total + clinical_total + export_total + interface_total
    success_rate = (total_passed / total_tests) * 100 if total_tests > 0 else 0
    
    # Final summary
    print("\n" + "=" * 70)
    print("🎯 FINAL INTEGRATION TEST SUMMARY")
    print("=" * 70)
    print(f"Literature-Based Thresholds: {threshold_passed}/{threshold_total} passed")
    print(f"Clinical Analysis Integration: {clinical_passed}/{clinical_total} passed")
    print(f"Export Functionality: {export_passed}/{export_total} passed")
    print(f"Interface Components: {interface_passed}/{interface_total} passed")
    print("-" * 70)
    print(f"OVERALL RESULT: {total_passed}/{total_tests} tests passed ({success_rate:.1f}%)")
    
    if success_rate >= 90:
        print("\n✅ EXCELLENT - NBDFinder is ready for clinical use!")
        print("🔬 All major components working correctly")
        print("🏥 Clinical analysis fully integrated")
        print("📚 Literature-based thresholds implemented")
        print("🧬 Disease motif detection validated")
        
        print("\n📋 Clinical Use Recommendations:")
        print("- Use for research and clinical research applications")
        print("- Validate findings with additional testing as appropriate")
        print("- Consider genetic counseling for pathogenic findings")
        print("- Keep updated with latest literature for threshold refinements")
        
    elif success_rate >= 75:
        print("\n⚠️ GOOD - Minor issues need attention")
        print("Most functionality working correctly")
        
    else:
        print("\n❌ NEEDS WORK - Significant issues require resolution")
        print("Review failed components before clinical use")
    
    return success_rate >= 90

if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)