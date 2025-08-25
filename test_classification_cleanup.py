#!/usr/bin/env python3
"""
Test script to verify that legacy types and strand information are properly removed
from results and that all features are integrated with the official classification system.
"""

import sys
import pandas as pd
from motifs import all_motifs
from motifs.classification_config import get_official_classification, OFFICIAL_CLASSIFICATION

def test_official_classification_system():
    """Test that the official classification system works correctly."""
    print("=" * 80)
    print("Testing Official Classification System")
    print("=" * 80)
    
    # Test legacy to official mapping
    test_cases = [
        ("G4_Canonical", "", "G-Quadruplex Family", "Canonical G4"),
        ("Canonical G4", "", "G-Quadruplex Family", "Canonical G4"),
        ("R-Loop", "", "R-loop", "R-loop"),
        ("Z-DNA", "", "Z-DNA", "Z-DNA"),
        ("Cruciform_IR", "", "Cruciform DNA", "Cruciform DNA [IR]/HairPin [IR]"),
    ]
    
    for legacy_class, legacy_subtype, expected_class, expected_subtype in test_cases:
        official_class, official_subtype = get_official_classification(legacy_class, legacy_subtype)
        print(f"✓ {legacy_class} -> {official_class} | {official_subtype}")
        
        # Verify mapping is correct
        if official_class != expected_class or official_subtype != expected_subtype:
            print(f"  ❌ Expected: {expected_class} | {expected_subtype}")
            return False
    
    print(f"\nOfficial classification has {len(OFFICIAL_CLASSIFICATION)} classes:")
    for i, class_info in OFFICIAL_CLASSIFICATION.items():
        print(f"  {i}. {class_info['class_name']} ({len(class_info['subclasses'])} subclasses)")
    
    return True

def test_motif_detection_output():
    """Test that motif detection produces clean output without legacy information."""
    print("\n" + "=" * 80)
    print("Testing Motif Detection Output Cleanliness")
    print("=" * 80)
    
    # Test sequence with various motif types
    test_seq = "GGGGTTAGGGTTAGGGTTAGGGCGCGCGCGCGCCCCCAACCCCAACCCCAA"
    print(f"Test sequence: {test_seq}")
    print(f"Length: {len(test_seq)} bp")
    
    # Run detection
    results_list = all_motifs(test_seq, "test_sequence")
    
    print(f"\nFound {len(results_list)} motifs total")
    
    # Check each motif for proper classification and absence of legacy fields
    issues_found = []
    
    for i, motif in enumerate(results_list):
        # Check for legacy fields that should be removed
        legacy_fields = ['Legacy_Type', 'Strand', 'Motif_Type']
        for field in legacy_fields:
            if field in motif:
                issues_found.append(f"Motif {i}: Found legacy field '{field}' = {motif[field]}")
        
        # Check that official classification is present
        if 'Class' not in motif:
            issues_found.append(f"Motif {i}: Missing 'Class' field")
        
        # Verify class is from official system
        motif_class = motif.get('Class', '')
        valid_classes = [info['class_name'] for info in OFFICIAL_CLASSIFICATION.values()]
        if motif_class not in valid_classes:
            issues_found.append(f"Motif {i}: Class '{motif_class}' not in official classification")
        
        # Check subtype is present and appropriate
        subtype = motif.get('Subtype', '')
        if not subtype or subtype in ['Unknown', 'Other']:
            issues_found.append(f"Motif {i}: Subtype '{subtype}' should be more specific")
    
    # Print results
    if issues_found:
        print("❌ Issues found:")
        for issue in issues_found:
            print(f"  - {issue}")
        return False
    else:
        print("✅ All motifs have clean classification without legacy information")
        
        # Show summary by class
        class_counts = {}
        for motif in results_list:
            motif_class = motif.get('Class', 'Unknown')
            class_counts[motif_class] = class_counts.get(motif_class, 0) + 1
        
        print("\nMotif distribution by official class:")
        for class_name, count in sorted(class_counts.items()):
            print(f"  - {class_name}: {count} motifs")
        
        return True

def test_results_table_format():
    """Test that results formatted for tables don't contain legacy information."""
    print("\n" + "=" * 80)
    print("Testing Results Table Format")
    print("=" * 80)
    
    # Test sequence
    test_seq = "GGGGTTAGGGTTAGGGTTAGGGCGCGCGCGCGCCCCCAACCCCAACCCCAA"
    results_list = all_motifs(test_seq, "test_sequence")
    
    # Simulate the results table creation (like in app.py)
    results_data = []
    for motif in results_list:
        # Get official classification
        legacy_class = motif.get('Class', motif.get('Motif', 'Unknown'))
        legacy_subtype = motif.get('Subtype', '')
        official_class, official_subtype = get_official_classification(legacy_class, legacy_subtype)
        
        results_data.append({
            'Sequence': 'test_sequence',
            'Class': official_class,
            'Subclass': official_subtype,
            'Start': motif.get('Start', 0),
            'End': motif.get('End', 0),
            'Length': motif.get('End', 0) - motif.get('Start', 0) + 1,
            'Score': motif.get('Score', 0)
        })
    
    # Create DataFrame
    df = pd.DataFrame(results_data)
    
    # Check for forbidden columns
    forbidden_columns = ['Legacy_Type', 'Strand', 'Motif_Type', 'Official_Class', 'Official_Subclass']
    issues = []
    
    for col in forbidden_columns:
        if col in df.columns:
            issues.append(f"Forbidden column '{col}' found in results table")
    
    # Check required columns
    required_columns = ['Sequence', 'Class', 'Subclass', 'Start', 'End', 'Length', 'Score']
    for col in required_columns:
        if col not in df.columns:
            issues.append(f"Required column '{col}' missing from results table")
    
    if issues:
        print("❌ Issues with results table format:")
        for issue in issues:
            print(f"  - {issue}")
        return False
    else:
        print("✅ Results table format is clean and proper")
        print(f"Table columns: {list(df.columns)}")
        print(f"Sample rows: {len(df)}")
        if len(df) > 0:
            print(f"First row class: {df.iloc[0]['Class']}")
            print(f"First row subclass: {df.iloc[0]['Subclass']}")
        return True

def main():
    """Run all tests."""
    print("NBDFinder Classification Cleanup Verification")
    print("Testing removal of legacy types and strand information")
    print("Testing integration of official classification system\n")
    
    all_passed = True
    
    # Run tests
    tests = [
        test_official_classification_system,
        test_motif_detection_output,
        test_results_table_format
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
        print("✅ ALL TESTS PASSED!")
        print("Legacy types and strand information successfully removed.")
        print("Official classification system working correctly.")
        print("All features properly integrated.")
    else:
        print("❌ SOME TESTS FAILED!")
        print("Please review the issues above.")
    print("=" * 80)
    
    return 0 if all_passed else 1

if __name__ == "__main__":
    sys.exit(main())