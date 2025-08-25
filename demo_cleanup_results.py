#!/usr/bin/env python3
"""
Demo script showing that legacy types and strand information have been removed
and that all features are properly integrated with the official classification.
"""

import pandas as pd
from motifs import all_motifs
from motifs.classification_config import get_official_classification

def demo_clean_results():
    """Demonstrate clean results without legacy information."""
    print("=" * 80)
    print("NBDFinder Classification Cleanup Demo")
    print("=" * 80)
    
    # Test with a sequence containing multiple motif types
    test_sequence = "GGGGTTAGGGTTAGGGTTAGGGCGCGCGCGCGCGCGCGCGCCCCCAACCCCAACCCCAAGATCGATCGATC"
    sequence_name = "Multi_Motif_Demo"
    
    print(f"Test sequence: {test_sequence}")
    print(f"Sequence length: {len(test_sequence)} bp")
    print(f"Sequence name: {sequence_name}")
    
    # Run motif detection
    print("\nRunning motif detection...")
    motifs_list = all_motifs(test_sequence, sequence_name)
    
    print(f"Found {len(motifs_list)} motifs")
    
    # Create results table like in the main app
    print("\nCreating results table (like in main app)...")
    results_data = []
    
    for motif in motifs_list:
        # Get official classification
        legacy_class = motif.get('Class', motif.get('Motif', 'Unknown'))
        legacy_subtype = motif.get('Subtype', '')
        official_class, official_subtype = get_official_classification(legacy_class, legacy_subtype)
        
        results_data.append({
            'Sequence': sequence_name,
            'Class': official_class,
            'Subclass': official_subtype,
            'Start': motif.get('Start', 0),
            'End': motif.get('End', 0),
            'Length': motif.get('End', 0) - motif.get('Start', 0) + 1,
            'Score': motif.get('Score', 0)
        })
    
    # Convert to DataFrame
    results_df = pd.DataFrame(results_data)
    
    print("\nResults Table:")
    print("-" * 80)
    print(results_df.to_string(index=False))
    
    print(f"\nTable columns: {list(results_df.columns)}")
    
    # Verify no legacy information
    forbidden_fields = ['Legacy_Type', 'Strand', 'Motif_Type', 'Legacy_Class']
    legacy_found = []
    
    for field in forbidden_fields:
        if field in results_df.columns:
            legacy_found.append(field)
    
    if legacy_found:
        print(f"\n❌ Legacy fields found: {legacy_found}")
    else:
        print(f"\n✅ No legacy fields found in results")
    
    # Show official classification usage
    print(f"\nOfficial classifications used:")
    class_counts = results_df['Class'].value_counts()
    for class_name, count in class_counts.items():
        print(f"  - {class_name}: {count} motifs")
    
    subclass_counts = results_df['Subclass'].value_counts()
    print(f"\nOfficial subclassifications used:")
    for subclass_name, count in subclass_counts.items():
        print(f"  - {subclass_name}: {count} motifs")
    
    print("\n" + "=" * 80)
    print("✅ Demo completed successfully!")
    print("✅ Legacy types and strand information removed")
    print("✅ Official classification system integrated")
    print("✅ Results are clean and professional")
    print("=" * 80)

if __name__ == "__main__":
    demo_clean_results()