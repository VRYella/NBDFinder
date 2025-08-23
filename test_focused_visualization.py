#!/usr/bin/env python3
"""
Focused Visualization Test for NBDFinder
=======================================

A streamlined test to identify and fix visualization issues quickly.
"""

import sys
from motifs import all_motifs
from nbdfinder_viz_integration import NBDFinderVisualizationHub

def test_simple_visualization():
    """Test visualization with a simple sequence to identify issues quickly."""
    print("=" * 60)
    print("FOCUSED VISUALIZATION TEST")
    print("=" * 60)
    
    # Simple test sequence
    test_seq = "GGGTTAGGGTTAGGGTTAGGG"  # Simple G4 sequence
    print(f"Testing with simple G4 sequence: {test_seq}")
    print(f"Length: {len(test_seq)} bp")
    
    # Run motif analysis with reduced complexity
    print("\n1. Running motif analysis...")
    try:
        motifs = all_motifs(test_seq, nonoverlap=True, report_hotspots=False)
        print(f"✅ Found {len(motifs)} motifs")
        
        if motifs:
            # Show first few motifs for inspection
            print("\n2. Sample motifs:")
            for i, motif in enumerate(motifs[:3]):
                print(f"   Motif {i+1}: {motif}")
            
            # Test visualization processing
            print("\n3. Testing visualization processing...")
            viz_hub = NBDFinderVisualizationHub()
            
            # Process motif data
            motifs_df = viz_hub.process_motif_data(motifs, len(test_seq))
            
            if not motifs_df.empty:
                print(f"✅ Visualization DataFrame created: {len(motifs_df)} rows")
                print(f"   Columns: {list(motifs_df.columns)}")
                
                # Check for N/A values
                na_counts = motifs_df.isnull().sum()
                total_na = na_counts.sum()
                
                if total_na == 0:
                    print("✅ No N/A values found in processed data")
                else:
                    print(f"⚠️  Found {total_na} N/A values:")
                    for col, count in na_counts.items():
                        if count > 0:
                            print(f"      {col}: {count} N/A values")
                
                print("\n4. Sample processed data:")
                print(motifs_df.head())
                
            else:
                print("❌ Empty visualization DataFrame")
                
        else:
            print("⚠️  No motifs detected")
            
    except Exception as e:
        print(f"❌ Error: {str(e)}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test_simple_visualization()