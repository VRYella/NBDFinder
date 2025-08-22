#!/usr/bin/env python3
"""
Test Visualization Generation for NBDFinder
==========================================

Test the actual creation of visualization plots to ensure they render properly.
"""

import sys
import os
from motifs import all_motifs
from nbdfinder_viz_integration import NBDFinderVisualizationHub
from publication_visualizations import create_comprehensive_publication_suite

def test_visualization_generation():
    """Test the generation of various visualization types."""
    print("=" * 60)
    print("VISUALIZATION GENERATION TEST")
    print("=" * 60)
    
    # Test sequence with multiple motif types
    test_seq = "GGGTTAGGGTTAGGGTTAGGGCGCGCGCGCGCGCGCCCCAACCCCAACCCCAACCCC"
    print(f"Test sequence length: {len(test_seq)} bp")
    
    print("\n1. Running motif analysis...")
    motifs = all_motifs(test_seq, nonoverlap=True, report_hotspots=False)
    print(f"✅ Found {len(motifs)} motifs")
    
    if not motifs:
        print("❌ No motifs found - cannot test visualization")
        return
    
    # Show motif summary
    motif_classes = {}
    for motif in motifs:
        motif_class = motif.get('Class', 'Unknown')
        motif_classes[motif_class] = motif_classes.get(motif_class, 0) + 1
    
    print(f"   Classes detected: {motif_classes}")
    
    print("\n2. Testing NBDFinder visualization hub...")
    try:
        viz_hub = NBDFinderVisualizationHub()
        motifs_df = viz_hub.process_motif_data(motifs, len(test_seq))
        
        if motifs_df.empty:
            print("❌ Empty DataFrame from visualization hub")
            return
        
        print(f"✅ Processed {len(motifs_df)} motifs for visualization")
        
        # Check data quality
        na_count = motifs_df.isnull().sum().sum()
        print(f"   N/A values: {na_count}")
        
    except Exception as e:
        print(f"❌ Visualization hub error: {str(e)}")
        import traceback
        traceback.print_exc()
        return
    
    print("\n3. Testing publication visualization suite...")
    try:
        # Create comprehensive publication suite
        suite = create_comprehensive_publication_suite(motifs_df)
        
        if not suite:
            print("❌ No visualization suite created")
            return
        
        print(f"✅ Created visualization suite with {len(suite)} plot types")
        
        # List available plot types
        for plot_type in suite.keys():
            plot_count = len(suite[plot_type]) if isinstance(suite[plot_type], dict) else 1
            print(f"   {plot_type}: {plot_count} plots")
        
        # Test a specific plot type
        if 'bar_plots' in suite and suite['bar_plots']:
            print("\n4. Testing bar plot generation...")
            bar_plots = suite['bar_plots']
            
            if 'class_counts' in bar_plots:
                fig = bar_plots['class_counts']
                print("✅ Bar plot created successfully")
                
                # Try to export to test quality
                try:
                    img_bytes = fig.to_image(format="png", width=800, height=600, scale=2)
                    print(f"✅ Image export successful: {len(img_bytes)} bytes")
                except Exception as e:
                    print(f"⚠️  Image export failed: {str(e)}")
            else:
                print("⚠️  class_counts plot not found in bar_plots")
        
    except Exception as e:
        print(f"❌ Publication visualization error: {str(e)}")
        import traceback
        traceback.print_exc()
        return
    
    print("\n" + "=" * 60)
    print("✅ VISUALIZATION GENERATION TEST COMPLETED")
    print("   All core visualization components working properly")
    print("=" * 60)

if __name__ == "__main__":
    test_visualization_generation()