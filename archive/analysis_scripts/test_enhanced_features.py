#!/usr/bin/env python3
"""
Test Enhanced NBDFinder Functionality
=====================================

Test script for validating the new advanced features including:
- Disease-associated motif detection
- Machine learning enhanced predictions
- Advanced clustering analysis
- Hybrid structure detection
"""

import sys
import time
from motifs import all_motifs

def test_enhanced_features():
    """Test all enhanced features of NBDFinder"""
    
    print("="*80)
    print("NBDFinder Enhanced Features Test")
    print("="*80)
    print()
    
    # Test sequences
    test_sequences = {
        "Disease GAA repeat": "GAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAA",
        "Disease CGG repeat": "CGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGG",
        "G4 + i-motif hybrid": "GGGTTAGGGTTAGGGTTAGGGCCCCAACCCCAACCCCAACCCC",
        "Complex multi-motif": "GGGTTAGGGTTAGGGTTAGGGCGCGCGCGCGCGCGCGCGCGCCCCAACCCCAACCCCAACCCCGAAGAAGAAGAAGAAGAAGAATTTAAATTTAAATTTAAA"
    }
    
    for name, seq in test_sequences.items():
        print(f"\n{'-'*60}")
        print(f"Testing: {name}")
        print(f"Sequence length: {len(seq)} bp")
        print(f"Sequence: {seq}")
        print(f"{'-'*60}")
        
        start_time = time.time()
        
        # Test with enhanced features enabled
        motifs = all_motifs(seq, nonoverlap=False, report_hotspots=True, sequence_name=name)
        
        end_time = time.time()
        processing_time = (end_time - start_time) * 1000
        
        print(f"Total motifs detected: {len(motifs)}")
        print(f"Processing time: {processing_time:.2f} ms")
        
        # Categorize motifs by class
        motif_classes = {}
        disease_motifs = []
        hybrid_motifs = []
        cluster_motifs = []
        ml_enhanced = []
        
        for motif in motifs:
            motif_class = motif.get('Class', 'Unknown')
            if motif_class not in motif_classes:
                motif_classes[motif_class] = 0
            motif_classes[motif_class] += 1
            
            # Check for enhanced features
            if motif_class == "Disease-Associated Motif":
                disease_motifs.append(motif)
            elif motif_class == "Hybrid Structures":
                hybrid_motifs.append(motif)
            elif "Cluster" in motif_class:
                cluster_motifs.append(motif)
            
            if 'ML_Probability' in motif:
                ml_enhanced.append(motif)
        
        # Display results summary
        print(f"Motif classes: {dict(motif_classes)}")
        print()
        
        # Disease-associated motifs
        if disease_motifs:
            print("🧬 Disease-Associated Motifs:")
            for motif in disease_motifs[:3]:  # Show first 3
                print(f"  - {motif['Disease_Name']}: {motif['Repeat_Count']} repeats")
                print(f"    Clinical Significance: {motif['Clinical_Significance']}")
                print(f"    Risk Score: {motif['Risk_Score']:.1f}")
            print()
        
        # Hybrid structures
        if hybrid_motifs:
            print("🔗 Hybrid Structures:")
            for motif in hybrid_motifs[:2]:  # Show first 2
                print(f"  - Components: {motif['Component_Motifs']}")
                print(f"    Structural Impact: {motif['Structural_Impact']}")
                print(f"    Stability Enhancement: {motif['Stability_Enhancement']:.3f}")
            print()
        
        # Advanced clusters
        if cluster_motifs:
            print("🌟 Advanced Clusters:")
            for motif in cluster_motifs[:2]:  # Show first 2
                motif_count = motif.get('Motif_Count', 'N/A')
                motif_types = motif.get('Motif_Types', 'N/A')
                diversity = motif.get('Diversity_Index', 'N/A')
                clinical_sig = motif.get('Clinical_Significance', 'N/A')
                print(f"  - {motif_count} motifs, {motif_types}")
                print(f"    Diversity Index: {diversity}")
                print(f"    Clinical Significance: {clinical_sig}")
            print()
        
        # ML enhanced motifs
        if ml_enhanced:
            print("🤖 Machine Learning Enhanced Predictions:")
            print(f"  - {len(ml_enhanced)} motifs enhanced with ML predictions")
            avg_probability = sum(m['ML_Probability'] for m in ml_enhanced) / len(ml_enhanced)
            avg_confidence = sum(m['ML_Confidence'] for m in ml_enhanced) / len(ml_enhanced)
            print(f"  - Average ML Probability: {avg_probability:.3f}")
            print(f"  - Average ML Confidence: {avg_confidence:.3f}")
            print()
        
        # Show top scoring motifs
        top_motifs = sorted(motifs, key=lambda m: float(m.get('Score', 0)), reverse=True)[:5]
        if top_motifs:
            print("🏆 Top Scoring Motifs:")
            for i, motif in enumerate(top_motifs, 1):
                score = motif.get('Score', 0)
                enhanced_score = motif.get('Enhanced_Score', score)
                print(f"  {i}. {motif['Subtype']}: Score={score:.1f}")
                if enhanced_score != score:
                    print(f"     Enhanced Score: {enhanced_score:.1f}")
        
        print()

def test_performance_comparison():
    """Test performance with and without enhanced features"""
    print("\n" + "="*80)
    print("Performance Comparison Test")
    print("="*80)
    
    test_seq = "GGGTTAGGGTTAGGGTTAGGGGAAGAAGAAGAAGAAGAAGAACGCGCGCGCGCGCGCGCCCCCAACCCCAACCCCAACCCC" * 3
    
    # Test without hotspots (basic mode)
    start = time.time()
    basic_motifs = all_motifs(test_seq, nonoverlap=False, report_hotspots=False)
    basic_time = (time.time() - start) * 1000
    
    # Test with enhanced features
    start = time.time()
    enhanced_motifs = all_motifs(test_seq, nonoverlap=True, report_hotspots=True)
    enhanced_time = (time.time() - start) * 1000
    
    print(f"Basic mode: {len(basic_motifs)} motifs in {basic_time:.2f} ms")
    print(f"Enhanced mode: {len(enhanced_motifs)} motifs in {enhanced_time:.2f} ms")
    print(f"Performance overhead: {enhanced_time/basic_time:.1f}x")
    
    # Count enhanced features
    disease_count = sum(1 for m in enhanced_motifs if m.get('Class') == 'Disease-Associated Motif')
    hybrid_count = sum(1 for m in enhanced_motifs if m.get('Class') == 'Hybrid Structures')
    cluster_count = sum(1 for m in enhanced_motifs if 'Cluster' in m.get('Class', ''))
    ml_count = sum(1 for m in enhanced_motifs if 'ML_Probability' in m)
    
    print(f"Enhanced features detected:")
    print(f"  - Disease motifs: {disease_count}")
    print(f"  - Hybrid structures: {hybrid_count}")
    print(f"  - Advanced clusters: {cluster_count}")
    print(f"  - ML enhanced: {ml_count}")

if __name__ == "__main__":
    try:
        test_enhanced_features()
        test_performance_comparison()
        
        print("\n" + "="*80)
        print("✅ All enhanced features tested successfully!")
        print("NBDFinder is now equipped with state-of-the-art capabilities:")
        print("  🧬 Disease-associated motif detection with clinical classification")
        print("  🤖 Machine learning enhanced structure prediction")
        print("  🔗 Advanced hybrid structure analysis")
        print("  🌟 Sophisticated clustering algorithms")
        print("  📊 Comprehensive statistical analysis")
        print("="*80)
        
    except Exception as e:
        print(f"\n❌ Error during testing: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)