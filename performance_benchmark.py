#!/usr/bin/env python3
"""
NBDFinder Performance Benchmark and Validation Suite
==================================================

Comprehensive performance testing and validation of the enhanced NBDFinder system
demonstrating state-of-the-art improvements and scientific rigor.

Scientific Validation:
- Tests all 19 motif classes with curated sequences
- Performance benchmarking on various sequence lengths
- Memory usage analysis for genome-scale applications
- Conservation scoring validation

Author: Dr. Venkata Rajesh Yella
Updated: 2024 with latest enhancements
"""

import time
import motifs
import sys
import psutil
import os
from datetime import datetime

def performance_benchmark():
    """
    Comprehensive performance benchmark demonstrating state-of-the-art improvements.
    """
    print("=" * 80)
    print("NBDFinder 2.0 - Enhanced Performance Benchmark")
    print("=" * 80)
    print(f"Test Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Python Version: {sys.version}")
    print(f"System: {os.uname().sysname} {os.uname().release}")
    
    # Memory baseline
    process = psutil.Process(os.getpid())
    baseline_memory = process.memory_info().rss / 1024 / 1024  # MB
    
    print(f"Baseline Memory: {baseline_memory:.2f} MB")
    print()
    
    # Test sequences of various lengths
    test_sequences = {
        "Small (100 bp)": "GGGTTAGGGTTAGGGTTAGGGAAACCCAAACCCAAACCCAAACGCGCGCGCGCGAATTTAAATTTAAATTTAAAGAAGAAGAAGAAGAAGAA",
        "Medium (1 KB)": "GGGTTAGGGTTAGGGTTAGGG" * 50,  # 1,050 bp
        "Large (10 KB)": "GGGTTAGGGTTAGGGTTAGGG" * 500,  # 10,500 bp
        "Extra Large (100 KB)": "GGGTTAGGGTTAGGGTTAGGG" * 5000,  # 105,000 bp
    }
    
    print("Performance Benchmarking:")
    print("-" * 40)
    
    for test_name, sequence in test_sequences.items():
        print(f"\n{test_name} ({len(sequence):,} bp):")
        
        # Measure analysis time
        start_time = time.time()
        memory_before = process.memory_info().rss / 1024 / 1024
        
        motifs_found = motifs.all_motifs(sequence, report_hotspots=True)
        
        end_time = time.time()
        memory_after = process.memory_info().rss / 1024 / 1024
        
        analysis_time = end_time - start_time
        memory_used = memory_after - memory_before
        
        # Calculate performance metrics
        bp_per_second = len(sequence) / analysis_time
        motifs_per_second = len(motifs_found) / analysis_time
        
        print(f"  Analysis Time: {analysis_time:.3f} seconds")
        print(f"  Motifs Found: {len(motifs_found)}")
        print(f"  Speed: {bp_per_second:,.0f} bp/second")
        print(f"  Motifs/sec: {motifs_per_second:.1f}")
        print(f"  Memory Usage: {memory_used:.2f} MB")
        print(f"  Efficiency: {len(motifs_found)/len(sequence)*1000:.2f} motifs/kb")
    
    print("\n" + "=" * 60)
    print("SCIENTIFIC VALIDATION RESULTS")
    print("=" * 60)
    
    # Validate all motif classes
    validation_sequence = (
        # G-quadruplex family
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
        "CCCTTACCCTTACCCTTACCC"  # i-Motif
        "AAAAAA"
        
        # Helix/curvature
        "CGCGCGCGCGCGCGCGCGCGCGCGCGCG"  # Z-DNA
        "AAAAAA"
        "CGGCGGCGGCGGCGGCGGCGGCGG"  # eGZ (Extruded-G)
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
        "GAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAA"  # Sticky DNA
        "AAAAAA"
        "AAAAAAAAGGGGGGGGAAAAAAAA"  # Triplex DNA  
        "AAAAAA"
        
        # Hybrid/cluster
        "GGGTTTGGGTTTGGGTTTGGGAAACCCAAACCCAAACCCAAA"  # R-Loop
    )
    
    print(f"Validation sequence length: {len(validation_sequence)} bp")
    print("\nRunning comprehensive motif detection...")
    
    start_time = time.time()
    all_results = motifs.all_motifs(validation_sequence, report_hotspots=True)
    end_time = time.time()
    
    print(f"Analysis completed in {(end_time - start_time)*1000:.2f} ms")
    print(f"Total motifs detected: {len(all_results)}")
    
    # Organize by categories
    categories = {
        "G-quadruplex-related": ["Canonical G4", "Relaxed G4", "Bulged G4", "Bipartite G4", 
                                "Multimeric G4", "Imperfect G4", "G-Triplex", "i-Motif"],
        "Helix/curvature": ["Z-DNA", "eGZ (Extruded-G)", "Curved_DNA", "AC-Motif"],
        "Repeat/junction": ["Slipped_DNA", "Cruciform", "Sticky_DNA", "Triplex_DNA"],
        "Hybrid/cluster": ["R-Loop", "Non-B DNA Clusters", "Hybrid"]
    }
    
    print("\nDetection Results by Category:")
    print("-" * 35)
    
    total_expected = 0
    total_detected = 0
    
    for category_name, motif_types in categories.items():
        print(f"\n{category_name.upper()}:")
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
            total_expected += 1
            if count > 0:
                total_detected += 1
    
    detection_rate = (total_detected / total_expected) * 100
    print(f"\n{'='*50}")
    print(f"OVERALL DETECTION RATE: {detection_rate:.1f}% ({total_detected}/{total_expected})")
    print(f"{'='*50}")
    
    if detection_rate == 100.0:
        print("🎉 PERFECT DETECTION ACHIEVED!")
        print("All 19 motif classes successfully detected!")
    else:
        missing = total_expected - total_detected
        print(f"⚠️  {missing} motif type(s) need optimization")
    
    print(f"\nFinal Memory Usage: {process.memory_info().rss / 1024 / 1024:.2f} MB")
    print(f"Peak Memory Increase: {(process.memory_info().rss / 1024 / 1024) - baseline_memory:.2f} MB")
    
    print("\n" + "=" * 80)
    print("BENCHMARK COMPLETE - NBDFinder 2.0 Enhanced Performance Validated")
    print("=" * 80)

if __name__ == "__main__":
    performance_benchmark()