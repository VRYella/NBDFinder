#!/usr/bin/env python3
"""
NBDFinder Performance Test - Fast Mode vs Complete Mode
======================================================

Tests the performance improvements achieved through optimization.
"""

import time
import sys
import os
sys.path.append('.')

from motifs import all_motifs

def performance_test():
    """Test complete mode performance (fast mode is now obsolete)"""
    
    print("=" * 80)
    print("NBDFinder Performance Test - Complete Mode Only")
    print("=" * 80)
    
    # Test sequences of different sizes
    test_sequences = {
        "Small (100bp)": "GGGTTAGGGTTAGGGTTAGGGAAACCCAAACCCAAACCCAAACGCGCGCGCGCGAATTTAAATTTAAATTTAAAGAAGAAGAAGAAGAAGAA",
        "Medium (500bp)": "GGGTTAGGGTTAGGGTTAGGG" * 25,
        "Large (2000bp)": "GGGTTAGGGTTAGGGTTAGGG" * 100,
    }
    
    results = {}
    
    for seq_name, sequence in test_sequences.items():
        print(f"\nTesting {seq_name} ({len(sequence)} bp):")
        print("-" * 50)
        
        # Test Complete Mode Only
        start_time = time.time()
        complete_motifs = all_motifs(sequence, report_hotspots=True)
        complete_time = (time.time() - start_time) * 1000  # Convert to ms
        
        print(f"Complete Mode: {len(complete_motifs):3d} motifs in {complete_time:7.2f} ms") 
        print(f"Speed:         {len(sequence)/complete_time*1000:.0f} bp/second")
        
        results[seq_name] = {
            'complete_time': complete_time,
            'complete_motifs': len(complete_motifs),
            'speed_bps': len(sequence)/complete_time*1000
        }
    
    print("\n" + "=" * 80)
    print("PERFORMANCE SUMMARY")
    print("=" * 80)
    
    for seq_name, data in results.items():
        print(f"{seq_name}:")
        print(f"  Time: {data['complete_time']:.2f} ms")
        print(f"  Motifs: {data['complete_motifs']}")
        print(f"  Throughput: {data['speed_bps']:.0f} bp/second")
        print()
    
    # Overall performance metrics
    avg_speed = sum(r['speed_bps'] for r in results.values()) / len(results)
    max_speed = max(r['speed_bps'] for r in results.values())
    
    print(f"Average throughput: {avg_speed:.0f} bp/second")
    print(f"Maximum throughput: {max_speed:.0f} bp/second")
    
    if avg_speed >= 10000:
        print("✅ EXCELLENT PERFORMANCE: >10,000 bp/second")
    elif avg_speed >= 5000:
        print("✅ GOOD PERFORMANCE: >5,000 bp/second")
    else:
        print("⚠️  PERFORMANCE BELOW TARGET: <5,000 bp/second")
    
    print("=" * 80)

if __name__ == "__main__":
    performance_test()