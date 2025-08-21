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
    """Test performance improvements between fast and complete modes"""
    
    print("=" * 80)
    print("NBDFinder Performance Test - Fast vs Complete Mode")
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
        
        # Test Fast Mode
        start_time = time.time()
        fast_motifs = all_motifs(sequence, fast_mode=True)
        fast_time = (time.time() - start_time) * 1000  # Convert to ms
        
        # Test Complete Mode  
        start_time = time.time()
        complete_motifs = all_motifs(sequence, fast_mode=False, report_hotspots=True)
        complete_time = (time.time() - start_time) * 1000  # Convert to ms
        
        # Calculate speedup
        speedup = complete_time / fast_time if fast_time > 0 else 1.0
        
        print(f"Fast Mode:     {len(fast_motifs):3d} motifs in {fast_time:7.2f} ms")
        print(f"Complete Mode: {len(complete_motifs):3d} motifs in {complete_time:7.2f} ms") 
        print(f"Speedup:       {speedup:.1f}x faster")
        print(f"Speed:         {len(sequence)/fast_time*1000:.0f} bp/second (fast mode)")
        
        results[seq_name] = {
            'fast_time': fast_time,
            'complete_time': complete_time,
            'fast_motifs': len(fast_motifs),
            'complete_motifs': len(complete_motifs),
            'speedup': speedup,
            'speed_bps': len(sequence)/fast_time*1000
        }
    
    print("\n" + "=" * 80)
    print("PERFORMANCE SUMMARY")
    print("=" * 80)
    
    for seq_name, data in results.items():
        print(f"{seq_name}:")
        print(f"  Speedup: {data['speedup']:.1f}x")
        print(f"  Throughput: {data['speed_bps']:.0f} bp/second")
        print(f"  Motif detection efficiency: {data['fast_motifs']/data['complete_motifs']*100:.1f}% of complete analysis")
        print()
    
    # Overall performance metrics
    avg_speedup = sum(r['speedup'] for r in results.values()) / len(results)
    max_speed = max(r['speed_bps'] for r in results.values())
    
    print(f"Average speedup: {avg_speedup:.1f}x")
    print(f"Maximum throughput: {max_speed:.0f} bp/second")
    
    if avg_speedup >= 10:
        print("✅ PERFORMANCE TARGET ACHIEVED: >10x speedup")
    elif avg_speedup >= 5:
        print("✅ GOOD PERFORMANCE: >5x speedup")
    else:
        print("⚠️  PERFORMANCE BELOW TARGET: <5x speedup")
    
    print("=" * 80)

if __name__ == "__main__":
    performance_test()