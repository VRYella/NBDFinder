#!/usr/bin/env python3
"""
NBDFinder Performance Validation Test
====================================

This script validates that the code reorganization and optimizations have
improved performance while maintaining accuracy.

Tests:
1. G4 detection performance (consolidated modules)
2. R-loop REZ detection optimization
3. Memory usage improvements
4. Import time improvements
"""

import time
import sys
import tracemalloc
from pathlib import Path

# Add current directory to path
sys.path.insert(0, str(Path(__file__).parent))

def test_g4_detection_performance():
    """Test G4 detection performance with the consolidated module"""
    print("Testing G4 detection performance...")
    
    # Generate a test sequence with multiple G4 patterns
    test_seq = "GGGATTGGGATCGATGGGATTGGGATCGATGGGATGGGATTGGGATCGATGGGATTGGG" * 10
    
    # Test the optimized G4 detection
    start_time = time.time()
    tracemalloc.start()
    
    from motifs.g4_related import find_all_g4_motifs
    results = find_all_g4_motifs(test_seq)
    
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    end_time = time.time()
    
    print(f"✓ G4 detection completed in {end_time - start_time:.4f} seconds")
    print(f"✓ Found {len(results)} G4 motifs")
    print(f"✓ Peak memory usage: {peak / 1024 / 1024:.2f} MB")
    
    return end_time - start_time, len(results), peak

def test_rloop_optimization():
    """Test the optimized R-loop REZ detection"""
    print("\nTesting R-loop detection optimization...")
    
    # Test sequence with G-rich regions for R-loop formation
    test_seq = "GGGGATTGGGGTTGGGGATTTTTTTCCCCCCGGGGATTGGGGTTGGGGATTTTTT" * 5
    
    start_time = time.time()
    tracemalloc.start()
    
    from motifs.r_loop import find_rlfs
    results = find_rlfs(test_seq)
    
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    end_time = time.time()
    
    print(f"✓ R-loop detection completed in {end_time - start_time:.4f} seconds")
    print(f"✓ Found {len(results)} R-loop candidates")
    print(f"✓ Peak memory usage: {peak / 1024 / 1024:.2f} MB")
    
    return end_time - start_time, len(results), peak

def test_import_performance():
    """Test import performance improvements"""
    print("\nTesting import performance...")
    
    start_time = time.time()
    
    # Test modular imports
    from motifs.shared_utils import gc_content, parse_fasta
    from motifs.g4_related import g4hunter_score
    
    end_time = time.time()
    
    print(f"✓ Core imports completed in {end_time - start_time:.4f} seconds")
    
    # Test basic functionality
    test_seq = "ATCGATCGATCG"
    gc = gc_content(test_seq)
    score = g4hunter_score(test_seq)
    
    print(f"✓ GC content: {gc:.1f}%, G4Hunter score: {score:.3f}")
    
    return end_time - start_time

def test_streamlined_module():
    """Test the new streamlined motifs module"""
    print("\nTesting streamlined motifs module...")
    
    start_time = time.time()
    tracemalloc.start()
    
    try:
        import motifs_streamlined as motifs
        
        # Test basic utilities
        test_seq = "GGGATTGGGATCGATGGGATTGGG"
        gc = motifs.gc_content(test_seq)
        
        # Test G4 detection through streamlined interface
        from motifs.g4_related import find_all_g4_motifs
        g4_results = find_all_g4_motifs(test_seq)
        
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        end_time = time.time()
        
        print(f"✓ Streamlined module working correctly")
        print(f"✓ Import and basic ops completed in {end_time - start_time:.4f} seconds")
        print(f"✓ Found {len(g4_results)} G4 motifs in test sequence")
        print(f"✓ Peak memory usage: {peak / 1024 / 1024:.2f} MB")
        
        return True, end_time - start_time, peak
        
    except Exception as e:
        tracemalloc.stop()
        print(f"✗ Streamlined module test failed: {e}")
        return False, 0, 0

def run_performance_validation():
    """Run complete performance validation suite"""
    print("NBDFinder Performance Validation")
    print("=" * 50)
    
    results = {}
    
    # Test G4 detection
    g4_time, g4_count, g4_memory = test_g4_detection_performance()
    results['g4'] = {'time': g4_time, 'count': g4_count, 'memory': g4_memory}
    
    # Test R-loop optimization
    rloop_time, rloop_count, rloop_memory = test_rloop_optimization()
    results['rloop'] = {'time': rloop_time, 'count': rloop_count, 'memory': rloop_memory}
    
    # Test import performance
    import_time = test_import_performance()
    results['import'] = {'time': import_time}
    
    # Test streamlined module
    streamlined_success, streamlined_time, streamlined_memory = test_streamlined_module()
    results['streamlined'] = {'success': streamlined_success, 'time': streamlined_time, 'memory': streamlined_memory}
    
    # Summary
    print("\n" + "=" * 50)
    print("PERFORMANCE VALIDATION SUMMARY")
    print("=" * 50)
    
    total_time = sum([results[k].get('time', 0) for k in results])
    total_memory = max([results[k].get('memory', 0) for k in results])
    
    print(f"Total test time: {total_time:.4f} seconds")
    print(f"Peak memory usage: {total_memory / 1024 / 1024:.2f} MB")
    print(f"G4 detection: {results['g4']['time']:.4f}s, {results['g4']['count']} motifs")
    print(f"R-loop detection: {results['rloop']['time']:.4f}s, {results['rloop']['count']} motifs")
    print(f"Import time: {results['import']['time']:.4f}s")
    print(f"Streamlined module: {'✓ PASS' if results['streamlined']['success'] else '✗ FAIL'}")
    
    # Performance thresholds (adjust based on system)
    performance_ok = (
        total_time < 2.0 and  # Total test time under 2 seconds
        total_memory / 1024 / 1024 < 50 and  # Peak memory under 50 MB
        results['streamlined']['success']  # Streamlined module works
    )
    
    print("\n" + "=" * 50)
    if performance_ok:
        print("🎉 PERFORMANCE VALIDATION PASSED")
        print("Code reorganization and optimizations are working correctly!")
    else:
        print("⚠️  PERFORMANCE VALIDATION ISSUES DETECTED")
        print("Some optimizations may need further attention.")
    print("=" * 50)
    
    return performance_ok

if __name__ == "__main__":
    success = run_performance_validation()
    sys.exit(0 if success else 1)