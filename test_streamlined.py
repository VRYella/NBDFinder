#!/usr/bin/env python3
"""
Streamlined NBDFinder Test Suite
==============================

Tests for the core functionality of the streamlined NBDFinder toolkit.
"""

import sys
import motifs

def test_basic_functionality():
    """Test basic motif detection functionality."""
    print("🧪 Testing basic functionality...")
    
    # Test sequences
    g4_seq = "GGGTTGGGTGGGTTGGG"  # Known G4
    z_dna_seq = "CGCGCGCGCGCGCGCG"  # Z-DNA sequence
    i_motif_seq = "CCCAACCCAACCCAACCC"  # i-motif sequence
    
    # Test G4 detection
    g4_results = motifs.find_gquadruplex(g4_seq)
    assert len(g4_results) == 1, f"Expected 1 G4, got {len(g4_results)}"
    print("  ✅ G4 detection works")
    
    # Test unified detection
    all_results = motifs.all_motifs(g4_seq, sequence_name="Test")
    assert len(all_results) >= 1, f"Expected ≥1 motifs, got {len(all_results)}"
    print("  ✅ Unified detection works")
    
    # Test field count (should be exactly 13 essential fields)
    if all_results:
        field_count = len(all_results[0])
        assert field_count == 13, f"Expected 13 fields, got {field_count}"
        print("  ✅ Output has exactly 13 essential fields")
    
    # Test required fields
    required_fields = {
        "Sequence Name", "Class", "Subtype", "Start", "End", "Length", 
        "Sequence", "Score", "Conservation_Score", "Conservation_P_Value", 
        "Conservation_Significance", "Arms/Repeat Unit/Copies", "Spacer"
    }
    
    if all_results:
        result_fields = set(all_results[0].keys())
        missing_fields = required_fields - result_fields
        extra_fields = result_fields - required_fields
        
        assert not missing_fields, f"Missing required fields: {missing_fields}"
        assert not extra_fields, f"Extra fields found: {extra_fields}"
        print("  ✅ All required fields present, no extra fields")

def test_g4_definitions():
    """Test exact G4 definitions as specified."""
    print("🧪 Testing G4 definitions...")
    
    # Test canonical G4
    canonical_seq = "GGGTTGGGTGGGTTGGG"
    canonical_results = motifs.find_gquadruplex(canonical_seq)
    assert len(canonical_results) >= 1, "Canonical G4 not detected"
    assert canonical_results[0]["Subtype"] == "Canonical"
    print("  ✅ Canonical G4 detection works")
    
    # Test that priorities work (no overlapping)
    test_seq = "GGGTTGGGTGGGTTTTTGGGTTGGG"  # Should detect canonical first
    canonical_results = motifs.find_gquadruplex(test_seq)
    relaxed_results = motifs.find_relaxed_gquadruplex(test_seq)
    
    # Check that relaxed doesn't overlap with canonical significantly
    if canonical_results and relaxed_results:
        canonical_pos = set(range(canonical_results[0]["Start"], canonical_results[0]["End"]))
        relaxed_pos = set(range(relaxed_results[0]["Start"], relaxed_results[0]["End"]))
        overlap = len(canonical_pos.intersection(relaxed_pos)) / len(relaxed_pos)
        assert overlap <= 0.5, f"Too much overlap between canonical and relaxed: {overlap}"
    
    print("  ✅ G4 priority system works")

def test_performance():
    """Test performance with pre-compiled patterns."""
    print("🧪 Testing performance...")
    
    import time
    
    # Test sequence
    long_seq = "GGGTTGGGTGGGTTGGG" * 100  # Repeat pattern
    
    start_time = time.time()
    results = motifs.all_motifs(long_seq)
    end_time = time.time()
    
    processing_time = end_time - start_time
    sequence_length = len(long_seq)
    
    print(f"  ✅ Processed {sequence_length} bp in {processing_time:.4f}s")
    print(f"  ✅ Rate: {sequence_length/processing_time:.0f} bp/second")
    
    # Should be much faster than 1 second for this small test
    assert processing_time < 1.0, f"Processing too slow: {processing_time}s"

def test_backward_compatibility():
    """Test that existing function names work."""
    print("🧪 Testing backward compatibility...")
    
    test_seq = "GGGTTGGGTGGGTTGGG"
    
    # Test all the preserved function names
    functions_to_test = [
        motifs.find_gquadruplex,
        motifs.find_relaxed_gquadruplex,
        motifs.find_bulged_gquadruplex,
        motifs.find_imperfect_gquadruplex,
        motifs.find_multimeric_gquadruplex,
        motifs.find_bipartite_gquadruplex,
        motifs.find_zdna,
        motifs.find_cruciform,
        motifs.find_gtriplex,
        motifs.find_imotif,
        motifs.find_slipped_dna
    ]
    
    for func in functions_to_test:
        try:
            result = func(test_seq)
            assert isinstance(result, list), f"{func.__name__} should return a list"
            print(f"  ✅ {func.__name__} works")
        except Exception as e:
            print(f"  ❌ {func.__name__} failed: {e}")
            raise

def test_utility_functions():
    """Test utility functions."""
    print("🧪 Testing utility functions...")
    
    # Test G4Hunter score
    score = motifs.g4hunter_score("GGGTTGGGTGGGTTGGG")
    assert isinstance(score, float), "G4Hunter should return float"
    assert score > 0, "G4Hunter score should be positive for G-rich sequence"
    print("  ✅ G4Hunter scoring works")
    
    # Test conservation scoring
    conservation = motifs.calculate_conservation_score("GGGTTGGGTGGGTTGGG", "G4")
    assert isinstance(conservation, dict), "Conservation should return dict"
    required_keys = {"enrichment_score", "p_value", "significance"}
    assert required_keys.issubset(conservation.keys()), "Missing conservation keys"
    print("  ✅ Conservation scoring works")
    
    # Test sequence utilities
    assert motifs.gc_content("GGCC") == 100.0, "GC content calculation error"
    assert motifs.reverse_complement("ATGC") == "GCAT", "Reverse complement error"
    print("  ✅ Sequence utilities work")

def run_all_tests():
    """Run all tests."""
    print("🚀 Starting Streamlined NBDFinder Test Suite")
    print("=" * 60)
    
    try:
        test_basic_functionality()
        test_g4_definitions()
        test_performance()
        test_backward_compatibility()
        test_utility_functions()
        
        print("=" * 60)
        print("✅ All tests passed! Streamlined NBDFinder is working correctly.")
        return True
        
    except Exception as e:
        print("=" * 60)
        print(f"❌ Test failed: {e}")
        return False

if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)