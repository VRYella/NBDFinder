#!/usr/bin/env python3
"""
Test script for Disease Expansion Analysis - runs on first 5 sequences only
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import time
from datetime import datetime
import json
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Add motifs directory to path for NBDFinder imports
sys.path.append(os.path.join(os.path.dirname(__file__), 'motifs'))
from motifs.shared_utils import all_motifs, gc_content, parse_fasta_multi
from disease_motifs import AdvancedDiseaseDetector

def test_analysis():
    """Test analysis on first 5 sequences"""
    print("🧬 Testing Disease Expansion Analysis on 5 sequences")
    print("=" * 50)
    
    # Load sequences
    with open('repeat_expansion_loci_annotated.fa', 'r') as f:
        content = f.read()
    
    sequences, names = parse_fasta_multi(content)
    print(f"Loaded {len(sequences)} sequences, testing first 5...")
    
    # Test first 5 sequences
    test_sequences = sequences[:5]
    test_names = names[:5]
    
    disease_detector = AdvancedDiseaseDetector()
    results = []
    
    for i, (sequence, name) in enumerate(zip(test_sequences, test_names)):
        print(f"\n🔬 Analyzing sequence {i+1}: {name[:60]}...")
        
        start_time = time.time()
        
        # Basic statistics
        seq_length = len(sequence)
        gc_percent = gc_content(sequence)
        
        # Extract gene info
        parts = name.split('|')
        gene_symbol = parts[1] if len(parts) > 1 else "Unknown"
        
        print(f"   📏 Length: {seq_length} bp, GC: {gc_percent:.1f}%")
        
        # Motif analysis
        motifs = all_motifs(sequence, sequence_name=gene_symbol)
        
        # Disease analysis
        disease_motifs = disease_detector.detect_pathogenic_repeats(sequence, gene_symbol)
        
        analysis_time = time.time() - start_time
        
        result = {
            'gene_symbol': gene_symbol,
            'sequence_length': seq_length,
            'gc_content': gc_percent,
            'total_motifs': len(motifs),
            'disease_motifs': len(disease_motifs),
            'analysis_time': analysis_time
        }
        
        results.append(result)
        
        print(f"   ✅ Found {len(motifs)} motifs, {len(disease_motifs)} disease-specific")
        print(f"   ⏱️  Completed in {analysis_time:.2f} seconds")
    
    # Summary
    print(f"\n📊 Test Summary:")
    print(f"   🧬 Sequences analyzed: {len(results)}")
    print(f"   🔍 Total motifs: {sum(r['total_motifs'] for r in results)}")
    print(f"   🦠 Disease motifs: {sum(r['disease_motifs'] for r in results)}")
    print(f"   ⏱️  Total time: {sum(r['analysis_time'] for r in results):.1f} seconds")
    print(f"   📈 Average motifs per gene: {np.mean([r['total_motifs'] for r in results]):.1f}")
    
    return results

if __name__ == "__main__":
    test_analysis()