#!/usr/bin/env python3
"""
NBDFinder Bacterial Genome Analysis Results Validation
======================================================

This script validates and summarizes the comprehensive bacterial genome analysis results.
"""

import pandas as pd
import json
import os
from pathlib import Path

def validate_analysis_results():
    """Validate the analysis results and provide summary statistics"""
    
    results_dir = Path("bacterial_genome_analysis_results")
    
    print("🧬 NBDFinder Bacterial Genome Analysis Results Validation")
    print("=" * 70)
    
    # Check directory structure
    required_dirs = ["representative_sequences", "motif_analysis_results", "comparative_plots", "summary_reports"]
    for dir_name in required_dirs:
        dir_path = results_dir / dir_name
        if dir_path.exists():
            file_count = len(list(dir_path.glob("*")))
            print(f"✅ {dir_name}: {file_count} files")
        else:
            print(f"❌ {dir_name}: Missing")
    
    print("\n📊 Analysis Summary:")
    print("-" * 30)
    
    # Load and validate motif results
    motif_files = list((results_dir / "motif_analysis_results").glob("*.json"))
    total_motifs = 0
    genome_data = []
    
    for motif_file in motif_files:
        try:
            with open(motif_file, 'r') as f:
                data = json.load(f)
                
            genome_data.append({
                'name': data['genome_name'],
                'accession': data['accession'],
                'gc_category': data['gc_category'],
                'actual_gc': data['actual_gc'],
                'genome_length': data['genome_length'],
                'total_motifs': data['total_motifs'],
                'motif_density': data['total_motifs'] / data['genome_length'] * 1000
            })
            
            total_motifs += data['total_motifs']
            
        except Exception as e:
            print(f"⚠️ Error reading {motif_file}: {e}")
    
    # Summary statistics
    df = pd.DataFrame(genome_data)
    
    print(f"📈 Total Genomes Analyzed: {len(genome_data)}")
    print(f"🧬 Total Motifs Detected: {total_motifs:,}")
    print(f"📏 Total Sequence Length: {df['genome_length'].sum():,} bp")
    print(f"🎯 Average Motif Density: {df['motif_density'].mean():.2f} per kb")
    
    print("\n🔬 Results by GC Category:")
    print("-" * 35)
    
    for category in ['Low GC', 'Medium GC', 'High GC']:
        cat_data = df[df['gc_category'] == category]
        if not cat_data.empty:
            print(f"{category}:")
            print(f"  - Genomes: {len(cat_data)}")
            print(f"  - GC Range: {cat_data['actual_gc'].min():.1f}% - {cat_data['actual_gc'].max():.1f}%")
            print(f"  - Avg Motif Density: {cat_data['motif_density'].mean():.2f} ± {cat_data['motif_density'].std():.2f} per kb")
            print(f"  - Total Motifs: {cat_data['total_motifs'].sum():,}")
    
    # Check for key files
    print("\n📁 Key Output Files:")
    print("-" * 25)
    
    key_files = [
        "summary_reports/comprehensive_motif_analysis.xlsx",
        "summary_reports/analysis_summary_report.txt", 
        "comparative_plots/comprehensive_motif_analysis.png",
        "README.md"
    ]
    
    for file_path in key_files:
        full_path = results_dir / file_path
        if full_path.exists():
            size_mb = full_path.stat().st_size / (1024 * 1024)
            print(f"✅ {file_path} ({size_mb:.2f} MB)")
        else:
            print(f"❌ {file_path}: Missing")
    
    print("\n🎉 Analysis Results Validation Complete!")
    print("📍 All results are stored in the 'bacterial_genome_analysis_results' directory")
    print("📖 See README.md for detailed documentation and usage instructions")

if __name__ == "__main__":
    validate_analysis_results()