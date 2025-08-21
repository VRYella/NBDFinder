#!/usr/bin/env python3
"""
NBDFinder Complete S. cerevisiae Genome Analysis
===============================================

Final comprehensive analysis script for all S. cerevisiae chromosomes.
Generates all required NBDFinder outputs in the genome_data folder.

Author: NBDFinder Analysis Team
"""

import os
import sys
import time
import pandas as pd
import numpy as np
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp

# Import NBDFinder modules
try:
    from motifs import all_motifs
except ImportError as e:
    print(f"Error importing NBDFinder modules: {e}")
    sys.exit(1)

def analyze_single_chromosome(chr_data):
    """Analyze a single chromosome - designed for multiprocessing"""
    chr_name, sequence = chr_data
    print(f"🧬 Processing {chr_name} ({len(sequence):,} bp)")
    
    # Clean sequence
    sequence = sequence.upper().replace('N', '').replace(' ', '').replace('\n', '')
    
    if len(sequence) < 1000:
        return chr_name, []
    
    # Analyze in chunks for large chromosomes
    chunk_size = 50000
    overlap = 5000
    all_motifs_chr = []
    
    if len(sequence) <= chunk_size:
        motifs = all_motifs(sequence, nonoverlap=False, report_hotspots=False, sequence_name=chr_name)
        all_motifs_chr.extend(motifs)
    else:
        chunk_count = 0
        for i in range(0, len(sequence), chunk_size - overlap):
            chunk_start = i
            chunk_end = min(i + chunk_size, len(sequence))
            chunk_seq = sequence[chunk_start:chunk_end]
            
            # Analyze chunk
            chunk_motifs = all_motifs(chunk_seq, nonoverlap=False, report_hotspots=False, sequence_name=chr_name)
            
            # Adjust coordinates
            for motif in chunk_motifs:
                if isinstance(motif, dict):
                    motif['Start'] = motif.get('Start', 0) + chunk_start
                    motif['End'] = motif.get('End', 0) + chunk_start
                    motif['Sequence Name'] = chr_name
            
            all_motifs_chr.extend(chunk_motifs)
            chunk_count += 1
            
            if chunk_count % 10 == 0:
                print(f"  {chr_name}: Processed {chunk_count} chunks")
    
    print(f"✅ {chr_name}: Found {len(all_motifs_chr)} motifs")
    return chr_name, all_motifs_chr

def main():
    """Main analysis function"""
    print("🧬 NBDFinder Complete S. cerevisiae Genome Analysis")
    print("=" * 60)
    
    # Setup
    output_dir = "example_inputs/genome_data"
    genome_file = "example_inputs/genome_data/GCF_000146045.2_R64_genomic.fna"
    
    if not os.path.exists(genome_file):
        print(f"❌ Genome file not found: {genome_file}")
        return
    
    # Read all chromosomes
    print(f"📖 Reading S. cerevisiae genome...")
    sequences = {}
    for record in SeqIO.parse(genome_file, "fasta"):
        sequences[record.id] = str(record.seq).upper()
    
    print(f"📊 Found {len(sequences)} chromosomes")
    
    # Analyze chromosomes
    start_time = time.time()
    all_results = {}
    
    # Process chromosomes sequentially for reliability
    for i, (chr_name, sequence) in enumerate(sequences.items(), 1):
        print(f"\n🧬 [{i}/{len(sequences)}] Analyzing {chr_name}")
        chr_name, motifs = analyze_single_chromosome((chr_name, sequence))
        all_results[chr_name] = motifs
    
    total_time = time.time() - start_time
    print(f"\n⏱️  Total analysis time: {total_time:.1f} seconds")
    
    # Combine all results
    all_motifs_combined = []
    for chr_motifs in all_results.values():
        all_motifs_combined.extend(chr_motifs)
    
    print(f"🎯 Total motifs across genome: {len(all_motifs_combined)}")
    
    # Generate comprehensive outputs
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Detailed results
    if all_motifs_combined:
        df_motifs = pd.DataFrame(all_motifs_combined)
        motifs_file = os.path.join(output_dir, "S_cerevisiae_NBDFinder_Detailed_Results.csv")
        df_motifs.to_csv(motifs_file, index=False)
        print(f"✅ Detailed results: {motifs_file}")
    
    # 2. Summary statistics
    type_counts = {}
    total_coverage = 0
    for motif in all_motifs_combined:
        if isinstance(motif, dict):
            motif_type = motif.get('Class', 'Unknown')
            type_counts[motif_type] = type_counts.get(motif_type, 0) + 1
            total_coverage += motif.get('Length', 0)
    
    # Calculate scores
    scores = []
    for motif in all_motifs_combined:
        if isinstance(motif, dict) and motif.get('Score'):
            try:
                scores.append(float(motif.get('Score', 0)))
            except:
                pass
    
    summary_stats = {
        'genome_name': 'S_cerevisiae_complete',
        'total_chromosomes': len(sequences),
        'total_motifs_detected': len(all_motifs_combined),
        'unique_motif_types': len(type_counts),
        'total_coverage_bp': total_coverage,
        'mean_score': np.mean(scores) if scores else 0,
        'median_score': np.median(scores) if scores else 0,
        'analysis_time_seconds': total_time
    }
    
    stats_file = os.path.join(output_dir, "S_cerevisiae_NBDFinder_Summary_Statistics.csv")
    pd.DataFrame([summary_stats]).to_csv(stats_file, index=False)
    print(f"✅ Summary statistics: {stats_file}")
    
    # 3. Motif distribution
    distribution_df = pd.DataFrame(list(type_counts.items()), columns=['Motif_Type', 'Count'])
    distribution_df = distribution_df.sort_values('Count', ascending=False)
    dist_file = os.path.join(output_dir, "S_cerevisiae_NBDFinder_Motif_Distribution.csv")
    distribution_df.to_csv(dist_file, index=False)
    print(f"✅ Motif distribution: {dist_file}")
    
    # 4. Positional matrix
    position_data = []
    for motif in all_motifs_combined:
        if isinstance(motif, dict):
            position_data.append({
                'Chromosome': motif.get('Sequence Name', 'Unknown'),
                'Motif_Type': motif.get('Class', 'Unknown'),
                'Position': motif.get('Start', 0),
                'Length': motif.get('Length', 0),
                'Score': motif.get('Score', 0)
            })
    
    pos_file = os.path.join(output_dir, "S_cerevisiae_NBDFinder_Positional_Matrix.csv")
    pd.DataFrame(position_data).to_csv(pos_file, index=False)
    print(f"✅ Positional matrix: {pos_file}")
    
    # 5. Analysis report
    report_file = os.path.join(output_dir, "S_cerevisiae_NBDFinder_Analysis_Report.txt")
    with open(report_file, 'w') as f:
        f.write("="*80 + "\n")
        f.write("NBDFinder 2.0 - Complete S. cerevisiae Genome Analysis Report\n")
        f.write("="*80 + "\n\n")
        
        f.write("ANALYSIS OVERVIEW\n")
        f.write("-" * 40 + "\n")
        f.write(f"Genome: S. cerevisiae (Complete)\n")
        f.write(f"Chromosomes Analyzed: {len(sequences)}\n")
        f.write(f"Analysis Date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"Total Analysis Time: {total_time:.1f} seconds\n\n")
        
        f.write("SUMMARY STATISTICS\n")
        f.write("-" * 40 + "\n")
        f.write(f"Total Motifs Detected: {len(all_motifs_combined):,}\n")
        f.write(f"Unique Motif Types: {len(type_counts)}\n")
        f.write(f"Total Coverage: {total_coverage:,} bp\n\n")
        
        f.write("MOTIF TYPE DISTRIBUTION\n")
        f.write("-" * 40 + "\n")
        for motif_type, count in sorted(type_counts.items(), key=lambda x: x[1], reverse=True):
            percentage = (count / len(all_motifs_combined)) * 100
            f.write(f"{motif_type:25}: {count:6,} ({percentage:5.1f}%)\n")
        
        f.write("\n" + "="*80 + "\n")
        f.write("Analysis completed successfully with NBDFinder 2.0\n")
        f.write("="*80 + "\n")
    
    print(f"✅ Analysis report: {report_file}")
    
    print("\n🎉 Complete S. cerevisiae Analysis Finished!")
    print(f"📋 All results saved to: {output_dir}")
    print(f"📊 {len(all_motifs_combined):,} total Non-B DNA motifs detected")
    print(f"🧬 {len(type_counts)} different motif types found")

if __name__ == "__main__":
    main()