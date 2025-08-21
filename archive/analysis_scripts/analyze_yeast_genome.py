#!/usr/bin/env python3
"""
NBDFinder S. cerevisiae Genome Analysis Script
=============================================

Analyzes the S. cerevisiae genome for Non-B DNA structures and generates
comprehensive results in the genomic data folder.

Author: NBDFinder Analysis Script
Updated: 2024
"""

import os
import sys
import time
import pandas as pd
import numpy as np
from pathlib import Path
from Bio import SeqIO

# Import NBDFinder modules
try:
    from motifs import all_motifs, parse_fasta
    from disease_motifs import AdvancedDiseaseDetector
    from advanced_clustering import AdvancedClusterAnalyzer
except ImportError as e:
    print(f"Error importing NBDFinder modules: {e}")
    print("Please ensure NBDFinder is properly installed.")
    sys.exit(1)

class YeastGenomeAnalyzer:
    """Analyzer for S. cerevisiae genome NBD structures"""
    
    def __init__(self, output_dir="example_inputs/genome_data"):
        """Initialize the analyzer"""
        self.version = "2.0.0"
        self.output_dir = output_dir
        self.all_motif_types = [
            "Bipartite G4", "Bulged G4", "Canonical G4", "Imperfect G4", 
            "Multimeric G4", "Relaxed G4", "G-Triplex", "AC-Motif", "i-Motif",
            "Curved DNA", "Z-DNA", "eGZ (Extruded-G)", "Cruciform", "R-Loop",
            "Slipped DNA", "Sticky DNA", "Triplex DNA", "Hybrid", "Non-B DNA Clusters"
        ]
        
        # Initialize advanced analyzers
        try:
            self.disease_detector = AdvancedDiseaseDetector()
            self.cluster_analyzer = AdvancedClusterAnalyzer()
        except Exception as e:
            print(f"Warning: Advanced analyzers not available: {e}")
            self.disease_detector = None
            self.cluster_analyzer = None
        
        print(f"🧬 NBDFinder {self.version} Yeast Genome Analyzer Initialized")
        print(f"📁 Output directory: {self.output_dir}")
        print(f"✅ Ready to analyze all {len(self.all_motif_types)} motif types")
    
    def read_fasta_file(self, fasta_path):
        """Read and parse FASTA file using BioPython"""
        sequences = {}
        try:
            for record in SeqIO.parse(fasta_path, "fasta"):
                sequences[record.id] = str(record.seq).upper()
            print(f"📖 Successfully read {len(sequences)} sequences from {fasta_path}")
            return sequences
        except Exception as e:
            print(f"❌ Error reading FASTA file: {e}")
            return None
    
    def analyze_sequence_chunk(self, sequence, sequence_name, chunk_size=50000):
        """Analyze a large sequence in chunks to manage memory"""
        print(f"\n🔬 Analyzing: {sequence_name}")
        print(f"📏 Length: {len(sequence):,} bp")
        
        if len(sequence) <= chunk_size:
            # Analyze the whole sequence if it's small enough
            return self._analyze_single_sequence(sequence, sequence_name)
        
        # For large sequences, analyze in overlapping chunks
        print(f"🧩 Processing in chunks of {chunk_size:,} bp...")
        all_motifs = []
        chunk_count = 0
        overlap = 5000  # Overlap to catch motifs spanning chunk boundaries
        
        for i in range(0, len(sequence), chunk_size - overlap):
            chunk_start = i
            chunk_end = min(i + chunk_size, len(sequence))
            chunk_seq = sequence[chunk_start:chunk_end]
            chunk_name = f"{sequence_name}_chunk_{chunk_count}"
            
            print(f"  Processing chunk {chunk_count + 1}: {chunk_start:,}-{chunk_end:,}")
            
            # Analyze this chunk
            chunk_motifs = self._analyze_single_sequence(chunk_seq, chunk_name, verbose=False)
            
            # Adjust coordinates to global sequence coordinates
            for motif in chunk_motifs:
                motif['Start'] += chunk_start
                motif['End'] += chunk_start
                motif['Sequence Name'] = sequence_name  # Set back to original name
            
            all_motifs.extend(chunk_motifs)
            chunk_count += 1
            
            # Progress update
            if chunk_count % 10 == 0:
                print(f"    ✅ Completed {chunk_count} chunks, found {len(all_motifs)} motifs so far")
        
        print(f"🎯 Total motifs found: {len(all_motifs)}")
        return all_motifs
    
    def _analyze_single_sequence(self, sequence, sequence_name, verbose=True):
        """Analyze a single sequence for Non-B DNA structures"""
        start_time = time.time()
        
        # Clean sequence
        sequence = sequence.upper().replace('N', '').replace(' ', '').replace('\n', '')
        
        if len(sequence) < 10:
            print(f"⚠️ Sequence too short: {len(sequence)} bp")
            return []
        
        # Find all motifs
        motifs = all_motifs(sequence, nonoverlap=False, report_hotspots=True, sequence_name=sequence_name)
        
        # Clean subtypes
        for motif in motifs:
            if 'Subtype' in motif and motif['Subtype']:
                motif['Subtype'] = str(motif['Subtype']).replace('_', ' ').strip()
        
        analysis_time = time.time() - start_time
        
        if verbose:
            print(f"⏱️  Analysis completed in {analysis_time*1000:.1f} ms")
            print(f"🎯 Detected {len(motifs)} total motifs")
        
        return motifs
    
    def generate_comprehensive_report(self, all_results, genome_name="S_cerevisiae"):
        """Generate comprehensive analysis report"""
        os.makedirs(self.output_dir, exist_ok=True)
        
        print(f"\n📊 Generating comprehensive report in: {self.output_dir}")
        
        # Combine all motifs from all sequences
        all_motifs_combined = []
        sequence_stats = {}
        
        for seq_name, motifs in all_results.items():
            all_motifs_combined.extend(motifs)
            sequence_stats[seq_name] = {
                'motif_count': len(motifs),
                'sequence_name': seq_name
            }
        
        print(f"📈 Total motifs across all sequences: {len(all_motifs_combined)}")
        
        # 1. Save detailed motif results
        if all_motifs_combined:
            df_motifs = pd.DataFrame(all_motifs_combined)
            motifs_file = os.path.join(self.output_dir, f"{genome_name}_NBDFinder_Detailed_Results.csv")
            df_motifs.to_csv(motifs_file, index=False)
            print(f"✅ Detailed motifs saved: {motifs_file}")
            print(f"   📊 {len(all_motifs_combined)} motifs detected")
        
        # 2. Generate summary statistics
        stats_summary = self._calculate_genome_stats(all_motifs_combined, sequence_stats, genome_name)
        stats_file = os.path.join(self.output_dir, f"{genome_name}_NBDFinder_Summary_Statistics.csv")
        pd.DataFrame([stats_summary]).to_csv(stats_file, index=False)
        print(f"✅ Summary statistics saved: {stats_file}")
        
        # 3. Generate motif type distribution
        if all_motifs_combined:
            type_counts = {}
            for motif in all_motifs_combined:
                motif_type = motif.get('Class', 'Unknown')
                type_counts[motif_type] = type_counts.get(motif_type, 0) + 1
            
            distribution_df = pd.DataFrame(list(type_counts.items()), columns=['Motif_Type', 'Count'])
            distribution_df = distribution_df.sort_values('Count', ascending=False)
            dist_file = os.path.join(self.output_dir, f"{genome_name}_NBDFinder_Motif_Distribution.csv")
            distribution_df.to_csv(dist_file, index=False)
            print(f"✅ Motif distribution saved: {dist_file}")
        
        # 4. Generate comprehensive text report
        self._generate_text_report(all_motifs_combined, stats_summary, genome_name)
        
        print(f"\n🎉 Analysis Complete!")
        print(f"📋 Files generated in {self.output_dir}:")
        print(f"   • {genome_name}_NBDFinder_Detailed_Results.csv - Complete motif detection results")
        print(f"   • {genome_name}_NBDFinder_Summary_Statistics.csv - Summary statistics")
        print(f"   • {genome_name}_NBDFinder_Motif_Distribution.csv - Motif type distribution")
        print(f"   • {genome_name}_NBDFinder_Analysis_Report.txt - Comprehensive text report")
    
    def _calculate_genome_stats(self, all_motifs, sequence_stats, genome_name):
        """Calculate comprehensive genome statistics"""
        total_motifs = len(all_motifs)
        total_sequences = len(sequence_stats)
        
        # Count by type
        type_counts = {}
        total_coverage = 0
        
        for motif in all_motifs:
            motif_type = motif.get('Class', 'Unknown')
            type_counts[motif_type] = type_counts.get(motif_type, 0) + 1
            total_coverage += motif.get('Length', 0)
        
        unique_types = len(type_counts)
        detection_rate = (unique_types / len(self.all_motif_types)) * 100
        
        # Score statistics
        scores = [float(motif.get('Score', 0)) for motif in all_motifs if motif.get('Score')]
        score_stats = {
            'mean_score': np.mean(scores) if scores else 0,
            'median_score': np.median(scores) if scores else 0,
            'std_score': np.std(scores) if scores else 0,
            'min_score': np.min(scores) if scores else 0,
            'max_score': np.max(scores) if scores else 0
        }
        
        return {
            'genome_name': genome_name,
            'total_sequences_analyzed': total_sequences,
            'total_motifs_detected': total_motifs,
            'unique_motif_types': unique_types,
            'motif_type_detection_rate_percent': detection_rate,
            'total_coverage_bp': total_coverage,
            'average_motifs_per_sequence': total_motifs / total_sequences if total_sequences > 0 else 0,
            **score_stats,
            **{f"{k}_count": v for k, v in type_counts.items()}
        }
    
    def _generate_text_report(self, all_motifs, stats_summary, genome_name):
        """Generate comprehensive text report"""
        report_path = os.path.join(self.output_dir, f"{genome_name}_NBDFinder_Analysis_Report.txt")
        
        with open(report_path, 'w') as f:
            f.write("="*80 + "\n")
            f.write(f"NBDFinder {self.version} - S. cerevisiae Genome Analysis Report\n")
            f.write("="*80 + "\n\n")
            
            f.write(f"ANALYSIS OVERVIEW\n")
            f.write("-" * 40 + "\n")
            f.write(f"Genome: {genome_name}\n")
            f.write(f"Analysis Date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"NBDFinder Version: {self.version}\n\n")
            
            f.write(f"SUMMARY STATISTICS\n")
            f.write("-" * 40 + "\n")
            f.write(f"Total Sequences Analyzed: {stats_summary['total_sequences_analyzed']}\n")
            f.write(f"Total Motifs Detected: {stats_summary['total_motifs_detected']}\n")
            f.write(f"Unique Motif Types: {stats_summary['unique_motif_types']}\n")
            f.write(f"Detection Rate: {stats_summary['motif_type_detection_rate_percent']:.1f}%\n")
            f.write(f"Total Coverage: {stats_summary['total_coverage_bp']:,} bp\n\n")
            
            f.write(f"MOTIF TYPE DISTRIBUTION\n")
            f.write("-" * 40 + "\n")
            
            # Count motif types
            type_counts = {}
            for motif in all_motifs:
                motif_type = motif.get('Class', 'Unknown')
                type_counts[motif_type] = type_counts.get(motif_type, 0) + 1
            
            for motif_type in sorted(type_counts.keys()):
                count = type_counts[motif_type]
                percentage = (count / len(all_motifs)) * 100 if all_motifs else 0
                f.write(f"{motif_type:25}: {count:6,} ({percentage:5.1f}%)\n")
            
            f.write("\n" + "="*80 + "\n")
            f.write("Analysis completed successfully with NBDFinder 2.0\n")
            f.write("For more information: https://github.com/VRYella/NBDFinder\n")
            f.write("="*80 + "\n")
        
        print(f"📄 Detailed text report saved: {report_path}")

def main():
    """Main function for S. cerevisiae genome analysis"""
    print("🧬 NBDFinder S. cerevisiae Genome Analysis")
    print("=" * 50)
    
    # Initialize analyzer
    analyzer = YeastGenomeAnalyzer()
    
    # Path to the S. cerevisiae genome file
    genome_file = "example_inputs/genome_data/GCF_000146045.2_R64_genomic.fna"
    
    if not os.path.exists(genome_file):
        print(f"❌ Genome file not found: {genome_file}")
        print("Please ensure the S. cerevisiae genome file is in the correct location.")
        return None
    
    # Read the genome sequences
    print(f"📖 Reading S. cerevisiae genome from: {genome_file}")
    sequences = analyzer.read_fasta_file(genome_file)
    
    if not sequences:
        print("❌ Failed to read genome file")
        return None
    
    print(f"📊 Found {len(sequences)} chromosome sequences")
    
    # Analyze each chromosome
    all_results = {}
    total_start = time.time()
    
    for i, (seq_name, sequence) in enumerate(sequences.items(), 1):
        print(f"\n🧬 Processing chromosome {i}/{len(sequences)}: {seq_name}")
        print(f"📏 Length: {len(sequence):,} bp")
        
        # Analyze this chromosome
        motifs = analyzer.analyze_sequence_chunk(sequence, seq_name)
        all_results[seq_name] = motifs
        
        print(f"✅ Found {len(motifs)} motifs in {seq_name}")
    
    total_time = time.time() - total_start
    print(f"\n⏱️  Total analysis time: {total_time:.1f} seconds")
    
    # Generate comprehensive report
    analyzer.generate_comprehensive_report(all_results, "S_cerevisiae")
    
    return all_results

if __name__ == "__main__":
    results = main()