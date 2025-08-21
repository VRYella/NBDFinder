#!/usr/bin/env python3
"""
NBDFinder S. cerevisiae Single Chromosome Analysis Script
=========================================================

Analyzes the first chromosome of S. cerevisiae genome for Non-B DNA structures 
and generates comprehensive results in the genomic data folder.

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
except ImportError as e:
    print(f"Error importing NBDFinder modules: {e}")
    print("Please ensure NBDFinder is properly installed.")
    sys.exit(1)

class YeastChromosomeAnalyzer:
    """Analyzer for S. cerevisiae chromosome NBD structures"""
    
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
        
        print(f"🧬 NBDFinder {self.version} Yeast Chromosome Analyzer Initialized")
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
    
    def analyze_sequence_optimized(self, sequence, sequence_name, chunk_size=25000):
        """Analyze a sequence with optimized chunking"""
        print(f"\n🔬 Analyzing: {sequence_name}")
        print(f"📏 Length: {len(sequence):,} bp")
        
        # Clean sequence
        sequence = sequence.upper().replace('N', '').replace(' ', '').replace('\n', '')
        
        if len(sequence) < 10:
            print(f"⚠️ Sequence too short: {len(sequence)} bp")
            return []
        
        analyzed_motifs = []
        
        if len(sequence) <= chunk_size:
            # Analyze the whole sequence if it's small enough
            print("🔬 Analyzing as single sequence...")
            motifs = all_motifs(sequence, nonoverlap=False, report_hotspots=True, sequence_name=sequence_name)
            analyzed_motifs.extend(motifs)
        else:
            # For large sequences, analyze in overlapping chunks
            print(f"🧩 Processing in chunks of {chunk_size:,} bp...")
            chunk_count = 0
            overlap = 2000  # Smaller overlap for speed
            
            for i in range(0, len(sequence), chunk_size - overlap):
                chunk_start = i
                chunk_end = min(i + chunk_size, len(sequence))
                chunk_seq = sequence[chunk_start:chunk_end]
                
                print(f"  Processing chunk {chunk_count + 1}: {chunk_start:,}-{chunk_end:,}")
                
                # Analyze this chunk with basic settings for speed
                chunk_motifs = all_motifs(chunk_seq, nonoverlap=False, report_hotspots=False, sequence_name=sequence_name)
                
                # Adjust coordinates to global sequence coordinates
                for motif in chunk_motifs:
                    if isinstance(motif, dict):
                        motif['Start'] = motif.get('Start', 0) + chunk_start
                        motif['End'] = motif.get('End', 0) + chunk_start
                        motif['Sequence Name'] = sequence_name
                
                analyzed_motifs.extend(chunk_motifs)
                chunk_count += 1
                
                # Progress update every 5 chunks
                if chunk_count % 5 == 0:
                    print(f"    ✅ Completed {chunk_count} chunks, found {len(analyzed_motifs)} motifs so far")
        
        # Clean subtypes
        for motif in analyzed_motifs:
            if isinstance(motif, dict) and 'Subtype' in motif and motif['Subtype']:
                motif['Subtype'] = str(motif['Subtype']).replace('_', ' ').strip()
        
        print(f"🎯 Total motifs found: {len(analyzed_motifs)}")
        return analyzed_motifs
    
    def generate_results(self, motifs, sequence_name, genome_name="S_cerevisiae_Chr1"):
        """Generate analysis results"""
        os.makedirs(self.output_dir, exist_ok=True)
        
        print(f"\n📊 Generating results in: {self.output_dir}")
        
        # 1. Save detailed motif results
        if motifs:
            df_motifs = pd.DataFrame(motifs)
            motifs_file = os.path.join(self.output_dir, f"{genome_name}_NBDFinder_Detailed_Results.csv")
            df_motifs.to_csv(motifs_file, index=False)
            print(f"✅ Detailed motifs saved: {motifs_file}")
            print(f"   📊 {len(motifs)} motifs detected")
        
        # 2. Generate summary statistics
        stats_summary = self._calculate_stats(motifs, sequence_name, genome_name)
        stats_file = os.path.join(self.output_dir, f"{genome_name}_NBDFinder_Summary_Statistics.csv")
        pd.DataFrame([stats_summary]).to_csv(stats_file, index=False)
        print(f"✅ Summary statistics saved: {stats_file}")
        
        # 3. Generate motif type distribution
        if motifs:
            type_counts = {}
            for motif in motifs:
                if isinstance(motif, dict):
                    motif_type = motif.get('Class', 'Unknown')
                    type_counts[motif_type] = type_counts.get(motif_type, 0) + 1
            
            distribution_df = pd.DataFrame(list(type_counts.items()), columns=['Motif_Type', 'Count'])
            distribution_df = distribution_df.sort_values('Count', ascending=False)
            dist_file = os.path.join(self.output_dir, f"{genome_name}_NBDFinder_Motif_Distribution.csv")
            distribution_df.to_csv(dist_file, index=False)
            print(f"✅ Motif distribution saved: {dist_file}")
        
        # 4. Generate positional matrix (simplified)
        position_matrix = self._create_position_matrix(motifs)
        if position_matrix is not None and not position_matrix.empty:
            pos_file = os.path.join(self.output_dir, f"{genome_name}_NBDFinder_Positional_Matrix.csv")
            position_matrix.to_csv(pos_file, index=False)
            print(f"✅ Positional matrix saved: {pos_file}")
        
        # 5. Generate comprehensive text report
        self._generate_text_report(motifs, stats_summary, genome_name)
        
        print(f"\n🎉 Analysis Complete!")
        print(f"📋 Files generated in {self.output_dir}:")
        print(f"   • {genome_name}_NBDFinder_Detailed_Results.csv - Complete motif detection results")
        print(f"   • {genome_name}_NBDFinder_Summary_Statistics.csv - Summary statistics")
        print(f"   • {genome_name}_NBDFinder_Motif_Distribution.csv - Motif type distribution")
        print(f"   • {genome_name}_NBDFinder_Positional_Matrix.csv - Nucleotide positional data")
        print(f"   • {genome_name}_NBDFinder_Analysis_Report.txt - Comprehensive text report")
    
    def _create_position_matrix(self, motifs):
        """Create a simplified positional matrix"""
        if not motifs:
            return None
        
        # Create basic positional data
        position_data = []
        for motif in motifs:
            if isinstance(motif, dict):
                position_data.append({
                    'Motif_Type': motif.get('Class', 'Unknown'),
                    'Position': motif.get('Start', 0),
                    'Length': motif.get('Length', 0),
                    'Score': motif.get('Score', 0)
                })
        
        return pd.DataFrame(position_data)
    
    def _calculate_stats(self, motifs, sequence_name, genome_name):
        """Calculate comprehensive statistics"""
        total_motifs = len(motifs)
        
        # Count by type
        type_counts = {}
        total_coverage = 0
        
        for motif in motifs:
            if isinstance(motif, dict):
                motif_type = motif.get('Class', 'Unknown')
                type_counts[motif_type] = type_counts.get(motif_type, 0) + 1
                total_coverage += motif.get('Length', 0)
        
        unique_types = len(type_counts)
        detection_rate = (unique_types / len(self.all_motif_types)) * 100
        
        # Score statistics
        scores = []
        for motif in motifs:
            if isinstance(motif, dict) and motif.get('Score'):
                try:
                    scores.append(float(motif.get('Score', 0)))
                except (ValueError, TypeError):
                    pass
        
        score_stats = {
            'mean_score': np.mean(scores) if scores else 0,
            'median_score': np.median(scores) if scores else 0,
            'std_score': np.std(scores) if scores else 0,
            'min_score': np.min(scores) if scores else 0,
            'max_score': np.max(scores) if scores else 0
        }
        
        return {
            'genome_name': genome_name,
            'sequence_name': sequence_name,
            'total_motifs_detected': total_motifs,
            'unique_motif_types': unique_types,
            'motif_type_detection_rate_percent': detection_rate,
            'total_coverage_bp': total_coverage,
            **score_stats,
            **{f"{k}_count": v for k, v in type_counts.items()}
        }
    
    def _generate_text_report(self, motifs, stats_summary, genome_name):
        """Generate comprehensive text report"""
        report_path = os.path.join(self.output_dir, f"{genome_name}_NBDFinder_Analysis_Report.txt")
        
        with open(report_path, 'w') as f:
            f.write("="*80 + "\n")
            f.write(f"NBDFinder {self.version} - S. cerevisiae Chromosome Analysis Report\n")
            f.write("="*80 + "\n\n")
            
            f.write(f"ANALYSIS OVERVIEW\n")
            f.write("-" * 40 + "\n")
            f.write(f"Genome: {genome_name}\n")
            f.write(f"Sequence: {stats_summary.get('sequence_name', 'Unknown')}\n")
            f.write(f"Analysis Date: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"NBDFinder Version: {self.version}\n\n")
            
            f.write(f"SUMMARY STATISTICS\n")
            f.write("-" * 40 + "\n")
            f.write(f"Total Motifs Detected: {stats_summary['total_motifs_detected']}\n")
            f.write(f"Unique Motif Types: {stats_summary['unique_motif_types']}\n")
            f.write(f"Detection Rate: {stats_summary['motif_type_detection_rate_percent']:.1f}%\n")
            f.write(f"Total Coverage: {stats_summary['total_coverage_bp']:,} bp\n\n")
            
            if stats_summary.get('mean_score', 0) > 0:
                f.write(f"SCORE STATISTICS\n")
                f.write("-" * 40 + "\n")
                f.write(f"Mean Score: {stats_summary.get('mean_score', 0):.3f}\n")
                f.write(f"Median Score: {stats_summary.get('median_score', 0):.3f}\n")
                f.write(f"Score Range: {stats_summary.get('min_score', 0):.3f} - {stats_summary.get('max_score', 0):.3f}\n\n")
            
            f.write(f"MOTIF TYPE DISTRIBUTION\n")
            f.write("-" * 40 + "\n")
            
            # Count motif types
            type_counts = {}
            for motif in motifs:
                if isinstance(motif, dict):
                    motif_type = motif.get('Class', 'Unknown')
                    type_counts[motif_type] = type_counts.get(motif_type, 0) + 1
            
            for motif_type in sorted(type_counts.keys()):
                count = type_counts[motif_type]
                percentage = (count / len(motifs)) * 100 if motifs else 0
                f.write(f"{motif_type:25}: {count:6,} ({percentage:5.1f}%)\n")
            
            f.write("\n" + "="*80 + "\n")
            f.write("Analysis completed successfully with NBDFinder 2.0\n")
            f.write("For more information: https://github.com/VRYella/NBDFinder\n")
            f.write("="*80 + "\n")
        
        print(f"📄 Detailed text report saved: {report_path}")

def main():
    """Main function for S. cerevisiae chromosome analysis"""
    print("🧬 NBDFinder S. cerevisiae Chromosome Analysis")
    print("=" * 50)
    
    # Initialize analyzer
    analyzer = YeastChromosomeAnalyzer()
    
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
    
    # Analyze the first chromosome (chromosome I is typically smallest)
    first_chromosome = list(sequences.items())[0]
    seq_name, sequence = first_chromosome
    
    print(f"\n🧬 Analyzing first chromosome: {seq_name}")
    print(f"📏 Length: {len(sequence):,} bp")
    
    # Analyze this chromosome
    start_time = time.time()
    motifs = analyzer.analyze_sequence_optimized(sequence, seq_name)
    analysis_time = time.time() - start_time
    
    print(f"⏱️  Total analysis time: {analysis_time:.1f} seconds")
    print(f"✅ Found {len(motifs)} motifs in {seq_name}")
    
    # Generate comprehensive report
    analyzer.generate_results(motifs, seq_name, "S_cerevisiae_Chr1")
    
    return motifs

if __name__ == "__main__":
    results = main()