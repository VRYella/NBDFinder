#!/usr/bin/env python3
"""
Comprehensive Disease Expansion Loci Analysis with NBDFinder
==========================================================

This script performs rigorous comparative analysis of all genes in 
repeat_expansion_loci_annotated.fa using the NBDFinder platform.

Features:
- Complete motif detection across 10 classes and 22 subclasses
- Disease-specific repeat expansion analysis
- Comparative analysis across all 153 disease-related genes
- Statistical summaries and visualizations
- Comprehensive export in multiple formats

Author: Dr. Venkata Rajesh Yella
Date: 2024
License: Academic Use
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

class DiseaseExpansionAnalyzer:
    """
    Comprehensive analyzer for disease expansion loci using NBDFinder
    """
    
    def __init__(self, output_dir="disease_expansion_analysis_results"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(exist_ok=True)
        
        # Create subdirectories
        self.individual_dir = self.output_dir / "individual_sequences"
        self.individual_dir.mkdir(exist_ok=True)
        
        self.plots_dir = self.output_dir / "visualizations"
        self.plots_dir.mkdir(exist_ok=True)
        
        self.summary_dir = self.output_dir / "summary_reports"
        self.summary_dir.mkdir(exist_ok=True)
        
        self.analysis_results = []
        self.disease_detector = AdvancedDiseaseDetector()
        
        print(f"🧬 Disease Expansion Analyzer initialized")
        print(f"📁 Output directory: {self.output_dir}")
        
    def load_sequences(self, fasta_file="repeat_expansion_loci_annotated.fa"):
        """Load sequences from the disease expansion FASTA file"""
        print(f"\n📖 Loading sequences from {fasta_file}...")
        
        with open(fasta_file, 'r') as f:
            content = f.read()
            
        self.sequences, self.names = parse_fasta_multi(content)
        
        print(f"✅ Loaded {len(self.sequences)} disease-related sequences")
        print(f"🧬 Average sequence length: {np.mean([len(seq) for seq in self.sequences]):.0f} bp")
        
        return len(self.sequences)
    
    def analyze_single_sequence(self, sequence, sequence_name, index):
        """Perform comprehensive analysis on a single sequence"""
        print(f"\n🔬 Analyzing sequence {index + 1}/{len(self.sequences)}: {sequence_name[:60]}...")
        
        start_time = time.time()
        
        # Basic sequence statistics
        seq_length = len(sequence)
        gc_percent = gc_content(sequence)
        
        # Extract gene information from sequence name
        parts = sequence_name.split('|')
        gene_symbol = parts[1] if len(parts) > 1 else "Unknown"
        disease_info = parts[2] if len(parts) > 2 else "Unknown"
        
        # Comprehensive motif analysis
        motifs = all_motifs(sequence, sequence_name=gene_symbol, report_hotspots=False)
        
        # Disease-specific analysis
        disease_motifs = self.disease_detector.detect_pathogenic_repeats(sequence, gene_symbol)
        
        # Count motifs by class
        motif_classes = {}
        for motif in motifs:
            motif_class = motif.get('Class', 'Unknown')
            motif_classes[motif_class] = motif_classes.get(motif_class, 0) + 1
        
        # Analysis timing
        analysis_time = time.time() - start_time
        
        # Store comprehensive results
        result = {
            'sequence_index': index + 1,
            'sequence_name': sequence_name,
            'gene_symbol': gene_symbol,
            'disease_info': disease_info,
            'sequence_length': seq_length,
            'gc_content': gc_percent,
            'total_motifs': len(motifs),
            'total_disease_motifs': len(disease_motifs),
            'motifs': motifs,
            'disease_motifs': disease_motifs,
            'motif_classes': motif_classes,
            'analysis_time': analysis_time,
            'timestamp': datetime.now().isoformat()
        }
        
        # Save individual results
        self._save_individual_results(result)
        
        print(f"   ✅ Found {len(motifs)} total motifs, {len(disease_motifs)} disease-specific")
        print(f"   ⏱️  Analysis completed in {analysis_time:.2f} seconds")
        
        return result
    
    def _save_individual_results(self, result):
        """Save individual sequence analysis results"""
        gene_symbol = result['gene_symbol']
        safe_filename = "".join(c for c in gene_symbol if c.isalnum() or c in (' ', '-', '_')).rstrip()
        
        # Save detailed JSON
        json_file = self.individual_dir / f"{safe_filename}_analysis.json"
        with open(json_file, 'w') as f:
            # Create JSON-serializable version
            json_result = {k: v for k, v in result.items() if k != 'motifs' and k != 'disease_motifs'}
            json_result['motif_count'] = len(result['motifs'])
            json_result['disease_motif_count'] = len(result['disease_motifs'])
            json.dump(json_result, f, indent=2)
        
        # Save motifs as CSV
        if result['motifs']:
            motifs_df = pd.DataFrame(result['motifs'])
            csv_file = self.individual_dir / f"{safe_filename}_motifs.csv"
            motifs_df.to_csv(csv_file, index=False)
        
        # Save disease motifs as CSV
        if result['disease_motifs']:
            disease_df = pd.DataFrame(result['disease_motifs'])
            disease_csv = self.individual_dir / f"{safe_filename}_disease_motifs.csv"
            disease_df.to_csv(disease_csv, index=False)
    
    def run_comprehensive_analysis(self, limit=20):
        """Run complete analysis on all sequences (or limited subset for demo)"""
        sequences_to_analyze = min(len(self.sequences), limit) if limit else len(self.sequences)
        print(f"\n🚀 Starting comprehensive analysis of {sequences_to_analyze} sequences...")
        
        start_time = time.time()
        
        for i in range(sequences_to_analyze):
            sequence = self.sequences[i]
            name = self.names[i]
            try:
                result = self.analyze_single_sequence(sequence, name, i)
                self.analysis_results.append(result)
                
                # Progress update every 5 sequences
                if (i + 1) % 5 == 0:
                    elapsed = time.time() - start_time
                    avg_time = elapsed / (i + 1)
                    remaining = avg_time * (sequences_to_analyze - i - 1)
                    print(f"📊 Progress: {i + 1}/{sequences_to_analyze} ({(i+1)/sequences_to_analyze*100:.1f}%)")
                    print(f"⏱️  Estimated time remaining: {remaining/60:.1f} minutes")
                    
            except Exception as e:
                print(f"❌ Error analyzing {name}: {str(e)}")
                continue
        
        total_time = time.time() - start_time
        print(f"\n🎉 Analysis completed in {total_time/60:.1f} minutes")
        print(f"✅ Successfully analyzed {len(self.analysis_results)} sequences")
        
        return self.analysis_results
    
    def generate_comparative_analysis(self):
        """Generate comprehensive comparative analysis"""
        print(f"\n📊 Generating comparative analysis...")
        
        # Create comprehensive summary DataFrame
        summary_data = []
        
        for result in self.analysis_results:
            base_data = {
                'Gene_Symbol': result['gene_symbol'],
                'Disease_Info': result['disease_info'],
                'Sequence_Length': result['sequence_length'],
                'GC_Content': result['gc_content'],
                'Total_Motifs': result['total_motifs'],
                'Disease_Motifs': result['total_disease_motifs'],
                'Motif_Density': result['total_motifs'] / result['sequence_length'] * 1000,  # per kb
                'Analysis_Time': result['analysis_time']
            }
            
            # Add motif class counts
            for class_name, count in result['motif_classes'].items():
                base_data[f'{class_name}_Count'] = count
                base_data[f'{class_name}_Density'] = count / result['sequence_length'] * 1000
            
            summary_data.append(base_data)
        
        self.summary_df = pd.DataFrame(summary_data)
        
        # Save comprehensive results
        excel_file = self.summary_dir / "comprehensive_disease_expansion_analysis.xlsx"
        with pd.ExcelWriter(excel_file, engine='openpyxl') as writer:
            self.summary_df.to_excel(writer, sheet_name='Summary', index=False)
            
            # Individual sequence details
            for result in self.analysis_results[:50]:  # Limit to first 50 for Excel
                if result['motifs']:
                    motif_df = pd.DataFrame(result['motifs'])
                    sheet_name = result['gene_symbol'][:31]  # Excel sheet name limit
                    motif_df.to_excel(writer, sheet_name=sheet_name, index=False)
        
        print(f"  ✅ Saved comprehensive analysis to {excel_file}")
        
        # Generate statistics
        self._generate_statistics()
        
        return self.summary_df
    
    def _generate_statistics(self):
        """Generate detailed statistics"""
        print(f"\n📈 Generating statistical analysis...")
        
        stats = {
            'total_sequences': len(self.analysis_results),
            'total_motifs': self.summary_df['Total_Motifs'].sum(),
            'total_disease_motifs': self.summary_df['Disease_Motifs'].sum(),
            'avg_motifs_per_gene': self.summary_df['Total_Motifs'].mean(),
            'avg_disease_motifs_per_gene': self.summary_df['Disease_Motifs'].mean(),
            'avg_sequence_length': self.summary_df['Sequence_Length'].mean(),
            'avg_gc_content': self.summary_df['GC_Content'].mean(),
            'genes_with_disease_motifs': (self.summary_df['Disease_Motifs'] > 0).sum(),
            'max_motifs_gene': self.summary_df.loc[self.summary_df['Total_Motifs'].idxmax(), 'Gene_Symbol'],
            'max_motifs_count': self.summary_df['Total_Motifs'].max(),
        }
        
        # Save statistics
        stats_file = self.summary_dir / "analysis_statistics.json"
        with open(stats_file, 'w') as f:
            json.dump(stats, f, indent=2)
        
        # Print summary
        print(f"📊 Analysis Statistics:")
        print(f"   🧬 Total sequences analyzed: {stats['total_sequences']}")
        print(f"   🔍 Total motifs found: {stats['total_motifs']}")
        print(f"   🦠 Disease-specific motifs: {stats['total_disease_motifs']}")
        print(f"   📏 Average motifs per gene: {stats['avg_motifs_per_gene']:.1f}")
        print(f"   🎯 Genes with disease motifs: {stats['genes_with_disease_motifs']}")
        print(f"   🏆 Gene with most motifs: {stats['max_motifs_gene']} ({stats['max_motifs_count']} motifs)")
        
        return stats
    
    def create_visualizations(self):
        """Create comprehensive visualizations"""
        print(f"\n🎨 Creating visualizations...")
        
        # Set up plotting style
        plt.style.use('default')
        sns.set_palette("husl")
        
        # 1. Distribution of motifs per gene
        plt.figure(figsize=(12, 8))
        plt.subplot(2, 2, 1)
        plt.hist(self.summary_df['Total_Motifs'], bins=30, alpha=0.7, edgecolor='black')
        plt.xlabel('Total Motifs per Gene')
        plt.ylabel('Number of Genes')
        plt.title('Distribution of Motifs per Gene')
        plt.grid(True, alpha=0.3)
        
        # 2. GC Content vs Motif Density
        plt.subplot(2, 2, 2)
        plt.scatter(self.summary_df['GC_Content'], self.summary_df['Motif_Density'], 
                   alpha=0.6, s=50)
        plt.xlabel('GC Content (%)')
        plt.ylabel('Motif Density (per kb)')
        plt.title('GC Content vs Motif Density')
        plt.grid(True, alpha=0.3)
        
        # 3. Disease motifs distribution
        plt.subplot(2, 2, 3)
        disease_genes = self.summary_df[self.summary_df['Disease_Motifs'] > 0]
        if not disease_genes.empty:
            plt.hist(disease_genes['Disease_Motifs'], bins=20, alpha=0.7, 
                    color='red', edgecolor='black')
            plt.xlabel('Disease Motifs per Gene')
            plt.ylabel('Number of Genes')
            plt.title('Distribution of Disease Motifs')
        else:
            plt.text(0.5, 0.5, 'No Disease Motifs Found', 
                    transform=plt.gca().transAxes, ha='center', va='center')
            plt.title('Distribution of Disease Motifs')
        plt.grid(True, alpha=0.3)
        
        # 4. Top genes by motif count
        plt.subplot(2, 2, 4)
        top_genes = self.summary_df.nlargest(10, 'Total_Motifs')
        plt.barh(range(len(top_genes)), top_genes['Total_Motifs'])
        plt.yticks(range(len(top_genes)), top_genes['Gene_Symbol'])
        plt.xlabel('Total Motifs')
        plt.title('Top 10 Genes by Motif Count')
        plt.grid(True, alpha=0.3)
        
        plt.tight_layout()
        plt.savefig(self.plots_dir / 'overview_analysis.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # 2. Motif class distribution
        self._create_motif_class_visualization()
        
        # 3. Heatmap of motif densities
        self._create_heatmap_visualization()
        
        print(f"  ✅ Visualizations saved to {self.plots_dir}")
    
    def _create_motif_class_visualization(self):
        """Create motif class distribution visualization"""
        # Collect all motif classes
        all_classes = set()
        for result in self.analysis_results:
            all_classes.update(result['motif_classes'].keys())
        
        if not all_classes:
            print("  ⚠️  No motif classes found for visualization")
            return
        
        # Create matrix for heatmap
        class_matrix = []
        gene_names = []
        
        for result in self.analysis_results[:50]:  # Limit for readability
            gene_names.append(result['gene_symbol'])
            row = []
            for class_name in sorted(all_classes):
                density = result['motif_classes'].get(class_name, 0) / result['sequence_length'] * 1000
                row.append(density)
            class_matrix.append(row)
        
        # Create heatmap
        plt.figure(figsize=(14, 10))
        sns.heatmap(class_matrix, 
                   xticklabels=sorted(all_classes),
                   yticklabels=gene_names,
                   cmap='viridis',
                   cbar_kws={'label': 'Motif Density (per kb)'})
        plt.title('Motif Class Distribution Across Genes')
        plt.xlabel('Motif Classes')
        plt.ylabel('Genes')
        plt.xticks(rotation=45)
        plt.yticks(rotation=0, fontsize=8)
        plt.tight_layout()
        plt.savefig(self.plots_dir / 'motif_class_heatmap.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    def _create_heatmap_visualization(self):
        """Create comprehensive heatmap of key metrics"""
        # Select top 30 genes by motif count for readability
        top_genes = self.summary_df.nlargest(30, 'Total_Motifs')
        
        # Prepare data for heatmap
        heatmap_data = top_genes[['Total_Motifs', 'Disease_Motifs', 'GC_Content', 
                                 'Sequence_Length', 'Motif_Density']].copy()
        
        # Normalize data for better visualization
        from sklearn.preprocessing import StandardScaler
        scaler = StandardScaler()
        heatmap_normalized = pd.DataFrame(
            scaler.fit_transform(heatmap_data),
            columns=heatmap_data.columns,
            index=top_genes['Gene_Symbol']
        )
        
        plt.figure(figsize=(10, 12))
        sns.heatmap(heatmap_normalized.T, 
                   cmap='RdYlBu_r',
                   center=0,
                   cbar_kws={'label': 'Normalized Value'})
        plt.title('Normalized Metrics for Top 30 Genes by Motif Count')
        plt.xlabel('Genes')
        plt.ylabel('Metrics')
        plt.xticks(rotation=90)
        plt.tight_layout()
        plt.savefig(self.plots_dir / 'metrics_heatmap.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    def generate_summary_report(self):
        """Generate comprehensive summary report"""
        print(f"\n📝 Generating summary report...")
        
        report_file = self.summary_dir / "analysis_report.md"
        
        with open(report_file, 'w') as f:
            f.write("# Disease Expansion Loci Analysis Report\n\n")
            f.write(f"**Analysis Date:** {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
            
            f.write("## Summary\n\n")
            f.write(f"- **Total sequences analyzed:** {len(self.analysis_results)}\n")
            f.write(f"- **Total motifs detected:** {self.summary_df['Total_Motifs'].sum()}\n")
            f.write(f"- **Disease-specific motifs:** {self.summary_df['Disease_Motifs'].sum()}\n")
            f.write(f"- **Average motifs per gene:** {self.summary_df['Total_Motifs'].mean():.1f}\n")
            f.write(f"- **Genes with disease motifs:** {(self.summary_df['Disease_Motifs'] > 0).sum()}\n\n")
            
            f.write("## Top Genes by Motif Count\n\n")
            top_10 = self.summary_df.nlargest(10, 'Total_Motifs')
            for _, row in top_10.iterrows():
                f.write(f"- **{row['Gene_Symbol']}**: {row['Total_Motifs']} motifs "
                       f"({row['Disease_Motifs']} disease-specific)\n")
            
            f.write("\n## Files Generated\n\n")
            f.write("- Individual sequence analyses in `individual_sequences/`\n")
            f.write("- Comprehensive Excel report: `comprehensive_disease_expansion_analysis.xlsx`\n")
            f.write("- Visualizations in `visualizations/`\n")
            f.write("- Analysis statistics: `analysis_statistics.json`\n")
        
        print(f"  ✅ Summary report saved to {report_file}")


def main():
    """Main execution function"""
    print("🧬 NBDFinder Disease Expansion Loci Analysis")
    print("=" * 50)
    
    # Initialize analyzer
    analyzer = DiseaseExpansionAnalyzer()
    
    # Load sequences
    if analyzer.load_sequences():
        # Run comprehensive analysis
        results = analyzer.run_comprehensive_analysis()
        
        if results:
            # Generate comparative analysis
            analyzer.generate_comparative_analysis()
            
            # Create visualizations
            analyzer.create_visualizations()
            
            # Generate summary report
            analyzer.generate_summary_report()
            
            print(f"\n🎉 Analysis completed successfully!")
            print(f"📁 Results saved to: {analyzer.output_dir}")
        else:
            print("❌ No results generated")
    else:
        print("❌ Failed to load sequences")


if __name__ == "__main__":
    main()