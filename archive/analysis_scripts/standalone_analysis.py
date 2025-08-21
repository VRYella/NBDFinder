#!/usr/bin/env python3
"""
NBDFinder 2.0: Standalone Comprehensive Analysis Script
======================================================

The most advanced non-B DNA structure detection platform.
Analyzes all 19 motif types with professional output and visualization.

Author: Dr. Venkata Rajesh Yella
Updated: 2024
"""

import os
import sys
import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# Import NBDFinder modules
try:
    from motifs import all_motifs
    from disease_motifs import AdvancedDiseaseDetector
    from advanced_clustering import AdvancedClusterAnalyzer
    from ml_predictor import AdvancedMLPredictor
except ImportError as e:
    print(f"Error importing NBDFinder modules: {e}")
    print("Please ensure NBDFinder is properly installed.")
    sys.exit(1)

class NBDFinderStandalone:
    """Standalone NBDFinder 2.0 analysis class"""
    
    def __init__(self):
        """Initialize the standalone analyzer"""
        self.version = "2.0.0"
        self.all_motif_types = [
            "Bipartite G4", "Bulged G4", "Canonical G4", "Imperfect G4", 
            "Multimeric G4", "Relaxed G4", "G-Triplex", "AC-Motif", "i-Motif",
            "Curved DNA", "Z-DNA", "eGZ (Extruded-G)", "Cruciform", "R-Loop",
            "Slipped DNA", "Sticky DNA", "Triplex DNA", "Hybrid", "Non-B DNA Clusters"
        ]
        
        # Initialize advanced analyzers
        self.disease_detector = AdvancedDiseaseDetector()
        self.cluster_analyzer = AdvancedClusterAnalyzer()
        self.ml_predictor = AdvancedMLPredictor()
        
        print(f"🧬 NBDFinder {self.version} Standalone Analyzer Initialized")
        print(f"✅ Ready to analyze all {len(self.all_motif_types)} motif types")
    
    def analyze_sequence(self, sequence, sequence_name="Sequence", verbose=True):
        """
        Comprehensive analysis of a single sequence
        
        Args:
            sequence (str): DNA sequence to analyze
            sequence_name (str): Name for the sequence
            verbose (bool): Print progress information
            
        Returns:
            dict: Complete analysis results
        """
        if verbose:
            print(f"\n🔬 Analyzing: {sequence_name}")
            print(f"📏 Length: {len(sequence):,} bp")
        
        start_time = time.time()
        
        # Core motif detection
        motifs = all_motifs(sequence)
        
        # Clean subtypes (remove underscores)
        for motif in motifs:
            if 'Subtype' in motif and motif['Subtype']:
                # Simple cleaning - remove underscores and normalize
                motif['Subtype'] = str(motif['Subtype']).replace('_', ' ').strip()
        
        # Disease association analysis
        disease_motifs = self.disease_detector.detect_pathogenic_repeats(sequence)
        
        # Advanced clustering
        cluster_regions = self.cluster_analyzer.detect_advanced_clusters(motifs, len(sequence))
        
        # Machine learning predictions
        ml_predictions = {}
        if motifs:
            for motif in motifs[:5]:  # Sample a few motifs for ML analysis
                motif_type = motif.get('Class', 'general')
                try:
                    pred = self.ml_predictor.predict_structure_probability(sequence, motif_type)
                    ml_predictions[motif_type] = pred
                except Exception as e:
                    if verbose:
                        print(f"⚠️ ML prediction failed for {motif_type}: {e}")
                    ml_predictions[motif_type] = {'probability': 0.5, 'confidence': 0.0}
        
        # Calculate comprehensive statistics
        stats = self._calculate_comprehensive_stats(motifs, sequence, sequence_name)
        
        analysis_time = time.time() - start_time
        
        if verbose:
            print(f"⏱️  Analysis completed in {analysis_time*1000:.1f} ms")
            print(f"🎯 Detected {len(motifs)} total motifs")
            print(f"🏥 Found {len(disease_motifs)} disease-associated motifs")
            print(f"🌟 Identified {len(cluster_regions)} structural clusters")
        
        return {
            'sequence_name': sequence_name,
            'sequence_length': len(sequence),
            'analysis_time_ms': analysis_time * 1000,
            'motifs': motifs,
            'disease_motifs': disease_motifs,
            'clusters': cluster_regions,
            'ml_predictions': ml_predictions,
            'statistics': stats,
            'coverage_analysis': self._analyze_coverage(motifs, len(sequence))
        }
    
    def _calculate_comprehensive_stats(self, motifs, sequence, sequence_name):
        """Calculate comprehensive statistics for the analysis"""
        # Basic counts
        total_motifs = len(motifs)
        unique_types = len(set(m.get('Class', 'Unknown') for m in motifs))
        
        # Coverage analysis
        total_coverage = sum(m.get('Length', 0) for m in motifs)
        coverage_percent = (total_coverage / len(sequence) * 100) if len(sequence) > 0 else 0
        
        # Motif type distribution
        type_counts = {}
        for motif in motifs:
            motif_type = motif.get('Class', 'Unknown')
            type_counts[motif_type] = type_counts.get(motif_type, 0) + 1
        
        # Identify missing types
        detected_types = set(type_counts.keys())
        missing_types = set(self.all_motif_types) - detected_types
        
        # Score distribution
        scores = [float(m.get('Score', 0)) for m in motifs if m.get('Score')]
        score_stats = {
            'mean': np.mean(scores) if scores else 0,
            'median': np.median(scores) if scores else 0,
            'std': np.std(scores) if scores else 0,
            'min': np.min(scores) if scores else 0,
            'max': np.max(scores) if scores else 0
        }
        
        return {
            'sequence_name': sequence_name,
            'sequence_length': len(sequence),
            'total_motifs': total_motifs,
            'unique_types_detected': unique_types,
            'total_possible_types': len(self.all_motif_types),
            'detection_rate_percent': (unique_types / len(self.all_motif_types)) * 100,
            'total_coverage_bp': total_coverage,
            'coverage_percent': coverage_percent,
            'type_counts': type_counts,
            'missing_types': list(missing_types),
            'score_statistics': score_stats
        }
    
    def _analyze_coverage(self, motifs, sequence_length):
        """Analyze sequence coverage by motif types"""
        coverage_by_type = {}
        
        for motif_type in self.all_motif_types:
            type_motifs = [m for m in motifs if m.get('Class') == motif_type]
            total_length = sum(m.get('Length', 0) for m in type_motifs)
            coverage_percent = (total_length / sequence_length * 100) if sequence_length > 0 else 0
            
            coverage_by_type[motif_type] = {
                'count': len(type_motifs),
                'total_length': total_length,
                'coverage_percent': coverage_percent,
                'status': 'Detected' if len(type_motifs) > 0 else 'Not Found'
            }
        
        return coverage_by_type
    
    def generate_comprehensive_report(self, results, output_dir="NBDFinder_Results"):
        """Generate comprehensive analysis report"""
        os.makedirs(output_dir, exist_ok=True)
        
        print(f"\n📊 Generating comprehensive report in: {output_dir}")
        
        # 1. Save detailed motif results
        if results['motifs']:
            df_motifs = pd.DataFrame(results['motifs'])
            df_motifs.to_csv(f"{output_dir}/detailed_motifs.csv", index=False)
            print(f"✅ Detailed motifs saved: {len(results['motifs'])} motifs")
        
        # 2. Save disease analysis
        if results['disease_motifs']:
            df_disease = pd.DataFrame(results['disease_motifs'])
            df_disease.to_csv(f"{output_dir}/disease_motifs.csv", index=False)
            print(f"🏥 Disease motifs saved: {len(results['disease_motifs'])} motifs")
        
        # 3. Save comprehensive statistics
        stats_df = pd.DataFrame([results['statistics']])
        stats_df.to_csv(f"{output_dir}/analysis_statistics.csv", index=False)
        
        # 4. Save coverage analysis
        coverage_df = pd.DataFrame.from_dict(results['coverage_analysis'], orient='index')
        coverage_df.to_csv(f"{output_dir}/coverage_analysis.csv")
        
        # 5. Generate visualizations
        self._generate_visualizations(results, output_dir)
        
        # 6. Generate text report
        self._generate_text_report(results, output_dir)
        
        print(f"📈 Report generation complete!")
        return output_dir
    
    def _generate_visualizations(self, results, output_dir):
        """Generate comprehensive visualizations"""
        # Set style for publication quality
        plt.style.use('seaborn-v0_8-whitegrid')
        sns.set_palette("husl")
        
        # 1. Motif type distribution
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        fig.suptitle(f"NBDFinder 2.0 Analysis: {results['sequence_name']}", fontsize=16, fontweight='bold')
        
        # Coverage by type
        coverage_data = results['coverage_analysis']
        types = list(coverage_data.keys())
        counts = [coverage_data[t]['count'] for t in types]
        colors = ['green' if c > 0 else 'lightgray' for c in counts]
        
        axes[0,0].barh(types, counts, color=colors)
        axes[0,0].set_title('Motif Detection Overview (All 19 Types)', fontweight='bold')
        axes[0,0].set_xlabel('Number of Motifs Detected')
        
        # Add "Not Found" annotations
        for i, (type_name, count) in enumerate(zip(types, counts)):
            if count == 0:
                axes[0,0].text(0.1, i, 'Not Found', va='center', ha='left', 
                             fontweight='bold', color='red', fontsize=8)
        
        # Coverage percentage
        coverage_pcts = [coverage_data[t]['coverage_percent'] for t in types]
        axes[0,1].barh(types, coverage_pcts, color=colors)
        axes[0,1].set_title('Sequence Coverage by Motif Type', fontweight='bold')
        axes[0,1].set_xlabel('Coverage Percentage (%)')
        
        # Statistics summary
        stats = results['statistics']
        summary_data = {
            'Total Motifs': stats['total_motifs'],
            'Types Detected': f"{stats['unique_types_detected']}/19",
            'Detection Rate': f"{stats['detection_rate_percent']:.1f}%",
            'Coverage': f"{stats['coverage_percent']:.1f}%"
        }
        
        axes[1,0].pie(summary_data.values(), labels=summary_data.keys(), autopct='%1.1f%%')
        axes[1,0].set_title('Analysis Summary', fontweight='bold')
        
        # Score distribution
        if results['motifs']:
            scores = [float(m.get('Score', 0)) for m in results['motifs'] if m.get('Score')]
            if scores:
                axes[1,1].hist(scores, bins=20, alpha=0.7, edgecolor='black')
                axes[1,1].set_title('Score Distribution', fontweight='bold')
                axes[1,1].set_xlabel('Motif Scores')
                axes[1,1].set_ylabel('Frequency')
        
        plt.tight_layout()
        plt.savefig(f"{output_dir}/comprehensive_analysis.png", dpi=300, bbox_inches='tight')
        plt.savefig(f"{output_dir}/comprehensive_analysis.pdf", bbox_inches='tight')
        plt.close()
        
        print("📊 Comprehensive visualization saved (PNG & PDF)")
    
    def _generate_text_report(self, results, output_dir):
        """Generate detailed text report"""
        report_path = f"{output_dir}/analysis_report.txt"
        
        with open(report_path, 'w') as f:
            f.write("="*80 + "\n")
            f.write("NBDFinder 2.0 Comprehensive Analysis Report\n")
            f.write("="*80 + "\n\n")
            
            # Basic information
            stats = results['statistics']
            f.write(f"Sequence Name: {results['sequence_name']}\n")
            f.write(f"Sequence Length: {results['sequence_length']:,} bp\n")
            f.write(f"Analysis Time: {results['analysis_time_ms']:.1f} ms\n")
            f.write(f"NBDFinder Version: {self.version}\n\n")
            
            # Summary statistics
            f.write("SUMMARY STATISTICS\n")
            f.write("-" * 40 + "\n")
            f.write(f"Total Motifs Detected: {stats['total_motifs']}\n")
            f.write(f"Motif Types Detected: {stats['unique_types_detected']}/19 ({stats['detection_rate_percent']:.1f}%)\n")
            f.write(f"Sequence Coverage: {stats['coverage_percent']:.2f}%\n")
            f.write(f"Total Coverage: {stats['total_coverage_bp']:,} bp\n\n")
            
            # Motif type breakdown
            f.write("MOTIF TYPE BREAKDOWN\n")
            f.write("-" * 40 + "\n")
            coverage_data = results['coverage_analysis']
            for motif_type in self.all_motif_types:
                data = coverage_data[motif_type]
                status = "✅ DETECTED" if data['count'] > 0 else "❌ NOT FOUND"
                f.write(f"{motif_type:20} {status:12} Count: {data['count']:3} Coverage: {data['coverage_percent']:5.2f}%\n")
            
            # Disease motifs
            if results['disease_motifs']:
                f.write(f"\nDISEASE-ASSOCIATED MOTIFS\n")
                f.write("-" * 40 + "\n")
                f.write(f"Found {len(results['disease_motifs'])} disease-associated motifs:\n")
                for motif in results['disease_motifs']:
                    f.write(f"- {motif.get('Disease_Name', 'Unknown')}: {motif.get('Clinical_Significance', 'Unknown')}\n")
            
            # Clusters
            if results['clusters']:
                f.write(f"\nSTRUCTURAL CLUSTERS\n")
                f.write("-" * 40 + "\n")
                f.write(f"Identified {len(results['clusters'])} structural clusters\n")
            
            f.write("\n" + "="*80 + "\n")
            f.write("Analysis completed successfully with NBDFinder 2.0\n")
            f.write("For more information: https://github.com/VRYella/NBDFinder\n")
            f.write("="*80 + "\n")
        
        print(f"📄 Detailed text report saved: {report_path}")

def main():
    """Main function for standalone execution"""
    print("🧬 NBDFinder 2.0 Standalone Analysis")
    print("=" * 50)
    
    # Initialize analyzer
    analyzer = NBDFinderStandalone()
    
    # Example sequences for demonstration
    example_sequences = {
        "G4_Rich_Example": "GGGTTAGGGTTAGGGTTAGGGCCCCAACCCCAACCCCAACCCC",
        "Disease_GAA_Repeat": "GAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAA",
        "Complex_Multi_Motif": "GGGTTAGGGTTAGGGTTAGGGCGCGCGCGCGCGCGCGCGCGCCCCAACCCCAACCCCAACCCCGAAGAAGAAGAAGAAGAAGAATTTAAATTTAAATTTAAA"
    }
    
    print(f"📝 Analyzing {len(example_sequences)} example sequences...")
    
    # Analyze each sequence
    all_results = []
    for seq_name, sequence in example_sequences.items():
        result = analyzer.analyze_sequence(sequence, seq_name)
        all_results.append(result)
        
        # Generate individual report
        output_dir = f"NBDFinder_Results_{seq_name}"
        analyzer.generate_comprehensive_report(result, output_dir)
    
    print("\n🎉 Analysis Complete!")
    print(f"✅ Analyzed {len(all_results)} sequences")
    print(f"📊 Generated comprehensive reports for each sequence")
    print(f"🔬 All 19 motif types assessed with professional visualization")
    
    return all_results

if __name__ == "__main__":
    results = main()