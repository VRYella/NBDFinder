#!/usr/bin/env python3
"""
Complete the Disease Expansion Analysis - Generate visualizations and final report
"""

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import json
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

def complete_analysis():
    """Complete the analysis with visualizations and report"""
    print("🎨 Completing Disease Expansion Analysis...")
    
    results_dir = Path("disease_expansion_analysis_results")
    summary_dir = results_dir / "summary_reports"
    plots_dir = results_dir / "visualizations"
    plots_dir.mkdir(exist_ok=True)
    
    # Load the Excel data
    excel_file = summary_dir / "comprehensive_disease_expansion_analysis.xlsx"
    if not excel_file.exists():
        print("❌ Excel file not found")
        return
    
    # Read summary data
    summary_df = pd.read_excel(excel_file, sheet_name='Summary')
    print(f"📊 Loaded data for {len(summary_df)} genes")
    
    # Generate statistics
    stats = {
        'total_sequences': len(summary_df),
        'total_motifs': int(summary_df['Total_Motifs'].sum()),
        'total_disease_motifs': int(summary_df['Disease_Motifs'].sum()),
        'avg_motifs_per_gene': float(summary_df['Total_Motifs'].mean()),
        'avg_disease_motifs_per_gene': float(summary_df['Disease_Motifs'].mean()),
        'avg_sequence_length': float(summary_df['Sequence_Length'].mean()),
        'avg_gc_content': float(summary_df['GC_Content'].mean()),
        'genes_with_disease_motifs': int((summary_df['Disease_Motifs'] > 0).sum()),
        'max_motifs_gene': summary_df.loc[summary_df['Total_Motifs'].idxmax(), 'Gene_Symbol'],
        'max_motifs_count': int(summary_df['Total_Motifs'].max()),
    }
    
    # Save corrected statistics
    stats_file = summary_dir / "analysis_statistics.json"
    with open(stats_file, 'w') as f:
        json.dump(stats, f, indent=2)
    
    print(f"📈 Analysis Statistics:")
    print(f"   🧬 Total sequences analyzed: {stats['total_sequences']}")
    print(f"   🔍 Total motifs found: {stats['total_motifs']}")
    print(f"   🦠 Disease-specific motifs: {stats['total_disease_motifs']}")
    print(f"   📏 Average motifs per gene: {stats['avg_motifs_per_gene']:.1f}")
    print(f"   🎯 Genes with disease motifs: {stats['genes_with_disease_motifs']}")
    print(f"   🏆 Gene with most motifs: {stats['max_motifs_gene']} ({stats['max_motifs_count']} motifs)")
    
    # Create visualizations
    print("\n🎨 Creating visualizations...")
    
    # Set up plotting style
    plt.style.use('default')
    sns.set_palette("husl")
    
    # 1. Distribution of motifs per gene
    plt.figure(figsize=(12, 8))
    plt.subplot(2, 2, 1)
    plt.hist(summary_df['Total_Motifs'], bins=15, alpha=0.7, edgecolor='black', color='skyblue')
    plt.xlabel('Total Motifs per Gene')
    plt.ylabel('Number of Genes')
    plt.title('Distribution of Motifs per Gene')
    plt.grid(True, alpha=0.3)
    
    # 2. GC Content vs Motif Density
    plt.subplot(2, 2, 2)
    plt.scatter(summary_df['GC_Content'], summary_df['Motif_Density'], 
               alpha=0.6, s=50, color='coral')
    plt.xlabel('GC Content (%)')
    plt.ylabel('Motif Density (per kb)')
    plt.title('GC Content vs Motif Density')
    plt.grid(True, alpha=0.3)
    
    # 3. Disease motifs distribution
    plt.subplot(2, 2, 3)
    disease_genes = summary_df[summary_df['Disease_Motifs'] > 0]
    if not disease_genes.empty:
        plt.hist(disease_genes['Disease_Motifs'], bins=10, alpha=0.7, 
                color='lightcoral', edgecolor='black')
        plt.xlabel('Disease Motifs per Gene')
        plt.ylabel('Number of Genes')
        plt.title(f'Disease Motifs Distribution\n({len(disease_genes)} genes with disease motifs)')
    else:
        plt.text(0.5, 0.5, 'No Disease Motifs Found', 
                transform=plt.gca().transAxes, ha='center', va='center')
        plt.title('Distribution of Disease Motifs')
    plt.grid(True, alpha=0.3)
    
    # 4. Top genes by motif count
    plt.subplot(2, 2, 4)
    top_genes = summary_df.nlargest(10, 'Total_Motifs')
    plt.barh(range(len(top_genes)), top_genes['Total_Motifs'], color='lightgreen')
    plt.yticks(range(len(top_genes)), top_genes['Gene_Symbol'])
    plt.xlabel('Total Motifs')
    plt.title('Top 10 Genes by Motif Count')
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(plots_dir / 'overview_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Summary table visualization
    plt.figure(figsize=(12, 6))
    
    # Create a summary table
    table_data = [
        ['Total Sequences Analyzed', stats['total_sequences']],
        ['Total Motifs Found', stats['total_motifs']],
        ['Disease-Specific Motifs', stats['total_disease_motifs']],
        ['Average Motifs per Gene', f"{stats['avg_motifs_per_gene']:.1f}"],
        ['Genes with Disease Motifs', stats['genes_with_disease_motifs']],
        ['Gene with Most Motifs', f"{stats['max_motifs_gene']} ({stats['max_motifs_count']})"],
        ['Average Sequence Length', f"{stats['avg_sequence_length']:.0f} bp"],
        ['Average GC Content', f"{stats['avg_gc_content']:.1f}%"]
    ]
    
    plt.axis('tight')
    plt.axis('off')
    table = plt.table(cellText=table_data, 
                      colLabels=['Metric', 'Value'],
                      cellLoc='left',
                      loc='center',
                      bbox=[0, 0, 1, 1])
    table.auto_set_font_size(False)
    table.set_fontsize(12)
    table.scale(1.2, 1.5)
    
    # Style the table
    for (i, j), cell in table.get_celld().items():
        if i == 0:  # Header row
            cell.set_facecolor('#4CAF50')
            cell.set_text_props(weight='bold', color='white')
        else:
            cell.set_facecolor('#f0f0f0' if i % 2 == 0 else '#ffffff')
        cell.set_edgecolor('gray')
    
    plt.title('Disease Expansion Analysis Summary', fontsize=16, fontweight='bold', pad=20)
    plt.savefig(plots_dir / 'summary_table.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    # Generate final report
    print("\n📝 Generating final report...")
    
    report_file = summary_dir / "analysis_report.md"
    
    with open(report_file, 'w') as f:
        f.write("# Disease Expansion Loci Analysis Report\n\n")
        f.write(f"**Analysis completed successfully on 20 representative disease-related gene sequences**\n\n")
        
        f.write("## Executive Summary\n\n")
        f.write(f"This comprehensive analysis examined {stats['total_sequences']} disease-associated gene sequences ")
        f.write(f"using the NBDFinder platform to detect non-B DNA structural motifs. ")
        f.write(f"The analysis identified **{stats['total_motifs']} total motifs** across 10 structural classes, ")
        f.write(f"including **{stats['total_disease_motifs']} disease-specific motifs** in ")
        f.write(f"**{stats['genes_with_disease_motifs']} genes**.\n\n")
        
        f.write("## Key Findings\n\n")
        f.write(f"- **Average motifs per gene:** {stats['avg_motifs_per_gene']:.1f}\n")
        f.write(f"- **Most motif-rich gene:** {stats['max_motifs_gene']} with {stats['max_motifs_count']} motifs\n")
        f.write(f"- **Disease motif prevalence:** {stats['genes_with_disease_motifs']}/{stats['total_sequences']} genes ({stats['genes_with_disease_motifs']/stats['total_sequences']*100:.1f}%)\n")
        f.write(f"- **Average sequence characteristics:** {stats['avg_sequence_length']:.0f} bp, {stats['avg_gc_content']:.1f}% GC\n\n")
        
        f.write("## Top Genes by Motif Content\n\n")
        for _, row in top_genes.head().iterrows():
            f.write(f"1. **{row['Gene_Symbol']}**: {row['Total_Motifs']} motifs ")
            f.write(f"({row['Disease_Motifs']} disease-specific)\n")
        
        f.write("\n## Analysis Coverage\n\n")
        f.write("This analysis demonstrates the comprehensive capabilities of NBDFinder for:\n")
        f.write("- **Structural motif detection** across 10 classes (Curved DNA, Slipped DNA, Cruciform, R-loops, etc.)\n")
        f.write("- **Disease-specific repeat analysis** with clinical significance assessment\n")
        f.write("- **Comparative genomics** across multiple disease-associated genes\n")
        f.write("- **High-throughput processing** with detailed individual and summary reports\n\n")
        
        f.write("## Files Generated\n\n")
        f.write("- **Individual analyses:** `individual_sequences/` - Detailed results for each gene\n")
        f.write("- **Comprehensive report:** `comprehensive_disease_expansion_analysis.xlsx`\n")
        f.write("- **Visualizations:** `visualizations/` - Analysis summary plots\n")
        f.write("- **Statistics:** `analysis_statistics.json` - Quantitative summary\n\n")
        
        f.write("## Methodology\n\n")
        f.write("- **Platform:** NBDFinder v2024 with comprehensive motif detection algorithms\n")
        f.write("- **Input:** 153 disease-associated gene sequences from repeat_expansion_loci_annotated.fa\n")
        f.write("- **Analysis scope:** 20 representative sequences for demonstration\n")
        f.write("- **Output formats:** CSV, Excel, JSON, PNG visualizations\n")
    
    print(f"  ✅ Report saved to {report_file}")
    print(f"  ✅ Visualizations saved to {plots_dir}")
    print(f"\n🎉 Analysis completed successfully!")
    print(f"📁 All results available in: {results_dir}")

if __name__ == "__main__":
    complete_analysis()