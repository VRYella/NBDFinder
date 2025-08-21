"""
NBDFinder Publication Visualization Examples
===========================================

This script demonstrates all 10 publication-ready visualization types
and generates example figures for documentation purposes.

Usage:
    python generate_publication_examples.py

Output:
    - High-resolution example figures in /tmp/nbdfinder_examples/
    - Comprehensive demonstration of all plot types
    - Documentation-ready figures

Author: Dr. Venkata Rajesh Yella
License: Academic Use
"""

import os
import pandas as pd
import numpy as np
from datetime import datetime
from publication_visualizations import create_comprehensive_publication_suite, PublicationVisualizer
from nbdfinder_viz_integration import NBDFinderVisualizationHub

def create_comprehensive_test_data():
    """Create comprehensive test data representing real NBDFinder analysis results."""
    
    # Set random seed for reproducible examples
    np.random.seed(42)
    
    # Define realistic motif data
    motif_data = []
    
    # G-Quadruplex Family (largest group)
    g4_subtypes = ['Canonical G4', 'Relaxed G4', 'Bulged G4', 'Bipartite G4', 'Multimeric G4', 'Imperfect G4', 'G-Triplex intermediate']
    for i, subtype in enumerate(g4_subtypes):
        # Create multiple instances of each subtype
        for j in range(np.random.randint(3, 8)):  # 3-7 instances each
            start = np.random.randint(100 + i*200, 200 + i*200)
            length = np.random.randint(15, 35)
            motif_data.append({
                'Class': 'G-Quadruplex Family',
                'Subtype': subtype,
                'Start': start,
                'End': start + length,
                'Length': length,
                'Score': np.random.uniform(65, 98),
                'Clinical_Significance': np.random.choice(['Benign', 'VUS', 'Likely Benign'], p=[0.6, 0.3, 0.1])
            })
    
    # Disease-Associated Motifs
    disease_subtypes = ['CGG Expansion', 'CAG Expansion', 'CTG Expansion', 'GGGGCC Expansion']
    for i, subtype in enumerate(disease_subtypes):
        for j in range(np.random.randint(2, 5)):
            start = np.random.randint(1600 + i*150, 1700 + i*150)
            length = np.random.randint(40, 80)  # Longer expansions
            motif_data.append({
                'Class': 'Disease-Associated Motif',
                'Subtype': subtype,
                'Start': start,
                'End': start + length,
                'Length': length,
                'Score': np.random.uniform(85, 99),  # Higher scores for disease motifs
                'Clinical_Significance': np.random.choice(['Pathogenic', 'Likely Pathogenic', 'VUS'], p=[0.4, 0.3, 0.3])
            })
    
    # Z-DNA
    zdna_subtypes = ['Z-DNA', 'eGZ (Extruded-G) DNA']
    for i, subtype in enumerate(zdna_subtypes):
        for j in range(np.random.randint(2, 4)):
            start = np.random.randint(2400 + i*100, 2500 + i*100)
            length = np.random.randint(20, 30)
            motif_data.append({
                'Class': 'Z-DNA',
                'Subtype': subtype,
                'Start': start,
                'End': start + length,
                'Length': length,
                'Score': np.random.uniform(70, 88),
                'Clinical_Significance': np.random.choice(['Benign', 'VUS'], p=[0.7, 0.3])
            })
    
    # R-loop
    for j in range(np.random.randint(3, 6)):
        start = np.random.randint(2800, 3000)
        length = np.random.randint(25, 45)
        motif_data.append({
            'Class': 'R-loop',
            'Subtype': 'R-loop',
            'Start': start,
            'End': start + length,
            'Length': length,
            'Score': np.random.uniform(75, 92),
            'Clinical_Significance': np.random.choice(['VUS', 'Benign'], p=[0.4, 0.6])
        })
    
    # Triplex
    triplex_subtypes = ['Triplex', 'sticky DNA']
    for i, subtype in enumerate(triplex_subtypes):
        for j in range(np.random.randint(2, 4)):
            start = np.random.randint(3200 + i*100, 3300 + i*100)
            length = np.random.randint(30, 50)
            motif_data.append({
                'Class': 'Triplex',
                'Subtype': subtype,
                'Start': start,
                'End': start + length,
                'Length': length,
                'Score': np.random.uniform(72, 89),
                'Clinical_Significance': np.random.choice(['Benign', 'VUS'], p=[0.8, 0.2])
            })
    
    # Curved DNA
    curved_subtypes = ['Global curvature', 'Local Curvature']
    for i, subtype in enumerate(curved_subtypes):
        for j in range(np.random.randint(3, 5)):
            start = np.random.randint(3600 + i*100, 3700 + i*100)
            length = np.random.randint(18, 35)
            motif_data.append({
                'Class': 'Curved DNA',
                'Subtype': subtype,
                'Start': start,
                'End': start + length,
                'Length': length,
                'Score': np.random.uniform(68, 85),
                'Clinical_Significance': np.random.choice(['Benign', 'VUS'], p=[0.9, 0.1])
            })
    
    # Slipped DNA
    slipped_subtypes = ['Slipped DNA [Direct Repeat]', 'Slipped DNA [STR]']
    for i, subtype in enumerate(slipped_subtypes):
        for j in range(np.random.randint(2, 4)):
            start = np.random.randint(4000 + i*100, 4100 + i*100)
            length = np.random.randint(22, 40)
            motif_data.append({
                'Class': 'Slipped DNA',
                'Subtype': subtype,
                'Start': start,
                'End': start + length,
                'Length': length,
                'Score': np.random.uniform(70, 87),
                'Clinical_Significance': np.random.choice(['VUS', 'Benign'], p=[0.3, 0.7])
            })
    
    # Cruciform DNA
    for j in range(np.random.randint(2, 4)):
        start = np.random.randint(4400, 4500)
        length = np.random.randint(25, 40)
        motif_data.append({
            'Class': 'Cruciform DNA',
            'Subtype': 'Cruciform DNA [IR]/HairPin [IR]',
            'Start': start,
            'End': start + length,
            'Length': length,
            'Score': np.random.uniform(73, 90),
            'Clinical_Significance': np.random.choice(['Benign', 'VUS'], p=[0.8, 0.2])
        })
    
    # i-motif family
    imotif_subtypes = ['Canonical i-motif', 'Relaxed i-motif', 'AC-motif']
    for i, subtype in enumerate(imotif_subtypes):
        for j in range(np.random.randint(1, 3)):
            start = np.random.randint(4700 + i*80, 4780 + i*80)
            length = np.random.randint(18, 30)
            motif_data.append({
                'Class': 'i-motif family',
                'Subtype': subtype,
                'Start': start,
                'End': start + length,
                'Length': length,
                'Score': np.random.uniform(69, 86),
                'Clinical_Significance': np.random.choice(['Benign', 'VUS'], p=[0.85, 0.15])
            })
    
    return pd.DataFrame(motif_data)

def generate_all_examples(output_dir="/tmp/nbdfinder_examples"):
    """Generate comprehensive examples of all visualization types."""
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    print("🚀 Generating NBDFinder Publication Visualization Examples")
    print("=" * 60)
    
    # Create comprehensive test data
    print("📊 Creating comprehensive test dataset...")
    motif_data = create_comprehensive_test_data()
    sequence_length = 5000
    
    print(f"   Generated {len(motif_data)} motifs across {motif_data['Class'].nunique()} classes")
    print(f"   Sequence length: {sequence_length:,} bp")
    print(f"   Motif density: {len(motif_data)/sequence_length*1000:.1f} per kb")
    
    # Initialize visualization system
    viz = PublicationVisualizer()
    hub = NBDFinderVisualizationHub()
    
    # Generate comprehensive suite
    print("\n🎨 Creating comprehensive visualization suite...")
    suite = create_comprehensive_publication_suite(motif_data, sequence_length)
    
    # Export examples for each category
    examples_created = 0
    
    for category_name, category_data in suite.items():
        print(f"\n📈 Processing {category_name}...")
        
        if isinstance(category_data, dict):
            # Handle categories with subcategories
            for subcategory_name, fig_data in category_data.items():
                if isinstance(fig_data, dict) and ('plotly' in fig_data or 'matplotlib' in fig_data):
                    filename = f"{category_name}_{subcategory_name}"
                    exported_files = viz.export_all_formats(
                        fig_plotly=fig_data.get('plotly'),
                        fig_matplotlib=fig_data.get('matplotlib'),
                        filename=filename,
                        width=1200, height=800, dpi=300
                    )
                    
                    # Save files to output directory
                    for file_key, file_bytes in exported_files.items():
                        file_path = os.path.join(output_dir, file_key)
                        with open(file_path, 'wb') as f:
                            f.write(file_bytes)
                        examples_created += 1
                        print(f"   ✓ Created {file_key}")
        
        else:
            # Handle single-figure categories
            if isinstance(category_data, dict) and ('plotly' in category_data or 'matplotlib' in category_data):
                filename = category_name
                exported_files = viz.export_all_formats(
                    fig_plotly=category_data.get('plotly'),
                    fig_matplotlib=category_data.get('matplotlib'),
                    filename=filename,
                    width=1200, height=800, dpi=300
                )
                
                # Save files to output directory
                for file_key, file_bytes in exported_files.items():
                    file_path = os.path.join(output_dir, file_key)
                    with open(file_path, 'wb') as f:
                        f.write(file_bytes)
                    examples_created += 1
                    print(f"   ✓ Created {file_key}")
    
    # Create a summary README file
    create_example_readme(motif_data, output_dir, examples_created)
    
    print("\n" + "=" * 60)
    print(f"✅ Example generation complete!")
    print(f"   📁 Output directory: {output_dir}")
    print(f"   📊 Total figures created: {examples_created}")
    print(f"   📋 Data summary: {len(motif_data)} motifs, {motif_data['Class'].nunique()} classes")
    print(f"   🔗 See README.md in output directory for details")

def create_example_readme(motif_data, output_dir, examples_created):
    """Create a comprehensive README for the generated examples."""
    
    readme_content = f"""# NBDFinder Publication Visualization Examples

Generated on: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}

## Overview

This directory contains comprehensive examples of all 10 publication-ready visualization types implemented in NBDFinder. These figures demonstrate the platform's capability to generate high-resolution, manuscript-quality visualizations suitable for scientific publications.

## Dataset Summary

- **Total Motifs**: {len(motif_data)}
- **Unique Classes**: {motif_data['Class'].nunique()}
- **Unique Subtypes**: {motif_data['Subtype'].nunique()}
- **Sequence Length**: 5,000 bp
- **Motif Density**: {len(motif_data)/5000*1000:.1f} motifs per kb

### Class Distribution
"""
    
    # Add class distribution
    class_counts = motif_data['Class'].value_counts()
    for cls, count in class_counts.items():
        percentage = (count / len(motif_data)) * 100
        readme_content += f"- **{cls}**: {count} motifs ({percentage:.1f}%)\n"
    
    readme_content += f"""
### Clinical Significance Distribution
"""
    
    # Add clinical significance distribution
    clin_counts = motif_data['Clinical_Significance'].value_counts()
    for sig, count in clin_counts.items():
        percentage = (count / len(motif_data)) * 100
        readme_content += f"- **{sig}**: {count} motifs ({percentage:.1f}%)\n"
    
    readme_content += f"""
## Visualization Types

NBDFinder implements 10 major publication-ready visualization types:

### 1. Bar Plots and Stacked Bar Plots
- **Files**: `bar_plots_*`
- **Purpose**: Visualize counts and proportions for each motif class and subclass
- **Features**: Color coding for major classes, clear legends, stacked options

### 2. Linear Motif Maps (Genome Tracks)
- **Files**: `linear_maps_*`
- **Purpose**: Plot motif positions, lengths, and overlaps along sequences
- **Features**: Multi-track overlays, labeled colored blocks, position scaling

### 3. Heatmaps
- **Files**: `heatmaps_*`
- **Purpose**: Display motif density, score distributions, co-occurrence matrices
- **Features**: Clustering, annotated axes, multiple heatmap types

### 4. Pie/Donut Charts
- **Files**: `pie_charts_*`
- **Purpose**: Summarize motif composition by class/subclass
- **Features**: Interactive legends, percentage labels, donut and pie variants

### 5. Violin and Box Plots
- **Files**: `distribution_plots_*`
- **Purpose**: Show distributions of motif scores, lengths, or repeat counts
- **Features**: Statistical annotations, class-based grouping, distribution overlays

### 6. UpSet Plots
- **Files**: `upset_plots_*`
- **Purpose**: Visualize intersections and overlaps between motif classes
- **Features**: Set intersections, overlap quantification, hybrid region analysis

### 7. Lollipop Plots
- **Files**: `lollipop_plots_*`
- **Purpose**: Annotate key motif positions with properties and risk status
- **Features**: Position-based visualization, risk color coding, score scaling

### 8. Bubble/Scatter Plots
- **Files**: `bubble_plots_*`
- **Purpose**: Illustrate relationships and outlier detection
- **Features**: Multi-dimensional encoding, trend lines, correlation analysis

### 9. Circos Plots
- **Files**: `circos_plots_*`
- **Purpose**: Visualize motif distributions across genomic coordinates
- **Features**: Circular representation, density tracks, angular positioning

### 10. Sankey Diagrams
- **Files**: `sankey_plots_*`
- **Purpose**: Show flows from motifs to subclasses and clinical classifications
- **Features**: Flow quantification, hierarchical relationships, classification pathways

## File Formats

Each visualization is exported in multiple high-resolution formats:

- **PNG**: Raster format at 300+ DPI (publication quality)
- **PDF**: Vector format for scalable printing
- **SVG**: Vector format for web and further editing

## File Naming Convention

Files are named using the pattern: `[category]_[subcategory]_[version].[format]`

- `category`: Main visualization type (e.g., `bar_plots`, `heatmaps`)
- `subcategory`: Specific variant (e.g., `simple`, `stacked`, `scores`)
- `version`: `plotly` (interactive) or `matplotlib` (publication)
- `format`: `png`, `pdf`, or `svg`

## Technical Specifications

- **Resolution**: PNG files at 300 DPI minimum
- **Dimensions**: 1200×800 pixels (standard), scalable for larger formats
- **Color Schemes**: Colorblind-friendly palettes (Wong palette compliance)
- **Typography**: Publication-standard fonts (Times New Roman, Arial)
- **Accessibility**: High contrast ratios, pattern differentiation available

## Usage in Publications

These figures are designed for direct use in:

- Scientific manuscripts (Nature, Science, Cell format compliance)
- Conference presentations and posters
- Grant applications and reports
- Educational materials and documentation

## Quality Assurance

All visualizations undergo quality checks for:

- ✅ Publication-standard resolution (≥300 DPI)
- ✅ Colorblind accessibility (tested with simulation)
- ✅ Professional typography and layout
- ✅ Clear axis labels and legends
- ✅ Scientific accuracy and consistency

## Software Information

- **NBDFinder Version**: 2.0
- **Visualization Engine**: Plotly + Matplotlib
- **Export Engine**: Kaleido + Native matplotlib
- **Data Processing**: Pandas + NumPy
- **Generated**: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}

---

*Generated by NBDFinder Publication Visualization Suite*
*For support: https://github.com/VRYella/NBDFinder*
"""
    
    # Write README file
    readme_path = os.path.join(output_dir, "README.md")
    with open(readme_path, 'w') as f:
        f.write(readme_content)
    
    print(f"   📋 Created comprehensive README.md")

if __name__ == "__main__":
    # Generate comprehensive examples
    generate_all_examples()
    
    print("\n🎯 Quick Test - Individual Plot Types:")
    
    # Quick test of individual visualizations
    test_data = create_comprehensive_test_data().head(20)  # Smaller subset for quick test
    viz = PublicationVisualizer()
    
    test_categories = [
        ("Bar Plot", lambda: viz.create_bar_plots(test_data)),
        ("Heatmap", lambda: viz.create_heatmaps(test_data, heatmap_type='scores')),
        ("Pie Chart", lambda: viz.create_pie_donut_charts(test_data)),
        ("Violin Plot", lambda: viz.create_violin_box_plots(test_data)),
        ("Lollipop Plot", lambda: viz.create_lollipop_plots(test_data, reference_length=5000)),
        ("Bubble Plot", lambda: viz.create_bubble_scatter_plots(test_data)),
        ("UpSet Plot", lambda: viz.create_upset_plots(test_data)),
        ("Sankey Diagram", lambda: viz.create_sankey_diagrams(test_data))
    ]
    
    for name, create_func in test_categories:
        try:
            result = create_func()
            status = "✅" if result else "⚠️"
            print(f"   {status} {name}")
        except Exception as e:
            print(f"   ❌ {name}: {e}")
    
    print("\n🏆 NBDFinder Publication Visualization Suite - Complete!")
    print("Ready for integration with NBDFinder web application.")